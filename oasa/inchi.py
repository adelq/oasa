#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003 Beda Kosata <beda@zirael.org>

#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     Complete text of GNU GPL can be found in the file gpl.txt in the
#     main directory of the program

#--------------------------------------------------------------------------

from plugin import plugin
from molecule import molecule
import periodic_table as pt
import molfile
import misc

import re
from sets import Set
import operator
import time
import xml.dom.minidom as dom
import dom_extensions
import string
import os.path
import tempfile
import popen2
import coords_generator
from oasa_exceptions import oasa_not_implemented_error
from tempfile import mkstemp


class inchi( plugin):

  name = "inchi"
  read = 1
  write = 0

  proton_donors = ['O','F','Cl','Br','I','N','S']
  proton_acceptors = ['P','N','S','O']
  

  def __init__( self, structure=None):
    self.structure = structure
    self._no_possibility_to_improve = False
    self.cleanup()

  def cleanup( self):
    self.forced_charge = 0
    self._protonation_dealt_with_already = 0
    self._added_hs = Set()

  def set_structure( self, structure):
    self.structure = structure

  def get_structure( self):
    return self.structure



  def split_layers( self, text):
    # try if we have a xml document in hand
    try:
      doc = dom.parseString( text)
    except:
      # it seems not to be the case
      return re.split( "/", text.strip())
    structs = doc.getElementsByTagName( 'structure')
    if not structs:
      raise "no structures found in xml string %s" % text
    struct = structs[0]
    layers = []
    for name in ('version','formula','connections','H'):
      try:
        layer = struct.getElementsByTagName( name)[0]
      except IndexError:
        raise "no '%s' tag found in xml string" % name
      layers.append( dom_extensions.getTextFromElement( layer))
    return layers




  def get_layer( self, prefix):
    for l in self.layers:
      if l.startswith( prefix):
        return l[1:]




  def read_inchi( self, text):
    self.structure = molecule()
    self.layers = self.split_layers( text)
    # version support (very crude)
    self.version = self._check_version( self.layers[0])
    if not self.version:
      raise ValueError, "this version of INChI is not supported - %s" % self.layers[0]
    
    self.hs_in_hydrogen_layer = self.get_number_of_hydrogens_in_hydrogen_layer()
    self.read_sum_layer()
    self.read_connectivity_layer()
    repeat = True
    run = 0
    # we have to repeat this step in order to find the right positioning of movable hydrogens
    while repeat and not self._no_possibility_to_improve:
      # cleanup
      for h in self._added_hs:
        self.structure.remove_vertex( h)
      for b in self.structure.edges:
        b.order = 1
      for v in self.structure.vertices:
        v.symbol = v.symbol
        v.charge = 0
      self.cleanup()

      # the code itself
      run += 1
      assert run < 20
      self.process_forced_charges()
      self.read_hydrogen_layer( run=run)
      self.read_charge_layer()
      self.read_p_layer()
      #[a.raise_valency_to_senseful_value() for a in self.structure.vertices]
      # if there is no possibility to improve via the hydrogen positioning we must try the retry here
      self.structure.add_missing_bond_orders()

      # here we check out if the molecule seems ok
      if not filter( None, [v.free_valency for v in self.structure.vertices]) and \
             not filter( None, [not v.order for v in self.structure.edges]):
        repeat = False
      else:
        pass
    if repeat and self._no_possibility_to_improve:
      pass
      #print "**error warning"




  def read_sum_layer( self):
    if "." in self.layers[1]:
      raise oasa_not_implemented_error( "INChI", "multiple compound systems are not supported by the library")

    form = pt.formula_dict( self.layers[1])
    processed_hs = 0 #for diborane and similar compounds we must process some Hs here
    for k in form.sorted_keys():
      for i in range( form[k]):
        if k == 'H':
          # we want to process only the Hs that are not in the h-layer
          if processed_hs >= form[k] - self.hs_in_hydrogen_layer:
            continue
          else:
            processed_hs += 1
        a = self.structure.create_vertex()
        a.symbol = k
        self.structure.add_vertex( a)



  def read_connectivity_layer( self):
    layer = self.get_layer( "c")
    chunks = re.split( "([0-9]*)", layer)
    chunks = filter( None, chunks)
    chunks = filter( lambda x: x!='-', chunks)
    last_atom = None
    bracket_openings = []
    for c in chunks:
      if c == '(':
        bracket_openings.append( last_atom)
      elif c == ')':
        last_atom = bracket_openings.pop(-1)
      elif c == ",":
        last_atom = bracket_openings.pop(-1)
        bracket_openings.append( last_atom)
      else:
        try:
          i = int( c)
        except:
          print "should not happen", c
          continue
        # atom
        if last_atom:
          self.structure.add_edge( last_atom-1, i-1)
        last_atom = i


  def read_hydrogen_layer( self, run=0):
    layer = self.get_layer( "h")
    # check presence of the layer
    if not layer:
      return

    re_for_brackets = "\([H\d,\-]+?\)"
    brackets = re.findall( re_for_brackets, layer)
    for bracket in brackets:
      self._process_moving_hydrogen( bracket[1:-1], run=run)

    # we improve only in the _process_moving_hydrogen so if it was not called there is no possibility for improvement
    if not self._added_hs:
      self._no_possibility_to_improve = True
      
    layer = re.sub( re_for_brackets, "", layer)  # clean the brackets out

    for vs, num in self._parse_h_layer( layer):
      for v in vs:
        for j in range( num):
          h = self.structure.create_vertex()
          h.symbol = 'H'
          self.structure.add_vertex( h)
          self.structure.add_edge( v-1, h)
          self._added_hs.add( h)



  def read_charge_layer( self):
    layer = self.get_layer( "q")
    if not layer:
      return

    charge = int( layer[1:]) - self.forced_charge

    change = True
    # at first we try to put the charge on atoms with negative free_valency
    steps = 0
    while charge and change:
      steps += 1
      assert steps < 10
      change = False
      for v in self.structure.vertices:
        if v.free_valency == -1:
          if charge > 0:
            v.charge = -v.free_valency
          else:
            v.charge = v.free_valency
          charge -= v.charge
          change = True
          break

    steps = 0
    change = True
    # then we try to put the charge on orphan atoms
    while charge and change:
      steps += 1
      assert steps < 10

      change = False
      for v in self.structure.vertices:
        if v.free_valency and sum( [n.free_valency for n in v.neighbors]) < v.free_valency:
          charge = self._valency_to_charge( v, charge)
          change = True
          break

    # this could be rewritten to put charges to places where a segment with odd number of free_valencies is
    # this would be more general than counting the free_valencies for the whole molecule
    free_valencies = sum( [a.free_valency for a in self.structure.vertices])
    if free_valencies % 2:
      # odd number of free_valencies means we have to put the charge somewhere to make a free_valency there
      change = True
      steps = 0

      while charge and change:
        change = False
        steps += 1
        assert steps < 10

        for v in self.structure.vertices:
          if v.symbol != 'C' and not v.free_valency and filter( None, [n.free_valency for n in v.neighbors]):
            v.charge = charge
            charge = 0
            change = True
            break

    # then we do the others, at first trying heteroatoms
    steps = 0
    change = True
    while charge and change:
      change = False
      steps += 1
      assert steps < 10

      for v in self.structure.vertices:
        if v.free_valency and v.symbol != 'C':
          charge = self._valency_to_charge( v, charge)
          change = True
          break


    # then we try the carbon atoms as well
    steps = 0
    change = True
    while charge and change:
      steps += 1
      assert steps < 10

      change = False
      for v in self.structure.vertices:
        if v.free_valency:
          charge = self._valency_to_charge( v, charge)
          change = True
          break

    if charge:
      raise oasa_exceptions.oasa_inchi_error( "The molecular charge could not be allocated to any atom (%d)." % charge)



  def read_p_layer( self):
    layer = self.get_layer( "p")
    if not layer:
      return

    p = int( layer[1:]) - self._protonation_dealt_with_already
    old_p = p

    # at first add protons to negative charges
    for v in self.structure.vertices:
      if p > 0:
        if v.charge < 0:
          h = self.structure.create_vertex()
          h.symbol = 'H'
          self.structure.add_vertex( h)
          self.structure.add_edge( v, h)
          p -= 1
          v.charge += 1
          self._protonation_dealt_with_already += 1

    if self._protonation_dealt_with_already:
      # p was solved above
      dp = self._protonation_dealt_with_already - sum( [v.charge for v in self.structure.vertices])
      if dp > 0:
        # there were no forced charges, so we are not in a zwittrion state, huh, who does not understand that?!
        for v in self.structure.vertices:
          if dp > 0 and v.symbol in self.proton_acceptors:
            v.charge += 1
            dp -= 1

    # then come the other additional protones
    while p:
      old_p = p
      for v in self.structure.vertices:
        if p > 0:
          if v.symbol in self.proton_acceptors:
            h = self.structure.create_vertex()
            h.symbol = 'H'
            self.structure.add_vertex( h)
            self.structure.add_edge( v, h)
            self._added_hs.add( h)
            p -= 1
            v.charge += 1
            break
        elif p < 0:
          if v.symbol in self.proton_donors:
            hs = [n for n in v.neighbors if n.symbol=='H']
            if hs:
              h = hs[0]
              self.structure.remove_vertex( h)
              p += 1
              v.charge -= 1
              break
      assert p < old_p




  def _valency_to_charge( self, v, charge):
    """sets charge of the atom so that it has minimal free_valency,
    returns the 'unused' charge in case it is not comsumed whole"""
    ch = misc.signum( charge) * min( abs( charge), v.free_valency)
    if charge < 0:
      charge += ch
      v.charge = -ch
    else:
      charge -= ch
      v.charge = ch
    return charge
        

  def get_number_of_hydrogens_in_hydrogen_layer( self):
    # version check
    layer = self.get_layer( "h")
    if not layer:
      return 0

    # check if we can handle it
    if "*" in layer or ";" in layer:
      raise oasa_not_implemented_error( "INChI", "multiple compound systems are not supported by the library")

    ret = 0

    re_for_brackets = "\([H\d,\-]+?\)"
    brackets = re.findall( re_for_brackets, layer)
    for bracket in brackets:
      ret += self._get_hs_in_moving_hydrogen( bracket[1:-1])
    layer = re.sub( re_for_brackets, "", layer)  # clean the brackets out

    for vs, num in self._parse_h_layer( layer):
      ret += len( vs) * num

    return ret


  def _get_hs_in_moving_hydrogen( self, chunk):
    chks = chunk.split( ',')
    if len( chks[0]) > 1:
      if chks[0][1:].endswith( "-"):
        hs = int( chks[0][1:-1] or 1)  # number of atoms 
      else:
        hs = int( chks[0][1:])  # number of atoms
    else:
      hs = 1
    return hs
    


  def _check_version( self, ver):
    """returns a tuple of two version numbers"""
    if "Beta" in ver:
      return None
    v = ver.strip(string.ascii_lowercase+string.ascii_uppercase)
    if "." in v:
      return v.split('.')
    else:
      return v, 0


  def _process_moving_hydrogen( self, chunk, run=0):
    chks = chunk.split( ',')
    if len( chks[0]) > 1:
      if chks[0][1:].endswith("-"):
        hs = int( chks[0][1:-1] or 1) + 1  # number of atoms 
        self._protonation_dealt_with_already += 1
      else:
        hs = int( chks[0][1:])  # number of atoms
    else:
      hs = 1

    vertices = [self.structure.vertices[ int( i)-1] for i in chks[1:]]
    vs = vertices
    take = hs

    # here we generate still shorter combinations of vertices that receive hydrogens
    # because of applying them in circles in the next step, we simulate the case
    # where more than one hydrogen is put to one atom
    while (run > misc.x_over_y( len( vs), take)) and take:
      run -= misc.x_over_y( len( vs), take)
      take -= 1
    if not take:
      self._no_possibility_to_improve = True
      return []
    variations = misc.gen_variations( vs, take)
    for i in range( run):
      vs = variations.next()
        
    while hs:
      for v in vs:
        h = self.structure.create_vertex()
        h.symbol = 'H'
        self.structure.add_vertex( h)
        self.structure.add_edge( v, h)
        self._added_hs.add( h)
        hs -= 1
        if not hs:
          break



  def _split_h_layer( self, layer):
    was_h = False
    chunks = []
    chunk = ""
    for c in layer:
      if c == 'H':
        was_h = True
        chunk += c
      elif c == "," and was_h:
        was_h = False
        chunks.append( chunk)
        chunk = ""
      else:
        chunk += c
    if chunk:
      chunks.append( chunk)
    return chunks


  def _parse_h_layer( self, layer):
    chunks = self._split_h_layer( layer)
    for chunk in chunks:
      head, tail = chunk.split( 'H')
      num_h = tail and int( tail) or 1
      vertices = []
      for p in head.split( ","):
        if "-" in p:
          a, b = map( int, p.split("-"))
          vertices.extend( range( a, b+1))
        else:
          vertices.append( int( p))
      yield vertices, num_h


  def process_forced_charges( self):
    """this marks the charges that are forced by the connectivity and thus helps
    process zwitrions"""
    forced_charge = 0
    for v in self.structure.vertices:
      if v.symbol in ("N",) and v.free_valency == -1:
        v.charge = 1
        forced_charge += 1
        
    self.forced_charge = forced_charge



def generate_inchi( m):
  program = "/home/beda/inchi/cInChI-1"

  f, name = mkstemp( text=True)
  os.close( f)
  
  file = open( name, "w")
  molfile.mol_to_file( m, file)
  file.close()

  f, in_name = mkstemp( text=True)
  os.close( f)

  if os.name == 'nt':
    options = "/AUXNONE"
  else:
    options = "-AUXNONE"
  command = ' '.join( (program, name, in_name, options))
  popen = popen2.Popen4( command)
  exit_code = popen.wait()
  #exit_code = os.spawnv( os.P_WAIT, program, (program, name, in_name, options))

  if exit_code == 0 or exit_code != 0:
    in_file = open( in_name, 'r')
    # go to the last line
    i = 0
    for line in in_file.readlines():
      i += 1
    if i > 1:
      out = ( line[6:].strip())
    else:
      # single line file is the only way how to determine it has crashed
      out = ""
    in_file.close()
  else:
    out = ''

  os.remove( in_name)
  os.remove( name)

  return out

    




##################################################
# MODULE INTERFACE

reads_text = 1
reads_files = 1
writes_text = 1
writes_files = 1

def text_to_mol( text, include_hydrogens=True, mark_aromatic_bonds=False, calc_coords=1):
  inc = inchi()
  inc.read_inchi( text)
  mol = inc.structure
  #mol.add_missing_bond_orders()
  if not include_hydrogens:
    mol.remove_all_hydrogens()
  if calc_coords:
    coords_generator.calculate_coords( mol)
  if mark_aromatic_bonds:
    mol.mark_aromatic_bonds()
  return mol


def file_to_mol( f):
  return text_to_mol( f.read())


def mol_to_text( mol):
  return generate_inchi( mol)


def mol_to_file( mol, f):
  f.write( mol_to_text( mol))



#
##################################################




##################################################
# DEMO

if __name__ == '__main__':

  import smiles

  def main( text, cycles):
    t1 = time.time()
    for jj in range( cycles):
      mol = text_to_mol( text)
      print map( str, [b for b in mol.bonds if b.order == 0])
      print "  smiles: ", smiles.mol_to_text( mol)
      print "  inchi:  ", mol_to_text( mol)
      print "  charge: ", sum( [a.charge for a in mol.vertices])
    t1 = time.time() - t1
    print 'time per cycle', round( 1000*t1/cycles, 2), 'ms'

  repeat = 3
  inch = "1/C4H6N2O3S/c1-5-2-3-6(4-5)10(7,8)9/h2-4H,1H3/p+1"
  print "oasa::INCHI DEMO"
  print "converting following inchi into smiles (%d times)" % repeat
  print "  inchi: %s" % inch
  
  #import profile
  #profile.run( 'main( inch, repeat)')
  main( inch, repeat)

# DEMO END
##################################################



##################################################
## TODO

## rename layer to sublayer


##################################################
