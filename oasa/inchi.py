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
from atom import atom
from bond import bond
import periodic_table as pt
import molfile

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


class inchi( plugin):

  name = "inchi"
  read = 1
  write = 0

  def __init__( self, structure=None):
    self.structure = structure

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


  def read_inchi( self, text):
    self.structure = molecule()
    layers = self.split_layers( text)
    # version support (very crude)
    self.version = self._get_version( layers[0])
    hs_in_hydrogen_layer = self.get_number_of_hydrogens_in_hydrogen_layer( layers[3])
    self.read_sum_layer( layers[1], hs_in_hydrogen_layer)
    self.read_connectivity_layer( layers[2])
    self.read_hydrogen_layer( layers[3])
    self.structure.add_missing_bond_orders()
    self.structure.remove_all_hydrogens()


  def read_sum_layer( self, layer, hs_in_hydrogen_layer):
    form = pt.formula_dict( layer)
    processed_hs = 0 #for diborane and similar compounds we must process some Hs here
    for k in form.sorted_keys():
      for i in range( form[k]):
        if k == 'H':
          # we want to process only the Hs that are not in the h-layer
          if processed_hs >= form[k] - hs_in_hydrogen_layer:
            continue
          else:
            processed_hs += 1
        self.structure.add_vertex( atom( symbol=k))



  def read_connectivity_layer( self, layer):
    # version check
    if self.version[0] >= 1 and self.version[1] > 11:
      if layer[0] != 'c':
        print "missing layers specifier at the begining of the string"
      else:
        layer = layer[1:]
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
      else:
        try:
          i = int( c)
        except:
          print "should not happen"
          continue
        # atom
        if last_atom:
          self.structure.add_edge( last_atom-1, i-1, e=bond())
        last_atom = i

  def read_hydrogen_layer( self, layer):
    # version check
    if self.version[0] >= 1 and self.version[1] > 11:
      if layer[0] != 'h':
        print "missing layers specifier at the begining of the string"
      else:
        layer = layer[1:]

    re_for_brackets = "\([H\d,]+?\)"
    brackets = re.findall( re_for_brackets, layer)
    for bracket in brackets:
      self._process_moving_hydrogen( bracket[1:-1])
    layer = re.sub( re_for_brackets, "", layer)  # clean the brackets out

    chunks = re.split( ",", layer)
    chunks = filter( None, chunks)
    hatoms = [a for a in self.structure.vertices if a.symbol=='H']
    for c in chunks:
      as, ns = re.split( "H", c)
      if '-' in as:
        a1, a2 = re.split( '-', as)
      else:
        a1 = as
        a2 = as
      a1 = int( a1)
      a2 = int( a2)
      ns = ns and int(ns) or 1
      for i in range( a1, a2+1):
        for j in range( ns):
          h = self.structure.add_vertex( atom( symbol='H'))
          self.structure.add_edge( i-1, h, e=bond())
        

  def get_number_of_hydrogens_in_hydrogen_layer( self, layer):
    # version check
    if self.version[0] >= 1 and self.version[1] > 11:
      if layer[0] != 'h':
        print "missing layers specifier at the begining of the string"
      else:
        layer = layer[1:]

    ret = 0

    re_for_brackets = "\([H\d,]+?\)"
    brackets = re.findall( re_for_brackets, layer)
    for bracket in brackets:
      ret += self._get_hs_in_moving_hydrogen( bracket[1:-1])
    layer = re.sub( re_for_brackets, "", layer)  # clean the brackets out

    chunks = re.split( ",", layer)
    chunks = filter( None, chunks)
    for c in chunks:
      as, ns = re.split( "H", c)
      if '-' in as:
        a1, a2 = re.split( '-', as)
      else:
        a1 = as
        a2 = as
      a1 = int( a1)
      a2 = int( a2)
      ns = ns and int(ns) or 1
      for i in range( a1, a2+1):
        ret += ns

    return ret


  def _get_hs_in_moving_hydrogen( self, chunk):
    chks = chunk.split( ',')
    if len( chks[0]) > 1:
      hs = int( chks[0][1:])  # number of atoms
    else:
      hs = 1
    return hs
    


  def _get_version( self, ver):
    """returns a tuple of two version numbers"""
    v = ver.strip(string.ascii_lowercase+string.ascii_uppercase)
    return v.split('.')


  def _process_moving_hydrogen( self, chunk):
    chks = chunk.split( ',')
    if len( chks[0]) > 1:
      hs = int( chks[0][1:])  # number of atoms
    else:
      hs = 1
    for i in range( hs):
      ai = int( chks[ i+1]) # atom index
      h = self.structure.add_vertex( atom( symbol='H'))
      self.structure.add_edge( ai, h, e=bond())
      
      

def generate_inchi( m):
  program = "/home/beda/inchi/cINChI11b"

  #try:
  name = os.path.join( tempfile.gettempdir(), "gen_inchi.mol")
  file = open( name, 'w')
  molfile.mol_to_file( m, file)
  file.close()
#  except:
#    print "It was imposible to write a temporary Molfile %s" % name
#    return

  in_name = os.path.join( tempfile.gettempdir(), "gen_inchi.temp")

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
    [line for line in in_file.readlines()]
    out = ( line[6:].strip())
    in_file.close()
  else:
    out = ''

  return out

    




##################################################
# MODULE INTERFACE

reads_text = 1
reads_files = 1
writes_text = 1
writes_files = 1

def text_to_mol( text, include_hydrogens=False, mark_aromatic_bonds=False, calc_coords=1):
  inc = inchi()
  inc.read_inchi( text)
  mol = inc.structure
  #mol.add_missing_bond_orders()
  #if not include_hydrogens:
  #  mol.remove_all_hydrogens()
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
      print "  smiles: ", smiles.mol_to_text( mol)
    t1 = time.time() - t1
    print 'time per cycle', round( 1000*t1/cycles, 2), 'ms'

  repeat = 3
  #inch = "1.12Beta/C18H36N2O6/c1-7-21-13-14-24-10-4-20-5-11-25-17-15-22-8-2-19(1)3-9-23-16-18-26-12-6-20/h1-18H2"
  inch = "1.12Beta/C8H8/c1-2-5-3(1)7-4(1)6(2)8(5)7/h1-8H"
  #inch = '1.12Beta/C3H4O2/c1-2-4-6-5-3-1/h1-3H' #1.12Beta/B2H6/c1-3-2-4-1/h1-2H2'
  #inch = "1.0Beta/C16H22NP/1(2-4-8-15-9-6-7-10-15)3-5-11-16-12-17-14-18-13-16/1-5H2,6-7H,8H2,9-10H,11H2,12-15H"
  print "oasa::INCHI DEMO"
  print "converting following inchi into smiles (%d times)" % repeat
  print "  inchi: %s" % inch
  
  main( inch, repeat)

# DEMO END
##################################################



##################################################
## TODO

## rename layer to sublayer


##################################################
