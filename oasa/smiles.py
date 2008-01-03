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
from molecule import molecule, equals
import periodic_table as PT
import oasa_exceptions

from config import Config

import re
from sets import Set
import operator
import time

class smiles( plugin):

  name = "smiles"
  read = 1
  write = 1

  smiles_to_oasa_bond_recode = {'-': 1, '=': 2, '#': 3, ':': 4, ".": 0}
  oasa_to_smiles_bond_recode = {1: '', 2: '=', 3: '#', 4:''}

  def __init__( self, structure=None):
    self.structure = structure

  def recode_oasa_to_smiles_bond( self, b):
    if b.aromatic:
      return ''
    else:
      if b.order == 1:
        a1, a2 = b.vertices
        if 'aromatic' in a1.properties_ and 'aromatic' in a2.properties_:
          # non-aromatic bond connecting two aromatic rings, we need to return -
          return '-'
      return self.oasa_to_smiles_bond_recode[ b.order]

  def set_structure( self, structure):
    self.structure = structure

  def get_structure( self):
    return self.structure

  def read_smiles( self, text):
    mol = Config.create_molecule()
    text = "".join( text.split())
    is_text = re.compile("^[A-Z][a-z]?$")
    chunks = re.split( "(\[.*?\]|[A-Z][a-z]?|[^A-Z]|[a-z])", text)
    chunks = filter( None, chunks)
    self._check_the_chunks( chunks)
    last_atom = None
    last_bond = None
    numbers = {}
    bracket_openings = []
    for c in chunks:
      # atom
      if is_text.match( c) or c.islower() or c[0] == "[":
        a = mol.create_vertex()
        if c[0] == "[":
          # atom spec in square brackets
          self._parse_atom_spec( c, a)
        else:
          # just atom symbol
          if c.islower():
            symbol = c.upper()
            a.properties_['aromatic'] = 1
          else:
            symbol = c
          a.symbol = symbol

        mol.add_vertex( a)
        if last_bond: # and not (not 'aromatic' in a.properties_ and last_bond.aromatic):
          mol.add_edge( last_atom, a, e=last_bond)
          last_bond = None
        elif last_atom:
          b = mol.add_edge( last_atom, a)
          if 'aromatic' in a.properties_:
            # aromatic bond
            b.order = 4
            b.type = 'n'
        last_atom = a
        last_bond = None
      # bond
      elif c in '-=#:.':
        order = self.smiles_to_oasa_bond_recode[ c]
        last_bond = mol.create_edge()
        last_bond.order = order
        last_bond.type = 'n'
      # ring closure
      elif c.isdigit():
        if c in numbers:
          if last_bond:
            b = last_bond
          else:
            b = mol.create_edge()
            if "aromatic" in numbers[c].properties_:
              b.order = 4
          mol.add_edge( last_atom, numbers[c], e=b)
          last_bond = None
          del numbers[ c]
        else:
          numbers[c] = last_atom
      elif c == '(':
        bracket_openings.append( last_atom)
      elif c == ')':
        last_atom = bracket_openings.pop(-1)

    ## FINISH
    for a in mol.vertices:
      if not 'explicit_hydrogens' in a.properties_:
        a.raise_valency_to_senseful_value()
      else:
        a.valency = a.occupied_valency + a.properties_['explicit_hydrogens']
        del a.properties_['explicit_hydrogens']
      try:
        del a.properties_['aromatic']
      except:
        pass
        
    self.structure = mol 


  def _parse_atom_spec( self, c, a):
    """c is the text spec,
    a is an empty prepared vertex (atom) instance"""
    bracketed_atom = re.compile("^\[(\d*)([A-z][a-z]?)(.*?)\]")
    m = bracketed_atom.match( c)
    if m:
      isotope, symbol, rest = m.groups()
    else:
      raise ValueError( "unparsable square bracket content '%s'" % c)
    if symbol.islower():
      symbol = symbol.upper()
      a.properties_['aromatic'] = 1
    a.symbol = symbol
    if isotope:
      a.isotope = int( isotope)
    # hydrogens
    _hydrogens = re.search( "H(\d*)", rest)
    h_count = 0
    if _hydrogens:
      if _hydrogens.group(1):
        h_count = int( _hydrogens.group(1))
      else:
        h_count = 1
    a.properties_['explicit_hydrogens'] = h_count
    # charge
    charge = 0
    # one possible spec of charge
    _charge = re.search( "[-+]{2,10}", rest)
    if _charge:
      charge = len( _charge.group(0))
      if _charge.group(0)[0] == "-":
        charge *= -1
    # second one, only if the first one failed
    else:
      _charge = re.search( "([-+])(\d?)", rest)
      if _charge:
        if _charge.group(2):
          charge = int( _charge.group(2))
        else:
          charge = 1
        if _charge.group(1) == "-":
          charge *= -1
    a.charge = charge
    # stereo
    _stereo = re.search( "@+", rest)
    if _stereo:
      stereo = _stereo.group(0)
      a.properties_['stereo'] = stereo

  def _check_the_chunks( self, chunks):
    is_text = re.compile("^[A-Z][a-z]?$")
    i = 0
    while i < len( chunks):
      c = chunks[i]
      if is_text.match( c):
        if not c in PT.periodic_table or c == "Sc": # Sc is S-c not scandium
          a,b = c
          del chunks[i]
          chunks.insert( i, b)
          chunks.insert( i, a)
          i += 1
      i += 1

  def get_smiles( self, mol):
    if not mol.is_connected():
      raise oasa_exceptions.oasa_not_implemented_error( "SMILES", "Cannot encode disconnected compounds, such as salts etc.")
    #mol = molec.copy()
    self.ring_joins = []
    self.branches = {}
    # at first we mark all the atoms with aromatic bonds
    # it is much simple to do it now when all the edges are present
    # we can make use of the properties attribute of the vertex
    for e in mol.edges:
      if e.aromatic:
        for v in e.vertices:
          v.properties_[ 'aromatic'] = 1
    ret = ''.join( [i for i in self._get_smiles( mol)])
    mol.reconnect_temporarily_disconnected_edges()
    # this is needed because the way temporarily_disconnected edges are handled is not compatible with the way smiles
    # generation failed - it splits the molecule while reusing the same atoms and bonds and thus disconnected bonds accounting fails
    for e in mol.edges:
      e.disconnected = False

    return ret



  def _get_smiles( self, mol, start_from=None):
    # single atoms
    if len( mol.vertices) == 1:
      v = mol.vertices[0]
      yield self._create_atom_smiles( v)
      for e in self.ring_joins:
        if v in e.get_vertices():
          yield self.recode_oasa_to_smiles_bond( e)
          yield str( self.ring_joins.index( e) +1)
      return
    while not (is_line( mol) and (not start_from or start_from.get_degree() <= 1)):
      if is_pure_ring( mol):
        if start_from:
          self.ring_joins.append( mol.temporarily_disconnect_edge( start_from.neighbor_edges[0]))
        else:
          self.ring_joins.append( mol.temporarily_disconnect_edge( list( mol.edges)[0]))
      else:
        e, mol, branch_vertex, branch = self.disconnect_something( mol, start_from=start_from)
        if branch_vertex:
          if branch_vertex in self.branches:
            self.branches[ branch_vertex].append((e, branch))
          else:
            self.branches[ branch_vertex] = [(e, branch)]
        else:
          self.ring_joins.append( e)
    try:
      start, end = filter( lambda x: x.get_degree() == 1, mol.vertices)
    except:
      #print filter( lambda x: x.get_degree() == 1, mol.vertices)
      raise "shit"
    if start_from == end:
      start, end = end, start
    v = start
    last = None
    while v != end:
      yield self._create_atom_smiles( v)
      # the atom
      for e in self.ring_joins:
        if v in e.get_vertices():
          yield self.recode_oasa_to_smiles_bond( e)
          yield str( self.ring_joins.index( e) +1)
      # branches
      if v in self.branches:
        for edg, branch in self.branches[ v]:
          yield '('
          yield self.recode_oasa_to_smiles_bond( edg)
          v1, v2 = edg.vertices
          vv = (v1 != v) and v1 or v2
          for i in self._get_smiles( branch, start_from=vv):
            yield i
          yield ')'
      # bond leading to the neighbor
      for e, neighbor in v.get_neighbor_edge_pairs():
        if neighbor != last:
          yield self.recode_oasa_to_smiles_bond( e)
          last = v
          v = neighbor
          break
    # the last atom - should make it somehow not to need this piece of code
    yield self._create_atom_smiles( v)
    for e in self.ring_joins:
      if v in e.get_vertices():
        yield self.recode_oasa_to_smiles_bond( e)
        yield str( self.ring_joins.index( e) +1)
      if v in self.branches:
        for edg, branch in self.branches[ v]:
          yield '('
          yield self.recode_oasa_to_smiles_bond( edg)
          v1, v2 = edg.vertices
          vv = (v1 != v) and v1 or v2
          for i in self._get_smiles( branch, start_from=vv):
            yield i
          yield ')'



  def _create_atom_smiles( self, v):
    if 'aromatic' in v.properties_.keys():
      symbol = v.symbol.lower()
    else:
      symbol = v.symbol

    if v.isotope or v.charge != 0 or v.valency != PT.periodic_table[ v.symbol]['valency'][0] or 'stereo' in v.properties_:
      # we must use square bracket
      isotope = v.isotope and str( v.isotope) or ""
      # charge
      if v.charge:
        sym = v.charge < 0 and "-" or "+"
        charge = sym + (abs( v.charge) > 1 and str( abs( v.charge)) or "")
      else:
        charge = ""
      # explicit hydrogens
      num_h = v.valency - v.occupied_valency
      h_spec = (num_h and "H" or "") + (num_h > 1 and str( num_h) or "")
      # stereo
      stereo = v.properties_.get( "stereo", "")
      return "[%s%s%s%s%s]" % (isotope, symbol, stereo, h_spec, charge)
    else:
      # no need to use square brackets
      return symbol


  def disconnect_something( self, mol, start_from=None):
    """returns (broken edge, resulting mol, atom where mol was disconnected, disconnected branch)"""
    # we cannot do much about this part
    if start_from and start_from.get_degree() != 1:
      for e,n in start_from.get_neighbor_edge_pairs():
        if n.get_degree() > 2:
          mol.temporarily_disconnect_edge( e)
          return e, mol, None, None
      mol.temporarily_disconnect_edge( e)
      return e, mol, None, None
    # at first try to find edges for which degree of neighbors is bigger
    # than [2,2] and at best they are not bridges
    # when no non-bridges are present use the other ones
    #
    # the edges with crowded atoms
    for e in mol.edges:
      d1, d2 = [x.get_degree() for x in e.get_vertices()]
      if d1 > 2 and d2 > 2 and not mol.is_edge_a_bridge_fast_and_dangerous( e):
        mol.temporarily_disconnect_edge( e)
        return e, mol, None, None
    # the other valuable non-bridge edges
    for e in mol.edges:
      d1, d2 = [x.get_degree() for x in e.get_vertices()]
      if (d1 > 2 or d2 > 2) and not mol.is_edge_a_bridge_fast_and_dangerous( e):
        mol.temporarily_disconnect_edge( e)
        return e, mol, None, None
    # there are no non-bridges
    # we want to split off the smallest possible chunks
    min_size = None
    the_right_edge = None
    the_right_mol = None
    the_right_branch = None
    the_right_branch_atom = None
    ring_joints_in_branch = 1000
    ring_join_vertices = Set( reduce( operator.add, [e.vertices for e in self.ring_joins], []))
    for e in mol.edges:
      d1, d2 = [x.get_degree() for x in e.get_vertices()]
      if d1 > 2 or d2 > 2:
        ps = mol.get_pieces_after_edge_removal( e)
        if len( ps) == 1:
          print "impossible"
          continue
        lenghts = map( len, ps)
        ms = min( lenghts)
        p1, p2 = ps
        the_mol = (len( p1) < len( p2)) and p2 or p1
        the_branch = (p1 == the_mol) and p2 or p1
        ring_joints = len( [i for i in the_branch if i in ring_join_vertices])
        if not min_size or ms < min_size or ring_joints_in_branch > ring_joints:
          min_size = ms
          the_right_edge = e
          the_right_mol = the_mol
          the_right_branch = the_branch
          ring_joints_in_branch = ring_joints
    if the_right_edge:
      # what is possible to make here instead in the loop is made here
      # it saves time
      v1, v2 = the_right_edge.vertices
      mol.temporarily_disconnect_edge( the_right_edge)
      the_right_branch_atom = (v1 in the_right_mol) and v1 or v2
      the_right_mol = mol.get_induced_subgraph_from_vertices( the_right_mol)
      the_right_branch = mol.get_induced_subgraph_from_vertices( the_right_branch)
      return (the_right_edge,
              the_right_mol,
              the_right_branch_atom,
              the_right_branch)
    #print mol, mol.is_connected()
    raise "fuck, how comes!?"




  def disconnect_something_simple( self, mol, start_from=None):
    """returns (broken edge, resulting mol, atom where mol was disconnected, disconnected branch)"""
    # we cannot do much about this part
    if not mol.is_connected():
      print "unconnected ", mol
    if start_from and start_from.get_degree() > 1:
      e = start_from._neighbors.keys()[0]
      mol.disconnect_edge( e)
      ps = [i for i in mol.get_connected_components()]
      if len( ps) == 1:
        return e, mol, None, None
      else:
        p1 = mol.get_induced_subgraph_from_vertices( ps[0])
        p2 = mol.get_induced_subgraph_from_vertices( ps[1])
        v1, v2 = e.get_vertices()
        p = (start_from in p1.vertices) and p1 or p2
        px = (p1 == p) and p2 or p1
        v = (v1 in p.vertices) and v1 or v2
        return e, p, v, px

    for e in mol.edges:
      d1, d2 = [x.get_degree() for x in e.get_vertices()]
      if (d1 > 2 or d2 > 2):
        mol.disconnect_edge( e)
        ps = [i for i in mol.get_connected_components()]
        if len( ps) == 1:
          return e, mol, None, None
        else:
          p1 = mol.get_induced_subgraph_from_vertices( ps[0])
          p2 = mol.get_induced_subgraph_from_vertices( ps[1])
          v1, v2 = e.get_vertices()
          if start_from:
            p = (start_from in p1.vertices) and p1 or p2
            px = (p1 == p) and p2 or p1
            v = (v1 in p.vertices) and v1 or v2
            return e, p, v, px

          v = (v1 in p1.vertices) and v1 or v2
          return e, p1, v, p2
    print mol, mol.is_connected(), ',', map( len, mol.get_connected_components()), ',', start_from
    raise "fuck, how comes!?"



def is_line( mol):
  """all degrees are 2 except of two with degree 1"""
  if len( mol.vertices) == 1:
    return True
  ones = 0
  for v in mol.vertices:
    d = v.get_degree() 
    if d == 1:
      if ones == 2:
        return False
      ones += 1
    elif d != 2:
      return False
  if ones == 2:
    return True
  return False

def is_pure_ring( mol):
  return filter( lambda x: x.get_degree() != 2, mol.vertices) == []


##################################################
## MODULE INTERFACE - newstyle

from converter_base import converter_base

class smiles_converter( converter_base):

  # standard converter attrs
  reads_text = True
  writes_text = True
  reads_files = True
  writes_files = True

  default_configuration = {"R_GENERATE_COORDS": True,
                           "R_BOND_LENGTH": 1,
                           "R_LOCALIZE_AROMATIC_BONDS": True,
                           
                           "W_AROMATIC_BOND_AUTODETECT": True,
                           "W_INDIVIDUAL_MOLECULE_SEPARATOR": ".",
                           }

  def __init__( self):
    converter_base.__init__( self)

  def mols_to_text( self, structures):
    converter_base.mols_to_text( self, structures)
    sm = smiles()
    ret = []
    for mol in structures:
      if self.configuration["W_AROMATIC_BOND_AUTODETECT"]:
        mol.mark_aromatic_bonds()
      ret.append( sm.get_smiles( mol))
    return self.configuration["W_INDIVIDUAL_MOLECULE_SEPARATOR"].join( ret)

  def text_to_mols( self, text):
    converter_base.text_to_mols( self, text)
    sm = smiles()
    sm.read_smiles( text)
    mol = sm.structure
    zero_bonds = [e for e in mol.edges if e.order == 0]
    for b in zero_bonds:
      mol.disconnect_edge( b)
    mols = mol.get_disconnected_subgraphs()
    for mol in mols:
      if self.configuration["R_LOCALIZE_AROMATIC_BONDS"]:
        mol.localize_aromatic_bonds()
        for b in mol.bonds:
          b.aromatic = 0
      if self.configuration["R_GENERATE_COORDS"]:
        coords_generator.calculate_coords( mol, bond_length=self.configuration['R_BOND_LENGTH'])
    return mols

  def mols_to_file( self, structures, f):
    converter_base.mols_to_file( self, structures, f)
    f.write( self.mols_to_text( structures))

  def file_to_mols( self, f):
    converter_base.file_to_mols( self, f)
    mols = []
    for line in f:
      mol = self.text_to_mols( line)
      mols.extend( mol)
    return mols
    
converter = smiles_converter

# END OF MODULE INTERFACE
##################################################

##################################################
## MODULE INTERFACE - oldstyle

import coords_generator

reads_text = True
writes_text = True
reads_files = True
writes_files = True

def mol_to_text( structure):
  sm = smiles()
  structure.mark_aromatic_bonds()
  return sm.get_smiles( structure)

def text_to_mol( text, calc_coords=1):
  sm = smiles()
  sm.read_smiles( text)
  mol = sm.structure
  mol.localize_aromatic_bonds()
  for b in mol.bonds:
    b.aromatic = 0
  if calc_coords:
    coords_generator.calculate_coords( mol, bond_length=calc_coords)
  return mol

def mol_to_file( mol, f):
  f.write( mol_to_text( mol))
  
def file_to_mol( f):
  return text_to_mol( f.read())


# END OF MODULE INTERFACE
##################################################



##################################################
# DEMO

if __name__ == '__main__':

  import sys

  def main( text, cycles):
    t = time.time()
    conv = converter()
    for j in range( cycles):
      mols = conv.text_to_mols( text)
      for mol in mols:
        mol.remove_all_hydrogens()
        print "  summary formula:   ", mol.get_formula_dict()        
      text = conv.mols_to_text( mols)
      print "  generated SMILES:   %s" % text
      print "  --"
    t = time.time()-t
    print 'time per cycle', round( 1000*t/cycles, 2), 'ms'

  repeat = 3

  if not len( sys.argv) > 1:
    text = "COc5ccc4c2sc(cc2nc4c5)-c(cc1nc3c6)sc1c3ccc6OC"  #"ccc4ccc2cc1cc3ccccc3cc1cc2c4"
  else:
    text = sys.argv[1]

  print "oasa::SMILES DEMO"
  print "converting following smiles to smiles (%d times)" % repeat
  print "  starting with:      %s" % text
  print "  --------------------"
  main( text, repeat)

# DEMO END
##################################################






##################################################
# TODO

# last branch does not need to be branch
# optimize for either speed or human-readability
# at first get rid of the edges common to two or more rings, not critical
# handling of the start_from in disconnect_something as in disconnect_something_simple
# could in disconnect_something happen that a start_from will be returned in a branch?
# if I do "print m.get_smiles( m.structure)" the ring counting does not work after that
## the transformation can be destructive!!!  - check

## THIS IS A PROBLEM : C=1ccC=2C=1C=CC=CC=2  (should be azulene)
