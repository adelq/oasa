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
from atom import atom
from bond import bond
import periodic_table as PT

import re
from sets import Set
import operator
import time

class smiles( plugin):

  name = "smiles"
  read = 1
  write = 1

  smiles_to_oasa_bond_recode = {'-': 1, '=': 2, '#': 3, ':': 4}
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
    mol = molecule()
    is_text = re.compile("^[A-Z][a-z]?$")
    is_small_text = re.compile( '^[a-z]')
    is_numer = re.compile("[0-9]")
    chunks = re.split( "([A-Z][a-z]?|[^A-Z]|[a-z])", text)
    chunks = filter( None, chunks)
    self._check_the_chunks( chunks)
    last_atom = None
    last_bond = None
    numbers = {}
    bracket_openings = []
    for c in chunks:
      # atom
      if is_text.match( c) or is_small_text.match( c):
        if is_small_text.match( c):
          symbol = c.upper()
        else:
          symbol = c
        a = atom( symbol=symbol)
        mol.add_vertex( a)
        if last_bond and not (not is_small_text.match( c) and last_bond.aromatic):
          mol.add_edge( last_atom, a, e=last_bond)
          last_bond = None
        elif last_atom:
          mol.add_edge( last_atom, a, e=bond())
        last_atom = a
        if is_small_text.match( c):
          # aromatic bond
          last_bond = bond( order=4, type='n')
          last_atom.properties_['aromatic'] = 1  #we need this bellow
        else:
          last_bond = None
      # bond
      elif c in '-=#:.':
        order = self.smiles_to_oasa_bond_recode[ c]
        type = 'n'
        last_bond = bond( order=order, type=type)
      # ring closure
      elif is_numer.match( c):
        b = last_bond or bond()
        if c in numbers:
          mol.add_edge( last_atom, numbers[c], e=b)
          if 'aromatic' in last_atom.properties_:
            last_bond = bond( order=4, type='n')
          else:
            last_bond = None
        else:
          numbers[c] = last_atom
      elif c == '(':
        bracket_openings.append( last_atom)
      elif c == ')':
        last_atom = bracket_openings.pop(-1)

    ## FINISH
    for a in mol.vertices:
      try:
        del a.properties_['aromatic']
      except:
        pass

    self.structure = mol 

  def _check_the_chunks( self, chunks):
    is_text = re.compile("^[A-Z][a-z]?$")
    i = 0
    while i < len( chunks):
      c = chunks[i]
      if is_text.match( c):
        if not c in PT.periodic_table:
          a,b = c
          del chunks[i]
          chunks.insert( i, b)
          chunks.insert( i, a)
          i += 1
      i += 1



  def get_smiles( self, molec):
    mol = molec.copy()
    self.ring_joins = []
    self.branches = {}
    # at first we mark all the atoms with aromatic bonds
    # it is much simple to do it now when all the edges are present
    # we can make use of the properties attribute of the vertex
    for e in mol.edges:
      if e.aromatic:
        for v in e.vertices:
          v.properties_[ 'aromatic'] = 1
    return ''.join( [i for i in self._get_smiles( mol)])



  def _get_smiles( self, mol, start_from=None):
    # single atoms
    if len( mol.vertices) == 1:
      v = mol.vertices[0]
      if 'aromatic' in v.properties_.keys():
        yield v.symbol.lower()
      else:
        yield v.symbol
      for e in self.ring_joins:
        if v in e.get_vertices():
          yield self.recode_oasa_to_smiles_bond( e)
          yield str( self.ring_joins.index( e) +1)
      return
    while not (is_line( mol) and (not start_from or start_from.get_degree() <= 1)):
      if is_pure_ring( mol):
        if start_from:
          self.ring_joins.append( mol.disconnect( start_from, start_from.get_neighbors()[0]))
        else:
          self.ring_joins.append( mol.disconnect( mol.vertices[0], mol.vertices[1])) # edge is returned from disconnect
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
      if 'aromatic' in v.properties_.keys():
        yield v.symbol.lower()
      else:
        yield v.symbol
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
    if 'aromatic' in v.properties_.keys():
      yield v.symbol.lower()
    else:
      yield v.symbol
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




  def disconnect_something( self, mol, start_from=None):
    """returns (broken edge, resulting mol, atom where mol was disconnected, disconnected branch)"""
    # we cannot do much about this part
    if start_from and start_from.get_degree() != 1:
      for e,n in start_from.get_neighbor_edge_pairs():
        if n.get_degree() > 2:
          mol.disconnect( start_from, n)
          return e, mol, None, None
      mol.disconnect( start_from, n)
      return e, mol, None, None
    # at first try to find edges for which degree of neighbors is bigger
    # than [2,2] and at best they are not bridges
    # when no non-bridges are present use the other ones
    #
    # the edges with crowded atoms
    for e in mol.edges:
      d1, d2 = [x.get_degree() for x in e.get_vertices()]
      if d1 > 2 and d2 > 2 and not mol.is_edge_a_bridge_fast_and_dangerous( e):
        mol.disconnect_edge( e)
        return e, mol, None, None
    # the other valuable non-bridge edges
    for e in mol.edges:
      d1, d2 = [x.get_degree() for x in e.get_vertices()]
      if (d1 > 2 or d2 > 2) and not mol.is_edge_a_bridge_fast_and_dangerous( e):
        mol.disconnect_edge( e)
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
        ring_joints = sum([1 for i in the_branch if i in ring_join_vertices])
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
      v1.remove_neighbor( v2)
      v2.remove_neighbor( v1)
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
## MODULE INTERFACE

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
    coords_generator.calculate_coords( mol)
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

  def main( text, cycles):
    t = time.time()
    #mol = molecule()
    #mol._read_file()
    for j in range( cycles):
      mol = text_to_mol( text)
      mol.remove_all_hydrogens()
      text = mol_to_text( mol)
      print "  generated: %s" % text
    t = time.time()-t
    print 'time per cycle', round( 1000*t/cycles, 2), 'ms'

  repeat = 3
  text = "COc5ccc4c2sc(cc2nc4c5)-c(cc1nc3c6)sc1c3ccc6OC"  #"ccc4ccc2cc1cc3ccccc3cc1cc2c4"
  #text = "c1ccc2c1ccccc2"
  #text = "C=1ccC=2C=1C=CC=CC=2"

  print "oasa::SMILES DEMO"
  print "converting following smiles to smiles (%d times)" % repeat
  print "  starting with: %s" % text
  main( text, repeat)
  #import profile
  #profile.run( 'main( text, repeat)')

  # test of equal function for comparison of molecules
##   a = text_to_mol('CC1=CCCCC1')
##   a.add_missing_hydrogens()
##   b = text_to_mol('CC1CC=CCC1')
##   b.add_missing_hydrogens()
##   for l in range( 1, 4):
##     print equals( a, b, level=l)

# DEMO END
##################################################






##################################################
# TODO

# last branch does not need to be branch
# optimize for either speed or human-readability
# at first get rid of the edges common to two or more rings, not critical
# handling of the start_from in disconnect_something as in disconnect_something_simple
# could in disconnect_something happen that a start_from will be returned in a branch?

# in the reading code the more possibilities will be there the more cleanup (last_bond=None,
# last_shit=None ...) will be need. Any way around this?
# if I do "print m.get_smiles( m.structure)" the ring counting does not work after that
## the transformation can be destructive!!!  - check

## THIS IS A PROBLEM : C=1ccC=2C=1C=CC=CC=2  (should be azulene)
