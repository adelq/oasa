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


import graph
from atom import atom
from bond import bond
from sets import Set, ImmutableSet
import copy
import common
import operator




class molecule( graph.graph):

  vertex_class = atom
  edge_class = bond


  def __init__( self, vertices =[]):
    graph.graph.__init__( self, vertices=vertices)
    # aliases
    self.atoms = self.vertices
    self.bonds = self.edges



  def __str__( self):
    return "molecule, %d atoms, %d bonds" % (len( self.vertices), len( self.edges))



  def add_missing_hydrogens( self):
    for v in copy.copy( self.vertices):
      for i in range( v.get_free_valency()):
        h = self.create_vertex()
        h.symbol = 'H'
        self.add_vertex( h)
        self.add_edge( h, v)



  def add_missing_bond_orders( self):
    #for b in self.edges:
    #  b.properties_['free_order'] = min( [a.get_free_valency() for a in b.vertices])
    processed = [1]
    while processed:
      processed = []
      for b in [bo for bo in self.edges if min( [a.get_free_valency() for a in bo.vertices])]:
        a1, a2 = b.get_vertices()
        as1 = [a for a in a1.get_neighbors() if a.get_free_valency() > 0]
        as2 = [a for a in a2.get_neighbors() if a.get_free_valency() > 0]
        if len( as1) == 1 or len( as2) == 1:
          b.order += min( [a.get_free_valency() for a in b.vertices])
          processed.append( b)
      if not processed:
        for b in self.edges:
          i = min( [a.get_free_valency() for a in b.vertices])
          if i:
            processed = [b]
            b.order += i
            break



  def mark_aromatic_bonds( self):
    rings = self.get_all_cycles()   #self.get_smallest_independent_cycles()
    solved = [1]
    # we need to repeat it to mark such things as 'badly' drawn naphtalene (no double bond in the centre) 
    while solved:
      solved = []
      for aring in rings:
        els = [self._get_atoms_possible_aromatic_electrons( a, aring) for a in aring]
        if () in els:
          continue
        for comb in common.gen_combinations_of_series( els):
          if sum( comb) % 4 == 2:
            bring = self.vertex_subgraph_to_edge_subgraph( aring)
            for b in bring:
              b.aromatic = 1
              #b.order = 4
            solved.append( aring)
            break
      for r in solved:
        rings.remove( r)




  def _get_atoms_possible_aromatic_electrons( self, at, ring):
    out = Set()
    if at.charge > 0:
      out.add( 0)
    elif at.charge < 0:
      out.add( 2)
    if at.symbol in ('N','S','O','Se','P'):
      out.add( 2)
    if at.symbol in ('B','Al'):
      out.add( 0)
    for b, a in at.get_neighbor_edge_pairs():
      if b.order > 1 and (a in ring or b.aromatic):
        out.add( 1)
      elif b.order > 1:
        out.add( 0)
    return tuple( out)




  def localize_aromatic_bonds( self):
    """localizes aromatic bonds (does not relocalize already localized ones),
    for those that are not aromatic but marked so
    (it is for instance possible to misuse 'cccc' in smiles to create butadiene)
    they will be properly localized but marked as non-aromatic"""
    rings = self.get_smallest_independent_cycles()
    # sort rings
    rings.sort( lambda x,y: len(y)%2 - len(x)%2) # odd size rings first
    last_rings = []
    while rings:
      # we have to continue with the neighbor rings of the last_ring (or the rings before)
      intersection = None
      if last_rings:
        aring = None
        while last_rings:
          last_ring = last_rings.pop( -1)
          for ring in rings:
            intersection = ring & last_ring
            if intersection:
              aring = ring
              last_rings.append( last_ring)
              break
          if aring:
            break
        if not aring:
          aring = rings.pop(0)
        else:
          rings.remove( aring)
      else:
        aring = rings.pop(0)
      last_rings.append( aring)
      # convert ring from set to list
      aring = self.sort_vertices_in_path( aring, start_from=intersection and intersection.pop() or None)
      # taken from mark_aromatic_bonds
      els = [self._get_atoms_possible_aromatic_electrons( a, aring) for a in aring]
      if () in els:
        continue  # misuse of aromatic bonds (e.g. by smiles) or e.g. tetrahydronaphtalene
      for comb in common.gen_combinations_of_series( els):
        if sum( comb) % 4 == 2:
          aring.append( aring[0])
          comb.append( comb[0])
          # at first we process the bonds that are surely single (like C-S-C in thiophene)
          already_set = None
          for i in range( len( aring) -1):
            if comb[i] + comb[i+1] != 2:
              # these bonds must be single
              b = aring[i].get_edge_leading_to( aring[i+1])
              b.order = 1
              already_set = i
          if already_set != None:
            # we reorder the comb and aring to start from the already set bonds
            j = already_set + 1
            aring = aring[j:len(aring)] + aring[1:j]
            aring.insert( 0, aring[-1])
            comb = comb[j:len(comb)] + comb[1:j]
            comb.insert( 0, comb[-1])
          i = 0
          while i+1 < len( aring):
            if comb[i] + comb[i+1] == 2:
              b = aring[i].get_edge_leading_to( aring[i+1])
              assert b != None # should be
              # to assure alternating bonds
              bs1 = [bo for bo in aring[i].get_neighbor_edges() if bo.order == 2 and bo.aromatic and bo!=b]
              bs2 = [bo for bo in aring[i+1].get_neighbor_edges() if bo.order == 2 and bo.aromatic and bo!=b]
              if len( bs1) == 0 and len( bs2) == 0:
                b.order = 2
              else:
                b.order = 1
            i += 1
          break
        
    self.localize_fake_aromatic_bonds()




  def localize_fake_aromatic_bonds( self):
    to_go = [b for b in self.bonds if b.order == 4]
    while to_go:
      # find the right bond
      b = None
      for bo in to_go:
        bs1, bs2 = bo.get_neighbor_edges2()
        if not bs1 or len( [e for e in bs1 if e.order == 1]) > 0 and len( [e for e in bs1 if e.order == 2]) == 0 \
           or not bs2 or len( [e for e in bs2 if e.order == 1]) > 0 and len( [e for e in bs2 if e.order == 2]) == 0:
          b = bo
          break
        # new start for iteration
      if not b:
        for bo in self.edges:
          if not [e for e in bo.get_neighbor_edges() if e.order != 4]:
            b = bo
            break
      if not b:
        b = to_go.pop(0)
      # the code itself
      b.order = 2
      b.aromatic = 0
      for bo in b.get_neighbor_edges():
        if bo.order == 4:
          bo.order = 1
          bo.aromatic = 0
      # next turn
      to_go = [b for b in self.bonds if b.order == 4]





  def remove_all_hydrogens( self):
    """removes all H atoms"""
    for v in copy.copy( self.vertices):
      if v.symbol == 'H' and v.degree <= 1:
        self.remove_vertex( v)



  def _get_atom_distance_matrix( self, a):
    self.mark_vertices_with_distance_from( a)
    big_out = []
    i = 0
    while 1:
      out = [a.symbol_number for a in self.vertices if a.properties_['d'] == i]
      if out:
        out.sort()
        # out.reverse()
        big_out.append( tuple( out))
      else:
        return big_out
      i += 1




  def get_symmetry_unique_atoms( self):
    out = []
    vs = copy.copy( self.vertices)
    for v in vs:
      out.append( [i for i in self._get_atom_distance_matrix( v)])
    while out:
      k = vs.pop( 0)
      v = out.pop( 0)
      o = [k]
      proc_vs = []
      proc_out = []
      while out:
        k2 = vs.pop( 0)
        v2 = out.pop( 0)
        if v2 == v:
          o.append( k2)
        else:
          proc_vs.append( k2)
          proc_out.append( v2)
      yield o
      out = proc_out
      vs = proc_vs




  def number_atoms_uniquely( self):
    out = {}
    for v in self.vertices:
      out[v] = [i for i in self._get_atom_distance_matrix( v)]
      v.properties_['distance_matrix'] = out[v]
    ms = out.values()
    #self.remove_all_hydrogens()
    ms.sort()
    i = 1
    ret = []
    for m in ms:
      ret.append( out.keys()[ out.values().index( m)])
      #if ret.symbol != 'H':
      #  print "%3d: " % i, ret, ret.get_free_valency()
      i += 1
    return ret




  def _read_file( self, name="/home/beda/oasa/oasa/mol.graph"):
    self.vertices = []
    self.edges = Set()
    f = file( name, 'r')
    vs = f.readline()
    for i in vs.split(' '):
      if i != '\n':
        v = self.create_vertex()
        v.symbol = i
        self.add_vertex( v)
    for l in f.readlines():
      o, a, b = l.split(' ')
      e = self.create_edge()
      e.order = int( o)
      self.add_edge( self.vertices[int(a)], self.vertices[int(b)], e=e)
    f.close()




  def find_longest_mostly_carbon_chain( self):
    if len( self.vertices) < 2:
      return copy.copy( self.vertices)
    ends = [v for v in self.vertices if v.get_degree() == 1]
    paths = []
    for e1 in ends:
      for e2 in ends:
        if e1 != e2:
          path = self.find_path_between( e1, e2)
          paths.append( path)  # could be optimized by not storing all the paths (once we know how to weight the chains)
    # for now we don't take care of chain composition
    l = 0
    p = None
    for path in paths:
      if len( path) > l:
        l = len( path)
        p = path
    return p


def the_right_sorting_function( t1, t2):
  for i,l in enumerate( t1):
    k = t2[i]
    for j in range( len( l)):
      if j > len( k)-1:
        return 1
      if l[j] > k[j]:
        return 1
      elif l[j] < k[j]:
        return -1
    if len( l) < len( k):
      return -1
  return 0




def equals( mol1, mol2, level=0):
  """don't forget to put all hydrogens and bond orders to the bonds
     level 1 - number of atoms and bonds,
     level 2 - the number of atoms with same symbols is the same,
     level 3 - the whole connecivity (no stereo)
     """
  # level 1
  if not level or level >= 1:
    if not len( mol1.vertices) == len( mol2.vertices):
      return False
    if not len( mol1.edges) == len( mol2.edges):
      return False
  # level 2
  if not level or level >= 2:
    symbols1 = [v.symbol for v in mol1.vertices]
    symbols1.sort()
    symbols2 = [v.symbol for v in mol2.vertices]
    symbols2.sort()
    if symbols2 != symbols1:
      return False
  # level 3
  if not level or level >= 3:
    vs1 = mol1.number_atoms_uniquely()
    vs2 = mol2.number_atoms_uniquely()
    for i, v1 in enumerate( vs1):
      v2 = vs2[ i]
      if v1.properties_['distance_matrix'] != v2.properties_['distance_matrix']:
        return False
  return True

#import psyco
#psyco.profile()

##################################################
# DEMO

if __name__ == '__main__':

  def main():
    g = molecule()
    g._read_file()

    for b in g.edges:
      print g.is_edge_a_bridge( b),

##     g.add_missing_hydrogens()

##     for v in g.vertices:
##   #    print v.symbol, v.is_chiral()
##       if v.symbol:
##         ## for a in v.get_neighbors_CIP_sorted():
##         ##   print a, [na.symbol for na in a.get_neighbors()]
##         if v.is_chiral():
##           print v.symbol, [na.symbol for na in v.get_neighbors_CIP_sorted()]



  import time
  t = time.time()
  main()
  print round( 1000*(time.time()-t), 2), 'ms'

# DEMO END
##################################################



##################################################
# TODO

# azulene is wrongly localized
# 0,2 must alternate in aromatic compound
# radical has 1 el in _get_atoms_possible_aromatic_electrons
# test the equals function

##################################################


