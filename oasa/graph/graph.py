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


"""this module contains a graph class that provides a minimalistic
graph implementation suitable for analysis of chemical problems"""

from sets import Set
from sets import ImmutableSet as ImSet
from edge import edge
from vertex import vertex
import warnings
import copy
from types import *
import operator
import time

class graph:
  """provides a minimalistic graph implementation suitable for analysis of chemical problems,
  even if some care was taken to make the graph work with nonsimple graphs, there are cases where it won't!"""

  def __init__( self, vertices = []):
    if vertices:
      self.vertices = vertices
    else:
      self.vertices = []
    self.edges = Set()
    self.connect = []
    for i in range( len( self.vertices)):
      self.connect.append( len( self.vertices) *[0])

  def __str__( self):
    str = "graph G(V,E), |V|=%d, |E|=%d" % ( len( self.vertices), len( self.edges))
    return str

  def copy( self):
    """provides a really shallow copy, the vertex and edge objects will remain the same,
    only the graph itself is different"""
    clazz = self.__class__
    c = clazz( vertices= copy.copy( self.vertices))
    for i, line in enumerate( self.connect):
      for j, e in enumerate( line):
        if e:
          c.add_edge( i, j, e)
    return c

  def deep_copy( self):
    """provides a deep copy of the graph. The result is an isomorphic graph,
    all the used objects are different"""
    c = graph()
    for v in self.vertices:
      c.add_vertex()
    for i, line in enumerate( self.connect):
      for j, e in enumerate( line):
        if j > i and e:
          c.add_edge( i, j)
    return c
    

  ## MODIFICATION METHODS

  def add_vertex( self, v=None):
    """adds a vertex to a graph, if v argument is not given creates a new one.
    returns None if vertex is already present or the vertex instance if successful"""
    if not v:
      v = vertex()
    if v not in self.vertices:
      self.vertices.append( v)
    else:
      warnings.warn( "Added vertex is already present in graph %s" % str( v), UserWarning, 2)
      return None
    for line in self.connect:
      line.append( 0)
    self.connect.append( len( self.vertices) *[0])
    return v



  def add_edge( self, v1, v2, e=None):
    """adds an edge to a graph connecting vertices v1 and v2, if e argument is not given creates a new one.
    returns None if operation fails or the edge instance if successful"""
    i1 = self._get_vertex_index( v1)
    i2 = self._get_vertex_index( v2)
    if i1 == None or i2 == None:
      warnings.warn( "Adding edge to a vertex not present in graph failed (of course)", UserWarning, 3)
      return None
    # to get the vertices if v1 and v2 were indexes
    v1 = self.vertices[ i1]
    v2 = self.vertices[ i2]
    if not e:
      e = edge()
    e.set_vertices( (v1,v2))
    self.edges.add( e)
    self.connect[i1][i2] = e
    self.connect[i2][i1] = e
    v1.add_neighbor( v2, e)
    v2.add_neighbor( v1, e)
    return e



  def disconnect( self, v1, v2):
    """disconnects vertices v1 and v2 (can be also indexes), on success returns the edge"""
    i1 = self._get_vertex_index( v1)
    i2 = self._get_vertex_index( v2)
    v1 = self.vertices[ i1]
    v2 = self.vertices[ i2]
    if i1 != None and i2 != None:
      e = self.connect[i1][i2]
      self.edges.remove( e)
      self.connect[i1][i2] = 0
      self.connect[i2][i1] = 0
      v1.remove_neighbor( v2)
      v2.remove_neighbor( v1)
      return e
    else:
      return None



  def disconnect_edge( self, e):
    v1, v2 = e.vertices
    self.disconnect( v1, v2)



  def remove_vertex( self, v):
    for neigh in v.get_neighbors():
      self.disconnect( v, neigh)
    i = self.vertices.index( v)
    del self.vertices[i]
    del self.connect[i]
    for line in self.connect:
      del line[i]
      


  def get_edge_between( self, v1, v2):
    """takes two vertices either as vertex instances or indexes in self.vertices"""
    i1 = self._get_vertex_index( v1)
    i2 = self._get_vertex_index( v2)
    return self.connect[i1][i2]



  def get_vertices_of_edge( self, e):
    """Info - available also trough the edge.get_vertices()"""
    for i, line in enumerate( self.connect):
      for j, ed in enumerate( line):
        if e == ed:
          return self.vertices[i], self.vertices[j]


  ## PROPERTIES METHODS
  ## BOOLEAN
    
  def is_connected( self):
    if len( self.edges) < len( self.vertices) -1:
      # in this case it cannot be connected
      return False
    i = 0
    for x in self.get_connected_components():
      i += 1
      if i > 1:
        return False
    return True
  #return len( self.get_connected_components()) == 1

  def is_tree( self):
    return self.is_connected() and len( self.vertices)-1 == len( self.edges)

  def is_cycle( self):
    for i in self.get_degrees():
      if i != 2:
        return False
    return True

  def is_euler( self):
    for i in self.get_degrees():
      if i % 2:
        return False
    return True

  def contains_cycle( self):
    """this assumes that the graph is connected"""
    if not self.is_connected():
      return None
    return not self.is_tree()

  def is_edge_a_bridge( self, e):
    v1, v2 = e.vertices
    self.disconnect( v1, v2)
    #x = len( [i for i in self.get_connected_components()])
    x = self.get_connected_components().next()
    if len( x) < len( self.vertices):
      x = 1
    else:
      x = 0
    self.add_edge( v1, v2, e=e)
    return x

  def is_edge_a_bridge_fast_and_dangerous( self, e):
    """should be used only in case of repetitive questions for the same edge in cases
    where no edges are added to the graph between the questions (if brigde==1 the value
    is stored and returned, which is safe only in case no edges are added)"""
    try:
      return e.properties_['bridge']
    except:
      if self.is_edge_a_bridge( e):
        e.properties_['bridge'] = 1
        return 1
      else:
        return 0

  def get_pieces_after_edge_removal( self, e):
    v1, v2 = e.vertices
    self.disconnect( v1, v2)
    ps = [i for i in self.get_connected_components()]
    self.add_edge( v1, v2, e=e)
    return ps

  def get_size_of_pieces_after_edge_removal( self, e):
    v1, v2 = e.vertices
    self.disconnect( v1, v2)
    ps = [i for i in self.get_connected_components()]
    self.add_edge( v1, v2, e=e)
    return map( len, ps)
    

##   def contains_vertex_path( self):
##     """can contain only two vertices with order different from 2"""
##     non2 = 0
##     for v in self.vertices:
##       if v.get_degree() != 2:
##         if non2 == 2:
##           return False
##         non2 += 1
##     return True

  ## ANALYSIS

  def get_connected_components( self):
    """returns the connected components of graph in a form o list of lists of vertices"""
    comp = Set() # just processed component 
    comps = []
    not_processed = Set( self.vertices)
    if not_processed:
      recent = Set( [not_processed.pop()])
      comp |= recent
    while not_processed:
      recent = Set( reduce( operator.add, [a.get_neighbors() for a in recent], [])) & not_processed
      if not recent:
        yield comp
        recent = Set( [not_processed.pop()])
        comp = recent
      else:
        comp |= recent
        not_processed -= recent


  def get_disconnected_subgraphs( self):
    vss = self.get_connected_components()
    out = []
    for vs in vss:
      out.append( self.get_induced_subgraph_from_vertices( vs))
    return out

  def get_induced_subgraph_from_vertices( self, vs):
    """it creates a new graph, however uses the old vertices and edges!"""
    g = self.__class__()
    for v in vs:
      g.add_vertex( v)
    for e in self.edges:
      v1, v2 = e.get_vertices()
      if v1 in vs and v2 in vs:
        g.add_edge( v1, v2, e)  # BUG - it should copy the edge?
    return g
    
      
  def get_vertex_degree( self, v):
    """accepts the vertex or its index.
    Info - available also trough the vertex.get_degree()"""
    i1 = self._get_vertex_index( v)
    return self._get_degree_of_index( i1)

  def get_degrees( self):
    """returns a generator of degrees, this is useful because for many properties
    the whole list is not important"""
    for i in range( len( self.vertices)):
      yield self._get_degree_of_index( i)
  
  def get_neighbors( self, v):
    """Info - available also trough the vertex.get_neighbors()"""
    for i in self.get_neighbors_indexes( v):
      yield self.vertices[ i]

  def get_neighbors_indexes( self, v):
    i1 = self._get_vertex_index( v)
    for i, x in enumerate( self.vertices[ i1]):
      if x and i!=i1:
        yield i

  def get_smallest_independent_cycles( self):
    """returns a set of smallest possible independent cycles,
    other cycles in graph are guaranteed to be combinations of them"""
    #t1 = time.time()
    ret = list( filter_off_dependent_cycles( self._get_some_cycles()))
    #print " %.5f - smallest rings search" % (time.time() - t1)
    return ret

  def get_smallest_cycles2( self):
    return list( filter_off_dependent_cycles( self.get_all_cycles()))

  def get_all_cycles( self):
    """takes all the smallest independent cycles found and combines them to provide
    the bigger dependent cycles. Warning - is relatively time consuming for fused cycles
    that give rise to many possible combinations"""
    #t1 = time.time()
    vcycles = self.get_smallest_independent_cycles() # vertex defined cycles
    ecycles = map( self.vertex_subgraph_to_edge_subgraph, vcycles) # edge defined cycles
    i = 0
    while i < len( ecycles):
      j = i+1
      while j < len( ecycles):
        if j != i:
          p = (ecycles[i] & ecycles[j]) 
          if p and not (vcycles[j] <= vcycles[i]) and not (vcycles[i] <= vcycles[j]) and self.is_path( p):
            new = ecycles[i] ^ ecycles[j]
            if new and new not in ecycles:
              vnew = self.edge_subgraph_to_vertex_subgraph( new)
              ecycles.append( new)
              vcycles.append( vnew)
        j += 1
      i += 1
    #print " %.2f ms - all rings search" % (1000*(time.time() - t1))
    return vcycles

  def mark_vertices_with_distance_from( self, v):
    self.clean_distance_from_vertices()
    v.properties_['d'] = 0
    self._mark_vertices_with_distance_from( v)

  def clean_distance_from_vertices( self):
    for i in self.vertices:
      try:
        del i.properties_['d']
      except KeyError:
        pass

  def vertex_subgraph_to_edge_subgraph( self, cycle):
    ret = Set()
    for v1 in cycle:
      for v2 in cycle:
        if v2 in v1.get_neighbors() and abs( v2.properties_['d']-v1.properties_['d']) <= 1:
          ret.add( self.get_edge_between( v1, v2))
    return ret

  def edge_subgraph_to_vertex_subgraph( self, cycle):
    ret = Set()
    for e in cycle:
      v1, v2 = e.get_vertices()
      ret.add( v1)
      ret.add( v2)
    return ret

  def get_new_induced_subgraph( self, vertices, edges):
    """returns a induced subgraph that is newly created and can be therefore freely
    changed without worry about the original."""
    sub = graph()
    r1, r2 = [], []
    for v in vertices:
      r1.append( v)
      sub.add_vertex()
    for e in edges:
      v1, v2 = e.get_vertices()
      i1 = r1.index( v1)
      i2 = r1.index( v2)
      sub.add_edge( i1, i2)
    return sub

  def is_path( self, edges):
    sub = self.get_new_induced_subgraph( self.edge_subgraph_to_vertex_subgraph( edges), edges)
    return sub.is_connected()

  def find_path_between( self, start, end, dont_go_through=[]):
    """finds path between two vertices, if dont_go_through is given (as a list of vertices and edges),
    only paths not containing these vertices will be given (or None is returned if such a path
    does not exist"""
    ### DOES NOT WORK WELL WITH RINGS, FOR THIS RECURSIVE DESIGN WILL BE NEEDED
    self.mark_vertices_with_distance_from( start)
    d = end.properties_['d']
    out = [end]
    rend = end
    for i in range( d, 0, -1):
      vs = [v for e,v in rend.get_neighbor_edge_pairs() if v.properties_['d'] == i-1 and (v not in dont_go_through or e not in dont_go_through)]
      if not vs:
        return None
      v = vs[0]
      out.append( v)
      rend = v
    return out

  def sort_vertices_in_path( self, path, start_from=None):
    """returns None if there is no path"""
    rng = copy.copy( path)
    if start_from:
      a = start_from
      rng.remove( a)
    else:
      a = None
      # for acyclic path we need to find one end
      for at in path:
        if len( [v for v in at.get_neighbors() if v in path]) == 1:
          a = at
          rng.remove( at)
          break
      if not a:
        a = rng.pop() # for rings
    out = [a]
    while rng:
      try:
        a = [i for i in a.get_neighbors() if i in rng][0]
      except IndexError:
        return None
      out.append( a)
      rng.remove( a)
    return out


  # PRIVATE METHODS

  def _get_vertex_index( self, v):
    """if v is already an index, return v, otherwise return index of v on None"""
    if type( v) == IntType and v < len( self.vertices):
      return v
    try:
      return self.vertices.index( v)
    except ValueError:
      return None

  def _get_degree_of_index( self, i1):
    """helps when index is already known"""
    deg = 0
    for i, b in enumerate( self.connect[i1]):
      if b:
        if i == i1:
          deg += 2
        else:
          deg += 1
    return deg

  def _read_file( self, name="/home/beda/oasa/oasa/mol.graph"):
    self.vertices = []
    self.edges = Set()
    self.connect = []
    f = file( name, 'r')
    vs = f.readline()
    [self.add_vertex() for i in vs.split(" ")]
    for l in f.readlines():
      c, a, b = l.split(' ')
      self.add_edge( self.vertices[int(a)], self.vertices[int(b)])
    f.close()

  def _mark_vertices_with_distance_from( self, v):
    d = v.properties_['d']
    to_mark = []
    for i in v.get_neighbors():
      if 'd' not in i.properties_ or i.properties_['d'] > d+1:
        i.properties_['d'] = d+1
        to_mark.append( i)
    for i in to_mark:
      self._mark_vertices_with_distance_from( i)

  def _get_some_cycles( self):
    self.mark_vertices_with_distance_from( self.vertices[0])
    for end in self._get_all_ring_end_points():
      for start in self._get_all_ring_start_points():
        ring = is_there_a_ring_between( start, end)
        if ring:
          yield ring

  def _get_all_ring_end_points( self):
    already_there = []
    for v in self.vertices:
      if v in already_there:
        continue
      vs_ed, vs_ver = is_ring_end_vertex( v)
      if vs_ed:
        yield vs_ed
      if vs_ver:
        already_there.extend( vs_ver)
        yield vs_ver

  def _get_all_ring_start_points( self):
    already_there = []
    for v in self.vertices:
      if v in already_there:
        continue
      vs_ed, vs_ver = is_ring_start_vertex( v)
      if vs_ed:
        yield vs_ed
      if vs_ver:
        already_there.extend( vs_ver)
        yield vs_ver



## common functions

def simplify_connectivity_matrix( mat):
  """makes a underlying simple graph connectivity matrix from mat, works on the mat itself"""
  for i, line in enumerate( mat):
    for j, x in enumerate( line):
      if mat[i][j] and i!=j:
        mat[i][j] = 1
      else:
        mat[i][j] = 0

def get_index_of_vertex_connected_to_first_vertex( mat):
  """returns an index of vertex that is connected to the vertex in first line
  of the matrix"""
  for i, x in enumerate( mat[0]):
    if x:
      return i

def fuse_vertices( i1, i2, mat):
  """fuses vertices with indexes i1, i2 to i1, works on the mat"""
  # sum the lines
  for i in range( len( mat)):
    if i != i1:
      mat[i1][i] = mat[i1][i] or mat[i2][i]
  # remove i2
  del mat[i2]
  for line in mat:
    del line[i2]
  # copy row i1 to column i1
  for i, x in enumerate( mat[i1]):
    mat[i][i1] = x

def is_ring_end_vertex( v):
# NEEDS NEW COMMENT
#  """tells if a vertex has two neighbors with 'd' one smaller and equal or
#  one neighbor with equal 'd'. These are the conditions for a cycle end.
#  Returns the Set([v]), in second case 'v' and the other with same 'd'"""
  ed, ver = None, None
  for x in v.get_neighbors():
    if v.properties_['d'] == x.properties_['d']:
      if len( v.get_neighbors_with_distance( v.properties_['d']-1)) and len( x.get_neighbors_with_distance( v.properties_['d']-1)):
        ed = Set([ x, v])
    for y in v.get_neighbors():
      if x != y:
        if (x.properties_['d'] == y.properties_['d']) and (x.properties_['d'] == v.properties_['d']-1):
          ver = Set([v])
  return ed, ver

def is_ring_start_vertex( v):
# NEEDS NEW COMMENT
#  """tells if a vertex has two neighbors with 'd' one higher and equal or one neighbor with
#  equal 'd' (then both have neighbors with one higher 'd'). These are the conditions for a cycle start.
#  Returns boolean"""
  ed, ver = None, None
  for x in v.get_neighbors():
    if v.properties_['d'] == x.properties_['d']:
      if len( v.get_neighbors_with_distance( v.properties_['d']+1)) and len( x.get_neighbors_with_distance( v.properties_['d']+1)):
        ed = Set([ x, v])
    for y in v.get_neighbors():
      if x != y:
        if (x.properties_['d'] == y.properties_['d']) and (x.properties_['d'] == v.properties_['d']+1):
          ver = Set([v])
          break
  return ed, ver

def get_first_closer_by_one( v):
  d = v.properties_['d']
  for x in v.get_neighbors():
    if x.properties_['d'] == d-1:
      return x
  return None

def is_there_a_ring_between( start, end):
  pths = []
  for e in end:
    for s in start:
      ps = get_paths_down_to( e, s)
      if ps:
        ## for end-edge there is only one path from each point important
        if len( end) == 2:
          ps = [ps[0]]
        ## for end-vertex the end vertex must not be part of the set in order to filter it
        if len( end) == 1:
          for p in ps:
            p.remove( e)
        pths.extend( ps)
  ## we use sets for filtering off the path that have common vertices 
  paths = [Set( p) for p in pths]
  final = []
  while len( paths) > 0:
    p1 = paths.pop(0)
    to_del= []
    for p2 in paths:
      if p1 & p2:
        to_del.append( p2)
    for p2 in to_del:
      paths.remove( p2)
    final.append( p1)
  if len( final) >= 2:
    final.append( end | start)
    ret = Set( reduce( operator.or_, final))
    return ret
  return False
                
def get_paths_down_to( end, start):
  paths = []
  if end == start:
    return None
  for x in end.get_neighbors():
    if x.properties_['d'] == end.properties_['d']-1:
      ps = get_path_down_to( x, start)
      if ps != None:
        ps.append( end)
        paths.append( ps)
  return paths

def get_path_down_to( end, start):
  if end == start:
    return []
  for x in end.get_neighbors():
    if x.properties_['d'] == end.properties_['d']-1:
      ps = get_path_down_to( x, start)
      if ps != None:
        ps.append( end)
        return ps
  return None

def filter_off_dependent_cycles( cycles):
  """this filtres off all cycles that are combinations of two smaller cycles.
  I don't want to go deeper (3 cycles) - it would be really slow and ugly"""
  # at first we use filter_off_supercycles - it is faster and catches almost all
  cs = list( filter_off_supercycles( cycles))
  i = 0
  while i < len( cs):
    j = i + 1
    while j < len( cs):
      k = j + 1
      while k < len( cs):
        if len( cs[k]) > len( cs[j]):
          if cs[k] <= (cs[i] | cs[j]):
            del cs[k]
            continue
        k += 1
      j += 1
    i += 1
  return cs
  
def filter_off_supercycles( cycles):
  """filtres off all the bigger cycles which completely contain any of smaller cycles;
  these should be dependent"""
  cs = list( cycles)
  cs.sort( lambda x, y: len(x)-len(y))
  while len( cs) > 0:
    to_del = []
    c0 = cs.pop(0)
    for c1 in cs:
      if c0 <= c1:
        to_del.append( c1)
    for c1 in to_del:
      cs.remove( c1)
    yield c0
  
  

## import profile

## g = graph()
## g._read_file()

## print [c for c in g.get_connected_components()]

## for e in g.edges:
##   print g.is_edge_a_bridge( e),
## print

## print len( [i for i in g.get_connected_components()])
## for e in g.edges:
##   d1, d2 = [x.get_degree() for x in e.get_vertices()]
##   if d1 != 2 or d2 != 2:
##     v1, v2 = e.vertices
##     break
## g.disconnect( v1, v2)
## print len( [i for i in g.get_connected_components()])


## g = g.deep_copy()
## print g
## #profile.run( 'g.get_all_cycles()')
## c = g.get_smallest_independent_cycles()
## print 'cycles %d, lengths %s' % (len( c), str( map( len, c)))
## c = g.get_all_cycles()
## print 'cycles %d, lengths %s' % (len( c), str( map( len, c)))




##################################################
# TODO


##################################################
