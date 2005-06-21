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

  vertex_class = vertex
  edge_class = edge


  def __init__( self, vertices = []):
    if vertices:
      self.vertices = vertices
    else:
      self.vertices = []
    self.edges = Set()

  def __str__( self):
    str = "graph G(V,E), |V|=%d, |E|=%d" % ( len( self.vertices), len( self.edges))
    return str

  def copy( self):
    """provides a really shallow copy, the vertex and edge objects will remain the same,
    only the graph itself is different"""
    c = self.create_graph()
    c.vertices = copy.copy( self.vertices)
    for e in self.edges:
      i, j = e.get_vertices()
      c.add_edge( i, j, e)
    return c

  def deep_copy( self):
    """provides a deep copy of the graph. The result is an isomorphic graph,
    all the used objects are different"""
    c = self.create_graph()
    for v in self.vertices:
      c.add_vertex()
    for e in self.edges:
      i, j = e.get_vertices()
      c.add_edge( i, j)
    return c
    
  def create_vertex( self):
    return self.vertex_class()

  def create_edge( self):
    return self.edge_class()

  def create_graph( self):
    return self.__class__()
  
  def delete_vertex( self, v):
    self.vertices.remove( v)


  ## MODIFICATION METHODS

  def add_vertex( self, v=None):
    """adds a vertex to a graph, if v argument is not given creates a new one.
    returns None if vertex is already present or the vertex instance if successful"""
    if not v:
      v = self.create_vertex()
    if v not in self.vertices:
      self.vertices.append( v)
    else:
      warnings.warn( "Added vertex is already present in graph %s" % str( v), UserWarning, 2)
      return None
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
      e = self.create_edge()
    e.set_vertices( (v1,v2))
    self.edges.add( e)
    v1.add_neighbor( v2, e)
    v2.add_neighbor( v1, e)
    return e



  def insert_a_graph( self, gr):
    """inserts all edges and vertices to the graph"""
    self.vertices.extend( gr.vertices)
    self.edges |= gr.edges


  # do we need it?
  def connect_a_graph( self, gr, v1, v2, e=None):
    """gr is a graph, v1 is vertex in self, v2 is vertex in gr, bond is what to use for connection"""
    self.insert_a_graph( gr)
    self.add_edge( v1, v2, e=e)
    


  def disconnect( self, v1, v2):
    """disconnects vertices v1 and v2, on success returns the edge"""
    if v1 != None and v2 != None:
      e = self.get_edge_between( v1, v2)
      if e:
        self.edges.remove( e)
        v1.remove_neighbor( v2)
        v2.remove_neighbor( v1)
      return e
    else:
      return None



  def disconnect_edge( self, e):
    self.edges.remove( e)
    v1, v2 = e.get_vertices()
    v1.remove_edge_and_neighbor( e)
    v2.remove_edge_and_neighbor( e)


  def remove_vertex( self, v):
    for neigh in v.get_neighbors():
      self.disconnect( v, neigh)
    self.delete_vertex( v)
      


  def get_edge_between( self, v1, v2):
    """takes two vertices"""
    for e in v1.get_neighbor_edges():
      if e in v2.get_neighbor_edges():
        return e
    return None



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
      recent = Set() # [not_processed.pop()])
    while not_processed:
      recent = Set( reduce( operator.add, [a.get_neighbors() for a in recent], [])) & not_processed
      if not recent:
        if comp:
          yield comp
        recent = Set( [not_processed.pop()])
        comp = recent
      else:
        comp |= recent
        not_processed -= recent
    # when there is only one atom in the last piece it is not yielded in the loop
    yield comp


  def get_disconnected_subgraphs( self):
    vss = self.get_connected_components()
    out = []
    for vs in vss:
      out.append( self.get_induced_subgraph_from_vertices( vs))
    return out

  def get_induced_subgraph_from_vertices( self, vs):
    """it creates a new graph, however uses the old vertices and edges!"""
    g = self.create_graph()
    for v in vs:
      g.add_vertex( v)
    for e in self.edges:
      v1, v2 = e.get_vertices()
      if v1 in vs and v2 in vs:
        g.add_edge( v1, v2, e)  # BUG - it should copy the edge?
    return g
    
      
  def get_degrees( self):
    """returns a generator of degrees, this is useful because for many properties
    the whole list is not important"""
    for v in self.vertices:
      yield v.get_degree()
  
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
    return map( self.edge_subgraph_to_vertex_subgraph, self.get_smallest_independent_cycles_e())
                


  def get_smallest_independent_cycles_e( self):
    """returns a set of smallest possible independent cycles as list of Sets of edges,
    other cycles in graph are guaranteed to be combinations of them"""
    ncycles = len( self.edges) - len( self.vertices) + 1
    if ncycles < 0:
      warnings.warn( "The number of edges is smaller than number of vertices -1, the molecule must be disconnected, which means there is something wrong with it.", UserWarning, 3)
      ncycles = 0
    cycles = self.get_all_cycles_e()
    while not len( cycles) <= ncycles:
      d = len( cycles) - ncycles
      cycles = filter_off_dependent_cycles( cycles)
      if len( cycles) - ncycles == d:
        # further cycling would not help
        break
    if len( cycles) < ncycles:
      warnings.warn( "The number of cycles found (%d) is smaller than the theoretical value %d (|E|-|V|+1)" % (len( cycles), ncycles), UserWarning, 2)
    elif len( cycles) > ncycles: 
      warnings.warn( "The number of independent cycles found (%d) is larger than the theoretical value %d (|E|-|V|+1), but I cannot improve it." % (len( cycles), ncycles), UserWarning, 2)

    return cycles



  def get_almost_all_cycles_e( self):
    """returns almost all cycles found in the graph as sets of edges
    this version is not perfect as it sometimes forgets a few rings"""
    all_cycles = []
    self.mark_vertices_with_distance_from( self.vertices[0])
    to_go = Set()
    for ps in self._get_all_ring_end_points():
      to_go |= ps
    for ps in self._get_all_ring_start_points():
      to_go |= ps
    while to_go:
      v = to_go.pop()
      cycles = self._get_cycles_for_vertex( v, to_reach=v)
      all_cycles += cycles
    all_cycles = Set( map( ImSet, all_cycles))
    self.clean_distance_from_vertices()
    return all_cycles



  def get_all_cycles_e( self):
    """returns all cycles found in the graph as sets of edges"""
    to_go = Set( self.vertices)
    for v in self.vertices:
      if v.degree == 1:
        to_go.remove( v)
    all_cycles = []
#    removed = Set( to_go)
    while to_go:
##       while removed:
##         new_removed = Set()
##         for v in removed:
##           for n in v.neighbors:
##             if n in to_go and n.degree == 2:
##               new_removed.add( n)
##               to_go.remove( n)
##         removed = new_removed
      if to_go:
        v = to_go.pop()
#        removed = Set([v])
        cycles = self._get_cycles_for_vertex( v, to_reach=v)
        all_cycles += cycles
    all_cycles = Set( map( ImSet, all_cycles))
    return all_cycles




  def _get_cycles_for_vertex( self, v, to_reach=None, processed=Set()):
    for e, neigh in v.get_neighbor_edge_pairs():
      if neigh == to_reach and len( processed) > 1:
        return [Set( [e])]
    all_cycles = []
    for e, neigh in v.get_neighbor_edge_pairs():
      if neigh not in processed:
        p = Set([v]) | processed
        cycles = self._get_cycles_for_vertex( neigh, to_reach=to_reach, processed=p)
        if cycles:
          for cycle in cycles:
            cycle.add( e)
        all_cycles += cycles
    return all_cycles




  def get_all_cycles( self):
    """returns all cycles found in the graph as sets of vertices,
    use get_all_cycles_e to get the edge variant, which is better for building new
    graphs as the mapping edges => vertices is unambiguous, while edges=>vertices=>edges might
    include some more edges"""
    return map( self.edge_subgraph_to_vertex_subgraph, self.get_all_cycles_e())


      


  def mark_vertices_with_distance_from( self, v):
    """returns the maximum d"""
    self.clean_distance_from_vertices()
    return self._mark_vertices_with_distance_from( v)

  def clean_distance_from_vertices( self):
    for i in self.vertices:
      try:
        del i.properties_['d']
      except KeyError:
        pass


  def get_diameter( self):
    diameter = 0
    best = None
    best_path = None
    for v in self.vertices:
      dist = self.mark_vertices_with_distance_from( v)
      if dist > diameter:
        diameter = dist
        best = v
        end = [x for x in self.vertices if x.properties_['d'] == dist][0]
        best_path = get_path_down_to( end, v)
    print best
    print "path"
    best_path.reverse()
    for v in best_path:
      print v
    return diameter



  def vertex_subgraph_to_edge_subgraph( self, cycle):
    ret = Set()
    for v1 in cycle:
      for (e,n) in v1.get_neighbor_edge_pairs():
        if n in cycle:
          ret.add( e)
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
    sub = self.create_graph()
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

  def defines_connected_subgraph_e( self, edges):
    sub = self.get_new_induced_subgraph( self.edge_subgraph_to_vertex_subgraph( edges), edges)
    return sub.is_connected()

  def defines_connected_subgraph_v( self, vertices):
    sub = self.get_new_induced_subgraph( vertices, self.vertex_subgraph_to_edge_subgraph( vertices))
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


  ## STATIC METHODS

##   def get_random_longest_path_numbered( self, start):
##     """vertices have to be freshly marked with distance"""
##     now = start
##     path = []
##     d = 0
##     while now:
##       d += 1
##       path.append( now)
##       ns = [v for v in now.neighbors if 'd' in v.properties_ and v.properties_['d'] == d]
##       if ns:
##         now = ns[0]
##       else:
##         now = None
##     return path

##   get_random_longest_path_numbered = classmethod( get_random_longest_path_numbered)


##   def get_random_longest_path( self, start):
##     self.mark_vertices_with_distance_from( start)
##     return self.get_random_longest_path_numbered( start)

##   get_random_longest_path = classmethod( get_random_longest_path)




  # PRIVATE METHODS

  def _get_vertex_index( self, v):
    """if v is already an index, return v, otherwise return index of v on None"""
    if type( v) == IntType and v < len( self.vertices):
      return v
    try:
      return self.vertices.index( v)
    except ValueError:
      return None


  def _read_file( self, name="/home/beda/oasa/oasa/mol.graph"):
    self.vertices = []
    self.edges = Set()
    f = file( name, 'r')
    vs = f.readline()
    [self.add_vertex() for i in vs.split(" ")]
    for l in f.readlines():
      c, a, b = l.split(' ')
      self.add_edge( self.vertices[int(a)], self.vertices[int(b)])
    f.close()


  def _mark_vertices_with_distance_from( self, v):
    """returns the maximum d"""
    d = 0
    to_mark = Set([v])
    while to_mark:
      to_mark_next = Set()
      for i in to_mark:
        i.properties_['d'] = d
        for j in i.get_neighbors():
          if 'd' not in j.properties_:
            to_mark_next.add( j)

      to_mark = to_mark_next
      d += 1

    return d-1


  def _get_some_cycles( self):
    if len( self.vertices) <= 2:
      raise StopIteration
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
  """this filtres off all cycles that are combinations of smaller cycles."""
  cs = list( cycles)
  level = 2
  cs.sort( lambda a,b: -len( a) + len(b))
  while level < len( cs):
    to_del = Set()
    for c in cs:
      combs = gen_variations( [x for x in cs if x != c and x not in to_del and x & c], level)
      for combl in combs:
        comb = Set( combl)
        if not comb & to_del and (c <= reduce( operator.or_, comb)) and reduce( operator.and_, comb):
          to_del.add( c)
          break
    [cs.remove( x) for x in to_del]
    level += 1
  return cs



## def filter_off_dependent_cycles( cycles):
##   cs = list( cycles)
##   level = 2
##   cs.sort( lambda a,b: -len( a) + len(b))
##   print len( cs), map( len, cs)
##   while level < len( cs):
##     to_del = Set()
##     for c in cs:
##       combs = gen_variations( [x for x in cs if x != c and x & c], level)
##       for combl in combs:
##         comb = Set( combl)
##         if not comb & to_del and (c <= reduce( operator.or_, comb)):
##           to_del.add( c)
##           print "remove", len( c), map( len, comb)
##           break
##     [cs.remove( x) for x in to_del]
##     level += 1
##   return cs




def gen_variations(items, n):
    if n==0:
      yield []
    else:
      for i in xrange( len(items)-n+1):
        for v in gen_variations(items[i+1:],n-1):
          yield [items[i]]+v



## import profile

## g = graph()
## g._read_file()
## print g
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
