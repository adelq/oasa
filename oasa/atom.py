#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003, 2004 Beda Kosata <beda@zirael.org>

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

import sys
sys.path.append( '../')

import graph
import periodic_table as PT
from common import is_uniquely_sorted

import copy
import itertools
from warnings import warn



class atom( graph.vertex):

  def __init__( self, symbol='C', charge=0, coords=None):
    graph.vertex.__init__( self)
    self.symbol = symbol
    self.charge = charge
    self.free_sites = 0
    # None means not set (used)
    if coords:
      self.x, self.y, self.z = coords
    else:
      self.x = None
      self.y = None 
      self.z = None
    self.multiplicity = 1


  def matches( self, other):
    # query atoms
    if self.query:
      # halogens
      if self.symbol == "X":
        if other.symbol in ("F","Cl","Br","I"):
          return True
        else:
          return False
      # Q (any except H,C)
      elif self.symbol == "Q":
        if other.symbol not in "HC":
          return True
        else:
          return False
      # A (any except H)
      elif self.symbol == "A":
        if other.symbol != "H":
          return True
        else:
          return False
      # R - anything
      elif self.symbol == "R":
        return True
      else:
        raise ValueError, "The atom does not specify a query but is marked as query %s" % self 

    # normal atoms
    elif self.symbol == other.symbol and self.valency == other.valency:
      return True
    return False


  ## PROPERTIES

  # coords
  def _set_coords( self, coords):
    if len( coords) == 2:
      self.x, self.y = coords
    elif len( coords) == 3:
      self.x, self.y, self.z = coords
    else:
      raise "wrong number of coordinates"

  def _get_coords( self):
    return self.x, self.y, self.z

  coords = property( _get_coords, _set_coords, None, "atom coords")




  # charge
  def _set_charge( self, charge):
    self._charge = charge

  def _get_charge( self):
    return self._charge

  charge = property( _get_charge, _set_charge, None, "atom charge")

  # multiplicity
  def _set_multiplicity( self, multiplicity):
    self._multiplicity = multiplicity

  def _get_multiplicity( self):
    return self._multiplicity

  multiplicity = property( _get_multiplicity, _set_multiplicity, None, "atom multiplicity")


  # symbol
  def _set_symbol( self, symbol):
    try:
      self.valency = PT.periodic_table[ symbol]['valency'][0]
      self.symbol_number = PT.periodic_table[ symbol]['ord']
      if "query" in PT.periodic_table[ symbol]:
        self._query = PT.periodic_table[ symbol]["query"]
      else:
        self._query = False
    except KeyError:
      raise ValueError, "wrong atom symbol %s" % symbol
    self._symbol = symbol

  def _get_symbol( self):
    return self._symbol

  symbol = property( _get_symbol, _set_symbol, None, "atom symbol")


  # query (does the atom specify a query?)
  def _get_query( self):
    return self._query

  query = property( _get_query, None, None, "does the atom specify a query?")



  # valency
  def _set_valency( self, valency):
    self._valency = valency

  def _get_valency( self):
    return self._valency

  valency = property( _get_valency, _set_valency, None, "atoms valency")


  # occupied_valency
  def _get_occupied_valency( self):
    i = 0
    for b in self._neighbors.keys():
      ord = b.order
      if ord == 4:
        ord = 1
      i += ord

    if self.charge:
      # now we deal with charge
      if abs( self.charge) > 1:
        # charges higher than one should always decrease valency
        charge = abs( self.charge)
      elif (self.symbol in PT.accept_cation) and (self.charge == 1) and (self.valency <= PT.accept_cation[self.symbol]):
        # elements that can accept cations to increase their valency (NH4+)
        charge = -1
      elif (self.symbol in PT.accept_anion) and (self.charge == -1) and (self.valency <= PT.accept_anion[self.symbol]):
        # elements that can accept anions to increase their valency (BH4-)
        charge = -1
      else:
        # otherwise charge reduces valency 
        charge = abs( self.charge)
    else:
      charge = 0

    return i+charge+self.multiplicity-1

  occupied_valency = property( _get_occupied_valency, None, None, "atoms occupied valency")


  # free_valency
  def _get_free_valency( self):
    return self.valency - self.occupied_valency

  free_valency = property( _get_free_valency, None, None, "atoms free valency")


  # free_sites
  def _set_free_sites( self, free_sites):
    self._free_sites = free_sites

  def _get_free_sites( self):
    return self._free_sites

  free_sites = property( _get_free_sites, _set_free_sites, None, "atoms free_sites")



  # electronegativity
  def _get_electronegativity( self):
    try:
      return PT.periodic_table[self.symbol]['en']
    except KeyError:
      return None

  electronegativity = property( _get_electronegativity, None, None, "atoms electronegativity")



  # oxidation number
  def _get_oxidation_number( self):
    en = self.charge
    for e, n in self.get_neighbor_edge_pairs():
      if n.symbol != self.symbol:
        en += e.order * (n.electronegativity > self.electronegativity and 1 or -1)
    hen = PT.periodic_table['H']['en']
    en += self.free_valency * (hen > self.electronegativity and 1 or -1) 
    return en

  oxidation_number = property( _get_oxidation_number, None, None, "atoms oxidation number as text")




      
  def __str__( self):
    return "atom '%s'" % str( self.symbol)



  def get_x( self):
    return self.x or 0

  def get_y( self):
    return self.y or 0

  def get_z( self):
    return self.z or 0

  def has_aromatic_bonds( self):
    for b in self._neighbors.keys():
      if b.aromatic:
        return 1
    return 0



  def raise_valency_to_senseful_value( self):
    """set atoms valency to the lowest possible, so that free_valency
    if non-negative (when possible) or highest possible,
    does not lower valency when set to higher then necessary value"""
    while self.free_valency < 0:
      if not self.raise_valency():
        return



  def raise_valency( self):
    """used in case where valency < occupied_valency to try to find a higher one"""
    for v in PT.periodic_table[ self.symbol]['valency']:
      if v > self.valency:
        self.valency = v
        return True
    return False

    

  def is_chiral( self):
    """this code is CIP (Cahn-Ingold-Prelog) based and therefore not necessarily
    the fastest for this job, however in newer versions it will be able to take care
    of the chirality of other centres and therefor will be universal,
    for now it takes only care of the connectivity"""
    if len( self._neighbors) < 4:
      # this is not true for N,P and similar !!!
      return False
    neighs = self.get_neighbors()
    cips = []
    for a in neighs:
      cips.append( [[], a, a.gen_CIP_sequence( came_from = self)])
    while not is_uniquely_sorted( cips, cip_sorting_function):
      for cip in cips:
        try:
          cip[0].append( cip[2].next().symbol_number)
        except StopIteration:
          cip[0].append( None)
        except AttributeError:
          cip[0].append( None)  # None was yielded and it has no symbol_number :)
      # to test if we have already reached the end
      for i in range( len( cips)-1):
        if cip_sorting_function( cips[i], cips[i+1]) == 0 and cips[i][0][-3:-1] == [None, None]:
          return False
      cips.sort( cip_sorting_function)
    return True



  def get_neighbors_CIP_sorted( self):
    """return neighbors sorted according to the CIP rules"""
    neighs = self.get_neighbors()
    cips = []
    for a in neighs:
      cips.append( [[], a, a.gen_CIP_sequence( came_from = self)])
    while not is_uniquely_sorted( cips, cip_sorting_function):
      for cip in cips:
        try:
          cip[0].append( cip[2].next().symbol_number)
        except StopIteration:
          cip[0].append( None)
        except AttributeError:
          cip[0].append( None)  # None was yielded and it has no symbol_number :)
      # to test if we have already reached the end
      all_finished = len( cip[0]) > 1
      for cip in cips:
        if len( cip[0]) > 1:
          if cip[0][-3:-1] != [None, None]:
            all_finished = 0
      if all_finished:
        break # it can't be uniquely sorted, we are out of atoms
      cips.sort( cip_sorting_function)
    return [cip[1] for cip in cips]





  def gen_CIP_sequence( self, iter_over=None, came_from=None):
    """generates the CIP (Cahn-Ingold-Prelog) stream of atoms suitable for
    comparison in searches for chiral centres, their configuration etc.,
    the values in different layers (with raising distance from self) are
    separated by Nones"""
    yield self
    yield None
    neighs = self.get_neighbors()
    if came_from:
      assert came_from in neighs
      neighs.remove( came_from)
    # end if there are no neighs
    if not neighs:
      raise StopIteration
    # generators of the neighs
    cips = []
    for a in neighs:
      cips.append( [[], [], a.gen_CIP_sequence( came_from = self), 0])
      #         number_list, atom_list, cip_generator, last yielded index
    while cips:
      # this returns the already generated cip of the new first atom in later cycles
      to_remove = []
      for cip in cips:
        try:
          cs = itertools.takewhile( lambda x: x!=None, cip[2])
        except StopIteration:
          to_remove.append( cip)
          continue
        cip[1] += cs
        cip[0] += [c.symbol_number for c in cs]
      [cips.remove( cip) for cip in to_remove]
      cips.sort( cip_sorting_function)
      for cip in cips:
        for i in range( cip[3], len( cip[1])):
          yield cip[1][i]
        cip[3] = len( cip[1])
      yield None




  def gen_CIP_sequence_does_something_else( self, iter_over=None, came_from=None):
    yield self
    neighs = self.get_neighbors()
    if came_from:
      assert came_from in neighs
      neighs.remove( came_from)
    # strip the neighbors that were already processed
    #neighs = [i for i in neighs if i in iter_over]
    # end if there are no neighs
    if not neighs:
      raise StopIteration
    # remove neighs from iter_over
##     for i in neighs:
##       try:
##         iter_over.remove( i)
##       except ValueError:
##         pass
    # generators of the neighs
    cips = []
    for a in neighs:
      cips.append( [0, [], a.gen_CIP_sequence( came_from = self)])
    while cips:
      reserve = [] # here we store whats ruled out to use it later
      # this returns the already generated cip of the new first atom in later cycles
      for a in cips[0][1]:
        if a != None:
          yield a
      while cips:
        for cip in cips:
          try:
            next_a = cip[2].next()
          except StopIteration:
            cip[1].append( None)
            continue
          cip[0] += next_a.symbol_number
          cip[1].append( next_a)
        if cips[0][1] == None:
          # here we remove the cip[0] that wa marked as finished
          del cips[0]
        if cips:
          cips.sort( cip_sorting_function)
          if cips[0][1][-1] != None:
            yield cips[0][1][-1]  # the first atom in the sequence
          else:
            del cips[0]
          # we filter of the already ruled out things
          ruled_out = [c for c in cips if c[0] < cips[0][0]]
          ruled_out.reverse() # to ensure sorting in the reserve
          reserve += ruled_out
          [cips.remove( c) for c in ruled_out]
          to_yield_now = [c for c in cips if c[1][-1] == None]  # e.g. 2 of the 3 H's in Me end here
          for c in to_yield_now:
            yield c[1][-2]
          [cips.remove( c) for c in to_yield_now]
      cips = reserve
      cips.reverse() # so we ensure that the 


def cip_sorting_function( a, b):
  return -cmp( a[0], b[0])



##################################################
# TODO

# chirality for those possesing a free electron pair
