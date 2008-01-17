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
#
#
#
#--------------------------------------------------------------------------

from __future__ import generators
#from generators import lazy_map, take
import re
import operator
import types


"""periodic table as a dictionary, plus functions for molecular
formula manipulation and computation"""

periodic_table = {
  "H"  : {'weight':   1.0079,  'ord':  1,  'en': 2.10,  'els':  1,  'valency': (1,)},
  "He" : {'weight':   4.0026,  'ord':  2,  'en': 0.00,  'els':  8,  'valency': (0, 2)},
  "Li" : {'weight':   6.9410,  'ord':  3,  'en': 0.98,  'els':  1,  'valency': (1,)},
  "Be" : {'weight':   9.0122,  'ord':  4,  'en': 1.57,  'els':  2,  'valency': (2,)},
  "B"  : {'weight':  10.8110,  'ord':  5,  'en': 2.04,  'els':  3,  'valency': (3,)},
  "C"  : {'weight':  12.0107,  'ord':  6,  'en': 2.55,  'els':  4,  'valency': (4, 2)},
  "N"  : {'weight':  14.0067,  'ord':  7,  'en': 3.04,  'els':  5,  'valency': (3, 5)},
  "O"  : {'weight':  15.9994,  'ord':  8,  'en': 3.44,  'els':  6,  'valency': (2,)},
  "F"  : {'weight':  18.9984,  'ord':  9,  'en': 3.98,  'els':  7,  'valency': (1,)},
  "Ne" : {'weight':  20.1797,  'ord': 10,  'en': 0.00,  'els':  8,  'valency': (0, 2)},
  "Na" : {'weight':  22.9898,  'ord': 11,  'en': 0.93,  'els':  1,  'valency': (1,)},
  "Mg" : {'weight':  24.3050,  'ord': 12,  'en': 1.31,  'els':  2,  'valency': (2,)},
  "Al" : {'weight':  26.9815,  'ord': 13,  'en': 1.61,  'els':  3,  'valency': (3,)},
  "Si" : {'weight':  28.0855,  'ord': 14,  'en': 1.91,  'els':  4,  'valency': (4,)},
  "P"  : {'weight':  30.9738,  'ord': 15,  'en': 2.19,  'els':  5,  'valency': (3, 5)},
  "S"  : {'weight':  32.0650,  'ord': 16,  'en': 2.58,  'els':  6,  'valency': (2, 4, 6)},
  "Cl" : {'weight':  35.4530,  'ord': 17,  'en': 3.16,  'els':  7,  'valency': (1, 3, 5, 7)},
  "Ar" : {'weight':  39.9480,  'ord': 18,  'en': 0.00,  'els':  8,  'valency': (0, 2)},
  "K"  : {'weight':  39.0983,  'ord': 19,  'en': 0.82,  'els':  1,  'valency': (1,)},
  "Ca" : {'weight':  40.0780,  'ord': 20,  'en': 1.00,  'els':  2,  'valency': (2,)},
  "Sc" : {'weight':  44.9559,  'ord': 21,  'en': 1.36,  'els':  3,  'valency': (3, 1)},
  "Ti" : {'weight':  47.8670,  'ord': 22,  'en': 1.54,  'els':  4,  'valency': (4, 3)},
  "V"  : {'weight':  50.9415,  'ord': 23,  'en': 1.63,  'els':  5,  'valency': (2, 4, 5)},
  "Cr" : {'weight':  51.9961,  'ord': 24,  'en': 1.66,  'els':  6,  'valency': (2, 3, 6)},
  "Mn" : {'weight':  54.9380,  'ord': 25,  'en': 1.55,  'els':  7,  'valency': (2, 3, 4, 6, 7)},
  "Fe" : {'weight':  55.8450,  'ord': 26,  'en': 1.83,  'els':  8,  'valency': (0, 2, 3)},
  "Co" : {'weight':  58.9332,  'ord': 27,  'en': 1.88,  'els':  8,  'valency': (2, 3)},
  "Ni" : {'weight':  58.6934,  'ord': 28,  'en': 1.91,  'els':  8,  'valency': (2, 3)},
  "Cu" : {'weight':  63.5460,  'ord': 29,  'en': 1.90,  'els':  1,  'valency': (0, 1, 2)},
  "Zn" : {'weight':  65.3900,  'ord': 30,  'en': 1.65,  'els':  2,  'valency': (2,)},
  "Ga" : {'weight':  69.7230,  'ord': 31,  'en': 1.81,  'els':  3,  'valency': (3,)},
  "Ge" : {'weight':  72.6400,  'ord': 32,  'en': 2.01,  'els':  4,  'valency': (4,)},
  "As" : {'weight':  74.9216,  'ord': 33,  'en': 2.18,  'els':  5,  'valency': (3, 5)},
  "Se" : {'weight':  78.9600,  'ord': 34,  'en': 2.55,  'els':  6,  'valency': (2, 4, 6)},
  "Br" : {'weight':  79.9040,  'ord': 35,  'en': 2.96,  'els':  7,  'valency': (1, 3, 5)},
  "Kr" : {'weight':  83.8000,  'ord': 36,  'en': 0.00,  'els':  8,  'valency': (0, 2)},
  "Rb" : {'weight':  85.4678,  'ord': 37,  'en': 0.82,  'els':  1,  'valency': (1,)},
  "Sr" : {'weight':  87.6200,  'ord': 38,  'en': 0.95,  'els':  2,  'valency': (2,)},
  "Y"  : {'weight':  88.9059,  'ord': 39,  'en': 1.22,  'els':  3,  'valency': (3,)},
  "Zr" : {'weight':  91.2240,  'ord': 40,  'en': 1.33,  'els':  4,  'valency': (4,)},
  "Nb" : {'weight':  92.9064,  'ord': 41,  'en': 1.60,  'els':  5,  'valency': (3, 5)},
  "Mo" : {'weight':  95.9400,  'ord': 42,  'en': 2.16,  'els':  6,  'valency': (3, 5, 6)},
  "Tc" : {'weight':  98.9063,  'ord': 43,  'en': 1.90,  'els':  7,  'valency': (5, 7)},
  "Ru" : {'weight': 101.0700,  'ord': 44,  'en': 2.20,  'els':  8,  'valency': (3, 4, 6, 8)},
  "Rh" : {'weight': 102.9055,  'ord': 45,  'en': 2.28,  'els':  8,  'valency': (3, 4)},
  "Pd" : {'weight': 106.4200,  'ord': 46,  'en': 2.20,  'els':  8,  'valency': (2, 4)},
  "Ag" : {'weight': 107.8682,  'ord': 47,  'en': 1.93,  'els':  1,  'valency': (1,)},
  "Cd" : {'weight': 112.4110,  'ord': 48,  'en': 1.69,  'els':  2,  'valency': (2,)},
  "In" : {'weight': 114.8180,  'ord': 49,  'en': 1.78,  'els':  3,  'valency': (3,)},
  "Sn" : {'weight': 118.7100,  'ord': 50,  'en': 1.96,  'els':  4,  'valency': (2, 4)},
  "Sb" : {'weight': 121.7600,  'ord': 51,  'en': 2.05,  'els':  5,  'valency': (3, 5)},
  "Te" : {'weight': 127.6000,  'ord': 52,  'en': 2.10,  'els':  6,  'valency': (2, 4, 6)},
  "I"  : {'weight': 126.9045,  'ord': 53,  'en': 2.66,  'els':  7,  'valency': (1, 3, 5, 7)},
  "Xe" : {'weight': 131.2930,  'ord': 54,  'en': 2.60,  'els':  8,  'valency': (0, 2)},
  "Cs" : {'weight': 132.9055,  'ord': 55,  'en': 0.79,  'els':  1,  'valency': (1,)},
  "Ba" : {'weight': 137.2370,  'ord': 56,  'en': 0.89,  'els':  2,  'valency': (2,)},
  "La" : {'weight': 138.9055,  'ord': 57,  'en': 1.10,  'els':  3,  'valency': (3,)},
  "Hf" : {'weight': 178.4900,  'ord': 72,  'en': 1.30,  'els':  4,  'valency': (4,)},
  "Ta" : {'weight': 180.9479,  'ord': 73,  'en': 1.50,  'els':  5,  'valency': (5,)},
  "W"  : {'weight': 183.8400,  'ord': 74,  'en': 2.36,  'els':  6,  'valency': (6,)},
  "Re" : {'weight': 186.2070,  'ord': 75,  'en': 1.90,  'els':  7,  'valency': (7,)},
  "Os" : {'weight': 190.2300,  'ord': 76,  'en': 2.20,  'els':  8,  'valency': (4, 6, 8)},
  "Ir" : {'weight': 192.2170,  'ord': 77,  'en': 2.20,  'els':  8,  'valency': (3, 4, 6)},
  "Pt" : {'weight': 195.0780,  'ord': 78,  'en': 2.28,  'els':  8,  'valency': (2, 4)},
  "Au" : {'weight': 196.9665,  'ord': 79,  'en': 2.54,  'els':  1,  'valency': (1, 3)},
  "Hg" : {'weight': 200.5900,  'ord': 80,  'en': 2.00,  'els':  2,  'valency': (1, 2)},
  "Tl" : {'weight': 204.3833,  'ord': 81,  'en': 2.04,  'els':  3,  'valency': (1, 3)},
  "Pb" : {'weight': 207.2000,  'ord': 82,  'en': 2.33,  'els':  4,  'valency': (2, 4)},
  "Bi" : {'weight': 208.9804,  'ord': 83,  'en': 2.02,  'els':  5,  'valency': (3, 5)},
  "Po" : {'weight': 208.9824,  'ord': 84,  'en': 2.00,  'els':  6,  'valency': (2, 4, 6)},
  "At" : {'weight': 209.9871,  'ord': 85,  'en': 2.20,  'els':  7,  'valency': (1, 7)},
  "Rn" : {'weight': 222.0176,  'ord': 86,  'en': 0.00,  'els':  8,  'valency': (0, 2)},
  "Ra" : {'weight': 226.0254,  'ord': 88,  'en': 0.90,  'els':  2,  'valency': (2,)},
  "U"  : {'weight': 238.0289,  'ord': 92,  'en': 1.38,  'els':  6,  'valency': (3, 4, 5, 6)},

  # query atoms
  "X":  {"weight": 0, 'ord': 300, "valency": (1,),       "query": True},  # halogen
  "Q":  {"weight": 0, 'ord': 301, "valency": (1,2,3,4),  "query": True},  # anything not H or C
  "A":  {"weight": 0, 'ord': 302, "valency": (1,2,3,4),  "query": True},  # anything not H
  "R":  {"weight": 0, 'ord': 303, "valency": (1,2,3,4),  "query": True},  # anything
  
  }

# elements that accept cations and raise their valency; for each element a valency is specified
# because this property is valency (oxidation state) specific 
accept_cation = {'N': 3, 'P': 3, 'O': 2, 'S': 2, 'Se': 2}

# elements that accept anions and raise their valency; for each element a valency is specified
# because this property is valency (oxidation state) specific 
accept_anion = {'B': 3, 'Al': 3, 'P': 5}


class composition_dict( dict):
  """special dict that automatically converts itself to human readable composition on str()"""
  def __str__( self):
    ret = ''
    for n in ('C','H'):
      if n in self:
        if ret:
          ret += ', '
        ret += "%s: %2.3f%%" % (n, self[n])
    k = self.keys()
    k.sort()
    for n in self:
      if n not in ('C','H'):
        if ret:
          ret += ', '
        ret += "%s: %2.3f%%" % (n, self[n])
    return ret

class formula_dict( dict):
  """special dict that automatically converts itself to human readable
  formula on str(). Implements += for convenient formula concatenation"""

  def __init__( self, form=None):
    dict.__init__( self)
    ## incomplete means that there were some problems to fully convert a formula to this dict
    self.incomplete = 0
    if type( form) in (types.StringType, types.UnicodeType):
      self.read_formula_string( form)
    elif type( form) == types.DictType:
      for key, val in form.iteritems():
        if key in periodic_table and type( val) == types.IntType:
          self[ key] = val
        else:
          raise ValueError, "some of the dictionary entries are not valid for formula_dict (%s => %s)" % (str( key), str( val))
  
  def __str__( self, reverse=0):
    sum = ''
    k = self.sorted_keys()
    if reverse:
      k.reverse()
    for s in k:
      if self[s] == 1:
        num = ''
      else:
        num = str( self[s])
      sum += s+num
    return sum

  def __iadd__( self, other):
    for s in other:
      if s in self:
        self[s] += other[s]
      else:
        self[s] = other[s]
    return self

  def __add__( self, other):
    ret = formula_dict()
    for form in (self, other):
      for s in form:
        if s in ret:
          ret[s] += form[s]
        else:
          ret[s] = form[s]
    return ret

  def __mul__( self, other):
    if not type( other) == types.IntType:
      raise TypeError, "formula_dict can be only multiplied by an integer"
    res = formula_dict()
    for key in self.keys():
      res[key] = other * self[key]
    return res


  def get_element_fraction( self, element):
    if element in self:
      return self[element]*periodic_table[element]['weight']/self.get_molecular_weight()
    return 0

  def get_molecular_weight( self):
    tot = 0
    for i in self:
      tot += self[i]* periodic_table[i]['weight']
    return tot
  
  def keys_in_order( self):
    return self.sorted_keys()

  def sorted_keys( self):
    k = self.keys()
    ret = []
    if 'C' in k:
      for a in ('C','H'):
        if a in k:
          ret.append( a)
          k.remove( a)
      k.sort()
      return ret+k
    else:
      k.sort()
      return k
    

  def read_formula_string( self, form):
    is_formula = re.compile("^([A-Z][a-z]?[0-9]*)*$")
    #form = "".join( form.split("."))
    form = form.replace( ".", "")
    if not is_formula.match( form):
      return None
    chunks = re.split( "([A-Z][a-z]*)", form)
    del chunks[0]
    for i in range( 0, len( chunks), 2):
      if chunks[i] in self:
        if chunks[i+1] == '':
          j = 1
        else:
          j = int( chunks[i+1])
        self[ chunks[i]] += j
      elif chunks[i] in periodic_table:
        if chunks[i+1] == '':
          j = 1
        else:
          j = int( chunks[i+1])
        self[ chunks[i]] = j
      else:
        self.incomplete = 1

  def get_html_repr_as_string( self, outer_element=None, reverse=0):
    sum = ''
    k = self.sorted_keys()
    if reverse:
      k.reverse()
    for s in k:
      if self[s] == 1:
        num = ''
      else:
        num = '<sub>%d</sub>' % self[s]
      sum += s+num
    if outer_element:
      return '<%s>%s</%s>' % (outer_element, sum, outer_element)
    return sum
    
  def is_saturated_alkyl_chain( self):
    if (self.sorted_keys() == ['C','H']) and (self['H'] == 2*self['C']+1):
      return 1
    else:
      return 0

  def to_tuple( self):
    return tuple( reduce( operator.add, [[k,self[k]] for k in self.sorted_keys()], []))


    
def dict_to_composition( form):
  w = form.get_molecular_weight()
  ret = composition_dict()
  for s in form:
    ret[ s] = form.get_element_fraction(s) * 100
  return ret

def formula_to_weight( formula):
  return formula_dict( formula).get_molecular_weight()

def formula_to_formula( formula):
  return str( formula_dict( formula))

def formula_to_composition( formula):
  return dict_to_composition( formula_to_dict( formula))


## other support functions

def text_to_hydrogenated_atom( text):
  a = re.match( '^([a-z]{1,2})(h)(\d*)$', text.lower())
  if a:
    atom = a.group( 1)
    hydrogens = a.group( 3)
  else:
    a = re.match( '^(h)(\d*)([a-z]{1,2})$', text.lower())
    if a:
      atom = a.group( 3)
      hydrogens = a.group( 2)
    else:
      return None

  if atom.capitalize() in periodic_table:
    ret = formula_dict()
    ret[ atom.capitalize()] = 1
    if hydrogens:
      ret[ 'H'] = int( hydrogens)
    else:
      ret[ 'H'] = 1
    return ret
  else:
    return None



def gen_bit_masks( length):
  ret = length * [0]
  yield ret
  for i in xrange( 2 ** length):
    ret[0] += 1
    for j in xrange( length):
      if ret[j] == 2:
        ret[j] = 0
        if j == length-1:
          raise StopIteration
        else:
          ret[j+1] += 1
      else:
        break
    yield ret
  

