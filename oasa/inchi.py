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

import re
from sets import Set
import operator
import time
import xml.dom.minidom as dom
import dom_extensions
import string


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
    self.read_sum_layer( layers[1])
    self.read_connectivity_layer( layers[2])
    self.read_hydrogen_layer( layers[3])
    self.structure.add_missing_bond_orders()
    self.structure.remove_all_hydrogens()



  def read_sum_layer( self, layer):
    form = pt.formula_dict( layer)
    for k in form.sorted_keys():
      if k == 'H':
        continue
      for i in range( form[k]):
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
    for c in chunks:
      as, ns = re.split( "H", c)
      if '-' in as:
        a1, a2 = re.split( '-', as)
      else:
        a1 = as
        a2 = as
      a1 = int( a1)
      a2 = int( a2)
      ns = ns or 1
      ns = int( ns)
      for i in range( a1, a2+1):
        for j in range( ns):
          h = self.structure.add_vertex( atom( symbol='H'))
          self.structure.add_edge( i-1, h, e=bond())
        


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
      
      
    




##################################################
# MODULE INTERFACE

reads_text = 1
reads_files = 1
writes_text = 0
writes_files = 0

def text_to_mol( text, include_hydrogens=False, mark_aromatic_bonds=False):
  inc = inchi()
  inc.read_inchi( text)
  mol = inc.structure
  #mol.add_missing_bond_orders()
  #if not include_hydrogens:
  #  mol.remove_all_hydrogens()
  if mark_aromatic_bonds:
    mol.mark_aromatic_bonds()

  return mol


def file_to_mol( f):
  return text_to_mol( f.read())


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
  inch = "1.0Beta/C16H22NP/1(2-4-8-15-9-6-7-10-15)3-5-11-16-12-17-14-18-13-16/1-5H2,6-7H,8H2,9-10H,11H2,12-15H"
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
