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


from molecule import molecule
from atom import atom
from bond import bond

import re
from periodic_table import periodic_table
import coords_generator
import misc


class linear_formula( object):

  def __init__( self, text="", valency=0, mol=None):
    """valency specifies the already occupied valency of the submited formula,
    it is usually used when parsing groups"""
    if text:
      self.molecule = self.parse_text( text, valency=valency, mol=mol)


  def parse_text( self, text, valency=0, mol=None):
    mol = self.parse_form(  text, valency=valency, mol=mol)
    if mol:
      self.molecule = mol

      # now we check if the structure is complete
      for v in self.molecule.vertices:
        if v.free_valency:
          return None

      if valency:
        self.molecule.remove_vertex( self.molecule.vertices[0]) # remove the dummy

      self.molecule.remove_all_hydrogens()
      return self.molecule



  def parse_form( self, text, valency=0, mol=None):
    is_formula = re.compile("^(([A-Z][a-z]?[0-9]*[+-]?)|\(|\)[0-9]?)*$")
    form = text
    if not is_formula.match( form):
      return None

    # the code itself
    if not mol:
      mol = molecule()

    # create the dummy atom
    if valency:
      last_atom = mol.create_vertex()
      last_atom.valency = valency
      mol.add_vertex( last_atom)
    else:
      last_atom = None

    # check if there are branches in the formula
    branching = re.split( "\(|\)([0-9]?)", form)
    if len( branching) == 1:
      # there are no subbranches
      chunks = re.split( "([A-Z][a-z]?[0-9]?[+-]?)", form)
      for chunk in chunks:
        if chunk:
          as = self.chunk_to_atoms( chunk, mol)
          for a in as:
            mol.add_vertex( a)
            if last_atom:
              max_val = min( last_atom.free_valency, a.free_valency, 3)
              b = mol.create_edge()
              b.order = max_val
              mol.add_edge( last_atom, a, b)
          last_atom = self.get_last_free_atom( mol)
    else:
      if len( branching) % 2:
        branching.append( '')
      for i in range( 0, len( branching), 2):
        chunk = branching[i]
        if chunk:
          rep = branching[i+1] and int( branching[i+1]) or 1
          for j in range( rep):
            m = self.parse_form( chunk, valency=1, mol=mol.create_graph())
            if not last_atom:
              mol.insert_a_graph( m) 
            else:
              m.remove_vertex( m.vertices[0]) # remove the dummy
              mol.insert_a_graph( m) 
              b = mol.create_edge()
              mol.add_edge( last_atom, m.vertices[0], b)
                          
        last_atom = self.get_last_free_atom( mol)

    return mol
        


  def chunk_to_atoms( self, chunk, mol):
    m = re.match( "([A-Z][a-z]?)([0-9])?([+-])?", chunk)    
    name = m.group( 1)
    number = m.group( 2) and int( m.group(2)) or 1
    sign = m.group( 3) and (int( m.group(3)+'1')) or 0
    ret = []
    for i in range( number):
      v = mol.create_vertex()
      v.symbol = name
      v.charge = sign
      ret.append( v)
    return ret


  def get_last_free_atom( self, mol):
    # check if there is something with a free valency
    atoms = [o for o in misc.reverse( mol.vertices)]
    for a in atoms:
      if a.get_free_valency() > 0:
        return a
    # if its not the case
    for i, a in enumerate( atoms):
      b = a.get_edge_leading_to( atoms[i+1])
      if b.order > 1:
        b.order -= 1
        return a
    # well, we cannot do anything else
    return a


## #form = 'OCH(COOCH3)2'
## form = "COO-"

## a = linear_formula( form , valency=1)
## m = a.molecule
## #coords_generator.calculate_coords( m)


## import smiles
## print form
## print smiles.mol_to_text( m)

## #coords_generator.show_mol( m)
