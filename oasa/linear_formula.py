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

from config import Config

import re
from periodic_table import periodic_table
import coords_generator
import misc
from known_groups import name_to_smiles
import smiles


class linear_formula( object):

  def __init__( self, text="", valency=0, mol=None):
    """valency specifies the already occupied valency of the submited formula,
    it is usually used when parsing groups"""
    if text:
      self.molecule = self.parse_text( text, valency=valency, mol=mol)


  def parse_text( self, text, valency=0, mol=None):
    text = self.expand_abbrevs( text)
    mol = self.parse_form( text, valency=valency, mol=mol, reverse=False)
    if mol:
      self.molecule = mol

      # are there any atoms?
      if not self.molecule.vertices:
        return None
      
      # now we check if the structure is complete
      for v in self.molecule.vertices:
        if v.free_valency:
          return None

      if valency:
        self.molecule.remove_vertex( self.molecule.vertices[0]) # remove the dummy

      # are there any atoms, again?
      if not self.molecule.vertices:
        return None

      self.molecule.remove_all_hydrogens()
      return self.molecule



  def parse_form( self, text, valency=0, mol=None, reverse=False):
    form = text

    # the code itself
    if not mol:
      mol = Config.create_molecule()

    # create the dummy atom
    if valency:
      last_atom = mol.create_vertex()
      last_atom.valency = valency
      mol.add_vertex( last_atom)
    else:
      last_atom = None

    # check if there are branches in the formula
    if "(" not in form:
      # there are no subbranches
      chunks = re.split( "([A-Z][a-z]?[0-9]?[+-]?)", form)
      if reverse:
        chunks = self.reverse_chunks( chunks)
      for chunk in chunks:
        if chunk:
          as = self.chunk_to_atoms( chunk, mol)
          if as == None:
            return None
          for a in as:
            last_atom = self.get_last_free_atom( mol)
            mol.add_vertex( a)
            if last_atom:
              max_val = min( last_atom.free_valency, a.free_valency, 3)
              b = mol.create_edge()
              b.order = max_val
              mol.add_edge( last_atom, a, b)
    else:
      for chunk, count in gen_formula_fragments( form, reverse=reverse):
        if chunk:
          last_atom = self.get_last_free_atom( mol)

          for j in range( count):
            if chunk[0] == "!":
              # the form should be a smiles
              m = smiles.text_to_mol( chunk[1:], calc_coords=0)
              m.add_missing_hydrogens()
              hs = [v for v in m.vertices[0].neighbors if v.symbol == 'H']
              m.disconnect( hs[0], m.vertices[0])
              m.remove_vertex( hs[0])
              smile = True
            else:
              val = last_atom and 1 or 0
              m = self.parse_form( chunk, valency=val, mol=mol.create_graph())
              smile = False
            if not last_atom:
              mol.insert_a_graph( m) 
            else:
              if not smile:
                m.remove_vertex( m.vertices[0]) # remove the dummy
              mol.insert_a_graph( m)
              b = mol.create_edge()
              mol.add_edge( last_atom, m.vertices[0], b)
                          

    return mol
        


  def chunk_to_atoms( self, chunk, mol):
    m = re.match( "([A-Z][a-z]?)([0-9])?([+-])?", chunk)    
    if m:
      name = m.group( 1)
      number = m.group( 2) and int( m.group(2)) or 1
      sign = m.group( 3) and (int( m.group(3)+'1')) or 0
      ret = []
      for i in range( number):
        v = mol.create_vertex()
        try:
          v.symbol = name
        except ValueError:
          # wrong atom symbol
          return None
        v.charge = sign
        ret.append( v)
      return ret


  def get_last_free_atom( self, mol):
    # check if there is something with a free valency
    if not mol.vertices:
      return None
    atoms = [o for o in misc.reverse( mol.vertices)]
    for a in atoms:
      if a.free_valency > 0:
        return a
    # if its not the case
    for i, a in enumerate( atoms[0:-1]):
      b = a.get_edge_leading_to( atoms[i+1])
      if b and b.order > 1:
        b.order -= 1
        return a
    # well, we cannot do anything else
    return atoms[-1]


  def expand_abbrevs( self, text):
    for key, val in name_to_smiles.iteritems():
      text = text.replace( key, "(!%s)" % val)
    return text


  def reverse_chunks( self, chunks):
    chunks.reverse()
    for i in range( 0, len( chunks), 2):
      chunks[i], chunks[i+1] = chunks[i+1],chunks[i]
    return chunks


def gen_formula_fragments( formula, reverse=False):
  chunks = list( gen_formula_fragments_helper( formula))
  if reverse:
    chunks.reverse()
  i = 0
  while i < len( chunks):
    chunk, brack = chunks[ i]
    if brack and i < len( chunks) - 1:
      next, nbrack = chunks[i+1]
      if not nbrack:
        count, rest = split_number_and_text( next)
        if count == None:
          yield chunk, 1
        else:
          chunks[i+1] = (rest, nbrack)
          if not rest:
            i += 1
          yield chunk, count
      else:
        yield chunk, 1
    else:
      yield chunk, 1
    i += 1


def gen_formula_fragments_helper( formula):
  opened_brackets = 0
  to_ret = []
  for ch in formula:
    if ch not in "()":
      to_ret.append( ch)
    elif ch == "(":
      if opened_brackets == 0:
        if to_ret:
          yield ''.join( to_ret), False
          to_ret = []
      else:
        to_ret.append( "(")
      opened_brackets += 1
    elif ch == ")":
      opened_brackets -= 1
      if opened_brackets == 0:
        if to_ret:
          yield ''.join( to_ret), True
          to_ret = []
      else:
        to_ret.append( ")")
  if to_ret:
    yield ''.join( to_ret), False



def split_number_and_text( txt):
  last = None
  for i in range( len( txt)):
    try:
      last = int( txt[0:i+1])
    except ValueError:
      return last, txt[i:]
  return last, ""
    
  


## form = 'CPh4'
## #form = "CH2(Cl)2"

## #print [i for i in gen_formula_fragments_helper( form)]
## #print [i for i in gen_formula_fragments( form)]

## a = linear_formula( form , valency=0)
## m = a.molecule
## #coords_generator.calculate_coords( m)

## print m

## import smiles
## print form

## if m:
##   print smiles.mol_to_text( m)

## #coords_generator.show_mol( m)


## #print [i for i in gen_formula_fragments( "CO(OH)2")]
