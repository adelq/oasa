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



class molfile( plugin):

  name = "molfile"
  read = 1
  write = 1

  def __init__( self, structure=None):
    self.structure = structure

  def set_structure( self, structure):
    self.structure = structure

  def get_structure( self):
    return self.structure

  def read_file( self, file):
    self._read_header( file)
    self._read_body( file)

  def _read_header( self, file):
    for i in range( 3):
      file.readline()

  def _read_body( self, file):
    atoms = read_molfile_value( file, 3, conversion=int)
    bonds = read_molfile_value( file, 3, conversion=int)
    # nothing more interesting
    file.readline()
    # read the structure
    self.structure = molecule()
    for i in range( atoms):
      a = self._read_atom( file)
      self.structure.add_vertex( v=a)
    for k in range( bonds):
      b, i, j = self._read_bond( file)
      self.structure.add_edge( i, j, e=b)

  def _read_atom( self, file):
    x = read_molfile_value( file, 10, conversion=float)
    y = read_molfile_value( file, 10, conversion=float)
    z = read_molfile_value( file, 10, conversion=float)
    read_molfile_value( file, 1) # empty space
    symbol = read_molfile_value( file, 3)
    mass_diff = read_molfile_value( file, 2)
    charge, multi = self._read_molfile_charge( read_molfile_value( file, 3, conversion=int))
    file.readline() # next line please
    ret = self.structure.create_vertex()
    ret.coords = (x,y,z)
    ret.charge = charge
    ret.symbol = symbol
    ret.multiplicity = multi
    return ret

  def _read_bond( self, file):
    a1 = read_molfile_value( file, 3, conversion=int) -1 # molfiles index from 1
    a2 = read_molfile_value( file, 3, conversion=int) -1 # how lamish :)
    order = read_molfile_value( file, 3, conversion=int)
    type = read_molfile_value( file, 3, conversion=int)
    type_remap = { 0: 'n', 1: 'w', 6: 'h', 4: 'a'}
    type = type_remap[ type]
    file.readline() # next line please
    b = self.structure.create_edge()
    b.order = order
    b.type = type
    return b, a1, a2

  def write_file( self, file):
    """file should be a writable file object"""
    if not self.structure:
      raise "No structure to write"
    self._write_header( file)
    self._write_counts_line( file)
    self._write_body( file)
    file.write( 'M  END\n')

  def get_text( self):
    return ''.join( (self._get_header(), self._get_counts_line(), self._get_body(), 'M  END'))

  def _write_header( self, file):
    file.write( self._get_header())

  def _get_header( self):
    return "\n\n\n"

  def _write_counts_line( self, file):
    file.write( self._get_counts_line())

  def _get_counts_line( self):
    atoms = len( self.structure.vertices)
    bonds = len( self.structure.edges)
    atom_lists = 0
    fff = 0 # obsolete
    chiral = 0
    stexts = 0
    obsolete = "  0  0  0  0"
    extras = 999
    mol_version = " V2000"
    #         1  2 3 4 5 6 7 8 9
    return "%3d%3d%3d%3d%3d%3d%s%s%s\n" % (atoms,bonds,atom_lists,fff,chiral,stexts,obsolete,extras,mol_version)

  def _write_body( self, file):
    file.write( self._get_body())

  def _get_body( self):
    lines = []
    for a in self.structure.vertices:
      lines.append( self._get_atom_line( a))
    for b in self.structure.edges:
      lines.append( self._get_bond_line( b))
    return '\n'.join( lines)+'\n'


  def _get_atom_line( self, a):
    x = a.get_x()
    y = a.get_y()
    z = a.get_z()
    symbol = a.symbol
    mass_diff = 0
    charge = self._get_molfile_charge( a.charge, a.multiplicity)
    rest = "  0  0  0  0  0  0  0  0  0  0"
    #            1    2     3     4  5  6 7
    return "%10.4f%10.4f%10.4f %-3s%2d%3d%s" % (x,y,z,symbol,mass_diff,charge,rest)

  def _get_bond_line( self, b):
    v1, v2 = b.get_vertices()
    a1 = self.structure.vertices.index( v1) +1
    a2 = self.structure.vertices.index( v2) +1
    order = b.order
    type_remap = { 'n': 0, 'w': 1, 'h': 6, 'a': 4, 'b': 0, 'd': 0}
    type = type_remap[ b.type]
    rest = "  0  0  0"
    #         1  2  3  4 5
    return "%3d%3d%3d%3d%s" % (a1,a2,order,type,rest)


  def _read_molfile_charge( self, value):
    if value == 0:
      return (0,1)
    ch = -(value-4)
    multi = ch == 0 and 2 or 1
    return (ch, multi)

  def _get_molfile_charge( self, charge, multiplicity):
    if charge == 0:
      if multiplicity == 1:
        return 0
      else:
        return 4
    else:
      return 4-charge


def read_molfile_value( file, length, strip=1, conversion=None):
  """reads specified number of characters, if strip strips whitespace,
  if conversion (a fuction taking one string argument) applies it;
  if empty string is obtained after stripping 0 is returned"""
  str = file.read( length)
  if strip:
    str = str.strip()
  if str == "":
    return 0
  if conversion:
    str = conversion( str)
  return str

  


##################################################
# MODULE INTERFACE

from StringIO import StringIO

reads_text = 1
reads_files = 1
writes_text = 1
writes_files = 1

def mol_to_text( mol):
  m = molfile()
  m.structure = mol
  return m.get_text()

def mol_to_file( mol, f):
  m = molfile()
  m.structure = mol
  m.write_file( f)

def file_to_mol( f):
  m = molfile()
  m.read_file( f)
  return m.structure

def text_to_mol( text):
  file_to_mol( StringIO( text))

#
##################################################
  

