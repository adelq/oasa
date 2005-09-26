#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2005 Beda Kosata <beda@zirael.org>

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



class oasa_error( Exception):

  def __init__( self):
    Exception.__init__( self)    




class oasa_periodic_table_error( oasa_error):
  """exception for reporting periodic_table related error"""

  def __init__( self, id, value, symbol=None):
    self.id = id
    self.value = value


  def __str__( self):
    return "OASA periodic_table error, id=%s, value=%s" % (self.id, self.value)




class oasa_invalid_atom_symbol( oasa_error):
  """exception for reporting invalid atom symbol use"""

  def __init__( self, value, symbol):
    self.value = value
    self.symbol = symbol


  def __str__( self):
    return "Symbol '%s' not allowed (%s)" % (self.symbol, self.value)



class oasa_not_implemented_error( oasa_error):

  def __init__( self, where, what):
    self.where = where
    self.what = what

  def __str__( self):
    return "Not implemented error in %s: %s" % (self.where, self.what)


class oasa_inchi_error( oasa_error):

  def __init__( self, what):
    self.what = what

  def __str__( self):
    return "InChI error: %s" % self.what



class oasa_unsupported_inchi_version_error( oasa_error):
  
  def __init__( self, version):
    self.version = version

  def __str__( self):
    return "The InChI has an unsupported version: %s" % self.version



