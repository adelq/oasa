#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003-2008 Beda Kosata <beda@zirael.org>

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

import oasa_exceptions


class stereochemistry( object):

  def __init__( self, center=None, references=None, value=None):
    self.center = center
    if not references:
      self.references = []
    else:
      self.references = references
    self.value = value

  # value property
  def _get_value( self):
    return self._value
  def _set_value( self, value):
    self._value = value
  value = property( _get_value, _set_value)

  # references property
  def _get_references( self):
    return self._references
  def _set_references( self, references):
    self._references = references
  references = property( _get_references, _set_references)

  # center property
  def _get_center( self):
    return self._center
  def _set_center( self, center):
    self._center = center
  center = property( _get_center, _set_center)



class cis_trans_stereochemistry( stereochemistry):

  SAME_SIDE = 0
  OPPOSITE_SIDE = 1

  def _set_value( self, value):
    if value not in (self.SAME_SIDE, self.OPPOSITE_SIDE):
      raise oasa_exceptions.oasa_stereochemistry_error( "invalid stereochemistry identifier '%s'" % value)
    super( self.__class__)._set_value( value)

  def _set_references( self, references):
    if len( references) != 2:
      raise oasa_exceptions.oasa_stereochemistry_error( "wrong number of references in stereochemistry specification '%s'" % len( references))
    super( self.__class__)._set_references( references)    

  def get_other( self, ref):
    if not ref in self.references:
      raise ValueError, "submitted object is not referenced in this stereochemistry object."
    ref1, ref2 = self.references
    return ref is ref1 and ref2 or ref1