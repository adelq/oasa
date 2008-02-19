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


"""Here resides the abstract parent class for all converters"""


class converter_base( object):

  # standard converter attrs
  reads_text = False
  writes_text = False
  reads_files = False
  writes_files = False

  default_configuration = {}

  def __init__( self):
    self.configuration = {}
    for k,v in self.default_configuration.iteritems():
      self.configuration[k] = v
    self.clean_logs()

  def clean_logs( self):
    # here all warnings for the last run should be stored
    self.warnings = []
    # here error messages for the last run will be stored
    self.errors = []

  def mols_to_text( self, structures):
    self.clean_logs()
    pass

  def text_to_mols( self, text):
    self.clean_logs()
    pass

  def mols_to_file( self, structures, f):
    self.clean_logs()
    pass

  def file_to_mols( self, f):
    self.clean_logs()
    pass

