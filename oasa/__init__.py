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

import sys
if not (sys.version_info[0] > 2 or (sys.version_info[0] == 2 and sys.version_info[1] >= 3)):
  raise ImportError, "system version %d.%d is lower than 2.3 which is needed by OASA" % sys.version_info[0:2]


import atom
import bond
import molecule
import smiles
import coords_generator
import coords_optimizer
import molfile
import inchi
import cdml
import graph
import linear_formula
import periodic_table
import config

atom = atom.atom
molecule = molecule.molecule
bond = bond.bond

__all__ = ['atom','bond','molecule','smiles','coords_generator','molfile','inchi','graph',"linear_formula",'periodic_table','config',
           'coords_optimizer']



  
