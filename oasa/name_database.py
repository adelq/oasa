#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2003-2007 Beda Kosata <beda@zirael.org>

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

import os, sys
import anydbm


class Config:
    database_file = os.path.abspath( os.path.join( os.path.dirname( __file__), "names.db"))


def normalize_inchi( inchi):
    if inchi.startswith("InChI="):
        return inchi[6:]
    else:
        return inchi
    

def compound_to_database_string( c):
    c['inchi'] = normalize_inchi( c['inchi'])
    return "%(inchi)s %(cid)s ### %(name)s\n" % c


def database_string_to_compound( line):
    a,name = line.split( "###")
    inchi, cid = a.split()
    return {'inchi':inchi.strip(), 'cid':cid.strip(), 'name':name.strip()}
    

def mydb_to_gdbm( infilename, outfilename):
    import gdbm
    infile = file( infilename, "r")
    base = gdbm.open( outfilename, "n")
    for line in infile:
        rec = database_string_to_compound( line)
        base[ rec['inchi']] = rec['cid'] + " " + rec['name']
        
    

def get_compound_from_database( inchi, database_file=None):
    inchi = normalize_inchi( inchi)
    for fname in (database_file, Config.database_file):
        if fname and os.path.exists(fname):
            break
    else:
        raise Exception("Name database not found")
    base = anydbm.open( fname)
    if base.has_key( inchi):
        cid, name = base[ inchi].split( " ", 1)
        return {'inchi': inchi, 'cid': cid, 'name': name}
    else:
        return None


def name_molecule( mol, database_file=None):
    """tries to find name for an OASA molecule in the database,
    it requires InChI generation to work"""
    import inchi
    inch = inchi.mol_to_text( mol)
    return get_compound_from_database( inch, database_file=database_file)


if __name__ == "__main__":
    import sys
    if len( sys.argv) > 1:
        fname = sys.argv[1]
        if os.path.exists( fname):
            try:
                mydb_to_gdbm( fname, Config.database_file)
            except:
                print "given file must be a text file with one compound per line in format 'InChI CID ### name'"
        else:
            print "you must supply a valid filename to update the database or no argument for a test to run"
    else:
        print get_compound_from_database( "1/C4H10/c1-3-4-2/h3-4H2,1-2H3")