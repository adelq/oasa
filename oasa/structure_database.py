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

import os, sys, re
try:
    from pysqlite2 import dbapi2 as sqlite
except ImportError:
    print "pysqlite module could not be loaded - please install the required package"
    raise


class Config:
    database_file = os.path.abspath( os.path.join( os.path.dirname( __file__), "structures.db"))


def create_database():
  connection = sqlite.connect( Config.database_file)
  c = connection.cursor()
  c.execute( "DROP TABLE IF EXISTS structures;")
  c.execute( """CREATE TABLE structures (
  id INTEGER PRIMARY KEY,
  name TEXT,
  inchikey TEXT,
  smiles TEXT);""")
  connection.commit()
  c.close()
  connection.close()


def normalize_inchi( inchi):
    if inchi.startswith("InChI="):
        return inchi[6:]
    else:
        return inchi

def _open_infile( infilename):
    if os.path.isfile( infilename):
        if infilename.endswith(".gz"):
            import gzip
            f = gzip.open( infilename, "r")
        else:
            f = file( infilename, "r")
        return f
    else:
        return None

def fill_database( infilename, name_cutoff=26, atom_count_cutoff=100):
    """compound is added if either name length of atom_count is below the
    corresponding cutoff value"""
    added = 0
    ignored = 0
    connection = sqlite.connect( Config.database_file)
    c = connection.cursor()
    f = _open_infile( infilename)
    if not f:
        raise ValueError( "File does not exist:", infilename)
    f.readline() # skip the first line
    i = 0
    for line in f:
        if i % 10000 == 0:
            print "done %8d, added %8d, ignored %8d" % (i,added,ignored)
        values = [x.strip() for x in line.strip().split("\t")]
        if len( values) != 4:
            print >> sys.stderr, "Ignoring line:", line,
            continue
        cid, inchikey, smiles, name = [x.strip() for x in line.strip().split("\t")]
        c.execute( "DELETE FROM structures WHERE id=?", (cid,))
        if len( name) <= name_cutoff and len( [x for x in smiles if x.isupper()]) <= atom_count_cutoff and _allow_molecule( name, smiles):
            c.execute( "INSERT INTO structures (id,name,inchikey,smiles) VALUES (?,?,?,?);", (cid, name, inchikey, smiles))
            added += 1
        else:
            ignored += 1
        if i % 100 == 0:
            connection.commit()
        i += 1
    f.close()
    connection.commit()
    connection.close()
    return added, ignored


def get_compounds_from_database( **kw):
    if 'inchi' in kw:
        if not 'inchikey' in kw:
            import inchi_key
            kw['inchikey'] = inchi_key.key_from_inchi( kw['inchi'])
        del kw['inchi']
    search = ["%s=?" % k for k in kw.keys()]
    values = kw.values()
    if search:
        sql = "SELECT * FROM structures WHERE %s" % (" AND ".join( search))
    else:
        sql = "SELECT * FROM structures"
    connection = sqlite.connect( Config.database_file)
    c = connection.cursor()
    c.execute( sql, values)
    ret = []
    for row in c:
       ret.append( row)
    c.close()
    connection.close()
    return ret

def _allow_molecule( name, smile):
    if smile.count("-]") > 2:
        # more than 2 negative charges
        return False
    if smile.count("+]") > 2:
        # more than 2 positive charges
        return False
    if re.search( "\[\d", smile):
        return False
    if smile.count(".") > 2:
        return False
    disallowed = ["Tl","U","Ra","Rn","At","Po","Bi","Ta","W","Hf","La","Y","Sr","Nb","Zr","Eu","Lu","Ho","Ce","Nd","Dy","Th","Gd","Yb","Sm",'Es','Np','Pu','Pr','Md','Pm','Pa','Tb','Bk','Cm','Er','Cf','Fm','No',"Am","Ac","Tm","Lr","Fr"]
    for ch in disallowed:
        if ch in smile:
            return False
    # try smiles
##     try:
##         smiles.text_to_mol( smile, calc_coords=False)
##     except oasa_exceptions.oasa_invalid_atom_symbol, e:
##         return False
    return True
    
def filter_src_file( infilename, name_cutoff=26, atom_count_cutoff=1000):
    """print only those lines in infilename that would be allowed
    in structure_database using current settings"""
    added = 0
    ignored = 0
    f = _open_infile( infilename)
    if not f:
        raise ValueError
    f.readline() # skip the first line
    for line in f:
        values = [x.strip() for x in line.strip().split("\t")]
        if len( values) != 4:
            print >> sys.stderr, "Ignoring line:", line,
            continue
        cid, inchikey, smiles, name = [x.strip() for x in line.strip().split("\t")]
        if len( name) <= name_cutoff and len( [x for x in smiles if x.isupper()]) <= atom_count_cutoff and _allow_molecule( name, smiles):
            print line[:-1]
            added += 1
        else:
            ignored += 1
    f.close()
    return added, ignored



if __name__ == "__main__":
    import sys
##     a,i = filter_src_file( 'selected_names_name50_inchi100.clean.txt')
##     print >> sys.stderr, a, i
##     sys.exit()

    if len( sys.argv) > 1:
        fname = sys.argv[1]
        if os.path.exists( fname):
            if not os.path.exists( Config.database_file):
                create_database()
            #try:
            added, ignored = fill_database( fname)
            print "Added: %d, Ignored: %d" % (added, ignored)
            #except:
            #    print "given file must be a text file with one compound per line in format 'InChI CID ### name'"
        else:
            print "you must supply a valid filename to update the database or no argument for a test to run"
    else:
        print get_compounds_from_database( inchi="1/C4H10/c1-3-4-2/h3-4H2,1-2H3")

