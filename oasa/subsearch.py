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

import smiles
from graph.digraph import digraph

class substructure_search_manager( object):

  def __init__( self):
    self.structures = digraph()  # graph describing the relations between individual structures
    self.search_trees = []  # list of substructure instances that for a tree

  def read_structure_file( self, name="subsearch_data.txt"):
    f = file( name, "r")
    for line in f:
      l = line.strip()
      if l and not l.startswith( "#"):
        parts = l.split(";")
        if len( parts) < 3:
          raise ValueError( "wrong line in data file: '%s'" % l)
        compound_type, name, smiles_string = parts[:3]
        if len( parts) > 3:
          to_ignore = parts[3].strip()
        else:
          to_ignore = ""
        to_ignore = map( int, filter( None, to_ignore.split(",")))
        if not name.strip():
          name = compound_type
        sub = substructure( name, compound_type, smiles=smiles_string, atoms_to_ignore=to_ignore)
        v = self.structures.create_vertex()
        v.value = sub
        self.structures.add_vertex( v)
    self._analyze_structure_dependencies()

  def _analyze_structure_dependencies( self):
    for v1 in self.structures.vertices:
      for v2 in self.structures.vertices:
        if v1 is not v2:
          sub1 = v1.value
          sub2 = v2.value
          if sub1.structure.contains_substructure( sub2.structure):
            self.structures.add_edge( v2, v1)

  def _find_head_structures( self):
    for v1 in self.structures.vertices:
      v1.properties_['in_links'] = []
      for v2 in self.structures.vertices:
        if v1 in v2.neighbors:
          v1.properties_['in_links'].append( v2.value)
    heads = []
    for v in self.structures.vertices:
      if len( v.properties_['in_links']) == 0:
        # we have got a top of the tree
        heads.append( v)
    return heads
  
  def _compute_search_trees( self):
    heads = self._find_head_structures()
    assert heads
    for head in heads:
      d = self.structures.mark_vertices_with_distance_from( head)
      dv = [(v.properties_['d'],v) for v in self.structures.vertices if 'd' in v.properties_]
      dv.sort( reverse=True)
      for d,v in dv:
        for d,parent in dv:
          if parent.value in v.properties_['in_links']:
            if v.value not in parent.value.children:
              # this ensures only one copy of a substructure in a tree
              parent.value.children.append( v.value)
            break
    return [h.value for h in heads]

  def find_substructures_in_mol( self, mol):
    heads = self._find_head_structures()
    structs = dict.fromkeys( self.structures.vertices, None)
    for head in heads:
      self._find_matching_substructure( head, mol, structs)
    return [k.value for k,v in structs.iteritems() if v==True]

      
  def _find_matching_substructure( self, v, mol, structs):
    if mol.contains_substructure( v.value.structure):
      structs[v] = True
      for n in v.neighbors:
        if structs[n] == None:
          self._find_matching_substructure( n, mol, structs)
          # if the child matches we set v to false - we only want the most
          # closely matching to be caught
          if structs[n] == True:
            structs[v] = False
    else:
      structs[v] = False



class substructure( object):

  def __init__( self, name, compound_type, smiles="", atoms_to_ignore=None):
    self.name = name
    self.compound_type = compound_type
    self.structure = None
    if smiles:
      self.read_smiles( smiles, atoms_to_ignore=atoms_to_ignore)
    self.children = []

  def __str__( self):
    return self.smiles_string

  def read_smiles( self, smiles_string, atoms_to_ignore=None):
    self.atoms_to_ignore = atoms_to_ignore
    self.smiles_string = smiles_string.strip()
    self.structure = smiles.text_to_mol( smiles_string)


  
if __name__ == "__main__":
  ssm = substructure_search_manager()
  ssm.read_structure_file()
  print ssm.structures
  print ssm.structures.is_connected()
  dump = ssm.structures.get_graphviz_text_dump()
  f = file( "dump.dot", "w")
  f.write( dump)
  f.close()

  def print_tree( x, l):
    print l*" ", x #, x.children
    for ch in x.children:
      print_tree( ch, l+2)
  
  #for tree in ssm._compute_search_trees():
  #  print_tree( tree, 0)

  mol = smiles.text_to_mol( "COCCC(=O)OC")
  subs = ssm.find_substructures_in_mol( mol)
  print map( str, subs)
