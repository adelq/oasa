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
from sets import Set

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

  def which_substructure_is_more_specific( self, s1, s2):
    """returns 1 if s1 is more specific,
               2 if s2 is more specific,
               3 if both are the same,
               4 if the relationship cannot be established
               (there is not connection between s1 and s2 in the
               substructure interdependency graph - this signifies
               a problem, probably in the matching algorithm
               """
    if s1 is s2:
      return 3
    # just find vertices corresponding to s1 and s2
    v1 = None
    v2 = None
    for v in self.structures.vertices:
      if v.value == s1:
        v1 = v
      elif v.value == s2:
        v2 = v
    assert v1
    assert v2
    p1 = self.structures.path_exists( v1, v2)
    p2 = self.structures.path_exists( v2, v1)
    print p1, p2, v1.value, v2.value
    if (p1 and p2) or (not p1 and not p2):
      # both paths exist or none does - this should not happen
      return 4
    if p1:
      return 2
    if p2:
      return 1

  def find_substructures_in_mol( self, mol):
    hits = []
    for v in self.structures.vertices:
      struct = v.value
      ms = struct.find_matches( mol)
      hits += ms
    hit_num = len( hits)
    to_delete = True
    while to_delete:
      to_delete = []
      for hit1 in hits:
        for hit2 in hits:
          if hit1 is not hit2:
            if Set(hit1.get_significant_atoms()) & Set(hit2.get_significant_atoms()):
              winner = self.which_substructure_is_more_specific( hit1.substructure, hit2.substructure)
              if winner == 1:
                to_delete.append( hit2)
              elif winner == 2:
                to_delete.append( hit1)
              elif winner == 3:
                pass # we preserve both hits
              else:
                raise ValueError( "Relationship between competing fragments could not be established,\nthere is probaly and error in the substructure matching code.")
        if to_delete:
          break
      hits = [hit for hit in hits if not hit in to_delete]
    return hits


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
    self.smiles_string = smiles_string.strip()
    self.structure = smiles.text_to_mol( smiles_string)
    if atoms_to_ignore:
      self.atoms_to_ignore = [self.structure.vertices[x-1] for x in atoms_to_ignore]
    else:
      atoms_to_ignore = []

  def find_matches( self, mol):
    ret = []
    ms = list( mol.select_matching_substructures( self.structure, implicit_freesites=True, auto_cleanup=False))
    for atoms in ms:
      num = min( atoms[0].properties_['subsearch'].keys())
      atoms_in_fragment = [a.properties_['subsearch'][num] for a in atoms]
      ret.append( substructure_match( atoms, atoms_in_fragment, self))
    mol.clean_after_search( self.structure)
    return ret
  



class substructure_match( object):

  def __init__( self, atoms_found, atoms_searched, substruct):
    """atoms_found are atoms in the molecule we searched in,
    atoms_searched are atoms in the fragment we used for search,
     (atoms_searched and atoms_found are guaranteed to have the same order, thus
      allowing matching between the two structures),
    substruct is the substructure instance"""
    self.substructure = substruct
    self.atoms_found = atoms_found
    self.atoms_searched = atoms_searched

  def __str__( self):
    return "<Match of %s with %d atoms, %d significant atoms>" % (self.substructure, len( self.atoms_found), len( self.get_significant_atoms()))

  def get_significant_atoms( self):
    ret = []
    for i,af in enumerate( self.atoms_found):
      as = self.atoms_searched[i]
      if as not in self.substructure.atoms_to_ignore:
        ret.append( af)
    return ret

  
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

  mol = smiles.text_to_mol( "COCC(=O)OC")
  subs = ssm.find_substructures_in_mol( mol)
  print map( str, subs)
