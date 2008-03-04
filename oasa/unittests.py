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

import unittest
import linear_formula
import smiles
import molecule

# helper functions

def create_test(i,name):
  def test(obj):
    getattr(obj, name)(i)
  return test


## linear formula testing

class TestLinearFormula(unittest.TestCase):

  formulas = [("CH3COOH","CC(=O)O",0,0),
              ("CH3C(CH3)3","CC(C)(C)C",0,0),
              ("CH2C(CH3)3","CC(C)(C)C",1,0),
              ("(CH2)7Cl","CCCCCCCCl",1,0),
              ]
    
  def _testformula(self, num):
    l = linear_formula.linear_formula()
    linear, smile, start_valency, end_valency = self.formulas[num]
    m1 = l.parse_text( linear, start_valency=start_valency, end_valency=end_valency) 
    m2 = smiles.text_to_mol( smile)
    self.assert_(molecule.equals(m1,m2,level=3))

# this creates individual test for formulas
for i in range( len( TestLinearFormula.formulas)):
  setattr( TestLinearFormula, "testformula"+str(i+1), create_test(i,"_testformula"))

## // linear formula testing



## substructure testing

class TestSubstructure(unittest.TestCase):

  formulas = [("CC(=O)C","CC(=O)C",True),
              ("CCCC","CCC",True),
              ("CC(=O)C","C=O",True),
              ("C=O","CC(=O)C",False),
              ("C(=O)H","C=O",True),
              ("C=O","C(=O)H",True),      # implicit hydrogens work
              ("C(=O)O","C(=O)OH",True),  # implicit hydrogens work
              ("C(=O)OH","C(=O)O",True),
              ("C(=O)OC","C(=O)OH",False),  # explicit hydrogens work
              ("C(=O)OC","C(=O)[OH]",False),  # explicit hydrogens in atom specs work
              ("C(=O)OC","C(=O)O",True),
              #
              ("CC(=O)H","C(=O)H",True),
              ]
    
  def _testformula(self, num):
    smile1, smile2, result = self.formulas[num]
    m1 = smiles.text_to_mol( smile1)
    m2 = smiles.text_to_mol( smile2)
    self.assert_(m1.contains_substructure(m2)==result)

# this creates individual test for substructures
for i in range( len( TestSubstructure.formulas)):
  setattr( TestSubstructure, "testformula"+str(i+1), create_test(i,"_testformula"))

## // substructure testing


## SMILES equality testing

class TestEqualSMILES(unittest.TestCase):

  formulas = [("Sc1ccccc1","S-c1ccccc1",True),  # check Sc (in PT scandium) bug
              ("Oc1ccccc1","O-c1ccccc1",True),
              ("c1ccccc1","C:1:C:C:C:C:C:1", True),
              ("c1ccccc1","[CH]:1:[CH]:[CH]:[CH]:[CH]:[CH]:1", True),
              ("c1cscc1","C1=C-S-C=C1", True),
              ]
    
  def _testformula(self, num):
    smile1, smile2, result = self.formulas[num]
    m1 = smiles.text_to_mol( smile1)
    m2 = smiles.text_to_mol( smile2)
    self.assert_(molecule.equals(m1,m2,level=3)==result)

# this creates individual test for substructures
for i in range( len( TestEqualSMILES.formulas)):
  setattr( TestEqualSMILES, "testformula"+str(i+1), create_test(i,"_testformula"))

## // SMILES equality testing


## SMILES reading testing

class TestSMILESReading(unittest.TestCase):

  formulas = [("Sc1ccccc1",("C6H6S",)),
              ("Oc1ccccc1",("C6H6O",)),
              ("[Na+].[Cl-]", ("Na","Cl")),
              ("[O-]c1ccccc1.[Na+]",("C6H5O","Na")),
              ("O=C[O-].[NH4+]",("CHO2","H4N")),
              ("c1ccccc1-c1ccccc1",("C12H10",)),
              ("c1cscc1",("C4H4S",)),
              ]
    
  def _testformula(self, num):
    smile1, sum_forms = self.formulas[num]
    conv = smiles.converter()
    mols = conv.read_text( smile1)
    for i,mol in enumerate( mols):
      assert i < len( sum_forms)
      self.assertEqual( str( mol.get_formula_dict()), sum_forms[i]) 

  def test_empty_smiles( self):
    conv = smiles.converter()
    for text in ("", "  "):
      mols = conv.read_text( "")
      self.assertEqual( mols, [])


# this creates individual test for substructures
for i in range( len( TestSMILESReading.formulas)):
  setattr( TestSMILESReading, "testformula"+str(i+1), create_test(i,"_testformula"))

## // SMILES equality testing


## SMILES Reaction support

class TestSMILESReactionSupport(unittest.TestCase):
  
  def test1(self):
    """tests handling of reactions by the SMILES reader on a preparation of methyl-formate"""
    c = smiles.converter()
    reacts = c.read_text( "O=CO.CO>[H+]>O=COC.O")
    self.assertEqual( reacts, c.result)
    self.assertEqual( c.last_status, c.STATUS_OK)
    self.assertEqual( len( reacts), 1)
    react = reacts[0]
    self.assertEqual( len( react.reactants), 2)
    self.assertEqual( len( react.reactants[0].molecule.atoms), 3)
    self.assertEqual( len( react.reactants[1].molecule.atoms), 2)
    self.assertEqual( len( react.reagents[0].molecule.atoms), 1)
    self.assertEqual( len( react.products[0].molecule.atoms), 4)
    self.assertEqual( len( react.products[1].molecule.atoms), 1)
    self.assertEqual( str( react.products[0].molecule.get_formula_dict()), "C2H4O2")

  def test2(self):
    """test reactions with some empty parts"""
    c = smiles.converter()
    reacts = c.read_text( "C=C.[H][H]>>CC")
    self.assertEqual( len( reacts), 1)
    react = reacts[0]
    self.assertEqual( len( react.reagents), 0)
    self.assertEqual( len( react.reactants), 2)
    self.assertEqual( len( react.products), 1)




## // SMILES Reaction support

## Reaction test

import reaction

class TestReactionComponent(unittest.TestCase):

  def test1(self):
    mol = smiles.text_to_mol( "CCCO")
    rc = reaction.reaction_component( mol, 2)
    self.assertEqual( rc.stoichiometry, 2)
    self.assertRaises( Exception, reaction.reaction_component, mol, "x")
    self.assertRaises( Exception, rc._set_stoichiometry, "x")
    self.assertRaises( Exception, reaction.reaction_component, 2, 2)
    self.assertRaises( Exception, rc._set_molecule, "x")    



if __name__ == '__main__':
  unittest.main()
