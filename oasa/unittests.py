#--------------------------------------------------------------------------
#     This file is part of OASA - a free chemical python library
#     Copyright (C) 2007 Beda Kosata <beda@zirael.org>

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
    
  def setUp(self):
    pass
    
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
              ]
    
  def setUp(self):
    pass
    
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
              ]
    
  def setUp(self):
    pass
    
  def _testformula(self, num):
    smile1, smile2, result = self.formulas[num]
    m1 = smiles.text_to_mol( smile1)
    m2 = smiles.text_to_mol( smile2)
    self.assert_(molecule.equals(m1,m2,level=3)==result)

# this creates individual test for substructures
for i in range( len( TestEqualSMILES.formulas)):
  setattr( TestEqualSMILES, "testformula"+str(i+1), create_test(i,"_testformula"))

## // substructure testing




if __name__ == '__main__':
  unittest.main()
