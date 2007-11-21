import unittest
import linear_formula
import smiles
import molecule

class TestLinearFormula(unittest.TestCase):

  formulas = [("CH3COOH","CC(=O)O",0,0),
              ("CH3C(CH3)3","CC(C)(C)C",0,0),
              ("CH2C(CH3)3","CC(C)(C)C",1,0),
              ("(CH2)7Cl","CCCCCCCCl",1,0),
              ]
    
  def setUp(self):
    self.seq = range(10)
    
  def _testformula(self, num):
    l = linear_formula.linear_formula()
    linear, smile, start_valency, end_valency = self.formulas[num]
    m1 = l.parse_text( linear, start_valency=start_valency, end_valency=end_valency) 
    m2 = smiles.text_to_mol( smile)
    self.assert_(molecule.equals(m1,m2,level=3))


def create_test(i):
  def test(obj):
    obj._testformula(i)
  return test

# this creates individual test for formulas
for i in range( len( TestLinearFormula.formulas)):
  setattr( TestLinearFormula, "testformula"+str(i+1), create_test(i))
    
if __name__ == '__main__':
  unittest.main()
