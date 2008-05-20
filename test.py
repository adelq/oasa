import oasa
#print oasa.CAIRO_AVAILABLE

def cairo_out_test2():
    #mol = smiles.text_to_mol( "COC(=O)CNC(C1=CC=CC=C1)C2=C(C=CC(=C2)Br)NC(=O)C3=CC(=CC=C3)Cl")
    #mol = smiles.text_to_mol( "c1nccc2c1cncn2")
    mol = oasa.smiles.text_to_mol( "OC(O)C(O)C(O)CO")
    mol.normalize_bond_length( 30)
    mol.remove_unimportant_hydrogens()
    c = oasa.cairo_out.cairo_out( color_bonds=True, color_atoms=True)
    c.show_hydrogens_on_hetero = True
    c.font_size = 20
    c.mol_to_cairo( mol, "test.png")

cairo_out_test2()
