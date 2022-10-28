from http.client import ImproperConnectionState
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
IPythonConsole.ipython_useSVG=False  #< set this to False if you want PNGs instead of SVGs
import simpangambarmodul

def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

# Test in a kinase inhibitor
mol = Chem.MolFromSmiles("C1CC2=C3C(=CC=C2)C(=CN3C1)[C@H]4[C@@H](C(=O)NC4=O)C5=CNC6=CC=CC=C65")
simpangambarmodul.simpangambar(mol, "original")

mol_with_atom_index(mol)
simpangambarmodul.simpangambar(mol_with_atom_index(mol), "denan indexing")