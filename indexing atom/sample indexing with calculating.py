import simpangambarmodul
from types import SimpleNamespace
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
IPythonConsole.drawOptions.addAtomIndices = True
IPythonConsole.molSize = 300,300


m = Chem.MolFromSmiles('c1ncncc1C(=O)[O-]')
AllChem.ComputeGasteigerCharges(m)
simpangambarmodul.simpangambar(m, "")

# with calculation
m2 = Chem.Mol(m)
for at in m2.GetAtoms():
    lbl = '%.2f'%(at.GetDoubleProp("_GasteigerCharge"))
    at.SetProp('atomNote',lbl)
simpangambarmodul.simpangambar(m2,"calculation")