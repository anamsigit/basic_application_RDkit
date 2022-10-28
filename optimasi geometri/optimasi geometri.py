from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from rdkit.Chem import Draw
import simpangambarmodul
from rdkit.Chem import PyMol

ibu = Chem.MolFromSmiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O')
ibuH = Chem.AddHs(ibu)
AllChem.EmbedMolecule(ibuH)
v = PyMol.MolViewer()
v.ShowMol(ibuH)
v.server.do('ray')
v.GetPNG()
AllChem.MMFFOptimizeMolecule(ibuH)
v.ShowMol(ibuH,name='optimized',showOnly=True);
v.server.do('ray')
v.GetPNG()
simpangambarmodul.simpangambar(v, "original")

