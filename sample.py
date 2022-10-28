from rdkit import Chem
from rdkit.Chem import Draw

m = Chem.MolFromSmiles('Cc1ccccc1')

print(Chem.MolToMolBlock(m))    

smiles = 'C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O'
n = Chem.MolFromSmiles(smiles)
p = Draw.MolToImage(n)
p.save('D:/dokumen/Academic/Pemograman/python/Library RDkit/rdkit output.png','PNG')
print(type(p))