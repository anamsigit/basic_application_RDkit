from rdkit import Chem
resulting= Chem.MolFromSmiles('CCO').HasSubstructMatch(Chem.MolFromSmiles('CCO'))
print(resulting)