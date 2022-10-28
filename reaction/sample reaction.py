import simpangambarmodul
from rdkit import Chem
from rdkit.Chem import rdChemReactions

rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]') #formatnya adalah reaction ::= reactants ">>" products
reacts = (Chem.MolFromSmiles('[Li+]'),Chem.MolFromSmiles('O'))
products = rxn.RunReactants(reacts)

sumproduk = len(products)
print(sumproduk)
# hasil = Chem.MolToSmiles(products[0][0])
# print(hasil)
# simpangambarmodul.simpangambar(products[0][0], "produk pertama reaksi")