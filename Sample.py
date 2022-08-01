#representasi dari smiles
from RDkit import Chem
m = Chem.MolFromSmiler(xxx=xx-x=xx)
m

#representasi dari mol file
#mol file dapat diperoleh dari drugbase
import requests
morphine_url = 'https://'
morphine_mol = requests.get(morphine_url).text

morphine = Chem.MolFromMolBlock(morphine_mol)
morphine
