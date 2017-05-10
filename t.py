#test lets grab some chemspider object. we can test common_name ==> CAS and CAS ==> SMILES
#also test obtaining density from CAS

from rdkit import Chem
from chemspipy import ChemSpider
cs = ChemSpider("ccab768d-273b-41cb-a63e-2b8ebf2508ad")

r = raw_input ("name >")

for result in cs.simple_search(r):
	info = cs.get_extended_compound_info(result.csid)
	for key in info:
		print key, info[key]
	#mol = Chem.MolFromSmiles(info["smiles"])
	#print mol.GetNumHeavyAtoms()
	#print info["smiles"], 
	print "\n"