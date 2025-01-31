import logging

from src.predict_pka import *


class CTSMolgpka:

	def __init__(self):
		self.pka_dec= 2

	def convert_floats(self, pka_list):
		"""
		Python's json serializer doesn't like np.float32 types, so
		converting them into python floats.
		"""
		return [round(float(i), self.pka_dec) for i in pka_list]

	def convert_to_json_serializable(self, obj):
	    if isinstance(obj, dict):
	        # Convert dictionary with potential tuple keys to string keys
	        return {str(k): self.convert_to_json_serializable(v) for k, v in obj.items()}
	    elif isinstance(obj, list):
	        return [self.convert_to_json_serializable(item) for item in obj]
	    elif isinstance(obj, tuple):
	        return str(obj)
	    elif isinstance(obj, np.floating):
	        return float(obj)
	    elif isinstance(obj, np.integer):
	        return int(obj)
	    elif isinstance(obj, (str, int, float, bool, type(None))):
	        return obj
	    return str(obj)

	def run_molgpka(self, smiles):
		mol = Chem.MolFromSmiles(smiles)
		molgpka_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol))) #ammended to return smile string to be used as input for chem axon
		base_dict, acid_dict, new_idx, mol = predict(mol)
		atom_idx = list(base_dict.keys()) + list(new_idx)
		pkas = list(base_dict.values()) + list(acid_dict.values())
		sites = len(pkas)
		yield sites, pkas, atom_idx, molgpka_smiles, acid_dict, base_dict, new_idx
		
	def main(self, smiles):
		
		"""
		Main function for returning pkas and/or microspecies.
		Examples=['CC(O)=O','CC(C)C(N)C(O)=O','C(O)1=CC=C(N)C=C1','NC(CCS)C(O)=O','NC(CC1=CN=CN1)C(O)=O']
		"""
		pka_sites, pka_list = None, None

		data = self.run_molgpka(smiles)
		
		for n,p,idx,smiles,a,b,new_idx in data:

			new_acid_dict = dict(zip(new_idx, a.values()))

			mg_dict = {}
			mg_tuples = []

			for k,v in new_acid_dict.items():
				cat='acid'
				atom=k
				pka=round(v,2)
				mg_tuples.append((cat,atom))
				mg_dict.update({(cat,atom):pka})
			
			for k,v in b.items():
				cat='base'
				atom=k
				pka=round(v,2)
				mg_tuples.append((cat,atom))
				mg_dict.update({(cat,atom):pka})

			pka_sites = n
			pka_list = p
			molgpka_smiles = smiles
			molgpka_index = idx

		pka_list = self.convert_floats(pka_list)

		molgpka_dict = {}
		for key, value in zip(molgpka_index, pka_list):
			if key in molgpka_dict:
				molgpka_dict[key].append(value)
			else:
				molgpka_dict[key] = [value]

		for key, val in molgpka_dict.items():
			molgpka_dict[key] = ', '.join(map(str, val))

		# TODO: Get rid of molgpka_dict once updates are working.

		return smiles, pka_sites, pka_list, molgpka_smiles, molgpka_dict, molgpka_index, mg_dict, mg_tuples



if __name__ == "__main__":
	# #expected output
	# SMILES: CC(O)=O 
	# 	# of sites: 1 
	# 	pKa: [8.337572]
	# SMILES: CC(C)C(N)C(O)=O 
	# 	# of sites: 2 
	# 	pKa: [9.808557, 8.171646]
	# SMILES: C(O)1=CC=C(N)C=C1 
	# 	# of sites: 3 
	# 	pKa: [5.06671, 14.116651, 12.194197]
	# SMILES: NC(CCS)C(O)=O 
	# 	# of sites: 3 
	# 	pKa: [8.882975, 10.914888, 7.571865]
	# SMILES: NC(CC1=CN=CN1)C(O)=O 
	# 	# of sites: 4 
	# 	pKa: [8.179463, 6.96088, 7.530961, 14.086538]

	examples = ['CC(O)=O','CC(C)C(N)C(O)=O','C(O)1=CC=C(N)C=C1','NC(CCS)C(O)=O','NC(CC1=CN=CN1)C(O)=O']

	molgpka = CTSMolgpka()

	for smi in examples:
		data = molgpka.main(smi)
