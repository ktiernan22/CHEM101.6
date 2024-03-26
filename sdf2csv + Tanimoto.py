## Converts a SDF file of small molecules to Smiles codes with their Tanimoto
## Similarity appended on the end
##Reference : https://projects.volkamerlab.org/teachopencadd/talktorials/T004_compound_similarity.html

import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdmolfiles
from rdkit.Chem import SDMolSupplier
from rdkit.Chem.rdmolfiles import SmilesWriter
from rdkit.Chem import (AllChem, rdchem, Draw, Descriptors, MACCSkeys, rdFingerprintGenerator)

#generate smiles file of all molecules from SDF file
i = 0
with Chem.SDMolSupplier('Dartmouth College_send-001.sdf') as sdf:
    for mol in sdf:
        if mol is not None: 
            i += 1
        if i == 2070:
            B5 = mol

#generate a dataframe    
    molecules = pd.DataFrame(sdf)

#calculate and generate molecular weight in dataframe
    MolWt = [Descriptors.MolWt(mol) for mol in sdf]
    Names = [mol.GetProp('IDNUMBER') for mol in sdf]
    smiles = [Chem.MolToSmiles(mol) for mol in sdf]

    molecules['SMILES'] = smiles
    molecules['Name'] = Names
    molecules['molecule_weight'] = MolWt

#add plate and row and column values to dataframe
    molecules['Plate'] = [mol.GetProp('PLATE') for mol in sdf]
    molecules['Column'] = [mol.GetProp('COLUMN') for mol in sdf]
    molecules['Row'] = [mol.GetProp('ROW') for mol in sdf]

#generate MACCS fingerprints
    maccs = [MACCSkeys.GenMACCSKeys(mol) for mol in sdf]

#generate Morgan fingerprints
    fpg = Chem.AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    morgan = [fpg.GetFingerprint(mol) for mol in sdf]

#Copy smiles code into last column to generate molecule drawing in chemdraw
    molecules['structure'] = smiles

#generate Tanimoto similarity comparing MACCS fingerprints
molecule_query_MACCS = maccs[2070]
#molecule_list_MACCS = maccs.to_list()
molecules['tanimoto_maccs'] = DataStructs.BulkTanimotoSimilarity(molecule_query_MACCS, maccs)

#generating Tanimoto similarity comparing Morgan fingerprints
molecule_query_morgan = morgan[2070]
#molecule_list_morgan = morgan.to_list()
molecules['tanimoto_morgan'] = DataStructs.BulkTanimotoSimilarity(molecule_query_morgan, morgan)

#convert data into csv file of SMILES code molecules with info
molecules.to_csv('tanimoto.csv', sep=',', index=False)
