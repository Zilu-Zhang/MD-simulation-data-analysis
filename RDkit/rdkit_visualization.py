from rdkit import Chem

## partial charge loading
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps

## logP loading
from rdkit.Chem import rdMolDescriptors

## index figure
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

IPythonConsole.drawOptions.addAtomIndices = True
IPythonConsole.molSize = 1000,1000
mol = Chem.MolFromSmiles(
    "CCOC1=NC2=CC=CC(=C2N1CC3=CC=C(C=C3)C4=CC=CC=C4C5=NN=NN5)C(=O)OC(C)OC(=O)OC6CCCCC6"
)

## H bonds analysis
import numpy as np
from rdkit.Chem import Descriptors

## partial charge
mol = Chem.MolFromSmiles('O=C(C1=CC(O)=C(C(O)=C1)O)OC2C(O)C(OC(COC(C3=CC(O)=C(C(O)=C3)O)=O)C2O)OC(C4=CC(O)=C(O)C(O)=C4)=O')
AllChem.ComputeGasteigerCharges(mol)
contribs = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(mol.GetNumAtoms())]
#fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)

## LogP
contribs = rdMolDescriptors._CalcCrippenContribs(mol)
#fig = SimilarityMaps.GetSimilarityMapFromWeights(mol,[x for x,y in contribs], colorMap='jet', contourLines=10)

mol = Chem.MolFromSmiles(
    'CNC(=O)C1=NC=CC(=C1)OC2=CC=C(C=C2)NC(=O)NC3=CC(=C(C=C3)Cl)C(F)(F)F'
)
contribs = np.zeros(mol.GetNumAtoms())
# N
contribs[20] = 9
contribs[17] = 8
contribs[1] = 2
contribs[5] = 0
# O
contribs[10] = 0
contribs[3] = 
contribs[19] = 0
fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)
