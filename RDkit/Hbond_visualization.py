### create atom-index map as reference (a nested dictionary)
import pandas as pd
import numpy as np
import math
import os
from rdkit import Chem
from rdkit.Chem.Draw import SimilarityMaps
import matplotlib.pyplot as plt

maps = pd.read_excel('map.xlsx', sheet_name = None)
ref = {}
for key, value in maps.items():
    sub = {}
    for i in range(len(value)):
        atom = value['atom'][i]
        index = value['index'][i]
        sub[atom] = index
    ref[key] = sub


### smiles
smile = {}
smile['sorafenib'] = 'CNC(=O)C1=NC=CC(=C1)OC2=CC=C(C=C2)NC(=O)NC3=CC(=C(C=C3)Cl)C(F)(F)F'
smile['budesonide'] = 'CCCC1O[C@@H]2C[C@H]3[C@@H]4CCC5=CC(=O)C=C[C@@]5([C@H]4[C@H](C[C@@]3([C@@]2(O1)C(=O)CO)C)O)C'
smile['candesartan_cilexetil'] = 'CCOC1=NC2=CC=CC(=C2N1CC3=CC=C(C=C3)C4=CC=CC=C4C5=NN=NN5)C(=O)OC(C)OC(=O)OC6CCCCC6'
smile['cholic_acid'] = 'C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C'
smile['glycocholic_acid'] = 'C[C@H](CCC(=O)NCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C'
smile['glycyrrhizin'] = 'C[C@]12CC[C@](C[C@H]1C3=CC(=O)[C@@H]4[C@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC2)C)C)(C)C)O[C@@H]6[C@@H]([C@H]([C@@H]([C@H](O6)C(=O)O)O)O)O[C@H]7[C@@H]([C@H]([C@@H]([C@H](O7)C(=O)O)O)O)O)C)(C)C(=O)O'
smile['hydrocortisone'] = 'C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@]4(C(=O)CO)O)C)O'
smile['indomethacin'] = 'CC1=C(C2=C(N1C(=O)C3=CC=C(C=C3)Cl)C=CC(=C2)OC)CC(=O)O'
smile['meloxicam'] = 'CC1=CN=C(S1)NC(=O)C2=C(C3=CC=CC=C3S(=O)(=O)N2C)O'
smile['neohesperidin_DC'] = 'C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=CC(=C(C(=C3)O)C(=O)CCC4=CC(=C(C=C4)OC)O)O)CO)O)O)O)O)O'
smile['prednisone'] = 'C[C@]12CC(=O)[C@H]3[C@H]([C@@H]1CC[C@@]2(C(=O)CO)O)CCC4=CC(=O)C=C[C@]34C'
smile['ursodiol'] = 'C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2[C@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C'


### correspond Hbond results to the previously created map lists
Hbonds = pd.read_excel('Hbonds_individual.xlsx', sheet_name = None)
for key, value in Hbonds.items():
    path = './%s/' % key
    if os.path.isdir(path) is False:
        os.mkdir(path)
    match_d = {}
    match_a = {}
    once_more = True
    print(key)
    for i in range(len(value)):
        atom = value['residue'][i][value['residue'][i].find('-')+1:]
        resi_cur = value['residue'][i][:value['residue'][i].find('-')]
        resi_nxt = value['residue'][i+1][:value['residue'][i+1].find('-')] if i < len(value)-1 else 0
        molecule = 'sorafenib' if 'ZPE' in resi_cur else key
        index = ref[molecule][atom]
        if once_more:
            if resi_cur != resi_nxt or resi_nxt == 0:
                once_more = False
            match_d[index] = value['Hdonors'][i] if not math.isnan(value['Hdonors'][i]) else 0
            match_a[index] = value['Hacceptors'][i] if not math.isnan(value['Hacceptors'][i]) else 0
            if once_more == True:
                continue
            # real code
            mol = Chem.MolFromSmiles(smile[molecule])
            contribs = np.zeros(mol.GetNumAtoms())
            for index, density in match_d.items():
                contribs[index] = int(density)
        
            fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)
            fig.savefig('./%s/donor_%s.png' % (key, resi_cur), bbox_inches = 'tight')
            plt.close()

            contribs = np.zeros(mol.GetNumAtoms())
            for index, density in match_a.items():
                contribs[index] = - int(density)

            fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)
            fig.savefig('./%s/acceptor_%s.png' % (key, resi_cur), bbox_inches = 'tight')
            plt.close()

            match_d = {}
            match_a = {}
            once_more = True


