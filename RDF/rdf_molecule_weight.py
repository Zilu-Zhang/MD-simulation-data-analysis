import mdtraj as md
import numpy as np  
import os
import os.path
import pandas as pd
import openpyxl as pxl
from statistics import mean
from math import sqrt

def dis(ref, tag):
    x = ref[0] - tag[0]
    y = ref[1] - tag[1]
    z = ref[2] - tag[2]
    return sqrt(x**2 + y**2 + z**2)

n_frames = 200
mass = {
    'C': 12,
    'H': 1,
    'O': 16,
    'N': 14,
    'S': 32,
    'F': 19,
    'Cl': 35.45
}
for filename in os.listdir('./'):
    if filename.endswith('.pdb'):
        excipient_name = filename[17:-4]
        traj = md.load(filename)
        top = traj.topology
        ori = 0
        total = np.zeros(16 * n_frames)
        
        length_drug = top.residue(0).n_atoms
        length_excp = top.residue(12).n_atoms
        
        element_drug = []
        element_excp = []
        
        for i in range(length_drug):
            element_mass = str(top.residue(0).atom(i).name)
            element_drug.append(mass[element_mass])
        mass_drug = sum(element_drug)
        
        for i in range(length_excp):
            element_mass_ex = str(top.residue(12).atom(i).name)
            element_excp.append(mass[element_mass_ex])
        mass_excp = sum(element_excp)
        
        mass_part = mass_drug + mass_excp
        
        for i in range(n_frames):
            start = 0
            position = np.zeros((17,3))
            for j in range(16):
                res = top.residue(j)
                length = res.n_atoms
                
                if j <= 11:
                    x = sum(list(traj.xyz[i, start:start + length, 0])) * element_drug) / mass_drug
                    y = sum(list(traj.xyz[i, start:start + length, 1])) * element_drug) / mass_drug
                    z = sum(list(traj.xyz[i, start:start + length, 2])) * element_drug) / mass_drug
                else:
                    x = sum(list(traj.xyz[i, start:start + length, 0])) * element_excp) / mass_excp
                    y = sum(list(traj.xyz[i, start:start + length, 1])) * element_excp) / mass_excp
                    z = sum(list(traj.xyz[i, start:start + length, 2])) * element_excp) / mass_excp
                
                position[j][:] = x, y, z
                start += length

            position[-1][:] = mean(position[:-1][0]), mean(position[:-1][1]), mean(position[:-1][2])
            
            distance = np.zeros(16)
            for h in range(16):
                distance[h] = dis(position[-1], position[h])
            
            total[ori:ori + 16] = distance
            ori += 16
            
        r_range = np.array([0, 5])
        bin_width = 0.005
        n_bins = int((r_range[1] - r_range[0]) / bin_width)
        g_r, edges = np.histogram(total, range=r_range, bins=n_bins)
        g_r = g_r / (16 * n_frames)
        r = 0.5 * (edges[1:] + edges[:-1])

        df = pd.DataFrame({'r': r, 'g_r': g_r})
        
        if not os.path.isfile('rdf_molecule.xlsx'):
            df.to_excel('rdf_molecule.xlsx', '%s' % excipient_name, index = True)
        
        else:
            excel_book = pxl.load_workbook('rdf_molecule.xlsx')
            with pd.ExcelWriter('rdf_molecule.xlsx', engine = 'openpyxl') as writer:
                writer.book = excel_book
                writer.sheets = {worksheet.title: worksheet for worksheet in excel_book.worksheets}
                df.to_excel(writer, '%s' % excipient_name, index = True)
                writer.save()
