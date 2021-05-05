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

        length_drug = top.residue(0).n_atoms
        length_excp = top.residue(12).n_atoms
        total_length = length_drug * 12 + length_excp * 4
        total = np.zeros(total_length * n_frames)

        ele_mass_drug = []
        ele_mass_excp = []
        ele_mass_part = []
        
        for i in range(length_drug):
            element = str(top.residue(0).atom(i).element.symbol)
            ele_mass_drug.append(mass[element])
        mass_drug = sum(ele_mass_drug)

        for i in range(length_excp):
            element = str(top.residue(12).atom(i).element.symbol)
            ele_mass_excp.append(mass[element])
        mass_excp = sum(ele_mass_excp)
        
        mass_part = mass_drug * 12 + mass_excp * 4

        
        for i in range(n_frames):
            position = np.zeros((total_length + 1, 3))

            if i == 0:
                for j in range(total_length):
                    ele_mass_part.append(mass[str(top.atom(j).element.symbol)])
            
            x = traj.xyz[i, :total_length, 0]
            y = traj.xyz[i, :total_length, 1]
            z = traj.xyz[i, :total_length, 2]

            lst_x, lst_y, lst_z = [], [], []

            for j in range(total_length):
                position[j][:] = x[j], y[j], z[j]
                lst_x.append(x[j] * ele_mass_part[j] / mass_part)
                lst_y.append(y[j] * ele_mass_part[j] / mass_part)
                lst_z.append(z[j] * ele_mass_part[j] / mass_part)

            position[-1][:] = sum(lst_x), sum(lst_y), sum(lst_z)
            
            distance = np.zeros(total_length)
            for j in range(total_length):
                distance[j] = dis(position[-1], position[j])
            
            total[ori:ori + total_length] = distance
            ori += total_length
            
        r_range = np.array([0, 5])
        bin_width = 0.05
        n_bins = int((r_range[1] - r_range[0]) / bin_width)
        g_r, edges = np.histogram(total, range=r_range, bins=n_bins)
        g_r = g_r / (total_length * n_frames)
        r = 0.5 * (edges[1:] + edges[:-1])

        df = pd.DataFrame({'r': r, 'g_r': g_r})
        
        if not os.path.isfile('rdf_atom_mass.xlsx'):
            df.to_excel('rdf_atom_mass.xlsx', '%s' % excipient_name, index = True)
        
        else:
            excel_book = pxl.load_workbook('rdf_atom_mass.xlsx')
            with pd.ExcelWriter('rdf_atom_mass.xlsx', engine = 'openpyxl') as writer:
                writer.book = excel_book
                writer.sheets = {worksheet.title: worksheet for worksheet in excel_book.worksheets}
                df.to_excel(writer, '%s' % excipient_name, index = True)
                writer.save()