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
for filename in os.listdir('./'):
    if filename.endswith('.pdb'):
        excipient_name = filename[17:-4]
        traj = md.load(filename)
        top = traj.topology
        ori = 0
        res_d = top.residue(0)
        res_e = top.residue(12)
        d_length = res_d.n_atoms
        e_length = res_e.n_atoms
        total_length = d_length * 12 + e_length * 4
        total = np.zeros(total_length * n_frames)
        
        for i in range(n_frames):
            position = np.zeros((total_length + 1, 3))
            
            a = traj.xyz[i, :total_length, 0]
            b = traj.xyz[i, :total_length, 1]
            c = traj.xyz[i, :total_length, 2]
            position[:total_length][:] = a, b, c
            
            x = mean(traj.xyz[i, :total_length, 0])
            y = mean(traj.xyz[i, :total_length, 1])
            z = mean(traj.xyz[i, :total_length, 2])
            position[-1][:] = x, y, z
            
            distance = np.zeros(total_length)
            for j in range(total_length):
                distance[j] = dis(position[-1], position[j])
            
            total[ori:ori + total_length] = distance
            ori += total_length
            
        r_range = np.array([0, 5])
        bin_width = 0.005
        n_bins = int((r_range[1] - r_range[0]) / bin_width)
        g_r, edges = np.histogram(total, range=r_range, bins=n_bins)
        g_r = g_r / (total_length * n_frames)
        r = 0.5 * (edges[1:] + edges[:-1])

        df = pd.DataFrame({'r': r, 'g_r': g_r})
        
        if not os.path.isfile('rdf_atom.xlsx'):
            df.to_excel('rdf_atom.xlsx', '%s' % excipient_name, index = True)
        
        else:
            excel_book = pxl.load_workbook('rdf_atom.xlsx')
            with pd.ExcelWriter('rdf_atom.xlsx', engine = 'openpyxl') as writer:
                writer.book = excel_book
                writer.sheets = {worksheet.title: worksheet for worksheet in excel_book.worksheets}
                df.to_excel(writer, '%s' % excipient_name, index = True)
                writer.save()
