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
        total = np.zeros(3200)
        
        for i in range(n_frames):
            start = 0
            position = np.zeros((17,3))
            for j in range(16):
                res = top.residue(j)
                length = res.n_atoms
                start += length
                x = mean(traj.xyz[i, start:start + length - 1, 0])
                y = mean(traj.xyz[i, start:start + length - 1, 1])
                z = mean(traj.xyz[i, start:start + length - 1, 2])
                position[j][:] = x, y, z
            position[-1][:] = mean(array[:-1][0]), mean(array[:-1][1]), mean(array[:-1][2])
            
            distance = np.zeros(16)
            for i in range(16):
                distance[i] = dis(position[-1], position[i])
            
            total[ori:ori + 16] = distance
            ori += 17
            
        r_range = np.array([0.0, 2.5])
        bin_width = 0.005
        n_bins = int((r_range[1] - r_range[0]) / bin_width)
        g_r, edges = np.histogram(total, range=r_range, bins=n_bins)
        g_r = g_r / 3200
        r = 0.5 * (edges[1:] + edges[:-1])

        df = pd.DataFrame({'r': r, 'g_r': g_r}
        
        if not os.path.isfile('rdf.xlsx'):
            df.to_excel('rdf.xlsx', '%s' % excipient_name, index = True)
        
        else:
            excel_book = pxl.load_workbook('rdf.xlsx')
            with pd.ExcelWriter('rdf.xlsx', engine = 'openpyxl') as writer:
                writer.book = excel_book
                writer.sheets = {worksheet.title: worksheet for worksheet in excel_book.worksheets}
                df.to_excel(writer, '%s' % excipient_name, index = True)
                writer.save()
          
