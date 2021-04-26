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

          
