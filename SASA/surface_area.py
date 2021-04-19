import os
import os.path
import mdtraj as md
import numpy as np
import pandas as pd
import openpyxl as pxl
import statistics

n_frames = 200
for filename in os.listdir('./'):
    if filename.endswith('.pdb'):
        excipient_name = filename[17:-4]
        total_Hbonds = np.zeros(n_frames)
        interactive_Hbonds = np.zeros(n_frames)
        ratio = np.zeros(n_frames)

        for i in range(n_frames):
            traj = md.load_frame(filename, i)
            hbonds = md.baker_hubbard(traj)
            number = 0
            for hbond in hbonds:
                n = label(hbond)
                number += n
