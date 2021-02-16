import os
import os.path
import mdtraj as md
import numpy as np
import pandas as pd
import openpyxl as pxl
import statistics

n_frames = 20000
for filename in os.listdir('./'):
    if filename.endswith('.pdb'):
        excipient_name = filename[17:-4]
        RMSD_1 = np.zeros(n_frames)
        RMSD_2 = np.zeros(n_frames)
        AVERAGE = np.zeros(n_frames)

        traj_1 = md.load(filename, atom_indices = np.arange(0,48))
        traj_2 = md.load(filename, atom_indices = np.arange(48,96))

        stable_traj_1 = traj_1.slice(500)
        stable_traj_2 = traj_2.slice(500)

        RMSD_1 = md.rmsd(traj_1, stable_traj_1)
        RMSD_2 = md.rmsd(traj_2, stable_traj_2)

        for i in range(n_frames):
            AVERAGE[i] = (RMSD_1[i] + RMSD_2[i])/2

        mean = statistics.mean(AVERAGE)
        sd = statistics.stdev(AVERAGE)
        df = pd.DataFrame({'RMSD_1': RMSD_1, 'RMSD_2': RMSD_2, 'AVERAGE': AVERAGE, 'Mean': mean, 'SD': sd})

        if not os.path.isfile('stable_results.xlsx'):
            df.to_excel('stable_results.xlsx', '%s' % excipient_name, index = False)

        else:
            excel_book = pxl.load_workbook('stable_results.xlsx')
            with pd.ExcelWriter('stable_results.xlsx', engine = 'openpyxl') as writer:
                writer.book = excel_book
                writer.sheets = {worksheet.title: worksheet for worksheet in excel_book.worksheets}
                df.to_excel(writer, '%s' % excipient_name, index = False)
                writer.save()
