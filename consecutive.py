import os
import os.path
import mdtraj as md
import numpy as np
import pandas as pd
import openpyxl as pxl
import statistics

n_frames = 20000 - 1
for filename in os.listdir('./'):
    if filename.endswith('.pdb'):
        excipient_name = filename[17:-4]
        RMSD_1 = np.zeros(n_frames)
        RMSD_2 = np.zeros(n_frames)
        AVERAGE = np.zeros(n_frames)

        traj_1 = md.load(filename, atom_indices = np.arange(0,48))
        traj_2 = md.load(filename, atom_indices = np.arange(48,96))

        for i in range(n_frames):
            if i == 0:
                frame_residue1_ref = traj_1.slice(i)
                frame_residue2_ref = traj_2.slice(i)

            else:
                frame_residue1_ref = frame_residue1_tag
                frame_residue2_ref = frame_residue2_tag

            frame_residue1_tag = traj_1.slice(i+1)
            frame_residue2_tag = traj_2.slice(i+1)

            RMSD_1[i] = md.rmsd(frame_residue1_tag, frame_residue1_ref)
            RMSD_2[i] = md.rmsd(frame_residue2_tag, frame_residue2_ref)
            AVERAGE[i] = (RMSD_1[i] + RMSD_2[i])/2

        mean = statistics.mean(AVERAGE)
        sd = statistics.stdev(AVERAGE)
        df = pd.DataFrame({'RMSD_1': RMSD_1, 'RMSD_2': RMSD_2, 'AVERAGE': AVERAGE, 'Mean': mean, 'SD': sd})

        if not os.path.isfile('convergence_results.xlsx'):
            df.to_excel('consecutive_results.xlsx', '%s' % excipient_name, index = False)

        else:
            excel_book = pxl.load_workbook('consecutive_results.xlsx')
            with pd.ExcelWriter('consecutive_results.xlsx', engine = 'openpyxl') as writer:
                writer.book = excel_book
                writer.sheets = {worksheet.title: worksheet for worksheet in excel_book.worksheets}
                df.to_excel(writer, '%s' % excipient_name, index = False)
                writer.save()
