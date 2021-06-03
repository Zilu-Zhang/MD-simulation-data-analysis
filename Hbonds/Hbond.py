### This will generate a Hbonds.xlsx file containing information about total Hbonds, interactive Hbonds and their ratios of each excipient
import os
import os.path
import mdtraj as md
import numpy as np
import pandas as pd
import openpyxl as pxl

def label(hbond):
    r1 = traj.topology.atom(hbond[0])
    r2 = traj.topology.atom(hbond[2])
    n = 0
    if str(r1)[:3] != str(r2)[:3]:
        n = 1
    return n

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
                total_Hbonds[i] = len(hbonds)
                interactive_Hbonds[i] = number
                ratio[i] = interactive_Hbonds[i] / total_Hbonds[i]

        df = pd.DataFrame({'total_Hbonds': total_Hbonds, 'interactive_Hbonds': interactive_Hbonds, 'ratio': ratio})
		
        if not os.path.isfile('Hbonds.xlsx'):
            df.to_excel('Hbonds.xlsx', '%s' % excipient_name, index = True)

        else:
            excel_book = pxl.load_workbook('Hbonds.xlsx')
            with pd.ExcelWriter('Hbonds.xlsx', engine = 'openpyxl') as writer:
                writer.book = excel_book
                writer.sheets = {worksheet.title: worksheet for worksheet in excel_book.worksheets}
                df.to_excel(writer, '%s' % excipient_name, index = True)
                writer.save()
