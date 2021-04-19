import mdtraj as md
import pandas as pd
import openpyxl as pxl
import os
import os.path

n_frames = 200
for filename in os.listdir('./'):
    if filename.endswith('.pdb'):
        excipient_name = filename[17:-4]
        lst = []
        traj = md.load(filename)
        surface = md.shrake_rupley(traj, mode = 'residue')
        top = traj.topology
        n_residue = 16
        
        for i in range(n_frame):
            area = {}
            for j in range(n_residue):
                residues = top.residue(j)
                area[residues] = surface[i][j]
            lst.append(area)
            
        df = pd.DataFrame(lst)
        df.loc['total'] = df.sum()
        
        if not os.path.isfile('sasa.xlsx'):
            df.to_excel('sasa.xlsx', '%s' % excipient_name, index = True)
        
        else:
            excel_book = pxl.load_workbook('sasa.xlsx')
            with pd.ExcelWriter('sasa.xlsx', engine = 'openpyxl') as writer:
                writer.book = excel_book
                writer.sheets = {worksheet.title: worksheet for worksheet in excel_book.worksheets}
                df.to_excel(writer, '%s' % excipient_name, index = True)
                writer.save()
