'''
This will generate a Hbonds.xlsx file containing information about total Hbonds, interactive Hbonds and the corresponding ratios
'''

# import packages
import os
import os.path
import mdtraj as md
import numpy as np
import pandas as pd
import openpyxl as pxl

# recognize the interactive H bonds between drug and excipient molecules
def label(hbond):
    r1 = str(traj.topology.atom(hbond[0]))
    r2 = str(traj.topology.atom(hbond[2]))
    n = 1 if r1[:3] != r2[:3] else 0
    if n == 1:
	m = 1 if r1[:3] == drug_code else 0
    return n, m

# main code
for filename in os.listdir('./'): 
    if filename.endswith('.pdb'):
	full = md.load(filename)
	n_frames = len(full)
	drug_code = str(full.topology.residue(0))
        excipient_name = filename[17:-4]
        frame_index = total_Hbonds = interactive_Hbonds = interactive_ratio = drug_donor_Hbonds = drug_donor_ratio = np.zeros(n_frames + 1)
	
        for i in range(n_frames):
            traj = full[i]
            hbonds = md.baker_hubbard(traj)
            number = 0
            donor = 0
            for hbond in hbonds:
                n, m = label(hbond)
                number += n
		donor += m
			
            frame_index[i] = i
            total_Hbonds[i] = len(hbonds)
            interactive_Hbonds[i] = number
            interactive_ratio[i] = interactive_Hbonds[i] / total_Hbonds[i]
            drug_donor_Hbonds[i] = donor
            drug_donor_ratio[i] = drug_donor[i] / interactive_Hbonds[i]


	
	frame_index[-1] = 'summary' # row title
	total_Hbonds[-1] = sum(total_Hbonds[:-2]) # total number of all Hbonds
	interactive_Hbonds[-1] = sum(interactive_Hbonds[:-2]) # total number of all interactive Hbonds
	interactive_ratio[-1] = interactive_Hbonds[-1] / total_Hbonds[-1] # overall interactive Hbonds ratios
	drug_donor_Hbonds[-1] = sum(drug_donor_Hbonds[:-2]) # total number of all interactive Hbonds
	drug_donor_ratio[-1] = drug_donor_Hbonds[-1] / interactive_Hbonds[-1] # overall drug-as-donor Hbonds ratios
	
        df = pd.DataFrame({'frame_index': frame_index, 
			   'total_Hbonds': total_Hbonds,
			   'interactive_Hbonds': interactive_Hbonds, 
			   'interactive_ratio':interactive_ratio,
			   'drug_donor_Hbonds':drug_donor_Hbonds,
			   'drug_donor_ratio': drug_donor_ratio})
		
        if not os.path.isfile('Hbonds.xlsx'):
            df.to_excel('Hbonds.xlsx', '%s' % excipient_name)

        else:
            excel_book = pxl.load_workbook('Hbonds.xlsx')
            with pd.ExcelWriter('Hbonds.xlsx', engine = 'openpyxl') as writer:
                writer.book = excel_book
                writer.sheets = {worksheet.title: worksheet for worksheet in excel_book.worksheets}
                df.to_excel(writer, '%s' % excipient_name, index = True)
                writer.save()
		writer.close()
