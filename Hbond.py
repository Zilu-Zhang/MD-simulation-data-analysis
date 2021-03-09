import mdtraj as md
import sys

filepath = './openmm_sorafenib_glycyrrhizin.pdb'

def label(hbond):
    r1 = traj.topology.atom(hbond[0])
    r2 = traj.topology.atom(hbond[2])
    if ('ZPE' in str(r1)) is not ('ZPE' in str(r2)):
        n = 1
    else:
        n = 0
    return (n, '%s -- %s' % (r1, r2))

traj = md.load_frame(filepath, 0)
hbonds = md.baker_hubbard(traj)
number = 0
readout = open('readout.txt', 'a')
sys.stdout = readout
for hbond in hbonds:
    n, note = label(hbond)
    number += n
    print(note)
print('total Hbonds: ', len(hbonds))
print('drug-excipient Hbonds: ', number)
print('--------------------------------')
readout.close()

traj = md.load_frame(filepath, 199)
hbonds = md.baker_hubbard(traj)
number = 0
readout = open('readout.txt', 'a')
sys.stdout = readout
for hbond in hbonds:
    n, note = label(hbond)
    number += n
    print(note)
print('total Hbonds: ', len(hbonds))
print('drug-excipient Hbonds: ', number)
print('--------------------------------')
readout.close()
