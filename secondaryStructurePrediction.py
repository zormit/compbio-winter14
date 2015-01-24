import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
import pymol
pymol.finish_launching()
from pymol import cmd

def extractSecondaryStructure(filename):
    label = 'targetStructure'
    cmd.load(filename, label)

    pymol.stored_ss = []
    cmd.iterate(label, 'stored_ss.append(ss)')

    secondaryStructure = ''.join(pymol.stored_ss)

    return secondaryStructure

s = extractSecondaryStructure("Data/input/2h3jA.pdb")
print s
print len(s)
