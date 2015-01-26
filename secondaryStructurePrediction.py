import pymol
from pymol import cmd

def extractSecondaryStructure(filename):
    label = 'targetStructure'
    cmd.load(filename, label)

    pymol.stored_ss = []
    cmd.iterate(label, 'stored_ss.append(ss)')

    secondaryStructure = ''.join(pymol.stored_ss)

    return secondaryStructure
