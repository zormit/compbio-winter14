import pymol
from pymol import cmd


def extract_secondary_structure(filename):
    label = 'targetStructure'
    cmd.load(filename, label)

    pymol.stored_ss = []
    cmd.iterate(label, 'stored_ss.append(ss)')

    secondary_structure = ''.join(pymol.stored_ss)

    return secondary_structure
