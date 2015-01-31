from nose.tools import *
from nocap.Constraint import Constraint, ConstraintGroup

c = Constraint('CB', 'CB', 4, 12)

# Constraint
def test_to_string():
    eq_(str(c), "AtomPair CB 4 CB 12 HARMONIC 8.0 1.5")

def test_is_enabled_true():
    ok_(c.is_enabled(7.0), "should be enabled")

def test_is_enabled_false():
    ok_(not c.is_enabled(9.0), "should not be enabled")

def test_distance_on_chain():
    eq_(c.distance_on_chain(), 8, "distance computation failed")

# ConstraintGroup
def test_read_from_file():
    filename = 'tests/test_constraints.txt'
    cg = ConstraintGroup.read_from_file(filename)

    eq_(len(cg.constraints), 5)

def test_write_to_file():
    filename = 'tests/test_constraints_out.txt'

    #guarantee, that it doesn't exist
    import os.path
    from os import remove
    if os.path.isfile(filename):
        remove(filename)

    constraints = [c, c]
    cg = ConstraintGroup(constraints) # 2 constraints
    cg.write_to_file(filename)

    import subprocess
    eq_(int(subprocess.check_output(['wc', '-l', filename])[0]), len(constraints))
    ok_(os.path.isfile(filename))

    # clean up
    remove(filename)
