import numpy as np

class Constraint():
    ''' representation of a rosetta constraint '''

    def __init__(self, atomA, atomB, residual_numberA, residual_numberB, x=8.0, weight=1.5):
        ''' constructor
            Usage:
                Constraint('CB', 'CA', 4, 12, 8.0, 4.0)
            Parameters:
                atomA: first atom identifier
                atomB: second atom identifier
                residual_numberA: corresponding residual number
                residual_numberB: corresponding residual number
                x: the target distance of the constraint
                weight: weight of the constraint
        '''
        self.atoms = (atomA, atomB)
        self.residual_numbers = (residual_numberA, residual_numberB)
        self.x = x
        self.weight = weight

    def __str__(self):
        return "AtomPair {} {} {} {} HARMONIC {} {}".format(
                self.atoms[0],
                self.residual_numbers[0],
                self.atoms[1],
                self.residual_numbers[1],
                self.x,
                self.weight)

    def is_enabled(self, distance_in_decoy):
        ''' when the distance in the prediction is lower or equal,
            the constraint is defined as enabled.
        '''
        return distance_in_decoy <= self.x

    def distance_on_chain(self):
        return abs(self.residual_numbers[0]-self.residual_numbers[1])

class ConstraintGroup():
    ''' representation of a rosetta constraint set '''

    def __init__(self, constraints):
        ''' Parameters:
                constraints: a list of constraint objects
        '''
        self.constraints = constraints

    def write_to_file(self, filename):
        with open(filename, 'w') as new_file:
            for constraint in self.constraints:
                new_file.write(str(constraint)+'\n')

    @staticmethod
    def read_from_file(filename):
        ''' constraint group factory '''
        raw_constraints = np.genfromtxt(filename, dtype=None, usecols=(1,3,2,4,6,7))
        constraints = [Constraint(*constraint) for constraint in raw_constraints]
        return ConstraintGroup(constraints)
