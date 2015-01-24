from __future__ import division
import numpy as np
import pdb
import matplotlib.pylab as plt

def visualizeConstraintDistances(constraintsFilename, proteinID):
    ''' generates a simple graph, which shows the constraint distances
        in an ascending manner

    Keyword arguments:
        constraintsFilename -- TODO: describe expected format.
    '''
    constraints = np.genfromtxt(constraintsFilename, dtype=float, usecols=(2,4,8))
    plt.figure()
    sorter = np.argsort(constraints[:,2])
    plt.plot(constraints[:,2][sorter])
    plt.plot((0,len(constraints)), (8,8), color="red")
    plt.xlabel("constraints")
    plt.ylabel("distance in Angstrom")
    #TODO: hardcoded
    plt.savefig('plots/{0}/{0}_constraintTrueDistances.png'.format(proteinID))
    #plt.show()


if __name__ == "__main__":
    proteinID = "2h3jA"
    proteinID = "2krkA"
    constraints = "/home/neeb/Data/input/{0}/generated/{0}_contact_constraints_withNativeDistances.txt".format(proteinID)
    visualizeConstraintDistances(constraints, proteinID)
