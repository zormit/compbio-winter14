from __future__ import division
from pymol import cmd
import numpy as np
import pdb
import os


def nativeConstraintsInDecoy(decoyFilename, constraintsFilename):
    try:
        constraints = np.genfromtxt(constraintsFilename, dtype=None, usecols=(6, 7, 11))
    except IOError as e:
    # no constraints and old numpy -> circumvent with "one constraint disabled"
        return [False]

    # no constraints -> circumvent with "one constraint disabled"
    if constraints.shape[0] == 0:
        return [False]

    lowerBounds = constraints[:, 0]
    upperBounds = constraints[:, 1]
    nativeDistances = constraints[:, 2]

    constraintsEnabled = np.logical_and(lowerBounds <= nativeDistances,
                                        nativeDistances <= upperBounds)

    return constraintsEnabled



def constraintsEnabledInDecoy(decoyFilename, constraintsFilename, threshold = 8.0):
    cmd.reinitialize()
    decoyLabel = 'decoy_0'
    cmd.load(decoyFilename, decoyLabel)
    decoyDistances = list()
    try:
        constraints = np.genfromtxt(constraintsFilename, dtype=None, usecols=(1,2,3,4,6,7))
    except IOError as e:
    # no constraints and old numpy -> circumvent with "one constraint disabled"
        return [False]

    # no constraints -> circumvent with "one constraint disabled"
    if constraints.shape[0] == 0:
        return [False]

    lowerBounds = []
    upperBounds = []
    for atomA, posA, atomB, posB, lowerBound, upperBound in constraints:
        template = "{}///{}/{}"
        lowerBounds.append(lowerBound)
        upperBounds.append(upperBound)
        # cmd.distance would draw the constraints
        decoyDistances.append(cmd.get_distance(
                    template.format(decoyLabel,posA,atomA),
                    template.format(decoyLabel,posB,atomB)))
    decoyDistances = np.array(decoyDistances)
    lowerBounds = np.array(lowerBounds)
    upperBounds = np.array(upperBounds)

    # enabled means, it is inside the bounds.
    constraintsEnabled = np.logical_and(lowerBounds <= decoyDistances, decoyDistances <= upperBounds)

    return constraintsEnabled

def writeDistancesToConstraintFile(realProteinFilename, inputPath, proteinID, logger):
    groundTruthLabel = 'groundTruth'
    cmd.load(realProteinFilename, groundTruthLabel)

    #TODO pull suffix to config? (cf. runExperiment -> runRosetta())
    constraintsFilename= "{}_contact_constraints.txt".format( os.path.join(inputPath, proteinID) )

    groundTruthDistances = list()
    # load constraints
    constraints = np.genfromtxt(constraintsFilename, dtype=None, usecols=(1,2,3,4))
    for atomA, posA, atomB, posB in constraints:
        template = "{}///{}/{}"
        groundTruthDistances.append(cmd.get_distance(
                    template.format(groundTruthLabel,posA,atomA),
                    template.format(groundTruthLabel,posB,atomB)))

    # do some path & filename magic
    constraintWithDistFilename = os.path.join(
            inputPath, 'generated',
            '{}_contact_constraints_withNativeDistances.txt'.format(proteinID))

    logger.debug("writing file: {}".format(constraintWithDistFilename))

    # add distance to each line of the old file
    with open(constraintWithDistFilename, 'w') as newConstraintsFile:
        with open(constraintsFilename, 'r') as oldConstraintsFile:
            i = 0
            for line in oldConstraintsFile:
                newConstraintsFile.write( "{} {}\n".format(line[:-1], groundTruthDistances[i]) )
                i +=1

    return constraintWithDistFilename

def checkConstraints(realProteinFilename, decoyFilename, constraintsFilename, threshold = 8.0):
    groundTruthLabel = 'groundTruth'
    decoyLabel = 'decoy_0'
    cmd.load(realProteinFilename, groundTruthLabel)
    cmd.load(decoyFilename, decoyLabel)

    decoyDistances = list()
    groundTruthDistances = list()
    # load constraints
    constraints = np.genfromtxt(constraintsFilename, dtype=None, usecols=(1,2,3,4))
    for atomA, posA, atomB, posB in constraints:
        template = "{}///{}/{}"
        # cmd.distance would draw the constraints
        decoyDistances.append(cmd.get_distance(
                    template.format(decoyLabel,posA,atomA),
                    template.format(decoyLabel,posB,atomB)))
        groundTruthDistances.append(cmd.get_distance(
                    template.format(groundTruthLabel,posA,atomA),
                    template.format(groundTruthLabel,posB,atomB)))

    decoyDistances = np.array(decoyDistances)
    groundTruthDistances = np.array(groundTruthDistances)

    trueConstraints = groundTruthDistances<=threshold
    positive = decoyDistances<=threshold

    truePositive = np.sum(np.logical_and(
                positive,
                trueConstraints))
    falsePositive = np.sum(np.logical_and(
                positive,
                np.logical_not(trueConstraints)))
    falseNegatives = np.sum(np.logical_and(
                np.logical_not(positive),
                trueConstraints))

    precision = truePositive/(truePositive+falsePositive)
    recall = truePositive/(truePositive+falseNegatives)

    return precision, recall

if __name__ == "__main__":
    decoy = "/home/neeb/Data/input-bounded/2h3jA/2h3jA.pdb"
    constraints = "/home/neeb/Data/input-bounded/2h3jA/2h3jA_contact_constraints.txt"
    print constraintsEnabledInDecoy(decoy, constraints)
