import numpy as np
import scipy.stats
import random

def generateSSContraintSubsets(secondaryStructure, constraintsFilename, logger, numGroups=None):

    n = len(secondaryStructure)

    constraintPositions = np.genfromtxt(constraintsFilename, dtype=None, usecols=(2,4))
    constraintPositions -= 1 # shift, to let it start by 0 instead 1

    prev = secondaryStructure[0]
    groups = np.zeros(n,dtype=int)
    group = 0
    for pos, sstype in enumerate(secondaryStructure):
        if prev != sstype:
            group+=1
            prev=sstype
        groups[pos]=group
        if sstype == 'L':
            groups[pos] = -1
            constraintPositions = np.delete(constraintPositions,
                                            np.where(constraintPositions==pos)[0], 0)

    adjacency = np.zeros((n,n),dtype=int)-1 # initialize with -1
    adjacency2 = np.zeros((n,n),dtype=int)-1 # initialize with -1

    # compute a unique number for the combination of two groups
    adjacency[constraintPositions[:,0],constraintPositions[:,1]] = (
        groups[constraintPositions[:,0]] + groups[constraintPositions[:,1]] * (groups.max()+1))

    frequencies = scipy.stats.itemfreq(adjacency.flatten())
    unique,counts = frequencies[:,0],frequencies[:,1]
    # zeroth position is no constraint
    unique, counts = unique[1:], counts[1:]

    maxNumGroups = len(unique)
    logger.info("there are {} different SS groups".format(maxNumGroups))

    if numGroups==None:
        numGroups=maxNumGroups
    chosenGroups = np.argsort(counts)[-numGroups:]

    for group in chosenGroups[-numGroups:]:
        # this is not totally sane (basefile gets opened and iterated three times), but it works
        with open(constraintsFilename) as baseConstraints:
            # prefix, suffix = constraintsFilename.rsplit('.',1)
            # subConstraintsFilename = "{}_ss{:03}.{}".format(prefix,group,suffix)
            # with open(subConstraintsFilename, 'w') as subConstraints:
            for i, l in enumerate(baseConstraints):
                row = l.split(' ')
                a, b = int(row[2])-1,int(row[4])-1
                if adjacency[a,b] == unique[group]:
                    adjacency2[a,b] = unique[group]
                    # subConstraints.write(l)

    unique = unique[chosenGroups]
    counts = counts[chosenGroups]

    return adjacency2, counts, unique

def generateRandomGroups(proteinLength, constraintsFilename, counts, unique):
    constraintPositions = np.genfromtxt(constraintsFilename, dtype=None, usecols=(2,4))
    constraintPositions -= 1 # shift, to let it start by 0 instead 1

    adjacency = np.zeros((proteinLength,proteinLength),dtype=int)-1 # initialize with -1

    takenPositions = set()
    allPositions = range(len(constraintPositions))
    for group,size in enumerate(counts):
        possiblePositions = [pos for pos in allPositions if not pos in takenPositions]
        sample = random.sample(possiblePositions,int(size))
        takenPositions.update(sample)
        for constraintPosition in constraintPositions[sample]:
            adjacency[constraintPosition[0],constraintPosition[1]] = unique[group]

    #TODO: this redundancy is necessary because of our group-system via adjacency,
    #      but not neccessarily what we want.
    return adjacency



def generatePSContraintSubsets(constraintsFilename, maxDistance):

    constraintPositions = np.genfromtxt(constraintsFilename, dtype=None, usecols=(2,4))
    constraintPositions -= 1 # shift, to let it start by 0 instead 1
    distances = np.abs(constraintPositions[:,1] - constraintPositions[:,0])

    with open(constraintsFilename) as baseConstraints:
        prefix, suffix = constraintsFilename.rsplit('.',1)
        subConstraintsFilename = "{}_ps{:03}.{}".format(prefix,maxDistance,suffix)
        with open(subConstraintsFilename, 'w') as subConstraints:
            for i, l in enumerate(baseConstraints):
                row = l.split(' ')
                a, b = int(row[2])-1,int(row[4])-1
                if  np.abs(a-b) <= maxDistance:
                    subConstraints.write(l)

    return distances,constraintPositions,distances<=maxDistance

def generateOccurencyConstraintScores(constraintsFilename, length, window):

    constraintPositions = np.genfromtxt(constraintsFilename, dtype=None, usecols=(2,4))
    constraintPositions -= 1 # shift, to let it start by 0 instead 1

    constraintScore = np.zeros(np.max(constraintPositions)+1)
    occurrences = np.bincount(constraintPositions.flatten())
    slices = (occurrences[i:i+window] for i in xrange(0, len(occurrences), 1))
    for i,s in enumerate(slices):
        constraintScore[i:i+window] += np.sum(s)

    return constraintScore


def generateNativeContraints(constraintsFilename, n=75):
    constraintPositions = np.genfromtxt(constraintsFilename, dtype=float, usecols=(2,4,11))
    constraintPositions[:,:2] -= 1 # shift, to let it start by 0 instead 1
    adjacency = np.zeros((n,n),dtype=int) # initialize with 0

    adjacency[constraintPositions[:,0].astype(int),constraintPositions[:,1].astype(int)] = constraintPositions[:,2] <= 8.0
    adjacency -= 1
    return adjacency
