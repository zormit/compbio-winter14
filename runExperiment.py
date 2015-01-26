from __future__ import division
from __future__ import print_function
import sys
import os
from os.path import isfile, isdir, join
import shutil
import pdb
import subprocess
import numpy as np
import random
import argparse
import matplotlib as mpl
mpl.use('Agg')
# http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pylab as plt
import logging
from time import strftime

# pymol has to be launched at the root-file-level, to not break stdout.
import pymol
pymol.pymol_argv = ['pymol', '-qc'] # Pymol: quiet and no GUI
pymol.finish_launching()

# own imports (after pymol launch!)
import constraintSubsets
import checkConstraints
import secondaryStructurePrediction

def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('inputPath', help='path to input directory')
    parser.add_argument('outputPath', help='path to output directory')
    parser.add_argument('proteinID', help='ID of protein (has to correspond with path&files)')
    parser.add_argument('-d', '--debug', help='run a light version of rosetta', action='store_true')
    parser.add_argument('-v','--verbose', help='be talkative', action="store_true")

    return parser.parse_args()

def setup_logger(args):
    # according to https://docs.python.org/2/howto/logging-cookbook.html

    logger = logging.getLogger('nocap')
    logger.setLevel(logging.DEBUG)

    # create file handler which logs even debug messages
    file_handler = logging.FileHandler('nocap.log')
    file_handler.setLevel(logging.DEBUG)

    # create console handler with a higher log level dependent on verbosity
    console_handler = logging.StreamHandler()
    if args.verbose:
        console_handler.setLevel(logging.INFO)
    else:
        console_handler.setLevel(logging.ERROR)


    # create formatter and add it to the handlers
    file_formatter = logging.Formatter('%(asctime)s (%(name)s - %(levelname)s): %(message)s',
                                       datefmt='%y-%m-%d %H:%M:%S')
    file_handler.setFormatter(file_formatter)
    console_formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console_handler.setFormatter(console_formatter)

    # add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger

def main(argv=None):

    # append argv to system argv if existing
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    try:
        # Process arguments
        args = parse_args()

        # start logging
        logger = setup_logger(args)
        logger.debug("start program with the following args: {}".format(args))

    except KeyboardInterrupt:
        ### handle keyboard interrupt silently ###
        return 0

if __name__ == "__main__":
    sys.exit(main())

# parameters
probabilityThreshold= 1.5

# the imports seem not to work, when loaded as module,
# that's why this function is copied here. TODO: fix somehow?
def extractSecondaryStructure(filename):
    label = 'targetStructure'
    #TODO: some more reasonable error message, when file doesn't exist
    cmd.load(filename, label)

    pymol.stored_ss = []
    cmd.iterate(label, 'stored_ss.append(ss)')

    secondaryStructure = ''.join(pymol.stored_ss)

    return secondaryStructure

# parse commandline arguments
inputPath = join(args.inputPath, args.proteinID)
outputPath = join(args.outputPath, args.proteinID)


# STEP 1: generate constraint groups, i.e. subsets of all constraints
logging.info("starting STEP 1: generate constraint groups")

nativeFile = "{}.pdb".format( join(inputPath, args.proteinID) )

secondaryStructure = extractSecondaryStructure(nativeFile)
sequenceLength = len(secondaryStructure)
logging.debug("secondary structure:{}".format(secondaryStructure))
logging.debug("sequence length:{}".format(sequenceLength))

constraintFile = checkConstraints.writeDistancesToConstraintFile(nativeFile, inputPath, args.proteinID)

adjacency, counts, unique = constraintSubsets.generateSSContraintSubsets(secondaryStructure, constraintFile)
numGroups = len(unique)

adjacencyRand = constraintSubsets.generateRandomGroups(len(secondaryStructure),constraintFile,counts,unique)

adjacency_nat = constraintSubsets.generateNativeContraints(constraintFile, sequenceLength)

# STEP 2: Run with all groups of constraints individually and compare with baseline

def runRosetta(constraintFile, inputPath, outputPath, proteinID, numGroups, groupSize, unique, adjacency, identifier=""):
    # Assign Probability
    constraints = []

    # load constraint file into dictionary
    with open(constraintFile) as baseConstraints:
        for i, l in enumerate(baseConstraints):
            row = l.split(' ')
            a, b = int(row[2])-1,int(row[4])-1
            c = {"positions":[a,b], "row":l, "probability":0}
            constraints.append(c)

    for group in range(numGroups):
        numConstraintsCurrent = 0
        for c in constraints:
            a,b = c["positions"]
            if adjacency[a,b] == unique[group]:
                # current constraint belongs to current group
                # => set to some weight
                c["probability"] = 1.5 # TODO: hardcoded! -> but never used as value: it just has to be >= threshold.
                numConstraintsCurrent += 1
            else:
                # reset previous weighting
                c["probability"] = 0

        sampleSet = range(numConstraintsCurrent)
        if not groupSize is None and groupSize <= numConstraintsCurrent:
            sampleSet = np.array(random.sample(sampleSet, groupSize))

        constraintIndex = 0
        # Create Constraint File with previously chosen constraints (according to group)
        prefix, suffix = constraintFile.rsplit('.',1)
        subConstraintsFilename = "{}_exp.{}".format(prefix,suffix)
        logging.debug("writing sub-constraints to:{}".format(subConstraintsFilename))
        with open(subConstraintsFilename, 'w') as subConstraints:
            for i,c in enumerate(constraints):
                if c["probability"] >= probabilityThreshold:
                    if constraintIndex in sampleSet:
                        if identifier == "nativesTrueConstraints":
                            newRow = c["row"].rstrip().split(" ")
                            newRow[7] = newRow[11]
                            newRow = " ".join(newRow[:-1])+"\n"
                            subConstraints.write(newRow)
                        else:
                            subConstraints.write(c["row"])
                    constraintIndex += 1
        # TODO: the previous three blocks might fit into one.

        if not args.debug:
            relaxFlags = ['-abinitio:relax',
                         '-relax:fast']
            nStruct = 20
        else:
            # debug mode => not relaxing yet ;)
            relaxFlags = []
            nStruct = 2

        # Standard Filenames.
        # TODO: put in some config?
        filePrefix = join(inputPath, proteinID)
        sequenceFile = "{}.fasta".format(filePrefix)
        frag3File = "{}.200.3mers".format(filePrefix)
        frag9File = "{}.200.9mers".format(filePrefix)
        nativeFile = "{}.pdb".format(filePrefix)

        FNULL = open(os.devnull, 'w')

        ## Run Prediction
        #TODO: hardcoded path
        logging.info("starting rosetta-run for group{}{} at {}".format(group, identifier,
                                                                    strftime('%H:%M:%S')))
        rosetta_cmd = ['/home/lassner/rosetta-3.5/rosetta_source/bin/AbinitioRelax.linuxgccrelease',
            '-in:file:fasta', sequenceFile,
            '-in:file:frag3', frag3File,
            '-in:file:frag9', frag9File,
            '-in:file:native', nativeFile,
            '-constraints:cst_file', subConstraintsFilename,
            '-out:nstruct', str(nStruct),
            '-out:pdb',
            '-mute core.io.database',
            '-database', '/home/lassner/rosetta-3.5/rosetta_database/']+relaxFlags
        logging.debug("executing rosetta with:{}".format(" ".join(rosetta_cmd)))
        subprocess.call(rosetta_cmd,
            stdout=FNULL, stderr=subprocess.STDOUT)

        ## Copy Results
        outputDir = join(outputPath, 'group{}{}'.format(group, identifier))
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        # backup subConstraints file with results
        subprocess.call(['cp',
            subConstraintsFilename,
            outputDir], stdout=FNULL, stderr=subprocess.STDOUT)

        inputDir = os.getcwd()
        #TODO: hardcoded file-matchers
        files = [ f for f in os.listdir(inputDir) if isfile(join(inputDir,f)) ]
        for f in files:
            if f.endswith('.pdb') or f.endswith('score.fsc'):
                shutil.move(join(inputDir,f),join(outputDir,f))

            if f.endswith('default.out'):
                os.remove(join(inputDir,f))
        logging.debug("moved results to {}. finished a rosetta-run".format(outputDir))

logging.info("starting STEP 2: run with all groups individually")

# Example Run for Rosetta native constraints only
runRosetta(constraintFile, inputPath, outputPath, args.proteinID, 1, None, [0], adjacency_nat, "natives")

# nativesTrueConstraints means, that the distances from the native are set as penalty
runRosetta(constraintFile, inputPath, outputPath, args.proteinID, 1, None, [0], adjacency_nat, "nativesTrueConstraints")


# Run for Rosetta Baseline:
runRosetta(constraintFile, inputPath, outputPath, args.proteinID, 1, None, [0], np.zeros((sequenceLength,sequenceLength),dtype=int)-1, "baseline")

# Run for Rosetta all groups of constraints individually
runRosetta(constraintFile, inputPath, outputPath, args.proteinID, numGroups, None, unique, adjacency, "individual")

runRosetta(constraintFile, inputPath, outputPath, args.proteinID, numGroups, None, unique, adjacencyRand, "random")

runRosetta(constraintFile, inputPath, outputPath, args.proteinID, numGroups, 8, unique, adjacency, "sizeCap")

# STEP 3: Compare scores and fulfillments of constraints
logging.info("starting STEP 3: compare scores and fulfillments of constraints")

dirIds= ['native', 'group0baseline', 'group0natives', 'group0nativesTrueConstraints']
dirs = []
for i in range(numGroups):
    dirIds.append('group'+str(i)+"individual")
    dirIds.append('group'+str(i)+"random")
    dirIds.append('group'+str(i)+"sizeCap")
for d in dirIds:
    if d in os.listdir(outputPath) and isdir(join(outputPath,d)):
        dirs.append(d)

for d in dirs:
    pdbfiles = [ join(outputPath,d,f) for f in os.listdir(join(outputPath,d)) if isfile(join(outputPath, d,f)) and f.endswith('.pdb')]
    pdbfiles.sort()
    ## Run different scoring
    logging.debug("rescoring {}".format(d))
    FNULL = open(os.devnull, 'w')
    subprocess.call(['/home/lassner/rosetta-3.5/rosetta_source/bin/score.linuxgccrelease',
        '-in:file:s']+pdbfiles+[
        '-in:file:fullatom',
        '-out:file:scorefile', join(outputPath, d, 'scoreV3.fsc'), # TODO: hardcoded.
        '-native', nativeFile,
        '-database', '/home/lassner/rosetta-3.5/rosetta_database/'],
        stdout=FNULL, stderr=subprocess.STDOUT)

scores = []
gdts = []
for d in dirs:
    # TODO: hardcoded.
    data = np.genfromtxt(join(outputPath,d)+'/scoreV3.fsc', dtype=None, usecols=(1,22), skip_header=1)
    scores.append(data[:,0])
    gdts.append(data[:,1])


scores = np.array(scores)
gdts = np.array(gdts)

def groupBoxplots(data, ylabel, groupnames, filename=None, ylim = None):
    fig, ax1 = plt.subplots()
    ax1.boxplot(data.T)
    xtickNames = plt.setp(ax1, xticklabels=groupnames)
    plt.setp(xtickNames, rotation=45, fontsize=8)
    ax1.set_xlabel("group identifiers")
    ax1.set_ylabel(ylabel)
    if not ylim is None:
        ax1.set_ylim(ylim)
    if filename is None:
        plt.show()
    else:
        plt.savefig(join("plots",args.proteinID,filename+".png"), bbox_inches='tight')

def gdtScoreScatterplot(gdts, scores, gdtsBaseline=None, scoresBaseline=None, filename=None):
    plt.figure()
    plt.title("{} {}".format(args.proteinID, filename))
    plt.scatter(gdts, scores)
    if scoresBaseline != None and gdtsBaseline!=None:
        plt.scatter(gdtsBaseline, scoresBaseline, c='grey', alpha=0.5)
    plt.xlabel('gdtmm')
    plt.xlim((0,1))
    plt.ylabel('score')
    if filename is None:
        plt.show()
    else:
        plt.savefig(join("plots",args.proteinID, "scatter",filename+".png"), bbox_inches='tight')

for group,label in enumerate(dirs):
    if "individual" in label:
        gdtScoreScatterplot(gdts[group], scores[group],
                            gdts[1], scores[1], label)

sizeCapGroupsIdx = []
randomGroupsIdx = []
individualGroupsIdx = []

sizeCapGroupsLabels = []
randomGroupsLabels = []
individualLabels = []

for k,label in enumerate(dirs):
    if "individual" in label:
        individualGroupsIdx.append(k)
        individualLabels.append(label)
    elif "random" in label:
        randomGroupsIdx.append(k)
        randomGroupsLabels.append(label)
    elif "sizeCap" in label:
        sizeCapGroupsIdx.append(k)
        sizeCapGroupsLabels.append(label)
    else:
        #add baseDirs to each of the groups
        individualGroupsIdx.append(k)
        individualLabels.append(label)
        randomGroupsIdx.append(k)
        randomGroupsLabels.append(label)
        sizeCapGroupsIdx.append(k)
        sizeCapGroupsLabels.append(label)


groupBoxplots(scores[individualGroupsIdx], "rosetta energy score", individualLabels, "scoresIndividual")
groupBoxplots(scores[randomGroupsIdx], "rosetta energy score", randomGroupsLabels, "scoresRandom")
groupBoxplots(scores[sizeCapGroupsIdx], "rosetta energy score", sizeCapGroupsLabels, "scoresSizeCap")

groupBoxplots(gdts[individualGroupsIdx], "gdt", individualLabels, "gdtsIndividual",ylim=(0,1))
groupBoxplots(gdts[randomGroupsIdx], "gdt", randomGroupsLabels, "gdtsRandom",ylim=(0,1))
groupBoxplots(gdts[sizeCapGroupsIdx], "gdt", sizeCapGroupsLabels, "gdtsSizeCap",ylim=(0,1))

prefix, suffix = constraintFile.rsplit('/',1)[1].rsplit('.',1)
subConstraintsFilename = "{}_exp.{}".format(prefix,suffix)

enabledConstraints = np.zeros(scores.shape)
nativeConstraints = np.zeros(scores.shape)
for group,inputDir in enumerate(dirs):
    files = [join(outputPath,inputDir,f) for f in os.listdir(join(outputPath,inputDir)) if isfile(join(outputPath,inputDir,f)) and f.endswith('.pdb') ]
    for decoy,f in enumerate(files):
        constraintsFilepath = join(join(outputPath,inputDir),subConstraintsFilename)
        constraintsEnabled = checkConstraints.constraintsEnabledInDecoy(f, constraintsFilepath)
        constraintsNative = checkConstraints.nativeConstraintsInDecoy(f, constraintsFilepath)
        enabledConstraints[group,decoy] = float(np.sum(constraintsEnabled))/len(constraintsEnabled)
        nativeConstraints[group,decoy] = float(np.sum(constraintsNative))/len(constraintsNative)

groupBoxplots(enabledConstraints[individualGroupsIdx], "ratio constraints enabled", individualLabels, "enabledConstraintsIndividual",ylim=(0,1))
groupBoxplots(enabledConstraints[randomGroupsIdx], "ratio constraints enabled", randomGroupsLabels, "enabledConstraintsRandom",ylim=(0,1))
groupBoxplots(enabledConstraints[sizeCapGroupsIdx], "ratio constraints enabled", sizeCapGroupsLabels, "enabledConstraintsSizeCap",ylim=(0,1))

# scatterplot enabled constraints vs native constraints
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.scatter(enabledConstraints.flatten(),nativeConstraints.flatten())
ax.set_xlabel("enabled constraints")
ax.set_ylabel("native constraints")
ax.set_xlim(0,1)
ax.set_ylim(0,1)
plt.savefig(join("plots",args.proteinID,"enabledVsNative.png"), bbox_inches='tight')

# STEP 3.5: get contact map
def contactMap(contactMatrix, contactMatrixNative, filename = None):
    truePositiveMatrix = np.logical_and(
                            contactMatrix,
                            contactMatrixNative)
    falsePositiveMatrix = np.logical_and(
                            contactMatrix,
                            np.logical_not(contactMatrixNative))
    truePositive = np.sum(truePositiveMatrix)
    falsePositive = np.sum(falsePositiveMatrix)
    falseNegatives = np.sum(np.logical_and(
                np.logical_not(contactMatrix),
                contactMatrixNative))

    precision = truePositive/(truePositive+falsePositive)
    recall = truePositive/(truePositive+falseNegatives)

    plt.figure()
    plt.scatter(*np.where(contactMatrixNative+contactMatrixNative.T), c='g', alpha=0.2)
    plt.scatter(*np.where(truePositiveMatrix.T), c='b', alpha=0.5)
    plt.scatter(*np.where(falsePositiveMatrix), c='r', alpha=0.5)
    plt.title("contact map")
    plt.xlabel("constraint position on backbone")
    plt.ylabel("constraint position on backbone")

    plt.text(5, 10, "precision:{:.2%}\nrecall:{:.2%}".format(precision,recall))

    # Set the ticks
    n = contactMatrix.shape[0]
    ticks = np.arange(0, n, 5)
    labels = range(ticks.size)
    plt.xticks(ticks)
    plt.yticks(ticks)
    plt.xlim(0, n)
    plt.ylim(0, n)

    if filename is None:
        plt.show()
    else:
        plt.savefig(join("plots",args.proteinID,"contactmaps", filename+".png"), bbox_inches='tight')

# STEP 3.6: get contact map for groups
def contactMapGroups(adjacency, adjacency_nat, filename = None):
    plt.figure()
    contactMatrixNative = adjacency_nat>=0
    plt.scatter(*np.where(contactMatrixNative+contactMatrixNative.T), c='g', alpha=0.2)

    def drawScatter(contactMatrix, contactMatrixNative, group, color='b'):
        truePositiveMatrix = np.logical_and(
                                contactMatrix,
                                contactMatrixNative)
        falsePositiveMatrix = np.logical_and(
                                contactMatrix,
                                np.logical_not(contactMatrixNative))
        truePositive = np.sum(truePositiveMatrix)
        falsePositive = np.sum(falsePositiveMatrix)
        falseNegatives = np.sum(np.logical_and(
                    np.logical_not(contactMatrix),
                    contactMatrixNative))

        precision = truePositive/(truePositive+falsePositive)
        recall = truePositive/(truePositive+falseNegatives)

        # http://matplotlib.org/examples/pylab_examples/scatter_star_poly.html
        if color == 'k':
            alpha=0.5
        else:
            alpha=1.0
        tp=plt.scatter(*np.where(truePositiveMatrix.T), c=color, alpha=alpha, marker='s',label='TP-group{}'.format(group))
        fp=plt.scatter(*np.where(falsePositiveMatrix), c=color, alpha=alpha, marker='^',label='FP-group{}'.format(group))

    plt.title("contact map")
    plt.xlabel("constraint position on backbone")
    plt.ylabel("constraint position on backbone")

    colors=['r','g','b','y','c','m', 'k']
    for group, groupid in enumerate(unique):
        drawScatter(adjacency==groupid, contactMatrixNative, group, color=colors[group])

    '''
    plt.legend(
        scatterpoints=1,
        loc='lower left',
        ncol=3,
        fontsize=8)
    '''

    # Set the ticks
    n = adjacency.shape[0]
    ticks = np.arange(0, n, 5)
    labels = range(ticks.size)
    plt.xticks(ticks)
    plt.yticks(ticks)
    plt.xlim(0, n)
    plt.ylim(0, n)

    if filename is None:
        plt.show()
    else:
        plt.savefig(join("plots",args.proteinID,"contactmaps", filename+".png"), bbox_inches='tight')
contactMapGroups(adjacency, adjacency_nat, 'allGroups')

for d in dirs:
    contactMatrix = np.zeros(adjacency.shape, dtype=bool)
    try:
        constraintPositions = np.genfromtxt(join(outputPath, d, subConstraintsFilename), dtype=None, usecols=(2,4))
    except IOError as e:
        constraintPositions = np.array([])

    if len(constraintPositions.shape) != 2:
        constraintPositions = constraintPositions[np.newaxis,:]
    if constraintPositions.shape[1] > 1:
        constraintPositions -= 1 # shift, to let it start by 0 instead 1
        contactMatrix[constraintPositions[:,0], constraintPositions[:,1]] = True
        contactMap(contactMatrix, adjacency_nat>=0, d)


def plotScoresAndEnabledConstraints(enabledConstraints, minScores):
    # [1:] to discard baseline
    baseline = minScores[0]
    sorter = np.argsort(enabledConstraints[1:])
    enabledConstraints = enabledConstraints[1:][sorter]
    minScores = minScores[1:][sorter]

    fig, ax1 = plt.subplots()
    ax1.plot(enabledConstraints, color="red")
    ax2 = ax1.twinx()
    ax2.plot(minScores, color="blue")
    ax2.plot((0,len(minScores)-1),(baseline,baseline),"b--")

    ax1.tick_params(axis='y', colors='red')
    ax2.tick_params(axis='y', colors='blue')

    ax1.set_xlabel("group identifier")
    ax1.set_ylabel("ratio of enabled constraints")
    ax2.set_ylabel("minimum energy score of the group")

    plt.show()

# TODO: this function fails at the moment anyway.
#plotScoresAndEnabledConstraints(enabledConstraints, minScores)


# STEP 4: Run predicition with combinations of promising groups
# TODO ^

