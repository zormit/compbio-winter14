from os.path import join
import matplotlib as mpl
mpl.use('Agg')
# http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pylab as plt
import numpy as np

plot_dir = './plots/2krkA'
def groupBoxplots(data, ylabel, groupnames, filename=None, ylim=None):
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
        plt.savefig(join(plot_dir, filename + ".png"), bbox_inches='tight')

def gdtScoreScatterplot(gdts, scores, gdtsBaseline=None, scoresBaseline=None, filename=None):
    plt.figure()
    plt.title(filename)
    plt.scatter(gdts, scores)
    if scoresBaseline != None and gdtsBaseline != None:
        plt.scatter(gdtsBaseline, scoresBaseline, c='grey', alpha=0.5)
    plt.xlabel('gdtmm')
    plt.xlim((0, 1))
    plt.ylabel('score')
    if filename is None:
        plt.show()
    else:
        plt.savefig(join(plot_dir, "scatter", filename + ".png"), bbox_inches='tight')

def enabled_native_scatterplot(enabledConstraints, nativeConstraints, filename=None):
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.scatter(enabledConstraints.flatten(), nativeConstraints.flatten())
    ax.set_xlabel("enabled constraints")
    ax.set_ylabel("native constraints")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if filename is None:
        plt.show()
    else:
        plt.savefig(join(plot_dir, filename + ".png"), bbox_inches='tight')

def contactMap(contactMatrix, contactMatrixNative, filename=None):
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

    precision = truePositive / (truePositive + falsePositive)
    recall = truePositive / (truePositive + falseNegatives)

    plt.figure()
    plt.scatter(*np.where(contactMatrixNative + contactMatrixNative.T), c='g', alpha=0.2)
    plt.scatter(*np.where(truePositiveMatrix.T), c='b', alpha=0.5)
    plt.scatter(*np.where(falsePositiveMatrix), c='r', alpha=0.5)
    plt.title("contact map")
    plt.xlabel("constraint position on backbone")
    plt.ylabel("constraint position on backbone")

    plt.text(5, 10, "precision:{:.2%}\nrecall:{:.2%}".format(precision, recall))

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
        plt.savefig(join(plot_dir, "contactmaps", filename + ".png"), bbox_inches='tight')

# STEP 3.6: get contact map for groups
def contactMapGroups(adjacency, adjacency_nat, unique, filename=None):
    plt.figure()
    contactMatrixNative = adjacency_nat >= 0
    plt.scatter(*np.where(contactMatrixNative + contactMatrixNative.T), c='g', alpha=0.2)

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

        precision = truePositive / (truePositive + falsePositive)
        recall = truePositive / (truePositive + falseNegatives)

        # http://matplotlib.org/examples/pylab_examples/scatter_star_poly.html
        if color == 'k':
            alpha = 0.5
        else:
            alpha = 1.0
        tp = plt.scatter(*np.where(truePositiveMatrix.T), c=color, alpha=alpha, marker='s',
                         label='TP-group{}'.format(group))
        fp = plt.scatter(*np.where(falsePositiveMatrix), c=color, alpha=alpha, marker='^',
                         label='FP-group{}'.format(group))

    plt.title("contact map")
    plt.xlabel("constraint position on backbone")
    plt.ylabel("constraint position on backbone")

    colors = ['r', 'g', 'b', 'y', 'c', 'm', 'k']
    for group, groupid in enumerate(unique):
        drawScatter(adjacency == groupid, contactMatrixNative, group, color=colors[group])

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
        plt.savefig(join(plot_dir, "contactmaps", filename + ".png"), bbox_inches='tight')
