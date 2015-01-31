from __future__ import division
from os import makedirs
from os.path import join, exists
import matplotlib as mpl

mpl.use('Agg')
# http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pylab as plt
import numpy as np

# if not overwritten by config, this is used
plot_dir = './plots/'

def custom_savefig(filename, subdir = '', extension = '.png'):
    save_dir = join(plot_dir, subdir)
    if not exists(save_dir):
        makedirs(save_dir)
    plt.savefig(join(save_dir, filename + extension), bbox_inches='tight')


def subset_boxplots(data, ylabel, groupnames, filename=None, ylim=None):
    fig, ax1 = plt.subplots()
    ax1.boxplot(data.T)
    xticknames = plt.setp(ax1, xticklabels=groupnames)
    plt.setp(xticknames, rotation=45, fontsize=8)
    ax1.set_xlabel("group identifiers")
    ax1.set_ylabel(ylabel)
    if ylim is not None:
        ax1.set_ylim(ylim)
    if filename is None:
        plt.show()
    else:
        custom_savefig(filename)


def gdt_score_scatter(gdts, scores, gdts_baseline=None, scores_baseline=None, filename=None):
    plt.figure()
    plt.title(filename)
    plt.scatter(gdts, scores)
    if scores_baseline is not None and gdts_baseline is not None:
        plt.scatter(gdts_baseline, scores_baseline, c='grey', alpha=0.5)
    plt.xlabel('gdtmm')
    plt.xlim((0, 1))
    plt.ylabel('score')
    if filename is None:
        plt.show()
    else:
        custom_savefig(filename, 'scatter')


def fulfilled_native_scatterplot(enabledConstraints, nativeConstraints, filename=None):
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
        custom_savefig(filename)


def contactmap(contactMatrix, contactMatrixNative, filename=None):
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
        custom_savefig(filename, 'contactmaps')


# STEP 3.6: get contact map for groups
def contactmap_subsets(subset_graph, subset_graph_native, unique, filename=None):
    plt.figure()
    contactmatrix_native = subset_graph_native >= 0
    plt.scatter(*np.where(contactmatrix_native + contactmatrix_native.T), c='g', alpha=0.2)

    def draw_scatter(contactmatrix, contactmatrix_native, group, color='b'):
        true_positive_matrix = np.logical_and(
            contactmatrix,
            contactmatrix_native)
        false_positive_matrix = np.logical_and(
            contactmatrix,
            np.logical_not(contactmatrix_native))
        true_positive = np.sum(true_positive_matrix)
        false_positive = np.sum(false_positive_matrix)
        false_negatives = np.sum(np.logical_and(
            np.logical_not(contactmatrix),
            contactmatrix_native))

        precision = true_positive / (true_positive + false_positive)
        recall = true_positive / (true_positive + false_negatives)

        # http://matplotlib.org/examples/pylab_examples/scatter_star_poly.html
        if color == 'k':
            alpha = 0.5
        else:
            alpha = 1.0
        tp = plt.scatter(*np.where(true_positive_matrix.T), c=color, alpha=alpha, marker='s',
                         label='TP-group{}'.format(group))
        fp = plt.scatter(*np.where(false_positive_matrix), c=color, alpha=alpha, marker='^',
                         label='FP-group{}'.format(group))

    plt.title("contact map")
    plt.xlabel("constraint position on backbone")
    plt.ylabel("constraint position on backbone")

    colors = ['r', 'g', 'b', 'y', 'c', 'm', 'k']
    for group, groupid in enumerate(unique):
        draw_scatter(subset_graph == groupid, contactmatrix_native, group, color=colors[group])

    '''
    plt.legend(
        scatterpoints=1,
        loc='lower left',
        ncol=3,
        fontsize=8)
    '''

    # Set the ticks
    n = subset_graph.shape[0]
    ticks = np.arange(0, n, 5)
    labels = range(ticks.size)
    plt.xticks(ticks)
    plt.yticks(ticks)
    plt.xlim(0, n)
    plt.ylim(0, n)

    if filename is None:
        plt.show()
    else:
        custom_savefig(filename, 'contactmaps')
