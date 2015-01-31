#!/usr/bin/python
from __future__ import division
from __future__ import print_function
import sys
import os
from os.path import isfile, join
import subprocess
import argparse
from ConfigParser import ConfigParser
import numpy as np
import logging
from time import strftime

# pymol has to be launched at the root-file-level, to not break stdout.
import pymol

pymol.pymol_argv = ['pymol', '-qc']  # Pymol: quiet and no GUI
pymol.finish_launching()

# own imports (after pymol launch!)
import constraintSubsets
import checkConstraints
import secondaryStructurePrediction
import plotting


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('config_filename', help='path to config')
    parser.add_argument('protein_ID',
                        help='ID of protein (has to correspond with path&files)')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='run a light version of rosetta')
    parser.add_argument('-v', '--verbose', action="store_true",
                        help='be talkative')

    return parser.parse_args()


def parse_config(config_filename):
    config = ConfigParser()
    config.read(config_filename)

    return config


def setup_logger(args, config):
    # according to https://docs.python.org/2/howto/logging-cookbook.html

    logger = logging.getLogger('nocap')
    logger.setLevel(logging.DEBUG)

    # create file handler which logs even debug messages
    file_handler = logging.FileHandler(config.get('filename', 'log'))
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


def generate_constraint_subsets(input_path, protein_ID, logger):
    # STEP 1: generate constraint groups, i.e. subsets of all constraints
    logger.info("starting: generate constraint groups")

    native_filename = "{}.pdb".format(join(input_path, protein_ID))

    secondary_structure = secondaryStructurePrediction.extract_secondary_structure(native_filename)
    sequenceLength = len(secondary_structure)
    logger.info("secondary structure:{}".format(secondary_structure))
    logger.debug("sequence length:{}".format(sequenceLength))

    constraint_filename = checkConstraints.writeDistancesToConstraintFile(
        native_filename, input_path, protein_ID, logger)

    subset_graph_SS, counts, subset_IDs = constraintSubsets.generateSSContraintSubsets(
        secondary_structure, constraint_filename, logger)

    subset_graph_rand = constraintSubsets.generateRandomGroups(sequenceLength,
                                                               constraint_filename, counts, subset_IDs)
    subset_graph_native = constraintSubsets.generateNativeContraints(constraint_filename, sequenceLength)

    subset_graph_all = np.zeros(subset_graph_native.shape, dtype=int)
    subset_graph_baseline = subset_graph_all - 1

    subset_IDs_all = [subset_IDs, subset_IDs, [0], [0], [0]]
    subset_labels = ["secstruct", "random", "native", "all", "baseline"]
    subset_graphs = [subset_graph_SS, subset_graph_rand,
                     subset_graph_native, subset_graph_all, subset_graph_baseline]

    logger.info("finished: generate constraint groups")

    return constraint_filename, subset_graphs, subset_IDs_all, subset_labels


def write_constraint_subset_files(output_path, constraint_filename,
                                  subset_graph, subset_IDs, dir_prefix,
                                  logger, config):
    constraint_residues = np.genfromtxt(constraint_filename,
                                        dtype=None, usecols=(2, 4))
    constraint_residues -= 1  # shift, to let it start by 0 instead 1

    # fill subset ID back to the each line
    subset_ID_mask = subset_graph[constraint_residues[:, 0], constraint_residues[:, 1]]

    # load file as array of lines
    with open(constraint_filename) as constraint_file:
        constraints = np.array(constraint_file.readlines())

    output_dirs = []

    for group, subset_ID in enumerate(subset_IDs):

        # create directory for results and constraint subset-file
        output_dir = join(output_path, '{}_{}'.format(dir_prefix, group))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_dirs.append(output_dir)

        # copy subset to its own file. name generated from base constraint file
        subset_filename = join(output_dir,
                               config.get('filename', 'constraint_subset'))
        logger.debug("writing sub-constraints to:{}".format(subset_filename))
        with open(subset_filename, 'w') as subset:
            # TODO allow for changed weights
            subset_indices = np.where(subset_ID_mask == subset_ID)
            subset.writelines(
                constraints[subset_indices])

    return output_dirs


def protein_structure_prediction(input_path, output_paths, protein_ID,
                                 logger, debug, config):
    logger.info("starting: PSP for all constraint subsets")

    for output_dir in output_paths:

        if not debug:
            relax_flags = ['-abinitio:relax',
                          '-relax:fast']
            number_of_decoys = 20
        else:
            # debug mode => not relaxing yet ;)
            relax_flags = []
            number_of_decoys = 2

        # Standard Filenames.
        # TODO: put in some config?
        file_prefix = join(input_path, protein_ID)
        sequence_filename = "{}.fasta".format(file_prefix)
        frag3_filename = "{}.200.3mers".format(file_prefix)
        frag9_filename = "{}.200.9mers".format(file_prefix)
        native_filename = "{}.pdb".format(file_prefix)
        subset_filename = join(output_dir,
                               config.get('filename', 'constraint_subset'))

        # Run Prediction
        logger.info(("starting rosetta-run on {} at {}").format(
            subset_filename, strftime('%H:%M:%S')))
        rosetta_cmd = [config.get('rosetta', 'abinitio'),
                       '-in:file:fasta', sequence_filename,
                       '-in:file:frag3', frag3_filename,
                       '-in:file:frag9', frag9_filename,
                       '-in:file:native', native_filename,
                       '-constraints:cst_file', subset_filename,
                       '-out:nstruct', str(number_of_decoys),
                       '-out:pdb',
                       '-out:path', output_dir,
                       '-out:sf', join(output_dir, 'score.fsc'),
                       '-out:file:silent', join(output_dir, 'default.out'),
                       '-mute core.io.database',
                       '-database', '/home/lassner/rosetta-3.5/rosetta_database/'] + relax_flags
        logger.debug("executing rosetta with:{}".format(" ".join(rosetta_cmd)))
        try:
            subprocess.call(rosetta_cmd, stdout=open(os.devnull, 'w'),
                            stderr=subprocess.STDOUT)
        except OSError:
            logger.exception("error calling rosetta. skipping.")
            break


def rescore_prediction(input_path, output_paths, protein_ID, logger, config):
    logger.info("starting: rescore all predictions")

    native_file = "{}.pdb".format(join(input_path, protein_ID))

    for output_dir in output_paths:
        pdbfiles = [join(output_dir, f) for f in os.listdir(output_dir)
                    if isfile(join(output_dir, f)) and f.endswith('.pdb')]
        pdbfiles.sort()

        logger.debug("rescoring {}".format(output_dir))

        scoring_cmd = [config.get('rosetta', 'score'),
                       '-in:file:s'] + pdbfiles + [
                       '-in:file:fullatom',
                       '-out:file:scorefile', join(output_dir,
                                    config.get('filename', 'rescore')),
                       '-native', native_file,
                       '-database', config.get('rosetta', 'database')]
        try:
            subprocess.call(scoring_cmd, stdout=open(os.devnull, 'w'),
                            stderr=subprocess.STDOUT)
        except OSError:
            logger.exception("error calling rosetta. skipping.")
            break


def plot(output_dirs_grouped, config):
    # STEP 3: Compare scores and fulfillments of constraints
    logging.info("starting STEP 3: compare scores and fulfillments of constraints")

    scores = {}
    gdts = {}
    subset_names = {}

    for group in output_dirs_grouped:
        scores[group] = []
        gdts[group] = []
        subset_names[group] = []

        for output_dir in output_dirs_grouped[group]:
            data = np.genfromtxt(
                join(output_dir, config.get('filename', 'rescore')),
                dtype=None, usecols=(1, 22), skip_header=1)
            scores[group].append(data[:, 0])
            gdts[group].append(data[:, 1])

            subset_names[group].append(os.path.basename(output_dir))

        scores[group] = np.array(scores[group])
        gdts[group] = np.array(gdts[group])

    target = 'secstruct'
    baseline = 'baseline'  # could also be 'all' or 'native'
    for i, name in enumerate(subset_names[target]):
        plotting.gdt_score_scatter(gdts[target][i], scores[target][i],
                                     gdts[baseline][0], scores[baseline][0],
                                     name)

    # individually selected baseline
    scores['combined_baseline'] = np.r_[scores['baseline'], scores['native']]

    gdts['combined_baseline'] = np.r_[gdts['baseline'], gdts['native']]
    subset_names['combined_baseline'] = (subset_names['baseline'] +
                                         subset_names['native'])

    baseline = 'combined_baseline'
    for target in ['secstruct', 'random']:
        partial_scores = np.r_[scores[baseline], scores[target]]
        partial_gdts = np.r_[gdts[baseline], gdts[target]]
        partial_subset_names = subset_names[baseline] + subset_names[target]

        plotting.subset_boxplots(partial_scores, "rosetta energy score",
                               partial_subset_names, "scores_" + target)
        plotting.subset_boxplots(partial_gdts, "gdt",
                               partial_subset_names, "gdts_" + target, ylim=(0, 1))

    ratio_fulfilled_constraints = {}
    ratio_native_constraints = {}
    for group in output_dirs_grouped:
        ratio_fulfilled_constraints[group] = np.zeros(scores[group].shape)
        ratio_native_constraints[group] = np.zeros(scores[group].shape)
        for i, output_dir in enumerate(output_dirs_grouped[group]):
            constraints_filename = join(output_dir,
                                        config.get('filename', 'constraint_subset'))

            native_constraints_subset = checkConstraints.native_constraints_subset(constraints_filename)
            ratio_native_constraints_subset = np.sum(native_constraints_subset) / len(native_constraints_subset)
            ratio_native_constraints[group][i, :] = ratio_native_constraints_subset

            pdbfiles = [join(output_dir, f) for f in os.listdir(output_dir)
                        if isfile(join(output_dir, f)) and f.endswith('.pdb')]
            for decoy, f in enumerate(pdbfiles):
                fulfilled_constraints_decoy = checkConstraints.constraints_fulfilled_decoy(f, constraints_filename)
                ratio_fulfilled_constraints[group][i, decoy] = np.sum(fulfilled_constraints_decoy) / len(fulfilled_constraints_decoy)

    baseline = 'native'
    for target in ['secstruct', 'random']:
        ratio_fulfilled_constraints_merged = np.r_[ratio_fulfilled_constraints[baseline],
                                           ratio_fulfilled_constraints[target]]
        partial_subset_names = subset_names[baseline] + subset_names[target]
        plotting.subset_boxplots(ratio_fulfilled_constraints_merged,
                               "ratio constraints fulfilled", partial_subset_names,
                               "constraints_fulfilled_" + target, ylim=(0, 1))

        # scatterplot enabled constraints vs native constraints
        plotting.fulfilled_native_scatterplot(ratio_fulfilled_constraints[target],
                                            ratio_native_constraints[target], "enabledVsNative_" + target)


def plot_contact_maps(subset_graphs, subset_labels, subset_IDs, output_dirs_grouped, config):
    idx_target = subset_labels.index('secstruct')
    idx_native = subset_labels.index('native')

    plotting.contactmap_subsets(subset_graphs[idx_target],
                              subset_graphs[idx_native],
                              subset_IDs[idx_target], 'secstruct_all')

    for group in output_dirs_grouped:
        for output_dir in output_dirs_grouped[group]:
            subset_graph_target = np.zeros(subset_graphs[0].shape, dtype=bool)
            try:
                constraint_residues = np.genfromtxt(
                    join(output_dir, config.get('filename', 'constraint_subset')),
                    dtype=None, usecols=(2, 4))
            except IOError as e:
                constraint_residues = np.array([])

            if len(constraint_residues.shape) != 2:
                constraint_residues = constraint_residues[np.newaxis, :]
            if constraint_residues.shape[1] > 1:
                constraint_residues -= 1  # shift, to let it start by 0 instead 1
                subset_graph_target[constraint_residues[:, 0], constraint_residues[:, 1]] = True
                plotting.contactmap(subset_graph_target,
                                    subset_graphs[idx_native] >= 0,
                                    os.path.basename(output_dir))


def main(argv=None):
    # append argv to system argv if existing
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    try:
        # Process arguments and config
        args = parse_args()
        config = parse_config(args.config_filename)

        # start logging
        logger = setup_logger(args, config)
        logger.debug("start program with the following args: {}".format(args))

        # append protein ID to paths
        protein_input_path = join(config.get('dir', 'input'), args.protein_ID)
        protein_output_path = join(config.get('dir', 'output'), args.protein_ID)
        plot_dir = join(config.get('dir', 'plot'), args.protein_ID)

        # generate constraint subsets based on secondary structure
        constraint_filename, subset_graphs, subset_IDs, subset_labels = (
            generate_constraint_subsets(protein_input_path,
                                        args.protein_ID, logger))

        prediction_paths_grouped = {}
        for i, label in enumerate(subset_labels):
            # write file for every subset of constraints
            prediction_paths = write_constraint_subset_files(
                protein_output_path, constraint_filename,
                subset_graphs[i], subset_IDs[i], label, logger, config)
            prediction_paths_grouped[label] = prediction_paths

        # reduce grouped paths to a consecutive list
        prediction_paths_all = sum(prediction_paths_grouped.values(), [])

        protein_structure_prediction(
            protein_input_path, prediction_paths_all, args.protein_ID,
            logger, args.debug, config)

        rescore_prediction(protein_input_path, prediction_paths_all,
                           args.protein_ID, logger, config)

        # plot-a-lot
        plotting.plot_dir = plot_dir
        plot(prediction_paths_grouped, config)
        plot_contact_maps(subset_graphs, subset_labels, subset_IDs,
                          prediction_paths_grouped, config)
        plotting.constraint_distances_graph(constraint_filename, "distances")

    except KeyboardInterrupt:
        logger.info("terminating...")
        return 0


if __name__ == "__main__":
    sys.exit(main())



