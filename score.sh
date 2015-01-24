#!/bin/bash
#
# Example calls of relaxing & scoring,
# to get a reasonable score native protein
/home/lassner/rosetta-3.5/rosetta_source/bin/relax.linuxgccrelease \
-in:file:s ~/Data/input/$1/$1.pdb \
-relax:fast \
-native ~/Data/input/$1/$1.pdb \
-database /home/lassner/rosetta-3.5/rosetta_database/ \
-out:file:silent $1_relaxed.pdb \
-out:nstruct 100

/home/lassner/rosetta-3.5/rosetta_source/bin/score.linuxgccrelease \
-in:file:silent $1_relaxed.pdb \
-in:file:fullatom \
-out:file:scorefile $1/scoreNative.fsc \
-native ~/Data/input/$1/$1.pdb \
-database /home/lassner/rosetta-3.5/rosetta_database/
