#!/bin/bash
dir_list=$(ssh compbio1 ls -d '/home/neeb/Data/output-sanitytest/{2h3jA,2krkA}/*')

for dir in $dir_list; do
    ldir='.'${dir:33}
    mkdir -p $ldir
    scp compbio1:$dir/*.fsc $ldir
done
