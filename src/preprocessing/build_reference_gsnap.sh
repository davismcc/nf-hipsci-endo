#!/bin/sh
reference = /nfs/research/stegle/datasets/references/human/gsnap-GRHC37_ERCC/GRCh37.p13.genome.chr_only_ERCC.fa
gmap_build -d human -k 15 -D /nfs/research/stegle/datasets/references/human/gsnap-GRHC37_ERCC $reference

