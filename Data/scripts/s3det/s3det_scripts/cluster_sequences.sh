#!/bin/bash
for value in orf3a orf3b orf6 orf7a orf7b orf8 orf9b orf9c orf10 E M N S nsp{1..16}
do
	cd-hit -i msas/$value.fasta -o cd_hit_clusters/$value.c98 -c 0.98 -s 0.9
done
