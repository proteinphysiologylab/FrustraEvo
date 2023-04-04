#!/bin/bash
for value in N S M E orf3a orf3b orf6 orf7a orf7b orf8 orf9b orf9c orf10 nsp{1..16}
do
	./s3det_linux64 -i ../../msas/sub_$value.fasta -o ../../S3_det_results/2nd_round_$value.S3det
done
