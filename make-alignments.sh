#!/bin/bash
relatives="ChimericAncestor.fasta \
	BtSY2.fasta \
	BANAL-20-236.fasta \
	BANAL-20-52.fasta \
	BANAL-20-103.fasta \
	RaTG13.fasta"

for rel in $relatives
do
	echo "Aligning WH1 with $rel..."
	cat WH1.fasta $rel | clustalo -i - > WH1-${rel%.*}.fasta
done
