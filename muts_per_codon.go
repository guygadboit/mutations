package main

type MutsPerCodonResult struct {
	codonPos int // This is the n'th codon
	muts     int // How many muts (0, 1 or 2)
	expected int // "Expected" number based on possible alternatives
}

/*
	Given two aligned genomes find how many silent muts they have between them
	per codon. Ignore non-silent muts
*/
func MutsPerCodon(genome *Genomes, results chan *MutsPerCodonResult) {
}
