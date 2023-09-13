package main

import (
	"fmt"
)

func Test() {
	genome := LoadGenomes("BANAL-20-52.fasta", "BANAL-20-52.orfs")

	fmt.Printf("Loaded %d genomes length %d\n",
		genome.NumGenomes(), genome.Length())

	genome.Save("B52", "B52-test.fasta", 0)
	fmt.Printf("Saved as B52-test.fasta\n")

	mutant := genome.Clone()

	nd := NewNucDistro(genome)
	MutateSilent(mutant, nd, 153)

	mutant.Save("Mutant", "B52-mutated.fasta", 0)
	fmt.Printf("Saved as B52-mutated.fasta\n")
}
