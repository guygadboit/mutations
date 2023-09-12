package main

import (
	"fmt"
)

func Trials(n int) int {
	genomes := LoadGenomes("./BANAL-20-52.fasta")
	orfs := LoadOrfs("./BANAL-20-52.orfs")
	b52 := genomes[0]
	nd := NewNucDistro(b52)
	mutant := make(Genome, len(b52))
	good := 0

	count, maxLength, unique := FindRestrictionMap(b52)
	fmt.Printf("Original: %d, %d, %t\n", count, maxLength, unique)

	reportProgress := func(n int) {
		fmt.Printf("Tested %d. Found %d/%d good mutants (%.2f%%)\n", n,
			good, n, float64(good*100)/float64(n))
	}

	for i := 0; i < n; i++ {
		copy(mutant, b52)
		MutateSilent(mutant, orfs, nd, 763)
		count, maxLength, unique = FindRestrictionMap(mutant)

		if unique && maxLength < 8000 {
			fmt.Printf("Mutant %d: %d, %d, %t\n", i, count, maxLength, unique)
			good += 1
		}

		if i%100 == 0 {
			reportProgress(i)
		}
	}

	reportProgress(n)
	return good
}