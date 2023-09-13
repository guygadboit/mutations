package main

import (
	"flag"
	"fmt"
)

func Trials(genome *Genome, nd *NucDistro, n int) int {
	good := 0

	count, maxLength, unique := FindRestrictionMap(genome)
	fmt.Printf("Original: %d, %d, %t\n", count, maxLength, unique)

	reportProgress := func(n int) {
		fmt.Printf("Tested %d. Found %d/%d good mutants (%.2f%%)\n", n,
			good, n, float64(good*100)/float64(n))
	}

	for i := 0; i < n; i++ {
		mutant := genome.Clone()
		MutateSilent(mutant, nd, 763)
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

func main() {
	var nTrials int
	flag.IntVar(&nTrials, "n", 10000, "Number of trials")
	flag.Parse()

	nd := NewNucDistro(nil)

	fnames := []string{
		"BtSY2",
		"BANAL-20-236",
		"BANAL-20-52",
		"RaTG13",
	}

	genomes := make([]*Genome, len(fnames))
	for i := 0; i < len(fnames); i++ {
		genomes[i] = LoadGenome(
			fnames[i]+".fasta",
			fnames[i]+".orfs",
		)
	}

	for i := 0; i < len(genomes); i++ {
		nd.Count(genomes[i])
	}
	nd.Show()

	results := make([]int, len(genomes))

	for i := 0; i < len(genomes); i++ {
		results[i] = Trials(genomes[i], nd, nTrials)
	}

	for i := 0; i < len(genomes); i++ {
		fmt.Printf("%s: %d/%d\n", fnames[i], results[i], nTrials)
	}
}
