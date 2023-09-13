package main

import (
	"flag"
	"fmt"
)

func Trials(genome *Genomes, nd *NucDistro, numTrials int, numMuts int) int {
	good := 0

	count, maxLength, unique := FindRestrictionMap(genome)
	fmt.Printf("Original: %d, %d, %t\n", count, maxLength, unique)

	reportProgress := func(n int) {
		fmt.Printf("Tested %d. Found %d/%d good mutants (%.2f%%)\n", n,
			good, n, float64(good*100)/float64(n))
	}

	for i := 0; i < numTrials; i++ {
		mutant := genome.Clone()
		MutateSilent(mutant, nd, numMuts)
		count, maxLength, unique = FindRestrictionMap(mutant)

		if unique && maxLength < 8000 {
			fmt.Printf("Mutant %d: %d, %d, %t\n", i, count, maxLength, unique)
			good += 1
		}

		if i%100 == 0 {
			reportProgress(i)
		}
	}

	reportProgress(numTrials)
	return good
}

func main() {
	var nTrials, nMuts int
	var test bool

	flag.IntVar(&nTrials, "n", 10000, "Number of trials")
	flag.IntVar(&nMuts, "m", 763, "Number of mutations per mutant")
	flag.BoolVar(&test, "t", false, "Just do some self-tests")
	flag.Parse()

	if test {
		Test()
		return
	}

	nd := NewNucDistro(nil)

	fnames := []string{
		"BtSY2",
		"BANAL-20-236",
		"BANAL-20-52",
		"RaTG13",
	}

	genomes := make([]*Genomes, len(fnames))
	for i := 0; i < len(fnames); i++ {
		genomes[i] = LoadGenomes(
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
		results[i] = Trials(genomes[i], nd, nTrials, nMuts)
	}

	for i := 0; i < len(genomes); i++ {
		fmt.Printf("%s: %d/%d %.2f%%\n", fnames[i], results[i], nTrials,
			float64(100.0*results[i])/float64(nTrials))
	}
}
