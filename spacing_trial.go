package main

import (
	"fmt"
	"io"
)

type SpacingTrial struct {
	runFunc func(genome *Genomes, numMuts int, results chan interface{})
}

func (t *SpacingTrial) Run(genome *Genomes, numMuts int,
	results chan interface{}) {
	t.runFunc(genome, numMuts, results)
}

func (t *SpacingTrial) WriteHeadings(w io.Writer) {
	fmt.Fprintln(w, "# Results from a Spacing Trial")
	fmt.Fprintln(w, "name count max_length unique acceptable"+
		" interleaved muts_in_sites total_sites total_singles")
}

type SpacingTrialResult struct {
	name         string // genome name
	count        int    // number of sites
	maxLength    int    // length of longest segment
	unique       bool   // unique sticky ends?
	acceptable   bool   // longest segment < 8kb and unique sticky?
	interleaved  bool   // BsaI interleaved with BsmBI?
	mutsInSites  int    // Number of silent muts in sites
	totalSites   int    // Total number of silently mutated sites
	totalSingles int    // Total number sites silently mutated with 1 mut
}

/*
func WriteHeadings(w io.Writer) {
	fmt.Fprintln(w, "name count max_length unique acceptable"+
		" interleaved muts_in_sites total_sites total_singles")
}
*/

func (r *SpacingTrialResult) Write(w io.Writer) {
	fmt.Fprintln(w, r.name,
		r.count,
		r.maxLength, r.unique, r.acceptable, r.interleaved,
		r.mutsInSites, r.totalSites, r.totalSites)
}

func SpacingTrials(genome *Genomes, nd *NucDistro,
	numTrials int, numMuts int, countSites bool,
	results chan interface{}) {
	good := 0

	count, maxLength, unique, interleaved := FindRestrictionMap(genome)
	fmt.Printf("Original: %d, %d, %t, %t\n", count,
		maxLength, unique, interleaved)

	reportProgress := func(n int) {
		fmt.Printf("Tested %d. Found %d/%d good mutants (%.2f%%)\n", n,
			good, n, float64(good*100)/float64(n))
	}

	for i := 0; i < numTrials; i++ {
		mutant := genome.Clone()
		MutateSilent(mutant, nd, numMuts)
		count, maxLength, unique, interleaved = FindRestrictionMap(mutant)

		acceptable := unique && maxLength < 8000
		if acceptable {
			/*
				fmt.Printf("Mutant %d: %d, %d, %t, %t\n", i, count,
					maxLength, unique, interleaved)
			*/
			good += 1
		}

		var sis SilentInSites
		if countSites {
			mutant.Combine(genome)
			sis = CountSilentInSites(mutant, RE_SITES, true)
		}

		results <- &SpacingTrialResult{genome.names[0],
			count, maxLength, unique, acceptable, interleaved,
			sis.totalMuts, sis.totalSites, sis.totalSites}

		if i%100 == 0 {
			reportProgress(i)
		}
	}

	reportProgress(numTrials)
}
