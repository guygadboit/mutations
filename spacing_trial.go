package main

import (
	"fmt"
	"io"
	"strings"
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
		" interleaved muts_in_sites total_sites total_singles"+
		" num_muts added removed genome_len positions")
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
	numMuts		 int	// How many muts did we do
	added        int    // How many sites were added?
	removed      int    // How many sites were removed?
	genomeLen    int    // length of the whole genome
	positions    []int  // the actual positions of the sites
}

func (r *SpacingTrialResult) Write(w io.Writer) {

	strPositions := make([]string, len(r.positions))
	for i, pos := range r.positions {
		strPositions[i] = fmt.Sprintf("%d", pos)
	}
	positions := "[" + strings.Join(strPositions, ",") + "]"

	fmt.Fprintln(w, r.name, r.count,
		r.maxLength, r.unique, r.acceptable, r.interleaved,
		r.mutsInSites, r.totalSites, r.totalSingles,
		r.numMuts, r.added, r.removed, r.genomeLen, positions)
}

func toSet(a []int) map[int]bool {
	ret := make(map[int]bool)
	for _, v := range a {
		ret[v] = true
	}
	return ret
}

// Return the number of sites added and removed
func addedRemoved(before map[int]bool, after []int) (int, int) {
	var added, removed int

	for _, v := range after {
		if !before[v] {
			added++
		}
	}

	afterSet := toSet(after)
	for k, _ := range before {
		if !afterSet[k] {
			removed++
		}
	}
	return added, removed
}

func SpacingTrials(genome *Genomes, nd *NucDistro,
	numTrials int, numMuts int, countSites bool,
	results chan interface{}) {
	good := 0

	count, maxLength, unique, interleaved, positions :=
		FindRestrictionMap(genome)
	originalPositions := toSet(positions)

	fmt.Printf("Original: %d, %d, %t, %t\n", count,
		maxLength, unique, interleaved)

	reportProgress := func(n int) {
		fmt.Printf("Tested %d. Found %d/%d good mutants (%.2f%%)\n", n,
			good, n, float64(good*100)/float64(n))
	}

	for i := 0; i < numTrials; i++ {
		mutant := genome.Clone()
		MutateSilent(mutant, nd, numMuts)
		count, maxLength, unique, interleaved, positions =
			FindRestrictionMap(mutant)

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

		added, removed := addedRemoved(originalPositions, positions)

		results <- &SpacingTrialResult{genome.names[0],
			count, maxLength, unique, acceptable, interleaved,
			sis.totalMuts, sis.totalSites,
			sis.totalSites, numMuts, added, removed,
			genome.Length(), positions}

		if i%100 == 0 {
			reportProgress(i)
		}
	}

	reportProgress(numTrials)
}
