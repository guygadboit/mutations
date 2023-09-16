package main

import (
	"fmt"
	"io"
	"math/rand"
)

type TamperTrial struct {
	runFunc func(genome *Genomes, results chan interface{})
}

func (t *TamperTrial) Run(genome *Genomes, results chan interface{}) {
	t.runFunc(genome, results)
}

func (t *TamperTrial) WriteHeadings(w io.Writer) {
	fmt.Fprintln(w, "name tampered muts_in_sites total_sites total_singles")
}

type TamperTrialResult struct {
	SilentInSites
	name     string
	tampered bool
}

func (r *TamperTrialResult) Write(w io.Writer) {
	fmt.Fprintln(w, r.name, r.tampered,
		r.totalMuts, r.totalSites, r.totalSingleSites)
}

func TamperTrials(genome *Genomes, nd *NucDistro,
	numTrials int, numMuts int, results chan interface{}) {

	reportProgress := func(n int) {
		fmt.Printf("Run %d/%d trials\n", n, numTrials)
	}

	for i := 0; i < numTrials; i++ {
		mutant := genome.Clone()
		MutateSilent(mutant, nd, numMuts)

		tampered := rand.Intn(2) == 1
		if tampered {
			Tamper(mutant, RE_SITES, 3, 3)
		}

		var result TamperTrialResult
		mutant.Combine(genome)
		result.SilentInSites = CountSilentInSites(mutant, RE_SITES, true)
		result.name = genome.names[0]
		result.tampered = tampered

		results <- &result

		if i%100 == 0 {
			reportProgress(i)
		}
	}
}
