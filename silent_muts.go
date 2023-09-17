/*
	Explore various scoring systems for how suspiciously engineered a genome
	might look based on whether it has more silent mutations, or more of
	particular kinds, in the sites of interest rather than outside them.
*/
package main

import (
	"fmt"
)

type SilentInSites struct {
	totalMuts, totalSites, totalSingleSites int
}

func (sis SilentInSites) Show() {
	fmt.Printf("Total muts, total sites, total singles: %d %d %d\n",
		sis.totalMuts, sis.totalSites, sis.totalSingleSites)
}

/*
	Assume there are two aligned genomes in Genomes. Return if they code for
	the same protein for count nts at pos
*/
func isSilent(genomes *Genomes, pos int, count int) bool {
	var aEnv, bEnv Environment

	aEnv.Init(genomes, pos, count, 0)
	bEnv.Init(genomes, pos, count, 1)

	aProt := aEnv.Protein()
	bProt := bEnv.Protein()

	for k := 0; k < count; k++ {
		if aProt[k] != bProt[k] {
			return false
		}
	}
	return true
}

/*
	Return the total number of mutations in the sites, the number of sites, and
	the number of sites with only one mutation. If assumeSilent don't bother
	checking if the mutations were silent (because in our common use case we
	only introduce silent mutations in the first place, so this saves time)
*/
func CountSilentInSites(genomes *Genomes,
	sites []ReSite, assumeSilent bool) SilentInSites {
	var ret SilentInSites
	var s Search
	m := len(sites[0].pattern)

	for s.Init(genomes, sites); ; {
		pos, _ := s.Iter()
		if s.End() {
			break
		}

		// First count how many muts
		numMuts := 0
		for k := 0; k < m; k++ {
			if genomes.nts[0][pos+k] != genomes.nts[1][pos+k] {
				numMuts++
			}
		}

		// Site is the same in both genomes. Not interested.
		if numMuts == 0 {
			continue
		}

		// Now given that it is mutated, check that it's silent (we might
		// already know that it is)
		if !assumeSilent {
			if !isSilent(genomes, pos, m) {
				continue
			}
		}

		/*
			first := string(genomes.nts[0][pos:pos+m])
			second := string(genomes.nts[1][pos:pos+m])
			fmt.Printf("%s silent mut %s (%s %s) at %d\n", genomes.names[0],
				string(site.pattern),
				first, second, pos)
		*/

		ret.totalMuts += numMuts
		ret.totalSites++
		if numMuts == 1 {
			ret.totalSingleSites++
		}
	}
	return ret
}

/*
	For each of our alignments of WH1 with various relatives, count the silent
	in sites. We will compare these to the simulated figures
*/
func CountSilentInSitesReference(name string, sites []ReSite,
	results chan interface{}) {

	// WH1 is the first genome in each of the alignments, so we use its
	// ORFS
	baseName := fmt.Sprintf("WH1-%s", name)
	fname := fmt.Sprintf("%s.fasta", baseName)
	genomes := LoadGenomes(fname, "WH1.orfs")

	var result TamperTrialResult
	result.SilentInSites = CountSilentInSites(genomes, RE_SITES, false)
	result.name = baseName

	results <- &result
}
