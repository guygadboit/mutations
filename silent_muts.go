/*
	Explore various scoring systems for how suspiciously engineered a genome
	might look based on whether it has more silent mutations, or more of
	particular kinds, in the sites of interest rather than outside them.
*/
package main

type MutationScore interface {
	Score(genomes *Genomes) float64
}

/*
	Return the total number of mutations in the sites, the number of sites, and
	the number of sites with only one mutation.
*/
func countSites(genomes *Genomes, sites []ReSite) (int, int, int) {
	var totalMuts, totalSites, totalSingleSites int
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

		// Now given that it is mutated, check that it's silent
		var aEnv, bEnv Environment

		aEnv.Init(genomes, pos, m, 0)
		bEnv.Init(genomes, pos, m, 1)

		iProt := aEnv.Protein()
		jProt := bEnv.Protein()

		for k := 0; k < m; k++ {
			if iProt[k] != jProt[k] {
				continue // Not silent -> not interested
			}
		}

		totalMuts += numMuts
		totalSites++
		if numMuts == 1 {
			totalSingleSites++
		}
	}
	return totalMuts, totalSites, totalSingleSites
}
