/*
	Explore various scoring systems for how suspiciously engineered a genome
	might look based on whether it has more silent mutations, or more of
	particular kinds, in the sites of interest rather than outside them.
*/
package main

type MutationScore interface{
	Score(genomes *Genomes)	float64
}

/*
	Return the total number of mutations in the sites, the number of sites, and
	the number of sites with only one mutation.
*/
func countSites(genomes *Genomes, sites []ReSite) (int, int, int) {
	var totalMuts, totalSites, totalSingleSites int
	var s Search
	m := len(sites[0].pattern)
	individualGenomes := genomes.Split()

	other := func(which int) int {
		return (which + 1) % 2
	}

	for i := 0; i < len(genomes.nts); i++ {
		for s.Init(genomes.nts[i], sites); ; {
			pos, _ := s.Iter()
			if s.End() {
				break
			}

			j := other(i)

			// First count how many muts
			numMuts := 0
			for k := 0; k < m; k++ {
				if genomes.nts[i][pos+k] != genomes.nts[j][pos+k] {
					numMuts++
				}
			}

			// Site is the same in both genomes. Not interested.
			if numMuts == 0 {
				continue
			}

			// Now given that it is mutated, check that it's silent
			var iEnv, jEnv Environment

			iEnv.Init(individualGenomes[i], pos, m)
			jEnv.Init(individualGenomes[j], pos, m)

			iProt := iEnv.Protein()
			jProt := iEnv.Protein()

			for k := 0; k < m; k++ {
				if iProt[k] != jProt[k] {
					continue	// Not silent -> not interested
				}
			}

			totalMuts += numMuts
			totalSites++
			if numMuts == 1 {
				totalSingleSites++
			}
		}
	}
	return totalMuts, totalSites, totalSingleSites
}
