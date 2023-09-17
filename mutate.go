package main

import (
	"math/rand"
)

/*
	Introduce num silent mutations into genome (the first one), selecting nts
	randomly from nucDist. Return the number of mutations
*/
func MutateSilent(genome *Genomes, nucDist *NucDistro, num int) int {
	numMuts := 0
	alreadyDone := make(map[int]int)
	nts := genome.nts[0]

	// Try to mutate silently at pos. Return true if we succeeded.
	tryMutate := func(pos int) bool {
		done, _ := alreadyDone[pos]
		if done != 0 {
			return false
		}

		var env Environment
		err := env.Init(genome, pos, 1, 0)
		if err != nil {
			// You get an error if pos is not in an ORF
			return false
		}

		existing := nts[pos]
		var replacement byte
		for {
			replacement = nucDist.Random()
			if replacement != existing {
				break
			}
		}

		silent, _ := env.Replace([]byte{replacement})
		if silent {
			nts[pos] = replacement
			alreadyDone[pos] = 1
			numMuts++
		}
		return silent
	}

mutations:
	for i := 0; i < num; {
		start := rand.Intn(genome.Length())

		for j := start; j < genome.Length(); j++ {
			if tryMutate(j) {
				i++
				continue mutations
			}
		}

		for j := 0; j < start; j++ {
			if tryMutate(j) {
				i++
				continue mutations
			}
		}

		// If we get here it means we've wrapped around the whole genome and
		// not found anywhere to put another silent mutation in. This shouldn't
		// ever happen.
		break
	}
	return numMuts
}

/*
	Returns the number of silent and non-silent mutations in an alignment of
	two genomes. Ignores indels.
*/
func CountMutations(genomes *Genomes) (int, int) {
	var nonSilent, silent int
	var env Environment
	a_nts := genomes.nts[0]
	b_nts := genomes.nts[1]
	n := genomes.Length()

	for i := 0; i < n; i++ {
		a := a_nts[i]
		b := b_nts[i]

		if a == b {
			continue
		}

		if a == '-' || b == '-' {
			continue
		}

		err := env.Init(genomes, i, 1, 0)
		if err != nil {
			// Ignore anything not in an ORF
			continue
		}
		isSilent, _ := env.Replace(b_nts[i : i+1])

		if isSilent {
			silent++
		} else {
			nonSilent++
		}
	}
	return silent, nonSilent
}
