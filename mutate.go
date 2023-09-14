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
