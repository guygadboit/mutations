package main

import (
	"fmt"
	"math/rand"
	"time"
)

// Counts for each nucleotide in a genome
type NucDistro struct {
	nts   map[byte]int
	total int
}

func (nd *NucDistro) Count(g Genome) {
	for i := 0; i < len(g); i++ {
		nt := g[i]
		count, _ := nd.nts[nt]
		nd.nts[nt] = count + 1
	}

	for k := range nd.nts {
		nd.total += nd.nts[k]
	}
}

func NewNucDistro(g Genome) *NucDistro {
	ret := NucDistro{nts: make(map[byte]int)}
	ret.Count(g)
	return &ret
}

func (nd *NucDistro) Show() {
	for k := range nd.nts {
		fmt.Printf("%c: %d\n", k, nd.nts[k])
	}
	fmt.Printf("Total: %d\n", nd.total)
}

var randGenerator *rand.Rand

/*
	Pick a nucleotide randomly from the distribution represented by nd
*/
func (nd *NucDistro) Random() byte {
	if randGenerator == nil {
		randGenerator = rand.New(rand.NewSource(time.Now().UnixNano()))
	}

	r := randGenerator.Intn(nd.total)

	var k byte
	for k = range nd.nts {
		if r < nd.nts[k] {
			break
		}
	}
	return k
}
