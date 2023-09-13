package main

import (
	"reflect"
)

type Search struct {
	genomes *Genomes
	reSites []ReSite
	i       int
}

func (s *Search) Init(genomes *Genomes, reSites []ReSite) {
	s.genomes = genomes
	s.reSites = reSites
	s.i = 0
}

/*
	Returns the position in the nts and the ReSite that matched. If there are
	multiple genomes we look for matches in any of them.
*/
func (s *Search) Iter() (int, *ReSite) {
	genomes := s.genomes
	nts := s.genomes.nts
	n := genomes.Length()
	m := len(s.reSites[0].pattern)

	for ; s.i < n-m; s.i++ {
		for j := 0; j < len(s.reSites); j++ {
			site := &s.reSites[j]
			for k := 0; k < genomes.NumGenomes(); k++ {
				if reflect.DeepEqual(nts[k][s.i:s.i+m], site.pattern) {
					retVal := s.i
					s.i++
					return retVal, site
				}
			}
		}
	}

	s.i = n
	return n, nil
}

func (s *Search) End() bool {
	return s.i == s.genomes.Length()
}
