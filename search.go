package main

import (
	"reflect"
)

type Search struct {
	genome  []byte
	reSites []ReSite
	i       int
}

func (s *Search) Init(genome []byte, reSites []ReSite) {
	s.genome = genome
	s.reSites = reSites
	s.i = 0
}

/*
	Returns the position in the genome and the ReSite that matched
*/
func (s *Search) Iter() (int, *ReSite) {
	n := len(s.genome)
	m := len(s.reSites[0].pattern)

	for ; s.i < n-m; s.i++ {
		for j := 0; j < len(s.reSites); j++ {
			site := &s.reSites[j]
			if reflect.DeepEqual(s.genome[s.i:s.i+m], site.pattern) {
				retVal := s.i
				s.i++
				return retVal, site
			}
		}
	}

	s.i = n
	return n, nil
}

func (s *Search) End() bool {
	return s.i == len(s.genome)
}
