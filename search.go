package main

import (
	"reflect"
)

type Search struct {
	genome   []byte
	patterns [][]byte
	i        int
}

func (s *Search) Init(genome []byte, patterns [][]byte) {
	s.genome = genome
	s.patterns = patterns
	s.i = 0
}

func (s *Search) Iter() int {
	n := len(s.genome)
	m := len(s.patterns[0])

	for ; s.i < n-m; s.i++ {
		for j := 0; j < len(s.patterns); j++ {
			if reflect.DeepEqual(s.genome[s.i:s.i+m], s.patterns[j]) {
				ret := s.i
				s.i++
				return ret
			}
		}
	}

	s.i = n
	return n
}

func (s *Search) End() bool {
	return s.i == len(s.genome)
}
