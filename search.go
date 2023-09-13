package main

import (
	"reflect"
)

type Search struct {
	nts     []byte
	reSites []ReSite
	i       int
}

func (s *Search) Init(nts []byte, reSites []ReSite) {
	s.nts = nts
	s.reSites = reSites
	s.i = 0
}

/*
	Returns the position in the nts and the ReSite that matched
*/
func (s *Search) Iter() (int, *ReSite) {
	n := len(s.nts)
	m := len(s.reSites[0].pattern)

	for ; s.i < n-m; s.i++ {
		for j := 0; j < len(s.reSites); j++ {
			site := &s.reSites[j]
			if reflect.DeepEqual(s.nts[s.i:s.i+m], site.pattern) {
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
	return s.i == len(s.nts)
}
