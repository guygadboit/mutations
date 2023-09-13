package main

import (
	"errors"
	//"fmt"
)

type ReSite struct {
	pattern     []byte
	stickyStart int
	stickyEnd   int  // I mean the end of the sticky end
	reverse     bool // Whether the sticky end needs to be reversed
}

var RE_SITES = []ReSite{
	{[]byte("CGTCTC"), 7, 11, false},
	{[]byte("GAGACG"), -5, -1, true},
	{[]byte("GGTCTC"), 7, 11, false},
	{[]byte("GAGACC"), -5, -1, true},
}

func reverse(b []byte) []byte {
	n := len(b)
	ret := make([]byte, n)
	for i := 0; i < n; i++ {
		ret[i] = b[n-i-1]
	}
	return ret
}

func getStickyEnd(genome *Genomes, pos int, site *ReSite) (string, error) {
	start := pos + site.stickyStart
	end := pos + site.stickyEnd

	if start < 0 || end > genome.Length() {
		return "", errors.New("Out of bounds")
	}

	s := genome.nts[0][start:end]
	if site.reverse {
		s = reverse(s)
	}

	/*
		fmt.Printf("pos: %d %s sticky %d->%d reverse %t: %s\n", pos,
		site.pattern, start, end, site.reverse, string(s))
	*/

	return string(s), nil
}

/*
	Returns the number of segments the length of the longest one, and whether
	the sticky ends are all unique
*/
func FindRestrictionMap(genome *Genomes) (int, int, bool) {
	var s Search
	prev, maxLength, count := 0, 0, 0
	stickyEnds := make(map[string]int)
	unique := true

	for s.Init(genome, RE_SITES); ; {
		pos, site := s.Iter()
		if s.End() {
			break
		}

		stickyEnd, err := getStickyEnd(genome, pos, site)
		if err == nil {
			n, _ := stickyEnds[stickyEnd]
			if n > 0 {
				unique = false
			}
			stickyEnds[stickyEnd] = n + 1
		}

		count++
		length := pos - prev
		if length > maxLength {
			maxLength = length
		}
		prev = pos
	}

	// One more segment from last position found to the end
	length := genome.Length() - prev
	if length > maxLength {
		maxLength = length
	}
	count++

	return count, maxLength, unique
}
