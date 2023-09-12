package main

type ReSite struct {
	pattern     []byte
	stickyStart int
	stickyEnd   int // I mean the end of the sticky end
}

var RE_SITES = []ReSite{
	{[]byte("CGTCTC"), 1, 5},
	{[]byte("GAGACG"), -11, -7},
	{[]byte("GGTCTC"), 1, 5},
	{[]byte("GAGACC"), -11, -7},
}

/*
	Returns the number of segments the length of the longest one, and whether
	the sticky ends are all unique
*/
func FindMap(genome Genome) (int, int, bool) {
	var s Search
	prev, maxLength, count := 0, 0, 0
	stickyEnds := make(map[string]int)
	unique := true

	for s.Init(genome, RE_SITES); ; {
		pos, site := s.Iter()
		if s.End() {
			break
		}

		stickyEnd := string(genome[pos+site.stickyStart : pos+site.stickyEnd])
		n, _ := stickyEnds[stickyEnd]
		if n > 0 {
			unique = false
		}
		stickyEnds[stickyEnd] = n + 1

		count++
		length := pos - prev
		if length > maxLength {
			maxLength = length
		}
		prev = pos
	}

	// One more segment from last position found to the end
	length := len(genome) - prev
	if length > maxLength {
		maxLength = length
	}
	count++

	return count, maxLength, unique
}
