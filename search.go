package main

import (
	"reflect"
)

type Search struct {
	genomes *Genomes // Where we're looking
	reSites []ReSite // What we're looking for

	i     int         // Where we got to looking for sites
	cache SearchCache // Cached results if we did this before
}

type CachedSearch struct {
	Search
	cache     SearchCache
	searchI   int  // Where we are in the cache
	cacheFull bool // Whether it's full of valid data yet
}

type SearchCacheResult struct {
	pos  int
	site *ReSite
}

type SearchCache []SearchCacheResult

func (s *Search) Init(genomes *Genomes, reSites []ReSite) {
	s.genomes = genomes
	s.reSites = reSites
	s.i = 0
}

func (s *CachedSearch) Init(genomes *Genomes, reSites []ReSite) {
	// If we have cached data, we just need to rewind so that we start
	// retrieving it.
	if s.cacheFull {
		// This iterator weirdly starts at -1. We increment it before returning
		// results so that End() behaves in the same way as when we aren't
		// cached (when we need to search the last segment of the genome and
		// maybe won't find anything).
		s.searchI = -1
		return
	}

	s.Search.Init(genomes, reSites)
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

func (s *Search) GetSites() []ReSite {
	return s.reSites
}

func (s *CachedSearch) cacheIter() (int, *ReSite) {
	s.searchI++
	if s.searchI < len(s.cache) {
		cachedVal := s.cache[s.searchI]
		pos, site := cachedVal.pos, cachedVal.site
		return pos, site
	}
	return len(s.cache), nil
}

func (s *CachedSearch) Iter() (int, *ReSite) {
	if s.cacheFull {
		return s.cacheIter()
	}

	if s.cache == nil {
		s.cache = make(SearchCache, 0)
	}

	pos, site := s.Search.Iter()
	if site != nil {
		s.cache = append(s.cache, SearchCacheResult{pos, site})
	}

	// If we found all the results the cache is full. But set searchI to the
	// end, because we shouldn't retrieve these results until the Search is
	// reinited.
	if s.Search.End() {
		s.cacheFull = true
		s.searchI = len(s.cache)
	}

	return pos, site
}

func (s *CachedSearch) End() bool {
	if s.cacheFull {
		return s.searchI == len(s.cache)
	}
	return false
}
