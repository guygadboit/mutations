package main

import (
	"errors"
	"math/rand"
)

/*
	Add one of the sites in sites somewhere randomly but not in notAt using a
	maximum of maxMuts mutations. Return where it was added or an error in the
	unlikely event that it couldn't be.
*/
func AddSite(genome *Genomes, sites []ReSite,
	notAt map[int]bool, maxMuts int) (int, error) {
	site := sites[rand.Intn(len(sites))]
	m := len(site.pattern)

	var tryAdd = func(pos int) bool {
		_, there := notAt[pos]
		if there {
			return false
		}

		var env Environment
		err := env.Init(genome, pos, m, 0)
		if err != nil {
			return false
		}

		silent, numMuts := env.Replace(site.pattern)
		if silent && numMuts <= maxMuts {
			copy(genome.nts[0][pos:pos+m], site.pattern)
		}
		return silent && numMuts <= maxMuts
	}

	start := rand.Intn(genome.Length())
	for i := start; i < genome.Length(); i++ {
		if tryAdd(i) {
			return i, nil
		}
	}

	for i := 0; i < start; i++ {
		if tryAdd(i) {
			return i, nil
		}
	}
	return 0, errors.New("Can't add site")
}

/*
	Remove a site from somewhere random, but not in notAt. Return the position
	it was removed from.
*/
func RemoveSite(genome *Genomes,
	search *CachedSearch, notAt map[int]bool) (int, error) {
	n := genome.Length()
	sites := search.GetSites()
	m := len(sites[0].pattern)
	nts := genome.nts[0]

	genomeStart := rand.Intn(n)

	var tryRemove = func(pos int) bool {
		_, there := notAt[pos]
		if there {
			return false
		}

		var env Environment
		err := env.Init(genome, pos, m, 0)
		if err != nil {
			return false
		}
		alternatives := env.FindAlternatives(1)

		if len(alternatives) == 0 {
			return false
		}

		alt := alternatives[rand.Intn(len(alternatives))]

		/*
			fmt.Printf("Replacing %s <- %s at %d\n",
				string(nts[pos:pos+m]), string(alt.nts), pos)
		*/

		copy(nts[pos:pos+m], alt.nts)
		return true
	}

	// First look for sites after our random starting point
	for search.Init(genome, sites); ; {
		pos, site := search.Iter()
		if site == nil {
			break
		}

		if pos >= genomeStart {
			if tryRemove(pos) {
				return pos, nil
			}
		}
	}

	// If we didn't find any, search before the random starting point
	for search.Init(genome, sites); ; {
		pos, site := search.Iter()
		if site == nil {
			break
		}

		if pos < genomeStart {
			if tryRemove(pos) {
				return pos, nil
			}
		}
	}

	return 0, errors.New("Can't find a site to remove")
}

/*
	Try to silently remove the specified numbers of sites. Return the actual number
	modified
*/
func Tamper(genome *Genomes, sites []ReSite, remove, add int) int {
	removed := make(map[int]bool)
	count := 0

	var search CachedSearch
	search.Init(genome, sites)

	for i := 0; i < remove; i++ {
		pos, err := RemoveSite(genome, &search, removed)
		if err == nil {
			removed[pos] = true
			count++
		} else {
			break
		}
	}

	for i := 0; i < add; i++ {
		_, err := AddSite(genome, search.GetSites(), removed, 1)
		if err == nil {
			count++
		} else {
			break
		}
	}

	return count
}
