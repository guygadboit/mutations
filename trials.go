package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"sync"
)

type TrialResult struct {
	name         string // genome name
	count        int    // number of sites
	maxLength    int    // length of longest segment
	unique       bool   // unique sticky ends?
	acceptable   bool   // longest segment < 8kb and unique sticky?
	interleaved  bool   // BsaI interleaved with BsmBI?
	mutsInSites  int    // Number of silent muts in sites
	totalSites   int    // Total number of silently mutated sites
	totalSingles int    // Total number sites silently mutated with 1 mut
}

func WriteHeadings(w io.Writer) {
	fmt.Fprintln(w, "name count max_length unique acceptable"+
		" interleaved muts_in_sites total_sites total_singles")
}

func (r *TrialResult) Write(w io.Writer) {
	fmt.Fprintln(w, r.name,
		r.count,
		r.maxLength, r.unique, r.acceptable, r.interleaved,
		r.mutsInSites, r.totalSites, r.totalSites)
}

func Trials(genome *Genomes, nd *NucDistro,
	numTrials int, numMuts int, countSites bool,
	results chan *TrialResult) {
	good := 0

	count, maxLength, unique, interleaved := FindRestrictionMap(genome)
	fmt.Printf("Original: %d, %d, %t, %t\n", count,
		maxLength, unique, interleaved)

	reportProgress := func(n int) {
		fmt.Printf("Tested %d. Found %d/%d good mutants (%.2f%%)\n", n,
			good, n, float64(good*100)/float64(n))
	}

	for i := 0; i < numTrials; i++ {
		mutant := genome.Clone()
		MutateSilent(mutant, nd, numMuts)
		count, maxLength, unique, interleaved = FindRestrictionMap(mutant)

		mutant.Combine(genome)

		acceptable := unique && maxLength < 8000
		if acceptable {
			/*
				fmt.Printf("Mutant %d: %d, %d, %t, %t\n", i, count,
					maxLength, unique, interleaved)
			*/
			good += 1
		}

		var sis SilentInSites
		if countSites {
			sis = CountSilentInSites(mutant, RE_SITES, true)
		}

		results <- &TrialResult{genome.names[0],
			count, maxLength, unique, acceptable, interleaved,
			sis.totalMuts, sis.totalSites, sis.totalSites}

		if i%100 == 0 {
			reportProgress(i)
		}
	}

	reportProgress(numTrials)
}

func main() {
	var nTrials, nMuts, nThreads int
	var test, countSites bool

	flag.IntVar(&nTrials, "n", 10000, "Number of trials")
	flag.IntVar(&nMuts, "m", 763, "Number of mutations per mutant")
	flag.IntVar(&nThreads, "p", 1, "Number of threads")
	flag.BoolVar(&test, "t", false, "Just do some self-tests")
	flag.BoolVar(&countSites, "c", false, "Count mutations per site etc.")
	flag.Parse()

	if test {
		Test()
		return
	}

	nd := NewNucDistro(nil)

	fnames := []string{
		"ChimericAncestor",
		"BtSY2",
		"BANAL-20-236",
		"BANAL-20-52",
		"BANAL-20-103",
		"RaTG13",
	}

	genomes := make([]*Genomes, len(fnames))
	for i := 0; i < len(fnames); i++ {
		genomes[i] = LoadGenomes(
			fnames[i]+".fasta",
			fnames[i]+".orfs",
		)
	}

	for i := 0; i < len(genomes); i++ {
		nd.Count(genomes[i])
	}
	nd.Show()

	fd, err := os.Create("results.txt")
	if err != nil {
		log.Fatal("Can't create results file")
	}
	defer fd.Close()

	resultsWriter := bufio.NewWriter(fd)
	WriteHeadings(resultsWriter)
	results := make(chan *TrialResult, 1000)
	var wg sync.WaitGroup

	// Cut the work up unto nThreads pieces, all writing their results to a
	// single channel. Each thread will do a portion of the tests but for all
	// genomes
	for i := 0; i < nThreads; i++ {
		wg.Add(len(genomes))

		go func() {
			for j := 0; j < len(genomes); j++ {
				Trials(genomes[j], nd, nTrials/nThreads,
					nMuts, countSites, results)
				wg.Done()
			}
		}()
	}

	// Keep reading out of the results channel and writing to the results file
	// until everyone has finished (we will write to stop once wg has
	// completed)
	stop := make(chan bool)
	go func(stop chan bool) {
	loop:
		for {
			select {
			case r := <-results:
				r.Write(resultsWriter)
			case <-stop:
				break loop
			}
		}
		resultsWriter.Flush()
		fmt.Println("Wrote results.txt")
	}(stop)

	wg.Wait()
	stop <- true
}
