package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"io"
)

func WriteHeadings(results io.Writer) {
	fmt.Fprintln(results, "name count max_length unique acceptable"+
		" muts_in_sites total_sites total_singles")
}

func Trials(genome *Genomes, nd *NucDistro,
	numTrials int, numMuts int, countSites bool, results io.Writer) int {
	good := 0

	count, maxLength, unique := FindRestrictionMap(genome)
	fmt.Printf("Original: %d, %d, %t\n", count, maxLength, unique)

	reportProgress := func(n int) {
		fmt.Printf("Tested %d. Found %d/%d good mutants (%.2f%%)\n", n,
			good, n, float64(good*100)/float64(n))
	}

	for i := 0; i < numTrials; i++ {
		mutant := genome.Clone()
		MutateSilent(mutant, nd, numMuts)
		count, maxLength, unique = FindRestrictionMap(mutant)

		mutant.Combine(genome)

		acceptable := unique && maxLength < 8000
		if acceptable {
			fmt.Printf("Mutant %d: %d, %d, %t\n", i, count, maxLength, unique)
			good += 1
		}

		var sis SilentInSites
		if countSites {
			sis = CountSilentInSites(mutant, RE_SITES, true)
		}
		fmt.Fprintln(results, genome.names[0], count,
			maxLength, unique, acceptable,
			sis.totalMuts, sis.totalSites, sis.totalSites)

		if i%100 == 0 {
			reportProgress(i)
		}
	}

	reportProgress(numTrials)
	return good
}

func main() {
	var nTrials, nMuts int
	var test, countSites bool

	flag.IntVar(&nTrials, "n", 10000, "Number of trials")
	flag.IntVar(&nMuts, "m", 763, "Number of mutations per mutant")
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
		/*
		"BtSY2",
		"BANAL-20-236",
		"BANAL-20-52",
		"BANAL-20-103",
		"RaTG13",
		*/
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

	results := make([]int, len(genomes))

	fd, err := os.Create("results.txt")
	if err != nil {
		log.Fatal("Can't create results file")
	}
	defer fd.Close()

	resultsWriter := bufio.NewWriter(fd)
	WriteHeadings(resultsWriter)

	for i := 0; i < len(genomes); i++ {
		results[i] = Trials(genomes[i], nd, nTrials, nMuts,
			countSites, resultsWriter)
	}

	for i := 0; i < len(genomes); i++ {
		fmt.Printf("%s: %d/%d %.2f%%\n", genomes[i].names[0],
			results[i], nTrials, float64(100.0*results[i])/float64(nTrials))
	}

	resultsWriter.Flush()
	fmt.Println("Wrote results.txt")
}
