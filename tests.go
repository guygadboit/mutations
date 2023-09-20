package main

import (
	"fmt"
)

func testMutations(genome *Genomes) {
	fmt.Printf("Loaded %d genomes length %d\n",
		genome.NumGenomes(), genome.Length())

	genome.Save("B52", "B52-test.fasta", 0)
	fmt.Printf("Saved as B52-test.fasta\n")

	nd := NewNucDistro(genome)

	var mutant *Genomes
	for {
		mutant = genome.Clone()
		MutateSilent(mutant, nd, 700)
		count, maxLength, unique, interleaved := FindRestrictionMap(mutant)
		if unique && maxLength < 8000 {
			fmt.Println(count, maxLength, unique, interleaved)
			break
		}
	}
	mutant.Save("Mutant", "B52-mutated.fasta", 0)
	fmt.Printf("Saved as B52-mutated.fasta\n")
}

func testTamper(genome *Genomes) {
	num := Tamper(genome, RE_SITES, 10, 10)
	fmt.Printf("Tampered with %d sites\n", num)

	genome.Save("Mutant", "B52-mutated.fasta", 0)
	fmt.Printf("Saved as B52-mutated.fasta\n")
}

func testAlternatives(genome *Genomes) {
	var env Environment
	env.Init(genome, 300, 6, 0)

	alternatives := env.FindAlternatives(6)
	fmt.Println(alternatives)
}

func testCachedSearch(genome *Genomes) {
	var cs CachedSearch

	for i := 0; i < 3; i++ {
		fmt.Printf("Starting search\n")
		cs.Init(genome, RE_SITES)
		for {
			pos, site := cs.Iter()
			if cs.End() {
				break
			}

			fmt.Printf("%s at %d\n", string(site.pattern), pos)
		}
	}
}

func testTranslate(genome *Genomes) {
	var it CodonIter
	it.Init(genome, 0)

	for {
		pos, codon, aa, err := it.Next()
		if err != nil {
			break
		}
		fmt.Printf("%d: %s %c\n", pos, codon, aa)
	}
	fmt.Printf("\n")
}

func Test() {
	genome := LoadGenomes("BANAL-20-52.fasta", "BANAL-20-52.orfs")
	// testCachedSearch(genome)
	// testMutations(genome)
	// testAlternatives(genome)
	// testTamper(genome)
	testTranslate(genome)
}
