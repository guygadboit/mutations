package main

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"log"
	"os"
	"reflect"
	"sort"
	"strconv"
	"strings"
)

var CodonTable = map[string]byte{
	"TTT": 'F', // Phenylalanine
	"TTC": 'F',

	"TTA": 'L', // Leucine
	"TTG": 'L',
	"CTT": 'L',
	"CTC": 'L',
	"CTA": 'L',
	"CTG": 'L',

	"ATT": 'I', // Isoleucine
	"ATC": 'I',
	"ATA": 'I',

	"ATG": 'M', // Methionine

	"GTT": 'V', // Valine
	"GTC": 'V',
	"GTA": 'V',
	"GTG": 'V',

	"TCT": 'S', // Serine
	"TCC": 'S',
	"TCA": 'S',
	"TCG": 'S',

	"CCT": 'P', // Proline
	"CCC": 'P',
	"CCA": 'P',
	"CCG": 'P',

	"ACT": 'T', // Threonine
	"ACC": 'T',
	"ACA": 'T',
	"ACG": 'T',

	"GCT": 'A', // Alanine
	"GCC": 'A',
	"GCA": 'A',
	"GCG": 'A',

	"TAT": 'Y', // Tyrosine
	"TAC": 'Y',

	"TAA": '*', // Stop
	"TAG": '*',

	"CAT": 'H', // Histidine
	"CAC": 'H',

	"CAA": 'Q', // Glutadine
	"CAG": 'Q',

	"AAT": 'N', // Asparagine
	"AAC": 'N',

	"AAA": 'K', // Lysine
	"AAG": 'K',

	"GAT": 'D', // Aspartic acid
	"GAC": 'D',

	"GAA": 'E', // Glutamic acid
	"GAG": 'E',

	"TGT": 'C', // Cysteine
	"TGC": 'C',

	"TGA": '*', // Stop
	"TGG": 'W', // Tryptophan

	"CGT": 'R', // Arginine
	"CGC": 'R',
	"CGA": 'R',
	"CGG": 'R',

	"AGT": 'S', // Serine
	"AGC": 'S',

	"AGA": 'R', // Arginine (again)
	"AGG": 'R',

	"GGT": 'G', // Glycine
	"GGC": 'G',
	"GGA": 'G',
	"GGG": 'G',
}

var ReverseCodonTable map[byte][]string

type Orf struct {
	start, end int
}

type Orfs []Orf

func LoadOrfs(fname string) *Orfs {
	ret := make(Orfs, 0)

	fd, err := os.Open(fname)
	if err != nil {
		log.Fatal("Can't open file")
	}
	defer fd.Close()

	fp := bufio.NewReader(fd)

loop:
	for {
		line, err := fp.ReadString('\n')
		switch err {
		case io.EOF:
			break loop
		case nil:
			break
		default:
			log.Fatal("Can't read file")
		}

		line = strings.TrimSpace(line)
		fields := strings.Fields(line)

		start, err := strconv.Atoi(fields[0])
		if err != nil {
			log.Fatal("Parse error in ORFs")
		}

		// ORFs seem to be conventionally 1-based
		start -= 1

		end, err := strconv.Atoi(fields[1])
		if err != nil {
			log.Fatal("Parse error in ORFs")
		}

		ret = append(ret, Orf{start, end})
	}

	return &ret
}

/*
	Return the start of the codon and where pos is in it. We just do a linear
	search since there aren't usually that many ORFs and this is probably as
	fast as anything else
*/
func (orfs Orfs) GetCodonOffset(pos int) (int, int, error) {
	for i := 0; i < len(orfs); i++ {
		orf := &orfs[i]
		if pos >= orf.start && pos < orf.end {
			orfPos := pos - orf.start // pos relative to start of ORF
			return orf.start + (orfPos/3)*3, orfPos % 3, nil
		}
	}
	return 0, 0, errors.New("Not in ORF")

}

// The "Environment" of a subsequence is the codon-aligned section that
// completely contains it.
type Environment struct {
	start int // Index into the original genome
	len   int // How many nts in the subsequence this represents

	window  []byte // The whole aligned section
	offset  int    // The offset to the start of the subsequence
	protein []byte // Its translation
}

// rounded up to the nearest multiple of 3
func ceil3(n int) int {
	return n + 3 - n%3
}

/*
	Assume nts are codon aligned and return a translation, with one amino-acid
	letter per nt, so something like LLLRRRIII
*/
func TranslateAligned(nts []byte) []byte {
	ret := make([]byte, len(nts))

	for i := 0; i < len(nts); i += 3 {
		aa := CodonTable[string(nts[i:i+3])]
		for j := 0; j < 3; j++ {
			ret[i+j] = aa
		}
	}
	return ret
}

func (env *Environment) Init(genome *Genomes,
	pos int, n int, which int) error {
	env.start = pos
	env.len = n

	windowStart, codonOffset, err := genome.orfs.GetCodonOffset(pos)
	if err != nil {
		return err
	}

	windowLen := ceil3(codonOffset + n)
	windowEnd := windowStart + windowLen

	env.offset = codonOffset
	env.window = genome.nts[which][windowStart:windowEnd]
	env.protein = TranslateAligned(env.window)
	return nil
}

func (env *Environment) Subsequence() []byte {
	return env.window[env.offset : env.offset+env.len]
}

func (env *Environment) Protein() []byte {
	return env.protein[env.offset : env.offset+env.len]
}

func (env *Environment) Print() {
	fmt.Println(string(env.Subsequence()))
	fmt.Println(string(env.Protein()))
}

/*
	If we were to replace the subsequence this is the environment of, would
	that be silent, and how many mutations would it contain?
*/
func (env *Environment) Replace(replacement []byte) (bool, int) {
	altWindow := make([]byte, len(env.window))
	copy(altWindow, env.window)
	copy(altWindow[env.offset:env.offset+env.len], replacement)

	protein := env.protein
	altProtein := TranslateAligned(altWindow)

	silent := true
	for i := 0; i < len(protein); i += 3 {
		if altProtein[i] != protein[i] {
			silent = false
			break
		}
	}

	subseq := env.Subsequence()
	differences := 0
	for i := 0; i < env.len; i++ {
		if replacement[i] != subseq[i] {
			differences++
		}
	}

	return silent, differences
}

// Silent alternative to a sequence of nts and how many muts that would require
type Alternative struct {
	numMuts int
	nts     []byte
}

type Alternatives []Alternative

// Iterator for finding the alternatives to a given subsequence
type altIter struct {
	protein  []byte // the protein we're finding nts for, as RL not RRRLLL
	odometer []int  // tracks the codon combinations as we iterate them
}

func (it *altIter) Init(protein []byte) {
	it.protein = protein
	it.odometer = make([]int, len(protein))
}

/*
	Return the next alternative and whether there are any more to come after
	it.
*/
func (it *altIter) Next() ([]byte, bool) {
	prot := it.protein
	ret := make([]byte, 0, len(prot)*3)

	for i := 0; i < len(prot); i++ {
		codons := ReverseCodonTable[prot[i]]
		ret = append(ret, []byte(codons[it.odometer[i]])...)
	}

	// Increment the odometer like a sort of odometer
	for j := 0; j < len(prot); j++ {
		codons := ReverseCodonTable[prot[j]]
		if it.odometer[j]+1 < len(codons) {
			it.odometer[j]++
			for k := 0; k < j; k++ {
				it.odometer[k] = 0
			}
			return ret, true
		}
	}

	return ret, false
}

func TestAlternatives() {
	protein := []byte("LF")
	var it altIter
	it.Init(protein)

	for {
		alt, more := it.Next()
		fmt.Println(string(alt))
		if !more {
			break
		}
	}
}

/*
	Newer versions of Go have a more "ergonomic" slices.SortFunc which saves
	you doing all this.
*/
func (a Alternatives) Len() int {
	return len(a)
}

func (a Alternatives) Less(i, j int) bool {
	return a[i].numMuts < a[j].numMuts
}

func (a Alternatives) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

/*
	Find the alternative nt sequences that would not change the protein here,
	ordered by fewest muts first.
*/
func (env *Environment) FindAlternatives(maxMuts int) Alternatives {
	var it altIter
	ret := make(Alternatives, 0)

	// The protein stored in env is like LLLRRRIII. We want just LRI.
	windowLen := len(env.window)
	protein := make([]byte, windowLen/3)
	for i := 0; i < len(protein); i++ {
		protein[i] = env.protein[i*3]
	}

	it.Init(protein)
	existing := env.Subsequence()

	for more := true; more; {
		var alt []byte
		alt, more = it.Next()
		start, end := env.offset, env.offset+env.len

		// The alternative is no good if it differs outside the subsequence
		if !reflect.DeepEqual(alt[:start], env.window[:start]) {
			continue
		}

		if !reflect.DeepEqual(alt[end:], env.window[end:]) {
			continue
		}

		numMuts := 0
		for i := 0; i < env.len; i++ {
			if alt[start+i] != existing[i] {
				numMuts++
			}
		}

		if numMuts > 0 && numMuts <= maxMuts {
			ret = append(ret, Alternative{numMuts,
				alt[start:end]})
		}

		if !more {
			break
		}
	}

	sort.Sort(ret)
	return ret
}

func init() {
	ReverseCodonTable = make(map[byte][]string)
	for k := range CodonTable {
		v, _ := CodonTable[k]

		codons, there := ReverseCodonTable[v]
		if !there {
			codons = make([]string, 0)
		}
		ReverseCodonTable[v] = append(codons, k)
	}
}
