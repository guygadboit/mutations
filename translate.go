package main

import (
	"bufio"
	//"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
	"errors"
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

type Orf struct {
	start, end int
}

type Orfs []Orf

func LoadOrfs(fname string) Orfs {
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

		end, err := strconv.Atoi(fields[0])
		if err != nil {
			log.Fatal("Parse error in ORFs")
		}

		ret = append(ret, Orf{start, end})
	}

	return ret
}

// Return the start of the codon and where pos is in it
func (orfs Orfs) GetCodonOffset(pos int) (int, int, error) {
	for i := 0; i < len(orfs); i++ {
		orf := &orfs[i]
		if pos >= orf.start && pos < orf.end {
			orf_rel := pos - orf.start
			return (orf_rel / 3) * 3, orf_rel % 3, nil
		}
	}
	return 0, 0, errors.New("Not in ORF")

}

// The "Environment" of a subsequence is the codon-aligned section that
// completely contains it.
type Environment struct {
	start, end	int			// Indices into the original genome
	env			[]byte		// The whole aligned section
	offset		int			// The offset to the start of the subsequence
	protein		[]byte		// Its translation 
}

func (env *Environment) Init(genome []byte, start int, n int) {
	env.start = start
	env.end = start + n



}
