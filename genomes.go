package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
)

/*
	Represents a collection of aligned genomes (usually two) with one genome
	per row. The orfs "belong" to the first one in the set. We also use this
	for a single genome.
*/
type Genomes struct {
	nts  [][]byte
	orfs *Orfs
}

func NewGenomes(orfs *Orfs, numGenomes int) *Genomes {
	return &Genomes{make([][]byte, numGenomes), orfs}
}

/*
	Load genomes, which might be a fasta file containing a single genome, or
	one containing a few of them in an alignment. Be a bit careful when working
	with alignments since there may be '-' in there
*/
func LoadGenomes(fname string, orfsName string) *Genomes {
	ret := NewGenomes(LoadOrfs(orfsName), 0)

	fd, err := os.Open(fname)
	if err != nil {
		log.Fatal("Can't open file")
	}
	defer fd.Close()

	fp := bufio.NewReader(fd)

	currentRow := make([]byte, 0)
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

		if strings.HasPrefix(line, ">") {
			if len(currentRow) > 0 {
				ret.nts = append(ret.nts, currentRow)
				currentRow = make([]byte, 0)
			}
			continue
		}

		currentRow = append(currentRow, []byte(line)...)
	}
	ret.nts = append(ret.nts, currentRow)
	return ret
}

func (g *Genomes) Save(name, fname string, which int) error {
	fd, err := os.Create(fname)
	if err != nil {
		return err
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	fmt.Fprintf(fp, ">%s\n", name)

	nts := g.nts[which]
	ll := 60
	var i int
	for i = 0; i < len(nts)-ll; i += ll {
		fmt.Fprintf(fp, "%s\n", string(nts[i:i+ll]))
	}

	fmt.Fprintf(fp, "%s\n", string(nts[i:]))
	fp.Flush()
	return nil
}

func (g *Genomes) Clone() *Genomes {
	ret := NewGenomes(g.orfs, g.NumGenomes())
	for i := 0; i < len(g.nts); i++ {
		ret.nts[i] = make([]byte, len(g.nts[i]))
		copy(ret.nts[i], g.nts[i])
	}
	return ret
}

/*
	Assuming align is aligned with g, add it to g's own nts array, just doing a
	shallow copy
*/
func (g *Genomes) Combine(other *Genomes) {
	for i := 0; i < other.NumGenomes(); i++ {
		g.nts = append(g.nts, other.nts[i])
	}
}

/*
	These little functions make it a bit easier not to get confused about
	which dimension is which.
*/
func (g *Genomes) NumGenomes() int {
	return len(g.nts)
}

func (g *Genomes) Length() int {
	return len(g.nts[0])
}
