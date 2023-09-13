package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
)

type Genome struct {
	nts  []byte
	orfs *Orfs
}

/*
	Represents a collection of aligned genomes (usually two) with one genome
	per row
*/
type Genomes struct {
	nts  [][]byte
	orfs *Orfs
}

/*
	Load genomes, which might be a fasta file containing a single genome, or
	one containing a few of them in an alignment. Be a bit careful when working
	with alignments since there may be '-' in there
*/
func LoadGenomes(fname string, orfsName string) *Genomes {
	ret := &Genomes{make([][]byte, 0), LoadOrfs(orfsName)}

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

func LoadGenome(fname string, orfsName string) *Genome {
	genomes := LoadGenomes(fname, orfsName)
	return &Genome{genomes.nts[0], genomes.orfs}
}

func (g *Genome) Save(name, fname string) error {
	fd, err := os.Create(fname)
	if err != nil {
		return err
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	fmt.Fprintf(fp, ">%s\n", name)

	ll := 60
	var i int
	for i = 0; i < len(g.nts)-ll; i += ll {
		fmt.Fprintf(fp, "%s\n", string(g.nts[i:i+ll]))
	}

	fmt.Fprintf(fp, "%s\n", string(g.nts[i:]))
	fp.Flush()
	return nil
}

func (g *Genome) Clone() *Genome {
	ret := &Genome{make([]byte, len(g.nts)), g.orfs}
	copy(ret.nts, g.nts)
	return ret
}
