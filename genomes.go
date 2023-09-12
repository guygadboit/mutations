package main

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
)

type Genome []byte

/*
	Represents a collection of aligned genomes (usually two) with one genome
	per row
*/
type Genomes []Genome

/*
	Load genomes, which might be a fasta file containing a single genome, or
	one containing a few of them in an alignment. Be a bit careful when working
	with alignments since there may be '-' in there
*/
func LoadGenomes(fname string) Genomes {
	ret := make(Genomes, 0, 2)

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
				ret = append(ret, currentRow)
				currentRow = make([]byte, 0)
			}
			continue
		}

		currentRow = append(currentRow, []byte(line)...)
	}
	ret = append(ret, currentRow)
	return ret
}

func (g Genome) Save(name, fname string) error {
	fd, err := os.Create(fname)
	if err != nil {
		return err
	}
	defer fd.Close()

	fp := bufio.NewWriter(fd)
	fmt.Fprintf(fp, ">%s\n", name)

	var i int
	for i = 0; i < len(g)-60; i += 60 {
		fmt.Fprintf(fp, "%s\n", string(g[i:i+60]))
	}

	fmt.Fprintf(fp, "%s\n", string(g[i:]))
	fp.Flush()
	return nil
}

func main() {
	Trials("BANAL-20-52.fasta", "BANAL-20-52.orfs", 100)
}
