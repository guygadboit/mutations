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
				fmt.Println("New row", len(ret))
				currentRow = make([]byte, 0)
			}
			continue
		}

		currentRow = append(currentRow, []byte(line)...)
	}
	ret = append(ret, currentRow)
	return ret
}

func main() {
	genomes := LoadGenomes("./WH1-RaT.fasta")
	orfs := LoadOrfs("./WH1.orfs")

	var env Environment
	env.Init(genomes[0], orfs, 716, 6)
	env.Print()

	env.Init(genomes[0], orfs, 718, 1)
	env.Print()

	var patterns = [][]byte{
		[]byte("CGTCTC"),
		[]byte("GAGACC"),
		[]byte("GGTCTC"),
		[]byte("GAGACG"),
	}

	var s Search

	for s.Init(genomes[0], patterns); ; {
		pos := s.Iter()
		if s.End() {
			break
		}
		fmt.Println(pos)
	}

	genomes = LoadGenomes("./BANAL-20-52.fasta")
	fmt.Println(len(genomes), len(genomes[0]))
	orfs = LoadOrfs("./BANAL-20-52.orfs")
	b52 := genomes[0]

	env.Init(b52, orfs, 716, 6)
	env.Print()

	nd := NewNucDistro(b52)
	nd.Show()

	fmt.Printf("%c\n", nd.Random())
	fmt.Printf("%c\n", nd.Random())
}
