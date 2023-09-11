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
	per row
*/
type Genomes [][]byte

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
}
