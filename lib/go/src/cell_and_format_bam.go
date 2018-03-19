package main

import "flag"
import "fmt"
import "os"
import "bufio"
import "github.com/brentp/bigly/bamat"
//import "strconv"
import "github.com/biogo/hts/bam"
//import "code.google.com/p/biogo.bam"
import "github.com/biogo/hts/sam"
import "log"
import "strings"
import "runtime/pprof"
import _ "net/http/pprof"

var bamin = flag.String("bam", "", "input bam from 10x single cell data")
var cells = flag.String("cells", "", "input csv file with cell barcodes as first column")
var chrom = flag.String("chrom", "", "chromosome of interest")
var start = flag.Int("start", 0, "start of fetch in bam")
var end = flag.Int("end", 10000000000, "end of fetch in bam")
var bamout = flag.String("bamout", "", "output bam filename")
var readGroupsOut = flag.String("readGroupsOut", "", "output readgroup file")
var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")

func main() {
    flag.Parse()
    //profiling code
    if *cpuprofile != "" {
        f, err := os.Create(*cpuprofile)
        if err != nil {
            log.Fatal("could not create CPU profile: ", err)
        }
        if err := pprof.StartCPUProfile(f); err != nil {
            log.Fatal("could not start CPU profile: ", err)
        }
        defer pprof.StopCPUProfile()
    }

    // make map of cell barcodes
    cellMap := GetCellBarcodes()
    readGroups := map[string]bool{}
    
    // open bam reader
    bamReader, bamReaderErr := bamat.New(*bamin)
    if bamReaderErr != nil { log.Fatal(bamReaderErr) } 
    defer bamReader.Close()
    bamIter, bamErr := bamReader.Query(*chrom, *start, *end)
    if bamErr != nil { log.Fatal(bamErr) }

    // open bam writer
    bamFile, bamFileErr := os.Open(*bamin)
    if bamFileErr != nil { log.Fatal(bamFileErr) }
    bamHeaderReader, bamHeaderReaderErr := bam.NewReader(bamFile, 1)
    if bamHeaderReaderErr != nil { log.Fatal(bamHeaderReaderErr) }
    bamFileOut, bamFileOutErr := os.Create(*bamout)
    if bamFileOutErr != nil { log.Fatal(bamFileOutErr) }
    bamWriter, bamWriterErr := bam.NewWriter(bamFileOut, bamHeaderReader.Header(), 1)
    if bamWriterErr != nil { log.Fatal(bamWriterErr) }
    defer bamWriter.Close()
        
    var err error
    var extant bool
    readGroupTag := sam.NewTag("RG")
    cellBxTag := sam.NewTag("CB")
    hasNs := false
    for {
        extant = bamIter.Next()
        if !extant {
            break
        }
        read := bamIter.Record()
        cellBarcode := read.AuxFields.Get(cellBxTag)
        if cellBarcode == nil {
            continue
        }
        _, extant = cellMap[cellBarcode.Value().(string)]
        if !extant {
            continue
        }
        readGroups[cellBarcode.Value().(string)] = true //add to list of read groups
        cigar := read.Cigar
        cigarOps := []sam.CigarOp(cigar)
        tagIndex := -1
        for index, tag := range read.AuxFields {
            if tag.Tag() == readGroupTag {
                tagIndex = index
                break
            }
        }
        //if tagIndex == -1 {
        //    fmt.Fprintln(os.Stderr, "RG tag didn't exist?")
        //    fmt.Fprintln(os.Stderr, read.String())
        //    os.Exit(-1)
        //}
        read.AuxFields[tagIndex], _ = sam.NewAux(readGroupTag, cellBarcode.Value()) //sam.ParseAux([]byte("RG:Z:"+cellBarcode.Value().(string)))
        
        hasNs = false
        for cigarOp := range cigarOps {
            if sam.CigarOp(cigarOp).Type() == sam.CigarSkipped {
                hasNs = true
            }
        } 
        if !hasNs {
            bamWriter.Write(read)
            //if err != nil { log.Fatal(err) }
        } else {
            seq := []byte(read.Seq.Expand())
            readOffsetStart := 0
            readOffsetEnd := 0
            refOffsetStart := read.Pos
            refOffsetEnd := read.Pos
            qual := make([]byte, len(read.Qual))
            copy(qual, read.Qual)
            cigarBuild := []sam.CigarOp{}
            for cigarOp := range cigarOps {
                cigarOpTyped := sam.CigarOp(cigarOp)
                cigarType := cigarOpTyped.Type()
                if cigarType == sam.CigarSkipped {
                    read.Seq = sam.NewSeq(seq[readOffsetStart:readOffsetEnd])
                    read.Qual = qual[readOffsetStart:readOffsetEnd]
                    read.Cigar = sam.Cigar(cigarBuild)
                    read.Pos = refOffsetStart
                    if !sam.IsValidRecord(read) {
                        fmt.Fprintln(os.Stderr, "record not valid")
                        fmt.Fprintln(os.Stderr, read.String())
                        os.Exit(-1)
                    }
                    err = bamWriter.Write(read)
                    if err != nil { log.Fatal(err) }
                    readOffsetStart = readOffsetEnd
                    refOffsetStart = refOffsetEnd
                    cigarBuild = []sam.CigarOp{}
                } else {
                    cigarBuild = append(cigarBuild, cigarOpTyped)
                }
                cigarLength := cigarOpTyped.Len()
                consumption := cigarType.Consumes()
                readOffsetEnd += consumption.Query * cigarLength
                refOffsetEnd += consumption.Reference * cigarLength
            }
            if len(cigarBuild) > 0 {
                read.Seq = sam.NewSeq(seq[readOffsetStart:readOffsetEnd])
                read.Qual = qual[readOffsetStart:readOffsetEnd]
                read.Cigar = sam.Cigar(cigarBuild)
                read.Pos = refOffsetStart
                if !sam.IsValidRecord(read) {
                    fmt.Fprintln(os.Stderr, "record not valid")
                    fmt.Fprintln(os.Stderr, read.String())
                    os.Exit(-1)
                }
                err = bamWriter.Write(read)
                if err != nil { log.Fatal(err) }
            }
        }
    }

    //output read groups
    readGroupWriter, readGroupWriterErr := os.Create(*readGroupsOut)
    if readGroupWriterErr != nil { log.Fatal(readGroupWriterErr) }
    defer readGroupWriter.Close()
    readGroupBufferedWriter := bufio.NewWriter(readGroupWriter)
    defer readGroupBufferedWriter.Flush()
    for barcode, _ := range readGroups {
        readGroupBufferedWriter.WriteString(barcode)
        readGroupBufferedWriter.WriteString("\n")
    }
}

func GetCellBarcodes() map[string]bool {
    cellMap := map[string]bool{}
    cellFile, cellFileErr := os.Open(*cells)
    if cellFileErr != nil { log.Fatal(cellFileErr) }
    defer cellFile.Close()
    scanner := bufio.NewScanner(cellFile)
    header := true
    for scanner.Scan() {
        if header {
            header = false
            continue
        }
        cellBc := strings.Split(scanner.Text(),",")[0]
        cellMap[cellBc] = true
    }
    return cellMap 
}
