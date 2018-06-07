package main

import "flag"
import "fmt"
import "os"
import "bufio"
import "github.com/brentp/bigly/bamat"
import "strconv"
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
var stats = flag.String("statsFile","","write statistics out to this file")
var readsPerCell = flag.String("readsPerCell","","file to write reads per cell info to")

func main() {
    flag.Parse()
    if *readsPerCell == "" {
        fmt.Fprintln(os.Stderr, "please provide -readsPerCell")
        os.Exit(-1)
    }
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

    statsFile, statsFileErr := os.Create(*stats)
    if statsFileErr != nil { log.Fatal(statsFileErr) }
    defer statsFile.Close()
    statsFileWriter := bufio.NewWriter(statsFile)
    defer statsFileWriter.Flush()
    noCellBarcode := 0
    noCell := 0
    readsTotal := 0
    
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
    readsPerCellFile, _ := os.Create(*readsPerCell)
    readsPerCellWriter := bufio.NewWriter(readsPerCellFile)
    defer readsPerCellFile.Close()
    readsPerCellMap := map[string]int{}
    for {
        extant = bamIter.Next()
        if !extant {
            break
        }
        readsTotal += 1
        read := bamIter.Record()
        cellBarcode := read.AuxFields.Get(cellBxTag)
        if cellBarcode == nil {
            noCellBarcode += 1
            continue
        }
        _, extant = cellMap[cellBarcode.Value().(string)]
        if !extant {
            noCell += 1
            continue
        }
        readGroups[cellBarcode.Value().(string)] = true //add to list of read groups
        tagIndex := -1
        for index, tag := range read.AuxFields {
            if tag.Tag() == readGroupTag {
                tagIndex = index
                break
            }
        }
        if read.Flags & sam.Duplicate == sam.Duplicate {
            continue
        }
        _, visited := readsPerCellMap[cellBarcode.Value().(string)]
        if visited {
            readsPerCellMap[cellBarcode.Value().(string)]++
        } else {
            readsPerCellMap[cellBarcode.Value().(string)] = 1
        }
        if tagIndex != -1 {
            read.AuxFields[tagIndex], _ = sam.NewAux(readGroupTag, cellBarcode.Value()) //sam.ParseAux([]byte("RG:Z:"+cellBarcode.Value().(string)))
        } else {
            auxes := []sam.Aux(read.AuxFields)
            newAux, _ := sam.NewAux(readGroupTag, cellBarcode.Value())
            auxes = append(auxes,newAux)
            read.AuxFields = sam.AuxFields(auxes)
        }
        hasNs = false
        for _, cigarOp := range read.Cigar {
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
            for _, cigarOp := range read.Cigar {
                cigarOpTyped := sam.CigarOp(cigarOp)
                cigarType := cigarOpTyped.Type()
                if cigarType == sam.CigarSkipped {
                    read.Seq = sam.NewSeq(seq[readOffsetStart:readOffsetEnd])
                    read.Qual = qual[readOffsetStart:readOffsetEnd]
                    read.Cigar = sam.Cigar(cigarBuild)
                    read.Pos = refOffsetStart
                    if !IsValidRecord2(read) {
                        fmt.Fprintln(os.Stderr, "record not valid")
                        fmt.Fprintln(os.Stderr, read.String())
                        os.Exit(-1)
                    }
                    if readOffsetEnd-readOffsetStart > 29 {
                        err = bamWriter.Write(read)
                        if err != nil { log.Fatal(err) }
                    }
                    refOffsetEnd += cigarOpTyped.Len()
                    readOffsetStart = readOffsetEnd 
                    refOffsetStart = refOffsetEnd
                    cigarBuild = []sam.CigarOp{}
                } else {
                    cigarBuild = append(cigarBuild, cigarOpTyped)
                    consumption := cigarType.Consumes()
                    cigarLength := cigarOpTyped.Len()
                    readOffsetEnd += consumption.Query * cigarLength
                    refOffsetEnd += consumption.Reference * cigarLength
                }
            }
            if len(cigarBuild) > 0 {
                read.Seq = sam.NewSeq(seq[readOffsetStart:readOffsetEnd])
                read.Qual = qual[readOffsetStart:readOffsetEnd]
                read.Cigar = sam.Cigar(cigarBuild)
                read.Pos = refOffsetStart
                if !IsValidRecord2(read) {
                    fmt.Fprintln(os.Stderr, "record not valid")
                    fmt.Fprintln(os.Stderr, read.String())
                    os.Exit(-1)
                }
                if readOffsetEnd - readOffsetStart > 29 {
                    err = bamWriter.Write(read)
                    if err != nil { log.Fatal(err) }
                }
            }
        }
    }

    for bx, count := range readsPerCellMap {
        readsPerCellWriter.WriteString(bx+","+strconv.Itoa(count)+"\n")
    }
    readsPerCellWriter.Flush()

    //output metrics
    statsFileWriter.WriteString("no_cell_barcode,cell_barcode_no_cell_called\n")
    statsFileWriter.WriteString(strconv.Itoa(noCellBarcode)+","+strconv.Itoa(noCell)+","+strconv.Itoa(readsTotal)+"\n")
    statsFileWriter.Flush()
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


func IsValidRecord2(r *sam.Record) bool {
        if (r.Ref == nil || r.Pos == -1) && r.Flags&sam.Unmapped == 0 {
                fmt.Println("r.Ref == nil || r.Pos == -1) && r.Flags&Unmapped == 0")
                return false
        }
        if r.Flags&sam.Paired != 0 && (r.MateRef == nil || r.MatePos == -1) && r.Flags&sam.MateUnmapped == 0 {
                fmt.Println("r.Flags&Paired != 0 && (r.MateRef == nil || r.MatePos == -1) && r.Flags&MateUnmapped == 0")
                return false
        }
        if r.Flags&(sam.Unmapped|sam.ProperPair) == sam.Unmapped|sam.ProperPair {
                fmt.Println("r.Flags&(Unmapped|ProperPair) == Unmapped|ProperPair")
                return false
        }
        if r.Flags&(sam.Paired|sam.MateUnmapped|sam.ProperPair) == sam.Paired|sam.MateUnmapped|sam.ProperPair {
                fmt.Println("r.Flags&(Paired|MateUnmapped|ProperPair) == Paired|MateUnmapped|ProperPair")
                return false
        }
        if len(r.Qual) != 0 && r.Seq.Length != len(r.Qual) {

               fmt.Println("len(r.Qual) != 0 && r.Seq.Length != len(r.Qual)")
                return false
        }
        return true
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
