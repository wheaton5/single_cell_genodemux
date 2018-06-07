#!/software/python-2.7.6/bin/python
import subprocess
import pyfasta
import pysam

def split(args):
    chunks = []
    fasta = pyfasta.Fasta(args.fasta)

    for chrom in sorted(fasta.keys()):
        chunk_size = 1000000000
        start = 0
        while start < len(fasta[chrom]):
            end = min(len(fasta[chrom]), start+chunk_size)
            chunks.append({'chrom':chrom, 'start':start, 'end':end ,'__mem_gb':8,'__threads':1})
            start += chunk_size
    return {'chunks':chunks,'join':{'mem_bg':8}}


def main(args, outs):
    outs.coerce_strings()
    args.coerce_strings()
    command = ['cell_and_format_bam','-bam='+args.bam,'-cells='+args.cells,'-chrom='+args.chrom,'-start='+str(args.start),'-end='+str(args.end),'-bamout='+outs.bam+"tmp.bam", '-readGroupsOut='+outs.read_groups,'-statsFile='+outs.bam+"stats.csv", '-readsPerCell='+outs.reads_per_cell]
    print " ".join(command)
    subprocess.check_call(command)
    
    subprocess.check_call(['samtools','sort',outs.bam+"tmp.bam",'-o',outs.bam])
    
    with open(outs.bam+"stats.csv") as stats: 
        stats.readline()
        tokens = stats.readline().strip().split(",")
        outs.no_cell_barcode = int(tokens[0])
        outs.no_cell_call = int(tokens[1])
        outs.reads_total = int(tokens[2])

def join(args, outs, chunk_defs, chunk_outs):
    read_groups = set()
    for chunk_out in chunk_outs:
        with open(chunk_out.read_groups) as rgs:
            for line in rgs:
                read_groups.add(line.strip())
    #outs.bams = [chunk_out.bam for chunk_out in chunk_outs]
    print len(read_groups)
    reads_per_cell = {}
    for chunk_out in chunk_outs:
        with open(chunk_out.reads_per_cell) as reads:
            for line in reads:
                tokens = line.strip().split(",")
                reads_per_cell.setdefault(tokens[0],0)
                reads_per_cell[tokens[0]] += int(tokens[1])
    with open(outs.reads_per_cell,'w') as rpc:
        for key, val in reads_per_cell.iteritems():
            rpc.write(key+","+str(val)+"\n")
    
    outs.no_cell_barcode = 0
    outs.no_cell_call = 0
    outs.reads_total = 0
    for chunk_out in chunk_outs:
        outs.no_cell_barcode += chunk_out.no_cell_barcode
        outs.no_cell_call += chunk_out.no_cell_call
        outs.reads_total += chunk_out.reads_total
    bamin = pysam.AlignmentFile(args.bam)
    outs.unmapped = bamin.unmapped
    outs.mapped = bamin.mapped
    header = bamin.header
    read_group_header = []
    for read_group in read_groups:
        single_read_group = {'ID':read_group,'SM':read_group}
        if 'RG' in bamin.header and len(bamin.header['RG']) > 0:
            for key, value in bamin.header['RG'][0].iteritems():
                if key == 'ID' or key == 'SM':
                    continue
                else:
                    single_read_group[key] = value
        read_group_header.append(single_read_group)
    header['RG'] = read_group_header
    bamout = pysam.AlignmentFile(outs.header,'wh', header=header)
    bamout.close()
    
    bams = [chunk.bam for chunk in chunk_outs]
    command = ['samtools','cat','-h',outs.header]
    command.extend(bams)
    print command
    with open(outs.bam,'w') as outbam:
        subprocess.check_call(command,stdout=outbam)
    #with open(outs.bam,'w') as outbam:
    #    subprocess.check_call(['samtools','sort',outs.bam+"unsorted.bam"],stdout=outbam)
    #subprocess.check_call(['rm',outs.bam+"unsorted.bam"])
    subprocess.check_call(['samtools','index',outs.bam])

def oldmain(args,outs):
    cells = set()
    with open(args.cells) as cellsfile:
        cellsfile.readline() #header
        for line in cellsfile:
            tokens = line.strip().split(',')
            cells.add(tokens[0]) 
    
    bamin = pysam.AlignmentFile(args.bam)
    bamout = pysam.AlignmentFile(outs.bam+"tmp.bam", 'wb', template=bamin)
    read_groups = {}
    for read in bamin.fetch(args.chrom, int(args.start), int(args.end)):
        if read.has_tag("CB"):
            if read.get_tag("CB") in cells:
                read_groups[read.get_tag("CB")] = True
    
                read.set_tag("RG", read.get_tag("CB"))
                r1 = read
                reference_offset_start = r1.pos
                reference_offset_end = r1.pos
                read_offset_start = 0
                read_offset_end = 0
                cigar_tuples = [tup for tup in r1.cigartuples]
                seq = r1.query_sequence
                quals = r1.query_qualities
                cigar_build = []
                for cigar in cigar_tuples:
                    cig = cigar[0]
                    run = cigar[1]
                    if cig == 0 or cig == 4:
                        cigar_build.append(cigar)
                        reference_offset_end += run
                        read_offset_end += run
                    elif cig == 1:
                        cigar_build.append(cigar)
                        read_offset_end += run
                    elif cig == 2:
                        cigar_build.append(cigar)
                        reference_offset_end += run
                    elif cig == 5:
                        cigar_build.append(cigar) 
                    elif cig == 3:
                        r1.reference_start = reference_offset_start
                        r1.query_sequence = seq[read_offset_start:read_offset_end]
                        r1.cigartuples = cigar_build
                        r1.query_qualities = quals[read_offset_start:read_offset_end]
                        bamout.write(r1)
                        reference_offset_start = reference_offset_end
                        read_offset_start = read_offset_end
                        cigar_build = []
                    else:
                        assert(False)
                if len(cigar_build) > 0:
                    r1.reference_start = reference_offset_start
                    r1.query_sequence = seq[read_offset_start:read_offset_end]
                    r1.cigartuples = cigar_build
                    r1.query_qualities = quals[read_offset_start:read_offset_end]
                    bamout.write(r1)
    bamin.close()
    bamout.close()
    subprocess.check_call(['samtools','sort',outs.bam+"tmp.bam",outs.bam[:-4]])
    #gatk = args.gatk
    #command = ['java','-jar','-Xmx8g', '-XX:ParallelGCThreads=1',gatk,'SplitNCigarReads','-R',args.fasta,'--input',outs.bam+"tmp.bam",'--output',outs.bam,'--process-secondary-alignments','true','--refactor-cigar-string','true']
    #subprocess.check_call(command)
    outs.read_groups = read_groups
 
