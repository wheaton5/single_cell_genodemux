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
            chunks.append({'chrom':chrom, 'start':start, 'end':end ,'__mem_gb':24,'__threads':1})
            start += chunk_size
    return {'chunks':chunks,'join':{'__mem_gb':16}}


def main(args, outs):
    outs.coerce_strings()
    args.coerce_strings()
    
    cells = set()
    with open(args.cells) as cellsfile:
        for line in cellsfile:
            barcode = line.strip().split(",")[0]
            cells.add(barcode)
    bam = pysam.AlignmentFile(args.bam)
    recent_umis = {} # map from bc+umi to last known location
    with open(outs.bam+".fq",'w') as fastq:
        for (index, read) in enumerate(bam.fetch(args.chrom, args.start, args.end)):
            if not read.has_tag("CB"):
                continue
            cell_barcode = read.get_tag("CB")
            if not cell_barcode in cells:
                continue
            if not read.has_tag("UB"):
                continue
            UMI = read.get_tag("UB")
            if read.is_secondary or read.is_supplementary:
                continue
            pos = read.pos
            full_umi = cell_barcode + UMI + str(pos)
            if full_umi in recent_umis:
                continue
            recent_umis[full_umi] = pos
            if index % 1000000 == 0:
                keys_to_remove = []
                for key, val in recent_umis.iteritems():
                    if val - pos > 20000:
                        keys_to_remove.append(key)
                for key in keys_to_remove:
                    del recent_umis[key]
            fastq.write("@"+read.qname+";"+cell_barcode+";"+UMI+"\n")
            fastq.write(read.seq+"\n")
            fastq.write("+\n")
            fastq.write(read.qual+"\n")
    print "bwa running"
    with open(outs.bam+"tmp.bam",'w') as tmpbam:
        #ps = subprocess.Popen(['bwa','mem',args.fasta,outs.bam+".fq"], stdout=subprocess.PIPE)
        command = ['minimap2', '-ax', 'splice', '-t','1','-G'+str(args.max_intron_length/1000)+'k', '-k','21','-w','11','--sr','-A2','-B8','-O12,32', '-E2,1', '-r200', '-p.5', '-N20', '-f1000,5000', '-n2', '-m20', '-s40', '-g2000', '-2K50m', '--secondary=no',args.fasta,outs.bam+".fq"]
        #command = ['hisat2','-x',args.fasta[:-3],'-U',outs.bam+".fq"]
        print " ".join(command)
        ps = subprocess.Popen(command,stdout=subprocess.PIPE)
        subprocess.check_call(['samtools','view','-1','-S'],stdin=ps.stdout,stdout=tmpbam)
        ps.wait()
    print "bwa done"
    bamout = pysam.AlignmentFile(outs.bam,'wb', template=bam)
    bwabam = pysam.AlignmentFile(outs.bam+"tmp.bam")
    read_groups = set()
    reads_per_cell = {}
    for read in bwabam:
        name = read.qname
        tokens = name.split(";")
        cb = tokens[1]
        ub = tokens[2]
        read_groups.add(cb)
        reads_per_cell.setdefault(cb,0)
        reads_per_cell[cb] += 1
        read.set_tag("CB",cb)
        read.set_tag("RG",cb)
        read.set_tag("UB",ub)
        bamout.write(read)
        
    bamout.close()
    with open(outs.read_groups,'w') as rgs:
        for cb in read_groups:
            rgs.write(cb+"\n")
    with open(outs.reads_per_cell,'w') as rps:
        for cb, count in reads_per_cell.iteritems():
            rps.write(cb+","+str(count)+"\n")
    #command = ['cell_and_format_bam','-bam='+args.bam,'-cells='+args.cells,'-chrom='+args.chrom,'-start='+str(args.start),'-end='+str(args.end),'-bamout='+outs.bam+"tmp.bam", '-readGroupsOut='+outs.read_groups,'-statsFile='+outs.bam+"stats.csv", '-readsPerCell='+outs.reads_per_cell]
    #print " ".join(command)
    #subprocess.check_call(command)
    
    #subprocess.check_call(['samtools','sort',outs.bam+"tmp.bam",'-o',outs.bam])
    

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
    #for chunk_out in chunk_outs:
    ##    #outs.no_cell_barcode += chunk_out.no_cell_barcode
    #    outs.no_cell_call += chunk_out.no_cell_call
    #    outs.reads_total += chunk_out.reads_total
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
    bamout = pysam.AlignmentFile(outs.bam+"header.bam",'wh', header=header)
    bamout.close()
    
    bams = [chunk.bam for chunk in chunk_outs]
    command = ['samtools','cat','-h',outs.bam+"header.bam"]
    command.extend(bams)
    print command
    with open(outs.bam+"unsorted.bam",'w') as outbam:
        subprocess.check_call(command,stdout=outbam)
    with open(outs.bam,'w') as outbam:
        subprocess.check_call(['samtools','sort',outs.bam+"unsorted.bam"],stdout=outbam)
    subprocess.check_call(['rm',outs.bam+"unsorted.bam"])
    subprocess.check_call(['samtools','index',outs.bam])

 
