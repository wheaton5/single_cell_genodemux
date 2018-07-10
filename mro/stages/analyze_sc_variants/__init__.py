#!/software/python-2.7.6/bin/python
import subprocess
import pyfasta
import vcf
import pysam
from collections import Counter
import sys

def split(args):
    chunks = []
    fasta = pyfasta.Fasta(args.fasta)

    for key in fasta.keys():
        chunk_size = 10000000
        start = 0
        while start < len(fasta[key]):
            end = min(len(fasta[key]), start+chunk_size)
            chunks.append({'chrom':key,'start':start,'end':end,'__mem_gb':8,'__threads':1})
            start += chunk_size
    return {'chunks':chunks,'join':{'__mem_gb':2}}


def main(args, outs):
    outs.coerce_strings()
    args.coerce_strings()
    fasta = pyfasta.Fasta(args.fasta)
    gtf_filename = args.gtf
    sc_calls_file = args.sc_calls
    gtf_regions = {}
    with open(gtf_filename) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            else:
                tokens = line.strip().split("\t")
                chrom_regions = gtf_regions.setdefault(tokens[0],[])
                chrom_regions.append((int(tokens[3]),int(tokens[4])))
    def in_regions(regions, chrom, pos, window=0):
        if not chrom in regions:
            return False
        regs = regions[chrom]
        for reg in regs:
            if pos >= reg[0]-window and pos <= reg[1]+window:
                return True
        return False
    bam_files = args.bam_files
    bam_files.append(args.sc_bam)
    samples = args.bam_names
    samples.append("sc_bam")
    bams = [pysam.AlignmentFile(bam) for bam in bam_files]
    vcfreader = vcf.Reader(open(args.ground_truth))
    def pair_iter(i1, i2, key1, key2):
        v1 = None
        v2 = None
        while True:
            if v1 is None:
                try:
                    v1 = i1.next()
                except StopIteration:
                    if v2 is not None:
                        yield (None, v2)
                    for x2 in i2:
                        yield (None, x2)
                    break
            if v2 is None:
                try:
                    v2 = i2.next()
                except StopIteration:
                    if v1 is not None:
                        yield (v1, None)
                    for x1 in i1:
                        yield (x1, None)
                    break
            k1 = key1(v1)
            k2 = key2(v2)
            sys.stdout.flush()
            if k1 == k2:
                yield (v1, v2)
                v1 = None
                v2 = None
            elif k1 < k2:
                yield (v1, None)
                v1 = None
            else:
                yield (None, v2)
                v2 = None

    def support(bams, chrom, pos, ref, alt):
        bam_support = []
        pos = int(pos)-1
        for bam in bams:
            pile = [[pileupread.alignment.query_sequence[pileupread.query_position] for pileupread in pile0.pileups if not pileupread.is_del and not pileupread.is_refskip] for pile0 in bam.pileup(chrom,int(pos),int(pos)+1) if pile0.pos == int(pos)]
            if len(pile) > 0:       
                allele_counts = Counter(pile[0])
                if ref in allele_counts:
                    bam_support.append(allele_counts[ref])
                else:
                    bam_support.append(0)
                if alt in allele_counts:
                    bam_support.append(allele_counts[alt])
                else:
                    bam_support.append(0)
            else:
                bam_support.append(0)
                bam_support.append(0)
        return bam_support
    def getGC(fasta, chrom, pos, window):
        seq = fasta[chrom][max(0,pos-window/2):min(len(fasta[chrom])-1,pos+window/2)].upper()
        gc = 0
        for char in seq:
            if char == 'G' or char == 'C':
                gc += 1
        return str(gc/float(min(1,len(seq))))
    with open(outs.analyze_sc_variants,'w') as cells_per_snp:
        vcf_iter = vcfreader.fetch(args.chrom, args.start, args.end)
        with open(outs.analyze_sc_variants+"tmpcalls.tsv",'w') as tmpcalls:
            subprocess.check_call(['tabix',args.sc_calls,args.chrom+":"+str(args.start)+"-"+str(args.end)],stdout=tmpcalls)
        with open(outs.analyze_sc_variants+"tmpcalls.tsv") as calls_iter:
            variant_pair_iter = pair_iter(calls_iter, vcf_iter, lambda x: (x.split("\t")[0], int(x.split("\t")[1]),x.split("\t")[4]), lambda x: (x.CHROM, x.POS, x.ALT[0]))
            for (call_var, vcf_var) in variant_pair_iter:
                if call_var:
                    tokens = call_var.strip().split("\t")
                    chrom = tokens[0]
                    pos = int(tokens[1])
                    ref = tokens[3]
                    alt = tokens[4].split(",")[0]
                    cells = tokens[7].split(";")
                    sc_qual = tokens[5]
                    ref_cells = str(len(cells[0].split(",")))
                    alt_cells = str(len(cells[1].split(",")))
                    in_single_cell = "1"
                else:
                    ref_cells = "0"
                    sc_qual = "-1"
                    alt_cells = "0"
                    in_single_cell = "0"

                if vcf_var:
                    chrom = vcf_var.CHROM
                    pos = vcf_var.POS
                    ref = vcf_var.REF
                    alt = str(vcf_var.ALT[0])
                    qual = str(vcf_var.QUAL)
                    in_vcf = "1"
                else:
                    qual = "-1"
                    in_vcf = "0"
                if len(ref) > 1 or len(alt) > 1 or not len(ref) == len(alt):
                    continue
                in_exon = "1" if in_regions(gtf_regions, chrom, int(pos)) else "0"
                near_exon = "1" if in_regions(gtf_regions, chrom, int(pos), args.near_exon_window) else "0"
                gc30 = getGC(fasta,chrom,pos,30)
                gc100 = getGC(fasta,chrom, pos, 100)
                gc200 = getGC(fasta,chrom, pos,200)
                refbefore = fasta[chrom][pos-2]
                refafter = fasta[chrom][pos]
                to_write = [chrom, str(pos), ref, alt,refbefore, refafter, ref_cells, alt_cells,qual, sc_qual, in_single_cell, in_vcf, gc30, gc100,gc200,str(in_exon),near_exon]
                to_write.extend([str(x) for x in support(bams, chrom, pos, ref, alt)])
                cells_per_snp.write(",".join(to_write)+"\n")
                
                

def join(args, outs, chunk_defs, chunk_outs):
    samples_header = []
    samples = args.bam_names
    for sample in samples:
        samples_header.append(sample+"_ref")
        samples_header.append(sample+"_alt")

    with open(outs.analyze_sc_variants+"header.csv",'w') as header:
        header.write("chrom,pos,ref,alt,refbefore,refafter,ref_cells,alt_cells,wgs_vcf_qual,sc_qual,in_single_cell_vcf,in_wgs_vcf,gc30,gc100,gc200,in_exon,near_exon,"+",".join(samples_header)+"\n")
    with open(outs.analyze_sc_variants,'w') as analyze:
        command = ['cat']
        command.append(outs.analyze_sc_variants+"header.csv")
        command.extend([chunk_out.analyze_sc_variants for chunk_out in chunk_outs])
        subprocess.check_call(command,stdout=analyze)
