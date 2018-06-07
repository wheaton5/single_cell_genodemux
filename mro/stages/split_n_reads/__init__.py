#!/software/python-2.7.6/bin/python
import subprocess
import pyfasta
import pysam

def split(args):
    chunks = []

    for bam in args.bams:
        chunks.append({'bam':bam,'__mem_gb':8,'__threads':4})
    return {'chunks':chunks,'join':{'mem_bg':8}}


def main(args, outs):
    outs.coerce_strings()
    args.coerce_strings()
    gatk = args.gatk
    command = ['java','-jar','-Xmx8g', '-XX:+UseSerialGC',gatk,'SplitNCigarReads','-R',args.fasta,'--input',args.bam,'--output',outs.bam,'--process-secondary-alignments','true','--refactor-cigar-string','true', '--reads-per-cell',outs.reads_per_cell]
    subprocess.check_call(command)
    
    

def join(args, outs, chunk_defs, chunk_outs):
    header = args.header
    bams = [chunk.bam for chunk in chunk_outs]
    command = ['samtools','cat','-h',header]
    command.extend(bams)
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
    print command
    with open(outs.bam,'w') as outbam:
        subprocess.check_call(command,stdout=outbam)
    subprocess.check_call(['samtools','index',outs.bam])
    #gatk = '/lustre/scratch118/malaria/team222/hh5/software/gatk-4.0.2.1/gatk-package-4.0.2.1-local.jar'

