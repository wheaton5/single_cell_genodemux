#!/software/python-2.7.6/bin/python
import subprocess
import pyfasta
import pysam

def split(args):
    chunks = []
    fasta = pyfasta.Fasta(args.fasta)

    for chrom in sorted(fasta.keys()):
        chunk_size = 20000000
        start = 0
        while start < len(fasta[chrom]):
            end = min(len(fasta[chrom]), start+chunk_size)
            chunks.append({'chrom':chrom, 'start':start, 'end':end ,'__mem_gb':8,'__threads':1})
            start += chunk_size
    return {'chunks':chunks,'join':{'mem_bg':8}}


def main(args, outs):
    outs.coerce_strings()
    args.coerce_strings()
    region = args.chrom+":"+str(args.start)+"-"+str(args.end)
    command = ['freebayes','-0','--region', region,'--fasta-reference',args.fasta, '--ploidy',str(args.ploidy)]
    command.extend(args.bams)
    print " ".join(command)
    with open(outs.vcf[:-3],'w') as vcf:
        subprocess.check_call(command, stdout = vcf)
    subprocess.check_call(['bgzip',outs.vcf[:-3]])
    subprocess.check_call(['tabix', '-p','vcf', outs.vcf])
    
def join(args, outs, chunk_defs, chunk_outs):
    vcfs = []
    for chunk_out in chunk_outs:
        vcfs.append(chunk_out.vcf)
    command = ['bcftools','concat']
    command.extend(vcfs)
    with open(outs.vcf[:-3],'w') as tmp:
        subprocess.check_call(command, stdout=tmp)
    subprocess.check_call(['bgzip',outs.vcf[:-3]])
    subprocess.check_call(['tabix','-p','vcf',outs.vcf])

 
