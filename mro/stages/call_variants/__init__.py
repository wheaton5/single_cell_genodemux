#!/software/python-2.7.6/bin/python
import subprocess
import pyfasta

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
    return {'chunks':chunks,'join':{'mem_bg':2}}


def main(args, outs):
    outs.coerce_strings()
    args.coerce_strings()
    region = args.chrom+":"+str(args.start)+"-"+str(args.end)
    num_regions = 0
    with open(args.gtf) as gtf:
        with open(outs.bed,'w') as bed:
            for line in gtf:
                if line.startswith("#"):
                    continue
                tokens = line.strip().split()
                chrom = tokens[0]
                if not chrom == args.chrom:
                    continue
                reg_type = tokens[2]
                if not reg_type == "exon":
                    continue
                start = int(tokens[3])
                end = int(tokens[4])
                if end > start:
                    if start >= args.start and start <= args.end:
                        bed.write(chrom+"\t"+str(start)+"\t"+str(end)+"\n")
                        num_regions += 1
            if num_regions == 0:
                bed.write(args.chrom+"\t"+str(args.start)+"\t"+str(int(args.start)+2)+"\n")
            

    with open(outs.vcf[:-7]+".bcf",'w') as vcf:
        #subprocess.check_call(['freebayes','-0','--pooled-continuous', '--use-best-n-alleles','4','--targets',outs.bed,'--fasta-reference',args.fasta,args.bam],stdout = vcf)
        #subprocess.check_call(['freebayes','-0','--dont-left-align-indels','--no-indels','--no-mnps','--no-complex','-p',str(args.ploidy), '--use-best-n-alleles','4','--targets',outs.bed,'--fasta-reference',args.fasta,args.bam],stdout = vcf)
        command = ['samtools','mpileup','-Bg','-m','3','-r',region,'-f',args.fasta,'-q','200',args.bam]

        subprocess.check_call(command,stdout=vcf)
    with open(outs.vcf[:-3],'w') as vcf:
        subprocess.check_call(['bcftools','view',outs.vcf[:-7]+".bcf"],stdout=vcf)
    subprocess.check_call(['bgzip',outs.vcf[:-3]])
    subprocess.check_call(['tabix','-p','vcf',outs.vcf])

def join(args, outs, chunk_defs, chunk_outs):
    assert(False)
    outs.vcfs = {}
    outs.vcf_indexes = {}
    sample_vcfs = {}
    for chunk_out in chunk_outs:
        vcfs = sample_vcfs.setdefault(chunk_out.sample_name,[])
        vcfs.append((chunk_out.region, chunk_out.vcf))
    for sample, sample_out in sample_vcfs.iteritems():
        out = sorted(sample_out)
        sample_calls = [x[1] for x in out]
        command = ['bcftools','concat']
        command.extend(sample_calls)
        with open(outs.vcf[0:-6]+sample+".vcf",'w') as tmp:
            subprocess.check_call(command,stdout=tmp)
        subprocess.check_call(['bgzip',outs.vcf[0:-6]+sample+".vcf"])
        subprocess.check_call(['tabix','-p','vcf', outs.vcf[0:-6]+sample+".vcf.gz"])
        outs.vcfs[sample] = outs.vcf[0:-6]+sample+".vcf.gz"
        outs.vcf_indexes[sample] = outs.vcf[0:-6]+sample+".vcf.gz.tbi"
