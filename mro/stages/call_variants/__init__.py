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
    return {'chunks':chunks,'join':{'__mem_gb':2}}


def main(args, outs):
    outs.coerce_strings()
    args.coerce_strings()
    region = args.chrom+":"+str(args.start)+"-"+str(args.end)
    num_regions = 0
    fasta = pyfasta.Fasta(args.fasta)
    with open(args.gtf) as gtf:
        with open(outs.bed+"tmp.bed",'w') as bed:
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
                start = max(0, int(tokens[3]) - args.call_near_exon_window)
                end = min(len(fasta[chrom])-1,int(tokens[4])+args.call_near_exon_window)
                if end > start:
                    if start >= args.start and start <= args.end:
                        bed.write(chrom+"\t"+str(start)+"\t"+str(end)+"\n")
                        num_regions += 1
            if num_regions == 0:
                bed.write(args.chrom+"\t"+str(args.start)+"\t"+str(int(args.start)+2)+"\n")
    with open(outs.bed,'w') as bed:
        subprocess.check_call(['bedtools','merge','-i',outs.bed+"tmp.bed"],stdout=bed)
    with open(outs.vcf+"ploidy.txt",'w') as ploidy_file:
        for key in fasta.keys():
            ploidy_file.write(key+"\t1\t"+str(len(fasta[key]))+"\tF\t"+str(args.ploidy)+"\n")

    with open(outs.vcf[:-7]+"tmp.bcf",'w') as vcf:
        #subprocess.check_call(['freebayes','-0','--pooled-continuous', '--use-best-n-alleles','4','--targets',outs.bed,'--fasta-reference',args.fasta,args.bam],stdout = vcf)
        #subprocess.check_call(['freebayes','-0','--no-mnps','--no-complex','-p',str(args.ploidy), '--use-best-n-alleles','4','--region',region,'--fasta-reference',args.fasta,args.bam],stdout = vcf)
        command = ['samtools','mpileup','-Bg','-m','3','-r',region,'-f',args.fasta,'-q',str(args.min_mapq)]
        if args.input_variants:
            command.extend(['-l',args.input_variants])
        elif args.call_exome_only:
            command.extend(['-l',outs.bed])
        command.append(args.bam)
        #proc = subprocess.Popen(command,stdout=subprocess.PIPE)
        #subprocess.check_call(['bcftools','call','-m','--variants-only','--ploidy',args.ploidy],stdin=proc.stdout,stdout=vcf, shell=True)
        #proc.wait()
        print " ".join(command)
        subprocess.check_call(command, stdout = vcf)
    with open(outs.vcf[:-7]+".bcf",'w') as vcf:
        subprocess.check_call(['bcftools','call','-m','--variants-only','--ploidy-file',outs.vcf+"ploidy.txt",outs.vcf[:-7]+"tmp.bcf"],stdout=vcf)
    with open(outs.vcf[:-3],'w') as vcf:
        subprocess.check_call(['bcftools','view',outs.vcf[:-7]+".bcf"],stdout=vcf)
    subprocess.check_call(['bgzip',outs.vcf[:-3]])
    subprocess.check_call(['tabix','-p','vcf',outs.vcf])

def join(args, outs, chunk_defs, chunk_outs):
    vcfs = []
    for chunk_out in chunk_outs:
        vcfs.append(chunk_out.vcf)
    outs.vcfs = vcfs
    command = ['bcftools','concat','-O','v']
    command.extend(vcfs)
    with open(outs.vcf[:-3]+"tmp.vcf",'w') as tmp:
        subprocess.check_call(command,stdout=tmp)
    with open(outs.vcf[:-3],'w') as tmp:
        subprocess.check_call(['sort', '-k1,1','-k2,2n',outs.vcf[:-3]+"tmp.vcf"],stdout=tmp)
        
        subprocess.check_call(['bgzip',outs.vcf[:-3]])
        subprocess.check_call(['tabix','-p','vcf', outs.vcf])
