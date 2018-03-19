#!/software/python-2.7.6/bin/python
import subprocess
import pyfasta

def split(args):
    chunks = []
    fasta = pyfasta.Fasta("/lustre/scratch118/malaria/team222/hh5/ref/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa")

    for (sample, bam) in args.bams.iteritems():
        for key in fasta.keys():
            chunk_size = 10000000
            start = 0
            while start < len(fasta[key]):
                end = min(len(fasta[key]), start+chunk_size)
                region = key+":"+str(start)+"-"+str(end)
                chunks.append({'bam':bam,'sample_name':sample, 'region':region,'__mem_gb':8,'__threads':1})
                start += chunk_size
    return {'chunks':chunks,'join':{'mem_bg':2}}


def main(args, outs):
    outs.coerce_strings()
    args.coerce_strings()
    with open(outs.vcf.strip(".gz"),'w') as vcf:
        subprocess.check_call(['freebayes','-0','--region',args.region,'--fasta-reference','/lustre/scratch118/malaria/team222/hh5/ref/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa',args.bam],stdout = vcf)
    subprocess.check_call(['bgzip',outs.vcf.strip(".gz")])
    subprocess.check_call(['tabix','-p','vcf',outs.vcf])
    outs.region = args.region
    outs.sample_name = args.sample_name

def join(args, outs, chunk_defs, chunk_outs):
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
