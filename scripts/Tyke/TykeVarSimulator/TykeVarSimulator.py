## THIS ALL WORKS WITH XINCHANG HG002 ALIGNMENT TO GRCh37/Hg37 REF !!!!!!!!
import sys
import pysam

if len(sys.argv) < 2:
    print("Usage: python vcfgen.py <path_to_bam> <path_to_ref> <output_path_prefix> <OPTIONAL:seed>")
    print("")
    print("e.g. python vcfgen.py chr22.bam hs37d5.fa chr22")
    print("The above generates a chr22SV.vcf and chr22SNV.vcf file")
    exit(0)

for a in sys.argv:
    if a.startswith("AF="):
        minAF=int(a.lstrip("AF="))
        maxAF=int(a.lstrip("AF="))
    else:
        minAF=0.05 # minimum allele frequency SV   
        maxAF=0.25  # maximum allele frequency SV
minAFsnp=0.05   # minimum allele frequency SNV   
maxAFsnp=0.25   # maximum allele frequency SNV
file=sys.argv[1]   # source bam to simulate SVs from
nosv=50 # no of SVs to simulate
nosnv=200 # no of SNVs to simulate
minsvl=50 # minimum SV length
maxsvl=10000 # maximum SV length
maxsnp=10 # maximum size of indel SNV
sub=0.7 # probability of producing a SNP vs indel in SNV
insdelsnv=0.5   # probability of producing INS vs DEL in SNV
insdel=0.7 # probability of producing INS vs DEL in SV
ref_path=sys.argv[2]
out_prefix=sys.argv[3]
SVvcf=f"{out_prefix}SV.vcf" # name of output vcf file for SV
SNVvcf=f"{out_prefix}SNV.vcf" # name of output vcf file for SNV
seed=sys.argv[4]

CHROMS = [str(i) for i in range(1,23)] + ["X","Y"]
CHROMS += [f'chr{c}' for c in CHROMS]

def get_chrom_lengths(bam_file, ref_chroms=CHROMS):
    with pysam.AlignmentFile(bam_file) as bam:
        chromol = dict((c,l) for c,l in zip(bam.references, bam.lengths) if c in ref_chroms)
    return chromol, tuple(chromol.keys())

def genloc(no,file,mincov=20):
    from numpy import random as nran
    from pysam import depth
    
    chromol, chrom = get_chrom_lengths(file)
    locations=[]
    for i in range(no):
        while True:
            ranchrom=nran.choice(chrom)
            loc=str(nran.randint(0,chromol[ranchrom]))
            try:
                cover=int(depth(file,'-r',ranchrom+":"+loc+"-"+loc).rstrip("\n").split("\t")[-1])
            except ValueError:
                continue
            if cover>=mincov:
                break
        locations.append(tuple([ranchrom,loc,cover]))
    return tuple(locations)
def genlocSV(no,file,mincov=20):
    from numpy import random as nran
    from pysam import depth

    chromol, chrom = get_chrom_lengths(file)
    locations=[]
    for i in range(no):
        while True:
            ranchrom=nran.choice(chrom)
            loc=str(nran.randint(0,chromol[ranchrom]))
            for i in locations:
                if i[0]==ranchrom and abs(int(i[1]))-int(loc)<50000:
                    continue
            try:
                cover=int(depth(file,'-r',ranchrom+":"+loc+"-"+loc).rstrip("\n").split("\t")[-1])
            except ValueError:
                continue
            if cover>=mincov:
                break
        locations.append(tuple([ranchrom,loc,cover]))
    return tuple(locations)
def genseq(minl,maxl):
    nucl=tuple(["A","T","C","G"])
    from numpy.random import choice
    res=''
    for i in range(choice(range(minl,maxl+1))):
        res+=choice(nucl)
    return res
def gensnps(maxsnp=maxsnp,sub=sub,insdelsnv=insdelsnv, snplist=[]):
    from numpy.random import choice
    result=[]
    nucl=["A","T","C","G"]
    for i in range(len(snplist)):
        draw = choice(tuple(['snp','indel']), 1, p= [sub,1-sub])    
        if draw=='snp':
            nucl=["A","T","C","G"]
            nucl.remove(snplist[i][3][0].upper())
            result.append(tuple(list(snplist[i][0:3])+[snplist[i][3][0].upper()]+[choice(nucl,1)[0]]))
        else:
            draw = choice(tuple(['in','del'],), 1, p= [insdelsnv,1-insdelsnv])    
            if draw=='in':
                result.append(tuple(list(snplist[i][0:3])+[snplist[i][3][0].upper()]+[snplist[i][3][0].upper()+genseq(1,maxsnp-1).upper()]))
            else:
                rmlen=choice(range(maxsnp+1)[1:])
                result.append(tuple(list(snplist[i][0:3])+[snplist[i][3].upper()[:rmlen]]+[snplist[i][3].upper()[0]]))

    return tuple(result)
def getrefsnp(reffile,snplist=1):
    with open(reffile,'r') as f:
        data=[]
        cta=[]
        for i in snplist:
            cta.append(i[0])
        cta=set(cta)
        start=1
        for line in f:
            if line.startswith(">"):
                if start==1:
                    start=0
                else:
                    if analyse==1:
                        data.append(tuple([chromosome,seq]))
                if "chr"==snplist[0][0][0:3]:
                    chromosome=line.lstrip(">").split(" ")[0].rstrip("\n")
                else:
                    chromosome=line.lstrip(">").split(" ")[0].lstrip("chr").rstrip("\n")
                seq=''
                analyse=0
                if chromosome in cta:
                    analyse=1
            else:
                if analyse==0:
                    continue
                seq+=line.rstrip("\n") 
        if analyse==1:
            data.append(tuple([chromosome,seq]))
        data=tuple(data)

    f.close()
    snplistmod=[]
    for i in snplist:
        for j in data:
            if i[0]==j[0]:
                snplistmod.append(tuple(list(i)+[j[1][int(i[1])-1:int(i[1])+maxsnp]]))
                break
    return tuple(snplistmod)
    


def main():
    from math import ceil
    from numpy.random import choice
    from numpy.random import uniform
    from numpy.random import seed as npseed
    if seed.isdigit()==1:
        npseed(int(seed))
    svloc=genlocSV(nosv,file,ceil(1/minAF))
    snvloc=genloc(nosnv,file,ceil(1/minAF))
    
    #print(snvloc)
    vcfsv=[
        '##fileformat=VCFv4.2',
        '##ALT=<ID=INS,Description="Insertion">',
        '##ALT=<ID=DEL,Description="Deletion">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])
    ]
    insertno=1
    delno=1
    for i in svloc:
        draw = choice(tuple(['in','del']), 1, p=[insdel,1-insdel])    
        if draw=='in':
            seq=genseq(minsvl,maxsvl)
            vcfsv.append(str(i[0])+"\t"+str(i[1])+"\tHackIns"+str(insertno)+"\tN\t"+seq+"\t60\tPASS\tPRECISE;SVTYPE=INS;SVLEN="+str(len(seq))+";END="+str(int(i[1])+1)+";AF="+str(round(uniform(minAF,maxAF),2))+"\tGT:GQ\t0/0:60")
            #PRECISE;SVTYPE=INS;SVLEN=333;END=748218 AF \t GT:GQ:DR:DV \t	0/0:28:28:5
            insertno+=1
        else:
            dellen=choice(range(minsvl,maxsvl))
            vcfsv.append(str(i[0])+"\t"+str(i[1])+"\tHackDel"+str(delno)+"\tN\t<DEL>\t60\tPASS\tPRECISE;SVTYPE=DEL;SVLEN=-"+str(dellen)+";END="+str(int(i[1])+dellen)+";AF="+str(round(uniform(minAF,maxAF),2))+"\tGT:GQ\t0/0:60")
            delno+=1
    with open(SVvcf,"w") as f:
        for i in tuple(vcfsv)[:-1]:
            f.write(i+'\n')
        f.write(tuple(vcfsv)[-1])
    f.close()
    ##print(snvloc)
    vcfsnv=[
        '##fileformat=VCFv4.2',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])
    ]
    snps=getrefsnp(ref_path,snvloc)
    
    from math import log
    
    for i in gensnps(maxsnp=maxsnp,sub=sub, snplist=snps):
        #TODO
        AFno=(round(uniform(minAF,maxAF),2))
        readno=ceil(AFno*i[2])
        vcfsnv.append(str(i[0])+'\t'+str(i[1])+'\t.\t'+str(i[3])+'\t'+str(i[4])+"\t1500\tPASS\tAF="+str(AFno)+"\tGT:AD\t0/0:"+str(i[2]-readno)+":"+str(readno))
    with open(SNVvcf,"w") as f:
        for i in tuple(vcfsnv)[:-1]:
            f.write(i+'\n')
        f.write(tuple(vcfsnv)[-1])
    f.close()
    #for i in gensnps():
    #    print(i) 
if __name__=="__main__":
    main()
