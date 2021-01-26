import datetime

def sortallsv(a):
	return int(a.split('\t')[1])


def writevcf(readpath,writepath,chromosomes,prefix,poa,genotype):
	if prefix:
		prefix=prefix='.'
	else:
		prefix=''
	g=open(writepath+prefix+'debreak.vcf','w')
	filedate=str(datetime.datetime.now()).split(' ')[0]
	g.write('##fileformat=VCFv4.2\n##source=DeBreak\n##fileDate='+filedate+'\n')
	g.write('##ALT=<ID=DEL,Description="Deletion">\n##ALT=<ID=DUP,Description="Duplication">\n##ALT=<ID=INV,Description="Inversion">\n##ALT=<ID=TRA,Description="Translocation">\n##ALT=<ID=INS,Description="Insertion">\n')
	g.write('##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">\n##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Mean mapping quality of supporting reads">\n##INFO=<ID=SUPPREAD,Number=1,Type=Integer,Description="number of supporting reads">\n##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">\n')
	g.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">\n##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+readpath+'filelist\n')


	f=open(writepath+'debreak-allsv-merged-final','r')
	allsv=f.read().split('\n')[:-1]

	foutype=[c for c in allsv if 'Translocation' not in c]
	tra=[c for c in allsv if 'Translocation' in c]


	if not poa:
		poainfo=''
	if not genotype:
		gtinfo=''
	svid=1
	for c in chromosomes:
		samechr=[]
		for d in foutype:
			if d.split('\t')[0]==c:
				samechr+=[d]
		samechr.sort(key=sortallsv)
		for d in samechr:
			if poa:
				poainfo='.;'
				if 'Precise' in d:
					poainfo='PRECISE;'
				if 'Imprecise' in d:
					poainfo='IMPRECISE;'

			if genotype:
				gtinfo='\tGT\t.'
				if 'GT=1/1' in d:
					gtinfo='\tGT\t1/1'
				if 'GT=1/0' in d or 'GT=0/1' in d:
					gtinfo='\tGT\t0/1'
			d=d.split('\t')
			if 'Insertion' in d:
				g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<INS>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+1)+';SVTYPE=INS;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[4])))+gtinfo+'\n')
			if 'Deletion' in d:
				g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<DEL>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+int(d[2]))+';SVTYPE=DEL;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[4])))+gtinfo+'\n')
			if 'Duplication' in d:
				g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<DUP>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+int(d[2]))+';SVTYPE=DUP;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[4])))+gtinfo+'\n')
			if 'Inversion' in d:
				g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<INV>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+int(d[2]))+';SVTYPE=INV;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[4])))+gtinfo+'\n')
			svid+=1

	for c in chromosomes:
		samechr=[]
		for d in tra:
			if d.split('\t')[0] == c:
				samechr+=[d]
		samechr.sort(key=sortallsv)

		for d in samechr:
			if poa:
				poainfo='.;'
				if 'Precise' in d:
					poainfo='PRECISE;'
				if 'Imprecise' in d:
					poainfo='IMPRECISE;'
			if genotype:
				gtinfo='\tGT\t.'
				if 'GT=1/1' in d:
					gtinfo='\tGT\t1/1'
				if 'GT=1/0' in d or 'GT=0/1' in d:
					gtinfo='\tGT\t0/1'
			d=d.split('\t')
			g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<TRA>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[2]+';END='+d[3]+';SVTYPE=TRA;SVLEN=NULL;SUPPREAD='+d[4]+';MAPQ='+d[5]+gtinfo+'\n')
			svid+=1	
	g.close()
	print 'VCF file is written.'
	return True


if __name__ == "__main__":
	readpath='/data/scratch/maggic/simulation/svsimulation/all_type/minimap2-pacbio/merged.sort.bam'
	writepath='/data/scratch/maggic/debreak_result_11.15.19/simulation_pacbio_allfunction/'
	chromosomes=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
	prefix=''
	writevcf(readpath,writepath,chromosomes,prefix,not True,True,True)

