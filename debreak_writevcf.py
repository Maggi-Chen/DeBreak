import datetime
#from pysam import VariantFile
import pysam
import sys

def sortallsv(a):
	return int(a.split('\t')[1])

def sortallsv2(a):
	return [a.split('\t')[1],int(a.split('\t')[1])]

def writevcf(readpath,writepath,chromosomes,prefix,poa,genotype):
	if prefix:
		prefix=prefix+'.'
	else:
		prefix=''
	g=open(writepath+prefix+'debreak.vcf','w')
	filedate=str(datetime.datetime.now()).split(' ')[0]
	g.write('##fileformat=VCFv4.2\n##source=DeBreak\n##fileDate='+filedate+'\n')
	g.write('##ALT=<ID=DEL,Description="Deletion">\n##ALT=<ID=INS,Description="Insertion">\n##ALT=<ID=DUP,Description="Duplication">\n##ALT=<ID=INV,Description="Inversion">\n##ALT=<ID=TRA,Description="Translocation">\n')
	g.write('##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">\n##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Mean mapping quality of supporting reads">\n##INFO=<ID=SUPPREAD,Number=1,Type=Integer,Description="number of supporting reads">\n##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Variant with precise breakpoint position from POA">\n')
	g.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">\n##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
	g.write('##INFO=<ID=MULTI,Number=0,Type=Flag,Description="If the SV is multi-allelic SV">\n')
	g.write('##INFO=<ID=LARGEINS,Number=0,Type=Flag,Description="Large insertion indentified from local assembly">\n')
	g.write('##INFO=<ID=START2,Number=1,Type=Integer,Description="SV start position on the second haplotype of multi-allilic SV">\n')
	g.write('##INFO=<ID=END2,Number=1,Type=Integer,Description="SV end position on the second haplotype of multi-allilic SV">\n')
	g.write('##INFO=<ID=SVLEN2,Number=1,Type=Integer,Description="SV length on the second haplotype of multi-allilic SV">\n')

	g.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+readpath+'\n')


	f=open(writepath+'debreak-allsv-merged-final','r')
	allsv=f.read().split('\n')[:-1]

	foutype=[c for c in allsv if 'Translocation' not in c]
	tra=[c for c in allsv if 'Translocation' in c]

	uniquefour=[c for c in foutype if 'CompoundSV' not in c]
	compoundfour=[c for c in foutype if 'CompoundSV' in c]

	compoundfour.sort(key=sortallsv2)
	ifsecond=0
	compoundsv=[]
	for i in range(len(compoundfour)-1):
		if ifsecond:
			ifsecond=0
			continue
		if compoundfour[i].split('\t')[0]==compoundfour[i+1].split('\t')[0] and  min(int(compoundfour[i].split('\t')[1])+int(compoundfour[i].split('\t')[2]),int(compoundfour[i+1].split('\t')[1])+int(compoundfour[i+1].split('\t')[2]))-max(int(compoundfour[i].split('\t')[1]),int(compoundfour[i+1].split('\t')[1]))>0:
			ifsecond=1
			svinfo1=compoundfour[i].split('\t')
			svinfo2=compoundfour[i+1].split('\t')
			if int(svinfo1[2])/float(svinfo2[2]) >=1.5 or int(svinfo1[2])/float(svinfo2[2]) <=1/1.5:
				if int(svinfo1[3])>=int(svinfo2[3]):
					try:
						compoundsv+=[svinfo1[0]+'\t'+svinfo1[1]+'\t'+svinfo1[2]+'\t'+svinfo1[3]+'\t'+svinfo1[4]+'\t'+svinfo1[5]+'\t'+svinfo1[6]+'\tCompoundSV:'+svinfo2[1]+':'+svinfo2[2]+':'+svinfo2[3]+':'+svinfo2[5]+'\t'+svinfo1[8]+'\t'+svinfo1[9]+'\t'+svinfo1[10]]
					except:
						compoundsv+=[svinfo1[0]+'\t'+svinfo1[1]+'\t'+svinfo1[2]+'\t'+svinfo1[3]+'\t'+svinfo1[4]+'\t'+svinfo1[5]+'\t'+svinfo1[6]+'\tCompoundSV:'+svinfo2[1]+':'+svinfo2[2]+':'+svinfo2[3]+':'+svinfo2[5]+'\t'+svinfo1[8]+'\t'+svinfo1[9]]
				else:
					try:
						compoundsv+=[svinfo2[0]+'\t'+svinfo2[1]+'\t'+svinfo2[2]+'\t'+svinfo2[3]+'\t'+svinfo2[4]+'\t'+svinfo2[5]+'\t'+svinfo2[6]+'\tCompoundSV:'+svinfo1[1]+':'+svinfo1[2]+':'+svinfo1[3]+':'+svinfo1[5]+'\t'+svinfo2[8]+'\t'+svinfo2[9]+'\t'+svinfo2[10]]
					except:
						compoundsv+=[svinfo2[0]+'\t'+svinfo2[1]+'\t'+svinfo2[2]+'\t'+svinfo2[3]+'\t'+svinfo2[4]+'\t'+svinfo2[5]+'\t'+svinfo2[6]+'\tCompoundSV:'+svinfo1[1]+':'+svinfo1[2]+':'+svinfo1[3]+':'+svinfo1[5]+'\t'+svinfo2[8]+'\t'+svinfo2[9]]
						print (len(svinfo2),len(svinfo1))
			else:
				if int(svinfo1[3])>=int(svinfo2[3]):
					compoundsv+=[svinfo1[0]+'\t'+svinfo1[1]+'\t'+svinfo1[2]+'\t'+str(int(svinfo1[3])+int(svinfo2[3]))+'\t'+str(int(svinfo1[4])+int(svinfo2[4]))+'\t'+svinfo1[5]+'\t'+svinfo1[6]+'\tUnique\t'+svinfo1[8]+'\t'+svinfo1[9]+'\t'+svinfo1[10]]
				else:
					compoundsv+=[svinfo2[0]+'\t'+svinfo2[1]+'\t'+svinfo2[2]+'\t'+str(int(svinfo1[3])+int(svinfo2[3]))+'\t'+str(int(svinfo1[4])+int(svinfo2[4]))+'\t'+svinfo2[5]+'\t'+svinfo2[6]+'\tUnique\t'+svinfo2[8]+'\t'+svinfo2[9]+'\t'+svinfo2[10]]
		else:
			svinfo=compoundfour[i].split('\t')
			try:
				compoundsv+=[svinfo[0]+'\t'+svinfo[1]+'\t'+svinfo[2]+'\t'+svinfo[3]+'\t'+svinfo[4]+'\t'+svinfo[5]+'\t'+svinfo[6]+'\tUnique\t'+svinfo[8]+'\t'+svinfo[9]+'\t'+svinfo[10]]
			except:
				compoundsv+=[svinfo[0]+'\t'+svinfo[1]+'\t'+svinfo[2]+'\t'+svinfo[3]+'\t'+svinfo[4]+'\t'+svinfo[5]+'\t'+svinfo[6]+'\tUnique\t'+svinfo[8]+'\t'+svinfo[9]]
			continue

	foutype=uniquefour+compoundsv
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
				poainfo=''
				if 'Precise' in d:
					poainfo='PRECISE;'
			if genotype:
				gtinfo='\tGT\t./.'
				if 'GT=1/1' in d:
					gtinfo='\tGT\t1/1'
				if 'GT=1/0' in d or 'GT=0/1' in d:
					gtinfo='\tGT\t0/1'
			d=d.split('\t')
			if 'Unique' in d or 'rescue_largeins' in d[6]:
				if 'Insertion' in d:
					if 'rescue_largeins' in d[6]:
						g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<INS>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+1)+';SVTYPE=INS;SVLEN='+d[2]+';LARGEINS'+gtinfo+'\n')
					else:
						g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<INS>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+1)+';SVTYPE=INS;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[5])))+gtinfo+'\n')
				if 'Deletion' in d:
					g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<DEL>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+int(d[2]))+';SVTYPE=DEL;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[5])))+gtinfo+'\n')
				if 'Duplication' in d:
					g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<DUP>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+int(d[2]))+';SVTYPE=DUP;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[5])))+gtinfo+'\n')
				if 'Inversion' in d:
					g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<INV>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+int(d[2]))+';SVTYPE=INV;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[5])))+gtinfo+'\n')
				svid+=1
			else:
				if 'Insertion' in d:
					g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<INS>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+1)+';SVTYPE=INS;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[5])))+';MULTI;START2='+d[7].split(':')[1]+';END2='+str(int(d[7].split(':')[1])+1)+';SVLEN2='+d[7].split(':')[2]+gtinfo+'\n')
				if 'Deletion' in d:
					g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<DEL>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+int(d[2]))+';SVTYPE=DEL;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[5])))+';MULTI;START2='+d[7].split(':')[1]+';END2='+str(int(d[7].split(':')[1])+int(d[7].split(':')[2]))+';SVLEN2='+d[7].split(':')[2]+gtinfo+'\n')
				if 'Duplication' in d:
					g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<DUP>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+int(d[2]))+';SVTYPE=DUP;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[5])))+';MULTI;START2='+d[7].split(':')[1]+';END2='+str(int(d[7].split(':')[1])+int(d[7].split(':')[2]))+';SVLEN2='+d[7].split(':')[2]+gtinfo+'\n')
				if 'Inversion' in d:
					g.write(d[0]+'\t'+d[1]+'\tdb'+str(svid)+'\tN\t<INV>\t.\tPASS\t'+poainfo+'SVMETHOD=DeBreak_1.0;CHR2='+d[0]+';END='+str(int(d[1])+int(d[2]))+';SVTYPE=INV;SVLEN='+d[2]+';SUPPREAD='+d[3]+';MAPQ='+str(round(float(d[5])))+';MULTI;START2='+d[7].split(':')[1]+';END2='+str(int(d[7].split(':')[1])+int(d[7].split(':')[2]))+';SVLEN2='+d[7].split(':')[2]+gtinfo+'\n')

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
	print ('VCF file is written.')

	return True



def writevcf_pysam(outpath,bampath):
	bamfile=pysam.AlignmentFile(bampath,'r')
	bamhead=bamfile.header

	filedate=str(datetime.datetime.now()).split(' ')[0]
	vcfheader=pysam.VariantHeader()
	vcfheader.add_sample(bampath)
	vcfheader.add_line('##source=DeBreak\n')
	vcfheader.add_line('##fileDate='+filedate+'\n')

	for i in range(len(bamhead.references)):
		vcfheader.contigs.add(bamhead.references[i],length=bamhead.lengths[i])

	vcfheader.add_line('##ALT=<ID=DEL,Description="Deletion">\n')
	vcfheader.add_line('##ALT=<ID=INS,Description="Insertion">\n')
	vcfheader.add_line('##ALT=<ID=DUP,Description="Duplication">\n')
	vcfheader.add_line('##ALT=<ID=INV,Description="Inversion">\n')
	vcfheader.add_line('##ALT=<ID=TRA,Description="Translocation">\n')

	vcfheader.info.add('CHR2',1,'String','Chromosome for END')
	vcfheader.info.add('END',1,'Integer','End position of the structural variant')
	vcfheader.info.add('MAPQ',1,'Integer','Mean mapping quality of supporting reads')
	vcfheader.info.add('SUPPREAD',1,'Integer','Number of supporting reads')
	vcfheader.info.add('SVLEN',1,'Integer','Length of the SV')
	vcfheader.info.add('SVMETHOD',1,'String','Type of approach used to detect SV')
	vcfheader.info.add('SVTYPE',1,'String','Type of structural variant')
	vcfheader.info.add('PRECISE',0,'Flag','Variant with precise breakpoint position from POA')

	vcfheader.info.add('MULTI',0,'Flag','If the SV is multi-allelic SV')
	vcfheader.info.add('LARGEINS',0,'Flag','Large insertion indentified from local assembly')
	vcfheader.info.add('START2',1,'Integer','SV start position on the second haplotype of multi-allilic SV')
	vcfheader.info.add('END2',1,'Integer','SV end position on the second haplotype of multi-allilic SV')
	vcfheader.info.add('SVLEN2',1,'Integer','SV length on the second haplotype of multi-allilic SV')

	vcfheader.add_meta('FORMAT', items=[('ID',"GT"), ('Number',1), ('Type','String'),('Description','Genotype')])
	vcfheader.add_line('##CommandLine=debreak '+" ".join(sys.argv[1:])+'\n')

	f=pysam.VariantFile(outpath+'debreak.vcf','w',header=vcfheader)
	allsv=open(outpath+'debreak-allsv-merged-final','r').read().split('\n')[:-1]

	svid=1

	for sv in allsv:
		sv=sv.split('\t')
		try:
			newrec=f.new_record(contig=sv[0],start=int(sv[1]),filter='PASS')
		except:
			continue
		newrec.id='DB'+str(svid)
		newrec.ref='N'
		svid+=1
		if 'Translocation' in sv:
			newrec.info['CHR2']=sv[2]
		else:
			newrec.info['CHR2']=sv[0]
		if 'Insertion' in sv:
			newrec.stop=newrec.start+1
		elif 'Translocation' in sv:
			newrec.stop=int(sv[3])
		else:
			newrec.stop=newrec.start+int(sv[2])
		if  'Translocation' in sv:
			newrec.info['SVLEN']=0
			newrec.info['SUPPREAD']=int(sv[4])
			newrec.info['MAPQ']=int(float(sv[5]))
			
		else:
			newrec.info['SVLEN']=int(sv[2])
			if 'rescue_largeins_' in ''.join(sv):
				newrec.info['SUPPREAD']=0
			else:
				newrec.info['SUPPREAD']=int(sv[3])
			newrec.info['MAPQ']=int(float(sv[5]))
		newrec.info['SVMETHOD']='DeBreak'
		if 'Insertion' in sv and 'rescue_largeins_' in ''.join(sv):
			newrec.info['LARGEINS']=True

		if 'Precise' in sv:
			newrec.info['PRECISE']=True
		newrec.info['SVTYPE']='DEL' if 'Deletion' in sv else ('INS' if 'Insertion' in sv else ('INV' if 'Inversion' in sv else ('DUP' if 'Duplication' in sv else 'TRA')))

		if 'GT=1/0' in sv:
			newrec.samples[bampath]['GT']=(0,1)
		if 'GT=1/1' in sv:
			newrec.samples[bampath]['GT']=(1,1)
		if 'GT=./.' in sv:
			newrec.samples[bampath]['GT']=(0,1)
		if 'CompoundSV' in sv:
			newrec.info['MULTI']=True
			
		f.write(newrec)

	f.close()
	return 0



