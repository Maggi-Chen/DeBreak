import datetime
#from pysam import VariantFile
import pysam
import sys

def sortallsv(a):
	return [a.split('\t')[0],int(a.split('\t')[1])]

def sortallsv2(a):
	return [a.split('\t')[1],int(a.split('\t')[1])]


def writevcf_pysam(outpath,bampath,ifpoa,reffile):
	bamfile=pysam.AlignmentFile(bampath,'r')
	bamhead=bamfile.header

	if ifpoa:
		reference=open(reffile,'r').read().split('>')[1:]
		refseq={}
		rawins={}
		for chrom in reference:
			chromseq=''.join(chrom.split('\n')[1:])
			refseq[chrom.split('\n')[0].split(' ')[0]]=chromseq
			rawins[chrom.split('\n')[0].split(' ')[0]]=[]
			try:
				samechrom=open(outpath+'sv_raw_calls/'+bampath.split('/')[-1]+'-'+chrom.split('\n')[0].split(' ')[0]+'.debreak.temp','r').read().split('\n')[:-1]
				samechrom=[mm.split('\t')[4]+'\t'+mm.split('\t')[8]  for mm in samechrom if 'I-cigar' in mm ]
				samechrom+=[mm.split('\t')[4]+'\t'+mm.split('\t')[7] for mm in samechrom if 'I-segment' in mm]
				rawins[chrom.split('\n')[0].split(' ')[0]]=samechrom
			except:
				pass

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
	allsv.sort(key=sortallsv)
	svid=1

	for sv in allsv:
		sv=sv.split('\t')
		try:
			newrec=f.new_record(contig=sv[0],start=int(sv[1]),filter='PASS')
		except:
			continue
		newrec.id='DB'+str(svid)

		if 'Duplication' in sv:
			newrec.alleles=('N','<DUP>')
		if 'Inversion' in sv:
			newrec.alleles=('N','<INV>')
		if 'Translocation' in sv:
			newrec.alleles=('N','<TRA>')

		if not ifpoa:
			if 'Deletion' in sv:
				newrec.alleles=('N','<DEL>')
			if 'Insertion' in sv:
				newrec.alleles=('N','<INS>')

		else:
			if 'Deletion' in sv:
				delseq=refseq[sv[0]][int(sv[1]):int(sv[1])+int(sv[2])]
				newrec.alleles=(delseq,'N')
			if 'Insertion' in sv:
				if 'Precise' in sv:
					newrec.alleles=('N',sv[10])
				else:
					goodraw=sv[6].split(';')
					candi=[mm for mm in rawins[sv[0]] if mm.split('\t')[0] in goodraw]
					samelen=[mm for mm in candi if len(mm.split('\t')[1]) == int(sv[2])]
					if samelen!=[]:
						newrec.alleles=('N',samelen[0].split('\t')[1])
					else:
						largelen=[mm for mm in candi if len(mm.split('\t')[1]) > int(sv[2])]
						try:
							newrec.alleles=('N',largelen[0].split('\t')[1][:int(sv[2])])
						except:
							newrec.alleles=('N','<INS>')
				
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




