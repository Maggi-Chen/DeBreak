import os
import pysam
import debreak_detect
import multiprocessing
import time

def sort_clip_pos(a):
	return int(a.split('\t')[1])

def find_candi_ins_bp(writefilepath,chrom,minsupp):
	allreadinfo=open(writefilepath+'debreak_ins_workspace/readinfo_start_end_'+chrom,'r').read().split('\n')[:-1]
	startc=[]
	stopc=[]
	for c in allreadinfo:
		c=c.split('\t')
		if c[2]!='0':
			startc+=[int(c[2])]
		if c[4]!='0':
			stopc+=[int(c[4])]

	inscandi=[]

	startc.sort()
	stopc.sort()
	startcandi=[]
	stopcandi=[]


	while startc!=[]:
		start=startc[0]
		allcandi=[]
		for c in startc:
			if c-start<=2000:
				allcandi+=[c]
			else:
				break

		if len(allcandi)<minsupp:
			startc=startc[len(allcandi):]
		else:
			startc=startc[len(allcandi):]
			medean=allcandi[len(allcandi)//2]
			allcandi=[c for c in allcandi if medean-1000<=c<=medean+1000]
			if 1*minsupp<=len(allcandi)<=minsupp*20:
				medean=sum(allcandi)//len(allcandi)
				startcandi+=[[medean,len(allcandi)]]
	
	print ('clean start of '+chrom)
	while stopc!=[]:
		start=stopc[0]
		allcandi=[]
		for c in stopc:
			if c-start<=2000:
				allcandi+=[c]
			else:
				break
		if len(allcandi)<minsupp:
			stopc=stopc[len(allcandi):]
		else:
			stopc=stopc[len(allcandi):]
			medean=allcandi[len(allcandi)//2]
			allcandi=[c for c in allcandi if medean-1000<=c<=medean+1000]
			if minsupp*1<=len(allcandi)<=minsupp*20:
				medean=sum(allcandi)//len(allcandi)
				stopcandi+=[[medean,len(allcandi)]]
	print ('clean stop of '+chrom)

	startc=startcandi;   stopc=stopcandi
	
	f=open(writefilepath+'debreak_ins_workspace/startcan_'+chrom,'w')
	for c in startc:
		f.write(chrom+'\t'+str(c[0])+'\t'+str(c[1])+'\n')
	f=open(writefilepath+'debreak_ins_workspace/stopcan_'+chrom,'w')
	for c in stopc:
		f.write(chrom+'\t'+str(c[0])+'\t'+str(c[1])+'\n')
	f.close()
	
	inscandi=[]

	while startc!=[] and stopc!=[]:
		if abs(startc[0][0]-stopc[0][0])<=200:
			if startc[0][1]+stopc[0][1]>=2*minsupp:
				inscandi+=[chrom+'\t'+str((int(startc[0][0])+int(stopc[0][0]))//2)]
			startc=startc[1:]
			stopc=stopc[1:]
		else:
			if startc[0][0]-stopc[0][0]>0:
				stopc=stopc[1:]
			else:
				startc=startc[1:]
	f=open(writefilepath+'debreak_ins_workspace/debreak_rescueins_chrom_'+chrom,'w')
	for c in inscandi:
		f.write(c+'\n')
	f.close()
	print ('finish '+chrom)

	return True



def extract_good_reads(writepath,bampath,inscandi):
	bamfile=pysam.AlignmentFile(bampath,'rb')
	
	for candi in inscandi:
		chrom=candi.split('\t')[0]
		pos=int(candi.split('\t')[1])
		allalign=bamfile.fetch(chrom,max(0,pos-500),pos+500)
		f=open(writepath+'debreak_ins_workspace/'+chrom+'_'+str(pos)+'.fa','w')
		for align in allalign:
			testif=0
			if align.is_secondary:
				continue
			cigar=align.cigar
			if 4<=cigar[0][0]<=5 and cigar[0][1]>=200 and abs(align.reference_start-pos)<=300:
				testif=1
			if 4<=cigar[-1][0]<=5 and cigar[-1][1]>=200 and abs(align.reference_end-pos)<=300:
				testif=1
			if testif!=1:
				for c in cigar[1:-1]:
					if c[0]==1 and int(c[1])>=1000:
						testif=1;  break
			if testif==1:
				f.write('>'+align.qname+'\n'+align.query_sequence+'\n')
		f.close()
	
	return True




def rescue_ins_bam(bampath,chromosomes,writepath,threads,refpath,min_supp,min_size,max_size):
	
	#os.system("cat "+writepath+"debreak_ins_workspace/readinfo*chr* > "+writepath+"debreak_ins_workspace/readinfoall")
	print (threads)
	debreak_resins=multiprocessing.Pool(threads)
	for i in range(len(chromosomes)):
		debreak_resins.apply_async(find_candi_ins_bp,args=(writepath,chromosomes[i],min_supp,))

	debreak_resins.close()
	debreak_resins.join()
	os.system("cat "+writepath+"debreak_ins_workspace/debreak_rescueins_chrom_*  > "+writepath+"debreak_ins_workspace/debreak_rescueins_allchr")

	f=open(writepath+"debreak_ins_workspace/debreak_rescueins_allchr",'r')
	inscandi=f.read().split('\n')[:-1]
	f.close()
	extract_good_reads(writepath,bampath,inscandi)
	os.system("ls "+writepath+"debreak_ins_workspace/chr*fa  >  "+writepath+"debreak_ins_workspace/falist")
	
	f=open(writepath+"debreak_ins_workspace/falist",'r')
	falist=f.read().split('\n')[:-1]
	f.close()
	os.system('mkdir '+writepath+"debreak_ins_workspace/workspace")
	os.system('mkdir '+writepath+"debreak_ins_workspace/doneassembly/")
	for c in falist:
		#print c
		os.system("wtdbg2 -q -i "+c+" -fo "+writepath+'debreak_ins_workspace/workspace/'+c.split('/')[-1][:-3]+'.wtdbg')
		os.system("wtpoa-cns  -q  -i "+writepath+'debreak_ins_workspace/workspace/'+c.split('/')[-1][:-3]+'.wtdbg.ctg.lay.gz  -fo '+writepath+'debreak_ins_workspace/workspace/'+c.split('/')[-1][:-3]+'.wtdbg.ctg.fa')
		os.system("mv "+writepath+'debreak_ins_workspace/workspace/'+c.split('/')[-1][:-3]+'.wtdbg.ctg.fa  '+writepath+'debreak_ins_workspace/doneassembly/'+c.split('/')[-1][:-3]+'.wtdbg.ctg.fa')
		os.system("rm "+writepath+'debreak_ins_workspace/workspace/'+c.split('/')[-1][:-3]+'.wtdbg*')
	
	os.system('ls '+writepath+"debreak_ins_workspace/doneassembly/chr*wtdbg.ctg.fa > "+writepath+'debreak_ins_workspace/doneassembly/ctgfalist')
	ctgfalist=open(writepath+'debreak_ins_workspace/doneassembly/ctgfalist','r').read().split('\n')[:-1]
	f=open(writepath+'debreak_ins_workspace/debreak_rescuelargeins_merged.fa','w')
	for ctgf in ctgfalist:
		ctg1=open(ctgf,'r').read()
		if ctg1!='':
			ctg1=ctg1.split('>')[1].split('\n')[1:-1]
			seq=''
			for mm in ctg1:
				seq+=mm
			f.write('>'+ctgf.split('/')[-1][:-3]+'.ctg1\n'+seq+'\n')
	f.close()

	
	os.system('minimap2 -a -t '+str(threads)+' '+refpath+' --secondary=no '+writepath+'debreak_ins_workspace/debreak_rescuelargeins_merged.fa  > '+writepath+'debreak_ins_workspace/debreak_rescuelargeins_merged.sam ')
	
	debreak_detect.detect_sam('debreak_rescuelargeins_merged.sam',writepath+'debreak_ins_workspace/',writepath+'debreak_ins_workspace/',chromosomes,min_size,max_size,False,False)
	
	allsv=open(writepath+'debreak_ins_workspace/debreak_rescuelargeins_merged.debreak.temp','r').read().split('\n')[:-1]
	allins=[c for c in allsv if 'I-' in c and int(c.split('\t')[2])>=1000]
	oldsvcalls=open(writepath+'debreak-allsv-merged-final','r').read().split('\n')[:-1]
	#oldsvcalls=open('/data/scratch/maggic/skbr3/newseq/debreak-hg19/debreak-allsv-merged-final-new','r').read().split('\n')[:-1]
	oldinscalls=[c for c in oldsvcalls if 'Insertion' in c or 'Duplication' in c]
	newins=[]
	for c in allins:
		testif=0
		for d in oldinscalls:
			if c.split('\t')[0]==d.split('\t')[0] and abs(int(c.split('\t')[1])-int(d.split('\t')[1]))<=500 and 0.7<=float(c.split('\t')[2])/float(d.split('\t')[2])<=1.43:
				testif=1; break
		if testif==0:
			c=c.split('\t')
			newins+=[c[0]+'\t'+c[1]+'\t'+c[2]+'\tcontig\t0\t'+c[6]+'\trescue_largeins_'+c[4]+'\tInsertion\tPrecise\tGT=./.']
	newins.sort(key=sort_ins)
	testif=0
	while testif==0:
		testif=1
		for i in range(len(newins)-1):
			if newins[i].split('\t')[0]==newins[i+1].split('\t')[0] and abs(int(newins[i].split('\t')[1])-int(newins[i+1].split('\t')[1]))<=500:
				testif=0
				if int(newins[i].split('\t')[2])>int(newins[i+1].split('\t')[2]):
					newins.remove(newins[i+1])
				else:
					newins.remove(newins[i])
				break

	cleanins=[]
	bamf=pysam.AlignmentFile(bampath,'rb')
	for c in newins:
		cov=bamf.count(c.split('\t')[0],int(c.split('\t')[1])-50,int(c.split('\t')[1])+50)
		if cov<=25*min_supp:
			cleanins+=[c]

	f=open(writepath+'debreak_ins_workspace/debreak_rescuelargeins_ins','w')
	for c in cleanins:
		f.write(c+'\n')
	f.close()
	f=open(writepath+'debreak-allsv-merged-final','a')
	for c in cleanins:
		f.write(c+'\n')
	f.close()
	return True


def sort_ins(a):
	return [a.split('\t')[0],int(a.split('\t')[1])]


