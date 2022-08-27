import os
import multiprocessing
import debreak_detect
import pysam


def extra_readname(writepath):
	readinfo={}
	insertion_poa=open(writepath+'insertion-merged','r').read().split('\n')[:-1]
	for c in insertion_poa:
		cinfo='ins__'+c.split('\t')[0]+'__'+c.split('\t')[1]+'__'+c.split('\t')[2]
		readnames=c.split('\t')[6].split(';')
		for readname in readnames:
			readname=readname.split('_cigar')[0].split('_seg')[0]
			if readname not in readinfo:
				readinfo[readname]=[cinfo]
			else:
				readinfo[readname]+=[cinfo]

	deletion_poa=open(writepath+'deletion-merged','r').read().split('\n')[:-1]
	for c in deletion_poa:
		cinfo='del__'+c.split('\t')[0]+'__'+c.split('\t')[1]+'__'+c.split('\t')[2]
		readnames=c.split('\t')[6].split(';')
		for readname in readnames:
			readname=readname.split('_cigar')[0].split('_seg')[0]
			if readname not in readinfo:
				readinfo[readname]=[cinfo]
			else:
				readinfo[readname]+=[cinfo]

	inversion_poa=open(writepath+'inversion-merged','r').read().split('\n')[:-1]
	for c in inversion_poa:
		cinfo='inv__'+c.split('\t')[0]+'__'+c.split('\t')[1]+'__'+c.split('\t')[2]
		readnames=c.split('\t')[6].split(';')
		for readname in readnames:
			readname=readname.split('_cigar')[0].split('_seg')[0]
			if readname not in readinfo:
				readinfo[readname]=[cinfo]
			else:
				readinfo[readname]+=[cinfo]

	duplication_poa=open(writepath+'duplication-merged','r').read().split('\n')[:-1]
	for c in duplication_poa:
		cinfo='dup__'+c.split('\t')[0]+'__'+c.split('\t')[1]+'__'+c.split('\t')[2]
		readnames=c.split('\t')[6].split(';')
		for readname in readnames:
			readname=readname.split('_cigar')[0].split('_seg')[0]
			if readname not in readinfo:
				readinfo[readname]=[cinfo]
			else:
				readinfo[readname]+=[cinfo]

	translocation_poa=open(writepath+'translocation-merged','r').read().split('\n')[:-1]
	for c in translocation_poa:
		cinfo='tra__'+c.split('\t')[0]+'__'+c.split('\t')[1]+'__'+c.split('\t')[2]+'__'+c.split('\t')[3]
		readnames=c.split('\t')[8].split(';')
		for readname in readnames:
			readname=readname.split('_cigar')[0].split('_seg')[0]
			if readname not in readinfo:
				readinfo[readname]=[cinfo]
			else:
				readinfo[readname]+=[cinfo]

	return readinfo


def extra_readseq_bam(readinfo,bampath,writepath,chrom):
	bamfile=pysam.AlignmentFile(bampath,"rb")
	allalignment=bamfile.fetch(chrom,)
	for align in allalignment:
		if align.flag <=16 and align.qname in readinfo:
			for svsupp in readinfo[align.qname]:
				svfile=open(writepath+"debreak_poa_workspace/"+svsupp+".readseq.fa",'a')
				svfile.write('>'+align.qname+'\n'+align.seq+'\n')
				svfile.close()

	return True

def extra_readseq_sam(readinfo,readpath,samfile,writepath):
	sam=open(readpath+samfile,'r')
	a=sam.readline()
	while a!='':
		if a[0]=='@' or int(a.split('\t')[1])>16:
			a=sam.readline()
			continue
		readname=a.split('\t')[0]
		if readname in readinfo:
			for svsupp in readinfo[readname]:
				svfile=open(writepath+"debreak_poa_workspace/"+svsupp+".fa",'a')
				svfile.write('>'+readname+'\n'+a.split('\t')[9]+'\n')
				svfile.close()
		a=sam.readline()
	sam.close()
	return True

def merge_ctgfa(writepath):
	os.system("ls "+writepath+"debreak_poa_workspace/*keepfile.ctg.fa > "+writepath+"debreak_poa_workspace/ctgfalist")
	ctgfalist=open(writepath+"debreak_poa_workspace/ctgfalist",'r').read().split('\n')[:-1]
	allfasta=open(writepath+"debreak_poa_workspace/allfasta",'w')
	for ctgfile in ctgfalist:
		allctg=open(ctgfile,'r').read().split('>')[1:]
		contigname=ctgfile.split('/')[-1].split('.')[0]
		for iii in range(len(allctg)):
			ctg1=allctg[iii].split('\n')[:-1]
			name1=contigname+'.'+str(iii+1)
			seq1=''
			for cc in ctg1[1:]:
				seq1+=cc
			allfasta.write('>'+name1+'\n'+seq1+'\n')
		os.system("rm "+ctgfile)
	allfasta.close()
	return True


def correct_bp(writepath):
	f=open(writepath+'debreak_poa_workspace/allsvpoa.debreak.temp','r')
	allsv_poa=f.read().split('\n')[:-1]
	f.close()
	del_org=open(writepath+'deletion-merged','r').read().split('\n')[:-1]
	ins_org=open(writepath+'insertion-merged','r').read().split('\n')[:-1]
	inv_org=open(writepath+'inversion-merged','r').read().split('\n')[:-1]
	dup_org=open(writepath+'duplication-merged','r').read().split('\n')[:-1]
	tra_org=open(writepath+'translocation-merged','r').read().split('\n')[:-1]
	allsv=del_org+ins_org+inv_org+dup_org+tra_org
	accuratesv=[]
	alldel=[c for c in allsv if 'Del' in c]
	allins=[c for c in allsv if 'Ins' in c]
	alldup=[c for c in allsv if 'Dup' in c]
	allinv=[c for c in allsv if 'Inv' in c]
	alltra=[c for c in allsv if 'Tra' in c]
	for c in allsv_poa:
		allsvs=[]
		c=c.split('\t')
		if 'D-' in c[3] and 'del__' in c[4]:
			allsvs=alldel
		if 'I-' in c[3] and 'ins__' in c[4]:
			allsvs=allins
		if 'DUP' in c[3] and 'dup__' in c[4]:
			allsvs=alldup
		if 'INV' in c[3] and 'inv__' in c[4]:
			allsvs=allinv
		if 'TRA' in c[4] and 'trai__' in c[5]:
			allsvs=alltra
		if allsvs==[] :
			continue
		if 'TRA-segment' not in c :
			try:
				if c[0]==c[4].split('__')[1] and abs(int(c[1])-int(c[4].split('__')[2]))<=500 and 0.5<=float(c[2])/float(c[4].split('__')[3].split('.')[0])<=2:
					pass
			except:
				print (c);quit()
			if c[0]==c[4].split('__')[1] and abs(int(c[1])-int(c[4].split('__')[2]))<=500 and 0.5<=float(c[2])/float(c[4].split('__')[3].split('.')[0])<=2:
				for d in allsvs:
					if d.split('\t')[0]==c[0] and d.split('\t')[1]==c[4].split('__')[2] and c[4].split('__')[3].split('.')[0]==d.split('\t')[2]:
						if 'I-cigar' in c or 'I-segment' in c:
							accuratesv+=[c[0]+'\t'+c[1]+'\t'+c[2]+'\t'+d.split('\t')[3]+'\t'+d.split('\t')[4]+'\t'+d.split('\t')[5]+'\t'+d.split('\t')[6]+'\t'+d.split('\t')[7]+'\t'+d.split('\t')[8]+'\tPrecise\t'+c[7]]
						else:
							accuratesv+=[c[0]+'\t'+c[1]+'\t'+c[2]+'\t'+d.split('\t')[3]+'\t'+d.split('\t')[4]+'\t'+d.split('\t')[5]+'\t'+d.split('\t')[6]+'\t'+d.split('\t')[7]+'\t'+d.split('\t')[8]+'\tPrecise\t<'+c[3].split('-')[0]+'>']
						allsvs.remove(d)
						break
				
		else:			
			if c[0]== c[5].split('__')[1] and c[2]==c[5].split('__')[3] and abs(int(c[1])-int(c[5].split('_')[2]))<=500 and abs(int(c[3])-int(c[5].split('__')[4].split('.')[0]))<=500:
				for d in allsvs:
					if  d.split('\t')[0]==c[0] and d.split('\t')[1]==c[5].split('__')[2] and d.split('\t')[2]==c[2] and d.split('\t')[3]==c[5].split('__')[4].split('.')[0]:
						accuratesv+=[c[0]+'\t'+c[1]+'\t'+c[2]+'\t'+c[3]+'\t'+d.split('\t')[4]+'\t'+d.split('\t')[5]+'\t'+d.split('\t')[6]+'\t'+d.split('\t')[7]+'\t'+d.split('\t')[8]+'\t'+d.split('\t')[9]+'\tPrecise\t<TRA>']
						allsvs.remove(d)
						break


	allsv=alldel+allins+alldup+allinv
	accuratesv+=[c.split('\t')[0]+'\t'+c.split('\t')[1]+'\t'+c.split('\t')[2]+'\t'+c.split('\t')[3]+'\t'+c.split('\t')[4]+'\t'+c.split('\t')[5]+'\t'+c.split('\t')[6]+'\t'+c.split('\t')[7]+'\t'+c.split('\t')[8]+'\tImprecise\t<'+c.split('\t')[8][:3].upper()+'>' for c in allsv]
	accuratesv+=[c+'\tImprecise\t<TRA>' for c in alltra]

	accuratesv.sort(key=sortallsv)

	f_del=open(writepath+'deletion-merged','w')
	f_ins=open(writepath+'insertion-merged','w')
	f_inv=open(writepath+'inversion-merged','w')
	f_dup=open(writepath+'duplication-merged','w')
	f_tra=open(writepath+'translocation-merged','w')

	for c in accuratesv:
		if 'Del' in c:
			f_del.write(c+'\n')
		if 'Ins' in c:
			f_ins.write(c+'\n')
		if 'Inv' in c:
			f_inv.write(c+'\n')
		if 'Dup' in c:
			f_dup.write(c+'\n')
		if 'Tra' in c:
			f_tra.write(c+'\n')
	f_del.close()
	f_ins.close()
	f_dup.close()
	f_inv.close()
	f_tra.close()
	return True


def sortallsv(a):
	return [a.split('\t')[0],int(a.split('\t')[1])]


def poa_bam(bampath,writepath,chromosomes,threads,ref,min_size,max_size):
	
	readinfo=extra_readname(writepath)
	debreak_poa=multiprocessing.Pool(threads)
	os.system("mkdir "+writepath+"debreak_poa_workspace/")
	for chrom in chromosomes:
		debreak_poa.apply_async(extra_readseq_bam,args=(readinfo,bampath,writepath,chrom,))
	debreak_poa.close()
	debreak_poa.join()
	
	os.system("ls "+writepath+"debreak_poa_workspace/*_*.readseq.fa > "+writepath+"debreak_poa_workspace/fastalist")
	
	fastalist=open(writepath+"debreak_poa_workspace/fastalist",'r').read().split('\n')[:-1]
	debreak_poa=multiprocessing.Pool(max(1,threads))
	for svevent in fastalist:
		debreak_poa.apply_async(call_wtdbg2,args=(svevent,))
		#abortable_call_wtdbg2(svevent)
	debreak_poa.close()
	debreak_poa.join()
	
	merge_ctgfa(writepath)
	os.system("minimap2 -a "+ref+"  "+writepath+"debreak_poa_workspace/allfasta -t "+str(threads)+" > "+writepath+"debreak_poa_workspace/allsvpoa.sam ")
		
	debreak_detect.detect_sam("debreak_poa_workspace/allsvpoa.sam",writepath,writepath,chromosomes,min_size,max_size,False,False)
		
	correct_bp(writepath)
	os.system("rm -r "+writepath+"debreak_poa_workspace/")	
	return True
	



def abortable_call_wtdbg2(svevent):
	timeout=1500
	testp=multiprocessing.Process(target=call_wtdbg2,args=(svevent,))
	testp.start()
	testp.join(timeout)
	if testp.is_alive():
		testp.terminate()
	return True
	'''
	testp = multiprocessing.Pool(1)
	testres = testp.apply_async(call_wtdbg2, args=(svevent,))
	try:
		testout = testres.get(timeout)  # Wait timeout seconds for func to complete.
		return testout
	except multiprocessing.TimeoutError:
		print("Aborting due to timeout")
		raise
	'''

def call_wtdbg2(svevent):
	#os.system("wtdbg2 -i "+svevent+" -fo "+svevent[:-11]+".wtdbg -q ")
	#os.system("wtpoa-cns -i "+svevent[:-11]+".wtdbg.ctg.lay.gz -fo "+svevent[:-11]+".wtdbg.ctg.fa -q ")
	#os.system("mv "+svevent[:-11]+".wtdbg.ctg.fa "+svevent[:-11]+".keepfile.ctg.fa ")
	#print ('6')
	#print ('bsalign poa '+svevent+' -o '+svevent[:-11]+".keepfile.ctg.fa")
	os.system('bsalign poa '+svevent+' -o '+svevent[:-11]+".keepfile.ctg.fa > /dev/null 2>&1 ")
	#quit()
	os.system("rm "+svevent[:-11]+".readseq.fa ")
	return True


def poa_sam(readpath,samlist,writepath,threads,ref,chromosomes,min_size,max_size):
	readinfo=extra_readname(writepath)
	debreak_poa=multiprocessing.Pool(threads)
	os.system("mkdir "+writepath+"debreak_poa_workspace/")
	for samfile in samlist:
		debreak_poa.apply_async(extra_readseq_sam,args=(readinfo,readpath,samfile,writepath,))
	debreak_poa.close()
	debreak_poa.join()
	os.system("ls "+writepath+"debreak_poa_workspace/*_*.readseq.fa > "+writepath+"debreak_poa_workspace/fastalist")
	
	fastalist=open(writepath+"debreak_poa_workspace/fastalist",'r').read().split('\n')[:-1]
	debreak_poa=multiprocessing.Pool(threads//4)	
	for svevent in fastalist:
		debreak_poa.apply_async(call_wtdbg2,args=(svevent,))
	debreak_poa.close()
	debreak_poa.join()

	merge_ctgfa(writepath)
	os.system("minimap2 -a "+ref+"  "+writepath+"debreak_poa_workspace/allfasta --secondary=no -t "+str(threads)+" > "+writepath+"debreak_poa_workspace/allsvpoa.sam ")
	
	debreak_detect.detect_sam("debreak_poa_workspace/allsvpoa.sam",writepath,writepath,chromosomes,min_size,max_size,False,False)
	correct_bp(writepath)
	os.system("rm -r "+writepath+"debreak_poa_workspace/")

	return True



