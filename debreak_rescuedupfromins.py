import os
import time

def sortsv(a):
	return [a.split('\t')[0],int(a.split('\t')[1])]


def identify_duplication(vcflist,writepath,min_support,refpath):
	f=open(writepath+'insertion-merged','r')
	ins=f.read().split('\n')[:-1]
	f.close()
	insread={}
	for c in ins:
		insinfo=c.split('\t')[0]+'_'+c.split('\t')[1]+'_'+c.split('\t')[2]
		readnames=c.split('\t')[6].split(';')
		for d in readnames:
			insread[d]=insinfo

	os.system("mkdir "+writepath+"debreak_resdup_insertseq/")
	os.system("mkdir "+writepath+"debreak_resdup_refseq/")
	os.system("mkdir "+writepath+"debreak_resdup_map_space/")

	#extrace insert sequence
	for vcffile in vcflist:
		f=open(writepath+vcffile,'r')
		a=f.readline()
		while a!='':
			if 'I-cigar' not in a:
				a=f.readline()
				continue
			readname=a.split('\t')[4]
			if readname in insread :
				insinfo=insread[readname]
				try:
					if a.split('\t')[0]==insinfo.split('_')[0] and abs(int(a.split('\t')[1])-int(insinfo.split('_')[1]))<=500 and 0.7<=int(a.split('\t')[2])/int(insinfo.split('_')[2])<=1.43:
						g=open(writepath+'debreak_resdup_insertseq/'+insinfo+'.fa','a')
						insseq=a.split('\t')[7]
						g.write('>'+insinfo+'_'+readname+'\n'+insseq+'\n')
						g.close()
				except:
					pass
			a=f.readline()

	#extra ref sequence
	allins=[]
	for c in ins:
	        insinfo=[c.split('\t')[0],int(c.split('\t')[1]),int(c.split('\t')[2])]
		allins+=[insinfo]

	f=open(refpath,'r')
	allref=f.read().split('>')[1:24]
	iii=0
	for i in range(23):
		chrom='chr'+str(i+1)
		if chrom=='chr23':
			chrom='chrX'
		ref=allref[i].split('\n')
		if ref[0].split(' ')[0]!=chrom:
			print 'chrom not correct '+chrom+'\t'+ref[0]
			return 0
		seq=''
		for c in ref[1:]:
			seq+=c
	
		for c in [m for m in allins if m[0]==chrom]:
			refseq=seq[c[1]-2000:c[1]+2000]
			f=open(writepath+'debreak_resdup_refseq/'+c[0]+'_'+str(c[1])+'_'+str(c[2])+'.ref.fa','w')
			f.write('>'+c[0]+'_'+str(c[1])+'_'+str(c[2])+'_local_ref\n'+refseq+'\n')
			f.close()
			iii+=1
	
	#write map.sh

	
	f=open(writepath+'debreak_mapdup.sh','w')
	for c in ins:
        	insinfo=c.split('\t')[0]+'_'+c.split('\t')[1]+'_'+c.split('\t')[2]
		f.write('minimap2 -a -k 9  '+writepath+'debreak_resdup_refseq/'+insinfo+'.ref.fa  '+writepath+'debreak_resdup_insertseq/'+insinfo+'.fa  > '+writepath+'debreak_resdup_map_space/'+insinfo+'.sam\n')
	f.close()


	os.system("bash "+writepath+"debreak_mapdup.sh")
	os.system("ls "+writepath+"debreak_resdup_map_space/chr*|rev |cut -d '/' -f 1 |rev > "+writepath+"debreak_resdup_map_space/filelist")
	# identify duplication
	f=open(writepath+'debreak_resdup_map_space/filelist','r')
	files=f.read().split('\n')[:-1]
	f.close()
	rescueddup=[]
	h=open(writepath+'duplication_local_map','w')
	for c in files:
		f=open(writepath+'debreak_resdup_map_space/'+c,'r')
		allf=f.read().split('\n')[:-1]
		f.close()
		allf=[m for m in allf if m[0]!='@' and m.split('\t')[1] in ['0','4','16']]
		if len(allf)==0:
			continue
		dup=len([m for m in allf if m.split('\t')[1]!='4'])
		if dup>=min_support/2 and float(dup)/len(allf) >0.5:
			insinfo=c.split('_')
			rescueddup+=[insinfo[0]+'\t'+insinfo[1]+'\t'+insinfo[2].split('.')[0]+'\t'+str(dup)+'\t'+str(float(dup)/len(allf))]
			h.write(insinfo[0]+'\t'+insinfo[1]+'\t'+insinfo[2].split('.')[0]+'\t'+str(dup)+'\t'+str(float(dup)/len(allf))+'\n')
	h.close()
	

	f=open(writepath+'duplication-merged','r')
	debreakdup=f.read().split('\n')[:-1]
	f.close()

	alldup=[]
	for c in rescueddup:
		chr1=c.split('\t')[0]
		pos1=int(c.split('\t')[1])
		len1=int(c.split('\t')[2])
		num=int(c.split('\t')[3])
		testif=0
		for d in debreakdup:
			if d.split('\t')[0]==chr1 and abs(int(d.split('\t')[1])-pos1)<=500 and 0.7<=len1/float(d.split('\t')[2]) <=1.43:
				testif=1; break
		if testif==0:
			alldup+=[c+'\t.\t.\tDuplication']
	vcf=alldup+debreakdup
	vcf.sort(key=sortsv)
	f=open(writepath+'duplication-merged','w')
	for c in vcf:
		f.write(c+'\n')
	f.close()

	f=open(writepath+'insertion-merged','r')
	debreakins=f.read().split('\n')[:-1]
	f.close()
	allins=[]
	for c in debreakins:
		chr1=c.split('\t')[0]
		pos1=int(c.split('\t')[1])
		len1=int(c.split('\t')[2])
		testif=0
		for d in rescueddup:
			if d.split('\t')[0]==chr1 and abs(int(d.split('\t')[1])-pos1)<=500 and 0.7<=len1/float(d.split('\t')[2]) <=1.43:
				testif=1;  break
		if testif==0:	
			allins+=[c]

	f=open(writepath+'insertion-merged','w')
	for c in allins:
		f.write(c+'\n')
	f.close()
	return True

	
