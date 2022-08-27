import os
import time


def check_dup(insinfo,writepath):
	allinfo=open(writepath+'debreak_resdup_map_space/'+insinfo+'.paf','r').read().split('\n')[:-1]
	if len(allinfo)==0:
		return (0,0)
	numsupp=0.0
	suppread=[]
	window=300
	leftmost=[]
	pos=int(insinfo.split('__')[2])+2000

	for align in allinfo:
		readname=align.split('\t')[0]
		if readname in suppread or int(align.split('\t')[10])<0.5*int(align.split('\t')[1])  or align.split('\t')[4]!='+':
			continue
		start=int(align.split('\t')[7])
		end=int(align.split('\t')[8])
		if start<=pos<=end or abs(end-pos)<=window or abs(start-pos)<=window:
			numsupp+=1
			suppread+=[readname]
			leftmost+=[int(align.split('\t')[7])]

	leftmost=leftmost[max(3,len(leftmost)//4)]
	bpshift=max(0,pos-leftmost)
	return (numsupp/int(insinfo.split('__')[3]),bpshift)


def sortdup(a):
	return [a.split('\t')[0],int(a.split('\t')[1])]

def identify_duplication(vcflist,writepath,refpath):
	ins=open(writepath+'insertion-merged','r').read().split('\n')[:-1]

	reffile=open(refpath,'r').read().split('>')[1:]
	refseq={}
	for chrom in reffile:
		chromname=chrom.split('\n')[0].split(' ')[0]
		chromseq=''.join(chrom.split('\n')[1:])
		refseq[chromname]=chromseq

	os.system("mkdir "+writepath+"debreak_resdup_insertseq/")
	os.system("mkdir "+writepath+"debreak_resdup_refseq/")
	os.system("mkdir "+writepath+"debreak_resdup_map_space/")
	
	insread={}
	for sv in ins:
		insinfo=sv.split('\t')[0]+'__'+sv.split('\t')[1]+'__'+sv.split('\t')[2]+'__'+sv.split('\t')[3]
		readnames=sv.split('\t')[6].split(';')
		for d in readnames:
			insread[d]=insinfo

		f=open(writepath+"debreak_resdup_refseq/"+insinfo+'.fa','w')
		f.write('>'+insinfo+'\n')
		f.write(refseq[sv.split('\t')[0]][int(sv.split('\t')[1])-2000-int(sv.split('\t')[2]):int(sv.split('\t')[1])+2000+int(sv.split('\t')[2])]+'\n')
		f.close()
	
	
	for vcffile in vcflist:
		f=open(writepath+'sv_raw_calls/'+vcffile,'r')
		a=f.readline()
		while a!='':
			if 'I-cigar' not in a or a.split('\t')[4] not in insread:
				a=f.readline(); continue
			readname=a.split('\t')[0]+'_'+a.split('\t')[1]+'_'+a.split('\t')[2]+'_'+a.split('\t')[4]
			insinfo=insread[a.split('\t')[4]]
			g=open(writepath+'debreak_resdup_insertseq/'+insinfo+'.fa','a')
			insseq=a.split('\t')[8]
			g.write('>'+readname+'\n'+insseq+'\n')
			g.close()
			a=f.readline()
		f.close()

	rescueddup=[]
	trueins=[]
	ff=open(writepath+'debreak_resdup_map_space/maprateinfo','w')
	for sv in ins:
		insinfo=sv.split('\t')[0]+'__'+sv.split('\t')[1]+'__'+sv.split('\t')[2]+'__'+sv.split('\t')[3]
		try:
			os.system('minimap2 -k 9 '+writepath+'debreak_resdup_refseq/'+insinfo+'.fa  '+writepath+'debreak_resdup_insertseq/'+insinfo+'.fa -N 100 > '+writepath+'debreak_resdup_map_space/'+insinfo+'.paf')
			(maprate,bpcorrection)=check_dup(insinfo,writepath)
		except:
			(maprate,bpcorrection)=(0.0,0)
		ff.write(insinfo+'\t'+str(maprate)+'\n')

		if maprate>=0.5:
			sv=sv.split('\t')
			sv[1]=str(int(sv[1])-bpcorrection)
			sv='\t'.join(sv)
			rescueddup+=[sv]
		else:
			trueins+=[sv]
	ff.close()

	f=open(writepath+'insertion-merged','w')
	for sv in trueins:
		f.write(sv+'\n')
	f.close()
	dup=open(writepath+'duplication-merged','r').read().split('\n')[:-1]
	rescueddup=[cc.replace('Insertion','Duplication') for cc in rescueddup]
	dup+=rescueddup
	dup.sort(key=sortdup)
	f=open(writepath+'duplication-merged','w')
	for sv in dup:
		f.write(sv+'\n')
	f.close()
	
	return 0


def clean_ins_dup(writepath):
	alldup=open(writepath+'duplication-merged','r').read().split('\n')[:-1]
	allins=open(writepath+'insertion-merged','r').read().split('\n')[:-1]
	toremove=[]
	for ins in allins:
		
		for dup in alldup:
			if ins.split('\t')[0]==dup.split('\t')[0] and 0.8<=int(ins.split('\t')[2])/float(dup.split('\t')[2])<=1/0.8 and abs(int(ins.split('\t')[1])-int(dup.split('\t')[1]))<=max(500,int(ins.split('\t')[2]),int(dup.split('\t')[2])):
				toremove+=[ins]
				alldup.remove(dup)
				dup=dup.split('\t')
				dup[3]=str(int(dup[3])+int(ins.split('\t')[3]))
				dup[4]=str(int(dup[4])+int(ins.split('\t')[4]))
				dup='\t'.join(dup)
				alldup+=[dup]
				break

	allins=[c for c in allins if c not in toremove]

	f=open(writepath+'insertion-merged','w')
	for c in allins:
		f.write(c+'\n')
	f.close()
	f=open(writepath+'duplication-merged','w')
	for c in alldup:
		f.write(c+'\n')
	f.close()
	return 0




