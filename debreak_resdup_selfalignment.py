import os
import time

def check_dup(insinfo,writepath):
	allinfo=open(writepath+'debreak_resdup_map_space/'+insinfo+'.paf','r').read().split('\n')[:-1]
	numsupp=0.0
	suppread=[]
	window=300
	for align in allinfo:
		readname=align.split('\t')[0]
		if readname in suppread or int(align.split('\t')[10])<0.5*int(align.split('\t')[1])  or align.split('\t')[4]!='+':
			continue
		pos=int(align.split('\t')[5].split('DeBreak')[-1])
		start=int(align.split('\t')[7])
		end=int(align.split('\t')[8])
		if start<=pos<=end or abs(end-pos)<=window or abs(start-pos)<=window:
			numsupp+=1
			suppread+=[readname]

	return numsupp/int(insinfo.split('__')[3])


def sortdup(a):
	return [a.split('\t')[0],int(a.split('\t')[1])]

def identify_duplication(vcflist,writepath):
	ins=open(writepath+'insertion-merged','r').read().split('\n')[:-1]
	
	insread={}
	for sv in ins:
		insinfo=sv.split('\t')[0]+'__'+sv.split('\t')[1]+'__'+sv.split('\t')[2]+'__'+sv.split('\t')[3]
		readnames=sv.split('\t')[6].split(';')
		for d in readnames:
			insread[d]=insinfo
	os.system("mkdir "+writepath+"debreak_resdup_insertseq/")
	os.system("mkdir "+writepath+"debreak_resdup_refseq/")
	os.system("mkdir "+writepath+"debreak_resdup_map_space/")
	
	for vcffile in vcflist:
		f=open(writepath+'sv_raw_calls/'+vcffile,'r')
		a=f.readline()
		while a!='':
			if 'I-cigar' not in a or a.split('\t')[4] not in insread:
				a=f.readline(); continue
			readname=a.split('\t')[4]
			insinfo=insread[readname]
			g=open(writepath+'debreak_resdup_insertseq/'+insinfo+'.fa','a')
			insseq=a.split('\t')[8]
			g.write('>'+readname+'\n'+insseq+'\n')
			g.close()
			g=open(writepath+'debreak_resdup_refseq/'+insinfo+'.fa','a')
			readseq=a.split('\t')[9]
			bppos=a.split('\t')[7]
			g.write('>'+readname+'DeBreak'+bppos+'\n'+readseq+'\n')
			a=f.readline()
		f.close()
	rescueddup=[]
	trueins=[]
	ff=open(writepath+'debreak_resdup_map_space/maprateinfo','w')
	for sv in ins:
		insinfo=sv.split('\t')[0]+'__'+sv.split('\t')[1]+'__'+sv.split('\t')[2]+'__'+sv.split('\t')[3]
		os.system('minimap2 -k 9 '+writepath+'debreak_resdup_refseq/'+insinfo+'.fa  '+writepath+'debreak_resdup_insertseq/'+insinfo+'.fa -N 100 > '+writepath+'debreak_resdup_map_space/'+insinfo+'.paf')
		try:
			maprate=check_dup(insinfo,writepath)
		except:
			maprate=0.0
		ff.write(insinfo+'\t'+str(maprate)+'\n')

		if maprate>=0.5:
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

