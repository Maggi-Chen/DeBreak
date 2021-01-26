import time
import os

def mergerpossort(a):
	return int(a.split('\t')[1])

def mergerlensort(a):
	return int(a.split('\t')[2])

def m_samechr_insertion(samechrom):
	ins=[]
	samechrom.sort(key=mergerpossort)
	samechrom+=['last_end\t999999999999\t999999999999\t1\t60\t0']
	candi=[]
	last=samechrom[0]
	for event in samechrom:
		maxlen=max(int(event.split('\t')[2]),int(last.split('\t')[2]))
		window=max(300,maxlen+10)
		window=min(800,window)
		if int(event.split('\t')[1]) < int(last.split('\t')[1])+window :
			candi+=[event]
			last=sorted(candi,key=mergerlensort)[-1]
			continue
		if len(candi)==1:
			ins+=[candi[0]]
		else:
			position=0; length=0; count=0; quality=0; sd=0; readnames=''
			candi.sort(key=mergerlensort,reverse=True)
			for can in candi:
				if int(can.split('\t')[2])>=0.5*length:
					can_count=int(can.split('\t')[3])
					position=int((position*count+int(can.split('\t')[1])*can_count)/(count+can_count))
					length=int((length*count+int(can.split('\t')[2])*can_count)/(count+can_count))
					quality=(quality*count+float(can.split('\t')[4])*can_count)/(count+can_count)
					sd=(sd*count+float(can.split('\t')[5])*can_count)/(count+can_count)
					readnames+=can.split('\t')[6]+';'
					count+=can_count
			readnames=readnames[:-1]
			ins+=[candi[0].split('\t')[0]+'\t'+str(position)+'\t'+str(length)+'\t'+str(count)+'\t'+str(quality)+'\t'+str(sd)+'\t'+readnames]
		candi=[event]
		last=event
	return ins

def sort_mostspupport(a):
	return [int(a.split('\t')[3]),int(a.split('\t')[2])]

def m_samechr_deletion(samechrom):
	dels=[]
	samechrom.sort(key=mergerpossort)
	samechrom+=['last_end\t999999999999\t999999999999\t1\t60\t0']
	candi=[]
	last=samechrom[0]
	for event in samechrom:
		maxlen=max(int(event.split('\t')[2]),int(last.split('\t')[2]))
		window=max(300,maxlen+10)
		window=min(800,window)
		if int(event.split('\t')[1]) < int(last.split('\t')[1])+int(last.split('\t')[2])+200 and int(event.split('\t')[1])-int(last.split('\t')[1]) < window:
			candi+=[event]
			last=event
			continue
		if len(candi)==1:
			dels+=[candi[0]]
		else:
			position=0; length=0; count=0; quality=0; sd=0; readnames=''
			candi.sort(key=sort_mostspupport,reverse=True)
			for can in candi:
				if int(can.split('\t')[2])>=0.5*length:
					can_count=int(can.split('\t')[3])
					position=int((position*count+int(can.split('\t')[1])*can_count)/(count+can_count))
					length=int((length*count+int(can.split('\t')[2])*can_count)/(count+can_count))
					quality=(quality*count+float(can.split('\t')[4])*can_count)/(count+can_count)
					sd=(sd*count+float(can.split('\t')[5])*can_count)/(count+can_count)
					readnames+=can.split('\t')[6]+';'
					count+=can_count
			readnames=readnames[:-1]
			dels+=[candi[0].split('\t')[0]+'\t'+str(position)+'\t'+str(length)+'\t'+str(count)+'\t'+str(quality)+'\t'+str(sd)+'\t'+readnames]
		candi=[event]
		last=event
	return dels

def mergertra(a):
	return 	[a.split('\t')[2],int(a.split('\t')[3])]

def m_samechr_translocation(samechrom):
	samechrom.sort(key=mergertra)
	iftrue=0
	while iftrue==0:
		iftrue=1
		for i in range(len(samechrom)-1):
			if samechrom[i].split('\t')[2]==samechrom[i+1].split('\t')[2] and abs(int(samechrom[i].split('\t')[3])-int(samechrom[i+1].split('\t')[3]))<=800:
				iftrue=0
				if samechrom[i].split('\t')[0]==samechrom[i+1].split('\t')[0] and abs(int(samechrom[i].split('\t')[1])-int(samechrom[i+1].split('\t')[1]))<=1000:
					count1=int(samechrom[i].split('\t')[4]); count2=int(samechrom[i+1].split('\t')[4])
					pos1=(int(samechrom[i].split('\t')[1])*count1+int(samechrom[i+1].split('\t')[1])*count2)/(count1+count2)
					pos2=(int(samechrom[i].split('\t')[3])*count1+int(samechrom[i+1].split('\t')[3])*count2)/(count1+count2)
					quality=(float(samechrom[i].split('\t')[5])*count1+float(samechrom[i+1].split('\t')[5])*count2)/(count1+count2)
					sd1=(float(samechrom[i].split('\t')[6])*count1+float(samechrom[i+1].split('\t')[6])*count2)/(count1+count2)
					sd2=(float(samechrom[i].split('\t')[7])*count1+float(samechrom[i+1].split('\t')[7])*count2)/(count1+count2)
					readname=samechrom[i].split('\t')[8]+';'+samechrom[i+1].split('\t')[8]
					mergedtra=samechrom[i].split('\t')[0]+'\t'+str(pos1)+'\t'+samechrom[i].split('\t')[2]+'\t'+str(pos2)+'\t'+str(count1+count2)+'\t'+str(quality)+'\t'+str(sd1)+'\t'+str(sd2)+'\t'+readname+'\tTranslocation'
					samechrom=samechrom[:i]+[mergedtra]+samechrom[i+2:]
				else:
					count1=int(samechrom[i].split('\t')[4]); count2=int(samechrom[i+1].split('\t')[4])
					if count1>=count2:
						samechrom.remove(samechrom[i+1])
					else:
						samechrom.remove(samechrom[i])
				break
	return samechrom

def standerd_varition(length):
	avelen=sum(length)/float(len(length))
	s=0.0
	for c in length:
		s+=(c-avelen)**2/float(len(length))
	s=s**0.5
	return s


def mergeinfosecpos(a):
	return int(a.split('\t')[3])

def mergeinfo_translocation(candi,min_support):
	chrom=candi[0].split('\t')[0]
	secchr=[c.split('\t')[2] for c in candi]
	secchr=max(set(secchr),key=secchr.count)
	candi=[c for c in candi if c.split('\t')[2]==secchr]

	candi.sort(key=mergeinfosecpos)
	if len(candi)%2==0:
		median=(int(candi[len(candi)/2-1].split('\t')[3])+int(candi[len(candi)/2].split('\t')[3]))/2
	else:
		median=int(candi[len(candi)/2].split('\t')[3])
	candi=[ c for c in candi if abs(int(c.split('\t')[3])-median)<=800]
	
	if len(candi)>=min_support:
		pos1=[int(c.split('\t')[1]) for c in candi]
		pos2=[int(c.split('\t')[3]) for c in candi]
		qual=[int(c.split('\t')[4]) for c in candi]
		sd1=standerd_varition(pos1)
		sd2=standerd_varition(pos2)
		readnames=''
		for c in candi:
			readnames+=c.split('\t')[5]+';'
		readnames=readnames[:-1]
		return [chrom+'\t'+str(sum(pos1)/len(pos1))+'\t'+secchr+'\t'+str(sum(pos2)/len(pos2))+'\t'+str(len(candi))+'\t'+str(sum(qual)/len(qual))+'\t'+str(sd1)+'\t'+str(sd2)+'\t'+readnames+'\tTranslocation']
	else:
		return []



def mergeinfolengthsort(a):
	return int(a.split('\t')[2])

def mergeinfo_insertion(candi):
	ins=[]
	candi.sort(key=mergeinfolengthsort)
	while len(candi)>2:
		if int(candi[-1].split('\t')[2]) > 1.5* int(candi[len(candi)/2].split('\t')[2]):
			candi.remove(candi[-1])
			continue
		if int(candi[len(candi)/2].split('\t')[2]) >  1.5*int(candi[0].split('\t')[2]):
			candi.remove(candi[0])
			continue
		break
	if len(candi)>2:
		chrom=candi[0].split('\t')[0]
		position=[int(c.split('\t')[1]) for c in candi]
		length=[int(c.split('\t')[2]) for c in candi]
		quality=[float(c.split('\t')[3]) for c in candi]
		position=sum(position)/len(position)
		quality=sum(quality)/float(len(quality))
		stand=standerd_varition(length)
		length=sum(length)/len(length)
		readnames=''
		for c in candi:
			readnames+=c.split('\t')[4]+';'
		readnames=readnames[:-1]
		return [chrom+'\t'+str(position)+'\t'+str(length)+'\t'+str(len(candi))+'\t'+str(quality)+'\t'+str(stand)+'\t'+readnames]
	else:
		return []

def mergeinfo_deletion(candi):
	dels=[]
	candi.sort(key=mergeinfolengthsort)
	while len(candi)>2:
		if int(candi[-1].split('\t')[2]) > 1.5* int(candi[len(candi)/2].split('\t')[2]):
			candi.remove(candi[-1])
			continue
		if int(candi[len(candi)/2].split('\t')[2]) >  1.5*int(candi[0].split('\t')[2]):
			candi.remove(candi[0])
			continue
		break
	if len(candi)>2:
		chrom=candi[0].split('\t')[0]
		position=[int(c.split('\t')[1]) for c in candi]
		length=[int(c.split('\t')[2]) for c in candi]
		quality=[float(c.split('\t')[3]) for c in candi]
		position=sum(position)/len(position)
		quality=sum(quality)/float(len(quality))
		stand=standerd_varition(length)
		length=sum(length)/len(length)
		return[chrom+'\t'+str(position)+'\t'+str(length)+'\t'+str(len(candi))+'\t'+str(quality)+'\t'+str(stand)]
	else:
		return []

def counttimesort_tra(a):
	return [int(a.split('\t')[1]),int(a.split('\t')[3])]

def counttime_translocation(samechrom,min_support):
	samechrom.sort(key=counttimesort_tra)
	samechrtra=[]
	start=int(samechrom[0].split('\t')[1])
	candi=[]
	window=800
	for event in samechrom:
		if int(event.split('\t')[1])<=start+window:
			candi+=[event]
			continue
		if len(candi)>=min_support:
			samechrtra+=mergeinfo_translocation(candi,min_support)
		candi=[event]
		start=int(event.split('\t')[1])
	if len(candi)>=min_support:
		samechrtra+=mergeinfo_translocation(candi,min_support)
		candi=[]
	return samechrtra


def counttimesort(a):
	return [int(a.split('\t')[1]),int(a.split('\t')[2])]


def counttime_insertion(samechrom,min_support):
	if samechrom==[]:
		return []
	samechrom.sort(key=counttimesort)
	samechrins=[]
	start=int(samechrom[0].split('\t')[1])
	candi=[]
	inslength=[]
	window=100
	for event in samechrom:
		if int(event.split('\t')[1])<=start+window:
			candi+=[event]
			inslength+=[int(event.split('\t')[2])]
			continue
		if window==100:
			length=sum(inslength)/len(inslength)
			if length<=100:
				window=200
			if 100<length<=500:
				window=400
			if length>500:
				window=800
			if int(event.split('\t')[1])<=start+window:
				candi+=[event]
				inslength+=[int(event.split('\t')[2])]
				continue
		if len(candi)>=min_support:
			samechrins+=mergeinfo_insertion(candi)
		candi=[event]
		inslength=[int(event.split('\t')[2])]
		start=int(event.split('\t')[1])
		window=100
	if len(candi)>=min_support:
		samechrins+=mergeinfo_insertion(candi)
		candi=[]
	return samechrins



def counttime_deletion(samechrom,min_support):
	if samechrom==[]:
		return []
	samechrom.sort(key=counttimesort)
	samechrdel=[]
	start=int(samechrom[0].split('\t')[1])
	candi=[]
	dellength=[]
	window=100
	for event in samechrom:
		if int(event.split('\t')[1])<=start+window:
			candi+=[event]
			dellength+=[int(event.split('\t')[2])]
			continue
		if window==100:
			length=sum(dellength)/len(dellength)
			if length<=100:
				window=200
			if 100<length<=500:
				window=400
			if length>500:
				window=800
			if int(event.split('\t')[1])<=start+window:
				candi+=[event]
				dellength+=[int(event.split('\t')[2])]
				continue
		if len(candi)>=min_support:
			samechrdel+=mergeinfo_insertion(candi)
		candi=[event]
		dellength=[int(event.split('\t')[2])]
		start=int(event.split('\t')[1])
		window=100
	if len(candi)>=min_support:
		samechrdel+=mergeinfo_insertion(candi)
		candi=[]
	return samechrdel

def merge_deletion(min_support,min_quality,readpath,samechrom_deletion,chrom,svtype,upper_bound):
	delt1=time.time()
	samechrom_deletion=[c for c in samechrom_deletion if float(c.split('\t')[3])>=min_quality]
	if samechrom_deletion==[]:
		return True
	tt1=time.time()
	deletions=counttime_deletion(samechrom_deletion,3)
	f=open(readpath+svtype+'-info-'+chrom,'w')
	for d in deletions:
		f.write(d+'\n')
	f.close()
	real=deletions
	'''
	for d in deletions:
		co=int(d.split('\t')[3])
		if co>=min_support and (float(d.split('\t')[4])>= or int(d.split('\t')[3])>=min_support+5 or int(d.split('\t')[2])>=1000):
			real+=[d]
	'''
	cleaneddels=m_samechr_deletion(real)
	if upper_bound:
		cleaneddels=[c for c in cleaneddels if min_support<=int(c.split('\t')[3])<=min_support*30]
	else:
		cleaneddels=[c for c in cleaneddels if min_support<=int(c.split('\t')[3])]
	f=open(readpath+svtype+'-merged-'+chrom,'w')
	if svtype=='del':
		sv_type='Deletion'
	if svtype=='dup':
		sv_type='Duplication'
	if svtype=='inv':
		sv_type='Inversion'
	merged_result=[]
	for d in cleaneddels:
		f.write(d+'\t'+sv_type+'\n')
		merged_result+=[d+'\t'+sv_type]
	f.close()
	delt2=time.time()
	print(chrom+'-'+svtype+' done, time costed: '+str(delt2-delt1))
	return merged_result

def  merge_insertion(min_support,min_quality,readpath,samechrom_insertion,chrom,svtype,upper_bound):
	samechrom_insertion=[c for c in samechrom_insertion if float(c.split('\t')[3])>=min_quality]
	if samechrom_insertion==[]:
		return True

	inst1=time.time()
	insertions=counttime_insertion(samechrom_insertion,3)
	f=open(readpath+svtype+'-info-'+chrom,'w')
	for d in insertions:
		f.write(d+'\n')
	f.close()
	real=[]
	for d in insertions:
		co=int(d.split('\t')[3])
		if (co>=min_support or (co>=max(min_support-3,3) and int(d.split('\t')[2])>=2000 )) and (float(d.split('\t')[4])>10 or co>=min_support+5 or int(d.split('\t')[2])>=500):
			real+=[d]
	cleanedins=m_samechr_insertion(insertions)
	if upper_bound:
		cleanedins=[c for c in cleanedins if min_support<=int(c.split('\t')[3])<=30*min_support]
	else:
		cleanedins=[c for c in cleanedins if min_support<=int(c.split('\t')[3])]
	f=open(readpath+svtype+'-merged-'+chrom,'w')
	if svtype=='ins':
		sv_type='Insertion'
	if svtype=='inv':
		sv_type='Inversion'
	merged_result=[]
	for d in cleanedins:
		f.write(d+'\t'+sv_type+'\n')
		merged_result+=[d+'\t'+sv_type]
	f.close()
	inst2=time.time()

	print(chrom+'-'+svtype+' done, time costed: '+str(inst2-inst1))
	return merged_result

def finalsorttra(a):
	return [a.split('\t')[0],int(a.split('\t')[1])]

def merge_translocation(min_support,min_qual,readpath,samechrom_translocation,chrom,upper_bound):
	samechrom_translocation=[c for c in samechrom_translocation if int(c.split('\t')[4])>=min_qual]
	if samechrom_translocation==[]:
		return True
	trat1=time.time()
	translocations=counttime_translocation(samechrom_translocation,min_support)
	translocations=m_samechr_translocation(translocations)
	if upper_bound:
		translocations=[c for c in translocations if int(c.split('\t')[4])<=30*min_support]

	translocations.sort(key=finalsorttra)
	merged_result=[]
	f=open(readpath+'tra-merged-'+chrom,'w')
	for d in translocations:
		f.write(d+'\n')
		merged_result+=[d]
	f.close()
	trat2=time.time()
	print(chrom+'-tra done, time costed: '+str(trat2-trat1))
	return merged_result


if __name__ == "__main__":
	t1=time.time()
	readpath="/data/scratch/maggic/DeBreak_manuscript/simulation/debreak/"
	f=open(readpath+'filelist','r')
	filelist=f.read().split('\n')[:-1]
	f.close()
	chromosomes=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
	#chromosomes=['chr17']
	rawinscall={}
	rawdelcall={}
	rawdupcall={}
	rawinvcall={}
	rawtracall={}
	for c in chromosomes:
		rawinscall[c]=[]
		rawdelcall[c]=[]
		rawdupcall[c]=[]
		rawinvcall[c]=[]
		rawtracall[c]=[]
	for d in filelist:
		f=open(readpath+d,'r')
		allsvcall=f.read().split('\n')[:-1]
		
		for m in allsvcall:
			if 'I-' in m:
				m=m.split('\t')
				rawinscall[m[0]]+=[m[0]+'\t'+m[1]+'\t'+m[2]+'\t'+m[6]+'\t'+m[4]]
			if 'D-' in m:
				m=m.split('\t')
				rawdelcall[m[0]]+=[m[0]+'\t'+m[1]+'\t'+m[2]+'\t'+m[6]+'\t'+m[4]]
			if 'DUP-segment' in m:
				m=m.split('\t')
				rawdupcall[m[0]]+=[m[0]+'\t'+m[1]+'\t'+m[2]+'\t'+m[6]+'\t'+m[4]]
			if 'INV-segment' in m:
				m=m.split('\t')
				rawinvcall[m[0]]+=[m[0]+'\t'+m[1]+'\t'+m[2]+'\t'+m[6]+'\t'+m[4]]
			if 'TRA-segment' in m:
				m=m.split('\t')
				rawtracall[m[0]]+=[m[0]+'\t'+m[1]+'\t'+m[2]+'\t'+m[3]+'\t'+m[7]+'\t'+m[5]]
	supp=1
	for c in chromosomes:
		merge_insertion(supp,0,readpath,rawinscall[c],c,'ins',False)
		merge_deletion(supp,0,readpath,rawdelcall[c],c,'del',False)
		#merge_deletion(supp,0,readpath,rawdupcall[c],c,'dup',True)
		#merge_insertion(supp,0,readpath,rawinvcall[c],c,'inv',True)
		#merge_translocation(supp,0,readpath,rawtracall[c],c,True)
	
	os.system("cat "+readpath+"ins-merged-* > "+readpath+"insertion-merged")
	#os.system("rm "+ readpath +"ins-*")
	os.system("cat "+readpath+"del-merged-* > "+readpath+"deletion-merged")
	#os.system("rm "+ readpath +"del-*")
	'''
	os.system("cat "+readpath+"dup-merged-* > "+readpath+"duplication-merged")
	os.system("rm "+ readpath +"dup-*")
	os.system("cat "+readpath+"inv-merged-* > "+readpath+"inversion-merged")
	os.system("rm "+ readpath +"inv-*")
	os.system("cat "+readpath+"tra-merged-* > "+readpath+"translocation-merged")
	os.system("rm "+ readpath +"tra-*")
	'''
	t2=time.time()
	print 'Total time used:'+str(t2-t1)


