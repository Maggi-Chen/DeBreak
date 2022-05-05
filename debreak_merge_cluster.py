import os
import time
import sys

def mergeinfolengthsort(a):
	return int(a.split('\t')[2])

def mergeinfo_insertion(candi,min_support,svtypeinfo):
	candi.sort(key=mergeinfolengthsort)

	if len(candi)>=1.5*min_support:
		upper=int(candi[len(candi)*3//4].split('\t')[2])
		lower=int(candi[len(candi)//4].split('\t')[2])
		if upper>1.75*lower and upper-lower>50:
			svgroups=assign_candi_insertion(candi,upper,lower)
			svgroups=assign_candi_insertion(candi,svgroups[2],svgroups[3])
			svgroups=assign_candi_insertion(candi,svgroups[2],svgroups[3])
			mergedsv=[]
			if len(svgroups[0])>=min_support:
				mergedsv+=mergeinfo_insertion_oneevent(svgroups[0],min_support)
			if len(svgroups[1])>=min_support:
				mergedsv+=mergeinfo_insertion_oneevent(svgroups[1],min_support)
			if len(mergedsv)==2:
				mergedsv=[c+'\tCompoundSV\t'+svtypeinfo for c in mergedsv]
			if len(mergedsv)==1:
				mergedsv=[mergedsv[0]+'\tUnique\t'+svtypeinfo]
			return mergedsv
	mergedsv=mergeinfo_insertion_oneevent(candi,min_support)
	if len(mergedsv)==1:
		return [mergedsv[0]+'\tUnique\t'+svtypeinfo]
	else:
		return []

def assign_candi_insertion(candi,mean1,mean2):
	group1=[]
	group2=[]
	for c in candi:
		if abs(int(c.split('\t')[2])-mean1)<=abs(mean2-int(c.split('\t')[2])):
			group1+=[c]
		else:
			group2+=[c]
	mean1_new=sum([int(c.split('\t')[2]) for c in group1])//len(group1)
	mean2_new=sum([int(c.split('\t')[2]) for c in group2])//len(group2)
	return [group1,group2,mean1_new,mean2_new]


def mergeinfo_insertion_oneevent(candi,min_support):
	candi.sort(key=mergeinfolengthsort)
	while len(candi)>max(2,min_support-2):
		if int(candi[-1].split('\t')[2]) > 2* int(candi[len(candi)//2].split('\t')[2]) and  int(candi[-1].split('\t')[2]) -int(candi[len(candi)//2].split('\t')[2]) >30:
			candi.remove(candi[-1])
			continue
		if int(candi[len(candi)//2].split('\t')[2]) >  2*int(candi[0].split('\t')[2]) and int(candi[len(candi)//2].split('\t')[2]) -int(candi[0].split('\t')[2]) >30:
			candi.remove(candi[0])
			continue
		break
	if len(candi)>=max(2,min_support):
		chrom=candi[0].split('\t')[0]
		position=[int(c.split('\t')[1]) for c in candi]
		length=[int(c.split('\t')[2]) for c in candi]
		quality=[float(c.split('\t')[6]) for c in candi]
		position=sum(position)//len(position)
		quality=sum(quality)//float(len(quality))
		length=sum(length)//len(length)
		readnames=''
		for c in candi:
			readnames+=c.split('\t')[4]+';'
		readnames=readnames[:-1]
		numread=len(readnames.split(';'))
		return[chrom+'\t'+str(position)+'\t'+str(length)+'\t'+str(len(candi))+'\t'+str(numread)+'\t'+str(quality)+'\t'+readnames]
	else:
		return []

def counttimesort(a):
	return [int(a.split('\t')[1]),int(a.split('\t')[2])]

def cluster(outpath,allsv,chrom,contiglength,mins,maxdepth,svtype,minmapq):
	if allsv==[]:
		return 0
	print (chrom,len(allsv))
	t1=time.time()
	# Large DEL
	svtypeinfo={}
	svtypeinfo['del']='Deletion';svtypeinfo['dup']='Duplication';svtypeinfo['inv']='Inversion'
	svtypeinfo=svtypeinfo[svtype]
	largesv=[c for c in allsv if  int(c.split('\t')[2])>2000]
	window=1600
	largesv.sort(key=counttimesort)
	largedel=[]
	start=0
	candi=[]
	for event in largesv:
		if int(event.split('\t')[1])<=start+window:
			candi+=[event]
			continue
		if len(candi)>=mins:
			largedel+=mergeinfo_insertion(candi,mins,svtypeinfo)
		candi=[event]
		start=int(event.split('\t')[1])
	if len(candi)>=mins:
		largedel+=mergeinfo_insertion(candi,mins,svtypeinfo)
		candi=[]



	#smaller DEL
	allsv=[c for c in allsv if int(c.split('\t')[2])<=3000]


	genomeposition=[0]*contiglength

	for c in allsv:
		start=int(c.split('\t')[1])
		end=int(c.split('\t')[1])+int(c.split('\t')[2])
		original=genomeposition[start-1:end-1]
		new=[mm+1 for mm in original]
		genomeposition[start-1:end-1]=new
	
	svregion=[]
	inblock=False
	threshold=3

	for i in range(len(genomeposition)):
		if inblock:
			if genomeposition[i]>=max(maxdep/10.0,threshold):
				localdep+=[genomeposition[i]]
				if genomeposition[i]>maxdep:
					maxdep=genomeposition[i]
			else:
				inblock=False
				end=i
				if maxdep<=maxdepth:
					peakpos=localdep.index(maxdep)
					peakleftsize=0
					for i in range(peakpos):
						if localdep[peakpos-i-1]>=maxdep/10.0:
							peakleftsize+=1
						else:
							break
					svregion+=[(start+peakpos-peakleftsize,end,maxdep)]
				start=0
				end=0
				maxdep=0

		else:
			if genomeposition[i] > threshold:
				inblock=True
				localdep=[genomeposition[i]]
				start=i
				maxdep=genomeposition[i]
	#return 0

	svregion=[c for c in svregion if c[2] < maxdepth]
	allsvinfo={}
	for c in svregion:
		allsvinfo[c]=[]

	for c in allsv:
		start=int(c.split('\t')[1])
		end=start+int(c.split('\t')[2])
		for d in svregion:
			if min(end,d[1])-max(d[0],start)>0:
				allsvinfo[d]+=[c]

	sv=[]
	for c in svregion:
		svinfo=allsvinfo[c]
		sv+=mergeinfo_insertion(svinfo,mins,svtypeinfo)

	newsv=[]
	for c in largedel:
		testif=0
		for d in sv:
			if min(int(c.split('\t')[1])+int(c.split('\t')[2]), int(d.split('\t')[1])+int(d.split('\t')[2])) - max(int(c.split('\t')[1]),int(d.split('\t')[1]))>0 and 0.8*int(d.split('\t')[2])<=int(c.split('\t')[2])<=int(d.split('\t')[2])/0.8:
				testif=1; break
		if testif==0:
			newsv+=[c]
	newsv+=sv
	newsv.sort(key=counttimesort)

	t2=time.time()
	print (chrom,'-',svtype,' finished, time cost: ',t2-t1)

	if newsv!=[]:
		f=open(outpath+svtype+'-merged-cluster-'+chrom,'w')
		for c in newsv:
			if minmapq==None:
				if float(c.split('\t')[5])<1 or float(c.split('\t')[5])>=10:
					f.write(c+'\n')
			else:
				if float(c.split('\t')[5])>=minmapq:
					f.write(c+'\n')
		f.close()


	return 0



def cluster_ins(outpath,allsv,chrom,contiglength,mins,maxdepth,minmapq):
	if allsv==[]:
		return 0
	print (chrom,len(allsv))
	t1=time.time()
	svtypeinfo='Insertion'
	# Large INS
	largesv=[c for c in allsv if  int(c.split('\t')[2])>2000]
	window=1600
	largesv.sort(key=counttimesort)
	largedel=[]
	start=0
	candi=[]
	for event in largesv:
		if int(event.split('\t')[1])<=start+window:
			candi+=[event]
			continue
		if len(candi)>=mins:
			largedel+=mergeinfo_insertion(candi,mins,svtypeinfo)
		candi=[event]
		start=int(event.split('\t')[1])
	if len(candi)>=mins:
		largedel+=mergeinfo_insertion(candi,mins,svtypeinfo)
		candi=[]
	
	# Small INS
	allsv=[c for c in allsv if  int(c.split('\t')[2])<=3000]

	print (len(allsv))

	genomeposition=[0]*contiglength

	for c in allsv:
		start=int(c.split('\t')[1])-100
		end=int(c.split('\t')[1])+100
		original=genomeposition[start-1:end-1]
		new=[mm+1 for mm in original]
		genomeposition[start-1:end-1]=new
	svregion=[]
	inblock=False
	threshold=3

	for i in range(len(genomeposition)):
		if inblock:
			if genomeposition[i]>=max(maxdep/10.0,threshold):
				localdep+=[genomeposition[i]]
				if genomeposition[i]>maxdep:
					maxdep=genomeposition[i]
			else:
				inblock=False
				end=i
				if maxdep<=maxdepth:
					peakpos=localdep.index(maxdep)
					peakleftsize=0
					for i in range(peakpos):
						if localdep[peakpos-i-1]>=maxdep/10.0:
							peakleftsize+=1
						else:
							break
					svregion+=[(start+peakpos-peakleftsize,end,maxdep)]
				start=0
				end=0
				maxdep=0

		else:
			if genomeposition[i] > threshold:
				inblock=True
				localdep=[genomeposition[i]]
				start=i
				maxdep=genomeposition[i]


	svregion=[c for c in svregion if c[2] < maxdepth]
	allsvinfo={}
	for c in svregion:
		allsvinfo[c]=[]

	for c in allsv:
		start=int(c.split('\t')[1])-50
		end=start+100
		for d in svregion:
			if min(end,d[1])-max(d[0],start)>0:
				allsvinfo[d]+=[c]

	sv=[]
	for c in svregion:
		svinfo=allsvinfo[c]
		mergedins=mergeinfo_insertion(svinfo,mins,svtypeinfo)
		for m in mergedins:
			sv+=[m+'\t'+chrom+'\t'+str(c[0])+'\t'+str(c[1])+'\t'+str(c[2])]

	newsv=[]
	for c in largedel:
		testif=0
		for d in sv:
			if min(int(c.split('\t')[1])+int(c.split('\t')[2]), int(d.split('\t')[1])+int(d.split('\t')[2])) - max(int(c.split('\t')[1]),int(d.split('\t')[1]))>0 and 0.8*int(d.split('\t')[2])<=int(c.split('\t')[2])<=int(d.split('\t')[2])/0.8:
				testif=1; break
		if testif==0:
			newsv+=[c]
	newsv+=sv
	newsv.sort(key=counttimesort)

	t2=time.time()
	print (chrom,'-ins finished, time cost: ',t2-t1)

	if newsv!=[]:
		f=open(outpath+'ins-merged-cluster-'+chrom,'a')
		for c in newsv:
			if minmapq==None:
				if float(c.split('\t')[5])<1 or float(c.split('\t')[5])>=10:
					f.write(c+'\n')
			else:
				if float(c.split('\t')[5])>=minmapq:
					f.write(c+'\n')
		f.close()



	return 0



