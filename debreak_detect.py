import pysam
import time
import os

def cigardeletion(flag,chrom,position,cigar,min_size,max_size):	#input a read line, return list of deletions
	flag=int(flag)
	if flag<=16:
		detect_cigar_sv=True
	else:
		detect_cigar_sv=False
	pos=int(position)
	numbers='1234567890'
	num=''
	reflen=0
	readlen=0
	leftclip=0
	rightclip=0
	deletions=[]
	insertions=[]
	for c in cigar:
		if c in numbers:
			num+=c
			continue
		if c in 'MNP=X':
			readlen+=int(num); reflen+=int(num);  num='';  continue
		if c=='I':
			if detect_cigar_sv and int(num)>=min_size and int(num)<=max_size:
				insertions+=[[chrom,pos+reflen,int(num),'I-cigar',readlen+leftclip]]
			readlen+=int(num)
			num=''; continue
		if c == 'D':
			if detect_cigar_sv and  int(num)>=min_size and int(num)<=max_size:
				deletions+=[[chrom,pos+reflen,int(num),'D-cigar']]
			reflen+=int(num);  num='';  continue
		if c in 'SH':
			if readlen==0:
				leftclip=int(num)
			else:
				rightclip=int(num)
			num=''; continue
	#merge deletions
	if detect_cigar_sv:
		testif=1
		window=500
		while testif==1:
			testif=0
			if len(deletions)==1:
				break
			i=len(deletions)-1
			while i>0:
				gaplength=deletions[i][1]-deletions[i-1][1]-deletions[i-1][2]
				if gaplength <= window:
					deletions=deletions[:i-1]+[[chrom,deletions[i-1][1],deletions[i-1][2]+deletions[i][2],'D-cigar']]+deletions[i+1:]
					testif=1
					break
				else:
					i-=1
		#merge insertions
		testif=1
		while testif==1:
			testif=0
			if len(insertions)==1:
				break
			i=len(insertions)-1
			while i>0:
				l1=insertions[i][2]
				l2=insertions[i-1][2]
				gaplength=insertions[i][1]-insertions[i-1][1]
				window=200 if max(l1,l2)<100 else 400
				window=400 if window==400 and max(l1,l2) <500 else 600
				if gaplength >window :
					i-=1
				else:
					insertions=insertions[:i-1]+[[chrom,insertions[i-1][1],l1+l2,'I-cigar',insertions[i-1][4]]]+insertions[i+1:]
					testif=1
					break
	
	svcallset=deletions+insertions
	
	return [svcallset,reflen,[leftclip,readlen,rightclip]]


def simplifycigar(cigar):
	numbers='1234567890'
	num=''
	readlen=0
	leftclip=0
	rightclip=0
	for c in cigar:
		if c in numbers:
			num+=c; continue
		if c in 'MNP=XI':
			readlen+=int(num);   num='';  continue
		if c == 'D':
			 num=''; continue
		if c in 'SH':
			if readlen==0:
				leftclip=int(num)
				num=''
			else:
				rightclip=int(num)
				num=''

	return [leftclip,readlen,rightclip]


def segmentdeletion(segments,min_size,max_size):  #input a list of segments,return list of deletions
	if len([c for c in segments if int(c[1])<=16])==0:
		return []
	if len(segments)<=1:
		return []
	svcallset=[]
	primary=[c for c in segments if int(c[1])<=16][0]
	segments=[c for c in segments if c != primary]
	chrom=primary[2]
	priflag=(int(primary[1])%32)>15
	samedirchr=[]
	samechr=[]
	diffchr=[]
	rawsegsvid=1
	for c in segments:
		ch=c[2]
		f=int(c[1])%32>15
		if c[5][1]<300:
			continue
		if ch!=chrom:
			diffchr+=[c]
		elif f!=priflag:
			samechr+=[c]
		else:
			samedirchr+=[c]
	for c in samedirchr:
		if c[3]>primary[3] :
			leftread=primary
			rightread=c
		elif c[3]<primary[3]:
			leftread=c
			rightread=primary
		else:
			continue
		leftinfo=leftread[5]
		rightinfo=rightread[5]
		#insertion:
		if abs(rightread[3]-leftread[4])<=300:
			overlap=rightread[3]-leftread[4]
			ins_size=rightinfo[0]-leftinfo[1]-leftinfo[0]-overlap
			if min_size<=ins_size<=max_size:
				svcallset+=[chrom+'\t'+str(min(rightread[3],leftread[4]))+'\t'+str(ins_size)+'\t'+'I-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(int(c[1])+int(primary[1]))+'\t'+str((int(c[6])+int(primary[6]))//2)]
				rawsegsvid+=1

		#deletion:
		overlapmap=leftinfo[0]+leftinfo[1]-rightinfo[0]
		window_max=1500
		if -200<overlapmap<window_max:
			del_size=rightread[3]-leftread[4]+overlapmap
			if min_size<=del_size<=max_size:
				svcallset+=[chrom+'\t'+str(leftread[4]-max(0,overlapmap))+'\t'+str(del_size)+'\t'+'D-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
				rawsegsvid+=1

		#duplication:
		overlapmap=leftinfo[0]+leftinfo[1]-rightinfo[0]
		window_max=500
		if -200<overlapmap<window_max and leftread[4]-rightread[3]>=max(50,overlapmap):
			dup_size=leftread[4]-rightread[3]-max(overlapmap,0)
			if min_size<=dup_size<=max_size:
				svcallset+=[chrom+'\t'+str(rightread[3])+'\t'+str(dup_size)+'\t'+'DUP-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
				rawsegsvid+=1
		overlapmap=rightinfo[0]+rightinfo[1]-leftinfo[0]
		if -200<overlapmap<window_max and (rightread[4]-leftread[3])>=max(1000,overlapmap):
			dup_size=rightread[4]-leftread[3]-overlapmap
			if min_size<=dup_size<=max_size:
				svcallset+=[chrom+'\t'+str(leftread[3])+'\t'+str(dup_size)+'\t'+'DUP-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
				rawsegsvid+=1
	#inversion:
	for c in samechr:
		if c[3]>primary[3] and c[4]-primary[4]>-200:
			leftread=primary
			rightread=c
		elif c[3]<primary[3] and primary[4]-c[4]>-200:
			leftread=c
			rightread=primary
		else:
			continue
		leftinfo=leftread[5]
		rightinfo=rightread[5]
		window_max=500
		overlapmap=rightinfo[0]+rightinfo[1]-leftinfo[2]
		if -200<overlapmap<window_max and (rightread[4]-leftread[4])>=max(100,overlapmap):
			inv_size=rightread[4]-leftread[4]-overlapmap
			if min_size<=inv_size<=max_size:
				svcallset+=[chrom+'\t'+str(leftread[4])+'\t'+str(inv_size)+'\t'+'INV-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
				rawsegsvid+=1
				continue
		overlapmap=rightinfo[1]+rightinfo[2]-leftinfo[0]
		if -200<overlapmap<window_max and (rightread[3]-leftread[3])>=max(100,overlapmap):
			inv_size=rightread[3]-leftread[3]-overlapmap
			if min_size<=inv_size<=max_size:
				svcallset+=[chrom+'\t'+str(leftread[3])+'\t'+str(inv_size)+'\t'+'INV-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
				rawsegsvid+=1
				continue
	#translocation:
	for c in diffchr:
		pinfo=primary[5]
		cinfo=c[5]
		window_max=200
		bp1=''
		bp2=''
		if abs(pinfo[0]-cinfo[0]-cinfo[1])<=window_max or abs(pinfo[0]-cinfo[1]-cinfo[2])<=window_max:
			chrom1=primary[2]
			bp1=primary[3]
		elif abs(pinfo[2]-cinfo[0]-cinfo[1])<=window_max or abs(pinfo[2]-cinfo[1]-cinfo[2])<=window_max:
			chrom1=primary[2]
			bp1=primary[4]
		if abs(cinfo[0]-pinfo[0]-pinfo[1])<=window_max or abs(cinfo[0]-pinfo[1]-pinfo[2])<=window_max:
			chrom2=c[2]
			bp2=c[3]
		elif abs(cinfo[2]-pinfo[0]-pinfo[1])<=window_max or abs(cinfo[2]-pinfo[1]-pinfo[2])<=window_max:
			chrom2=c[2]
			bp2=c[4]
		if bp1!='' and bp2!='':
			if chrom1 > chrom2:
				svcallset+=[chrom2+'\t'+str(bp2)+'\t'+chrom1+'\t'+str(bp1)+'\tTRA-segment\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
				rawsegsvid+=1
			if chrom1 < chrom2:
				svcallset+=[chrom1+'\t'+str(bp1)+'\t'+chrom2+'\t'+str(bp2)+'\tTRA-segment\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
				rawsegsvid+=1

	return svcallset




def segmentdeletion_tumor(segments,min_size,max_size):  #input a list of segments,return list of deletions
	svcallset=[]
	for i in range(len(segments)-1):
		primary=segments[i]
		others=segments[i+1:]

		chrom=primary[2]
		priflag=(int(primary[1])%32)>15
		samedirchr=[]
		samechr=[]
		diffchr=[]
		rawsegsvid=1
		for c in others:
			ch=c[2]
			f=int(c[1])%32>15
			if c[5][1]<300:
				continue
			if ch!=chrom:
				diffchr+=[c]
			elif f!=priflag:
				samechr+=[c]
			else:
				samedirchr+=[c]

		for c in samedirchr:
			if c[3]>primary[3] :
				leftread=primary
				rightread=c
			elif c[3]<primary[3]:
				leftread=c
				rightread=primary
			else:
				continue	
			leftinfo=leftread[5]
			rightinfo=rightread[5]
		#insertion:
			if abs(rightread[3]-leftread[4])<=300:
				overlap=rightread[3]-leftread[4]
				ins_size=rightinfo[0]-leftinfo[1]-leftinfo[0]-overlap
				if min_size<=ins_size<=max_size:
					svcallset+=[chrom+'\t'+str(min(rightread[3],leftread[4]))+'\t'+str(ins_size)+'\t'+'I-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(int(c[1])+int(primary[1]))+'\t'+str((int(c[6])+int(primary[6]))//2)]
					rawsegsvid+=1

		#deletion:
			overlapmap=leftinfo[0]+leftinfo[1]-rightinfo[0]
			window_max=1500
			if -200<overlapmap<window_max:
				del_size=rightread[3]-leftread[4]+overlapmap
				if min_size<=del_size<=max_size:
					svcallset+=[chrom+'\t'+str(leftread[4]-max(0,overlapmap))+'\t'+str(del_size)+'\t'+'D-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
					rawsegsvid+=1
		
		#duplication:
			overlapmap=leftinfo[0]+leftinfo[1]-rightinfo[0]
			window_max=500
			if -200<overlapmap<window_max and leftread[4]-rightread[3]>=max(50,overlapmap):
				dup_size=leftread[4]-rightread[3]-max(overlapmap,0)
				if min_size<=dup_size<=max_size:
					svcallset+=[chrom+'\t'+str(rightread[3])+'\t'+str(dup_size)+'\t'+'DUP-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
					rawsegsvid+=1
			overlapmap=rightinfo[0]+rightinfo[1]-leftinfo[0]
			if -200<overlapmap<window_max and (rightread[4]-leftread[3])>=max(1000,overlapmap):
				dup_size=rightread[4]-leftread[3]-overlapmap
				if min_size<=dup_size<=max_size:
					svcallset+=[chrom+'\t'+str(leftread[3])+'\t'+str(dup_size)+'\t'+'DUP-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
					rawsegsvid+=1
	#inversion:
		for c in samechr:
			if c[3]>primary[3] and c[4]-primary[4]>-200:
				leftread=primary
				rightread=c
			elif c[3]<primary[3] and primary[4]-c[4]>-200:
				leftread=c
				rightread=primary
			else:
				continue
			leftinfo=leftread[5]
			rightinfo=rightread[5]
			window_max=500
			overlapmap=rightinfo[0]+rightinfo[1]-leftinfo[2]
			if -200<overlapmap<window_max and (rightread[4]-leftread[4])>=max(100,overlapmap):
				inv_size=rightread[4]-leftread[4]-overlapmap
				if min_size<=inv_size<=max_size:
					svcallset+=[chrom+'\t'+str(leftread[4])+'\t'+str(inv_size)+'\t'+'INV-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
					rawsegsvid+=1
					continue
			overlapmap=rightinfo[1]+rightinfo[2]-leftinfo[0]
			if -200<overlapmap<window_max and (rightread[3]-leftread[3])>=max(100,overlapmap):
				inv_size=rightread[3]-leftread[3]-overlapmap
				if min_size<=inv_size<=max_size:
					svcallset+=[chrom+'\t'+str(leftread[3])+'\t'+str(inv_size)+'\t'+'INV-segment'+'\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
					rawsegsvid+=1
					continue
	#translocation:
		for c in diffchr:
			pinfo=primary[5]
			cinfo=c[5]
			window_max=200
			bp1=''
			bp2=''
			if abs(pinfo[0]-cinfo[0]-cinfo[1])<=window_max or abs(pinfo[0]-cinfo[1]-cinfo[2])<=window_max:
				chrom1=primary[2]
				bp1=primary[3]
			elif abs(pinfo[2]-cinfo[0]-cinfo[1])<=window_max or abs(pinfo[2]-cinfo[1]-cinfo[2])<=window_max:
				chrom1=primary[2]
				bp1=primary[4]
			if abs(cinfo[0]-pinfo[0]-pinfo[1])<=window_max or abs(cinfo[0]-pinfo[1]-pinfo[2])<=window_max:
				chrom2=c[2]
				bp2=c[3]
			elif abs(cinfo[2]-pinfo[0]-pinfo[1])<=window_max or abs(cinfo[2]-pinfo[1]-pinfo[2])<=window_max:
				chrom2=c[2]
				bp2=c[4]
			if bp1!='' and bp2!='':
				if chrom1 > chrom2:
					svcallset+=[chrom2+'\t'+str(bp2)+'\t'+chrom1+'\t'+str(bp1)+'\tTRA-segment\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
					rawsegsvid+=1
				if chrom1 < chrom2:
					svcallset+=[chrom1+'\t'+str(bp1)+'\t'+chrom2+'\t'+str(bp2)+'\tTRA-segment\t'+primary[0]+'_seg'+str(rawsegsvid)+'\t'+str(c[1])+'\t'+str((int(c[6])+int(primary[6]))//2)]
					rawsegsvid+=1


	return svcallset



def detect_sam(filename,readpath,writepath,chromosomes,min_size,max_size,record_depth,if_rawsvcall):
	print('Start to detect SVs from '+filename)
	f=open(readpath+filename,'r')
	c=f.readline()
	if if_rawsvcall:
		g=open(writepath+'sv_raw_calls/'+filename[:-4]+'.debreak.temp','w')
	else:
		g=open(writepath+filename[:-4]+'.debreak.temp','w')

	lastname=''
	segments=['']
	totalmaplength=0
	while c!='':
		#remove headerlines, secondary alignments, alignment on scallfolds
		if c[0]=='@' or int(c.split('\t')[1])%8>3 or int(c.split('\t')[1])%512>255:
			c=f.readline()
			continue
		if chromosomes!=[] and c.split('\t')[2] not in chromosomes:
			c=f.readline();continue
		#detect the deletion from cigar 
		readname=c.split('\t')[0]
		flag=c.split('\t')[1]
		chrom=c.split('\t')[2]
		position=int(c.split('\t')[3])
		mappingquality=c.split('\t')[4]
		readseq=c.split('\t')[9]
		cigar=c.split('\t')[5]
		cigarinfo=cigardeletion(flag,chrom,position,cigar,min_size,max_size)
		refend=position+cigarinfo[1]
		cimplecigar=str(cigarinfo[2][0])+'\t'+str(cigarinfo[2][1])+'\t'+str(cigarinfo[2][2])
		# if primary: write deletions from cigar string
		if int(c.split('\t')[1])%4096<2048:
			cigarsv=cigarinfo[0]
			totalmaplength+=len(c.split('\t')[9])
			rawsvid=0
			for d in cigarsv:
				rawsvid+=1
				if 'I-cigar' in d:
					insertseq=readseq[d[4]:d[4]+d[2]]
					g.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'_cigar'+str(rawsvid)+'\t'+flag+'\t'+mappingquality+'\t'+insertseq+'\n')
				else:
					g.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'_cigar'+str(rawsvid)+'\t'+flag+'\t'+mappingquality+'\n')
		readinfo=[readname,flag,chrom,position,refend,cigarinfo[2],mappingquality]
		if readname!=lastname:
			if 1<len(segments)<=20:
				segmentd=segmentdeletion(segments,min_size,max_size)
				for d in segmentd:
					g.write(d+'\n')
			lastname=readname
			segments=[readinfo]
		else:
			segments+=[readinfo]
		c=f.readline()
	if 1<len(segments)<=20:
		segmentd=segmentdeletion(segments,min_size,max_size)
		for d in segmentd:
			g.write(d+'\n')
	segments=[]

	f.close()
	g.close()
	if record_depth and totalmaplength>0:
		f=open(writepath+'map_depth/maplength_'+filename,'w')
		f.write(str(totalmaplength)+'\n')
		f.close()	
	return True

def detect_sortbam(filename,writepath,min_size,max_size,chrom,chromosomes,record_clip,ifrescuedup,iftumor):
	print ('Start to detect SV from '+filename)
	samfile=pysam.AlignmentFile(filename,"rb")
	tempfile=open(writepath+'sv_raw_calls/'+filename.split('/')[-1]+'-'+chrom+'.debreak.temp','w')
	allreads=samfile.fetch(chrom)
	segmentreads={}
	tracandidate=[]
	clipinfo=[]
	totalmaplength=0
	for align in allreads:
		if align.is_unmapped or align.is_secondary:
			continue
		readname=align.query_name
		flag=align.flag
		position=align.reference_start+1
		mappingquality=align.mapping_quality
		cigar=align.cigar
		readclipinfo=chrom
		readstart=align.reference_start
		readstop=align.reference_end
		if record_clip and ( cigar[0][0] in [4,5] and cigar[0][1]>=200 ) or (cigar[-1][0] in [4,5] and  cigar[-1][1]>=200):
			if (cigar[0][0]==4 or cigar[0][0]==5) and  cigar[0][1]>=200:
				readclipinfo+='\t0\t'+str(readstart)
			else:
				readclipinfo+='\t'+str(readstart)+'\t0'
			if (cigar[-1][0]==4 or cigar[-1][0]==5) and  cigar[-1][1]>=200:
				readclipinfo+='\t0\t'+str(readstop)
			else:
				readclipinfo+='\t'+str(readstop)+'\t0'		
			clipinfo+=[readclipinfo]

		if align.is_supplementary:
			cigarinfo=[0,0,0]
			cigarleft=align.cigar[0]
			cigarright=align.cigar[-1]
			if cigarleft[0]==4 or cigarleft[0]==5:
				cigarinfo[0]=cigarleft[1]
			if cigarright[0]==4 or cigarright[0]==5:
				cigarinfo[2]=cigarright[1]
			cigarinfo[1]=align.query_alignment_length
			refend=align.reference_end+1
			readinfo=[readname,flag,chrom,int(position),refend,cigarinfo,mappingquality]
			if readname in segmentreads:
				segmentreads[readname]+=[readinfo]
			else:
				segmentreads[readname]=[readinfo]
		else:
			cigar=align.cigarstring
			readseq=align.query_sequence
			cigarinfo=cigardeletion(flag,chrom,position,cigar,min_size,max_size)
			cigarsv=cigarinfo[0]
			totalmaplength+=align.query_length
			rawsvid=0
			for d in cigarsv:
				rawsvid+=1
				if 'I-cigar' in d:
					insertseq=readseq[d[4]:d[4]+d[2]]
					if  ifrescuedup:
						rawseq=readseq[:d[4]]+readseq[d[4]+d[2]:]
					else:
						rawseq=''
					tempfile.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'_cigar'+str(rawsvid)+'\t'+str(flag)+'\t'+str(mappingquality)+'\t'+str(d[4])+'\t'+insertseq+'\t'+rawseq+'\n')
				else:
					tempfile.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'_cigar'+str(rawsvid)+'\t'+str(flag)+'\t'+str(mappingquality)+'\n')
					
				
			if align.has_tag("SA"):
				refend=position+cigarinfo[1]
				readinfo=[readname,flag,chrom,position,refend,cigarinfo[2],mappingquality]
				if readname in segmentreads:
					segmentreads[readname]+=[readinfo]
				else:
					segmentreads[readname]=[readinfo]
				suppaligns=align.get_tag("SA").split(';')[:-1]
				diffchr=[]
				for c in suppaligns:
					if c.split(',')[0]==chrom or c.split(',')[0] not in chromosomes:
						continue
					if c.split(',')[2]=='+':
						suppflag=2048
					else:
						suppflag=2064
					cigarinfo=simplifycigar(c.split(',')[3])
					suppinfo=[readname,suppflag,c.split(',')[0],int(c.split(',')[1]),int(c.split(',')[1])+cigarinfo[1],cigarinfo,c.split(',')[4]]
					segmentreads[readname]+=[suppinfo]
	for readgroup in segmentreads:
		if len(segmentreads[readgroup])<2 or len(segmentreads[readgroup])>20:
			continue
		if iftumor:
			segmentsv=segmentdeletion_tumor(segmentreads[readgroup],min_size,max_size)
		else:
			segmentsv=segmentdeletion(segmentreads[readgroup],min_size,max_size)
		for d in segmentsv:
			tempfile.write(d+'\n')
	tempfile.close()
	if record_clip:
		f=open(writepath+'debreak_ins_workspace/readinfo_start_end_'+chrom,'w')
		for c in clipinfo:
			f.write(c+'\n')
		f.close()
	print (chrom+' is done with maplength: '+str(totalmaplength))
	if totalmaplength!=0:
		f=open(writepath+'map_depth/maplength_'+chrom,'w')
		f.write(str(totalmaplength)+'\n')
		f.close()
	return True

				

