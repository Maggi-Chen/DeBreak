import time
def if_overlap(a,b):
	window=100
	if 'Deletion' in a or 'Duplication' in a:
		chr1=a.split('\t')[0]
		pos1=int(a.split('\t')[1])	
		end1=pos1+int(a.split('\t')[2])

		if 'Deletion' in b or 'Duplication' in b:
			chr2=b.split('\t')[0]		
			pos2=int(b.split('\t')[1])
			end2=pos2+int(b.split('\t')[2])
			if chr1==chr2 and min(end1,end2)-max(pos1,pos2)>-window:
				return True
			else:
				return False


		if 'Insertion' in b:
			chr2=b.split('\t')[0]
			pos2=int(b.split('\t')[1])
			if chr1==chr2 and pos1-window<pos2<end1+window:
				return True
			else:
				return False

		if 'Inversion' in b:
			chr2=b.split('\t')[0]
			pos2=int(b.split('\t')[1])
			end2=pos2+int(b.split('\t')[2])
			if chr1==chr2 and (pos1-window<pos2<end1+window or pos1-window<end2<end1+window):
				return True
			else:
				return False


	if 'Insertion' in a:
		chr1=a.split('\t')[0]
		pos1=int(a.split('\t')[1])

		if 'Deletion' in b or 'Duplication' in b:
			chr2=b.split('\t')[0]
			pos2=int(b.split('\t')[1])
			end2=pos2+int(b.split('\t')[2])
			if chr1==chr2 and pos2-window<pos1<=end2+window:
				return True
			else:
				return False
		if 'Insertion' in b:
			chr2=b.split('\t')[0]
			pos2=int(b.split('\t')[1])
			if chr1==chr2 and abs(pos1-pos2)<window:
				return True
			else:
				return False
		if 'Inversion' in b:
			chr2=b.split('\t')[0]
			pos2=int(b.split('\t')[1])
			end2=pos2+int(b.split('\t')[2])
			if chr1==chr2 and (abs(pos1-pos2)<window or abs(pos1-end2)<window):
				return True
			else:
				return False


	if 'Inversion' in a:
		chr1=a.split('\t')[0]
		pos1=int(a.split('\t')[1])
		end1=pos1+int(a.split('\t')[2])
		if 'Deletion' in b or 'Duplication' in b:
			chr2=b.split('\t')[0]
			pos2=int(b.split('\t')[1])
			end2=pos2+int(b.split('\t')[2])
			if chr1==chr2 and  ( pos2-window<pos1<end2+window or pos2-window < end1 < end2+window):
				return True
			else:
				return False
		if 'Insertion' in b:
			chr2=b.split('\t')[0]
			pos2=int(b.split('\t')[1])
			if chr1==chr2 and ( abs(pos1-pos2)<window or abs(end1-pos2)<window ):
				return True
			else:
				return False
		if 'Inversion' in b:
			chr2=b.split('\t')[0]
			pos2=int(b.split('\t')[1])
			end2=pos2+int(b.split('\t')[2])
			if chr1==chr2 and (abs(pos1-pos2)<window or abs(end1-end2)<window):
				return True
			else:
				return False

def if_tra_not_overlap(a,allsv):
	window=100
	chr1=a.split('\t')[0]
	chr2=a.split('\t')[2]
	pos1=int(a.split('\t')[1])
	pos2=int(a.split('\t')[3])
	testif=0
	for c in allsv:
		if c.split('\t')[0]!=chr1 and c.split('\t')[0]!=chr2:
			continue
		chr3=c.split('\t')[0]
		pos3=int(c.split('\t')[1])
		end3=pos3+int(c.split('\t')[2])
		if 'Deletion' in c or 'Duplication' in c:
			if (chr1==chr3 and pos3-window<pos1<end3+window) or ( chr2==chr3 and pos3-window<pos2<end3+window):
				print a
				print c
				testif=1; break
		if 'Insertion' in c:
			if (chr1==chr3 and abs(pos3-pos1)<window) or(chr2==chr3 and  abs(pos3-pos2)<window):
				testif=1; break
		if 'Inversion' in c:
			if (chr1==chr3 and (abs(pos3-pos1)<window or abs(end3-pos1)<window)) or (chr2==chr3 and (abs(pos2-pos3)<window or abs(end3-pos2)<window)):
				testif=1; bireak
	if testif==1:
		return False
	else:
		return True


def sortpos(a):			
		return [a.split('\t')[0], int(a.split('\t')[1])]

def clean_all_new(alldel,allins,alldup,allinv,alltra):
	allsv=alldel+allins+alldup+allinv
	allsv.sort(key=sortpos)
	badsv=[]
	a=allsv
	for i in range(len(allsv)-1):
		if if_overlap(a[i],a[i+1]):
			if ('Insertion' in a[i] and 'Duplication' in a[i+1]) or int(a[i].split('\t')[3])<int(a[i+1].split('\t')[3]):
				badsv+=[a[i]]
				continue
			if ('Insertion' in a[i+1] and 'Duplication' in a[i]) or int(a[i].split('\t')[3])>=int(a[i+1].split('\t')[3]):
				badsv+=[a[i+1]]		
				continue
	allsv=[c for c in allsv if c not in badsv]
	goodtra=[]
	print len(alltra)
	for c in alltra:
		if if_tra_not_overlap(c,allsv):
			goodtra+=[c]
	allsv+=goodtra
	return allsv

if __name__ == "__main__":
	t1=time.time()
	readpath='/data/scratch/maggic/summary-results/debreak-1.0-3.10.19/simulation/all-type-new/lowerror/minimap2/'
	allins=open(readpath+'insertion-merged-new','r').read().split('\n')[:-1]
	alldel=open(readpath+'deletion-merged','r').read().split('\n')[:-1]
	allinv=open(readpath+'inversion-merged','r').read().split('\n')[:-1]
	alldup=open(readpath+'duplication-merged-new-add','r').read().split('\n')[:-1]
	alltra=open(readpath+'translocation-merged','r').read().split('\n')[:-1]

	allsv=clean_all_new(alldel,allins,alldup,allinv,alltra)

	f=open(readpath+'allsv','w')
	for c in allsv:
		f.write(c+'\n')
	f.close()
	t2=time.time()
	print t2-t1
