f=open('/data/scratch/maggic/skbr3/newseq/debreak-hg38/four-gt','r')
a=f.read().split('\n')[:-1]
f.close()


f=open('/data/scratch/maggic/skbr3/newseq/debreak-hg38/complex_candidate','w')
for i in range(len(a)-1):
	if a[i].split('\t')[0]!=a[i+1].split('\t')[0] or int(a[i].split('\t')[1])+int(a[i].split('\t')[2])<int(a[i+1].split('\t')[1]) or max(int(a[i].split('\t')[2]),int(a[i+1].split('\t')[2]))>10000:
		continue
	if a[i].split('\t')[7]==a[i+1].split('\t')[7]:
		continue
	if ('Insertion' in a[i] and 'Deletion' in a[i+1]) or ('Insertion' in a[i+1] and 'Deletion' in a[i]):
		continue
	if 'Insertion' not in a[i] and 'Insertion' not in a[i+1]:
		f.write(a[i].split('\t')[0]+'\t'+a[i].split('\t')[1]+'\t'+a[i].split('\t')[2]+'\t'+a[i].split('\t')[3]+'\t'+a[i].split('\t')[7]+'\t'+a[i].split('\t')[8]+'\n'+a[i+1].split('\t')[0]+'\t'+a[i+1].split('\t')[1]+'\t'+a[i+1].split('\t')[2]+'\t'+a[i+1].split('\t')[3]+'\t'+a[i+1].split('\t')[7]+'\t'+a[i+1].split('\t')[8]+'\n\n')
	if 'Insertion' in a[i] and 'Insertion' not in a[i+1]:
		continue
	if 'Insertion' not in a[i] and 'Insertion' in a[i+1]:
		f.write(a[i].split('\t')[0]+'\t'+a[i].split('\t')[1]+'\t'+a[i].split('\t')[2]+'\t'+a[i].split('\t')[3]+'\t'+a[i].split('\t')[7]+'\t'+a[i].split('\t')[8]+'\n'+a[i+1].split('\t')[0]+'\t'+a[i+1].split('\t')[1]+'\t'+a[i+1].split('\t')[2]+'\t'+a[i+1].split('\t')[3]+'\t'+a[i+1].split('\t')[7]+'\t'+a[i+1].split('\t')[8]+'\n\n')


f.close()
