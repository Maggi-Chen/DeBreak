def findgenename(chrom,pos,geneindex):
	regions=[]
	genes=[]
	for c in geneindex:
		if chrom==c.split('\t')[0] and int(c.split('\t')[1])<=pos<=int(c.split('\t')[2]):
			regions=geneindex[c]
	for c in regions:
		if int(c.split('\t')[3])<=pos<=int(c.split('\t')[4]) and c.split('\t')[8].split('"')[5] not in genes:
			genes+=[[c.split('\t')[8].split('"')[5],c.split('\t')[6]]]
	return genes

def findgenename_region(chrom,pos,end,geneindex):
	regions=[]
	genes=[]
	for c in geneindex:
		if chrom==c.split('\t')[0] and min(int(c.split('\t')[2]),end)-max(pos,int(c.split('\t')[1]))>=0:
			regions+=geneindex[c]
	for c in regions:
		if min(end,int(c.split('\t')[4]))-max(pos,int(c.split('\t')[3]))>=0 and c.split('\t')[8].split('"')[5] not in genes:
			genes+=[[c.split('\t')[8].split('"')[5],c.split('\t')[6]]]

	return genes


def annotate(svfile,annotafile,chromosomes,genefusion):
	alldel=open(svfile,'r').read().split('\n')[:-1]
	gene=open(annotafile,'r').read().split('\n')[:-1]
	gene=['chr'+c for c in gene if c[0]!='#' and c.split('\t')[2]=='gene']
	#index gene gtf
	geneindex={}
	for chrom in chromosomes: 
		genechr=[c for c in gene if c.split('\t')[0]==chrom]
		for i in range((len(genechr)+999)/1000):
			minim=genechr[1000*i].split('\t')[3]
			maxm=genechr[min(1000*(i+1),len(genechr))-1].split('\t')[4]
			info=chrom+'\t'+minim+'\t'+maxm
			regions=genechr[1000*i:min(1000*(i+1),len(genechr))]
			geneindex[info]=regions

	f=open(svfile+'-annotation','w')
	if genefusion:
		g=open(svfile+'-genefusion','w')
		genefucandi=[]
	for c in alldel:
		if 'Insertion' in c:
			chrom=c.split('\t')[0]
			pos=int(c.split('\t')[1])
			posgenename=findgenename(chrom,pos,geneindex)
			genenames=''
			if posgenename!=[]:
				for cc in posgenename:
					genenames+=cc[0]+','
				genenames=genenames[:-1]
				f.write(c+'\tgene_name='+genenames+'\n')
			else:
				f.write(c+'\t'+'gene_name=.\n')
			continue
		if 'Translocation' in c:
			chr1=c.split('\t')[0]
			pos1=int(c.split('\t')[1])
			pos1genename=findgenename(chr1,pos1,geneindex)
			chr2=c.split('\t')[2]
			pos2=int(c.split('\t')[3])
			pos2genename=findgenename(chr2,pos2,geneindex)
			genenames=''
			posgenename=pos1genename+pos2genename
			if posgenename!=[]:
				for cc in posgenename:
					genenames+=cc[0]+','
				genenames=genenames[:-1]
				f.write(c+'\tgene_name='+genenames+'\n')
			else:
				f.write(c+'\t'+'gene_name=.\n')
			
			if genefusion:
				pos1genename=[c for c in pos1genename if '.' not in c[0]]
				pos2genename=[c for c in pos2genename if '.' not in c[0]]
				if pos1genename!=[] and pos2genename!=[]:
					genefucandi+=[[chr1+'\t'+str(pos1)+'\t'+chr2+'\t'+str(pos2)+'\tTranslocation',pos1genename,pos2genename]]
			continue

		chrom=c.split('\t')[0]
		pos=int(c.split('\t')[1])
		end=int(c.split('\t')[2])+pos
		posgenename=findgenename_region(chrom,pos,end,geneindex)
		genenames=''
		if posgenename!=[]:
			for cc in posgenename:
				genenames+=cc[0]+','
			genenames=genenames[:-1]
			f.write(c+'\tgene_name='+genenames+'\n')
		else:
			f.write(c+'\t'+'gene_name=.\n')
		if genefusion and 'Deletion' in c and end-pos<=10000000:
			pos1genename=findgenename(chrom,pos,geneindex)
			pos1genename=[c for c in pos1genename if '.' not in c[0]]
			pos2genename=findgenename(chrom,end,geneindex)
			pos2genename=[c for c in pos2genename if '.' not in c[0]]
			if pos1genename!=[] and pos2genename!=[] and pos1genename!=pos2genename:
				genefucandi+=[[chrom+'\t'+str(pos)+'\t'+str(end-pos)+'\tDeletion',pos1genename,pos2genename]]

		if genefusion and 'Inversion' in c :
			pos1genename=findgenename(chrom,pos,geneindex)
			pos1genename=[c for c in pos1genename if '.' not in c[0]]
			pos2genename=findgenename(chrom,end,geneindex)
			pos2genename=[c for c in pos2genename if '.' not in c[0]]
			if pos1genename!=[] and pos2genename!=[] and pos1genename!=pos2genename:
				genefucandi+=[[chrom+'\t'+str(pos)+'\t'+str(end-pos)+'\tInversion',pos1genename,pos2genename]]

	f.close()
	if genefusion:
		for c in genefucandi:
			svcandi=c[0]
			pos1genename=c[1]
			pos2genename=c[2]
			pos1genename1=[c for c in pos1genename if c[1]=='+']
			pos1genename2=[c for c in pos1genename if c[1]=='-']
			pos2genename1=[c for c in pos2genename if c[1]=='+']
			pos2genename2=[c for c in pos2genename if c[1]=='-']
			if pos1genename1!=[] and pos2genename1!=[] and pos1genename1!=pos2genename1:
				if len(pos1genename1)!=1:
					pos1genename1=[c for c in pos1genename1 if c not in pos2genename1]
				if len(pos2genename1)!=1:
					pos2genename1=[c for c in pos2genename1 if c not in pos1genename1]
				genes1='';  genes2=''
				for c in pos1genename1:
					genes1+=c[0]+','
				for c in pos2genename1:
					genes2+=c[0]+','
				genes1=genes1[:-1];  genes2=genes2[:-1]
				g.write(svcandi+'\t'+genes1+';'+genes2+';+\n')
			if pos1genename2!=[] and pos2genename2!=[] and pos1genename2!=pos2genename2:
				if len(pos1genename2)!=1:
					pos1genename2=[c for c in pos1genename2 if c not in pos2genename2]
				if len(pos2genename2)!=1:
					pos2genename2=[c for c in pos2genename2 if c not in pos1genename2]
				genes1='';  genes2=''
				for c in pos1genename2:
					genes1+=c[0]+','
				for c in pos2genename2:
					genes2+=c[0]+','
				genes1=genes1[:-1];  genes2=genes2[:-1]
				g.write(svcandi+'\t'+genes1+';'+genes2+';-\n')
		g.close()


if __name__ =="__main__":
	svfile='/data/scratch/maggic/HG00733/debreak-pacbio/debreak-allsv-merged-final'
	annotafile='/data/user/maggic/svstudy/data/reference/Homo_sapiens.GRCh38.96.gtf'
	chromosomes=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
	annotate(svfile,annotafile,chromosomes,True)
	#genefusion(svfile,annotafile)

