#!/usr/mbu/python/x86_64/bin/python
# SPIKING DATA WITH A HOMOZYGOUS DELETION TO GET SOME IDEA OF THE DISTRIBUTION OF POSSIBLE HOMOZYGOUS DELETIONS
# NOT ACTUALLY LIKELY THAT THERE WILL BE MANY OF THESE IN THE TRUE DATA IN ANY INTERESTING GENES

from intermine.webservice import Service
import math, sys, random
sys.path.append('/data/bioinf/home/cls87/Modules')
import NGSStats

# ----------------------------------------------------------------------------------------------------------------
# OPEN INTERVAL SUMMARY DATA FOR HEADERS
intervalSummary = open('/data/bioinf/home/cls87/DOCAnalysis/pt1/pt1DOCAnalysis.sample_interval_summary').readlines()
headers = intervalSummary[0]

# ----------------------------------------------------------------------------------------------------------------
# IDENTIFY SAMPLES AND CREATE A METHOD OF ORDERING NUMERICALLY
sHeaders = headers.split('\t')
sampleHeaders = sHeaders[3:]

# LIST KEEPS THE ORIGINAL ORDERS OF THE SAMPLES
sampleList = []
# DICTIONARY ASSOCIATES SAMPLE WITH POSITION IN THE INTVERAL SUMMARY FILE
sampleDict = {}

for a in sampleHeaders:
	sampleNo = a[6:8]
	meanCovName = 'NGS000'+sampleNo+'_mean_cvg'
	
	if sampleNo not in sampleDict:
		indexSample = sHeaders.index(meanCovName)
		sampleDict[sampleNo] = indexSample
		sampleList.append(sampleNo)

sortSampleList = sorted(sampleList)

# INDEXES OF THE SAMPLE COLUMN POSITIONS, IN SAMPLE NUMERICAL ORDER
orderList = []
for b in sortSampleList:
	orderList.append(sampleDict[b])

# ----------------------------------------------------------------------------------------------------------------
# READ IN UNIVERSAL VARIATION DATA
uniData = open('/data/bioinf/home/cls87/OutputAnalysis/SeparateGroups/SeparateMedians.txt').readlines()

varDict = {}

for line in uniData:
	line = line[:-1]
	sLine = line.split()
	# Sample number
	sampleNumber = sLine[0]
	# Median, SD
	varDict[sampleNumber] = [sLine[1], sLine[2]]

# ----------------------------------------------------------------------------------------------------------------
# BUILD A DICTIONARY TO MAKE IT POSSIBLE TO EXTRACT GENE NAME FROM INTERVAL POSITION
geneRefFile = open('/data/processing/PostProcessing/Databases/refGene/refGene_sorted.txt').readlines()

intervalRefDict = {}

for e in geneRefFile:
	# EXTRACTING INFORMATION ABOUT EACH GENE
	e = e[:-1]
	se = e.split('\t')
	geneName = se[12]
	chromosome = se[2]
	exonStarts = se[9]
	sExonStarts = (exonStarts.split(','))[:-1]
	exonEnds = se[10]
	sExonEnds = (exonEnds.split(','))[:-1]

	if chromosome not in intervalRefDict:
		intervalRefDict[chromosome] = {}
	
	# FOR EACH EXON START/END, CHECK TO SEE IF ITS IN THE DICTIONARY
	# THEN ASSOCIATE A GENE NAME WITH THIS INFORMATION
	# THERE CAN BE MULTIPLE GENE NAMES ASSOCIATED WITH AN EXON
	for f in sExonStarts:
		if f not in intervalRefDict[chromosome]:
			intervalRefDict[chromosome][f] = [geneName]
		else:
			if geneName not in intervalRefDict[chromosome][f]:
				intervalRefDict[chromosome][f].append(geneName)

	for g in sExonEnds:
		if g not in intervalRefDict[chromosome]:
			intervalRefDict[chromosome][g] = [geneName]
		else:
			if geneName not in intervalRefDict[chromosome][g]:
				intervalRefDict[chromosome][g].append(geneName)

# ---------------------------------------------------------------------------------------------------------------
# DEFINE THE SAMPLE GROUPS
NGS2 = ['04','05','06','07','08','09']
NGS3 = ['11','12','13','14','15']

# ----------------------------------------------------------------------------------------------------------------
# CREATE OUTPUT FILES
NGS2Output = open('/data/bioinf/home/cls87/SpikingData/HomDel/NGS2HomDel.txt','w')
NGS3Output = open('/data/bioinf/home/cls87/SpikingData/HomDel/NGS3HomDel.txt','w')

NGS2Output.write('Mito Ranking	Protein	Protein Description	Gene	Chromosome	Interval Start	Interval End	Average Coverage	Identified Sample	Identified Value	Difference	NGS00004 Coverage	NGS00005 Coverage	NGS00006 Coverage	NGS00007 Coverage	NGS0008 Coverage	NGS00009 Coverage	NGS00011 Coverage	NGS00012 Coverage	NGS00013 Coverage	NGS00014 Coverage	NGS00015 Coverage	NGS00004 Score	NGS00005 Score	NGS00006 Score	NGS00007 Score	NGS00008 Score	NGS00009 Score	NGS00011 Score	NGS00012 Score	NGS00013 Score	NGS00014 Score	NGS00015 Score\n')

NGS3Output.write('Mito Ranking	Protein	Protein Description	Gene	Chromosome	Interval Start	Interval End	Average Coverage	Identified Sample	Identified Value	Difference	NGS00004 Coverage	NGS00005 Coverage	NGS00006 Coverage	NGS00007 Coverage	NGS0008 Coverage	NGS00009 Coverage	NGS00011 Coverage	NGS00012 Coverage	NGS00013 Coverage	NGS00014 Coverage	NGS00015 Coverage	NGS00004 Score	NGS00005 Score	NGS00006 Score	NGS00007 Score	NGS00008 Score	NGS00009 Score	NGS00011 Score	NGS00012 Score	NGS00013 Score	NGS00014 Score	NGS00015 Score\n')

# SET UP CONNECTION TO THE MITOMINER DATABASE
TOKEN = 'T1V4d694f3z7o301H0q9'
service = Service('http://mitominer.mrc-mbu.cam.ac.uk/release-3.1/service', token = TOKEN)
# ----------------------------------------------------------------------------------------------------------------
# SET UP SOEM LISTS FOR DATA STORAGE
NGS2CovList = []
NGS3CovList = []

NGS2SpikedList = []
NGS3SpikedList = []

NGS2VarList = []
NGS3VarList = []

# FOR EACH PART (WHICH WAS ANALYSED SEPARATELY IN DOC ANALYSIS)
for partNo in range(1,17):
	# OPEN INTERVAL SUMMARY DATA FOR THAT PART
	intervalSummary = open('/data/bioinf/home/cls87/DOCAnalysis/pt'+str(partNo)+'/pt'+str(partNo)+'DOCAnalysis.sample_interval_summary').readlines()
	intervalSummary = intervalSummary[1:]

	for d in intervalSummary:
		# EXTRACTING INFORMATION ABOUT THE CHOSEN INTERVAL
		d = d[:-1]
		sd = d.split('\t')
		interval = sd[0]
		sInterval = interval.split(':')
		chromosome = sInterval[0]
		startEnd = sInterval[1]
		sStartEnd = startEnd.split('-')
		start = sStartEnd[0]
		end = sStartEnd[1]
	
		# SET UP MORE DATA STORAGE LISTS, AT AN INTERVAL LEVEL
		NGS2intervalList = []
		NGS2intervalOrder = []
		NGS2covList = []
	
		NGS3intervalList = []
		NGS3intervalOrder = []
		NGS3covList = []

		# FOR EACH SAMPLE IN NUMERICAL ORDER
		for h in orderList:
			# IDENTIFYING THE POSITION OF THE SAMPLE, SO CAN OUTPUT IN NUMERICAL ORDER
			whichSample = orderList.index(h)
			numWhichSam = whichSample + 1
			strWhichSam = str(numWhichSam)
			if len(strWhichSam) == 1:
				strWhichSam = '0'+strWhichSam
			else:
				strWhichSam = str(int(strWhichSam) + 1)

			# COVERAGE FOR THE INTERVAL FOR THE SPECIFIC SAMPLE
			intSamMeanCov = float(sd[int(h)])

			# DISCARD NGS1
			if (strWhichSam in NGS2) or (strWhichSam in NGS3):
				# EXTRACTING MEDIAN AND SD FOR THE SPECIFIC SAMPLE
				samMed = float(varDict[strWhichSam][0])
				samDiff = float(varDict[strWhichSam][1])
		
				# CALCULATE THE STANDARDISED SCORE FOR THE INTERVAL RELATIVE TO SAMPLE COVERAGE
				score = round((intSamMeanCov - samMed)/samDiff, 4)	

			# ADD SOME DATA TO STORAGE
			if strWhichSam in NGS2:
				NGS2intervalList.append(score)
				NGS2covList.append(intSamMeanCov)
				NGS2intervalOrder.append(strWhichSam)

			elif strWhichSam in NGS3:
				NGS3intervalList.append(score)
				NGS3covList.append(intSamMeanCov)
				NGS3intervalOrder.append(strWhichSam)

		# CALCULATE THE SCORE MEDIAN FOR EACH SAMPLE GROUP
		NGS2median = NGSStats.median(NGS2intervalList)
		NGS3median = NGSStats.median(NGS3intervalList)
		
		# ONLY IF AVERAGE COVERAGE FOR THE INTERVAL FOR THE SAMPLE GROUP IS OVER 10
		if NGSStats.mean(NGS2covList) > 10:
			NGS2intDiff = []
			NGS3intDiff = []

			# CALCULATE THE AVERAGE COVERAGE FOR THE INTERVAL, WITHIN THE SAMPLE GROUP
			avgCov = NGSStats.mean(NGS2covList)
			NGS2CovList.append(avgCov)
	
			# WORK OUT THE SCORE DIFFERENCE FROM MEDIAN
			for each in NGS2intervalList:
				scoreDiff = (each - NGS2median)
				NGS2intDiff.append(scoreDiff)

			for each in NGS3intervalList:
				scoreDiff = (each - NGS3median)
				NGS3intDiff.append(scoreDiff)

			# FOR EACH SAMPLE IN NGS2 SAMPLE GROUP
			for sample in NGS2:	
				position = NGS2.index(sample)
				sampleNGS = 'NGS000'+sample
				samMed = float(varDict[sample][0])
				samDiff = float(varDict[sample][1])
	
				# CALCULATE VARIANCE IN SCORE DIFFERENCE FROM MEDIAN FOR THE INTERVAL + SAMPLE GROUP
				MADVar = NGSStats.MAD(NGS2intDiff, NGSStats.median(NGS2intDiff))
				NGS2VarList.append(MADVar)
		
				# USE THE SLOPES CALCULATED FROM SPIKED DATA PREVIOUSLY
				slope = -0.01208 
				trendline = slope * avgCov
				MADslope = -0.009202
				MADline = MADslope * avgCov
	
				# IF SAMPLE FALLS WITHIN THE SET PARAMETERS
				intDiff = NGS2intDiff[position]
				if (intDiff <= trendline + MADVar) or (intDiff <= MADline):	
					
					# IDENTIFY GENE OR GENES WHICH THE INTERVAL COVERS 
					startGenes = intervalRefDict[chromosome][str(start)]
					endGenes = intervalRefDict[chromosome][str(end)]
			
					geneList = startGenes[:]
					for j in endGenes:
						if j not in geneList:
							geneList.append(j)
			
					strGenes = ','.join(geneList)
	
					# USE A QUERY TO MITOMINER TO EXTRACT SOME INFORMATION ABOUT THE RELEVANCE OF THE
					# PROTEIN(S) TO MITOCHONDRIA
					query = service.new_query("Gene")
			
					# SPECIFY THE OUTPUT COLUMNS
					query.add_view("proteins.primaryAccession", "proteins.uniprotName", "proteins.name", "proteins.organism.shortName", "proteins.mitoEvidenceGFP", "proteins.mitoEvidenceMassSpec",  "proteins.mitoEvidenceGO", "proteins.mitoEvidenceHumanProteinAtlas", "proteins.mitoEvidenceIMPI", "proteins.mitoEvidenceMitoCarta", "proteins.mitoEvidenceMitoMinerRefSet", "proteins.MTSipsort",  "proteins.MTSmitoprot", "proteins.MTStargetP", "Gene.organism.name")

					# ADD CONSTRAINTS TO THE QUERY - GENES TO LOOK FOR AND ONLY IN HUMAN		
					query.add_constraint("proteins.organism.shortName", "=", "H. sapiens")
					query.add_constraint("Gene.symbol", "ONE OF", geneList)

					# CREATE A RANK OF HOW 'MITOCHONDRIAL' THE MOST MITOCHONDRIAL PROTEIN OF THE SELECTION IS
					mitoRank = 0
					mitoProtein = ""
					proteinDesc = ''

					for row in query.rows():
			
						newMitoRank = 0	
						newMitoProtein = ''
		
						#GFPEvi = row["proteins.mitoEvidenceGFP"]
						#MassSpecEvi = row["proteins.mitoEvidenceMassSpec"]
	
						GOEvi = row["proteins.mitoEvidenceGO"]
						if GOEvi == True:
							newMitoRank +=1

						HPAEvi = row["proteins.mitoEvidenceHumanProteinAtlas"]
						if HPAEvi == True:
							newMitoRank +=1
	
						IMPI = row["proteins.mitoEvidenceIMPI"]
						if IMPI == True:
							newMitoRank +=1

						MitoCarta = row["proteins.mitoEvidenceMitoCarta"]
						if MitoCarta == True:
							newMitoRank +=1

						MitoMinerRef = row["proteins.mitoEvidenceMitoMinerRefSet"]
						if MitoMinerRef == True:
							newMitoRank +=1
	
						ipsort = row["proteins.MTSipsort"]
						if ipsort == 1:
							newMitoRank +=1

						mitoprot = row["proteins.MTSmitoprot"]
						if mitoprot > 0.8:
							newMitoRank +=1
			
						targetP = row["proteins.MTStargetP"]
						if targetP > 0.8:
							newMitoRank +=1

						if newMitoRank > mitoRank:
							mitoRank = int(str(newMitoRank)[:])
							mitoProtein = row['proteins.primaryAccession']
							proteinDesc = row['proteins.name']

						elif mitoProtein == '':
							mitoProtein = row['proteins.primaryAccession']
							proteinDesc = row['proteins.name']
	
					# WRITE OUTPUT TO FILE
					chrm = chromosome[3:]			
					NGS2Output.write(str(mitoRank) +'\t'+ mitoProtein +'\t'+ proteinDesc +'\t'+ strGenes +'\t'+ chrm +'\t'+ str(start) +'\t'+ str(end) +'\t')
					roundDiff = round(intDiff,2)
					NGS2Output.write(str(avgCov) +'\t'+ sampleNGS + '\t'+ str(NGS2intervalList[position]) +'\t'+ str(roundDiff)+'\t')

					for entry in NGS2covList:
						NGS2Output.write(str(entry) + '\t')

					for entry in NGS3covList:
						NGS2Output.write(str(entry) + '\t')

					for entry in NGS2intDiff:
						NGS2Output.write(str(entry) + '\t')
	
					for entry in NGS3intDiff:
						NGS2Output.write(str(entry) + '\t')

					NGS2Output.write('\n')

					print 'Written '+strGenes+' to file'
		
		# ONLY IF AVERAGE COVERAGE FOR THE INTERVAL FOR THE SAMPLE GROUP IS OVER 10
		if NGSStats.mean(NGS3covList) > 10:
			NGS2intDiff = []
			NGS3intDiff = []

			# CALCULATE AVERAGE COVERAGE FOR THE INTERVAL, WITHIN THE SAMPLE GROUP
			avgCov = NGSStats.mean(NGS3covList)
			NGS3CovList.append(avgCov)

			# WORK OUT THE SCORE DIFFERENCE FROM THE MEDIAN
			for each in NGS2intervalList:
				scoreDiff = (each - NGS2median)
				NGS2intDiff.append(scoreDiff)

			for each in NGS3intervalList:
				scoreDiff = (each - NGS3median)
				NGS3intDiff.append(scoreDiff)

			# FOR EACH SAMPLE IN NGS3 SAMPLE GROUP
			for sample in NGS3:	
				position = NGS3.index(sample)
				sampleNGS = 'NGS000'+sample
				samMed = float(varDict[sample][0])
				samDiff = float(varDict[sample][1])
	
				# CALCULATE VARIANCE IN SCORE DIFFERENCE FROM MEDIAN FOR THE INTERVAL + SAMPLE GROUP
				MADVar = NGSStats.MAD(NGS3intDiff, NGSStats.median(NGS3intDiff))
				NGS3VarList.append(MADVar)

				# USE THE SLOPES CALCULATED FROM SPIKED DATA PREVIOUSLY
				slope = -0.01515
				trendline = slope * avgCov
				MADslope = -0.01094
				MADline = MADslope * avgCov
	
				# IF SAMPLE FALLS WITHIN THE SET PARAMETERS, ADD TO THE LIST
				intDiff = NGS3intDiff[position]
				if (intDiff <= trendline + MADVar) or (intDiff <= MADline):
	
					# IDENTIFY GENE OR GENES WHICH THE INTERVAL COVERS 
					startGenes = intervalRefDict[chromosome][str(start)]
					endGenes = intervalRefDict[chromosome][str(end)]
			
					geneList = startGenes[:]
					for j in endGenes:
						if j not in geneList:
							geneList.append(j)
			
					strGenes = ','.join(geneList)
	
					# USE A QUERY TO MITOMINER TO EXTRACT SOME INFORMATION ABOUT THE RELEVANCE OF THE
					# PROTEIN(S) TO MITOCHONDRIA
					query = service.new_query("Gene")
			
					# SPECIFY THE OUTPUT COLUMNS
					query.add_view("proteins.primaryAccession", "proteins.uniprotName", "proteins.name", "proteins.organism.shortName", "proteins.mitoEvidenceGFP", "proteins.mitoEvidenceMassSpec",  "proteins.mitoEvidenceGO", "proteins.mitoEvidenceHumanProteinAtlas", "proteins.mitoEvidenceIMPI", "proteins.mitoEvidenceMitoCarta", "proteins.mitoEvidenceMitoMinerRefSet", "proteins.MTSipsort",  "proteins.MTSmitoprot", "proteins.MTStargetP", "Gene.organism.name")

					# ADD CONSTRAINTS TO THE QUERY - GENES TO LOOK FOR AND ONLY IN HUMAN		
					query.add_constraint("proteins.organism.shortName", "=", "H. sapiens")
					query.add_constraint("Gene.symbol", "ONE OF", geneList)

					# CREATE A RANK OF HOW 'MITOCHONDRIAL' THE MOST MITOCHONDRIAL PROTEIN OF THE SELECTION IS
					mitoRank = 0
					mitoProtein = ""
					proteinDesc = ''

					for row in query.rows():
			
						newMitoRank = 0	
						newMitoProtein = ''
		
						#GFPEvi = row["proteins.mitoEvidenceGFP"]
						#MassSpecEvi = row["proteins.mitoEvidenceMassSpec"]
	
						GOEvi = row["proteins.mitoEvidenceGO"]
						if GOEvi == True:
							newMitoRank +=1

						HPAEvi = row["proteins.mitoEvidenceHumanProteinAtlas"]
						if HPAEvi == True:
							newMitoRank +=1
	
						IMPI = row["proteins.mitoEvidenceIMPI"]
						if IMPI == True:
							newMitoRank +=1

						MitoCarta = row["proteins.mitoEvidenceMitoCarta"]
						if MitoCarta == True:
							newMitoRank +=1

						MitoMinerRef = row["proteins.mitoEvidenceMitoMinerRefSet"]
						if MitoMinerRef == True:
							newMitoRank +=1
	
						ipsort = row["proteins.MTSipsort"]
						if ipsort == 1:
							newMitoRank +=1

						mitoprot = row["proteins.MTSmitoprot"]
						if mitoprot > 0.8:
							newMitoRank +=1
			
						targetP = row["proteins.MTStargetP"]
						if targetP > 0.8:
							newMitoRank +=1

						if newMitoRank > mitoRank:
							mitoRank = int(str(newMitoRank)[:])
							mitoProtein = row['proteins.primaryAccession']
							proteinDesc = row['proteins.name']

						elif mitoProtein == '':
							mitoProtein = row['proteins.primaryAccession']
							proteinDesc = row['proteins.name']
	
					# WRITE OUTPUT TO FILE	
					chrm = chromosome[3:]		
					NGS3Output.write(str(mitoRank) +'\t'+ mitoProtein +'\t'+ proteinDesc +'\t'+ strGenes +'\t'+ chrm +'\t'+ str(start) +'\t'+ str(end) +'\t')
					roundDiff = round(intDiff,2)
					NGS3Output.write(str(avgCov) +'\t'+ sampleNGS + '\t'+ str(NGS3intervalList[position]) +'\t'+ str(roundDiff)+'\t')

					for entry in NGS2covList:
						NGS3Output.write(str(entry) + '\t')

					for entry in NGS3covList:
						NGS3Output.write(str(entry) + '\t')

					for entry in NGS2intDiff:
						NGS3Output.write(str(entry) + '\t')
	
					for entry in NGS3intDiff:
						NGS3Output.write(str(entry) + '\t')

					NGS3Output.write('\n')

					print 'Written '+strGenes+' to file'

NGS2Output.close()
NGS3Output.close()
