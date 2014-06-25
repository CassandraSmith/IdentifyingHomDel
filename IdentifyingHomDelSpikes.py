#!/usr/mbu/python/x86_64/bin/python
# SPIKING A RANDOM SAMPLE OF INTERVALS TO LOOK LIKE A HOMOZYGOUS EXON DELETION
# THEN ATTEMPTING TO IDENTIFY THESE USING A METHOD BASED ON THE VARIANCE CALCULATED
# OUTPUTS A PRINTED SUMMARY OF HOW MANY OF THE SPIKED INTERVALS WERE IDENTIFIED IN NGS2 AND NGS3
# AND HOW MANY OTHER INTERVALS WERE ALSO IDENTIFIED

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

# ---------------------------------------------------------------------------------------------------------------
# DEFINE THE SAMPLE GROUPS
NGS2 = ['04','05','06','07','08','09']
NGS3 = ['11','12','13','14','15']

# ----------------------------------------------------------------------------------------------------------------
# CREATE OUTPUT FILES
NGS2Output = open('/data/bioinf/home/cls87/SpikingData/HomDel/NGS2MissedSpikes.txt','w')
NGS3Output = open('/data/bioinf/home/cls87/SpikingData/HomDel/NGS3MissedSpikes.txt','w')

# ----------------------------------------------------------------------------------------------------------------
# CHOOSE AT RANDOM A PART, INTERVAL ADN SAMPLE TO BE 'SPIKED' FOR HOMOZYGOUS DELETION
# CHOOSE A PART; THIS IS REPEATED FOR THE NUMBER OF TIMES IN RANGE()

randParts = []
for i in range(1,20001):
	randPart = random.randrange(1,17,1)
	randParts.append(randPart)

# SET UP SOME LISTS AND COUNTERS FOR STORING INFORMATION
NGS2CovList = []
NGS3CovList = []

NGS2SpikedList = []
NGS3SpikedList = []

NGS2VarList = []
NGS3VarList = []

identified = 0
NGS2identified = 0
NGS3identified = 0
missed = []
NGS2missed = []
NGS3missed = []

scoreident = 0
NGS2scoreident = 0
NGS3scoreident = 0
totaleach = 0
NGS2total = 0
NGS3total = 0

# FOR EACH PART THAT WAS RANDOMLY CHOSEN
for partNo in randParts:
	# OPEN INTERVAL SUMMARY DATA FOR THAT PART
	intervalSummary = open('/data/bioinf/home/cls87/DOCAnalysis/pt'+str(partNo)+'/pt'+str(partNo)+'DOCAnalysis.sample_interval_summary').readlines()
	intervalSummary = intervalSummary[1:]

	# RETRIEVE THE NUMBER OF INTERVALS IN THE PART AND RANDOMLY CHOOSE ONE OF THEM
	lengthPart = len(intervalSummary)
	randInt = random.randrange(0, lengthPart, 1)

	# EXTRACTING INFORMATION ABOUT THE CHOSEN INTERVAL
	d = intervalSummary[randInt]
	d = d[:-1]
	sd = d.split('\t')
	interval = sd[0]
	sInterval = interval.split(':')
	chromosome = sInterval[0]
	startEnd = sInterval[1]
	sStartEnd = startEnd.split('-')
	start = sStartEnd[0]
	end = sStartEnd[1]

	# SET UP SOME MORE LISTS FOR STORING DATA, AT AN INTERVAL LEVEL THIS TIME
	NGS2intervalList = []
	NGS2intervalOrder = []
	NGS2covList = []
	NGS2intDiff = []

	NGS3intervalList = []
	NGS3intervalOrder = []
	NGS3covList = []
	NGS3intDiff = []

	# CHOOSE GROUP AND SAMPLE TO 'SPIKE'
	# THE PRINT STATEMENTS ALLOW THE USER TO KEEP TRACK OF THE PROGRAMS PROGRESS
	randGroup = random.randrange(2,4,1)
	if randGroup == 2:
		randSample = random.choice(NGS2)
		print partNo, randInt, randSample
	elif randGroup == 3:
		randSample = random.choice(NGS3)
		print partNo, randInt, randSample

	# FOR EACH SAMPLE IN THE INTERVAL, IN NUMERICAL ORDER
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
		# AND SPIKING THE DATA TO ZERO FOR THE IDENTIFIED SAMPLE
		# (THIS WILL BE MODIFIED LATER ON TO CREATE SOME NOISE, 
		# BUT THIS ALLOWS AN APROX INTERVAL AVERAGE COVERAGE TO BE CALCULATED)
		if randSample == strWhichSam:
			intSamMeanCov = 0
		else:
			intSamMeanCov = float(sd[int(h)])

		# THERE WAS NO MEDIAN/SD CALCULATED FOR NGS00001
		if strWhichSam != '01':
			# EXTRACTING MEDIAN AND SD FOR THE SPECIFIC SAMPLE
			samMed = float(varDict[strWhichSam][0])
			samDiff = float(varDict[strWhichSam][1])
		
			# CALCULATE THE STANDARDISED SCORE COMPARATIVE TO SAMPLE COVERAGE
			score = round((intSamMeanCov - samMed)/samDiff, 4)	

		# APPEND SCORES AND COVERAGE TO LIST
		if strWhichSam in NGS2:
			NGS2intervalList.append(score)
			NGS2covList.append(intSamMeanCov)
			NGS2intervalOrder.append(strWhichSam)

		elif strWhichSam in NGS3:
			NGS3intervalList.append(score)
			NGS3covList.append(intSamMeanCov)
			NGS3intervalOrder.append(strWhichSam)

	# CALCULATE THE MEDIAN SCORE FOR EACH SAMPLE GROUP
	NGS2median = NGSStats.median(NGS2intervalList)
	NGS3median = NGSStats.median(NGS3intervalList)

	# ONLY CONTINUE FOR INTERVALS WITH AN AVERAGE COVERAGE OVER THE SAMPLE GROUP >10	
	if randGroup == 2 and NGSStats.mean(NGS2covList) > 10:
		spikePos = NGS2intervalOrder.index(randSample)
		samMed = float(varDict[randSample][0])
		samDiff = float(varDict[randSample][1])

		# CHANGE TO RANDOM SCORE CALCULATED FROM <10% AVERAGE COVERAGE OF THE INTERVAL
		avgCov = NGSStats.mean(NGS2covList)
		change = (float(random.uniform(0,10)))/100
		NGS2covList[spikePos] = avgCov * change

		NGS2intervalList[spikePos] = round((NGS2covList[spikePos] - samMed)/samDiff, 4)	
		NGS2CovList.append(avgCov)

		# CALCULATE DIFFERENCE FROM MEDIAN SCORE FOR BOTH SPIKED AND NON-SPIKED DATA
		# ONE SAMPLE SHOULD ALWAYS BE 0, DUE TO THE USE OF THE MEDIAN RATHER THAN THE MEAN
		spikeScore = NGS2intervalList[spikePos]
		spikeDiff = (spikeScore - NGS2median)
		NGS2intDiff.append(spikeDiff)
		NGS2SpikedList.append(spikeDiff)

		del NGS2intervalList[spikePos]
		for each in NGS2intervalList:
			scoreDiff = (each - NGS2median)
			NGS2intDiff.append(scoreDiff)

		# CALCULATE VARIANCE MEASUREMENT FOR RANGES OF SPIKED DATA
		# HERE, HAVE USED MAD - MEAN AVERAGE DISTANCE FROM THE MEDIAN
		# THIS IS LESS SUSCEPTABLE TO SKEWING BY OUTLIERS
		MADVar = NGSStats.MAD(NGS2intDiff, NGSStats.median(NGS2intDiff))
		NGS2VarList.append(MADVar)

		# THESE VALUES FOR THE FIT OF THE SPIKED DATA AND SPIKED DATA PLUS ONE MAD
		# WERE CALCULATED USING R IN THE DIFFVCOVERAGE PROGRAM, 20000 SAMPLE SIZE, FORCED THROUGH ORIGIN
		slope = -0.01208
		trendline = slope * avgCov
		MADslope = -0.009202
		MADline = MADslope * avgCov

		# TWO METHODS OF IDENTIFYING SPIKES TO REDUCE FALSE NEGATIVE RATE
		# 1. IF SPIKE DIFF SCORE IS LESS THAN THE SPIKED LEAST SQUARES + ONE MEASURE OF VARIANCE FOR THAT INTERVAL
		if spikeDiff <= trendline +MADVar:
			print 'spike identified'
			identified += 1
			NGS2identified +=1
		# 2. IF SPIKE DIFF SCORE IS LESS THAN THE LEAST SQUARES FIT FOR ALL INTERVALS OF PLUS 1 MEASURE OF VARIANCE
		elif spikeDiff <= MADline:
			print 'spike identified'
			identified += 1
			NGS2identified +=1
		# IF NOT IDENTIFIED, ADD INFORMATION TO THE MISSED LIST
		else:
			missed.append(avgCov)
			NGS2missed.append([spikeDiff,avgCov])

		# SEE IF THE PROCESS IDENTIFIES ANY MORE INTERVALS AS POSSIBLE HOMOZYGOUS DELETIONS
		for each in NGS2intDiff:
			totaleach += 1
			NGS2total += 1
			if each <= trendline + MADVar:			
				scoreident +=1
				NGS2scoreident +=1
			elif each <= MADline:
				scoreident +=1
				NGS2scoreident +=1

	# ONLY CONTINUE FOR INTERVALS WITH AN AVERAGE COVERAGE OVER THE SAMPLE GROUP >10	
	if randGroup == 3 and NGSStats.mean(NGS3covList) > 10:
		spikePos = NGS3intervalOrder.index(randSample)
		samMed = float(varDict[randSample][0])
		samDiff = float(varDict[randSample][1])

		# CHANGE TO RANDOM SCORE CALCULATED FROM <10% AVERAGE COVERAGE OF THE INTERVAL
		avgCov = NGSStats.mean(NGS3covList)
		change = (float(random.uniform(0,10)))/100
		NGS3covList[spikePos] = avgCov * change

		NGS3intervalList[spikePos] = round((NGS3covList[spikePos] - samMed)/samDiff, 4)
		NGS3CovList.append(avgCov)

		# CALCULATE DIFFERENCE FROM MEDIAN SCORE FOR BOTH SPIKED AND NON-SPIKED DATA
		# ONE SAMPLE SHOULD ALWAYS BE 0, DUE TO THE USE OF THE MEDIAN RATHER THAN THE MEAN
		spikeScore = NGS3intervalList[spikePos]
		spikeDiff = (spikeScore - NGS3median)
		NGS3intDiff.append(spikeDiff)
		NGS3SpikedList.append(spikeDiff)

		del NGS3intervalList[spikePos]
		for each in NGS3intervalList:
			scoreDiff = (each - NGS3median)
			NGS3intDiff.append(scoreDiff)

		# CALCULATE VARIANCE MEASUREMENT FOR RANGES OF SPIKED DATA
		# HERE, HAVE USED MAD - MEAN AVERAGE DISTANCE FROM THE MEDIAN
		# THIS IS LESS SUSCEPTABLE TO SKEWING BY OUTLIERS
		MADVar = NGSStats.MAD(NGS3intDiff, NGSStats.median(NGS3intDiff))
		NGS3VarList.append(MADVar)

		# THESE VALUES FOR THE FIT OF THE SPIKED DATA AND SPIKED DATA PLUS ONE MAD
		# WERE CALCULATED USING R IN THE DIFFVCOVERAGE PROGRAM, 20000 SAMPLE SIZE, FORCED THROUGH ORIGIN
		slope = -0.01515
		MADslope = -0.01094
		trendline = slope * avgCov
		MADline = MADslope * avgCov

		# TWO METHODS OF IDENTIFYING SPIKES TO REDUCE FALSE NEGATIVE RATE
		# 1. IF SPIKE DIFF SCORE IS LESS THAN THE SPIKED LEAST SQUARES + ONE MEASURE OF VARIANCE FOR THAT INTERVAL
		if spikeDiff <= trendline +MADVar:
			print 'spike identified'
			identified += 1
			NGS3identified +=1
		# 2. IF SPIKE DIFF SCORE IS LESS THAN THE LEAST SQUARES FIT FOR ALL INTERVALS OF PLUS 1 MEASURE OF VARIANCE
		elif spikeDiff <= MADline:
			print 'spike identified'
			identified += 1
			NGS3identified +=1
		# IF NOT IDENTIFIED, ADD INFORMATION TO THE MISSED LIST
		else:
			missed.append(avgCov)
			NGS3missed.append([spikeDiff,avgCov])

		# SEE IF THE PROCESS IDENTIFIES ANY MORE INTERVALS AS POSSIBLE HOMOZYGOUS DELETIONS
		for each in NGS3intDiff:
			totaleach += 1
			NGS3total += 1
			if each <= trendline + MADVar:		
				scoreident += 1
				NGS3scoreident +=1
			elif each <= MADline:
				scoreident += 1
				NGS3scoreident +=1	

# PRINT OUTPUT - SOME INFORMATION ABOUT HOW MANY SPIKED DELETIONS WERE IDENTIFIED AND MISSED
# PRINT OUTPUT - HOW MANY OTHER INTERVALS WERE IDENTIFIED AS POSSIBLE DELETIONS?
print 'NGS2: Identified -', NGS2identified, 'Missed -', len(NGS2missed), 'New -', (NGS2scoreident - NGS2identified), 'Total -', NGS2total
# FILE OUTPUT - INFORMATION ABOUT MISSED SPIKES AVERAGE COVERAGE AND DIFFERENCE FROM INTERVAL, SO CAN LOOK FOR TRENDS IN WHAT IS MISSED
for entry in NGS2missed:
	NGS2Output.write(str(entry[0]) + '\t'+str(entry[1]) + '\n')
print 'NGS3: Identified -', NGS3identified, 'Missed -', len(NGS3missed), 'New -', (NGS3scoreident - NGS3identified), 'Total -', NGS3total
for entry in NGS3missed:
	NGS3Output.write(str(entry[0]) + '\t'+str(entry[1])+ '\n')
