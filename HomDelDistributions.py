#!/usr/mbu/python/x86_64/bin/python
# SPIKING DATA WITH A HOMOZYGOUS DELETION TO GET SOME IDEA OF THE DISTRIBUTION OF POSSIBLE HOMOZYGOUS DELETIONS
# OUTPUTS INFORMATION FOR EXCEL TO PRODUCE A 'HISTOGRAM'

from intermine.webservice import Service
import math, sys, random
sys.path.append('/data/bioinf/home/cls87/Modules')
import NGSStats
 ----------------------------------------------------------------------------------------------------------------
# OPEN INTERVAL SUMMARY DATA FOR HEADERS
intervalSummary = open('/data/bioinf/home/cls87/DOCAnalysis/pt1/pt1DOCAnalysis.sample_interval_summary').readlines()
headers = intervalSummary[0]
 ----------------------------------------------------------------------------------------------------------------
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
NGS2Output = open('/data/bioinf/home/cls87/SpikingData/HomDel/Distribution/NGS2Spiked.txt','w')
NGS3Output = open('/data/bioinf/home/cls87/SpikingData/HomDel/Distribution/NGS3Spiked.txt','w')
NGS2ScoreOutput = open('/data/bioinf/home/cls87/SpikingData/HomDel/Distribution/NGS2Score.txt','w')
NGS3ScoreOutput = open('/data/bioinf/home/cls87/SpikingData/HomDel/Distribution/NGS3Score.txt','w')
# ----------------------------------------------------------------------------------------------------------------
# CHOOSE A PART, INTERVAL AND SAMPLE TO SPIKE WITH A 'HOMOZYGOUS DELETION'

randParts = []
for i in range(1,20001):
	randPart = random.randrange(1,17,1)
	randParts.append(randPart)

# SET UP SOME LISTS FOR STORING INFORMATION
NGS2ScoreList = []
NGS3ScoreList = []

NGS2SpikedList = []
NGS3SpikedList = []

# FOR EACH RANDOMLY CHOSEN PART
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

	# SET UP SOME MORE LISTS FOR STORING DATA, AT AN INTERVAL LEVEL
	NGS2intervalList = []
	NGS2intervalOrder = []
	NGS2covList = []

	NGS3intervalList = []
	NGS3intervalOrder = []
	NGS3covList = []

	# CHOOSE GROUP AND SAMPLE TO 'SPIKE'
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

	# CALCULATE THE MAD VARIATION FOR EACH SAMPLE GROUP OVER THE INTERVAL
	NGS2MAD = NGSStats.MAD(NGS2intervalList, NGS2median)
	NGS3MAD = NGSStats.MAD(NGS3intervalList, NGS3median)

	# IN THIS CASE LOOK FOR AVERAGE COVERAGE OVER THE INTERVAL OF MORE THAN 5
	# THIS MEANS WE ARE NOT CONSIDERING THOSE INTERVALS ALREADY DETERMINED TO BE COMPLETELY 'MISSING' PREVIOUSLY	
	if randGroup == 2 and NGSStats.mean(NGS2covList) > 5:
		# CALCULATE THE SCORE DIFFERENCE FROM THE MEAN FOR THE SPIKED SAMPLE
		spikePos = NGS2intervalOrder.index(randSample)
		spikeScore = NGS2intervalList[spikePos]
		spikeDiff = (spikeScore - NGS2median)
		NGS2SpikedList.append(spikeDiff)

		# CALCULATE THE SCORE DIFFERENCE FROM THE MEAN FOR THE OTHER SAMPLES
		del NGS2intervalList[spikePos]
		for each in NGS2intervalList:
			scoreDiff = (each - NGS2median)
			NGS2ScoreList.append(scoreDiff)		

	elif randGroup == 3 and NGSStats.mean(NGS3covList) > 5:
		# CALCUALTE THE SCORE DIFFERENCE FROM THE MEAN FOR THE SPIKED SAMPLE
		spikePos = NGS3intervalOrder.index(randSample)
		spikeScore = NGS3intervalList[spikePos]
		spikeDiff = (spikeScore - NGS3median)
		NGS3SpikedList.append(spikeDiff)

		# CALCULATE THE SCORE DIFFERENCE FROM THE MEAN FOR THE OTHER SAMPLES
		del NGS3intervalList[spikePos]
		for each in NGS3intervalList:
			scoreDiff = (each - NGS3median)
			NGS3ScoreList.append(scoreDiff)

# ------------------------------------------------------------------------------------------------
# OUTPUTTING DATA FOR EXCEL 'HISTOGRAM', WITH ALL BINS THE SAME SIZE
# ONE EACH OUTPUTTED FOR SPIKED AND NON-SPIKED DATA, FOR EACH SAMPLE GROUP

# -------------------------------------------------------------------------------------------------
# DISTRIBUTION FOR SPIKED NGS2
maximum = max(NGS2SpikedList)
minimum = min(NGS2SpikedList)

sortCovList = sorted(NGS2SpikedList)

binFloor = round(math.floor(minimum), 1)
binSize = 0.01

while binFloor <= math.ceil(maximum):
	frequencyCount = 0
	binCeil = binFloor + binSize

	while len(sortCovList) > 0 and binFloor < sortCovList[0] <= binCeil:
		frequencyCount +=1
		del sortCovList[0]
	else:
		frequencyDensity = frequencyCount/binSize
		print binCeil, frequencyDensity
		NGS2Output.write(str(round(binCeil, 3))+'\t'+str(frequencyDensity)+'\n')
		binFloor += binSize

NGS2Output.close()

# ------------------------------------------------------------------------------------------------
# DISTRIBUTION FOR SPIKED NGS3
maximum = max(NGS3SpikedList)
minimum = min(NGS3SpikedList)

sortCovList = sorted(NGS3SpikedList)

binFloor = round(math.floor(minimum), 1)
binSize = 0.01

while binFloor <= math.ceil(maximum):
	frequencyCount = 0
	binCeil = binFloor + binSize

	while len(sortCovList) > 0 and binFloor < sortCovList[0] <= binCeil:
		frequencyCount +=1
		del sortCovList[0]
	else:
		frequencyDensity = frequencyCount/binSize
		print binCeil, frequencyDensity
		NGS3Output.write(str(round(binCeil, 3))+'\t'+str(frequencyDensity)+'\n')
		binFloor += binSize

NGS3Output.close()

# ------------------------------------------------------------------------------------------------
# DISTRIBUTION FOR NON-SPIKED NGS2
maximum = max(NGS2ScoreList)
minimum = min(NGS2ScoreList)

sortCovList = sorted(NGS2ScoreList)

binFloor = round(math.floor(minimum), 1)
binSize = 0.01

while binFloor <= math.ceil(maximum):
	frequencyCount = 0
	binCeil = binFloor + binSize

	while len(sortCovList) > 0 and binFloor < sortCovList[0] <= binCeil:
		frequencyCount +=1
		del sortCovList[0]
	else:
		frequencyDensity = frequencyCount/binSize
		print binCeil, frequencyDensity
		NGS2ScoreOutput.write(str(round(binCeil, 3))+'\t'+str(frequencyDensity)+'\n')
		binFloor += binSize

NGS2ScoreOutput.close()

# ------------------------------------------------------------------------------------------------
# DISTRIBUTION FOR NON-SPIKED NGS3
maximum = max(NGS3ScoreList)
minimum = min(NGS3ScoreList)

sortCovList = sorted(NGS3ScoreList)

binFloor = round(math.floor(minimum), 1)
binSize = 0.01

while binFloor <= math.ceil(maximum):
	frequencyCount = 0
	binCeil = binFloor + binSize

	while len(sortCovList) > 0 and binFloor < sortCovList[0] <= binCeil:
		frequencyCount +=1
		del sortCovList[0]
	else:
		frequencyDensity = frequencyCount/binSize
		print binCeil, frequencyDensity
		NGS3ScoreOutput.write(str(round(binCeil, 3))+'\t'+str(frequencyDensity)+'\n')
		binFloor += binSize

NGS3ScoreOutput.close()

# THESE CHARACTERISTICS GIVE A VERY BRIEF IDEA OF THE OVERLAP BETWEEN SPIKED AND NON-SPIKED DATA WITHOUT BUILDING THE HISTOGRAM
print min(NGS2ScoreList), max(NGS2SpikedList)
print min(NGS3ScoreList), max(NGS3SpikedList)		
