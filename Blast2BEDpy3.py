'''
Blastn output -> BED parser. 
record bait hits against contigs

REQUIRES PYTHON3

HOW TO USE:
[$user]: python3 blastnC.py [path/to/input.txt] [path/to/output.bed] ([path/to.CoverageAverage.txt])
() = optional
blastnC.py needs to be in the current working directory to be used in this form, otherwise just ensure the full path to
the blastnC.py script is entered

jwarnold1995@outlook.com
'''
#---------------------IMPORT MODULES---------------------------
import sys
import os
import pprint
import statistics


class NewRead: #The main class
	#--------------DEFINE GLOBAL VARIABLES-------------------
	global rArray, cArray
	rArray = {}
	cArray = {}
	global red, blue, lightBlue, green, yellow
	red = "255,0,0" #WEAK HITS <80% COVERAGE
	blue = "0,0,255" #80% COVERAGE
	lightBlue = "0,255,239" #90% COVERAGE
	green = "0,255,0" #95% COVERAGE
	yellow = "255,255,0" #98% COVERAGE


	def __init__(self,newFile):
		if len(sys.argv) <2:
			print('Error: Incorrect Arguments Given: python3 blastnC.py [input.txt] [output.bed]')
		else:
			self.cFile = newFile
			self.sRead(self.cFile)
		
	def sRead(self, iFile):#Parse the file into variables which can be assessed and generate colour coded results
		self.nRead = iFile
		for line in self.nRead:
			(self.baitID,
			 self.contigID,
			 self.percIdentity,
			  self.alnLength,
			   self.mismatchCount,
			    self.gapOpenCount,
			     self.queryStart,
			      self.queryEnd,
			       self.subjectStart,
			        self.subjectEnd,
			         self.eVal,
			          self.bitScore) = line.split('\t')#Assign the blastn output to corresponding variables
			'''
			if float(self.percIdentity) >= 80.00 and int(self.alnLength) >= 100:
				self.addToArray(self.baitID, self.contigID, self.subjectStart, self.subjectEnd)
			else:
				print('Error lost')
			'''
			self.coverageAverage(self.contigID, self.percIdentity)
#--------------------------------------ASSESSING STRENGTH OF HIT--------------------------------------
			if float(self.percIdentity) >= 98.00 and int(self.alnLength) >= 100:
				self.addToArray(self.contigID, self.subjectStart, self.subjectEnd, yellow)
			elif float(self.percIdentity) >= 95.00 and int(self.alnLength) >= 100:
				self.addToArray(self.contigID, self.subjectStart, self.subjectEnd, green)
			elif float(self.percIdentity) >= 90.00 and int(self.alnLength) >= 100:
				self.addToArray(self.contigID, self.subjectStart, self.subjectEnd, lightBlue)
			elif float(self.percIdentity) >= 80.00 and int(self.alnLength) >= 100:
				self.addToArray(self.contigID, self.subjectStart, self.subjectEnd, blue)
			else:#WEAK HITS
				self.addToArray(self.contigID, self.subjectStart, self.subjectEnd, red)

		self.createBed(rArray) #Create the .BED file from the completed Dictionary
		#self.printResult(rArray)
		self.printCovAvg()

	def addToArray(self, iContig, iStart, iEnd, colour):#This assigns the input from sRead into the global array 
		self._contig = iContig
		self._start = iStart
		self._end = iEnd
		self._color = colour


		if self._contig not in rArray:#Initializing the contig header
			rArray.setdefault(self._contig, [self._start,self._end,self._color])
		else:
			rArray[self._contig] +=  [self._start,self._end,self._color]

		
	'''
	#A testing function to evaluate the success of the conversion into array
	def printResult(self, nArray):
		self.neArray = nArray	
		pprint.pprint(self.neArray)	
	'''

	def createBed(self, bArray):#Parses the global array 'rArray' into a BED file
		self.finalArray = bArray
		self.nString = ''
		self.nName = sys.argv[2]#Retrieve title of desired output file name from argument line
		self.bedFILE = open(self.nName, 'w+')#create and open the new file
		self.trackString = "track name=\""+self.nName+"\" itemRgb=\"On\"\n"
		self.bedFILE.write(self.trackString)
		for key in self.finalArray.keys():
			#self.leng = len(self.finalArray[key])
			#print(self.leng)
			for i in range(0, len(self.finalArray[key]),3):#allows for lists within the dictionary to be fed out in 3s
				
				self.nString = key + '\t' + self.finalArray[key][i] + '\t' + self.finalArray[key][i+1] + '\t' + self.finalArray[key][i+2] +'\n'
				self.bedFILE.write(self.nString)
				self.nString = ''
		self.bedFILE.close()
				
	def coverageAverage(self, con, cov):
		self._nContig = con
		self._nCov = cov

		if self._nContig not in cArray:
			cArray.setdefault(self._nContig, [float(self._nCov),])
		else:
			cArray[self._nContig] += [float(self._nCov),]

	def printCovAvg(self):
		self._avgCov = 0.0
		try:

			self.argv3 = sys.argv[3]
			self.anotherOne = open(self.argv3, 'w+')
			for k in cArray.keys():
				self.lenAr = len(cArray[k])
				self.lenArC = float(self.lenAr)
				self.sumAr = cArray[k]

				self.s = statistics.mean((self.sumAr))
				self.anotherOne.write(k + '\t' + str(self.s) + '\n') 
		except ValueError:
			print('?')
			pass
		
			

nRequest = sys.argv[1]
nFile = open(nRequest, 'r')

if __name__=='__main__':
	run = NewRead(nFile)

