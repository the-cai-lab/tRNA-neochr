import os
import Bio.SeqIO
import re
import subprocess
import json

strandDict = {-1:"-",0:".",1:"+"}
anticodonList = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"] 
codonDict = {"AAA":"Lys", "AAC":"Asn", "AAG":"Lys", "AAT":"Asn", "ACA":"Thr", "ACC":"Thr", "ACG":"Thr", "ACT":"Thr", "AGA":"Arg", "AGC":"Ser", "AGG":"Arg", "AGT":"Ser", "ATA":"Ile", "ATC":"Ile", "ATG":"Met", "ATT":"Ile", "CAA":"Gln", "CAC":"His", "CAG":"Gln", "CAT":"His", "CCA":"Pro", "CCC":"Pro", "CCG":"Pro", "CCT":"Pro", "CGA":"Arg", "CGC":"Arg", "CGG":"Arg", "CGT":"Arg", "CTA":"Leu", "CTC":"Leu", "CTG":"Leu", "CTT":"Leu", "GAA":"Glu", "GAC":"Asp", "GAG":"Glu", "GAT":"Asp", "GCA":"Ala", "GCC":"Ala", "GCG":"Ala", "GCT":"Ala", "GGA":"Gly", "GGC":"Gly", "GGG":"Gly", "GGT":"Gly", "GTA":"Val", "GTC":"Val", "GTG":"Val", "GTT":"Val", "TAA":"*", "TAC":"Tyr", "TAG":"*", "TAT":"Tyr", "TCA":"Ser", "TCC":"Ser", "TCG":"Ser", "TCT":"Ser", "TGA":"*", "TGC":"Cys", "TGG":"Trp", "TGT":"Cys", "TTA":"Leu", "TTC":"Phe", "TTG":"Leu", "TTT":"Phe",} 

anticodonDict = {} 




def extractAnticodon(specie,feature):
	if specie == "SC":
		return re.findall(r'\((.+)\)',feature.qualifiers['locus_tag'][0])[0]
	elif specie == "AG":
		return feature.qualifiers['codon_recognized'][0]
	elif specie == "KL":
		return re.findall(r'\((.+)\)',feature.qualifiers['note'][0])[0]

#this method is used for instead extractAnticodon, which not read anticodon from the file, it use tRNAscan to scan it
message = ""
def scanAnticodon(seq):
	f = open("test.fa","w")
	f.write(">test\n")
	f.write(seq+"\n")
	f.close()
	process = subprocess.Popen('tRNAscan-SE test.fa -q'.split(' '), stdout=subprocess.PIPE)
	try:
		out = process.communicate()[0]
		line = out.split('\n')[3].split('\t')
		anticodon = line[5]
		aa = line[4]
		return (anticodon,aa)
	except:
		return ("NNN","unknow")


def findfeature(startA,endA,start,end,strand,features,fType='gene',includeSameStrand=False):
	featureNames = []
	for f in features:
		if f.type == fType:
			startC = f.location.start
			endC = f.location.end
			if includeSameStrand:
				if ((endC > startA and startC < startA) or (startC < endA and startC>startA)) or ( startC>startA and endC < endA):
					featureNames.append({"name":f.qualifiers['locus_tag'][0],"start":startC,"end":endC,"strand":f.location.strand})
			else:
				if ((endC > startA and startC < startA) or (startC < endA and startC>startA)) and (strand != f.strand) and not(startC==start and endC==end):
					featureNames.append({"name":f.qualifiers['locus_tag'][0],"start":startC,"end":endC,"strand":f.location.strand})
	return featureNames

def find5(featureLoc,features,fType = 'gene',includeSameStrand=False):
	start = featureLoc.start	
	end = featureLoc.end
	strand = featureLoc.strand
	return findfeature(start,end,0,0,strand,features,fType,includeSameStrand)

def find3(feature,features,tail):
	start = feature.location.start	
	end = feature.location.end
	strand = feature.location.strand
	if strand > 0:
		startA = end
		endA = end+tail 
	else:
		startA = start-tail 
		endA = start 
	return findfeature(startA,endA,start,end,strand,features,True)

	

def findSolo(promoter,features,record):
	start = promoter.start
	end = promoter.end
	strand = promoter.strand

	length = end-start
	for f in features:
		if f.type == 'misc_feature' and len(re.findall('solo',f.qualifiers['note'][0]))>0:
			startC = f.location.start
			endC = f.location.end
			if (startC>start and startC<end) or (endC>start and endC<end):
				print "match %r "%(f)
				soloLength = endC-startC
				if strand>0:
					startN = start-soloLength
					endN = end
					seq = str(Bio.SeqFeature.FeatureLocation(startN,startC,strand).extract(record).seq) + str(Bio.SeqFeature.FeatureLocation(endC,endN,strand).extract(record).seq)
				else:
					startN = start
					endN = end + soloLength
					seq = str(Bio.SeqFeature.FeatureLocation(endC,endN,strand).extract(record).seq) + str(Bio.SeqFeature.FeatureLocation(startN,startC,strand).extract(record).seq)
				return seq 
	
	return ""

def extendIfSolo(promoter,features):
	start = promoter.start
	end = promoter.end
	strand = promoter.strand
	length = end-start
	solo = []
	for f in features:
		if f.type == 'misc_feature' and len(re.findall('solo',f.qualifiers['note'][0]))>0:
			startC = f.location.start
			endC = f.location.end
			if (startC>start and startC<end) or (endC>start and endC<end):
				#print "match %r "%(f)
				soloLength = endC-startC
				if strand>0:
					start = start-soloLength
				else:
					end = end + soloLength
				solo.append(Bio.SeqFeature.FeatureLocation(f.location.start,f.location.end,strand))
				print "ext %d + %d = %d"%(length,soloLength,end-start)
				print "%d->%d %d->%d"%(promoter.start,start,promoter.end,end)
				break
	
	return (Bio.SeqFeature.FeatureLocation(start,end,strand),solo)

def removeSolo(seq5,record,solo):
	seq = seq5
	if solo == []:
		return seq
	#print "before %s"%seq ,
	for s in solo:
		soloSeq = str(s.extract(record).seq)
		seq = seq.replace(soloSeq,"")
	#print "after %s"%seq
	print "%d - %d = %d"%(len(seq5),len(soloSeq),len(seq))
	return seq

		
				

def showDirtyType(dirtyType):
	s = ""
	if dirtyType & 1 >0:
		s+= "5 "
	if dirtyType & 2 >0:
		s+= "3 "
	if dirtyType & 4 >0:
		s+= "T "
	if dirtyType & 8 >0:
		s+= "S "
	if dirtyType & 16 >0:
		s+= "N "
	return s

			

def mergeStrArray(strArray):
	if strArray == []:
		return ""
	re = strArray[0]
	for i in range(1,len(strArray)):
		re+=" & "
		re+=strArray[i]
	return re

class Neochromosome:

	def __init__(self):
		self.species=["SC","AG","KL","EC"]
		self.tRnaData = {}
		self.anticodon = {}
		self.tRnaData["SC"] = {"gbPath":"gb/Saccharomyces cerevisiae"};
		self.tRnaData["AG"] = {"gbPath":"gb/Ashbya gossypii"};
		self.tRnaData["KL"] = {"gbPath":"gb/Kluyveromyces lactis"};
		self.tRnaData["EC"] = {"gbPath":"gb/Eremothecium Cymbalariae"};
		self.soloSet = {"AG":[],"SC":[],"KL":[],"EC":[]}
		self.matched = {} 


	def importDataFromJson(self):
		if os.path.isfile("anticodon.json"):
			with open("anticodon.json") as f:
				jsonStr = f.read()
				self.anticodon = json.loads(jsonStr)
		
		if os.path.isfile("neochromosome.json"):
			with open("neochromosome.json") as f:
				jsonStr = f.read()
				self.tRnaData = json.loads(jsonStr)
			return True 
		else:
			return False 
	def saveData(self):
		with open("anticodon.json","w") as f:
			s = json.dumps(self.anticodon)
			f.write(s)
		with open("soloSet.json","w") as f:
			s = json.dumps(self.soloSet)
			f.write(s)
		with open("matchedSet.json","w") as f:
			s = json.dumps(self.matched)
			f.write(s)
		with open("neochromosome.json","w") as f:
			s = json.dumps(self.tRnaData)
			f.write(s)
	def importData(self):
		if not self.importDataFromJson():
			self.importDataFromGB()
	def importDataFromGB(self):
		for specie in self.tRnaData.keys():
			fr = []
			for p,d,n in os.walk(self.tRnaData[specie]["gbPath"]):
				for fname in n:
					mainName, ext =  os.path.splitext(fname)
					lastFileName = fname
					if ext == ".gb":
						if specie =="SC":
							chromosomeName = mainName[3:]
						else:
							chromosomeName = mainName.split(" ")[-1]
						print os.path.join(p,fname)
						f = open(os.path.join(p,fname))
						count = 0
						for record in Bio.SeqIO.parse(f,"genbank"):
							for feature in record.features:	
								if feature.type == 'tRNA' and feature.qualifiers.has_key('pseudo')==False:
									lastFeature = feature
									count+=1
									fd={}
									fd["matched"] = {"AG":"","KL":"","EC":""}
									fd["id"] = "%s.t%02d.%02d"%(specie,int(chromosomeName),count)
									fd["chromosome"] = int(chromosomeName)
									fd["count"] = count
									fd["tRnaName"] = feature.qualifiers['locus_tag'][0]
									fd["original_sequence"] = str(Bio.SeqFeature.FeatureLocation(feature.location.start,feature.location.end,feature.strand).extract(record).seq)
									fd["sequence"] = str(feature.extract(record).seq)
									fd["orientation"] = strandDict[feature.strand]
									#fd["anticodon"] = extractAnticodon(specie,feature)
									if not self.anticodon.has_key(fd["tRnaName"]):
										self.anticodon[fd["tRnaName"]] = scanAnticodon(fd["sequence"])
									fd["anticodon"],fd["aa"] = self.anticodon[fd["tRnaName"]]

									fd["start"] = feature.location.start
									fd["end"] = feature.location.end
									fd["strand"] = feature.location.strand
									fd["position"] = "%d-%d"%(feature.location.start+1, feature.location.end)
									if len(feature.location.parts) >= 2:
										fd["intron"] = "%d-%d"%(feature.location.parts[0].end+1, feature.location.parts[1].start)
										fd["intronSeq"] = str(Bio.SeqFeature.FeatureLocation(feature.location.parts[0].end,feature.location.parts[1].start,feature.strand).extract(record).seq)
									else:
										fd["intron"] = ""
										fd["intronSeq"] = ""
									if feature.location.strand>0:
										seq5Loc = Bio.SeqFeature.FeatureLocation(feature.location.start-500,feature.location.start,feature.location.strand)
										fd["seq5"] = str(seq5Loc.extract(record).seq)
										fd["seq5P"] = fd["seq5"]
										seq3Loc = Bio.SeqFeature.FeatureLocation(feature.location.end,feature.location.end+20,feature.location.strand)
										fd["seq3"] = str(seq3Loc.extract(record).seq)
										fd["seq3_40"] = str(Bio.SeqFeature.FeatureLocation(feature.location.end,feature.location.end+40,feature.location.strand).extract(record).seq)
									else:
										seq5Loc = Bio.SeqFeature.FeatureLocation(feature.location.end,feature.location.end+500,feature.location.strand)
										fd["seq5"] = str(seq5Loc.extract(record).seq)
										fd["seq5P"] = fd["seq5"]
										fd["seq3"] = str(Bio.SeqFeature.FeatureLocation(feature.location.start-20,feature.location.start,feature.location.strand).extract(record).seq)
										fd["seq3_40"] = str(Bio.SeqFeature.FeatureLocation(feature.location.start-40,feature.location.start,feature.location.strand).extract(record).seq)
									
									#check solo at here, extend if nessary
									seq5LocExt,soloLocs = extendIfSolo(seq5Loc,record.features)

									#featuresIn5 = find5(feature.location, record.features)
									featuresIn5 = find5(seq5LocExt, record.features)
									fd["overlapped5"] = mergeStrArray([x['name'] for x in featuresIn5])
									featuresIn3 = find3(feature, record.features,20)
									fd["overlapped3"] = mergeStrArray([x['name'] for x in featuresIn3])
									featuresIn3_40 = find3(feature, record.features,40)
									fd["overlapped3_40"] = mergeStrArray([x['name'] for x in featuresIn3])
									t5 = re.findall(r'ttttt',str(fd["seq3"]),re.I)
									fd["ttttt"] = "true" if t5 else "false"
									t5_40 = re.findall(r'ttttt',str(fd["seq3_40"]),re.I)
									fd["ttttt_40"] = "true" if t5_40 else "false"


									tRnaNearby = find5(seq3Loc,record.features,'tRNA',True)
									tRnaNearby = tRnaNearby + find5(seq5LocExt,record.features,'tRNA',True)


									fd["dirty"] = 0


									seq5Final = str(seq5LocExt.extract(record).seq)
									fd["seq5_clean"] = ""

									if fd["overlapped5"]!="":
										fd["dirty"]+=1
										#generate clean seq5
										for fea in featuresIn5:
											if fd["strand"] != fea["strand"]:
												if fd["strand"]>0:
													startCodon = (fea["end"]-3,fea["end"])
													if str(Bio.SeqFeature.FeatureLocation(startCodon[0],startCodon[1],fd["strand"]).extract(record).seq) == "CAT" and startCodon[0]>=seq5LocExt.start and startCodon[1]<=seq5LocExt.end:
														fd["seq5_clean"] = str(Bio.SeqFeature.FeatureLocation(seq5LocExt.start,startCodon[0],feature.location.strand).extract(record).seq)+"TTA"+str(Bio.SeqFeature.FeatureLocation(startCodon[1],seq5LocExt.end,feature.location.strand).extract(record).seq)
														len1 = len(seq5Final)
														seq5Final = fd["seq5_clean"]
														if len1 != len(seq5Final):
															print "clean %d->%d"%(len1, len(seq5Final)) 
															print "%d %d %d %d"%(seq5LocExt.start,startCodon[0],startCodon[1],seq5LocExt.end)
												else:
													startCodon = (fea["start"],fea["start"]+3)
													if str(Bio.SeqFeature.FeatureLocation(startCodon[0],startCodon[1],fd["strand"]).extract(record).seq) == "CAT" and startCodon[0]>=seq5LocExt.start and startCodon[1]<=seq5LocExt.end:
														fd["seq5_clean"] = str(Bio.SeqFeature.FeatureLocation(startCodon[1],seq5LocExt.end,feature.location.strand).extract(record).seq)+"TTA"+str(Bio.SeqFeature.FeatureLocation(seq5LocExt.start,startCodon[0],feature.location.strand).extract(record).seq)
														seq5Final = fd["seq5_clean"]
														len1 = len(seq5Final)
														if len1 != len(seq5Final):
															print "clean %d->%d"%(len1, len(seq5Final)) 
															print "%d %d %d %d"%(seq5LocExt.start,startCodon[0],startCodon[1],seq5LocExt.end)
												fd["seq5P"] = fd["seq5_clean"]
											
									if fd["overlapped3"]!="":
										fd["dirty"]+=2
									if fd["ttttt"] == "false":
										fd["dirty"]+=4
									if soloLocs != []:
										fd["seq5_clean_solo"] = removeSolo(seq5Final,record,soloLocs)
										fd["seq5P"] = fd["seq5_clean_solo"]
										print fd["tRnaName"]
										fd["dirty"]+=8
									else:
										fd["seq5_clean_solo"] = ""
									if tRnaNearby != []:
										fd["dirty"]+=16
										print "dirty 16"
										



									fr.append(fd)
							lastRecord = record
						f.close()
			self.tRnaData[specie] = fr

	def loadUsedTable(self):
		self.matched = {}
		with open("matchedTRNA.json") as f:
			self.matched = json.loads(f.read())
	def exportCsvSC(self):
		specie = "SC"
		targetFile = open("%s.csv"%specie,"w")
		targetFile.write("id,tRNA Name,Original Sequence,Orintation,Anticodon,Amino Acid,Position,Sequence without intron,Intron,Intron Sequence,tRNA matched ,\n")
		for fd in self.tRnaData[specie]:	
			if self.matched.has_key(fd['tRnaName']):
				tRnaMatched = self.matched[fd['tRnaName']]
			else:
				tRnaMatched = ""
			targetFile.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n"%(fd['id'],fd['tRnaName'],fd['original_sequence'],fd['orientation'],fd['anticodon'],fd['aa'],fd['position'],fd['sequence'],fd['intron'],fd['intronSeq'],tRnaMatched))
		targetFile.close()
		
	def exportCsv(self):
		for specie in self.species:
			if specie == "SC":
				self.exportCsvSC()
			else:
				targetFile = open("%s.csv"%specie,"w")
				targetFile.write("id,tRNA Name,Original Sequence,Orintation,Anticodon,Amino Acid,Position,Sequence without intron,Intron,Intron Sequence,5,3,ttttt,feature at 5, feature at 3,3_40,ttttt_40,feature at 3_40,dirty type,seq5_clean,seq5_clean_solo\n")
				for fd in self.tRnaData[specie]:	
					targetFile.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n"%(fd['id'],fd['tRnaName'],fd['original_sequence'],fd['orientation'],fd['anticodon'],fd['aa'],fd['position'],fd['sequence'],fd['intron'],fd['intronSeq'],fd['seq5'],fd['seq3'],fd['ttttt'],fd['overlapped5'],fd['overlapped3'],fd['seq3_40'],fd['ttttt_40'],fd['overlapped3_40'],showDirtyType(fd["dirty"]),fd["seq5_clean"],fd["seq5_clean_solo"]))
				targetFile.close()

	#auto match tRNA from AG to SC:
	def autoMatch(self,specie):

		for fd in self.tRnaData[specie]:
			if fd['dirty'] & (4+16) == 0: #TTTTT and nearby is false, that means it is clean or cleared.
				anticodon = fd['anticodon']
				#filter same anticodon
				sameACinSC = [f for f in self.tRnaData['SC'] if f['anticodon'] == anticodon and not self.matched.has_key((f['tRnaName'])) and not self.matched.values().__contains__(fd['tRnaName'])]

				if len(sameACinSC)>0:
					self.matched[sameACinSC[0]['tRnaName']] = fd['tRnaName']
					sameACinSC[0]['matched'][specie] = fd['tRnaName']
					if fd['tRnaName'] == "AGOS_t0154":
						print "mat"

				elif self.matched.values().__contains__(fd['tRnaName']):
					print "XXXXXXXXXXXXXXXXXXXXX"
				else:
					self.soloSet[specie].append(fd['tRnaName'])
					if fd['tRnaName'] == "AGOS_t0154":
						print "Aaadd"

	def autoMatchAgain2(self,specie):
		for name in self.soloSet["SC"]:
			sc = [x for x in self.tRnaData['SC'] if x['tRnaName']== name][0]
			other = [x for x in self.tRnaData[specie] if x['aa'] == sc['aa'] and self.soloSet[specie].__contains__(x['tRnaName'])]
			if other != []:
				otherName = other[0]['tRnaName']
				sc['matched'][specie] = otherName
				self.matched[name] = otherName
				self.soloSet["SC"].remove(name)
				self.soloSet[specie].remove(otherName)
				print name, otherName
				
	def autoMatchFinal3(self,specie):
		print "final---"+specie
		print "sc = %d"%len(self.soloSet["SC"])
		print specie+" = %d"%len(self.soloSet[specie])
		scs = [x for x in self.tRnaData['SC'] if self.soloSet["SC"].__contains__(x['tRnaName'])]

		for sc in scs: 
			name = sc['tRnaName']
			if self.soloSet[specie] == []:
				print "target null"
				return
			otherName = self.soloSet[specie][0]
			sc['matched'][specie] = otherName
			self.matched[name] = otherName
			self.soloSet["SC"].remove(name)
			self.soloSet[specie].remove(otherName)
			print name, otherName

	
	def findSoloSC(self):
		a = [x['tRnaName'] for x in self.tRnaData['SC']]
		b = self.matched.keys()
		self.soloSet["SC"] = [x for x in a if x not in b]

	
	def statSolo(self,specie):
		targetFile = open("%sleft.csv"%specie,"w")
		targetFile.write("anticodon,aminoacid,count,\n")
		anticodonList = [x['anticodon'] for x in self.tRnaData[specie]]
		aaList = dict([(x['anticodon'],x['aa']) for x in self.tRnaData[specie]])
		count = dict.fromkeys(anticodonList,0)
		for name in self.soloSet[specie]:
			n = [x for x in self.tRnaData[specie] if x['tRnaName']==name][0]
			anticodon = n['anticodon']
			count[anticodon]+=1
		
		for ac in count:
			if count[ac] >0:
				targetFile.write("%s,%s,%d,\n"%(ac,aaList[ac],count[ac]))

		targetFile.write(",,,,\n");
		for name in self.soloSet[specie]:
			n = [x for x in self.tRnaData[specie] if x['tRnaName']==name][0]
			anticodon = n['anticodon']
			targetFile.write("%s,%s,,\n"%(anticodon,name))


		targetFile.close()

	def genDescription(self,sc,ag):
		s = "The 5` and 3` tRNA of %s have been replaced by %s. "%(sc['tRnaName'],ag['tRnaName'])
		if(sc['sequence'] != sc['original_sequence']):
			s+= "Intron has been removed. "
		d = ag['dirty']
		if d & 1 >0:
			s+="Unwanted gene has been deleted from 5`. "
		if d & 2 >0:
			s+="Unwanted gene has been deleted from 3`. "
		if d & 8 > 0:
			s+="SoloLTR has been deleted."
		return s

				
	def mergeAll(self):
		self.tRnaData['SC'].sort(key=lambda x: x['id'])
		for i in range(1,17):	
			base = 0
			roxStr = "taactttaaataatgccaattatttaaagtta"
			rox = Bio.Seq.Seq(roxStr,Bio.Alphabet.DNAAlphabet())
			roxLen = len(roxStr)
			seq = ""
			seq+=rox
			frLoc =  Bio.SeqFeature.FeatureLocation(base,base+roxLen,0)
			frox = Bio.SeqFeature.SeqFeature(frLoc,"misc_feature",strand=0,qualifiers={"note":"rox recombination site","color":"color: #ff0000"})
			features = [frox]
			base+=roxLen


			for sc in [x for x in self.tRnaData['SC'] if x['chromosome']==i]:
				if self.matched.has_key(sc['tRnaName']):
					matchedName = self.matched[sc['tRnaName']]
					targetSpecie = matchedName[0:2].upper()
					replaceTRna = [x for x in self.tRnaData[targetSpecie] if x['tRnaName']==matchedName][0]
					neoTRna = Bio.Seq.Seq(replaceTRna['seq5P']+sc['sequence']+replaceTRna['seq3'],Bio.Alphabet.DNAAlphabet())
					if sc["orientation"] == '-':
						strand = -1
						neoTRna = neoTRna.reverse_complement()
						f3Len = len(replaceTRna['seq3'])
						f3Loc = Bio.SeqFeature.FeatureLocation(base,base+f3Len,strand)
						base+= f3Len
						fsLen = len(sc['sequence'])
						fsLoc = Bio.SeqFeature.FeatureLocation(base,base+fsLen,strand)
						base+= fsLen
						f5Len = len(replaceTRna['seq5P'])
						f5Loc = Bio.SeqFeature.FeatureLocation(base,base+f5Len,strand)
						base+= f5Len
					else:
						strand = 1
						f5Len = len(replaceTRna['seq5P'])
						f5Loc = Bio.SeqFeature.FeatureLocation(base,base+f5Len,strand)
						base+= f5Len
						fsLen = len(sc['sequence'])
						fsLoc = Bio.SeqFeature.FeatureLocation(base,base+fsLen,strand)
						base+= fsLen
						f3Len = len(replaceTRna['seq3'])
						f3Loc = Bio.SeqFeature.FeatureLocation(base,base+f3Len,strand)
						base+= f3Len
					f3 = Bio.SeqFeature.SeqFeature(f3Loc,"3'UTR",strand=f3Loc.strand,qualifiers={"note":replaceTRna['id']+"_3"})
					f5 = Bio.SeqFeature.SeqFeature(f5Loc,"5'UTR",strand=f5Loc.strand,qualifiers={"note":replaceTRna['id']+"_5"})
					fs = Bio.SeqFeature.SeqFeature(fsLoc,"tRNA",strand=fsLoc.strand,qualifiers={"note":sc['tRnaName'],"desc":self.genDescription(sc,replaceTRna)})
					features.append(f5)
					features.append(fs)
					features.append(f3)

					frLoc =  Bio.SeqFeature.FeatureLocation(base,base+roxLen,0)
					frox = Bio.SeqFeature.SeqFeature(frLoc,"misc_feature",strand=0,qualifiers={"note":"rox recombination site","color":"color: #ff0000"})
					features.append(frox)

					seq+=neoTRna
					base+=roxLen
					seq+=rox

			record = Bio.SeqRecord.SeqRecord(seq,'17','Exported','the neochromosome',["Exported 8 Dec 2014 from SnapGene Viewer 2.6.0"],features=features)
			with open("neochromosome%02d.gb"%i,"w") as f:
				Bio.SeqIO.write(record,f,"genbank")
			
			#convert to gbk-snapgene format
			with open("neochromosome%02d.gb"%i,"r") as f,open("neochromosome%02d.gbk"%i,"w") as f2:
				a = f.readlines()
				for line in a:
					if line[0:6] == "DBLINK":
						f2.write("REFERENCE   1\n  AUTHORS   Isaac Luo\n  TITLE     Direct Submission\n  JOURNAL"+line[6:])
					else:
						line = line.replace(r"/color=",r"/note=")
						f2.write(line.replace(r"/desc=",r"/note="))
				


			
		
	
neo = Neochromosome()
neo.importData()
neo.loadUsedTable()
neo.autoMatch("AG")
neo.autoMatch("EC")
neo.findSoloSC()
#neo.autoMatch("KL")
print len(neo.soloSet['SC'])
print len(neo.soloSet['AG'])
print len(neo.soloSet['EC'])
neo.statSolo('SC')
neo.statSolo('AG')
neo.statSolo('EC')
neo.autoMatchAgain2('AG')
neo.autoMatchAgain2('EC')
print len(neo.soloSet['SC'])
print len(neo.soloSet['AG'])
print len(neo.soloSet['EC'])
neo.autoMatchFinal3('AG')
neo.autoMatchFinal3('EC')
print neo.soloSet['SC']
print neo.soloSet['AG']
print neo.soloSet['EC']
print len(neo.soloSet['SC'])
print len(neo.soloSet['AG'])
print len(neo.soloSet['EC'])
neo.mergeAll()

neo.saveData()
neo.exportCsv()
		

