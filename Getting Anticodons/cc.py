import sys
import os

codonDict = {"AAA":"Lys", "AAC":"Asn", "AAG":"Lys", "AAT":"Asn", "ACA":"Thr", "ACC":"Thr", "ACG":"Thr", "ACT":"Thr", "AGA":"Arg", "AGC":"Ser", "AGG":"Arg", "AGT":"Ser", "ATA":"Ile", "ATC":"Ile", "ATG":"Met", "ATT":"Ile", "CAA":"Gln", "CAC":"His", "CAG":"Gln", "CAT":"His", "CCA":"Pro", "CCC":"Pro", "CCG":"Pro", "CCT":"Pro", "CGA":"Arg", "CGC":"Arg", "CGG":"Arg", "CGT":"Arg", "CTA":"Leu", "CTC":"Leu", "CTG":"Leu", "CTT":"Leu", "GAA":"Glu", "GAC":"Asp", "GAG":"Glu", "GAT":"Asp", "GCA":"Ala", "GCC":"Ala", "GCG":"Ala", "GCT":"Ala", "GGA":"Gly", "GGC":"Gly", "GGG":"Gly", "GGT":"Gly", "GTA":"Val", "GTC":"Val", "GTG":"Val", "GTT":"Val", "TAA":"*", "TAC":"Tyr", "TAG":"*", "TAT":"Tyr", "TCA":"Ser", "TCC":"Ser", "TCG":"Ser", "TCT":"Ser", "TGA":"*", "TGC":"Cys", "TGG":"Trp", "TGT":"Cys", "TTA":"Leu", "TTC":"Phe", "TTG":"Leu", "TTT":"Phe",} 
anticodonDict = {"AAA":"Lys", "AAC":"Asn", "AAG":"Lys", "AAT":"Asn", "ACA":"Thr", "ACC":"Thr", "ACG":"Thr", "ACT":"Thr", "AGA":"Arg", "AGC":"Ser", "AGG":"Arg", "AGT":"Ser", "ATA":"Ile", "ATC":"Ile", "ATG":"Met", "ATT":"Ile", "CAA":"Gln", "CAC":"His", "CAG":"Gln", "CAT":"His", "CCA":"Pro", "CCC":"Pro", "CCG":"Pro", "CCT":"Pro", "CGA":"Arg", "CGC":"Arg", "CGG":"Arg", "CGT":"Arg", "CTA":"Leu", "CTC":"Leu", "CTG":"Leu", "CTT":"Leu", "GAA":"Glu", "GAC":"Asp", "GAG":"Glu", "GAT":"Asp", "GCA":"Ala", "GCC":"Ala", "GCG":"Ala", "GCT":"Ala", "GGA":"Gly", "GGC":"Gly", "GGG":"Gly", "GGT":"Gly", "GTA":"Val", "GTC":"Val", "GTG":"Val", "GTT":"Val", "TAA":"*", "TAC":"Tyr", "TAG":"*", "TAT":"Tyr", "TCA":"Ser", "TCC":"Ser", "TCG":"Ser", "TCT":"Ser", "TGA":"*", "TGC":"Cys", "TGG":"Trp", "TGT":"Cys", "TTA":"Leu", "TTC":"Phe", "TTG":"Leu", "TTT":"Phe",} 

with open('codon_stat.csv') as f, open('codon_stat_2.csv','w') as f2:
	
	lines = f.readlines()
	for line in lines[1:]:
		x = line.split(",")
		anticodon = x[0]
		aa = codonDict[anticodon]
		f2.write(aa)
		f2.write(",")
		f2.write(line)
