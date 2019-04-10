# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 14:42:08 2019

@author: jamesbutterworth
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 10:17:24 2019

@author: jamesbutterworth
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 09:57:16 2019

@author: jamesbutterworth
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 10:14:16 2019

@author: jamesbutterworth
"""


import pysam
import pandas
import pickle
from collections import Counter
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
from Bio.Seq import MutableSeq 
from Bio.Alphabet import generic_dna
#from Bio import PatternIO
from Bio.SeqRecord import SeqRecord
import re
import os
import time
import pprint
import sys
import csv
import operator
import csv
from Bio import Align
from Bio.pairwise2 import format_alignment
from Bio import SeqIO



import numpy as np
from Bio import SeqIO
import datetime
from operator import add
#import StringBuilder
import pandas
from collections import OrderedDict
from Bio import pairwise2
#samfile = pysam.AlignmentFile("referenceSY1.bam", "rb")
#for read in samfile.fetch('NC_003888', 0, 100):
#for pileupcolumn in samfile.pileup("NC_003888", 0, 100):
#    for pileupread in pileupcolumn.pileups:
#        base=pileupread.alignment.query_sequence      
#print base
#for read in samfile.fetch('NC_003888', 0, 100):
#     print read
#import pysamstats

bamFile = "referenceSY4.bam"  ## Reads in bam file for the strain of interest as well as the reference fasta file.
referenceFile = "reference4.fasta"
filename = "GCF_000203835.1_ASM20383v1_genomic.gbff"

nums=[] ## Defines parameters
Sym=[]
baseRef=0
locus_to_gene = dict()
loci={}
Chromosome_name= 'NC_003888'
num_index=0
AInsert=None
newDele=0
Upos=-1
start = time.time()
NoCov=[]
ran=1000
deletionsRefPos =[]
deletionsStrainPos=[]
PointPos=[]
PointBefore=[]
PointAfter=[]
PointCoverage=[]
PointStrainPos=[]
all_count={}
deletionsBase=[]
insertionsRefPos=[]
insertionsBase=[]
insertionsStrainPos=[]
query1=[]
query2=[]
query3=[]
count=[]
reference_positions = []
nsegments=[]
reference_id=[]
reference_name=[]
bases=[]
tmpI=[]
Sym2=[]
GeneName=[]
w=0
FrequencyDict={}
CoverageDict={}
DeletionPosition=None
InsertionPosition=None
i=93
n_dels=0
x=0
StrainPosition=0
reference = []
df=pandas.DataFrame()
dfINT=pandas.DataFrame()
dfTheirs=pandas.DataFrame()
dfTheirsLE=pandas.DataFrame()
dfSym=pandas.DataFrame()
deletion_symbol='*'
n_ins=0
SynonymousStatus=""
Gene_name=""
gene=""
geneD={}
proD={}
locD={}
product=[]
A=""
B=""

                        

#for record in SeqIO.parse(filename, "genbank"): ## reads genebank file for Streptomyces Coelicolor found on NCBI website.
 #   for f in record.features: ##Looks through file and strips to leave only information about genes.  Data about location and gene name is extracted and placed in a dictionary.
  #              if "product" in f.qualifiers:
   #                 locus_tags = f.qualifiers["locus_tag"]
    #                if locus_tags==['SCO0022']:
     #                  None
      #              else:
       #                 product_tag=f.qualifiers["product"]  
        #                if product_tag=="":
         #                   None
          #              else:
           #                 proc=str(product_tag)
            #                proc=proc[:-2]
             #               proc=proc[2:]
              #              loc=f.location
  #                          loc=loc[:-3]
#                            loc=loc[:-1]
 #                           loc=loc[1:]
   #                         loc=loc.replace(">","")
    #                        loc=loc.replace("<","")
     #                       loc=loc.split(':')
      #                      locr=range(int(loc[0]), int(loc[1]))
       #                     proD[proc]=locr
                        
                        
                        
                        ##if locus_tags==['SCO5824']:
                         ##   print locr
                        
                        
             #               lokL=str(locus_tags)
                            ##proc=str(product_tag)
              #              lokL=lokL[:-2]
               #             lokL=lokL[2:]
                            
                            ## for x,y in loci.items():
                            ##    if locus_tags in x:
                            ##       loci[loki]=locr, y
                            
                            ## loci[y].append(locr)
                          #  locD[lokL]=locr
                           # if locus_tags==['SCO0126']:
                            #   print product_tag
                            #product_tag=""
                #if "locus_tag" in f.qualifiers:
                 #   locus_tags = f.qualifiers["locus_tag"]
                    
                  #  if locus_tags==['SCO0022']:
                   #    None
                    #else:
                     #   z=""
                      #  loc=f.location
                       # loc=str(loc)
                        #loc=loc[:-3]
                        #loc=loc[:-1]
                        #loc=loc[1:]
                        #loc=loc.replace(">","")
                        #loc=loc.replace("<","")
                        #loc=loc.split(':')
                        #locr=range(int(loc[0]), int(loc[1]))
                        
                        ##if locus_tags==['SCO5824']:
                         ##   print locr
                        
                        
                        #loki=str(locus_tags)
                        ##proc=str(product_tag)
                        #loki=loki[:-2]
                        #loki=loki[2:]
                       
                       ## for x,y in loci.items():
                        ##    if locus_tags in x:
                         ##       loci[loki]=locr, y
                       
                       ## loci[y].append(locr)
                        #loci[loki]=locr
                        ##proD[proc]=locr
                        
                        
                        
                        
                        
                        
                        
                        
                        #locr=None
                        
                        #loc=[]
                        #z=""
                        #if "gene" in f.qualifiers:
                         #   genes=f.qualifiers["gene"]
                          #  genei=str(genes)
                           # genei=genei[:-2]
                            #genei=genei[2:]
                            #loc=f.location
                            #loc=str(loc)
                            #loc=loc[:-3]
                            #loc=loc[:-1]
                            #loc=loc[1:]
                            #loc=loc.split(':')
                            #locr=range(int(loc[0]), int(loc[1]))
                            #geneD[genei]=locr





#for record in SeqIO.parse(filename, "genbank"): ## reads genebank file for Streptomyces Coelicolor found on NCBI website.
 #   for f in record.features: ##Looks through file and strips to leave only information about genes.  Data about location and gene name is extracted and placed in a dictionary.
  #              if "locus_tag" in f.qualifiers:
   #                 locus_tags = f.qualifiers["locus_tag"]
    #                if locus_tags==['SCO0022']:
     #                  None
      #              else:
       #                 z=""
        #                loc=f.location
         #               loc=str(loc)
          #              loc=loc[:-3]
           #             loc=loc[:-1]
            #            loc=loc[1:]
             #           loc=loc.replace(">","")
              #          loc=loc.replace("<","")
               #         loc=loc.split(':')
                #        locr=range(int(loc[0]), int(loc[1]))
                        
                        ##if locus_tags==['SCO5824']:
                          ##  print locr
                        
                        
                 #       loki=str(locus_tags)
                  #      loki=loki[:-2]
                   #     loki=loki[2:]
                    #    for x,y in loci.items():
                     #       if loki in x and locr not in y:
                      #          h=locr + y                                
 #                               loci[loki]=h
#                        if "gene" in f.qualifiers:
  #                          genes=f.qualifiers["gene"]
   #                         genei=str(genes)
    #                        genei=genei[:-2]
     #                       genei=genei[2:]
      #                      loc=f.location
 #                           loc=loc[:-3]
#                            loc=loc[:-1]
  #                          loc=loc[1:]
   #                         loc=loc.split(':')
    #                        locr=range(int(loc[0]), int(loc[1]))
                            #geneD[genei]=locr
                        
     #                       for x,y in geneD.items():
      #                          if genei in x and locr not in y:
       #                             h=locr + y                                
        #                            geneD[genei]=""
         #                           geneD[genei]=h
          #              if "product" in f.qualifiers:
           #                 product_tag=f.qualifiers["product"]  
            #                proc=str(product_tag)
             #               proc=proc[:-2]
              #              proc=proc[2:]
               #             loc=f.location
                #            loc=str(loc)
                 #           loc=loc[:-3]
                  #          loc=loc[:-1]
 #                           loc=loc.replace(">","")
#                            loc=loc.replace("<","")
  #                          loc=loc.split(':')
   #                         locr=range(int(loc[0]), int(loc[1]))
    #                        lokL=str(locus_tags)
                        ##proc=str(product_tag)
     #                       lokL=lokL[:-2]
      #                      lokL=lokL[2:]
                                
       #                     for x,y in proD.items():
        #                        if proc == x and locr not in y:
         #                           h=locr + y                                
          #                          proD[proc]=""
           #                         proD[proc]=h
            #                for x,y in locD.items():
             #                   if lokL in x and locr not in y:
              #                      h=locr + y                                
               #                     locD[lokL]=""
                #                    locD[lokL]=h
#OnlyInTheirsSY1=[238465,303359,335691,363530,461702,480026,542144,567543,615966,617040,685256,685260,898700,954304,975289,1120938,1172685,1196433,1305955,1305958,1577555,1616264,1619309,1799094,1881782,1914409,1940745,1940747,1990806,2047555,2047808,2187131,2187145,2254778,2262762,2366366,2463719,2466884,2494047,2903373,3134658,3244356,3315219,3372057,3425798,3442231,3677481,3681457,3784285,3799532,3820495,3883759,4120143,4182372,4242460,4335562,4368010,4482704,4854072,4886043,5015343,5016697,5017102,5017391,5017659,5017856,5018196,5018258,5018303,5179038,5284407,5309494,5640821,5857376,5872646,5892942,5916923,5970641,6175175,6276084,6296134,6296845,6303814,6323795,6332825,6411145,6430907,6494953,6494957,6652183,6685111,6756626,6894463,6894868,6895425,6895622,6895962,6896621,6897032,7013124,7125988,7188923,7199898,7236674,7347280,7752588,7790124,7812884,7932765,8045464,8209958,8209987,8351167,8527898,8638333]

#OnlyInTheirsSY2=[615963,757255,804152,1120933,1734126,2254778,2262781,2463719,2476629,3000905,3677481,4335562,4854072,5015343,5872646,5906039,5916923,6008595,6235401,6319425,6946819,7188919,7188923,7236674,7790124,8388985,8493775,480026,1881782,2187145,3000902,5284407,5640821,88679,685256,1267083,1305955,1305958,1619309,2047555,3246641,3639579,3681457,3799532,3883779,4886043,5234366,5234371,5234379,5892942,7236677,7614529,5098324,6017682,3784275,1327295,3600287,5017391,5017856,6895157,6895355,6895622,542144,898700,1120938,1230280,1577555,1914409,1940745,2366366,2452318,2707231,3033245,3246623,3425798,3883759,4242460,4482704,4635208,5017102,5017129,5018258,5018303,5608083,5693587,5697263,6276084,6319416,6894868,6894895,6896024,7347280,8351167,8499568,8527898,4175173,6090133,2959823,5018196,6895962,2254780,5016697,6296134,6296845,6652183,7199898,102908,238465,335691,363530,685260,1940747,2466884,2494047,3244356,3315219,3393392,3960557,4182372,4407626,4630733,5309494,5443070,5536064,5560348,6895425,7483927,7939516,8638333,5479806]

#OnlyInTheirsSY3=[130,335691,363530,449721,480026,542144,552452,552461,617040,620758,685256,685260,804152,898700,954304,1110675,1120938,1172685,1267083,1305955,1305958,1530460,1592977,1619309,1658179,1658187,1734126,1799094,1881782,1914409,1940745,1940747,1990806,2047555,2138835,2187131,2187145,2254778,2254780,2354209,2403736,2458983,2463719,2466884,2494047,2608348,2707231,2889536,3054640,3244356,3246623,3315219,3425798,3666746,3677481,3681457,3799532,3820495,3926305,3960557,4182372,4242460,4335562,4347064,4368010,4377131,4407626,4482704,4527587,4638934,4638936,4692054,4854072,4886043,5015343,5017391,5017589,5017856,5018303,5169938,5234366,5284407,5309494,5312790,5535165,5560348,5640821,5693587,5697263,5762043,5857376,5892942,5916923,6017682,6175175,6223482,6276084,6296134,6296904,6296907,6393402,6430907,6494953,6494957,6652183,6756626,6792528,6850006,6894463,6895157,6895622,6895962,6896621,6897032,7155242,7188919,7188923,7199898,7236674,7236677,7313675,7347280,7483927,7614529,7752588,7790124,7929748,7939516,7958526,8060637,8140456,8209958,8209987,8264168,8308247,8351167,8388985,8493775,8527898,8562221,8638333,8666781]

#OnlyInTheirsSY4=[88679,259647,303359,449721,542144,564099,615966,615974,617040,620758,679307,685256,685260,804152,898700,954304,1110675,1120933,1120938,1196433,1230280,1267083,1305958,1353701,1353825,1577555,1616345,1616354,1619309,1649434,1696166,1734126,1799094,1881782,1914409,1940747,2042340,2047555,2187131,2187145,2222038,2262781,2339352,2354209,2366366,2370626,2458948,2463719,2466884,2494047,2521464,2707231,2889536,2889541,3033245,3147715,3246623,3246641,3315219,3334861,3425798,3677481,3698423,3784275,3784276,3784285,3799532,3820495,3883779,3960557,3987880,4182372,4242460,4335562,4347064,4377131,4407626,4527587,4630733,4854072,4886043,5015343,5016697,5017102,5017129,5017391,5017589,5017659,5017856,5018196,5082759,5284407,5293858,5309494,5443070,5535165,5575580,5640821,5693587,5706577,5769289,5857376,5872646,5873901,5873903,5892942,5916923,5970641,6017682,6175175,6276084,6296134,6296845,6323795,6652183,6792528,6894463,6894868,6895157,6895355,6895425,6895622,6895962,7155242,7188923,7199898,7236674,7347280,7483943,7564996,7614529,7790124,7939516,8045464,8060637,8194142,8209958,8209987,8303825,8388985,8527898,8638333,8666781]
##Creates 4 lists which contain the positions of polymorphisms only found in their results.

fasta_sequences=SeqIO.parse(open("reference4.fasta"), 'fasta') ##Opens fasta file and the strain of interest's bam file as a sam file. 
samfile = pysam.Samfile("referenceSY4.bam", "rb")


start=0   ## Intial start and stop sites for iterating through positions of the bam file.
stop=8667506

#deletion 4438967
#8667506
#Insertion 3054000-3055000

iterator = samfile.pileup( 'NC_003888',start,stop) # what are start and stop

for base in fasta_sequences: ##Reads through bases within the fasta file and converts them all to lower case.
    seq=base.seq 
    seql= seq.lower()
seql= seql.tomutable() ##Converts fasta sequence into mutable string

colnames = ['X', 'Base.After', 'Base.Before', 'Base.Coverage', 'Bases.After', 'Bases.Before', "Frequency.of.Mutation", "Overall Coverage", "Mutation Coverage", "Position.in.Strain", "RefPos", "Product.of.Rank", "Rank.of.Coverage", "Rank.of.Frequency", "Status", "Sum.of.Rank", "Amino Before", "Amino After" , "Synonymous Status", "Gene", "Locus Tag", "Chromosome Name"]
#data = pandas.read_csv('SY4OursPY.csv', names=colnames)  ##Retrieves and formats data from csv containg polymorphisms that are only observed in our results.
#pos = data.RefPos.tolist()
#del pos[0]

#SeqA=seq.translate(table=11) ##Translates fasta sequence using the genetic code for Streptomyces Coelicolor


#f=open('SY4SequenceTest.txt','r') ##Opens Sequence produced from script previously.
#for x in f: ##Converts sequence into String and formats it so is suitable for translation.
#    y=str(Seq(x))
#h=y.replace('*','')
#z=h.replace('-','N')
#dna=Seq(z)
#SeqRA=dna.translate(table=11)

seql[3521520:3600277]="-"*78757 ##Converts the base at these specific positions into "-" representing known deleted regions.
seql[5515933:5532917]="-"*16984
seql[6434878:6462139]="-"*27261
seql[6896194:6946387]="-"*50193


for pileupcolumn in iterator: ##Creates for loop for iterating through postions of reference file.
    
    ref_pos = pileupcolumn.reference_pos ##Positions of refernce file.
    #for x in OnlyInTheirs:
     #   if x !=pileupcolumn.reference_pos:
      #      print x
 #   if 3521520<=ref_pos<=3600277:
  #      seql[ref_pos]="-"
   # if 5515933 <=seql[ref_pos]<=5532917:
    #    seql[ref_pos]="-"
    #i#f 6434878 <=seql[ref_pos]<=6462139:
      #  seql[ref_pos]='-'
    #if 6896194 <= seql[ref_pos]<=6946387:
     #   seql[ref_pos]='-'
         
    C_ins_count=0 ##Defining parameters
    G_ins_count=0
    A_ins_count=0
    T_ins_count=0
    del_count=0
    C_count=0
    G_count=0
    A_count=0
    T_count=0
    Status=''
    
    for pileupread in pileupcolumn.pileups: ##Looks at readds for positioin in bam file
        flag_change=False
        if pileupread.is_del==1: ##Identifies if deletion has occured.
            del_count=del_count+1
        else: ##Count the coverage for indivdual bases
            if pileupread.alignment.query_sequence[pileupread.query_position]=='C':
              C_count=C_count+1  
            if pileupread.alignment.query_sequence[pileupread.query_position]=='G':
              G_count=G_count+1
            if pileupread.alignment.query_sequence[pileupread.query_position]=='A':
              A_count=A_count+1
            if pileupread.alignment.query_sequence[pileupread.query_position]=='T':
              T_count=T_count+1
        if pileupread.indel==1: ##Identifies if insertion has taken place and counts coverage for insertion
            if pileupread.alignment.query_sequence[pileupread.query_position+1]=='C':
              C_ins_count=C_ins_count+1  
            if pileupread.alignment.query_sequence[pileupread.query_position+1]=='G':
              G_ins_count=G_ins_count+1
            if pileupread.alignment.query_sequence[pileupread.query_position+1]=='A':
              A_ins_count=A_ins_count+1
            if pileupread.alignment.query_sequence[pileupread.query_position+1]=='T':
              T_ins_count=T_ins_count+1
          
    
   # for x in OnlyInTheirsSY4: ##Identifies if polymorphisms only found in theirs only occur in regions of zero coverage. 
    #    if x != ref_pos+1:
     #      None 
      #  else:
       #     NoCov.append(x)

        
    all_count = {'C':C_count,'G':G_count,'T':T_count,'A':A_count}   ##Places base coverage in dictionary. 
    cov= C_count+G_count+T_count+A_count 
    StrainBase=max(all_count.iteritems(),key=operator.itemgetter(1))[0] ##Identifes base with the heighest amount of coverage.
    Frequency=max(all_count.values()) ## Identifies the value associated with the base that has heighest coverage.
    
    if all_count.get(StrainBase)>.7*cov: ## Identifies Point Mutations
        seql[ref_pos+n_ins]=StrainBase.upper() ##If the coverage of the StrainBase is greater then 70% then write the StrainBase to the Sequence.
        
        if seq[ref_pos] is not StrainBase: ## If the StrainBase is not the same as the base found at that position within the orignal sequecne then mark it as a point mutation.
            flag_change=True  
            Status=Status + 'Point mutation'
            FrequencyDict.update({ref_pos: Frequency})
            CoverageDict.update({ref_pos: cov})
            BaseAfter=StrainBase
            BaseBefore=seq[ref_pos]
            OverallCov=cov
            
		

                    
                
    
    if all_count.get(StrainBase)<=.7*cov: ##If the StrainBase coverage is less than 70% then use the base located within that position of the reference.
        seql[ref_pos+n_ins]=seq[ref_pos].upper()
    
    if del_count > 1.5*cov: ## Identifies if deletions have taken place, but will only accept deletion if the count is greater then 1.5 * the coverage.
        seql[ref_pos+n_ins] = deletion_symbol
        flag_change=True
        Status=Status + 'Deletion'
        Frequency=del_count
        FrequencyDict.update({ref_pos: Frequency})
        CoverageDict.update({ref_pos: cov})
        BaseBefore=seq[ref_pos]
        BaseAfter=deletion_symbol
        OverallCov=cov+del_count
        n_dels=n_dels-1
        
		
    
    all_count_ins = {'C':C_ins_count,'G':G_ins_count,'T':T_ins_count,'A':A_ins_count}
    StrainBaseIns=max(all_count_ins.iteritems(),key=operator.itemgetter(1))[0]
    if all_count_ins.get(StrainBaseIns)>.7*cov: ## Identifies if Insertions have taken place, again will only accept if the insertion count is greater then 70% of the coverage.
        seql=seql[:(ref_pos+n_ins+1)] + StrainBaseIns + seql[(ref_pos+n_ins+1):]
        n_ins=n_ins+1
        flag_change=True
        Frequency=all_count_ins.get(StrainBaseIns)
        Status=Status +'Insertion'
        FrequencyDict.update({ref_pos: Frequency})
        CoverageDict.update({ref_pos: cov})
        BaseBefore="-"
        BaseAfter=StrainBaseIns
        OverallCov=cov+all_count_ins.get(StrainBaseIns)
        
		
    #word=all_count.get(StrainBase) ##Used to identify point of questionable polymorphisms found in their results, this particular loop is used to identify a polymorphism that is observed when the base coverage is 50% without the reference.
    #for word in all_count:
     #   w=w+1
   # pos2=[]
    #for x in pos: ## Converts the postions identifed from Our data set into useable integers.
     #   pos2.append(int(x))
    
        
        
        
   
    #if ref_pos+1 in pos2: ##This if staement is used to identify synonymous and non-synonymous changes 
     #   z= float(ref_pos+1) / 3
      #  p=float(ref_pos+1)
       # if (z).is_integer()==True:
        #    z=z-1
        #A=SeqA[int(ref_pos)/3]
        #B=SeqRA[int((ref_pos+n_ins+n_dels)/3)]
        #print ref_pos+1
        #print A
        #print B
        #dfSym=dfSym.append({"Reference point": ref_pos+1, "Strain position":ref_pos+1+n_ins+n_dels, "Amino before":A, "Amino after": B, "Number of insertions":n_ins, "Number of deletions": n_dels }, ignore_index=True)
        #dfSym.to_csv("SymTest1INTTESTOSY4GeneTEST.csv", sep='\t')
        
		
        #if A !=B: ##Creates two lists, one contiaing synonymous changes and one containing non-synonymous changes.
         #   Sym.append(ref_pos+1)
        #if A==B:
         #   Sym2.append(ref_pos+1)
    #if ref_pos+1 in OnlyInTheirsSY4: ## If statement used to split polymorphisms only found in theirs into 3 categories (Technically 4 with the ones found in zero coverage)
     #   ref_pos2=ref_pos + 1
      
      #  if all_count[seq[ref_pos]] >= all_count.get(StrainBase): ## 1st Category: Reference Base has the heighest coverage
         
       #     dfTheirs=dfTheirs.append({"Mutations where reference coverage is highest": ref_pos+1, "Coverage": cov, "Reference Coverage": all_count[seq[ref_pos]], "Base Frequnecy Count": all_count, "Reference Base": seq[ref_pos] }, ignore_index=True)
       # else: ## 2nd Category contains mutations of little evidence to support their existence.
           
        #    dfTheirsLE=dfTheirsLE.append({"Little evidence": ref_pos+1, "Reference Base": seq[ref_pos], "Little evidence Count": all_count}, ignore_index=True)
        
        #if all_count.get(StrainBase)==.5*cov and all_count[seq[ref_pos]]==0 and w>1: ##3rd Category contains Mutations where the base covergae of 2 StrainBases is 50% and neither are the reference base.
         #   dfTheirs=dfTheirs.append({"Mutation with only 50% coverage": ref_pos+1, "50% Coverage": cov, "50% Count": all_count, "Reference Base": BaseBefore }, ignore_index=True)
            
        
  
    
    if flag_change: ##Allows the updating of dataframes for insertions. point mutations and deletions.
                   
        None
        z= float(ref_pos+1) / 3
        p=float(ref_pos+1)
        if (z).is_integer()==True:
            z=z-1
        #A=SeqA[int(ref_pos)/3]
        #B=SeqRA[int((ref_pos+n_ins+n_dels)/3)]
        print ref_pos+1
       # print A
        #print B
        #if A==B:
        #    SynonymousStatus="2"
        #else:
         #   SynonymousStatus="1"
        #dfSym=dfSym.append({"Reference point": ref_pos+1, "Strain position":ref_pos+1+n_ins+n_dels, "Amino before":A, "Amino after": B, "Number of insertions":n_ins, "Number of deletions": n_dels }, ignore_index=True)
        #dfSym.to_csv("SymTest1INTTEST.csv", sep='\t')
        ref_pos2=ref_pos+1
        print ref_pos2
        print ref_pos2+n_ins
        print Status
        BasesAfter= seql[ref_pos+n_ins-5:ref_pos+n_ins+5]
        BasesBefore= seq[ref_pos-5:ref_pos+5]
        
        #for x, y in loci.items(): ##Used to identify wether the mutation takes place within a known gene or not.    
         #   if ref_pos in y:
          #      GeneName.append(x)
        #for x, y in geneD.items():
         #   if ref_pos in y:
          #      gene=x
        #for x,y in proD.items():
         #   if ref_pos in y:
          #      product.append(x)
       
        df=df.append({'Chromosome Name': Chromosome_name, 'Status': Status, "Position in reference": ref_pos2, "Position in Strain": ref_pos2+n_ins+n_dels, 'Status': Status, 'Bases Before': BasesBefore, "Bases After": BasesAfter, "Mutation Coverage": all_count, "Frequency of Mutation": Frequency, "Base Coverage": cov, "Overall Coverage": OverallCov,  "Base Before": BaseBefore, "Base After": BaseAfter, "Amino before":A, "Amino after": B, "Number of insertions":n_ins, "Number of deletions": n_dels, "Synonymous Status": SynonymousStatus, "Locus Tag": GeneName, "Gene": gene,"Product": product, "Frame Shift position": n_ins+n_dels }, ignore_index=True)
        df.to_csv("dfTestSY4GeneTESTFran.csv", sep='\t')       
        Gene_name=""
        gene=""
        product=[]
        GeneName=[]
        
                ##else: 
                  ##  df=df.append({"Gene": Gene_name}, ignore_index=True)
RE=OrderedDict(FrequencyDict) ##Creates an ordered Frequency Dictionary

        
        

RankedF={key: rank for rank, key in enumerate(sorted(FrequencyDict, key=FrequencyDict.get, reverse=True), 1)} ##Ranks the ordered Frequency dictionary 

listF=[]


for key in sorted(RankedF.iterkeys()): ##Appends dictionary to dataframe
    keyF=(RankedF[key])
    df=df.append({"Rank of Frequency": keyF }, ignore_index=True)
    listF.append(keyF)
    
listP=[]
list=[]
RankedC={key: rank for rank, key in enumerate(sorted(CoverageDict, key=CoverageDict.get, reverse=True), 1)} ##Ranks the coverage dictionary

for key in sorted(RankedC.iterkeys()):    
    keyC=(RankedC[key])
    df=df.append({"Rank of Coverage": keyC }, ignore_index=True)
    listP.append(keyC)

list=map(add,listF,listP) ##Allows the product of the two lists to be created.

FreqM = np.array(listF)
CovM = np.array(listP)
M=CovM*FreqM ##Allows the sum of two lists to be created.

df2=pandas.DataFrame(list,columns=['Product of Rank'])
df3=pandas.DataFrame(M, columns=['Sum of Rank'])

df=df.append(df2)
df=df.append(df3)
print dfINT

df1 = df.apply(lambda x: pandas.Series(x.dropna().values)) ##Removes all NaN values
dfINT=dfINT.fillna("")

#dfSym.to_csv("SymTestSY4INTTESTOGeneTEST.csv", sep='\t') ##Writes dataframe to csv file
#dfINT.to_csv("MutationInSY4GeneTESTOGeneTEST.csv", sep='\t')
df1.to_csv("TheirsTestSymIntSY4TESTOGeneTESTFran.csv", sep='\t')

f=open("PointTestSY4INTTESTOGeneFran.txt", 'w') ##Writes new sequence to text file
seqW=str(seql)
f.write(seqW)

#NoCover=set(OnlyInTheirsSY4)-set(NoCov)
#for x in NoCover:
#    dfTheirs=dfTheirs.append({"Mutations at positions of Zero coverage": x, "Reference Base at Zero Coverage ": seq[x]}, ignore_index=True)
#print (test[150])
#dfTheirs=dfTheirs.apply(lambda x: pandas.Series(x.dropna().values))      
#dfTheirs.to_csv("TheirsTestMSY150.csv", sep='\t')

#print Sym  
#print Sym2 


#dfTheirs.to_csv("TheirsTestMSY150.csv", sep='\t')
#dfTheirsLE.to_csv("TheirsTestMSY1LE50.csv", sep='\t')

#df.to_csv("dummyRank.csv", sep='\t')
#df.dropna(inplace=True)
#df2.to_csv("dummyRank.csv", mode='a', header=False, sep='\t')
#with open('dummyRank.csv', 'a') as f:
    #(df2).to_csv(f, header=False)
#df=df.append({"Rank of Coverage": keyF }, ignore_index=True)
#df.to_csv("dummyRank.csv", sep='\t')
#with open('dummyRank.csv', 'a') as f:  # Just use 'w' mode in 3.x
#    w = csv.DictWriter(f, RankedC.keys())
#    w.writeheader()
#    w.writerow(RankedC)
#print RankedC
#with open('dummyRank.csv', mode='a') as file:
 #           file_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
  #          file_writer.writerow(["Rank of Coverage", RankedF])
   #         file_writer.writerow(["Rank of Frequency", RankedC])
#with open('dummyRank.csv', 'wb') as outcsv:
#    writer = csv.DictWriter(outcsv, fieldnames = ["Status", "Reference Position", "Strain Position"])
#    writer.writeheader()
#for x in OnlyInTheirs:
 #   if x in NoCov:
  #      ZeroCoverageMutants=x
   #     dfTheirs=dfTheirs.append({"Mutations at positions of Zero coverage":ZeroCoverageMutants, "Reference Base": seq[x] }, ignore_index=True)
            
#with open("SequenceCovSY1.txt","wb") as handle:
#    pickle.dump(seql, handle)
#    handle.close() 
#infile= open("SequenceCovSY1.txt", 'rb')
#test= pickle.load(infile)
#infile.close()

#for x in NoCov:
    #if x in OnlyInTheirs:
     #   y=OnlyInTheirs.remove(x)
      #  print y

        #print seq[ref_pos]
  
        #remove seq.base[ref_pos]? and replace with StrainBase
        # append coverage and old/new base to meta file
    #if StrainBase==base[ref_pos]:
        #conver base[ref_pos] to upper case
        
    
        
        
    
              
