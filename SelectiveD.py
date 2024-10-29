#########################################################################################################
#########################################################################################################
#
# THIS SCRIPT CALCULATES REGULAR, SYNONYMOUS AND NONSYNONYMOUS D-TEST FROM VCF INPUT
# (WORKS WELL IF CODON POSITIONS PROPERLY MAPPED - TAKE CARE OF IT!)
# 
#########################################################################################################
#########################################################################################################


#   The workflow:
# - Read the whole multi-species VCF file line by line
# - Make all possible quartets form species, listed in VCF, with single outgroup (defined here in script)
#   (in selecting quartets with proper topology this script relies on user's tree)
# - Analyse all quartets for each line (huge time economy on file reading)
# - Calculate D (and F4 - to be impemented later), for syn and nonsyn sites independently


#########################################################################################################
#
# LIB IMPORTS
#
#########################################################################################################

#/data/antwerpen/207/vsc20750/AMPHIPODS/NEW_D_TEST
import os, yaml
import numpy as np
import pandas as pd
import pysam
from matplotlib import pyplot as plt
#from pypopgen3.modules import treetools as tt
import ete3
from ete3 import Tree
import itertools
jn = os.path.join
eu = os.path.expanduser

from collections import OrderedDict
import time


#########################################################################################################
#
# READING VCF INPUT
#
#########################################################################################################

# VCF with list of samples:

#from short tests to long ones:

#vcf_fn = '/data/antwerpen/207/vsc20750/AMPHIPODS/NEW_D_TEST/MERGED_test_quartet_subsample_50rows.vcf' # just string
#vcf_fn = '/data/antwerpen/207/vsc20750/AMPHIPODS/NEW_D_TEST/MERGED_test_quartet_subsample.vcf' # just string
#vcf_fn = '/data/antwerpen/207/vsc20750/AMPHIPODS/NEW_D_TEST/MERGED_test_quartet.vcf' # just string
#vcf_fn = '/data/antwerpen/207/vsc20750/AMPHIPODS/NEW_D_TEST/MERGED_five.vcf' # just string
#vcf_fn = '/data/antwerpen/207/vsc20750/AMPHIPODS/NEW_D_TEST/MERGED_eleven.vcf' # just string

#main vcf:

vcf_fn = '/data/antwerpen/207/vsc20750/AMPHIPODS/NEW_D_TEST/MERGED.vcf' # just string

bcf_in = pysam.VariantFile(vcf_fn)  # auto-detect input format
#bcf_in.header.samples
samples = [s for s in bcf_in.header.samples]  # iterator fetches names from object, puts them into a list
n_samples = len(samples)



#########################################################################################################
#
# READ THE TREE
#
#########################################################################################################

tree = Tree("/data/antwerpen/207/vsc20750/AMPHIPODS/NEW_D_TEST/Gammaridae_tree_scinames_nobrlen_46.newick")
#tree.render("%%inline") 


#########################################################################################################
#
# MAKING ALL POSSIBLE TRIOS WITH GIVEN OUTGROUP AND SPECIFIC TOPOLOGY
#
#########################################################################################################


#List of all species, which are not outgroup:
outgroup = 'Pandorites_podoceroides'
ingroups = [s for s in samples if s != outgroup]
#ingroups

#tree.prune([a,b,c])
#ingroups = [s for s in samples if s != outgroup]

trios = []
for a,b,c in itertools.combinations(ingroups,3):
    ab_LCA=tree.get_common_ancestor(a,b)
    ac_LCA=tree.get_common_ancestor(a,c)
    bc_LCA=tree.get_common_ancestor(b,c)
    
    if ac_LCA==bc_LCA and ab_LCA!=ac_LCA: #second check - for avoiding unresolved trees
        p1, p2, p3 = a, b, c
    elif ab_LCA==bc_LCA and ac_LCA!=ab_LCA: #second check - for avoiding unresolved trees
        p1, p2, p3 = a, c, b
    elif ab_LCA==ac_LCA and bc_LCA!=ab_LCA: #second check - for avoiding unresolved trees
        p1, p2, p3 = b, c, a
    else:
        print("exception - unresolved node in the tree, or other strange problem, that may be serious")
        break
    trios.append((p1, p2, p3))
#trios    
#len(trios)



#########################################################################################################
#########################################################################################################
#
# HASHES
#
#########################################################################################################
#########################################################################################################

#########################################################################################################
#
# MAKING CODON POSITION HASH FOR GIVEN VCF
#
#########################################################################################################

# Final hash CODON_POSITIONS defines codon position for each raw of VCF file
# CODON_POSITIONS key - merged VCF fields "CHROM" and "POS" (with underline in between)
# CODON_POSITIONS value - 1,2,3 if codon is defined, 0 if it's unclear (these positions should not be used)


###############################################
# FIRST PART: we make hash with codon positions
# - position 1, 2 or 3 is identified by division of position in CDS by 3
# - pos 2 must be precieded by pos 1
# - pos 3 must e precieded by pos 2 and 1
# - but! there is a problem - if middle of the gene is missing, codon may consist of only positions 1 and 2
# this case will be fixed in SECOND PART of this script
###############################################

#CODON_POSITIONS = {}
CODON_POSITIONS = OrderedDict()

# I expect that first codon in chromosome may be lost due to imperfect seed structure
# Can be changed later
Minus_1_row_CHROM="seed"
Minus_1_row_CHROM_POS=1
Minus_2_row_CHROM="seed"
Minus_2_row_CHROM_POS=1
    
for i, rec in enumerate(bcf_in.fetch(*())):
    
    #if i>10:
    #    break
    
    
    
    #print(rec)
    
    Current_row_CHROM = rec.chrom
    #print("Current row CHROM: ",Current_row_CHROM)
    Current_row_CHROM_POS = rec.pos
    dict_key=Current_row_CHROM+"_"+str(Current_row_CHROM_POS)
    #print("Current row CHROM POS: ",Current_row_CHROM_POS)
    Current_row_CHROM_POS=int(Current_row_CHROM_POS)
    
    
    #if Current_row_CHROM!="cog1004":
        #break
    
    
    if Current_row_CHROM_POS%3 == 1:
        #print ("FIRST POSITION OF THE CODON")
        # check of the codon contiguity:
        # no checks for the first position
        CODON_POSITIONS[dict_key]=1
        
    if Current_row_CHROM_POS%3 == 2:
        #print ("SECOND POSITION OF THE CODON")
        # check of the codon contiguity:
        # If we handle 2nd position - we should check that
        # - previous row (first position of same codon) of VCF corresponds to previous position in alignment
        # - previous row (first position of same codon) belong to same chromosome
        if Current_row_CHROM_POS-Minus_1_row_CHROM_POS==1 and Current_row_CHROM==Minus_1_row_CHROM:
            #print("GOOD CHROM CONTIGUITY")
            CODON_POSITIONS[dict_key]=2
        else:
            #print("CHROMOSOME CONTIGUITY DISRUPTION!")
            CODON_POSITIONS[dict_key]=0
            #break
            #continue

        
    if Current_row_CHROM_POS%3 == 0:
        # print ("THIRD POSITION OF THE CODON")
        # check of the codon contiguity:
        # If we handle 3rd position - we should check that
        # - previous row (second position of same codon) of VCF corresponds to previous position in alignment
        # - previous row (second position of same codon) belong to same chromosome
        # - The row before previous row (first position of same codon) of VCF corresponds to -2 position in alignment
        # - The row before previous row (first position of same codon) belong to same chromosome
        if (Current_row_CHROM_POS-Minus_1_row_CHROM_POS==1 and Current_row_CHROM==Minus_1_row_CHROM) and (Current_row_CHROM_POS-Minus_2_row_CHROM_POS==2 and Current_row_CHROM==Minus_2_row_CHROM):
            #print("GOOD CHROM CONTIGUITY")
            CODON_POSITIONS[dict_key]=3
        else:
            #print("CHROMOSOME CONTIGUITY DISRUPTION!")
            CODON_POSITIONS[dict_key]=0
   
    
    Minus_2_row_CHROM=Minus_1_row_CHROM
    Minus_2_row_CHROM_POS=Minus_1_row_CHROM_POS
    Minus_1_row_CHROM=Current_row_CHROM
    Minus_1_row_CHROM_POS=Current_row_CHROM_POS

    


#################################################################
# SECOND PART
# Here we go through dataframe, and check that each codon has position 1, 2 and 3
# If one of positions is missing - it substitutes all available positions of broken codon by 0s
#################################################################

Min2_key=None
Min2_val=None
Min1_key=None
Min1_val=None
for key in CODON_POSITIONS:

    #print(key, "->", CODON_POSITIONS[key])
    Current_key=key
    Current_val=CODON_POSITIONS[key]
    
    if Min2_val==1:
        #print("Min2_val==1")
        if Min1_val!=2 or Current_val!=3:
            #print("BAD CODON AROUND POS 1, CHANGE 1 TO 0")
            CODON_POSITIONS[Min2_key]=0
        #else:
            #print("GOOD CODON AROUND POS 1")
    if Min1_val==2:
        #print("Min1_val==1")
        if Min2_val!=1 or Current_val!=3:
            #print("BAD CODON AROUND POS 2, CHANGE 1 TO 0")
            CODON_POSITIONS[Min1_key]=0
        #else:
            #print("GOOD CODON AROUND POS 2")
    if Current_val==3:
        #print("Current_val==3")
        if Min2_val!=1 or Min1_val!=2:
            #print("BAD CODON AROUND POS 3, CHANGE 1 TO 0")
            CODON_POSITIONS[Current_key]=0
        #else:
            #print("GOOD CODON AROUND POS 3")

    Min2_key=Min1_key
    Min2_val=Min1_val
    Min1_key=Current_key
    Min1_val=Current_val


    
    
#########################################################################################################
#  
#  HASHES FOR DETECTION OF NON-DEGENERATE AND 4-FOLD-DEGENERATE CODONS
#
#########################################################################################################


# 1) Hash for detection of codons, which are non-degenerate by 1st position
# POS1_NONDEGEN_dict keys - dinucleiotides made of 2nd and 3rd codon position, lowercase
# POS1_NONDEGEN_dict values - True is codon is non-degenerate, False if not

POS1_NONDEGEN_dict={}

bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
   
    
for a, nuc1 in enumerate (bases):
    for b, nuc2 in enumerate (bases):
        
        #print(nuc1,nuc2)
        
        dinuc_key=nuc1+nuc2
        codon1="t"+dinuc_key
        codon2="c"+dinuc_key
        codon3="a"+dinuc_key
        codon4="g"+dinuc_key

        AA1=codon_table.get(codon1)
        AA2=codon_table.get(codon2)
        AA3=codon_table.get(codon3)
        AA4=codon_table.get(codon4)

        # if one of the codons is stop codon - we don't concider it 4-fold-degenerate (no matter what, we just ignore this codon)
        if (AA1=="*" or AA2=="*" or AA3=="*" or AA4=="*"):
            POS1_NONDEGEN_dict[dinuc_key]=False
        else:    
            # Looking at these comparisons I become worried about correctness of P-test. Take a look at GitHub!
            if AA1!=AA2 and AA2!=AA3 and AA3!=AA4 and AA1!=AA4 and AA1!=AA3 and AA2!=AA4:
                #print ("4-fold degenerate codon")
                POS1_NONDEGEN_dict[dinuc_key]=True
            else:
                POS1_NONDEGEN_dict[dinuc_key]=False   

                
                
# 2) Hash for detection of codons, which are non-degenerate by 2nd position
# POS2_NONDEGEN_dict keys - dinucleiotides made of 1st and 3rd codon position, lowercase
# POS2_NONDEGEN_dict values - True is codon is non-degenerate, False if not

                
POS2_NONDEGEN_dict={}

bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
   
    
for a, nuc1 in enumerate (bases):
    for b, nuc2 in enumerate (bases):
        
        #print(nuc1,nuc2)
        
        dinuc_key=nuc1+nuc2
        codon1=nuc1+"t"+nuc2
        codon2=nuc1+"c"+nuc2
        codon3=nuc1+"a"+nuc2
        codon4=nuc1+"g"+nuc2

        AA1=codon_table.get(codon1)
        AA2=codon_table.get(codon2)
        AA3=codon_table.get(codon3)
        AA4=codon_table.get(codon4)

        # if one of the codons is stop codon - we don't concider it 4-fold-degenerate (no matter what, we just ignore this codon)
        if (AA1=="*" or AA2=="*" or AA3=="*" or AA4=="*"):
            POS2_NONDEGEN_dict[dinuc_key]=False
        else:    
            # Looking at these comparisons I become worried about correctness of P-test. Take a look at GitHub!
            if AA1!=AA2 and AA2!=AA3 and AA3!=AA4 and AA1!=AA4 and AA1!=AA3 and AA2!=AA4:
                #print ("4-fold degenerate codon")
                POS2_NONDEGEN_dict[dinuc_key]=True
            else:
                POS2_NONDEGEN_dict[dinuc_key]=False                 
                
                


# 3) Hash for detection of codons, which are 4-fold-degenerate by 3rd position
# POS3_DEGEN_dict keys - dinucleiotides made of 1st and 2nd codon position, lowercase
# POS3_DEGEN_dict values - True is codon is 4-fold-degenerate, False if not

POS3_DEGEN_dict={}

bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))
   
    
for a, nuc1 in enumerate (bases):
    for b, nuc2 in enumerate (bases):
        
        #print(nuc1,nuc2)
        
        dinuc_key=nuc1+nuc2
        codon1=dinuc_key+"t"
        codon2=dinuc_key+"c"
        codon3=dinuc_key+"a"
        codon4=dinuc_key+"g"

        AA1=codon_table.get(codon1)
        AA2=codon_table.get(codon2)
        AA3=codon_table.get(codon3)
        AA4=codon_table.get(codon4)

        # if one of the codons is stop codon - we don't concider it 4-fold-degenerate (no matter what, we just ignore this codon)
        if (AA1=="*" or AA2=="*" or AA3=="*" or AA4=="*"):
            POS3_DEGEN_dict[dinuc_key]=False
        else:    
            # Looking at these comparisons I become worried about correctness of P-test. Take a look at GitHub!
            if AA1==AA2==AA3==AA4:
                #print ("4-fold degenerate codon")
                POS3_DEGEN_dict[dinuc_key]=True
            else:
                POS3_DEGEN_dict[dinuc_key]=False
                
                
#########################################################################################################
#########################################################################################################
#
# FUNCTIONS
#
#########################################################################################################
#########################################################################################################

#########################################################################################################
#
# FUNCTIONS FOR HYBRIDIZATION TESTS ON SINGLE LOCUS (D and F4)
#
#########################################################################################################                

# D test
# for input takes DIPLOID genotypes of p1, p2, p3 and HAPLOID genotype of outgroup
#
# outgroup genotype we pass as one letter (diploid converted to "haploid"),
# because homozygotes can be discribed by one letter, while
# heterozygous outgroups can't be used for calculation of D and are skipped in main script


def D_test(p1a, p2a, p3a, a):
    
    # check that all genotypes are biallelic
    # except the outgroup - because outgroup is haploid in my case (and it's often haploid, in general,
    # because it's made from distant sequence (often not from a VCF))

    allele_str = ''.join([str(i) for i in (a,) + p1a + p2a + p3a])
    alleles = [int(s) for s in set(allele_str)]
    #print(alleles)

    #print(allele_str) # string, which merged all alleles
    #print(alleles)  # a list of alleles (lenght should be =2, if site is biallelic)
    
    # skip more than biallelic and monomorphic genotypes (only for binary genotypes we can count abba and baba)
    if len(alleles) != 2:
        return(0,0)
    else:

        ##############################################################
        #
        # HERE WE RECODE BINARY SNP INTO 0 and 1
        # There are 4 options: (0, else)(1, else)(else,else)()
        #
        ##############################################################

        # (0, else) case

        if 0 in alleles and 1 not in alleles:
            other_allele = str([s for s in alleles if s!=0][0])  # find other allele
            allele_str1 = allele_str.replace(other_allele, '1') # substitute it by 1
            #cool, now we have string of allels, which only contains 0s and 1s

        # (1, else) case

        elif 0 not in alleles and 1 in alleles:
            other_allele = str([s for s in alleles if s!=1][0]) # find other allele
            allele_str1 = allele_str.replace(other_allele, '0')  # substitute it by 0
            #cool, now we have string of allels, which only contains 0s and 1s

        # (else, else) case

        elif 0 not in alleles and 1 not in alleles:

            #alleles=[2,3]
            #print(allele_str)
            allele_str1=allele_str
            #print(alleles)
            for new_allele, a in zip(['0','1'], alleles): # no idea how it works!
                a=str(a)
                new_allele=str(new_allele)
                allele_str1=allele_str1.replace(a, new_allele)
            #print(allele_str1)
            #cool, now we have string of allels, which only contains 0s and 1s (hopefully)

        # (0, 1) case
        else:
            allele_str1 = allele_str
            #cool, now we have string of allels, which only contains 0s and 1s

        # GOOD. NOW WE HANDLE STRING OF ALLELES, RECODED TO 0s and 1s, to get freqs for F4 test

        # Calculating mean frequencies of alleles.
        # Doesn't make much sense for one sample per species,
        # But we'll make it implemented with a population sample (later, TODO!)

        a, p11, p12, p21, p22, p31, p32 = allele_str1
        p1b = np.mean([int(p11), int(p12)]) # mean frequency: 1 or 0 if homo, 0.5 if hetero
        p2b = np.mean([int(p21), int(p22)]) # mean frequency: 1 or 0 if homo, 0.5 if hetero
        p3b = np.mean([int(p31), int(p32)]) # mean frequency: 1 or 0 if homo, 0.5 if hetero  
        a = int(a)



        # HERE WE COUNT ABBA-BABA score for the quartet (and F4 will be done here, later):  

        # p1b, p2b, p3b can be 0(homo), 0.5(hetero), 1(other homo, not necessarilly alt)

        # ABBA can be formed by either 0110, or 1001
        # so it's either 0*0*

        # equals 1 if pattern is 1001 or 0110
        # equals 0.5 if pattern is (0.5,0,0,1) or (0,0.5,1,0) or (0,1,1,0.5) or (1,0,0.5,1) or similar
        # equals 0.25 if pattern is (0.5,0,0,0.5) or (0,0.5,0.5,0) or (0.5,1,1,0.5) or (1,0.5,0.5,1)
        abba = p1b*(1-p2b)*(1-p3b)*a + (1-p1b)*p2b*p3b*(1-a)

        #works same as above for patterns 1010, 0101 and it's variations with heterozygous allele 0.5 
        baba = (1-p1b)*p2b*(1-p3b)*a + p1b*(1-p2b)*p3b*(1-a)

        #print(abba)
        #print(baba)
        return(abba,baba)
                    

  
        
        
        
#########################################################################################################
#########################################################################################################
#
# MAIN SCRIPT
#
#########################################################################################################                
#########################################################################################################
        
       
        
#%%prun -s cumulative


start_time = time.time()

# Let's Initialize arrays of len(trios) length:
# needed to store info for recent codon:
import pdb
#pdb.set_trace()

# POS1
POS1_trio_presense=[False] * len(trios)

POS1_conserv=[False] * len(trios)
POS1_consensus=[None] * len(trios)

POS1_D_arr_abba=[0] * len(trios)
POS1_D_arr_baba=[0] * len(trios)


# POS2

POS2_trio_presense=[False] * len(trios)

POS2_conserv=[False] * len(trios)
POS2_consensus=[None] * len(trios)

POS2_D_arr_abba=[0] * len(trios)
POS2_D_arr_baba=[0] * len(trios)

# POS3

POS3_trio_presense=[False] * len(trios)

POS3_conserv=[False] * len(trios)
POS3_consensus=[None] * len(trios)

POS3_D_arr_abba=[0] * len(trios)
POS3_D_arr_baba=[0] * len(trios)


# General vars (here we'll sum values for each trio along the whole VCF):

abba_raw = [0 for k in range(len(trios))]
baba_raw = [0 for k in range(len(trios))]
abba_s = [0 for k in range(len(trios))]
baba_s = [0 for k in range(len(trios))]
abba_n = [0 for k in range(len(trios))]
baba_n = [0 for k in range(len(trios))]



# marker of conservativeness:
#Conservative_POS1=False
#Conservative_POS2=False
#Conservative_POS3=False

# arrays of 4 species nucleotides for particular position handling
#CURRENT_CODON_POS1=[]
#CURRENT_CODON_POS2=[]
#CURRENT_CODON_POS3=[]



###############################################################
# HERE WE READ VCF LINE BY LINE
###############################################################


i=0

# i - counter, rec - informative (genotype) strings of VCF 
for i, rec in enumerate(bcf_in.fetch(*())):
    
    if i%1000000 == 0:
        print(i," rows handled")
        
    #if i<214176:
    #    continue
    #if i>214178:
    #    break
    
#    if i<8470 or i>8479:
#        continue
    #print("\nNEW ROW:")
    #print(rec)
    #print("i at the beggining: ",i)
    
    # Row info:
    
    Current_row_CHROM = rec.chrom
    Current_row_CHROM_POS = rec.pos
    dict_key=Current_row_CHROM+"_"+str(Current_row_CHROM_POS)
    ACTIVE_CODON_POSITION=CODON_POSITIONS.get(dict_key,"Error: no such chromosome or position in this VCF")
    
    #print("codon position: ", ACTIVE_CODON_POSITION)
    
    # Alleles:
    
    REF = rec.ref
    ALT = rec.alts
    ALL_ALLELES=rec.alleles
    #print("alt1: ",rec.alts[0])          # first of alt alleles!
    #print("alt2: ",rec.alts[1])          # second of alt alleles! (RETURNS ERROR IF OUT OF RANGE)
    #print("all_alleles: ",all_alleles)          #list of all alleles?

    
    # Outgroup
    outgroup_gt = rec.samples[outgroup]['GT']
    
    
    
    #if ACTIVE_CODON_POSITION==0:
        #print("THE CODON FRAME IS TEMPORARILY SHIFTED, IGNORE THIS POSITION")
        # This option is most painful - think about it later, may cause many errors!
    
    
    
    if ACTIVE_CODON_POSITION==1:
        
        #print("ACTIVE_CODON_POSITION=1")
        if None in outgroup_gt:
            POS1_trio_presense=[False] * len(trios)
            #if j>1:
            #    print("POS1 SKIPPED ROW - OUTGROUP NONE")
            continue # go to the next row

            
        # For each row we check all possible trios
        # One trio per cycle
        for j,(p1, p2, p3) in enumerate(trios):     # p1, p2, p3 - species names

            # p1a, p2a, p3a - genotypes, could be [0,0], [0,1], [2,2], etc
            p1a = rec.samples[p1]['GT']
            p2a = rec.samples[p2]['GT']
            p3a = rec.samples[p3]['GT']
            
            # if alleles contain a missing genotype - we skip this trio
            if None in p1a or None in p2a or None in p3a:
                POS1_trio_presense[j]=False
                #if j>1: # here it works very well, many rows with non empty 4 species (which still have False  POS3_trio_presense later)
                #    print("POS1, absent")
                continue # go to the next trio
            
            # here we finally handle case, where all 4 genotypes are present! let's mark it in hash:
            POS1_trio_presense[j]=True
            
            #if j>1: # here it works very well, many rows with non empty 4 species (which still have False  POS3_trio_presense later)
            #    print("POS1, present")
            #    print(p1, p2, p3)
            #    print(p1a,p2a,p3a,outgroup_gt)
            #    print(POS1_trio_presense[j])
            #    print("i: ",i)
            #    print("j: ",j)
           
            
            # Nucleotides, now we need nucleotides!
                
            # p1
            p11=p1a[0]
            p11_nuc=ALL_ALLELES[p11]
            p12=p1a[1]
            p12_nuc=ALL_ALLELES[p12]
            # p2
            p21=p2a[0]
            p21_nuc=ALL_ALLELES[p21]
            p22=p2a[1]
            p22_nuc=ALL_ALLELES[p22]
            # p3
            p31=p3a[0]
            p31_nuc=ALL_ALLELES[p31]
            p32=p3a[1]
            p32_nuc=ALL_ALLELES[p32]
            # outgroup_gt
            outgroup_gt1=outgroup_gt[0]
            outgroup_gt1_nuc=ALL_ALLELES[outgroup_gt1]
            outgroup_gt2=outgroup_gt[1]
            outgroup_gt2_nuc=ALL_ALLELES[outgroup_gt2]
            
            CURRENT_CODON_POS1=[p11_nuc,p12_nuc,p21_nuc,p22_nuc,p31_nuc,p32_nuc,outgroup_gt1_nuc,outgroup_gt2_nuc]
            
   
            
            if len(set(CURRENT_CODON_POS1)) == 1 and CURRENT_CODON_POS1[0]!='.':
                    POS1_conserv[j]=True
                    POS1_consensus[j]=CURRENT_CODON_POS1[0]
                    
                 
            
            
            # ABBA-BABA PART:
            
            # if outgroup is absent or is not homozygous, hybrid tests (D and F4) are not possible,
            if outgroup_gt[0] == outgroup_gt[1]:
                
                

                #set outgroup allele
                a = outgroup_gt[0]

                (abba,baba)=D_test(p1a,p2a,p3a,a)
                
                abba_raw[j] += abba
                baba_raw[j] += baba
                
                POS1_D_arr_abba[j]=abba
                POS1_D_arr_baba[j]=baba
                

                    
                    
  
    if ACTIVE_CODON_POSITION==2:
        
        #print("ACTIVE_CODON_POSITION=2")
        
        if None in outgroup_gt:
            POS2_trio_presense=[False] * len(trios)
            #if j>1:
            #    print("POS2 SKIPPED ROW - OUTGROUP NONE")
            continue # go to the next row
            
            
        # For each row we check all possible trios
        # One trio per cycle
        for j,(p1, p2, p3) in enumerate(trios):     # p1, p2, p3 - species names

            
            # p1a, p2a, p3a - genotypes, could be [0,0], [0,1], [2,2], etc
            p1a = rec.samples[p1]['GT']
            p2a = rec.samples[p2]['GT']
            p3a = rec.samples[p3]['GT']
            
            # if alleles contain a missing genotype - we skip this trio
            if None in p1a or None in p2a or None in p3a:
                POS2_trio_presense[j]=False
                #if j>1: # here it works very well, many rows with non empty 4 species (which still have False  POS3_trio_presense later)
                #    print("POS2, absent")
                continue # go to the next trio
                
                
        
            
            # here we finally handle case, where all 4 genotypes are present! let's mark it in hash:
            POS2_trio_presense[j]=True
            
            #if j>1: # here it works very well, many rows with non empty 4 species (which still have False  POS3_trio_presense later)
            #    print("POS2, present")
            #    print(p1, p2, p3)
            #    print(p1a,p2a,p3a,outgroup_gt)
            #    print(POS2_trio_presense[j])
            #    print("i: ",i)
            #    print("j: ",j)
            
            # Nucleotides, now we need nucleotides!
                
            # p1
            p11=p1a[0]
            p11_nuc=ALL_ALLELES[p11]
            p12=p1a[1]
            p12_nuc=ALL_ALLELES[p12]
            # p2
            p21=p2a[0]
            p21_nuc=ALL_ALLELES[p21]
            p22=p2a[1]
            p22_nuc=ALL_ALLELES[p22]
            # p3
            p31=p3a[0]
            p31_nuc=ALL_ALLELES[p31]
            p32=p3a[1]
            p32_nuc=ALL_ALLELES[p32]
            # outgroup_gt
            outgroup_gt1=outgroup_gt[0]
            outgroup_gt1_nuc=ALL_ALLELES[outgroup_gt1]
            outgroup_gt2=outgroup_gt[1]
            outgroup_gt2_nuc=ALL_ALLELES[outgroup_gt2]
            
            CURRENT_CODON_POS2=[p11_nuc,p12_nuc,p21_nuc,p22_nuc,p31_nuc,p32_nuc,outgroup_gt1_nuc,outgroup_gt2_nuc]
            
            if len(set(CURRENT_CODON_POS2)) == 1 and CURRENT_CODON_POS2[0]!='.':
                    POS2_conserv[j]=True
                    POS2_consensus[j]=CURRENT_CODON_POS2[0]
            
            
            # ABBA-BABA PART:
            
            # if outgroup is absent or is not homozygous, hybrid tests (D and F4) are not possible,
            if outgroup_gt[0] == outgroup_gt[1]:

                #set outgroup allele
                a = outgroup_gt[0]

                (abba,baba)=D_test(p1a,p2a,p3a,a)

                abba_raw[j] += abba
                baba_raw[j] += baba
                
                POS2_D_arr_abba[j]=abba
                POS2_D_arr_baba[j]=baba
  
    
    if ACTIVE_CODON_POSITION==3:
        
        #print("ACTIVE_CODON_POSITION=3")
        
        if None in outgroup_gt:
            POS3_trio_presense=[False] * len(trios)
            #print("NONE")
            #if j>1:
            #    print("POS3 SKIPPED ROW - OUTGROUP NONE")
            continue # go to the next row
            
        # For each row we check all possible trios
        # One trio per cycle
        for j,(p1, p2, p3) in enumerate(trios):     # p1, p2, p3 - species names
        

            
            # p1a, p2a, p3a - genotypes, could be [0,0], [0,1], [2,2], etc
            p1a = rec.samples[p1]['GT']
            p2a = rec.samples[p2]['GT']
            p3a = rec.samples[p3]['GT']
            

                
            
            # if alleles contain a missing genotype - we skip this trio
            if None in p1a or None in p2a or None in p3a:
                POS3_trio_presense[j]=False
                #if j>1: # here it works very well, many rows with non empty 4 species (which still have False  POS3_trio_presense later)
                #    print("POS3, absent")
                continue # go to the next trio
                

            
            # here we finally handle case, where all 4 genotypes are present! let's mark it in hash:
            POS3_trio_presense[j]=True
           
            #if j>1: # here it works very well, many rows with non empty 4 species (which still have False  POS3_trio_presense later)
            #    print("POS3, present")
            #    print(p1, p2, p3)
            #    print(p1a,p2a,p3a,outgroup_gt)
            #    print(POS3_trio_presense[j])
            #    print("i: ",i)
            #    print("j: ",j)
                #break
            # Nucleotides, now we need nucleotides!
                
            # p1
            p11=p1a[0]
            p11_nuc=ALL_ALLELES[p11]
            p12=p1a[1]
            p12_nuc=ALL_ALLELES[p12]
            # p2
            p21=p2a[0]
            p21_nuc=ALL_ALLELES[p21]
            p22=p2a[1]
            p22_nuc=ALL_ALLELES[p22]
            # p3
            p31=p3a[0]
            p31_nuc=ALL_ALLELES[p31]
            p32=p3a[1]
            p32_nuc=ALL_ALLELES[p32]
            # outgroup_gt
            outgroup_gt1=outgroup_gt[0]
            outgroup_gt1_nuc=ALL_ALLELES[outgroup_gt1]
            outgroup_gt2=outgroup_gt[1]
            outgroup_gt2_nuc=ALL_ALLELES[outgroup_gt2]
            
            CURRENT_CODON_POS3=[p11_nuc,p12_nuc,p21_nuc,p22_nuc,p31_nuc,p32_nuc,outgroup_gt1_nuc,outgroup_gt2_nuc]
            
            if len(set(CURRENT_CODON_POS3)) == 1 and CURRENT_CODON_POS3[0]!='.':
                    POS3_conserv[j]=True
                    POS3_consensus[j]=CURRENT_CODON_POS3[0]
            
            
            # ABBA-BABA PART:
            
            # if outgroup is absent or is not homozygous, hybrid tests (D and F4) are not possible,
            if outgroup_gt[0] == outgroup_gt[1]:

                #set outgroup allele
                a = outgroup_gt[0]

                (abba,baba)=D_test(p1a,p2a,p3a,a)
                
                abba_raw[j] += abba
                baba_raw[j] += baba
                
                POS3_D_arr_abba[j]=abba
                POS3_D_arr_baba[j]=baba
        
        #print("ARRS:")
        #print(POS1_trio_presense)
        #print(POS2_trio_presense)
        #print(POS3_trio_presense)
   
        #############################################################################################################
        # Codon handling part - happens after each third position of codon
        #############################################################################################################    
    
        for j,(p1, p2, p3) in enumerate(trios):     # p1, p2, p3 - species names
        

            
            #if j>1:
            #    print("CODON HANDLING")   #HERE works, all presense false
            #    #print(POS1_D_arr_abba[j])
            #    
            #    print("ARRS2:")
            #    print(POS1_trio_presense)
            #    print(POS2_trio_presense)
            #    print(POS3_trio_presense)
            #    
            #    print(p1, p2, p3)
            #    print(POS1_trio_presense[j])
            #    print(POS2_trio_presense[j])
            #    print(POS3_trio_presense[j])
            #    print("i: ",i)
            #    print("j: ",j)
            #    print("\n")

            
            # check that for this trio (tree specs and outgroup) all genotypes for codon are not None
            if POS1_trio_presense[j]==True and POS2_trio_presense[j]==True and POS3_trio_presense[j]==True:

                # aha, cool, now we can check if codons are non-degenerate of 4-fold-degenerate

                
                #if j>1:
                #    print("CODON HANDLING INSIDE THE CONS CHECK")   #NOTHING HERE!!
                    #print(POS1_D_arr_abba[j])
            
            
            
                ##########################################################################    

                # Here we handle 3 possible conservative codons
                # 1) Positions 2 and 3 are conservative + 1st position is NON-DEGENERATE
                # 2) Positions 1 and 3 are conservative + 2nd position is NON-DEGENERATE
                # 3) Positions 1 and 2 are conservative + 3rd position is 4 fold DEGENERATE

                ##########################################################################
                # here we get rid of conservative codons:
                if POS1_conserv[j] and POS2_conserv[j] and POS3_conserv[j]:
                    
                    #print(POS1_conserv)
                    #print(POS2_conserv)
                    #print(POS3_conserv)
                    
                    
                    # before skipping, we need to return default values of the variables
                    # (these vars were collected for this codon, but will not be used, because it's conservative)

                    # POS1
                    #POS1_trio_presense=[False] * len(trios)

                    #POS1_conserv=[False] * len(trios)
                    #POS1_consensus=[None] * len(trios)

                    #POS1_D_arr_abba=[0] * len(trios)
                    #POS1_D_arr_baba=[0] * len(trios)


                    # POS2

                    #POS2_trio_presense=[False] * len(trios)

                    #POS2_conserv=[False] * len(trios)
                    #POS2_consensus=[None] * len(trios)

                    #POS2_D_arr_abba=[0] * len(trios)
                    #POS2_D_arr_baba=[0] * len(trios)

                    # POS3

                    #POS3_trio_presense=[False] * len(trios)

                    #OS3_conserv=[False] * len(trios)
                    #OS3_consensus=[None] * len(trios)

                    #OS3_D_arr_abba=[0] * len(trios)
                    #OS3_D_arr_baba=[0] * len(trios)


                    # marker of conservativeness:
                    #Conservative_POS1=False
                    #Conservative_POS2=False
                    #Conservative_POS3=False

                    # arrays of 4 species nucleotides for particular position handling
                    #CURRENT_CODON_POS1=[]
                    #CURRENT_CODON_POS2=[]
                    #CURRENT_CODON_POS3=[]
                   
                    continue # we don't need to re-init vars here - it will be done after the cycle
                    
                
                # (1) codones, which are 4-fold-nondegenerate by 1st position:

                if POS2_conserv[j] and POS3_conserv[j]:

                    #print ("this codon is conservative by position 2 and 3, let's check if it's 4-fold-nondegenerate by 1st position")

                    # First we have to pick random nucleotide from conservative positions 2 and 3
                    # Than we'll be able to form the codon and test if it's non-degenerate

                    # POS2 NUCLEOTIDE:
                    NUC_POS2=POS2_consensus[j]  # random conservative position

                    # POS3 NUCLEOTIDE:
                    NUC_POS3=POS3_consensus[j]  # random conservative position

                    #Let's make a dinucleotide from positions 2 and 3:

                    #dinuc="AC"
                    #dinuc=NUC_POS2+NUC_POS3

                    # CHECKING IF THE CODON IS NON-DEGENERATE BY 1st POSITION:

                    #codon1=("A"+dinuc).lower()
                    #codon2=("C"+dinuc).lower()
                    #codon3=("T"+dinuc).lower()
                    #codon4=("G"+dinuc).lower()
                    #print(codon1,codon2,codon3,codon4)
                    #%debug

                    #decision_pos1=is_nondegenerate(codon1,codon2,codon3,codon4)
                    #print(decision_pos1)
                    
                    dinuc=(NUC_POS2+NUC_POS3).lower()
                    decision_pos1=POS1_NONDEGEN_dict[dinuc]
                    
                    # TEMP!!!
                    if(decision_pos1==True):
                        #print ("this codon is non-degenerate by 1st position")
                        
                        #########################################################
                        # Non-degenerate codon found - time to collect proper ABBA BABA patterns
                        #########################################################
                        
                        # TEMP!!
                        abba_n[j]+=POS1_D_arr_abba[j]
                        baba_n[j]+=POS1_D_arr_baba[j]
                        


                        ####################################################################################################    

                # (2) codones, which are 4-fold-nondegenerate by 2nd position:

                if POS1_conserv[j] and POS3_conserv[j]:

                    #print ("this codon is 4-fold-nondegenerate by 2nd position")

                    # First we have to pick random nucleotide from conservative positions 1 and 3
                    # Than we'll be able to form the codon and test if it's non-degenerate

                    # POS1 NUCLEOTIDE:
                    NUC_POS1=POS1_consensus[j]  # random conservative position

                    # POS3 NUCLEOTIDE:
                    NUC_POS3=POS3_consensus[j]  # random conservative position

                    # CHECKING IF THE CODON IS NON-DEGENERATE BY 1st POSITION:

                    #codon1=(NUC_POS1+"A"+NUC_POS3).lower()
                    #codon2=(NUC_POS1+"C"+NUC_POS3).lower()
                    #codon3=(NUC_POS1+"T"+NUC_POS3).lower()
                    #codon4=(NUC_POS1+"G"+NUC_POS3).lower()
                    #print(codon1,codon2,codon3,codon4)
                    #decision_pos2=is_nondegenerate(codon1,codon2,codon3,codon4)
                    #print(decision_pos2)
                    
                    dinuc=(NUC_POS1+NUC_POS3).lower()
                    decision_pos2=POS2_NONDEGEN_dict[dinuc]

                    # TEMP!!!
                    if(decision_pos2==True):
                        #print ("this codon is non-degenerate by 2nd position")
                        
                        #########################################################
                        # Non-degenerate codon found - time to collect proper ABBA BABA patterns
                        #########################################################
                        
                        #TEMP!!!
                        abba_n[j]+=POS2_D_arr_abba[j]
                        baba_n[j]+=POS2_D_arr_baba[j]

                ############################################################################################  


                # (3) codones, which are 4-fold-degenerate by 3rd position:

                if POS1_conserv[j] and POS2_conserv[j]:

                    #if j>1:
                    #    print(abba_n)
                    #    print("HERE")
                    
                    #print ("this codon is 4-fold-degenerate by 3rd position")

                    # First we have to pick random nucleotide from conservative positions 1 and 2
                    # Than we'll be able to form the codon and test if it's degenerate

                    # POS1 NUCLEOTIDE:
                    NUC_POS1=POS1_consensus[j]  # random conservative position
                    
                    # POS2 NUCLEOTIDE:
                    NUC_POS2=POS2_consensus[j]  # random conservative position


                    #Let's make dinucleotide from 2 and 3 position:
                    #dinuc=(NUC_POS1+NUC_POS2).lower

                    # CHECKING IF THE CODON IS NON-DEGENERATE BY 1st POSITION:

                    #codon1=(dinuc+"A").lower()
                    #codon2=(dinuc+"C").lower()
                    #codon3=(dinuc+"T").lower()
                    #codon4=(dinuc+"G").lower()
                    #print(codon1,codon2,codon3,codon4)
                    #dinuc=dinuc.lower()
                    #decision_pos3=is_4fold_degenerate(codon1,codon2,codon3,codon4)
                    
                    dinuc=(NUC_POS1+NUC_POS2).lower()
                    decision_pos3=POS3_DEGEN_dict[dinuc]
                    
                    

                    #if not POS3_conserv[j]:
                      #codon_pos_string="".join(CURRENT_CODON_POS1) + " " + "".join(CURRENT_CODON_POS2) + " " + "".join(CURRENT_CODON_POS3)
                      #print (f"4fd3rd: {i} {dict_key} {decision_pos3} {POS3_D_arr_abba} {POS3_D_arr_baba} {dinuc} {codon_pos_string}")

                    
                    if(decision_pos3==True):
                        #print ("this codon is degenerate by 3rd position")
                        #########################################################
                        # 4-fold-degenerate codon found - time to collect proper ABBA BABA patterns
                        #########################################################
                        abba_s[j]+=POS3_D_arr_abba[j]
                        baba_s[j]+=POS3_D_arr_baba[j]
                        
                        



        ################################################################################################## 



        # returning default values of the variables

        # POS1
        POS1_trio_presense=[False] * len(trios)

        POS1_conserv=[False] * len(trios)
        POS1_consensus=[None] * len(trios)

        POS1_D_arr_abba=[0] * len(trios)
        POS1_D_arr_baba=[0] * len(trios)


        # POS2

        POS2_trio_presense=[False] * len(trios)

        POS2_conserv=[False] * len(trios)
        POS2_consensus=[None] * len(trios)

        POS2_D_arr_abba=[0] * len(trios)
        POS2_D_arr_baba=[0] * len(trios)

        # POS3

        POS3_trio_presense=[False] * len(trios)

        POS3_conserv=[False] * len(trios)
        POS3_consensus=[None] * len(trios)

        POS3_D_arr_abba=[0] * len(trios)
        POS3_D_arr_baba=[0] * len(trios)


        # marker of conservativeness:
        #Conservative_POS1=False
        #Conservative_POS2=False
        #Conservative_POS3=False

        # arrays of 4 species nucleotides for particular position handling
        #CURRENT_CODON_POS1=[]
        #CURRENT_CODON_POS2=[]
        #CURRENT_CODON_POS3=[]

    
    #print(i)
    #if i==7775:
    #    break
    #if i%100 == 0:
    #    print(i," rows handled")
            
            
            
    # this output block will be moved outside this cycle, when debug finished:
    #if i>5000000:
    #if i>40:
        #print(abba_raw)
        #print(baba_raw)
D_out = open('/data/antwerpen/207/vsc20750/AMPHIPODS/NEW_D_TEST/D_raw_smaller_ns_14000.txt', 'w')

header="P1,P2,P3,ABBA,BABA,D,ABBA_n,BABA_n,D_n,ABBA_s,BABA_s,D_s"
print(header)
D_out.write(header)
for x in range(len(trios)):
    # species names:
    mystring = ' '.join(map(str, trios[x]))

    print(abba_n[x])

    # regular D
    if abba_raw[x]==0 or baba_raw[x]==0:
        D="NA"                
    else:
        D=str((abba_raw[x]-baba_raw[x])/(abba_raw[x]+baba_raw[x]))

    # nonsyn D
    if abba_n[x]==0 or baba_n[x]==0:
        D_n="NA"                
    else:
        D_n=str((abba_n[x]-baba_n[x])/(abba_n[x]+baba_n[x]))

    # syn D
    if abba_s[x]==0 or baba_s[x]==0:
        D_s="NA"                
    else:
        D_s=str((abba_s[x]-baba_s[x])/(abba_s[x]+baba_s[x]))

    ABBA_raw=str(abba_raw[x])
    BABA_raw=str(baba_raw[x])
    ABBA_n=str(abba_n[x])
    BABA_n=str(baba_n[x])
    ABBA_s=str(abba_s[x])
    BABA_s=str(baba_s[x])


    #out_line=mystring+" "+str(abba_raw[x])+" "+str(baba_raw[x])+" "+D+"\n"
    out_line=mystring+" "+ABBA_raw+" "+BABA_raw+" "+D+" "+ABBA_n+" "+BABA_n+" "+D_n+" "+ABBA_s+" "+BABA_s+" "+D_s+"\n"
    #print(out_line)

    # OUTPUT OF RAW ABBA AND BABA
    D_out.write(out_line)

    #print(p1,p2,p3)
    #exit

    # Gammarus_lacustris Eulimnogammarus_cruentus Asprogammarus_rhodophthalmus
    #if trios[x][0]=="Gammarus_lacustris" and trios[x][1]=="Eulimnogammarus_cruentus" and trios[x][2]=="Asprogammarus_rhodophthalmus":
    #        print("first_trio: ",out_line)
    #if trios[x][0]=="Eulimnogammarus_cruentus" and trios[x][1]=="Gammarus_lacustris" and trios[x][2]=="Asprogammarus_rhodophthalmus":
    #        print("second_trio: ",out_line)


    print(x," ",out_line)




    #if x>10:
    #    break
D_out.close()
#bcf_in.close()
#this break is a break for i>5000
        
  

print("--- %s seconds ---" % (time.time() - start_time))            
            
            
            
            
                    