import math
import argparse
import pandas as pd
import numpy as np
import subprocess
import os
import collections
import re

parser = argparse.ArgumentParser()
parser.add_argument("--min-coverage", dest="mc", help="minimum coverage threshold")
parser.add_argument("--min-count", dest="m", help="minimum count threshold")
parser.add_argument("--input", dest="input", help="input file")
parser.add_argument("--output", dest="o", help="output file(s)")
args = parser.parse_args()

# Synopsis
#python3.9 extract_wolbachia3.py --min-coverage 20 --input syncountswCer12.3allele.txt --min-count 20 --output output-file
# awk 'BEGIN{ FS=OFS="\t" }{print $1, $2, $3, $10, $15, $16}' syncounts2.txt > scounts.txt

# Quality strings:
# ________________

# r: target allele similar to reference
# a: target allele different from reference
# c: coverage too low
# p: reference is polymorphic
# /: alternative AF outside boundaries of BCI for the particular individual
# z: more than two alleles 
# m: target allele coverage too low
# M: positions is monomorphic - target allele same as reference  



#################### confidence interval ###########################

def binomial_ci(x, N, alpha=0.1):
    #x is number of cases, N is sample size, alpha is set at 95% confidence intervals

    from scipy import stats
    ### lower limit
    if x==0:
        c1 = 0
    else:
        c1 = stats.beta.interval(1-alpha, x,N-x+1)[0]
    ### upper limit
    if (x==N and N!=0):
        c2=1
    elif N==0:
        c2=0
    else:
        c2 = stats.beta.interval(1-alpha, x+1,N-x)[1]
    return c1, c2
    #result = {'Proportion':x, 'Lower CI': c1, 'Upper CI': c2}
    #return result


def multinomial_ci(x, alpha=0.1):

    from statsmodels.stats.proportion import multinomial_proportions_confint

    #X=x[0]
    m=multinomial_proportions_confint(x, alpha, method='sison-glaz')
    return m


def CI_ind():
    for pop,alleles in sorted(popalleles.items()): #for population number and respective alleles with counts in each gen positions for each population
        #print(alleles.items())
        if sum(alleles.values())<int(args.mc): #if sum of alleles for each pop in each position is lower than --min-coverage threshold
            #print(sum(alleles.values()))
            hapallele.append("N") #append "N" in hapallele columns of the consensus file
            qual.append("c") #append "c" for quality in the quality column of the consensus file
            freqal.append("NA") #append "NA" for list of allele frequncies of haplotype alleles in the *.af output file

        else:
            ## test if alleles found in individual is fullfilling the minimum count criteria else delete the allele
            for k,val in list(alleles.items()): #for alleles and the respective count in each genoomic pos for each population
                #print(val)
                if val<=2: #if the count is lower or equal 2 
                    #print(alleles[k])
                    del alleles[k] #delete the alleles
            if len(alleles)==2: #if length of alleles equals 2
                #print(alleles.keys())
                if refal not in alleles.keys(): #and if reference allele is not in alleles
                    hapallele.append("N")
                    qual.append("z")
                    freqal.append("NA")
                    continue
                total=sum(alleles.values())
                #print(min(alleles.values()))
                #print(max(alleles.values()))
                #print(min(alleles.values())/float(sum(alleles.values())))
                #print(binomial_ci(total/2,total/2)[0])
                #perc=max(alleles.values())/min(alleles.values())/100
                perc1=(list(alleles.values())[0]/total)
                perc2=(list(alleles.values())[1]/total)
                #print(perc)
                #print(binomial_ci(total*perc,total*(1-perc)))
                ## test if CI fits the data:
                if binomial_ci(total*perc1, total*perc2)[0]<=(min(alleles.values())/float(sum(alleles.values()))):
                    #print(alleles.keys())
                    #print(min(alleles.values()))
                    for k in alleles.keys():
                        if k!=refal[0]:
                            if k!=refal2[0]:
                             ## test if alternative allele occurs more than args.m in cummulative sample
                                if fullalleles[k]>=int(args.m):
                                    hapallele3.append(k)
                                    qual3.append("a3")
                                    freqal3.append(alleles[k]/float(sum(alleles.values())))
                                else:
                                    hapallele3.append("N")
                                    qual3.append("m3")
                                    freqal3.append("NA")
                            else:
                                if fullalleles[k]>=int(args.m):
                                    hapallele.append(k)
                                    #hapallele3.append(k)
                                    qual.append("a")
                                    #qual3.append("a")
                                    freqal.append(alleles[k]/float(sum(alleles.values())))
                                    #freqal3.append(alleles[k]/float(sum(alleles.values())))
                                else:
                                    hapallele.append("N")
                                    qual.append("m")
                                    freqal.append("NA")
                else:
                    hapallele.append("N")
                    qual.append("/")
                    freqal.append("NA")

            elif len(alleles)==3: # if the length of alleles is equal 3 (in triply-infected samples)
                if refal not in alleles.keys():
                    hapallele.append("N")
                    qual.append("z")
                    freqal.append("NA")
                    continue
                total=sum(alleles.values())
                #print(alleles.values())
                perc1=(list(alleles.values())[0]/total)
                perc2=(list(alleles.values())[1]/total)
                perc3=(list(alleles.values())[2]/total)
                print(1-(perc1+perc2))
                indx=list(alleles.values()).index(min(alleles.values()))
                print(indx)
                # test if CI fits the data -test whether the minimum allele frequency (the lowest frequency of the three alleles) is within the 90% confidence interval: 
                if min(multinomial_ci([total*perc1, total*perc2, total*perc3])[indx])<=(min(alleles.values())/float(sum(alleles.values()))):
                    print(alleles.keys()) 
                    for k in alleles.keys():
                        if k!=refal[0]:
                            if k!=refal2[0]:
                                if fullalleles[k]>=int(args.m):
                                    hapallele3.append(k)
                                    qual3.append("a3")
                                    freqal3.append(alleles[k]/float(sum(alleles.values())))
                                else:
                                    hapallele3.append("N")
                                    qual3.append("m3")
                                    freqal3.append("NA")
                            else:
                                if fullalleles[k]>=int(args.m):
                                    hapallele.append(k)
                                    qual.append("a")
                                    freqal.append(alleles[k]/float(sum(alleles.values())))
                                else:
                                    hapallele.append("N")
                                    qual.append("m")
                                    freqal.append("NA")
                else:
                    hapallele.append("N")
                    qual.append("/")
                    freqal.append("NA")

            elif len(alleles)==1:
                #print(alleles)
                hapallele.append(list(alleles.keys())[0])
                qual.append("r")
                freqal.append("NA")

            else:
                hapallele.append("N")
                qual.append("z")
                freqal.append("NA")


## add monomorphic positions
def monomorphic():
    for pop,alleles in sorted(popalleles.items()):
            ### test if coverage of individual is above threshold, else append "N"
        if sum(alleles.values())<int(args.mc):
            hapallele.append("N")
            qual.append("c")
        else:
            hapallele.append(refal)
            qual.append("M")



synccode=["A","C","G","T"]
chromo=["tig00000001_segment0_pilon_pilon_pilon_pilon"]
#popstotest=range(3,len(open(args.input).readline().split()[3:])+3)
popstotest=range(3, len(open(args.input).readline().split()[3:])+3)
print(popstotest)
############################# determine haplotype ##########################
out=open(args.o+".consensus","w")
out2=open(args.o+".af","w")
out3=open(args.o+"_ref.af","w")
out3rd=open(args.o+".consensus3","w")
out23rd=open(args.o+".af3","w")
out33rd=open(args.o+"_ref2.af","w")
count=1
for line in open(args.input, 'r'):
    a=line.split()
    chrom=a[0]
    #print(chrom)
    if count%1000000==0:
        print(count,"lines processed")
    count+=1
    #print(count)
    popalleles=collections.defaultdict(lambda:collections.defaultdict(lambda:0))
    fullalleles=collections.defaultdict(lambda:0)
    #print(popalleles)
    #print(fullalleles)

    for pop in popstotest:
        popalleles[pop]
        #print(popalleles[pop])

        ### go through all nucleotides per populations and test whether larger than one, if yes, keep!
        for i in range(len(synccode)):
            if int(a[pop].split(":")[i])>0:

                ## store counts cummulatively for all populations
                fullalleles[synccode[i]]+=int(a[pop].split(":")[i])

                ## store counts for each population separately
                popalleles[pop][synccode[i]]+=int(a[pop].split(":")[i])
                #print(popalleles[4])
                #print(fullalleles)
                #print(len(fullalleles))

    ## if overall coverage is too low just print N's for all pops and continue
    if len(fullalleles)==0:
        out.write("\t".join(a[:3])+"\t"+a[2]+"\t"+"N"*(len(popstotest)-1)+"\t"+"c"*(len(popstotest)-1)+"\t-\n")
        continue

    ########################## extract reference 1  ######################################

    reference1=popalleles[3]
    #print(fullalleles.items())
    #print(popalleles.items())
    for k,val in list(reference1.items()):
        #print(reference.items())
        #print(val)
        if val<=2:
            del reference1[k]

    ## test if reference is above minimum coverage
    if sum(reference1.values())<int(args.mc): #or sum(reference.values())>covdict[3][chrom]:
        refal=a[2]
        #print(refal)
        del popalleles[3]

    ## print "N's"  and exit, if reference is ambiguous
    elif len(reference1)>1:
        total=sum(reference1.values())
        #print(total)

        ## test if minor allele is within the range of 90% interval of a minimum AF of 0.05 because reference is a population of Wolbachia. If it passes, reference in considered polymorphic and therefore print "N's" and exit
        if binomial_ci(total*0.05,total*0.95)[0]<=min(reference1.values())/float(total):
            out.write("\t".join(a[:3])+"\t"+"/".join(reference1.keys())+"\t"+"N"*(len(popstotest)-1)+"\t"+"p"*(len(popstotest)-1)+"\t"+"/".join(fullalleles.keys())+"\n")
            out3.write("\t".join(a[:3])+"\t"+"/".join(reference1.keys())+"\t"+str(min(reference1.values())/float(total))+"\n")
            continue

        ## else if minor alleles does not pass use major allele as reference; reference is not considered polymorphic
        else:
            tar= dict((value,key) for key,value in reference1.items())
            del popalleles[3]
            refal=tar[max(tar.keys())]
    ## if (for whatever reason) there is no allele, take the old reference
    elif len(reference1)==0:
        refal=a[2]
        del popalleles[3]
    ## else use the only allele as the new reference allele
    else:
        refal=list(reference1.keys())[0]
        #print(refal)
        del popalleles[3]
        #print(popalleles)

    ########################## extract reference 2  ######################################
    reference2 = popalleles[4]

    for k,val in list(reference2.items()):
        if count<=2:
            del reference2[k]

    if sum(reference2.values())<int(args.mc): # if the strain refrence 2 is lower than the min coverage threshold
        refal2 = a[2] # ref alleles is taken from the mapping reference
        del popalleles[4]

    ## print "N's"  and exit, if reference is ambiguous
    elif len(reference2)>1:
        total2=sum(reference2.values())

    ## test if minor allele is within the range of 90% interval of a minimum AF of 0.05 because reference is a population of Wolbachia. If it passes, reference in considered polymorphic and therefore print "N's" and exit
        if binomial_ci(total2*0.05, total2*0.95)[0]<=min(reference2.values())/float(total2):
            out3rd.write("\t".join(a[:3])+"\t"+"/".join(reference2.keys())+"\t"+"N"*(len(popstotest)-1)+"\t"+"p"*(len(popstotest)-1)+"\t"+"/".join(fullalleles.keys())+"\n")
            out33rd.write("\t".join(a[:3])+"\t"+"/".join(reference2.keys())+"\t"+str(min(reference2.values())/float(total2))+"\n")

    ## else if minor alleles does not pass use major allele as reference; reference is not considered polymorphic in this case.
        else:
            tar2=dict((k,val) for val,k in reference2.items())
            del popalleles[4]
            refal2=tar2[max(tar2.keys())]

    ## if (for whatever reason) there is no allele, take the old reference.
    elif len(reference1)==0:
        refal2=a[2]
        del popalleles[4]

    ## else use the only allele as the new reference allele.
    else:
        refal2=list(reference2.keys())[0]
        del popalleles[4]



    ################# determine the paternal allele in the indivduals ###########################

    freqal=[]  ## list of allelesfreqs from the haplotypes
    hapallele=[] ## list of alleles from the haplotypes 
    qual=[] ## list of description symbols for each individual
    freqal3=[]
    hapallele3=[]
    qual3=[]

    if len(fullalleles)>1:
        #print(fullalleles)

        CI_ind()

    else:
        monomorphic()

    out.write("\t".join(a[:3])+"\t"+"/".join(refal)+"\t"+"".join(hapallele)+"\t"+"".join(qual)+"\t"+"/".join(fullalleles.keys())+"\n")

    out3rd.write("\t".join(a[:3])+"\t"+"/".join(refal2)+"\t"+"".join(hapallele3)+"\t"+"".join(qual3)+"\t"+"/".join(fullalleles.keys())+"\n")

    if len(freqal)!=0 and list(set(freqal))!=["NA"]:
    #if len(freqal)!=0:
        out2.write("\t".join(a[:2])+"\t"+"\t".join(map(str,freqal))+"\n")

    if len(freqal3)!=0 and list(set(freqal3))!=["NA"]:
        out23rd.write("\t".join(a[:2])+"\t"+"\t".join(map(str,freqal3))+"\n")
