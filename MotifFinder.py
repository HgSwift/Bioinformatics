from typing import Counter
from Bio import SeqIO
import sys
import os
import shutil
import random
import time

bases = ['A','C','G','T']
path = "/Users/ralhazmy/Documents/Spring 22/Bioinfo/AustinRaghidMotif-main/ParametersGenerated/"

resultsFile = "/Users/ralhazmy/Documents/Spring 22/Bioinfo/AustinRaghidMotif-main/results/"
if os.path.exists(resultsFile):
    shutil.rmtree(resultsFile)
os.makedirs(resultsFile)


#Finds the probability of background
def sequenceProbability(seq):
    count = [0, 0, 0, 0]
    prob = [0, 0, 0, 0]
    for base in seq:
        if base == 'A':
            count[0] += 1
        elif base == 'C':
            count[1] += 1
        elif base == 'G':
            count[2] += 1
        elif base == 'T':
            count[3] += 1
    total = sum(count)
    pribIndex =0
    for l in count:
        prob[pribIndex] = l/total
        pribIndex += 1
    return count

# creates a PWM
def makePWM(motifs, seqList, motifLen):
    print(motifs)
    print(len(seqList))
    pwm = []
    for i in range(0, motifLen):
        count = [0,0,0,0]
        prob = [0,0,0,0]
        for j in range(0,len(seqList)):
            if motifs[j] == '':# motif being replace
                continue
            for k in range(0,len(bases)):
                if (motifs[j][i] == bases[k]):
                    count[k] = count[k] + 1
                    break
        total = sum(count)
        pribIndex =0
        for l in count:
            prob[pribIndex] = l/total
            pribIndex += 1
        pwm.append(count)
    return pwm

def gibbsAlgo(seqsFile, MLFile, dir,index):
    #Open seqs file + Append sequences to a list
    start_time = time.time()
    seqList = []
    seqs = open(seqsFile,'r')
    for sequence in SeqIO.parse(seqs,'fasta'):
        seqList.append(str(sequence.seq))
    seqLength = len(seqList[0])
    
    #Open MLFile
    motifLen = 0
    with open(MLFile) as f:
        motifLen = int(f.read())
    
    #Initliaze gibbs algorithm by choosing a random site as the motif for each sequence
    currSolution = []
    sites = []
    motifs = []
    for sequence in seqList:
        site = random.randint(0, seqLength - motifLen - 1)
        sites.append(site)
        motif = sequence[site: site + motifLen]
        motifs.append(motif)

    loop = True
    idk=0
    while(loop):
        #choose one motif and replace
        motif2Replace = random.randint(0, len(motifs) - 1)
        sites[motif2Replace] = ''
        motifs[motif2Replace] = ''
        seqZ = seqList[motif2Replace]

        #make pwm from the remaining motifs
        pwm = makePWM(motifs, seqList, motifLen)
        #For each candidate site x  in sequence z, calculate Qx and Px:
        #Qx: the 	probability of generating x according to θt(pwm); 	
        #Px: the	probability of generating x according to the background model
        bestCandid = []
        for candSite in range(0, seqLength - motifLen):
            candMotif = seqZ[candSite : motifLen+candSite]
            i = 0
            Qx = 1
            for motifBase in candMotif:
                if motifBase == 'A':
                    Qx = Qx * pwm[i][0]
                elif motifBase == 'C':
                    Qx = Qx * pwm[i][1]
                elif motifBase == 'G':
                    Qx = Qx * pwm[i][2]
                elif motifBase == 'T':
                    Qx = Qx * pwm[i][3]
                i += 1
            
            bkgrndProb = sequenceProbability(seqZ)
            Px = 1
            for motifBase in candMotif:
                if motifBase == 'A':
                    Px = Px * bkgrndProb[0]
                elif motifBase == 'C':
                    Px = Px * bkgrndProb[1]
                elif motifBase == 'G':
                    Px = Px * bkgrndProb[2]
                elif motifBase == 'T':
                    Px = Px * bkgrndProb[3]

            #Among all possible candidates, choose one (say x) with probability proportional to Qx/Px
            proportional = Qx/Px
            #print(proportional, Qx, Px)
            if candSite == 0:
                bestCandid = [candSite, candMotif, proportional]
            elif proportional > bestCandid[2]:
                bestCandid = [candSite, candMotif, proportional]


        sites[motif2Replace] = bestCandid[0]
        motifs[motif2Replace] = bestCandid[1]

        """
        Local Optimal solution
        """
        idk+=1
        if idk == 500:
            loop = False

    bestMotifFound = makePWM(motifs, seqList, motifLen)

    dir = dir[87: ]
    f = open(resultsFile + str(index) + 'predictedmotif.txt', 'w')
    f.write('>'+ 'Motif' +'\t'+str(motifLen)+'\n')
    for i in range(0,len(bestMotifFound)):
        for j in range(0,len(bestMotifFound[i])):
            f.write(str(bestMotifFound[i][j])+'\t')
        f.write('\n')
    f.write('<')
    f.close()
#resultsFile + dir + "/"
    f=open(resultsFile + str(index) + 'predictedsites.txt','w')
    for i in range(0,len(seqList)):
        f.write(str(sites[i])+'\n')
    f.close()
    total_time = time.time()-start_time
    f=open(resultsFile + str(index) + 'runtime.txt','w')
    f.write(str(total_time))
    f.close()

def main():
    i = 0
    j = 1
    index=0
    for root, dirs, files in os.walk(path):
        seqsFile = ""
        motifLength = ""
        for file in root.splitlines():
            print()
            print(file)
            if i < 2 or j > 10:
                i += 1
                j = 1
            else:
                seqsFile = file+"/"+"sequences.fa"
                motifLength = file+"/"+"motiflength.txt"
                gibbsAlgo(seqsFile, motifLength, file, index)
                index+=1
                j += 1

if __name__ == '__main__':
    main()