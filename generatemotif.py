import sys
import random
import time
import os

## Constants
bases = ['A','C','G','T']
freq = [0.25,0.25,0.25,0.25]
filename='uniform_bg_ICPC_to_p.txt'
sequence_filename='sequences.fa'
site_filename='sites.txt'
folder='dataset'

def create_seq(SL, freq = [0.25, 0.25, 0.25, 0.25]):
    seq=[]
    for i in range(0,SL):
        x=random.randint(1, 4)
        y = x/4
        if (y<=freq[0]):
            seq.append(bases[0])
        elif (y<=sum(freq[0:2])):
            seq.append(bases[1])
        elif (y<=sum(freq[0:3])):
            seq.append(bases[2])
        else:
            seq.append(bases[3])
    return ''.join(seq)


def plant_motifs(seq,motifs,ML,SL, SC):
	locations=random.randrange(0,SL-ML)
	motif_count=random.randrange(0,SC)
	plant_seq=seq[:locations]+motifs[motif_count]+seq[locations+ML:]
	return plant_seq,locations
	
def write_to_fa(SC, seq, filename):
    ## Write plnted sequences to sequences.fa file
    f=open(filename,'w')
    for i in range(0,SC):
        f.write('>seq'+str(i+1)+'\n'+seq[i]+'\n')
    f.close()



def build_bench(ICPC, ML, SC, z, SL = 500):
    # Create SC sequences
    seqs=[]
    for i in range(0,SC):
        temp=create_seq(SL,freq)
        seqs.append(temp)
    # Create motif and add to motifs.txt
    #f_motif=open('motif.txt','wb')
    #f_motif.write('>'+gen_seq+'\t'+str(ML)+'\n')
    # Assign p given ICPC
    if ICPC == 1:
        p = 0.8105
    elif ICPC == 1.5:
        p = 0.9245
    else:
        p = 1
    q = (1-p)/3
    # make random motif
    gen_seq=create_seq(ML)
	## Create SC motifs of length ML
    motifs_temp=[]
    for i in range(0,ML):
        pq_freq=[q,q,q,q]
        for j in range(0,4):
            if(gen_seq[i]==bases[j]):
                pq_freq[j]=p
        temp=create_seq(SC, pq_freq)
        motifs_temp.append(temp)
    motifs=[]
    for i in range(0,SC):
        temp=[]
        for j in range(0,ML):
            temp.append(motifs_temp[j][i])
        motifs.append(''.join(temp))
    ## Planting motifs in seqs
    planted_seq=[]
    sites=[]
    for i in range(0,SC):
        temp,site=plant_motifs(seqs[i],motifs,ML,SL, SC)
        planted_seq.append(temp)
        sites.append(site)

    write_to_fa(SC, planted_seq, os.path.join(str(z), 'sequences.fa'))
    #Write 
    f=open(os.path.join(str(z), 'sites.txt'),'w')
    for i in range(0,SC):
        f.write(str(sites[i])+'\n')
    f.close()

    #Create PWM
    pwm=[]
    for i in range(0,ML):
        count=[0,0,0,0]
        for j in range(0,SC):
            for k in range(0,len(bases)):
                if (motifs[j][i]==bases[k]):
                    count[k]=count[k]+1
        pwm.append(count)
    #Write PWM
    f = open(os.path.join(str(z), 'motif.txt'), 'w')
    f.write('>'+gen_seq+'\t'+str(ML)+'\n')
    for i in range(0,len(pwm)):
        for j in range(0,len(pwm[0])):
            f.write(str(pwm[i][j])+'\t')
        f.write('\n')
    f.write('<')
    f.close()
	## Write motiflength
    f=open(os.path.join(str(z), 'motiflength.txt'),'w')
    f.write(str(ML))
    f.close()
for i in range(1, 11):
    build_bench(2, 8, 20, i)
    #print(i, end=' ')
# build_bench(1.5, 8, 10, 2)
# build_bench(2, 8, 10, 3)
# build_bench(2, 6, 10, 4)
# build_bench(2, 7, 10, 5)
# build_bench(2, 8, 5, 6)
# build_bench(2, 8, 10, 7)
# build_bench(2, 8, 20, 8)
