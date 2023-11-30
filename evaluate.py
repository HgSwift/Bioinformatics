import os, math
import matplotlib.pyplot as plt
import statistics



def readMotif(file):
	motif = []
	for line in file.readlines()[1:]:
		if line[0] != '<':
			sp = line.split()
			split = []
			for s in sp:
				split.append(float(s))
			motif.append(split)
	return motif



def getEntropy(directory):
	motif_file = open(os.path.join(directory,'motif.txt'),'r')
	predicted_file = open(os.path.join(directory,'predictedmotif.txt'),'r')
	motif = readMotif(motif_file)
	#print('ML: ', len(motif))
	predicted = readMotif(predicted_file)
	entropy = []
	sc = motif[0][0] + motif[0][1] + motif[0][2] + motif[0][3] 
	for p in range(len(motif)):
		e = 0.0
		for b in range(4):
			if motif[p][b] != 0 and predicted[p][b] != 0:
				e += (motif[p][b] / sc)* math.log(((motif[p][b]/sc)/(predicted[p][b]/sc)),2)
		e = e/len(motif)	
		entropy.append(e)
	return entropy

def getOverlap(directory):
	site_file = open(os.path.join(directory,'sites.txt'))
	predicted_file = open(os.path.join(directory,'predictedsites.txt'),'r')
	sites = site_file.readlines()
	predicted = predicted_file.readlines()
	overlaps = []
	for i in range(len(sites)):
		o = int(predicted[i]) - int(sites[i])
		overlaps.append(o)
	return statistics.mean(overlaps)

def getRuntime(directory):
	run_file = open(os.path.join(directory,'runtime.txt'))
	return float(run_file.readline())

defaultdir = 'ICPC=2,ML=8,SC=10'
icpc_1 = 'ICPC=1,ML=8,SC=10'
icpc_1_5 = 'ICPC=1.5,ML=8,SC=10'
ml_6 = 'ICPC=2,ML=6,SC=10'
ml_7 = 'ICPC=2,ML=7,SC=10'
sc_5 = 'ICPC=2,ML=8,SC=5'
sc_20 = 'ICPC=2,ML=8,SC=20'

directories = [defaultdir, icpc_1, icpc_1_5, ml_6, ml_7, sc_5, sc_20]


entropies = [[], [], [], [], [], [], []]
overlaps = [[], [], [], [], [], [], []]
runtimes = [[], [], [], [], [], [], []]
entropies_std = []
entropies_mean = []
j = 0
for dir in directories:
	for i in range(1,11):
		entropies[j].append(getEntropy(os.path.join(dir, (str(i)))))
		overlaps[j].append(getOverlap(os.path.join(dir, (str(i)))))
		runtimes[j].append(getRuntime(os.path.join(dir, (str(i)))))
	entropies_mean = [ sum(x) for x in zip(*entropies[j]) ]
	entropies_std = [ statistics.stdev(x) for x in zip(*entropies[j]) ]
	plt.errorbar(range(1,len(entropies_mean)+1), entropies_mean, entropies_std, linestyle='None', marker='^', color='b', ecolor='r')
	plt.xlabel('Position')
	plt.ylabel('Entropy w/ stddev')
	plt.title('Positional entropy for ' + directories[j])
	plt.show()
	plt.clf()
	#print(entropies_mean)
	#print(entropies_std)
	entropies_std = []
	entropies_mean = []
	j += 1
#print(entropies)

#print(entropies_mean)
#print(entropies_std)
#defaultdir, icpc_1, icpc_1_5, ml_6, ml_7, sc_5, sc_20
plt.plot( [1]*len(runtimes[1]), runtimes[1], 'r.', linestyle='None')
plt.plot(1, statistics.mean(runtimes[1]), 'g', linestyle='None', marker='^')
plt.plot( [1.5]*len(runtimes[2]), runtimes[2], 'r.', linestyle='None')
plt.plot(1.5, statistics.mean(runtimes[2]), 'g^', linestyle='None')
plt.plot( [2]*len(runtimes[0]), runtimes[0], 'r.', linestyle='None')
plt.plot(2, statistics.mean(runtimes[0]), 'g^', linestyle='None')
plt.title('Runtime vs ICPC')
plt.xlabel('ICPC')
plt.ylabel('Runtime')
plt.show()
plt.clf()

plt.plot( [6]*len(runtimes[3]), runtimes[3], 'r.', linestyle='None')
plt.plot(6, statistics.mean(runtimes[3]), 'g', linestyle='None', marker='^')
plt.plot( [7]*len(runtimes[4]), runtimes[4], 'r.', linestyle='None')
plt.plot(7, statistics.mean(runtimes[4]), 'g^', linestyle='None')
plt.plot( [8]*len(runtimes[0]), runtimes[0], 'r.', linestyle='None')
plt.plot(8, statistics.mean(runtimes[0]), 'g^', linestyle='None')
plt.title('Runtime vs Motif Length')
plt.xlabel('Motif Length')
plt.ylabel('Runtime')
plt.show()
plt.clf()

plt.plot( [5]*len(runtimes[5]), runtimes[5], 'r.', linestyle='None')
plt.plot(5, statistics.mean(runtimes[5]), 'g', linestyle='None', marker='^')
plt.plot( [20]*len(runtimes[6]), runtimes[6], 'r.', linestyle='None')
plt.plot(20, statistics.mean(runtimes[6]), 'g^', linestyle='None')
plt.plot( [10]*len(runtimes[0]), runtimes[0], 'r.', linestyle='None')
plt.plot(10, statistics.mean(runtimes[0]), 'g^', linestyle='None')
plt.title('Runtime vs Sequence Count')
plt.xlabel('Sequence Count')
plt.ylabel('Runtime')
plt.show()
plt.clf()


plt.plot( [1]*len(overlaps[1]), overlaps[1], 'r.', linestyle='None')
plt.plot(1, statistics.mean(overlaps[1]), 'g', linestyle='None', marker='^')
plt.plot( [1.5]*len(overlaps[2]), overlaps[2], 'r.', linestyle='None')
plt.plot(1.5, statistics.mean(overlaps[2]), 'g^', linestyle='None')
plt.plot( [2]*len(overlaps[0]), overlaps[0], 'r.', linestyle='None')
plt.plot(2, statistics.mean(overlaps[0]), 'g^', linestyle='None')
plt.title('Mean number of Overlapping Sites vs ICPC')
plt.xlabel('ICPC')
plt.ylabel('Overlapping Sites')
plt.show()
plt.clf()

plt.plot( [6]*len(overlaps[3]), overlaps[3], 'r.', linestyle='None')
plt.plot(6, statistics.mean(overlaps[3]), 'g', linestyle='None', marker='^')
plt.plot( [7]*len(overlaps[4]), overlaps[4], 'r.', linestyle='None')
plt.plot(7, statistics.mean(overlaps[4]), 'g^', linestyle='None')
plt.plot( [8]*len(overlaps[0]), overlaps[0], 'r.', linestyle='None')
plt.plot(8, statistics.mean(overlaps[0]), 'g^', linestyle='None')
plt.title('Mean number of Overlapping Sites vs Motif Length')
plt.xlabel('Motif Length')
plt.ylabel('Overlapping Sites')
plt.show()
plt.clf()

plt.plot( [5]*len(overlaps[5]), overlaps[5], 'r.', linestyle='None')
plt.plot(5, statistics.mean(overlaps[5]), 'g', linestyle='None', marker='^')
plt.plot( [20]*len(overlaps[6]), overlaps[6], 'r.', linestyle='None')
plt.plot(20, statistics.mean(overlaps[6]), 'g^', linestyle='None')
plt.plot( [10]*len(overlaps[0]), overlaps[0], 'r.', linestyle='None')
plt.plot(10, statistics.mean(overlaps[0]), 'g^', linestyle='None')
plt.title('Mean number of Overlapping Sites vs Sequence Count')
plt.xlabel('Sequence Count')
plt.ylabel('Overlapping Sites')
plt.show()
plt.clf()