import numpy as np
import matplotlib.pyplot as plt
 
def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)


a = np.array((1, 0.5))
b = np.array((2.5, 3))
c = np.array((2, 1))
d = np.array((3, 2))
e = np.array((3.5, 2))
f = np.array((-0.5, 0))
g = np.array((-0.5, 1))
h = np.array((-1, 0.5))
i = np.array((1, -1))
j = np.array((0.5, -1))
points = [a, b, c, d, e, f, g, h, i, j]
cen1 = np.array((0,0))
cen2 = np.array((2, 3))
cen3 = np.array((1.5, -1))
centers = [cen1, cen2, cen3]
K = 3
plt.plot(cen1[0], cen1[1], 'b^')
plt.plot(cen2[0], cen2[1], 'r^')
plt.plot(cen3[0], cen3[1], 'g^')
cmap = plt.cm.get_cmap('hsv', 4)
distances = np.zeros([K])
for x in range(3):
    clusters = [[], [], []]
    for point in points:
        for i in range(K):
            distances[i] = np.linalg.norm(point - centers[i])
        minid = np.argmin(distances)
        clusters[minid].append(point)

    print(clusters)
    # naming the x axis
    plt.xlabel('X')
    # naming the y axis
    plt.ylabel('y')
    
    # giving a title to my graph
    plt.title('Clusters')
    xcoords = [[], [], []]
    ycoords = [[], [], []]
    clusnum = 0
    for cluster in clusters:
        for point in cluster:
            xcoords[clusnum].append(point[0])
            ycoords[clusnum].append(point[1])
        clusnum += 1

    for i in range(K):
        plt.plot(xcoords[i], ycoords[i], 'o', c=cmap(i))
        plt.plot(centers[i][0], centers[i][1], '^', c=cmap(i))
    xc = [[], [], []]
    yc = [[], [], []]
    for i in range(K):
        xc[i]=(sum(xcoords[i])/len(xcoords[i]))
        yc[i]=(sum(ycoords[i])/len(ycoords[i]))
        centers[i] = np.array((xc[i], yc[i]))
        

    for i in range(K):
        print('Cluster: ', i+1, clusters[i])
        print('Center: ', i+1, centers[i])
    plt.show()
