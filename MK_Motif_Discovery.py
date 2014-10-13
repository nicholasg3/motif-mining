#!/usr/bin/env python2.7
import numpy as np
import random

def distance(x, y):
	sum = 0
	i = 0
	while (i < len(x)):
		sum = sum + (x[i] - y[i])**2
		i = i + 1
	return (sum ** 0.5)

# MK Motif Discovery
# @param R, int number of reference points
# @param D, list, time series
# @param tempLen, the length of the Template to Discover
# @param r, the minimum dsitance between two motifs to be considered
def MK_Motif(D,R, tempLen = 30):
	best_so_far = 999999999.0
	D = np.array(D)
	m = len(D)-tempLen # the number of time series to compare
	S = []
	Dist = []
	for i in range(R):
		r = random.randint(0,len(D)-tempLen)
		ref_i = D[r:r+tempLen] # a randomly chosen time series Dr from D
		Dist_i = []
		for j in range(m):
			Dj = D[j:j+tempLen]
			Dist_ij = distance(ref_i,Dj) # euclidian or DTW distance between two time series
			Dist_i.append(Dist_ij)

			if Dist_ij < best_so_far:
				best_so_far = Dist_ij
				# print "best_so_far: "+ `best_so_far`
				L1 = r
				L2 = j
		Si = np.std(Dist_i)
		S.append(Si)
		Dist.append(Dist_i)
	# find an ordering "Z" of the indicdes to the reference time series ref such that S_z(i) >= S_z(i+1)
	S_arr = np.array(S) # convert list to array
	Z = (-S_arr).argsort()
	# find an ordering "I" of the indices to the time series in D such that Dist_z(1),I(j) <= Dist_z(1),I(J+1)
	Dist = np.array(Dist)

	I = (-Dist[Z[0]]).argsort()
	offset = 0
	abandon = False

	while abandon == False:
		offset = offset + 1
		abandon = True

		for j in range(R):
			reject = False
			for i in range(R):
				# print 
				# print "i " + `i`
				# print "j" + `j`
				# print "offset " + `offset`
				# print "len[I]: " + `len(I)`
				lower_bound = abs(Dist[Z[i],I[j]] - Dist[Z[i],I[j+offset]])
				if lower_bound > best_so_far:
					reject = True
					break
				elif i == 1:
					abandon = False
			if reject == False:
				check_d = distance(D[I[j]:(I[j]+tempLen)],D[I[j+offset]:(I[j+offset]+tempLen)])
				# print "check_d: " + `check_d`
				if check_d < best_so_far:
					best_so_far = check_d
					L1 = I[j]
					L2 = I[j+offset]
	return [L1, L2]