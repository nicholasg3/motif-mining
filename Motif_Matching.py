#!/usr/bin/env python2.7
import numpy as np
import csv
import time


def zNormalize(my_time_series):
    data = np.array(my_time_series)
    y_normed = (data-np.mean(data))/np.std(data)
    return y_normed

# adapted from Mueen and Keogh 2012
# Q query: a list of numbers. The motif to be matched to.
# filePathT: A the path pointing to a csv containing one number per row. 
# sample_rate: an integer. If you wish to sample every nth value of the time series, set this to n. Note that the matching motif much also be collected at this same sampling rate.
# normalize: If true, compare distaces between z-normalized query and test-data. If false compare absolute distance
# returns the first index of the closest fit, nn
def find_matches(filePathT,Q,sample_rate=1, normalize = True):
    with open(filePathT) as tsFile:
        T = csv.reader(tsFile)

        bestSoFar = float("inf")
        count = 0 # the 'step' in the time series
        if normalize:
            Q = zNormalize(Q)
        m = len(Q)
        X = [0]*m
        ex = 0; ex2 = 0
        
        try:
            tnext = float(T.next()[0]) # This implemntation assumes the CSV stores one value per line
        except StopIteration:
            print "empty list"
            return 0
        else:
            hasNext = True
            while hasNext: #can choose to limit the depth you look into the TS with "and count < 100000:"
                count = count + 1
                if count % sample_rate == 0:
                    i = count/sample_rate % m
                    X[i] = tnext # circular buffer to store current subsequence. 
                    ex = ex + X[i] # iteratively sum up values to use for the mean
                    ex2 = ex2 + X[i] ** 2 # sum up squared values to use for sdev

                    if count >= m-1:
                        u = ex/m # u is mu, or the mean
                        sdev = abs(ex2/m - u**2) ** (0.5) 
                        j = 0
                        dist = 0

                        # compare Q and T[i]
                        while j < m and dist < bestSoFar:
                            if normalize:
                                dist = dist + (Q[j]-(X[(i+j)%m]-u)/sdev)**2
                            else:
                                dist = dist + abs(Q[j]-(X[(i+j+1)%m]))
                            j = j + 1

                        if dist < bestSoFar:
                            bestSoFar = dist
                            nn = count - m # count gives the end of the matched motif. Move m spaces back to finds its head.
                            
                        # keep the mean and sdev and moving averages.
                        ex = ex - X[(i+1)%m]
                        ex2 = ex2 - X[(i+1)%m]**2
                try:
                    tnext = float(T.next()[0])
                except StopIteration:
                    print "end of list"
                    hasNext = False
        return nn # closest match spot in time series T