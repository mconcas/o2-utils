#! /usr/bin/env python

import uproot
from sklearn.cluster import DBSCAN
from numpy import array

cent = uproot.open("vertexer_serial_data.root")["centroids"].arrays(["id","x","y","z"])
eventDataReady = array([ [x[1], x[2], x[3]] for x in array(filter(lambda x : x[0]==20, array((cent["id"], cent["x"], cent["y"], cent["z"])).transpose()))])
clustering = DBSCAN(eps=.1, min_samples=16).fit(array([ [x[1], x[2], x[3]] for x in array(filter(lambda x : x[0]==14, array((cent["id"], cent["x"], cent["y"], cent["z"])).transpose()))]))
# print(eventDataReady)

