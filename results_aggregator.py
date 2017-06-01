from __future__ import division
import cPickle as pickle
import os
import sys
from collections import Counter

if len(sys.argv) != 2:
    raise Exception("Call this file like so: python results_aggregator.py source_path")

pkl_path = sys.argv[1]
pkl_files = [f for f in os.listdir(pkl_path) if os.path.isfile(os.path.join(pkl_path, f)) and 'pkl' in f]
print 'Found the following files:', pkl_files

unconnected_paths = Counter()
counter_unconnected = 0
connected_paths = Counter()
counter_connected = 0
for f in pkl_files:
    if f[:4] == 'unco':
        unconnected_paths+=pickle.load(open(f,'rb'))        
        counter_unconnected+=1
    elif f[:4] == 'conn':
        connected_paths+=pickle.load(open(f,'rb'))        
        counter_connected+=1
        print pickle.load(open(f,'rb'))
    else:
       print 'Skipping unrecognized filename',f

for key, value in unconnected_paths.items():
    unconnected_paths[key] = value / counter_unconnected
for key, value in connected_paths.items():
    connected_paths[key] = value / counter_connected

print 'Unconnected paths on average'
print unconnected_paths
print '\nConnected paths on average'
print connected_paths


