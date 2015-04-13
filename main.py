# -*- coding: utf-8 -*-
"""
EdgeClipper Algorithm
Author: Andrew P. Hodges
Created on Mon Apr 25 16:54:46 2011
Note: translated from php code in MARIMBA
"""

#Open system variables
from sys import argv
import os
import re

#Major steps:
"""###Main function###"""
"""global variables"""
timestamp = "20101105124824";
maindir = "F:\\SomeDirectory\\";
scriptDir = maindir + "scripts\\";
workingDir = maindir + "reportFiles\\";
outputDir = maindir + "outputFiles\\";

"""1) first import dependencies"""
os.chdir(scriptDir);
#if 1:
import funsEC;
#else:
#    reload(funsEC);
os.chdir(workingDir);

"""#2) Next get list of nodes"""
#Note: may need additional function to reverse the index->gene to gene->index
file_nodes = timestamp+".probelist.txt";
nodes = funsEC.getNodes(file_nodes);

"""#3) Next get set of networks from the BANJO analyses"""
#directory = diry + 'testReports\\';
#Note: this section should be adjusted for your particular file naming strategy!
filenames = funsEC.generateFileNames(timestamp,range(0,20));
#networks = funsEC.getNetworks(filenames);

"""#4) generate the condensed network information file with sorted scores"""
netfile1 = "tempNetFile.txt";
netfile2 = "tempNetFile2.txt";

"""#5 rewrite the networks into a parsable format"""
networkScores = funsEC.getNetsAll(filenames,netfile1,os,re);

#Next reformat the testNetOut temporary file so that networks are ordered by score
scoresUniq = networkScores[2].keys();
scoresUniq.sort();
funsEC.rewriteScoresFile(netfile1,netfile2,scoresUniq,re);

#"""#6 get b-value information"""
outputBvals = funsEC.bvalCalculator(netfile2,e);

#"""#7 write b-value info and corresponding consensus networks to file"""
edgePercent = 1.0;
bvalFile = "bvalEdgeInfo.txt";
funsEC.bvalWrapper(outputBvals, bvalFile, netfile2,edgePercent);

#"""finally, generate c-values for edges and store to file"""
cvalFile = "cvalEdgeInfo.txt";
funsEC.cvalWrapper(bvalFile, cvalFile, nodes);
