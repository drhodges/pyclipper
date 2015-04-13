# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 18:44:51 2011
@author: Andrew Hodges
"""

def openReportFile(filename):    
    #opens & reads files.    
    fileHandle = open(filename,'r');
    fileInfo = fileHandle.readlines();
    fileHandle.close();
    return fileInfo;

def generateFileNames(root,indexes):
    #filename extraction (note: you may need to write your own!)    
    out = [];
    for i in indexes:
        out.append(str(root)+str(i)+".report.txt");
    return out;

def getNodes(filename):
    #get the nodes from a file with node names separated by "\n"
    fileHandle = open(filename,'r');
    nodesx = fileHandle.readlines();
    nodes = [];
    #genes = genes.rstrip();
    for x in nodesx:
        x = str(x)        
        y=x.rstrip();        
        nodes.append(y);
    fileHandle.close();    
    return nodes;

def getNetworks(fnames):
    #function will return all networks with no parsing
    networks = [];    
    for i in fnames: #change this later to use indexes array
        #indexes:
        fileHandle = open(i,'r');
        network = fileHandle.readlines();
        fileHandle.close();
        #print network;
        #print scores;
        networks.append(network);
    return networks;

def getNetsAll(files,outfile,os,re):
    #Purpose: get all of the networks from the BANJO files
    #Banjo files should all be saved in the same directory   
    scores = {};
    fileOutArray = [];
    i = -1;
    #fileOut = open(outfile,"wb");    
    netcount = 0;
    i=i+1;
    for fileinx in files:
        x = 0;
        edges = [];
        BanjoFile = openReportFile(fileinx);        
        for line in BanjoFile:
            linex = str(line);
            linex = linex.rstrip();        
            if x:
                #currently in the list of top-scoring networks
                m = re.search("Network .*, score: ([0-9\-\.]+),",linex); 
                if m:
                    score = m.group(1);
                    if len(edges) >= 1:
                        fileOutArray.append(";"+",".join(edges)+"\n");
                        #fileOut.write(";"+",".join(edges)+"\n");
                    fileOutArray.append(score);
                    #fileOut.write(score);
                    edges=[]
                    if not scores.has_key(str(score)):
                        scores[str(score)]=1;
                    else:
                        scores[str(score)]= scores[str(score)] + 1;
                elif re.search("Search Statistics",linex):
                    if len(edges) >= 1:
                        fileOutArray.append(";"+",".join(edges)+"\n");
                        #fileOut.write(";"+",".join(edges)+"\n");
                        x=0;
                        break;
                else:
                    m = re.search("([\d]+) ([\d]+)",linex);                    
                    if m and x:
                        if int(m.group(2)) > 0:
                            entries = linex.split();
                            child = entries.pop(0);
                            for parent in entries:
                                edges.append(parent+"->"+child);
            else:
                m = re.search("These are the ([0-9]+) ",linex);
                if m:
                    netcount = m.group(1);
                else:
                    m = re.search("Best ([\d]+) Structures",linex);
                    if m:
                        netcount = m.group(1);
                        x = 1;
                    else:
                        m = re.search("\- Final report",linex);
                        if m:
                            x = 1;
    #print edges
    scorex = scores.keys();
    fileOut = open(outfile,"wb");
    fileOut.write("".join(fileOutArray));
    fileOut.close();    
    return (scorex, netcount, scores)

def reorderNetworks(infile, currS, re):
    #pass in regex excpression... write each network only if it matches the score regex
    fileHandle = open(infile,'r');
    fileInfo = fileHandle.readlines();
    fileHandle.close();
    returnArray = [];
    for i in fileInfo:
        if re.match(currS, i):
            returnArray.append(i);
    return returnArray;

def rewriteScoresFile(infile,outfile,scores,re):
    #function to parse the network/score file... wrapper for reorderNetworks funtion    
    scoresSorted = scores;
    #scoresSorted.sort();
    #scoresSorted.reverse();
    allReturnedInfo = [];
    for currS in scoresSorted:
        currS = re.sub("[\.]","\.",currS)+";";
        R = reorderNetworks(infile,currS,re);
        allReturnedInfo.append("".join(R));
    fileOut = open(outfile,"wb");    
    fileOut.write("".join(allReturnedInfo));
    fileOut.close();
    return 1;

def grabScores(infile):
    #function gets the parsed score information and returns the scores     
    scores = [];    
    netInfo = openReportFile(infile);
    for i in netInfo:
        i = i.split(";")
        scores.append(i[0]);
    return scores;        

def getMaxScore(listx):
    #function grabs the max score and returns it    
    scores = [];
    for i in listx:
        scores.append(float(i));
    return max(scores);

def bvalCalculator(infile,e):
    #function calculates the b-values    
    #initialize variables    
    allScores = {};
    uniqScores = [];
    combined_probabilities = 0.0;
    scoreChanges = {};
    bval = {};    
    bval_mapper = {};
    #executable portion
    #first get the unique scores and max value  
    scores = grabScores(infile);
    for score in scores:
        allScores[score]=1;
    uniqScores = allScores.keys()
    uniqScores.sort(); #python uses ascending score size
    #uniqScores.reverse();
    #uniqScores = uniqScores[0:10];
    maxscore = getMaxScore(uniqScores);
    #Step one of computation... compute score change
    for key in range(len(uniqScores)):
        scoreChanges[key]=e**(float(uniqScores[key])-maxscore);
        combined_probabilities += scoreChanges[key];
        #x = key;
    #step 2: generate right-tail density for each score cutoff
    for key in range(len(uniqScores)):
        bval=0.0;
        for key2 in range(len(uniqScores)-1,-1,-1):
            if float(uniqScores[key])>float(uniqScores[key2]):
                bval += scoreChanges[key2]/combined_probabilities;
            else:
                break;
        bval_mapper[str(uniqScores[key])] = bval;
    return bval_mapper;
    
def percent_appear(networks={},count=0,edgePercent=1.0):
    #function generates the decisions for whether an edge should be included or not
    #major assumption is that an edge must be present with freq = 1.0 to include it    
    sortEdges = {};
    final_edges = [];
    count = edgePercent * count;
    for edge in networks.keys():
        nodes = edge.split("->");
        nodes.sort();
        edger = str(nodes[0])+"->"+str(nodes[1]);
        sortEdges[edger]=1;
    #sorted_edges = sortEdges.keys();
    for edge in sortEdges.keys():
        numOccur = 0;
        nodes=edge.split("->");
        edge2=str(nodes[1])+"->"+str(nodes[0]);
        if networks.has_key(edge):
            numOccur += float(networks[edge]);
        if networks.has_key(edge2):
            numOccur += float(networks[edge2]);
        if numOccur >= count:
            final_edges.append(edge);
    return final_edges;            

def bvalConsensus(score,infile,edgePercent):
    #function generates the consensus network given the b-values
    #note: function is a wrapper for the percent_appear function    
    consensusEdges = {};
    y = 0;
    include_edges = {};
    z = 0;
    fileInfo = openReportFile(infile);
    for string in fileInfo:
        y+=1;
        scorei = 0;
        edgestr = "";
        string = str(string)        
        string=string.rstrip();
        if not string == "":
            info = string.split(";");
            scorei = info[0];
            edgestr = info[1];
            if float(scorei)>=float(score):
                z+=1;
                edges = [];
                edges = edgestr.split(",");
                for edge in edges:
                    if not include_edges.has_key(edge):
                        include_edges[edge]=0;
                    include_edges[edge]+=1;
            else:
                consensusEdges[score]="";
                consensusEdges[score] = percent_appear(include_edges,z,edgePercent);
                #consensusEdges[score]=
        elif not consensusEdges.has_key(score):
            consensusEdges[score] = percent_appear(include_edges,z,edgePercent);
    return consensusEdges[score];

def readEdgeInformation(filename,nodes,writer):
    #function generates the formatted information for the c-value output    
    cutoffArray = {};    
    #counter = len(nodes);
    includeEdges = [];
    lines = openReportFile(filename);
    string = lines[1];  #line 0 has the header, while line 1 has the top-level network
    string = string.rstrip();
    info = string.split("\t");
    edgesx = info[2];
    edges = edgesx.split(",");
    for edge in edges:
        edger = edge.split("->");
        cutoffArray[edge]=[info[1],nodes[int(edger[0])],nodes[int(edger[1])],0,0]
        includeEdges.append(edge);
    print cutoffArray;
    for i in range(2,len(lines)):
        string = lines[i];
        string = string.rstrip();
        if string != "":
            vals = string.split("\t");
            if len(vals)==3:            
                bval = vals[1];
                currentEdges = vals[2];
                print str(bval)+"\t"+str(currentEdges)+"\n";
                currentEdges = currentEdges.split(",");
                for edge in includeEdges:
                        if cutoffArray[edge][3]==0 and not edge in currentEdges:
                            cutoffArray[edge][3]=1;
                            cutoffArray[edge][4]=bval;
            else:
                break;
    return cutoffArray;

def bvalWrapper(outputBvals, bvalFile, netfile2,edgePercent):
    #main wrapper for b-value computation and output    
    outputArray = [];
    scoresN = [];
    scoresN = outputBvals.keys();
    scoresN.sort();    
    for score in scoresN:
        edgeInfo = bvalConsensus(str(score),netfile2,edgePercent);
        stringx = str(score)+"\t"+str(outputBvals[str(score)])+"\t"+str(",".join(edgeInfo))+"\n";
        outputArray.append(stringx);
    fout = open(bvalFile,'wb');
    fout.write("Score\tBval\tEdges\n");
    fout.write("".join(outputArray)+"\n");
    fout.close();
    return 1;

def cvalWrapper(bvalFile, cvalFile, nodes):    
    #main wrapper for c-value computation and output    
    outputArray = [];
    outputCvalues = readEdgeInformation(bvalFile,nodes);
    for edge in outputCvalues.keys():
        outputArray.append(edge+"\t"+"\t".join(outputCvalues[edge])+"\n");
        #fHandleCval.write(edge+"\t"+"\t".join(outputCvalues[edge])+"\n");
    fHandleCval = open(cvalFile,"wb");
    fHandleCval.write("".join(outputArray)+"\n");
    fHandleCval.close();
    return 1;

if __name__ == "__main__":
    import sys
    from sys import argv
    import os
    