#
## This script in an internal Python script designed to collate all of the scripts previously generated from the pipeline.
## It is intended to be the 2nd to last script run in the WGS alignment pipeline.
## This script is designed to run in a non-cluster environment but NOT on the cluster login node. Open an interactive node or run this script on the file system nodes.
#
import sys
import glob
import fnmatch
import pathlib
import subprocess
import time
import socket
import shutil
import os
import gzip
import numpy as np
#
### SETUP ##
#
## Obtain the exported environment variables passed from the parent bash shell script that invoked this subscript
#
file_dir = os.getenv('file_dir') # This imports the exported variable from the parent shell script of the samples location
txt_dir = os.getenv('txt_dir')
bqsr_dir = os.getenv('bqsr_dir')
gvcf_dir = os.getenv('gvcf_dir')
## Helper function to run commands and return values.
#
def runCMD(cmd):
    val = subprocess.Popen(cmd, shell=True).wait()
    if val == 0:
        pass
    else:
        print('command failed')
        print(cmd)
        sys.exit(1)
#
## Change directories to get old file name and then in the dictionary, strip the extension to create the old file name for use in compression ratio stats.
#
os.chdir(bqsr_dir)
#
bqsrname = []
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*.BQSR.bam'):
        print(file)
        bqsrname = file
#
print(bqsrname)
#
oldname = []
#
oldname = bqsrname.removesuffix('.BQSR.bam') 
#os.system('cd $gvcf_dir')
#oldgvcf = []
#for file in glob.glob(gvcf_dir+'*_g.vcf.gz'):
#    oldgvcf.append(os.path.basename(file))
#
os.chdir(txt_dir)
#
## Parse the rename file to get the string of the internal lab IDs which also are the prefixes to the files required to combine stats.
#
with open('rename.txt', 'r') as file:
    samplename = file.read().rstrip()
#
## Create a dictionary in order to hold the samplename to be passed.
myData = {}
myData['file_dir'] = file_dir
myData['txt_dir'] = txt_dir
myData['bqsr_dir'] = bqsr_dir
myData['gvcf_dir'] = gvcf_dir
myData['oldname'] = oldname
#myData['oldgvcf'] = oldgvcf
myData['samplename'] = samplename
myData['knownSitesVCF'] = myData['txt_dir'] + '/' + myData['samplename'] + '.knownsites.vcf.gz'
myData['alignMetricsFileName'] = myData['txt_dir'] + '/' + myData['samplename'] + '.alignment_summary_metrics' 
myData['sizeMetricsFileName'] = myData['txt_dir'] + '/' + myData['samplename'] + '.insert_size_metrics.txt'
myData['CRAM'] = myData['file_dir'] + '/' + myData['samplename'] + '.BQSR.cram'
myData['bamsize'] = myData['file_dir'] + '/BQSR/' + myData['oldname'] + '.BQSR.bam'
myData['gvcf'] = myData['file_dir'] + '/' + myData['samplename'] + '.g.vcf.gz'
myData['gvcfsize'] = myData['file_dir'] + '/gVCF/' + myData['oldname'] + '_g.vcf.gz'
#
####
#print('Combining Stats for %s' % myData['samplename'], flush=True)
def calc_effective_depth(myData):
    print('Running calc effective depth', flush=True)
    myData['depthSummary'] = myData['samplename'] + '.knownsites.depth.txt'

    if os.path.isfile(myData['depthSummary']) is True:
        print('%s exists' % myData['depthSummary'],flush=True)
        return
    autoDp = []
    xDp = []
    inFile = gzip.open(myData['knownSitesVCF'],'rt')
    for line in inFile:
        if line[0] == '#':
            continue
        line = line.rstrip()
        line = line.split()

        infoField = line[7]
        infoField = infoField.split(';')
        dp = -1
        for i in infoField:
            if i[0:3] == 'DP=':
                dp = int(i.split('=')[1])
                break

        if dp == -1:
            dp = 0
        if line[0] == 'chrX':
            xDp.append(dp)
        else:
            autoDp.append(dp)
    inFile.close()

    outFile = open(myData['depthSummary'],'w')

    outFile.write('#chrom\ttotalSites\tMean\tStd\tMedian\n')
    outFile.write('Autos\t%i\t%.2f\t%.2f\t%.1f\n' % (len(autoDp),np.mean(autoDp),np.std(autoDp),np.median(autoDp) ))
    outFile.write('ChrX\t%i\t%.2f\t%.2f\t%.1f\n' % (len(xDp),np.mean(xDp),np.std(xDp),np.median(xDp) ))
    outFile.close()
#####
def summarize_stats(myData):
    myData['statsSummary'] = myData['samplename'] + '.stats.txt'

    outFile = open(myData['statsSummary'],'w')

    sn = myData['samplename'].split('/')[-1].split('.')[0]

    outFile.write('SampleName\t%s\n' % (sn))

    cramFileSize = os.path.getsize(myData['CRAM']) / (1024*1024*1024)

    outFile.write('CramSize\t%.2f Gb\n' % cramFileSize)

    bamFileSize = os.path.getsize(myData['bamsize']) / (1024*1024*1024)

    bamcramRatio = (int(cramFileSize) / int(bamFileSize) ) * 100

    outFile.write('Cram Compression Ratio\t%.2f%%\n' % bamcramRatio)

    oldgvcfFileSize = os.path.getsize(myData['gvcfsize']) / (1024*1024*1024)

    gvcfFileSize = os.path.getsize(myData['gvcf']) / (1024*1024*1024)

    outFile.write('GVCFSize\t%.2f Gb\n' % gvcfFileSize)

    gvcfRecompressRatio = (int(gvcfFileSize) / int(oldgvcfFileSize)) * 100

    outFile.write('gVCF Compression Ratio\t%.2f%%\n' % gvcfRecompressRatio)

    # now get coverage stats

    inFile = open(myData['depthSummary'],'r')
    lines=inFile.readlines()
    inFile.close()

    autoLine = lines[1].rstrip().split()
    xLine = lines[2].rstrip().split()

    outFile.write('effectiveAutoMean\t%s\neffectiveAutoMedian\t%s\n' % (autoLine[2],autoLine[4]) )
    outFile.write('effectiveXMean\t%s\neffectiveXMedian\t%s\n' % (xLine[2],xLine[4]) )

    xvsAutoMean =  float(xLine[2]) / float(autoLine[2])

    if float(autoLine[4]) == 0.0:
        xvsAutoMedian = 0.0
    else:
        xvsAutoMedian = float(xLine[4]) / float(autoLine[4])

    outFile.write('Mean(X/Auto)\t%.2f\n' % xvsAutoMean)
    outFile.write('Median(X/Auto)\t%.2f\n' % xvsAutoMedian)

    # dup mark stats

    dupMetFileName = myData['samplename'] + '.sort.md.metrics.txt'
    inFile = open(dupMetFileName,'r')
    lines = inFile.readlines()
    inFile.close()

    metLine = lines[7].rstrip().split()
    dupF = metLine[8]

    outFile.write('fractionDup\t%s\n' % dupF)

    # now get read alignment stats

    inFile = open(myData['alignMetricsFileName'],'r')
    lines = inFile.readlines()
    pairLine = lines[9].rstrip().split()
    inFile.close()

    outFile.write('totalPairedReads\t%s\n' % pairLine[1])
    outFile.write('fractionAligned\t%s\n' % pairLine[6])
    outFile.write('mismatchRate\t%s\n' % pairLine[12])
    outFile.write('indelRate\t%s\n' % pairLine[14])
    outFile.write('meanReadLen\t%s\n' % pairLine[15])
    outFile.write('fractionImproperPairs\t%s\n' % pairLine[19])
    outFile.write('fractionChimera\t%s\n' % pairLine[22])

    # now get insert len stats
    inFile = open(myData['sizeMetricsFileName'],'r')
    lines = inFile.readlines()
    inFile.close()

    statLine = lines[7].rstrip().split()

    numPairsFirstLine= int(statLine[7])
    totPairs = numPairsFirstLine
    i = 8
    while lines[i] != '\n':
       nl = lines[i].rstrip().split()
       totPairs += int(nl[7])
       i += 1
    fractionPairsAssigned = numPairsFirstLine/totPairs



    # check to see if


    outFile.write('pairOrientation\t%s\n' % statLine[8])
    outFile.write('fractionWithPairOrientation\t%.4f\n' % fractionPairsAssigned )
    outFile.write('meanInsertLen\t%s\n' % statLine[5])
    outFile.write('stdInsertLen\t%s\n' % statLine[6])
    outFile.write('medianInsertLen\t%s\n' % statLine[0])
    outFile.write('madInsertLen\t%s\n' % statLine[2])
    outFile.close()



    print('Summary written to',myData['statsSummary'])
######
#
# Make sure programs are available
calc_effective_depth(myData)
summarize_stats(myData)
