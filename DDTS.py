#!/usr/bin/python
'''
Identifies open chromatin regions (DNase Hypersensitive Sites, DHS).
It may be run using either bam files or the processed wig files.
Example run:

#Can be run on variable number of replicates. Whether you are running on 1 or 3 replicates.
'''

from scipy import stats
import sys
import numpy as np
import statsmodels.api as sm
import smtplib
from optparse import OptionParser
import string
import os
import os.path
from multiprocessing import Process
import multiprocessing
import subprocess


def readLines(file, window):
    """
    Reads in wig files ~100,000 bp at a time to keep memory low in arrays.

    :param file: file name
    :param window:  window size (in bp)
    
    :return: an array of part of the wig file, ~100,000bp long.
    """

    x = 100000 / window + 1
    end = window * x
    readLine = []
    for f in range(0, end):
        line = file.readline()
        if line == '':
            break  # Break if line is empty, end of file reached
        columns = line.split("\t")
        columns[-1] = columns[-1].replace("\n", "")  # remove new line character
        readLine.append(columns)
    return readLine


def readLinesSize(file, size):
    """
    Reads in a wig file up to a customized point depending on the respective wig files. Will read to a specified bp position 
    This is for the first read through each wig file, so that we scan through each wig file at the same bp position
   
    :param file: file name
    :param size: the bp position to read up to, or the size of the array to be made

    :return : returns an array of the wig file 
    """

    readLine = []
    for f in range(0, size):
        line = file.readline()
        if line == '':
            break  # Break if line is empty, end of file reached
        columns = line.split("\t")
        columns[-1] = columns[-1].replace("\n", "")  # remove new line character
        readLine.append(columns)
    return readLine


def findOp(ctrlList, expmList, zzcut, outName, distance, wind, shif):
    """
    Identifies DHS sites between control and a digested samples using a specified Z value cut of and a specified distance between control and digested.

    :param ctrlList: An array containing the names of the control files
    :param expmList: An array containing the names of the experimental or digested files
    :param zzcut: The z value cut off for the statistical test
    :param outName: The output file name
    :param distance: The cut off distance between the control and experimental. The mean control value is divided by the experimental value and to be an open region has to be greater than this value.

    :return : Writes a bed file to the current directory containing the DHS sites. 
    """

    window = wind  # Window size to scan around
    shift = shif  # Shift size to shift window to next location

    # Open each of the files
    for i in range(len(ctrlList)):
        print(ctrlList[i], expmList[i])
        ctrlList[i] = open(ctrlList[i], "r")
        expmList[i] = open(expmList[i], "r")
    output = outName
    zCut = zzcut
    dist = distance

    # For use later, endFile=the results BED file. statArr =holds any stats formed
    endFile = []

    # Reading info section of each file
    infoCtrlFile = []
    infoExpmFile = []
    for i in range(len(ctrlList)):
        infoCtrlFile.append(ctrlList[i].readline().split("\t"))
        infoExpmFile.append(expmList[i].readline().split("\t"))

    # Chromosome number is grabbed from one of the files for creation of bed file
    chrmNum = infoCtrlFile[0][0].split()[1].split("=")[1]

    # Grab start bp position of each wig file (it varies from wig to wig)
    startCtrlList = []
    startExpmList = []
    for i in range(len(ctrlList)):
        startCtrlList.append(int(infoCtrlFile[i][0].split()[2].split("=")[1]))
        startExpmList.append(int(infoExpmFile[i][0].split()[2].split("=")[1]))

    # Find the max start position to start the other files at.
    maxArray = [np.max(startCtrlList), np.max(startExpmList)]
    maxVal = np.max(maxArray)

    # Read in first bit of each file and put into an array
    x = 100000 / window + 1
    end = window * x

    # Read in first bit of each file up to the same end bp position (end - start accomplishes this)
    readCtrlArray = []
    readExpmArray = []
    for i in range(len(ctrlList)):
        readCtrlArray.append(readLinesSize(ctrlList[i], end - startCtrlList[i]))
        readExpmArray.append(readLinesSize(expmList[i], end - startExpmList[i]))

    # Slice first chunk to the starting position. so we are scanning through each file at the same start bp
    for i in range(len(readCtrlArray)):
        readCtrlArray[i] = readCtrlArray[i][(maxVal - startCtrlList[i]) + 1:]
        readExpmArray[i] = readExpmArray[i][(maxVal - startExpmList[i]) + 1:]


    # Cut all to multiple of window, for looping correctly to end of current slice
    slice = len(readCtrlArray[0]) - (len(readCtrlArray[0]) / window * window)
    for i in range(len(readCtrlArray)):
        readCtrlArray[i] = readCtrlArray[i][slice:]
        readExpmArray[i] = readExpmArray[i][slice:]

    # The overall start position in bp.
    start = int(maxVal)
    curPos = start + slice
    flag = True

    # loop until end of file is met
    while flag == True:
        # Identifying the file with the smallest length (the farthest bp we can scan)
        minLength = []
        for i in range(len(readCtrlArray)):
            minLength.append(len(readCtrlArray[i]))
            minLength.append(len(readExpmArray[i]))
        minLen = np.min(minLength)
        # If a file length is less than 100000 then turn flag to false as we have reached end of one of the files.
        if minLen < 100000:
            if curPos > 10000:
                flag = False

        # Set up arrays to store each replicate values for each bp in the window
        # As well as a running total to calculate average
        totalCtrl = []
        totalExpm = []
        totValCtrl = []
        totValExpm = []

        # Begin looping through the files one window at a time
        for i in range(0, ((minLen - window) / shift) + 1):  # shift
            totalCtrl = []  # Set up arrays to store each replicate values for each bp in the window
            totalExpm = []  # As well as a running total to calculate average

            # Create arrays and running totals for each separate replicate
            for h in range(len(readCtrlArray)):
                totalCtrl.append([])
                totalExpm.append([])

            # Loop through window of each file and grab x*window number of values and put into separate arrays.
            for j in range(0, window):  # window
                for k in range(len(readCtrlArray)):
                    tempCVal = float(readCtrlArray[k][i * shift + j][0])
                    tempEVal = float(readExpmArray[k][i * shift + j][0])
                    if tempCVal == 0:
                        tempCVal = 0.0001  # 0.0001
                    if tempEVal == 0:
                        tempEVal = 0.0001
                    totalCtrl[k].append(tempCVal)
                    totalExpm[k].append(tempEVal)

            # Grab means for each replicate
            controlArray = []
            expmermArray = []
            for k in range(len(totalCtrl)):
                controlArray.append(np.mean(totalCtrl[k]))
                expmermArray.append(np.mean(totalExpm[k]))

            # Perform statistical test between each ctrl and expm means as long as there is atleast 3 replicates
            if len(controlArray) >= 3:
                stats = sm.stats.ttest_ind(controlArray, expmermArray, alternative='larger', usevar='unequal')
                pVal = np.mean(stats[1])
                zVal = np.mean(stats[0])
            else:  # if < 3 replicates, can't do t-test between means, do t-test versus each window of each replicate.
                statArray = []
                for k in range(len(totalCtrl)):
                    statArray.append(sm.stats.ztest(totalCtrl[k], totalExpm[k]))
                pValArray = []
                zValArray = []
                for k in range(len(statArray)):
                    pValArray.append(statArray[k][1])
                    zValArray.append(statArray[k][0])
                pVal = np.mean(pValArray)
                zVal = np.mean(zValArray)

            # statistical values to be above a cut off.
            if zVal > zCut:
                flag2 = True
                for k in range(len(controlArray)):  # Require all ctrl file means > expm mean values by a certain amount
                    if (float(controlArray[k]) / expmermArray[k] < dist):
                        flag2 = False
                # If previous is true, create window entry and store in the output.
                if flag2 == True:
                    newFile = []
                    newFile.append(str(chrmNum))
                    newFile.append(str(curPos))
                    newFile.append(str(curPos + window))
                    newFile.append(str(zVal))
                    newFile.append(str(pVal))
                    for k in range(len(controlArray)):
                        newFile.append(str(controlArray[k]))
                    for k in range(len(expmermArray)):
                        newFile.append(str(expmermArray[k]))
                    endFile.append(newFile)
            curPos = curPos + shift

        # Slice next 100 000, if less than 100 000 then finish last set and break.
        for k in range(len(ctrlList)):
            readCtrlArray[k] = readLines(ctrlList[k], window)
            readExpmArray[k] = readLines(expmList[k], window)

        curPos = curPos + window - shift

    print('Number of Sites: ', len(endFile))

    # Write output to a file
    with open(output, "w") as file:
        file.writelines('\t'.join(i) + '\n' for i in endFile)
    return 0


def vararg_callback(option, opt_str, value, parser):
    """
    Helper function for OptionParser to grab multiple files with one flag instance 
    For use with the -c and -e flags, allow option parser to have multiple files after one flag instance

    :param option: option instance that is calling the callback
    :param opt_str: the option string, eg -c 
    :param value: the arguments to this option, eg. file1, file2
    :param parser: the OptionParser instance

    :return : None
    """
    assert value is None
    value = []

    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False

    for arg in parser.rargs:
        if arg[:2] == "--" and len(arg) > 2:  # stop on --check like options
            break
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):  # stop on -d, but not on -6 or -5.0
            break
        value.append(arg)

    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)


if __name__ == '__main__':
    """
    The main function that drives this whole program
    Specifies the flags you can use.
    Will convert bam to bed files and run F-seq if that is the specified route.
    Or will run findOp on the wig files to identify the open regions.
    #can be run on a variable number of replicates.
    """

    # Add flag options
    usage = "If BAM files run \n usage: %prog -c [undig1.bam] [undig2.bam] [undig3.bam] -e [expm1.bam] [expm2.bam] [expm3.bam] \n This will convert BAM files to BED file and run FSEQ on each BED file to report WIG files. Once that finishes run the next command on each chromosome wig file \n usage: %prog -c [undig1_1.wig] [undig2_1.wig] [undig3_1.wig] -e [expm1_1.wig] [expm2_1.wig] [expm3_1.wig] -z [zcut] -d [foldcut] -o [output] -n [chromosomes] -w"
    parser = OptionParser(usage=usage)
    parser.add_option('-c', '--ctrl', dest='ctrl', action="callback", callback=vararg_callback,
                      help='Control/Background DataSets')
    parser.add_option('-e', '--expm', dest='expm', action="callback", callback=vararg_callback,
                      help='Experimental/Digested Datasets')
    parser.add_option('-z', '--zcut', dest='zcut', type="int",
                      help='Cut off value for the z-test', default=3)
    parser.add_option('-o', '--output', dest='output',
                      help='Output file name')
    parser.add_option('-d', '--dist', dest='dist', type="float",
                      help='Fold change cut-off between undigested and digested', default=1.5)
    parser.add_option('-l', '--flength', dest='flength', type="int",
                      help='Fseq feature length', default=200)
    parser.add_option('-w', '--wig', action="store_true", dest="wig",
                      help='If you\'ve already run F-Seq or have wig files. Each chromosome must be run separately in this mode')
    parser.add_option('-W', '--window', dest='window', type="int",
                      help='Length of scan window (bp)', default=300)
    parser.add_option('-s', '--step', dest='step', type="int",
                      help='Length of step increase (bp)', default=25)
    options, args = parser.parse_args()
    if not options.ctrl:
        parser.error('Must provide control files')
    if not options.expm:
        parser.error('Must provide experimental files')
    if not options.wig:  # If wig option not selected, convert bam to bed, run Fseq, and find open regions

        # check program dependencies
        try:
            subprocess.call('type fseq', shell=True)  # use unix type builtin to check path
        except EnvironmentError:
            print('Error: dependecy \'fseq\' not found in your path (https://github.com/aboyle/F-seq). Stopping.')
            system.exit(1)
        try:
            subprocess.call('type bedtools', shell=True)
        except EnvironmentError:
            print('Error: dependancy \'bedtools\' not found in your path (https://github.com/arq5x/bedtools2). Stopping.')
            system.exit(1)

        # Convert all control Bam files to bed files
        fileC = []
        cmd = 'mkdir run'
        os.system(cmd)
        processes = []
        print('Converting control bam files to bed files.')
        for f in options.ctrl:
            fileTag = os.path.splitext(f)[0]
            fileC.append(fileTag)
            cmd = 'bedtools bamtobed -cigar -i {} > {}.bed'.format(f, fileTag)
            print(cmd)
            processes.append(subprocess.Popen(cmd, shell=True))

        # Convert all experimental Bam files to bed files
        fileE = []
        print('Converting experimental bam files to bed files.')
        for f in options.expm:
            fileTag = os.path.splitext(f)[0]
            fileE.append(fileTag)
            cmd = 'bedtools bamtobed -cigar -i {} > {}.bed'.format(f, fileTag)
            print(cmd)
            processes.append(subprocess.Popen(cmd, shell=True))

        # Wait for bamtobed to finish before running F-seq
        for p in processes:
            p.wait()

        # Run Fseq on all given control bed files to produce wig outputs
        processes = []
        print('Running fseq on control files.')
        for f in fileC:
            cmd = 'mkdir {}'.format(f)
            os.system(cmd)
            os.chdir(f)
            cmd = 'fseq -v -l {} ../{}.bed > {}_fseq.log'.format(options.flength, f, f)
            print(cmd)
            processes.append(subprocess.Popen(cmd, shell=True))
            os.chdir("..")

        # Run Fseq on all given experimental bed files to produce wig outputs
        print('Running fseq on experimental files.')
        for f in fileE:
            cmd = 'mkdir {}'.format(f)
            os.system(cmd)
            os.chdir(f)
            cmd = 'fseq -v -l {} ../{}.bed > {}_fseq.log'.format(options.flength, f, f)
            print(cmd)
            processes.append(subprocess.Popen(cmd, shell=True))
            os.chdir("..")

        # wait for f-seq to finish before moving files to desired location
        for p in processes:
            p.wait()

        # Move all wig files from each control replicate to combined folder adn label correctly
        for f in fileC:
            os.chdir(f)
            cmd = 'for filename in *.wig; do mv "$filename" "{}_$filename"; done'.format(f)
            os.system(cmd)
            os.chdir("../run")
            cmd = 'ln -s ../{}/*.wig ./'.format(f)
            os.system(cmd)
            os.chdir("..")

        # Move all wig files from each experimental replicate to combined folder adn label correctly
        for f in fileE:
            os.chdir(f)
            cmd = 'for filename in *.wig; do mv "$filename" "{}_$filename"; done'.format(f)
            os.system(cmd)
            os.chdir("../run")
            cmd = 'ln -s ../{}/*.wig ./'.format(f)
            os.system(cmd)
            os.chdir("..")

        os.chdir("run")

    else:  # if the wig file is specified, then we have wig files, Can identify open regions.
        if not options.output:
            parser.error('Must provide an output file name')

        print('Identifying open regions.')
        findOp(options.ctrl, options.expm, options.zcut, options.output, options.dist, options.window, options.step)
