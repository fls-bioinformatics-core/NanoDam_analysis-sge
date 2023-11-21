#!/usr/local/bin/python3

import argparse
import os
import sys
import re
import shlex
import time
import shutil
import subprocess
import pandas as pd
import numpy as np
import pybedtools
from difflib import SequenceMatcher

FDRs=(
    2000, 1900, 1800, 1700, 1600, 1500, 1400, 1300, \
    1200, 1100, 1000, 900, 800, 700, 600, 500, \
    475, 450, 425, 400, 375, 350, 325, 300, \
    275, 250, 225, 200, 175, 150, 125, 100, \
    75, 50, 25, 10, 5, 3, 2, 1, 0
    )

##-----------------##
##----Arguments----##
##-----------------##

def parse_args():
    '''Generate parser & define arguments'''
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-r", "--repos",
        nargs = '*',
        type = str,
        required = True,
        help = "List of repositories (i.e., directories)."
        )
    parser.add_argument(
        "-o", "--out",
        type = str,
        required = True,
        help = "Directory for output."
        )

    arguments = parser.parse_args()
    return arguments

##-----------------##
##----Functions----##
##-----------------##

def checkSl(dirLS, regex, verbose=False):
    '''
    Check presence of files with regex
    containting names in parsed set of dirs.
    '''

    # Create a dictionary where keys are the
    # input directory paths
    fir = dict()
    for fIN in dirLS:
        fir[fIN] = False

    sys.stdout.write('\tChecking for files.\n')
    rist = True
    while rist == True:
        for slIN in [k for k,v in fir.items() if v == False]:
            if verbose:
                print("checkSl: looking in %s:" % slIN)
                for f in os.listdir(slIN):
                    print("checkSl: -- %s" % f)
            fir[slIN] = any(
                r == True \
                    for r in [True\
                            if re.compile(regex).search(f) \
                            else False \
                                for f in os.listdir(slIN)\
                            ]\
                )
            if not all(value == True for value in fir.values()):
                time.sleep(5)
                continue
            else:
                rist = False
                sys.stdout.write('\tAll files present.\n')
                break

def evalDir(path):
    try:
        os.makedirs(path)
        return(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def createDir(ori,out,suf,files):
    '''Create dir & copy '*.broadPeak'-files.'''

    dirName = evalDir(os.path.join(ori, str(out + suf)))
    sys.stdout.write("\n>Copy '*.broadPeak'-files to " + dirName + '\n')
    for f in files:
        sys.stdout.write('\t' + os.path.basename(f) + '\n')
        shutil.copy2(f, dirName)

    return(dirName)

def renamer(curDIR):
    '''Rename macs2-derived sge files.'''

    ##Identify_sge-file
    ##-------------------
    try:
        sls = [ \
            f for f in os.listdir() \
            if re.compile('^sge-.*\.o.*').search(f) \
            ]
    except IndexError:
        sys.stdout.write("\t'sge-*.o*'-file already renamed.\n")

    ##Rename_sge-file
    ##-----------------
    date = time.strftime("%Y%m%d", time.localtime())
    for sl in sls:
        os.rename(
            str(curDIR + '/' + sl),
            str( \
                curDIR + '/' + \
                date + '_' + \
                re.sub('(\..*)$', '', sl) + \
                '.log'
                )
            )

    ##Identify_'*.broadPeak'-file
    ##---------------------------
    bP = os.path.abspath(
        [\
            f for f in os.listdir() \
            if re.compile('^(?=(.*-vs-)).*\.broadPeak$').search(f)\
        ][0]\
    )
    ##Identify_DamOnly_'*.broadPeak'-file
    ##-----------------------------------
    bPDN = os.path.abspath(
        [\
            f for f in os.listdir() \
            if re.compile('^(?!(.*-vs-)).*\.broadPeak$').search(f)\
        ][0]
    )
    return(bP,bPDN)

def reader(file, cols = [
        'chr',
        'start',
        'end',
        'pkID',
        'dis',
        'nd',
        'fc',
        'neglog10pval',
        'neglog10qval'
        ]
    ):

    types = {
        'chr': str,
        'start': int,
        'end': int,
        'pkID': str,
        'dis': int,
        'nd': str,
        'fc': float,
        'neglog10pval': float,
        'neglog10qval': float
        }

    types={k:types[k] for k in cols if k in types.keys()}

    try:
        df = pd.read_csv(
            file,
            names = cols,
            header = None,
            dtype = types,
            delim_whitespace = True
            )
    except pd.io.common.EmptyDataError:
        sys.exit("\nEmpty file:\t" + file + "\n")

    return(df)

def populater(ori,dir):

    os.chdir(dir)
    bPs = [f for f in os.listdir() if re.compile('^.*\.broadPeak').search(f)]

    for el in bPs:
        sys.stdout.write("\t" + el + "\n")
        df = reader(el)
        nam = re.sub('^(.*)\/(.*?)\.(.*)$', r'\2', el)
        for FDR in FDRs:
            rP = dir + "/" + str(FDR) + ".regionPeak"

            newDF = (
                df
                .query('neglog10qval >= @FDR')
                .iloc[:,0:3]
                .assign(sample = lambda x: nam)
                .copy()
                )
            newDF['chr'].replace(
                to_replace='^chr(.*)',
                value=r'\1',
                regex=True,
                inplace=True
                )

            try:
                with open(rP, 'a') as curFile:
                    newDF.to_csv(
                        path_or_buf = curFile,
                        sep = '\t',
                        header = False,
                        index = False
                    )
            except IOError as e:
                print("Error: cannot open file")
                if e.errno == errno.EACCES:
                    print("\tPermission denied.")
                    print("\tError message: {0}".format(e))
                    sys.exit()
                # Not a permission error.
                print("\tDoes file exist?")
                print("\tError message: {0}".format(e))
                sys.exit()

    os.chdir(ori)

def colorize(row, cut):
    '''Helper function for color assignment.'''
    if row["pkID"] > cut:
        return("48,8,177")
    else:
        return("213,24,14")

def writer(fdr,df,out,track):
    '''Write out current *.bed file.'''

    try:
        with open(out, 'w') as curFile:
            if track:
                curFile.write(
                    'track name="' + str(fdr) + \
                    '" description="' + str(fdr) + \
                    '" visibility=2 itemRgb="On"\n'
                )
            df.to_csv(
                path_or_buf = curFile,
                sep = '\t',
                header = False,
                index = False
            )
    except IOError as e:
        print("Error: cannot open file")
        if e.errno == errno.EACCES:
            print("\tPermission denied.")
            print("\tError message: {0}".format(e))
            sys.exit()
        # Not a permission error.
        print("\tDoes file exist?")
        print("\tError message: {0}".format(e))
        sys.exit()

def sorter(ori,dir):
    '''Sort populated *.regionPeak files.'''

    os.chdir(dir)

    for FDR in FDRs:
        rP = dir + "/" + str(FDR) + ".regionPeak"
        sys.stdout.write('\t' + str(FDR) + '.regionPeak\n')

        if os.path.getsize(rP) > 0:
            df = reader(rP, cols=['chr', 'start', 'end', 'pkID'])
        else:
            sys.stdout.write('\tEmpty:\t' + rP + '\n')
            continue

        # sys.stdout.write('\tSort:\t' + rP + '\n')
        df.sort_values(
            by = ['chr', 'start'],
            ascending = [1,1],
            axis = 0,
            inplace = True
            )

        writer(FDR,df,rP,False)

    os.chdir(ori)

def merger(ori,dir,id):
    '''Merge overlapping peaks in '*.regionPeak' files.'''

    os.chdir(dir)

    aFS = len([f for f in os.listdir() if re.compile('^.*\.broadPeak').search(f)])

    for FDR in FDRs:
        rP = dir + "/" + str(FDR) + ".regionPeak"
        mP = dir + "/" + str(FDR) + ".mergePeak"
        rpoP = dir + "/" + str(FDR) + ".reproPeak"
        sys.stdout.write('\t' + str(FDR) + '.regionPeak\n')

        if os.path.getsize(rP) > 0:
            df = reader(rP, cols=['chr', 'start', 'end', 'pkID'])
        else:
            sys.stdout.write('\tEmpty:\t' + rP + '\n')
            continue

        bdO = pybedtools.BedTool.from_dataframe(df)
        mbdO = bdO.merge(
                c = 4,
                o = 'count_distinct'
                )

        pdDF = pd.read_table(
            mbdO.fn,
            names=['chr', 'start', 'end','pkID'],
            dtype={
                'chr': str,
                'start': int,
                'end': int,
                'pkID': int
                }
            )

        ##''*.mergePeak'-file
        ##-------------------
        mergDF = (
            pdDF
            .assign(rep = lambda x: np.round_(((x.pkID/aFS)*100), decimals=2))
            .assign(score = lambda x: "0")
            .assign(strand = lambda x: ".")
            .assign(thickStart = lambda x: x.start)
            .assign(thickEnd = lambda x: x.end)
            .assign(rgb = pdDF.apply(colorize, axis=1, args = (int(aFS/2),)))
            .iloc[:, np.r_[0:3,4:10]]
            )

        writer(FDR,mergDF,mP,True)

        ##''*.reproPeak'-file
        ##-------------------
        rpoDF = (
            mergDF
            .loc[mergDF['rep'] > 50]
            .iloc[:, np.r_[0:3]]
            .assign(name = lambda x: str(id))
            )

        writer(FDR,rpoDF,rpoP,True)

    os.chdir(ori)

##---------------------##
##----Main_workflow----##
##---------------------##

def main():
    args = parse_args()
    oriDIR = os.getcwd()

    ##Check_presence_of_both_'.*\.broadPeak'-files_in_all_dirs
    ##--------------------------------------------------------
    ##Note:Alternative_is_to_search_for_slurm-file_by_jobIDs
    sys.stdout.write("\n>Checking presence of '*.broadPeak'-files\n")
    checkSl(args.repos, '^(?=(.*-vs-)).*\.broadPeak$')
    sys.stdout.write("\n>Checking presence of DamOnly-'*.broadPeak'-files\n")
    checkSl(args.repos, '^(?!(.*-vs-)).*\.broadPeak$', verbose=True)

    ##Rename_slurm_files
    ##------------------
    sys.stdout.write("\n>Rename files\n")
    BPs = list()
    damOBPs = list()
    for el in args.repos:

        absDIR = os.path.abspath(el)
        os.chdir(absDIR)
        sys.stdout.write('\t' + absDIR + '\n')

        BP,damOBP = renamer(absDIR)

        BPs.append(BP)
        damOBPs.append(damOBP)

        os.chdir(oriDIR)

    ##Deduplicate_damONs
    ##------------------
    subDic = dict()
    subSet = set()
    for el in damOBPs:
        subDic[el] = os.path.basename(el)
        subSet.add(os.path.basename(el))

    subDONs = list()
    for setEL in subSet:
        subDONs.append([k for k,v in subDic.items() if v == setEL][0])

    ##Create_dirs_&_copy_bedgraph-files
    ##---------------------------------
    damOBPDIR = createDir(oriDIR,args.out,"_DamOnly_peaks",subDONs)
    BPDIR = createDir(oriDIR,args.out,"_peaks", BPs)

    ##Populate_'FDR.regionPeak'-files
    ##-------------------------------
    sys.stdout.write("\n>Read in '*.broadPeak'-files\n")
    populater(oriDIR,BPDIR)
    sys.stdout.write("\n>Read in DamOnly '*.broadPeak'-files\n")
    populater(oriDIR,damOBPDIR)

    ##Sort_'{FDR}.regionPeak'-files
    ##-----------------------------
    sys.stdout.write("\n>Sort '*.regionPeak'-files.\n")
    sorter(oriDIR,BPDIR)
    sys.stdout.write("\n>Sort DamOnly '*.regionPeak'-files\n")
    sorter(oriDIR,damOBPDIR)

    ##Merge_'{FDR}.regionPeak'-files
    ##------------------------------
    sys.stdout.write("\n>Merge '*.regionPeak'-files.\n")
    merger(oriDIR,BPDIR,args.out)
    sys.stdout.write("\n>Merge DamOnly '*.regionPeak'-files\n")
    merger(oriDIR,damOBPDIR,args.out)

    ##Remove_'{FDR}.regionPeak'-files
    ##-------------------------------
    sys.stdout.write("\n>Remove '*.regionPeak'-files\n")
    for dir in [BPDIR,damOBPDIR]:
        for FDR in FDRs:
            rP = dir + "/" + str(FDR) + ".regionPeak"
            os.remove(rP)

    sys.stdout.write('\nAll done.\n')

if __name__ == '__main__':
    main()
