#!/usr/local/bin/python3

import argparse
import os
import sys
import re
import shlex
import time
import shutil
import subprocess
from difflib import SequenceMatcher

shItr = 1
tmpl = """\
#!/bin/bash
#!
#! Name of the job:
#$ -N sge-{name}
#$ -V
#$ -j y
#$ -pe smp.pe {cores}
#$ -cwd

echo -e "JobID: $JOBID\\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
echo -e "\\nExecuting command:\\n==================\\n{SGE_CMD}\\n"

eval {SGE_CMD}
"""

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
        "-d", "--defaults",
        type = str,
        default = "dm6",
        help = "Load defaults for species of interest."
        )
    parser.add_argument(
        "-f", "--feedback",
        type = str,
        required = False,
        default="noemail",
        help = "Complete mail address to receive slurm feedback."
        )
    parser.add_argument(
        "-o", "--out",
        type = str,
        required = True,
        help = "Directory for output."
        )
    parser.add_argument(
        "-p", "--exppre",
        type = str,
        required = True,
        help = "Common string in experimental samples."
        )
    parser.add_argument(
        "-c", "--ctrlpre",
        type = str,
        required = True,
        help = "Common string in control Dam-samples."
        )
    parser.add_argument(
        "-l", "--chrSize",
        type = str,
        default = None,
        help = "List of chromsome sizes."
        )
    parser.add_argument(
        "-m", "--macs2",
        type = str,
        default = "/usr/bin/macs2",
        help = "Path to 'MACS2'."
        )
    parser.add_argument(
        "-n", "--quantile",
        type = str,
        default = "./quantile_norm_bedgraph.pl",
        help = "Path to 'quantile_norm_bedgraph.pl'."
        )
    parser.add_argument(
        "-a", "--average",
        type = str,
        default = "./average_tracks.pl",
        help = "Path to 'average_tracks.pl'."
        )
    parser.add_argument(
        "-b", "--bgToBw",
        type = str,
        default = "/usr/bin/bedGraphToBigWig",
        help = "Path to 'bedGraphToBigWig'."
    )

    arguments = parser.parse_args()
    return arguments

##-----------------##
##----Functions----##
##-----------------##

def checkt(toolPath):
    '''Checking paths of used tools.'''
    tool = os.path.basename(toolPath)
    #logging.info('Checking: ' + tool)
    sys.stdout.write('\tChecking: '+ tool + '\n')

    if os.path.exists(toolPath):
        tcho = toolPath
    else:
        tcho = ""

    if shutil.which(tool):
        tdet = shutil.which(tool)
    else:
        tdet = ""

    if tcho == tdet and tdet != "" and tcho != "":
        ##Chosen_&_detected_dirs_exist_and_are_similar

        #logging.info(tool + ' checked: ' + tdet)
        tuse = tdet
        return tuse

    elif tcho != tdet and tdet != "" and tcho != "":
        ##Chosen_&_detected_dirs_exist_but_are_not_similar

        #logging.warning('"' + tcho + '" not default installation: "' + tdet + '"')
        sys.stdout.write(tcho + '" not default installation: "' + tdet + '"\n')

        while True:
            desc = input('Use: [1 - ' + tcho + '] or [2 - ' + tdet + ']?\n')
            if desc == "1" or desc == tcho:
                tuse = tcho
                break
            elif desc == "2" or desc == tdet:
                tuse = tdet
                break
            else:
                sys.stdout.write('\nNot a valid choice.\n')
                continue

        #logging.info(tool + ' used: ' + tuse)
        sys.stdout.write(tool + ' used: ' + tuse + '\n')
        return tuse

    elif tcho != tdet and tdet != "" and tcho == "":
        ##Detected_but_not_chosen_dir_exist
        #logging.warning('"' + toolPath + \
        #   '" not default installation: "' + tdet + '"')
        sys.stdout.write(toolPath + \
            '" not default installation: "' + tdet + '"\n')

        #logging.warning('Not found: ' + toolPath)
        tuse = tdet
        #logging.info(tool + ' used: ' + tuse)
        sys.stdout.write(tool + ' used: ' + tuse + '\n')
        return tuse

    elif tcho != tdet and tdet == "" and tcho != "":
        ##Chosen_but_not_detected_dir_exist
        tuse = tcho
        #logging.info(tool + ' used: ' + tuse)
        sys.stdout.write(tool + ' used: ' + tuse + '\n')
        return tuse

    elif tcho =="" and tdet == "":
        ##Neither_chosen_nor_detected_dir_exists
        #logging.error('Not found: ' + toolPath)
        #logging.error('Not found: ' + tool)
        sys.exit('Error: Not found: ' + tool + \
            '\nUse python3 bow2map.py --help.\n')

def checkSl(dirLS, regex):
    '''
    Check presence of files with regex
    containting names in parsed set of dirs.
    '''

    fir = dict()
    for fIN in dirLS:
        fir[fIN] = False

    sys.stdout.write('\tWaiting for cluster.\n')
    rist = True
    while rist == True:
        for slIN in [k for k,v in fir.items() if v == False]:
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

def checkFin(jobIDs):
    '''Check all provided jobIDs are no longer registered by sge.'''

    sys.stdout.write('\tWaiting for cluster.\n')
    rest = True
    interval = 30
    while rest == True:
        print("checkFin: getting output from 'qstat'")
        chk = subprocess.check_output(['qstat'])
        qstat_list = []
        for line in chk.decode("utf-8").split('\n'):
            print("QSTAT> %s" % line.rstrip())
            try:
                if line.split()[0].isdigit():
                    qstat_list.append(line.split()[0])
                    print("...... extracted job ID '%s'" % line.split()[0])
            except IndexError:
                pass
        chk = qstat_list
        print("checkFin: job list = %s" % chk)
        if not set(jobIDs).intersection(chk):
            rest = False
            break
        else:
            print("checkFin: wait for %ss..." % interval)
            time.sleep(interval)
            continue
    sys.stdout.write('\tJob(s) finished.\n')

def checkQue(jobIDs):
    '''Check if jobs are registered by sge.'''

    sys.stdout.write('\tWaiting for cluster.\n')
    print("checkQue: waiting for jobs: %s" % jobIDs)
    print("checkQue: getting output from 'qstat'")
    chk = subprocess.check_output(['qstat'])
    qstat_list = []
    for line in chk.decode("utf-8").split('\n'):
        print("QSTAT> %s" % line.rstrip())
        try:
            if line.split()[0].isdigit():
                qstat_list.append(line.split()[0])
                print("...... extracted job ID '%s'" % line.split()[0])
        except IndexError:
            pass
    chk = qstat_list
    print("checkQue: job list = %s" % chk)
    if set(jobIDs).issubset(chk):
        sys.stdout.write('\tJob(s) running.\n')
    else:
        sys.exit("\nOne or more job(s) not running.\n")

def evalDir(path):
    try:
        os.makedirs(path)
        return(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def createDir(ori,out,suf,files):
    '''Create dir & copy '*.bedgraph'-files.'''

    dirName = evalDir(os.path.join(ori, str(out + suf)))
    sys.stdout.write("\n>Copy '*.bedgraph'-files to " + dirName + '\n')
    for f in files:
        sys.stdout.write('\t' + os.path.basename(f) + '\n')
        shutil.copy2(f, dirName)

    return(dirName)

def create_sh(cmd,mailAc,cores=8):
    '''OBS! 'dir' as in 'damMer.py' changed to 'os.getcwd()'.'''

    global shItr
    cmdName = re.compile('\..*').sub('', os.path.basename(cmd.split(" ")[0]))
    #sys.stdout.write('\tcmdName:\t' + cmdName + "\n")
    fileName = os.getcwd() + "/" + str(shItr) + "_" + cmdName + ".sh"
    with open(fileName, 'w') as shOUT:
        shOUT.write(tmpl.format(name=cmdName, SGE_CMD=cmd, cores=cores, mail=mailAc))
    shItr += 1

    print("create_sh: created submission script '%s'" %
                     fileName)
    with open(fileName, 'rt') as shOUT:
        for line in shOUT:
            print("> %s" % line.rstrip('\n'))
    return(fileName)

def submit(cmdSH, dpdIDs=''):
    '''dpdIDs should have format: '<jobID1>,<jobID2>,etc.' '''

    Sub = "qsub -terse"
    if dpdIDs:
        Sub += " -hold_jid " + dpdIDs
    Sub += " " + cmdSH
    print("submit: created command '%s'" % cmdSH)
    #sys.stdout.write("\tSubmit:\t" + str(Sub) + "\n")
    try:
        prc = subprocess.Popen(
            shlex.split(Sub),
            shell = False,
            stderr = subprocess.PIPE,
            stdout = subprocess.PIPE
            )
    except Exception as e:
        sys.exit("\nERROR: 'devscript_damMer.py' aborted:\t" + type(e).__name__ + "\n")

    submit_stderr = str(prc.communicate()[1].decode("utf-8"))
    print("submit: stderr from qsub: %s" % submit_stderr)
    jobID = str(prc.communicate()[0].decode("utf-8")).rstrip()

    print("submit: returned job ID '%s'" % jobID)
    if not jobID:
        sys.exit("submit: failed to get job id\n")
    return(jobID)

def readlines_reverse(filename):
    '''Retrieve individual lines from file end.'''

    with open(filename) as qfile:
        qfile.seek(os.SEEK_SET, os.SEEK_END)
        position = qfile.tell()
        line = ''
        while position >= 0:
            qfile.seek(position)
            next_char = qfile.read(1)
            if next_char == "\n":
                yield line[::-1]
                line = ''
            else:
                line += next_char
            position -= 1
        yield line[::-1]

def screener(f):
    '''Check for 'All done.' in all 'sge-*'-files.'''

    i = 0
    for qline in readlines_reverse(f):
        if re.compile('.*All\s+done\..*').search(qline):
            return(True)
        else:
            i += 1
            if i <=3:
                continue
            else:
                return(False)

def extractor(path):
    pat = re.compile('.*Reading\sdata\sfiles.*')
    rec = re.compile('.*Using\s(.+?)\sas.*')

    tr = False
    sh = 0
    sams = dict()
    with open(path, 'r') as slIn:
        for lNr, l in enumerate(slIn,1):
            #print(str(l.split(" ")))
            if tr == True and sh in [0,1]:
                sam = l.split('\t')
                sams[sam[0]]=sam[1]
                sh += 1
                continue
            if tr == True and sh == 2:
                dam = rec.search(l).group(1)
                break
            if pat.search(l):
                tr = True

    pres = dict()
    for k,v in sams.items():
        pre = re.compile('^(.+?)(\..*)+').search(os.path.basename(v)).group(1)
        pres[k]=pre

    return(pres,dam)

def renamer(curDIR, ctrlpre, exppre):
    '''Rename output files of damidseq_pipeline_vR.1.'''

    ##Identify_sge-file
    ##-------------------
    sl = [
        f for f in os.listdir() \
        if re.compile('^sge-.*\.o.*', re.IGNORECASE).search(f)
        ][0]
    sys.stdout.write('\tsge file:\t' + sl + '\n')

    ##Extract_info_from_sge-file
    ##----------------------------
    fs, dam = extractor(os.path.join(curDIR,sl))
    # for k,v in fs.items():
    #     if k == dam:
    #         sys.stdout.write('\tDam:\t' + k + '\t' + v + '\n')
    #     else:
    #         sys.stdout.write('\tPro:\t' + k + '\t' + v + '\n')

    ##Rename_'dam-ext300.bam'-file
    ##----------------------------
    damFile = [
        f for f in os.listdir() \
        if re.search(dam, f, re.IGNORECASE) \
            and re.compile('.*-ext300.bam').search(f) \
        ][0]
    #sys.stdout.write('\tdamFile:\t' + damFile + '\n')
    damNew = str(curDIR + '/' + fs[dam] + ".ext300.bam")
    if re.search(ctrlpre, fs[dam], re.IGNORECASE):
        os.rename(
            str(curDIR + '/' + damFile),
            damNew
            )

    ##Rename 'dam-DamOnly.gatc.bedgraph'
    ##----------------------------------
    damOnlyFile = [
        f for f in os.listdir() \
        if re.search(dam, f, re.IGNORECASE) \
            and re.compile('.*-DamOnly.gatc.bedgraph').search(f) \
        ][0]
    #sys.stdout.write('\tdamOnlyFile:\t' + damOnlyFile + '\n')
    damOnlyNew = str(curDIR + '/' + fs[dam] + ".DamOnly.gatc.bedgraph")
    if re.search(ctrlpre, fs[dam], re.IGNORECASE):
        os.rename(
            str(curDIR + '/' + damOnlyFile),
            damOnlyNew
            )

    ##Rename_'exp-ext300.bam'-file
    ##----------------------------
    expFile = [
        f for f in os.listdir() \
        if not re.search(dam,f, re.IGNORECASE) \
            and re.compile('.*-ext300.bam').search(f)
        ][0]
    #sys.stdout.write('\texpFile:\t' + expFile + '\n')
    expFilePre = (
                re
                .compile('^(.*)-ext300(\..*)+', re.IGNORECASE)
                .search(expFile)
                .group(1)
                )
    #sys.stdout.write('\texpFilePre:\t' + expFilePre + '\n')
    expNew = str(curDIR + '/' + fs[expFilePre] + ".ext300.bam")
    if re.search(exppre, fs[expFilePre], re.IGNORECASE):
        os.rename(
            str(curDIR + '/' + expFile),
            expNew
            )

    ##Rename_bedgraph-file
    ##--------------------
    bGF = [
        f for f in os.listdir() \
        if re.compile('^.*-vs-.*\.gatc\.bedgraph').search(f)
        ][0]
    #sys.stdout.write('\tbGF:\t' + bGF + '\n')
    nbGF = str(curDIR + '/' + fs[expFilePre] + '-vs-' + fs[dam] + ".gatc.bedgraph")
    os.rename(
        str(curDIR + '/' + bGF),
        nbGF
        )

    ##Rename_sge-_&_pipeline-file
    ##-----------------------------
    date = time.strftime("%Y%m%d", time.localtime())
    os.rename(
        str(curDIR + '/' + sl),
        str(curDIR + '/' + date + '_' + re.sub('(\..*)$', '', sl) + '.log')
        )
    pi = [
        f for f in os.listdir() \
        if re.compile('^pipeline.*').search(f)
        ][0]
    #sys.stdout.write('\tpi:\t' + pi + '\n')
    os.rename(
        str(curDIR + '/' + pi),
        str(curDIR + '/' + date + '_pipeline_' + re.sub('(\..*)$', '', sl) + '.log')
        )

    return(nbGF, damOnlyNew, damNew, expNew)

def peakCalling(macs2,genomeSize,ctrl,*trt):
    '''Function to call peaks via macs2.'''

    pkc = macs2 + \
        " callpeak" + \
        " --format BAM" + \
        " --gsize " + str(genomeSize) + \
        " --keep-dup all" + \
        " --bw " + str(300) + \
        " --qvalue " + str(0.05) + \
        " --mfold " + str(5) + \
           " " + str(50) + \
        " --broad"

    if trt:
        pkc = pkc + \
        " --treatment " + str(trt[0]) + \
        " --control " + ctrl + \
        " --name " + re.sub('\.ext300\.bam$', '', os.path.basename(str(trt[0]))) + \
        "-vs-" + re.sub('\.ext300\.bam$', '', os.path.basename(ctrl))
    else:
        pkc = pkc + \
        " --treatment " + str(ctrl) + \
        " --name " + re.sub('\.ext300\.bam$', '', os.path.basename(ctrl))

    return(pkc)

def quantNorm(ori,dir,quant,mailAc):
    '''
    Perform Quantile normalization on provided set of *.bedgraph files.
    'quantile_norm_bedgraph.pl'-script erases all trailing 'chr'-indicator.
    '''

    sys.stdout.write("\n>Quantile normalization - '*.gatc.bedgraph' files\n")
    os.chdir(dir)
    fs = [f for f in os.listdir() if re.compile('.*\.gatc\.bedgraph').search(f)]

    qna = "perl" + \
        " " + quant + \
        " " + ' '.join(fs)
    qnaSH = create_sh(qna,mailAc)
    #sys.stdout.write("\t" + qnaSH + "\n")
    qnaID = submit(qnaSH)

    os.chdir(ori)
    return(qnaID)

def average(ori,dir,aver,mailAc):
    '''Average provided *.bedgraph files per GATC fragment.'''

    sys.stdout.write("\n>Averaging - '*.quant.norm.bedgraph' files\n")
    os.chdir(dir)
    qGFs = [f for f in os.listdir() if re.compile('.*quant\.norm\.bedgraph').search(f)]

    avg = "perl" + \
        " " + aver + \
        " " + ' '.join(qGFs)
    avgSH = create_sh(avg,mailAc)
    #sys.stdout.write("\t" + avgSH + "\n")

    avgID = submit(avgSH)

    os.chdir(ori)
    return(avgID)

def bwer(ori,dir,chroms,bGTBW,mailAc):
    '''Convert all '*.quant.norm.*' files into *.bw format.'''

    sys.stdout.write("\n>Convert '*.quant.norm.*'-files into '*.bw'\n")
    os.chdir(dir)
    qnGFs = [f for f in os.listdir() if re.compile('.*\.quant\.norm.*').search(f)]

    jobIDs = list()
    for qnGF in qnGFs:
        bw = bGTBW + \
            " " + qnGF + \
            " " + chroms + \
            " " + qnGF + ".bw"
        bwSH = create_sh(bw,mailAc)
        #sys.stdout.write("\t" + bwSH + "\n")

        bwID = submit(bwSH)
        jobIDs.append(bwID)

    os.chdir(ori)
    return(jobIDs)

##---------------------##
##----Main_workflow----##
##---------------------##

def main():
    args = parse_args()
    oriDIR = os.getcwd()

    ##Tester----------------------------------------------------------------------------
    # for arg in vars(args):
    #     sys.stdout.write('{arg}:\t{value}\n'.format(arg=arg,value=getattr(args,arg)))
    # sys.stdout.write('ctrlpre:\t' + args.ctrlpre + '\texppre:\t' + args.exppre + '\n')
    ##----------------------------------------------------------------------------------

    ##Calculate_genomeSize
    ##--------------------
    if args.chrSize == None:
        inDir = "/mnt/home1/brand/rk565/resources"
        if args.defaults == "dm6":
            args.chrSize = '/'.join(
                [inDir, "dm6.chrom.sizes.mod"]
            )
        elif args.defaults == "mm10":
            args.chrSize = '/'.join(
                [inDir, "mm10.chrom.sizes"]
            )
        else:
            sys.exit('Unsupported species: --defaults=[dm6/mm10].\n')

    with open(args.chrSize, 'r') as inFile:
        genSize = sum(int(l.split()[1]) for l in inFile)
    sys.stdout.write('\n>Checking defaults\n')
    sys.stdout.write(
        '\tSpecies:\t' + args.defaults + '\n'\
        '\tGenomesize:\t' + str(genSize) + '\n'
        )

    ##Checking_executables
    ##--------------------
    sys.stdout.write('\n>Checking executables\n')
    qnause = checkt(args.quantile)
    avguse = checkt(args.average)
    macuse = checkt(args.macs2)
    bwuse = checkt(args.bgToBw)

    ##Check_presence_of_'sge-.*\.o.*'-files
    ##---------------------------------------
    ##Note:Alternative_is_to_check_for_bedgraph-file_presence
    ##Note:Alternative_is_to_search_for_sge-file_by_jobIDs
    sys.stdout.write("\n>Checking presence of 'sge-.*\.o.*'-files\n")
    checkSl(args.repos, '^sge-.*\.o.*')

    ##Check_end_of_job_via_'sge-*'-files
    ##------------------------------------
    ##Note:alternative_is_to_check_for_last_coord_bedgraph-file/length
    sys.stdout.write("\n>Check complete 'sge-.*\.o.*'-files\n")
    nov = dict()
    for fIN in args.repos:
        sl = [f for f in os.listdir(fIN) if re.compile('^sge-.*\.o.*').search(f)][0]
        nov[os.path.join(fIN,sl)] = False

    sys.stdout.write('\tWaiting for jobs finishing.\n')
    rast = True
    while rast == True:
        for el in [k for k,v in nov.items() if v == False]:
            nov[el] = screener(el)
            if not all(value == True for value in nov.values()):
                time.sleep(5)
                continue
            else:
                rast = False
                sys.stdout.write('\tAll jobs finished.\n')
                break

    ##Rename_files_in_individual_dirs_&_initiate_peakcalling
    ##------------------------------------------------------
    sys.stdout.write('\n>Rename files & initiate peak calling\n')
    bGFs = list()
    damONs = list()
    jobIDs = list()
    for el in args.repos:

        absDIR = os.path.abspath(el)
        os.chdir(absDIR)
        sys.stdout.write('\t' + absDIR + '\n')

        nbGF, damOnlyNew, damNew, expNew = renamer(absDIR, args.ctrlpre, args.exppre)
        bGFs.append(nbGF)
        damONs.append(damOnlyNew)

        ##Create_peak_calling_commands_&_scripts
        ##--------------------------------------
        pcc = peakCalling(macuse,genSize,damNew,expNew)
        pccSH = create_sh(pcc,args.feedback)
        pccDO = peakCalling(macuse,genSize,damNew)
        pccDOSH = create_sh(pccDO,args.feedback)

        ##Submit_peak_calling_scripts
        ##---------------------------
        #sys.stdout.write("\tPC_trt-vs-ctrl:\t" + pccSH + "\n")
        jobID = submit(pccSH)
        jobIDs.append(jobID)
        #sys.stdout.write("\tPC_ctrl-alone:\t" + pccDOSH + "\n")
        jobID = submit(pccDOSH)
        jobIDs.append(jobID)

        os.chdir(oriDIR)

    #Check_all_peak_calling_jobs_are_queued
    ##-------------------------------------
    sys.stdout.write("\n>Check peak calling jobs\n")
    checkQue(jobIDs)

    ##Deduplicate_damONs
    ##------------------
    subDic = dict()
    subSet = set()
    for el in damONs:
        subDic[el] = os.path.basename(el)
        subSet.add(os.path.basename(el))

    subDONs = list()
    for setEL in subSet:
        subDONs.append([k for k,v in subDic.items() if v == setEL][0])

    ##Create_dirs_&_copy_bedgraph-files
    ##---------------------------------
    damONDIR = createDir(oriDIR,args.out,"_DamOnly_tracks",subDONs)
    bGFDIR = createDir(oriDIR,args.out,"_tracks", bGFs)

    ##Process_'*.bedgraph'_files
    ##--------------------------
    ##Quantile_normalize_all_bGFs
    jobID = quantNorm(oriDIR,bGFDIR,qnause,args.feedback)
    ##Check_normalization_job_finished
    checkFin([jobID])
    ##Average_all_normalized_bGFs
    jobID = average(oriDIR,bGFDIR,avguse,args.feedback)
    ##Ensure_all_jobs_are_finished
    checkFin([jobID])
    ##Convert_*.bedgraph_files_into_*.bw
    jobIDs = bwer(oriDIR,bGFDIR,args.chrSize,bwuse,args.feedback)
    ##Ensure_all_jobs_are_running
    checkQue(jobIDs)

    ##Process_DamOnly_'*.bedgraph'in_files
    ##------------------------------------
    ##Quantile_normalize_all_DamOnlyNew_files
    jobID = quantNorm(oriDIR,damONDIR,qnause,args.feedback)
    ##Check_normalization_job_finished
    checkFin([jobID])
    ##Average_all_normalized_damONs
    jobID = average(oriDIR,damONDIR,avguse,args.feedback)
    ##Ensure_all_jobs_are_running
    checkFin([jobID])
    ##Convert_*.bedgraph_files_into_*.bw
    jobIDs = bwer(oriDIR,damONDIR,args.chrSize,bwuse,args.feedback)
    ##Ensure_all_jobs_are_running
    checkQue(jobIDs)

    sys.stdout.write('\nAll done.\n')

if __name__ == '__main__':
    main()
