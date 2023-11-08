#!/usr/local/bin/python3
'''
#Parse_files_into_'damMer.py'-script:
dam=($(find . -type f -iname "*.fastq.gz" -and -iname "dam_*"))
exp=($(find . -type f -iname "*.fastq.gz" -and -iname "experiment_*"))
python3 ~/Desktop/damMer.py -e "${exp[@]}" -c "${dam[@]}"
'''

import argparse
import os
import sys
import logging
import shutil
import gzip
import re
import errno
import subprocess
import shlex
import time
from difflib import SequenceMatcher

shItr = 1
tmpl = """\
#!/bin/bash
#!
#! Name of the job:q
#$ -N {name}
#$ -V
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
    parser = argparse.ArgumentParser(description="Wrapper script for damidseq_pipeline.")

    parser.add_argument(
        "-e", "--experiment",
        nargs = '*',
        type = str,
        required = True,
        help = "List of experimental '*.fastq.gz'-files for the Dam-fusion samples."
        )
    parser.add_argument(
        "-c", "--control",
        nargs = '*',
        type = str,
        required = True,
        help = "List of control '*.fastq.gz'-files for the Dam-only samples."
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
        "-i", "--index",
        type = str,
        default = None,
        help = "'bowtie2_build'-derived genome index."
        )
    parser.add_argument(
        "-g", "--gatcfrag",
        type = str,
        default = None,
        help = "'*.GATC.gff'-file listing coordinates of GATC-fragments."
        )
    parser.add_argument(
        "-b", "--bow2dir",
        type = str,
        default = "/usr/bin/bowtie2",
        help = "Path to bowtie2 executables."
        )
    parser.add_argument(
        "-s", "--samdir",
        type = str,
        default = "usr/bin/samtools",
        help = "Path to samtools_mt executables."
        )
    parser.add_argument(
        "-q", "--damidseq",
        type = str,
        default = "./damidseq_pipeline_vR.1.pl",
        help = "Path to damidseq_pipeline executable."
        )

    arguments = parser.parse_args()
    return arguments

##-----------------##
##----Functions----##
##-----------------##

def filing(allargs):
    '''Generate ''*.log'-file & save arguments'''

    global prefixq
    prefix = os.path.basename(dir)
    global logname
    logname = dir + '/' + prefix + ".log"
    sys.stdout.write("\t" + logname + "\n")
    logging.basicConfig(
        filename = logname,
        level = logging.DEBUG,
        format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )

    for arg in vars(allargs):
        logging.info('%s: %s', arg, getattr(allargs,arg))
    return(prefix)

def checkt(toolPath):
    '''Checking paths of used tools.'''

    tool = os.path.basename(toolPath)
    logging.info('Checking: ' + tool)
    sys.stdout.write('\t'+ tool + '\n')

    if os.path.exists(toolPath):
        tcho = toolPath
    else:
        tcho = ""

    if shutil.which(tool):
        tdet = shutil.which(tool)
    else:
        tdet = ""

    #Tester
    #sys.stdout.write('\nChosen dir: ' + tcho + '\n')
    #sys.stdout.write('\nDetected dir: ' + tdet + '\n')

    if tcho == tdet and tdet != "" and tcho != "":
        ##Chosen_&_detected_dirs_exist_and_are_similar

        logging.info(tool + ' checked: ' + tdet)
        tuse = tdet
        return tuse

    elif tcho != tdet and tdet != "" and tcho != "":
        ##Chosen_&_detected_dirs_exist_but_are_not_similar

        logging.warning('"' + tcho + '" not default installation: "' + tdet + '"')
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

        logging.info(tool + ' used: ' + tuse)
        sys.stdout.write(tool + ' used: ' + tuse + '\n')
        return tuse

    elif tcho != tdet and tdet != "" and tcho == "":
        ##Detected_but_not_chosen_dir_exist
        logging.warning('"' + toolPath + \
            '" not default installation: "' + tdet + '"')
        sys.stdout.write(toolPath + \
            '" not default installation: "' + tdet + '"\n')

        logging.warning('Not found: ' + toolPath)
        tuse = tdet
        logging.info(tool + ' used: ' + tuse)
        sys.stdout.write(tool + ' used: ' + tuse + '\n')
        return tuse

    elif tcho != tdet and tdet == "" and tcho != "":
        ##Chosen_but_not_detected_dir_exist
        tuse = tcho
        logging.info(tool + ' used: ' + tuse)
        sys.stdout.write(tool + ' used: ' + tuse + '\n')
        return tuse

    elif tcho =="" and tdet == "":
        ##Neither_chosen_nor_detected_dir_exists
        logging.error('Not found: ' + toolPath)
        logging.error('Not found: ' + tool)
        sys.exit('Error: Not found: ' + tool + \
            '\nUse python3 bow2map.py --help.\n')

def checki(indices):
    '''Check path and validity of bowtie2-indices.'''

    logging.info('Checking: ' + indices)
    sys.stdout.write('\t' + indices + '\n')

    #default = "/mnt/home1/brand/rk565/resources/bowtie2_BDGP6.ensembl/Ensembl_BDGP6_genome"
    suffices=['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    for s in suffices:
        index = indices + s
        if not os.path.isfile(index):
            logging.error('Index missing: ' + index)
            sys.exit('Error: Index missing: ' + index + '\n')

    logging.info(os.path.basename(indices) + '-indices validated.')
    #sys.stdout.write(os.path.basename(indices) + '-indices validated.\n')

def binary_tester(fqFile):
    '''Test binary status of fastq-file.'''

    fqF = open(fqFile, 'rb')
    try:
        while True:
            block = fqF.read(1024)
            if b'\0' in block:
                return True
            if len(block) < 1024:
                break
    finally:
        fqF.close()

    return False

def checkf(fastq):
    '''Check path and validity of fastq-file.'''

    logging.info('Checking: ' + fastq)
    sys.stdout.write('\t'+ fastq + '\n')

    ##Check_path_of_fastq-file
    if os.path.isfile(fastq):
        logging.info('File exists: ' + fastq)
    else:
        logging.warning('File not found: ' + fastq)
        sys.stderr.write('WARNING: File not found: ' + fastq + \
            '\nUse python3 bow2map.py --help.\n')
        return

    ##Check_whether_fastq_is_compressed
    ##---------------------------------
    binary = binary_tester(fastq)
    logging.info(fastq + ' binary: ' + str(binary))

    ##Read_in_sample_of_first_200_lines
    ##---------------------------------
    if binary:
        try:
            with gzip.open(fastq, "r") as fq:
                header = [next(fq) for i in range(200)]
        except StopIteration:
            logging.warning('File has less than 200 lines.')
            sys.stderr.write('WARNING: File has less than 200 lines.\n')
            return

        header = [header[i].decode("utf-8") for i in range(200)]
    else:
        try:
            with open(fastq, "r") as fq:
                header = [next(fq) for i in range(200)]
        except StopIteration:
            logging.warning('File has less than 200 lines.')
            sys.stderr.write('WARNING: File has less than 200 lines.\n')
            return

    ##Check_beginning_of_every_fastq_entry
    ##------------------------------------
    j = 0
    while j < 200:
        if bool(re.compile('^@').search(header[j])):
            j += 4
            continue
        else:
            logging.warning('Read header not starting with "@". line: ' \
                + str(j) + ' in ' + fastq)
            sys.stderr.write('WARNING: Read header not starting with "@".\n \
                See: ' + logname + '\n')
            return

    logging.info("All lines start with '@'.")

    ##Check_sequence_for_DNA-bases_&_quality-scores
    ##---------------------------------------------
    k = 1
    while k < 200:
        if not re.compile('[^AGCTN]').search(header[k].replace('\n', '')):
            if len(header[k].replace('\n', '')) == len(header[k+2].replace('\n', '')):
                k += 4
                continue
            else:
                logging.warning('Length of sequence doesn\'t match quality \
                    value length. line: ' + str(k) + ' in ' + fastq)
                sys.stderr.write('WARNING: Length of sequence doesn\'t match quality \
                    value length.\nSee: ' + logname + '\n')
                k += 4
                continue
        else:
            logging.warning('Sequence contains non-\'GACT\'-letters. \
                line: ' + str(k) + ' in ' + fastq)
            sys.stderr.write('WARNING: Sequence contains non-\'GACT\'-letters.\n \
                See: ' + logname + '\n')
            k += 4

    logging.info(fastq + ' validated.')
    #sys.stdout.write(fastq + ' validated.\n')

    return(fastq)

def evalDir(path):
    '''Create actual path for dir.'''

    try:
        os.makedirs(path)
        return(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def matcher(strings):
    '''Identify common denominator string amongst filenames.'''

    strings = [
        re.compile('(.*/.*?)(\..*)+$').search(stri).group(1) \
        for stri in strings\
        ]

    match = os.path.basename(strings[0])
    for stri in strings[1:]:
        stri = os.path.basename(stri)
        res = SequenceMatcher(None, match, stri).find_longest_match(0, len(match), 0, len(stri))
        match = match[res.a: res.a + res.size]


    sys.stdout.write("\tPrefix: " + match + "\n")
    return(match)

def create_sh(cmd,mailAc,cores=2):
    '''Create submission script for current command.'''

    global shItr
    cmdName = re.compile('\..*').sub('', os.path.basename(cmd.split(" ")[0]))
    fileName = dir + "/" + str(shItr) + "_" + cmdName + ".sh"
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
    '''
    Submit the script for the current command.
    dpdIDs should have format: '<jobID1>,<jobID2>,etc.'
    '''

    Sub = "qsub -terse"
    if dpdIDs:
        Sub += " -hold_jid=" + dpdIDs
    Sub += " " + cmdSH
    print("submit: created command '%s'" % cmdSH)
    #sys.stdout.write("\nSub command: " + str(shlex.split(Sub)) + "\n")
    try:
        prc = subprocess.Popen(
            shlex.split(Sub),
            shell = False,
            stderr = subprocess.PIPE,
            stdout = subprocess.PIPE
            )
    except Exception as e:
        sys.exit("devscript_damMer.py aborted due to:\t" + type(e).__name__ + "\n")

    jobID = str(prc.communicate()[0].decode("utf-8")).rstrip()
    #sys.stdout.write("ID of current job: " + jobID + "\n")

    print("submit: returned job ID '%s'" % jobID)
    return(jobID)

def checkFin(jobIDs):
    '''Check all provided jobIDs are no longer registered by sge.'''

    sys.stdout.write('\nWaiting for cluster.\n')
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
                    print("...... extracted job ID '%s'" % line.split[0])
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
    sys.stdout.write('Job(s) finished.\n')

def checkQue(jobIDs):
    '''Check if jobs are registered by sge.'''

    sys.stdout.write('\nWaiting for cluster.\n')
    print("checkQue: getting output from 'qstat'")
    chk = subprocess.check_output(['qstat'])
    qstat_list = []
    for line in chk.decode("utf-8").split('\n'):
        print("QSTAT> %s" % line.rstrip())
        try:
            if line.split()[0].isdigit():
                qstat_list.append(line.split()[0])
                print("...... extracted job ID '%s'" % line.split[0])
        except IndexError:
            pass
    chk = qstat_list
    print("checkQue: job list = %s" % chk)
    if set(jobIDs).issubset(chk):
        sys.stdout.write('Job(s) running.\n')
    else:
        sys.exit("One or more job(s) not running.\n")

##---------------------##
##----Main_workflow----##
##---------------------##

def main():
    args = parse_args()

    ##Set_global_variable_'dir'
    ##-------------------------
    global dir
    dir = os.path.dirname(os.path.abspath(args.experiment[0]))

    ##Set_variable_for_index_directory
    ##--------------------------------
    '''In absence of specified 'defaults', 'index' needs to be provided.'''
    if args.index == None:
        inDir = "resources"
        if args.defaults == "dm6":
            args.index = '/'.join(
                [inDir,"bowtie2_BDGP6.ensembl/Ensembl_BDGP6_genome"]
            )
            args.gatcfrag = '/'.join(
                [inDir, "Ensembl_BDGP6.GATC.mod.gff"]
            )
        elif args.defaults == "mm10":
            args.index = '/'.join(
                [inDir, "bowtie2_GRCm38.ensembl/bowtie2_GRCm38.ensembl"]
                #[inDir, "bowtie2_GRCm38.exclContigs.masked/mm10.chrom.masked"]
            )
            args.gatcfrag = '/'.join(
                [inDir, "bowtie2_GRCm38.ensembl/Mus_musculus.GRCm38.dna.primary_assembly.GATC.gff"]
            )
        else:
            sys.exit('Unsupported species: --defaults=[dm6/mm10].\n')

    ##Generate_&_initiate_logfile
    ##---------------------------
    sys.stdout.write('\n>Logfile\n')
    preOut = filing(args)

    ##Checking_indices_&_executables
    ##------------------------------
    sys.stdout.write('\n>Checking executables\n')
    bowuse = checkt(args.bow2dir)
    samuse = checkt(args.samdir)
    damuse = checkt(args.damidseq)
    sys.stdout.write('\n>Checking indices\n')
    checki(args.index)

    ##Checking_all_fastq-files
    ##------------------------
    sys.stdout.write("\n>Checking '*.fastq.gz'-files\n")
    exps = list()
    for f in args.experiment:
        fckd = checkf(f)
        if fckd is not None:
            exps.append(fckd)

    expsPre = matcher(exps)
    #expsPre = re.compile('_|\.').sub('', expsPre)

    ctrls = list()
    for c in args.control:
        cckd = checkf(c)
        if cckd is not None:
            ctrls.append(cckd)

    ctrlsPre = matcher(ctrls)
    #ctrlsPre = re.compile('_|\.').sub('', ctrlsPre)

    ##Create_WDs_&_copy_*.fastq.gz-files
    ##----------------------------------
    sys.stdout.write('\n>Create directories & copy files\n')
    dirs = list()
    jobIDs = list()
    for e in exps:
        eb = os.path.basename(e)
        eb = re.compile('\..*\..*|\..*').sub('', eb)

        for d in ctrls:
            db = os.path.basename(d)
            db = re.compile('\..*\..*|\..*').sub('', db)

            dirName = dir + "/" + eb + "-vs-" + db + "/"
            sys.stdout.write('\t' + eb + '-vs-' + db + '/\n')
            evalDir(dirName)
            dirs.append(dirName)

            ##Create_copy_script
            cpy = "cp" + \
                " " + e + \
                " " + dirName + \
                "; cp " + d + \
                " " + dirName
            cpSH = create_sh(cpy,args.feedback)
            #sys.stdout.write("\nList script: " + cpSH + "\n")

            crnID = submit(cpSH)
            jobIDs.append(crnID)

    #jobIDs = [str(elem) for elem in jobIDs]
    cpJobs = (',').join(jobIDs)

    ##Tester-----------------------------------------------------
    #
    #sys.stdout.write('Number of dirs: ' + str(len(dirs)) + '\n')
    #sys.stdout.write('Pending jobs: ' + cpJobs + '\n')
    #
    ##-----------------------------------------------------------

    ##Wait_until_all_dirs_include_both_'*.fastq.gz'-files
    ##---------------------------------------------------
    sys.stdout.write("\n>Check for presence of all '*.fastq.gz'-files\n")
    rest = True
    interval = 10
    sys.stdout.write('\tWaiting for cluster.\n')
    while rest == True:
        cou = 0
        for dirCurr in dirs:
            if len(os.listdir(dirCurr))==2:
                cou += 1
                if not cou == len(dirs):
                    continue
                else:
                    rest = False
                    break
            else:
                time.sleep(interval)
                break
    sys.stdout.write('\tAll files copied.\n')

    ##Run_damid_for_all_combinations
    ##------------------------------
    sys.stdout.write("\t>Initialize damidseq_pipeline_vR.1 in all directories\n")
    jobIDs = list()
    for cwd in dirs:
        os.chdir(cwd)
        sys.stdout.write('\t' + cwd + '\n')
        fs = [f for f in os.listdir(cwd) if os.path.isfile(os.path.join(cwd, f))]
        dam = [f for f in fs if re.search(ctrlsPre, f)]
        #sys.stdout.write('dam: '+str(dam)+'\n')
        exp = [f for f in fs if re.search(expsPre, f)]
        #sys.stdout.write('exp: '+str(exp)+'\n')

        ##'damid'_command
        dsq = damuse + \
            " --bins=300" + \
            " --gatc_frag_file=" + args.gatcfrag + \
            " --bowtie2_genome_dir=" + args.index + \
            " --samtools_path=" + os.path.dirname(samuse) + "/" + \
            " --bowtie2_path=" + os.path.dirname(bowuse) + "/" + \
            " --dam=" + cwd + dam[0] + \
            " " + cwd + exp[0]
        dsqSH = create_sh(dsq,args.feedback,cores=8)
        #sys.stdout.write("\nList script:\t" + dsqSH + "\n")

        jobID = submit(dsqSH, dpdIDs=cpJobs)
        jobIDs.append(jobID)

    ##Ensure_all_jobs_are_running
    ##---------------------------
    time.sleep(10)
    sys.stdout.write('\n>Check all jobs are registered by sge\n')
    checkQue(jobIDs)

if __name__ == '__main__':
    main()
