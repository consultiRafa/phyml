#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author: Imran Fanaswala (imranf@student.ethz.ch)

import subprocess, shlex, sys, os, time, argparse    
from os.path import join

_PHYML = None
_RESULTS_DIR = "./datasets/results"

_DATASETS_NT = [
  "./datasets/toy.txt",
  "./datasets/16.codon.paml",
  "./datasets/10.codon.paml",
  "./datasets/5.codon.paml",
  #"./datasets/4.codon.paml",
  #"./datasets/1.codon.paml",
  ]

_DATASETS_AA = [
  "./datasets/medium_protein/proteic_M1992_46x68_2004.phy",
  "./datasets/medium_protein/proteic_M1990_35x153_2004.phy",
  "./datasets/medium_protein/proteic_M2476_30x719_2005.phy",
  ]
  
  
#These commands will be run for EACH dataset. IOW, len(_DATASETS)*len(_CLI_PARAMS) commands
_CLI_PARAMS_NT = [
  r"{0} -i {1} -d nt -q -c 2 -v 0 -t 4         -m JC69  -f '0.25,0.25,0.25,0.25' -o n   -b -5 --r_seed 1999",
  r"{0} -i {1} -d nt -q -c 3 -v 0 -t 2         -m TN93  -f '0.01,0.04,0.05,0.90' -o lr  -b -5 --r_seed 1999",
  r"{0} -i {1} -d nt -q -c 4 -v 0 -t e -s NNI  -m F84   -f e                     -o r   -b -5 --r_seed 1999",
  r"{0} -i {1} -d nt -q -c 4 -v 0 -t e -s SPR  -m GTR   -f e                     -o tlr -b -5 --r_seed 1999",
  r"{0} -i {1} -d nt -q -c 7 -v 0 -t e -s BEST -m HKY85 -f e                     -o tlr -b -5 --r_seed 1999",
  ]
  
_CLI_PARAMS_AA = [
  r"{0} -i {1} -d aa -c 2 -v 0         -m WAG        -o n   -b -5 --r_seed 1999",
  r"{0} -i {1} -d aa -c 4 -v 0 -s BEST -m JTT   -f e -o r   -b -5 --r_seed 1999",
  r"{0} -i {1} -d aa -c 5 -v 0 -s NNI  -m MtArt -f e -o tlr -b -5 --r_seed 1999",
  r"{0} -i {1} -d aa -c 7 -v 0 -s NNI  -m HIVw  -f e -o tlr -b -5 --r_seed 1999",  
  ]

def generate_commands():
  commands = []
  for ds in _DATASETS_NT:
    for p in _CLI_PARAMS_NT:
      commands.append(p.format(_PHYML, ds))
  for ds in _DATASETS_AA:
    for p in _CLI_PARAMS_AA:
      commands.append(p.format(_PHYML, ds))      
  return commands
  
def main():
  #parse arguments
  parser = argparse.ArgumentParser()
  parser.add_argument("--beagle", help="Use PhyML-BEAGLE", action="store_true")
  args = parser.parse_args()
  global _PHYML 
  _PHYML = "./src/phyml-beagle" if args.beagle else "./src/phyml"
  
  commands = generate_commands()
  fname=join(_RESULTS_DIR,time.strftime("%d_%m_%Y_%H:%M:%S[%a_%d_%b_%Y]", time.localtime())+"_result.txt")
  
  with open(fname,"w") as logfile:
    #First show what will be executed
    logfile.write("Commands executed:\n")
    for cmd in commands:
      print cmd
      logfile.write(cmd+"\n")
      #Force flushing..
      logfile.flush()
      os.fsync(logfile.fileno())
    #And now execute...
    sys.stdout.write("\n===============================================================\nStarting batch job\nAll results (stdout AND stderr) will be dumped into {0}\nYou can inspect this file while the job is running\n===============================================================\n".format(fname))      
    for cmd in commands:
      cmd_list = shlex.split(cmd)
      sys.stdout.write("\nExecuting:"+cmd+"\n")
      subp = subprocess.Popen(cmd_list, shell=False, stdout=logfile, stderr=subprocess.STDOUT)
      subp.wait()
      #Force flushing..
      logfile.flush()
      os.fsync(logfile.fileno())
  
  #Write the timings back into the file
  with open(fname,"r+a") as logfile:
    #Extract the timings
    timings = filter(None,[(line.partition("Time used")[2]).rstrip("\n") for line in logfile])
    logfile.seek(0,os.SEEK_END)
    logfile.write("\n\nTimings Summary:\n")
    for i,cmd in enumerate(commands):
      logfile.write(timings[i]+"  "+cmd+"\n")
    if len(commands)!=len(timings):
      logfile.write("\nWarning: The number of commands executed does not match the number of timings I retreived\n")
   
      
    

if __name__ == "__main__":
  main()
