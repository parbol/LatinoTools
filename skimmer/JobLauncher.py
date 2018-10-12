#//----------------------------------------------------------------------//
#// ___  ___                    _____           _                        //
#// |  \/  |                   /  ___|         | |                       //
#// | .  . |_   _  ___  _ __   \ `--. _   _ ___| |_ ___ _ __ ___  ___    //
#// | |\/| | | | |/ _ \| '_ \   `--. \ | | / __| __/ _ \ '_ ` _ \/ __|   //
#// | |  | | |_| | (_) | | | | /\__/ / |_| \__ \ ||  __/ | | | | \__ \   //
#// \_|  |_/\__,_|\___/|_| |_| \____/ \__, |___/\__\___|_| |_| |_|___/   //
#//                                    __/ |                             //
#//----------------------------------------------------------------------//
#// A project by: C. Diez, P. Gomez and P. Martinez                      //
#//----------------------------------------------------------------------//
#//----------------------------------------------------------------------//
#// muonProducer.py                                                      //
#//----------------------------------------------------------------------//
#// Running sequentially the production.                                 //
#//                                                                      //
#//----------------------------------------------------------------------//
#//----------------------------------------------------------------------//


#!/usr/bin/env python

import os, sys, subprocess
from optparse import OptionParser


#Coloured output
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


#Absolute paths
basepath = '/home/pablom/Documentos/DATA2/'
runpath = '/home/pablom/Documentos/work/MTP/'


#Functions
def showBanner():

    print ''
    print ''
    print ''
    
    print bcolors.HEADER + ' _____ _    _                       _                            _                    ' + bcolors.ENDC
    print bcolors.HEADER + '/  ___| |  (_)                     | |                          | |                   ' + bcolors.ENDC
    print bcolors.HEADER + '\ `--.| | ___ _ __ ___   ___ _ __  | |     __ _ _   _ _ __   ___| |__   ___ _ __      ' + bcolors.ENDC
    print bcolors.HEADER + ' `--. \ |/ / | \'_ ` _ \ / _ \ \'__| | |    / _` | | | | \'_ \ / __| \'_ \ / _ \ \'__|     ' + bcolors.ENDC
    print bcolors.HEADER + '/\__/ /   <| | | | | | |  __/ |    | |___| (_| | |_| | | | | (__| | | |  __/ |        ' + bcolors.ENDC 
    print bcolors.HEADER + '\____/|_|\_\_|_| |_| |_|\___|_|    \_____/\__,_|\__,_|_| |_|\___|_| |_|\___|_|        ' + bcolors.ENDC
    print ''
    print ''
    print ''
                                                                                                                      



if __name__ == "__main__":
    
    
    showBanner()	

    parser = OptionParser(usage="%prog --help")
    parser.add_option("-d", "--directory",    dest="directory",     type="string",   default='none',   help="Directory containing the files that should be skimmed.")
    parser.add_option("-o", "--output",       dest="output",        type="string",   default='none',   help="Directory where the skimmed files will be located.")
    parser.add_option("-t", "--tag",          dest="tag",           type="string",   default='all',    help="The name of the sample that will be processed.")
    parser.add_option("-n", "--name",         dest="skimName",      type="string",   default='none',   help="The name of the skim.")
    parser.add_option("-q", "--queue",        dest="queue",         type="string",   default='8nh',    help="Name of the queue.")
    parser.add_option("-c", "--cmssw",        dest="cmsswpath",     type="string",   default='none',   help="Name of the CMSSW path (src).")
    parser.add_option("-w", "--work",         dest="workdirectory", type="string",   default='./',     help="Name of the working directory.")
    (options, args) = parser.parse_args()


    ############################ Checking all the inputs are OK ##################################
    if not os.path.exists(options.directory):
        print bcolors.FAIL + 'The path ' + options.directory + ' does not exist. ' + bcolors.ENDC
        exit()
    if not os.path.exists(options.output):
        print bcolors.FAIL + 'The path ' + options.output + ' does not exist. ' + bcolors.ENDC
        exit()
    if not os.path.exists(options.cmsswpath):
        print bcolors.FAIL + 'The path ' + options.cmsswpath + ' does not exist. ' + bcolors.ENDC
        exit()
    if not os.path.exists(options.workdirectory):
        print bcolors.FAIL + 'The path ' + options.workdirectory + ' does not exist. ' + bcolors.ENDC
        exit()
    if not os.path.exists(options.workdirectory + '/runSkim'):
        print bcolors.FAIL + 'The runSkim executable cannot be found in' + options.workdirectory + bcolors.ENDC
        exit()


    ############################### Collecting the list of files ##################################
    filesToProcess = []
    for f in os.listdir(options.directory):
        if options.tag == 'all':
            filesToProcess.append(f)
        else:    
            if f.find("_" + options.tag + "__") != -1:
                filesToProcess.append(f)
    if len(filesToProcess) == 0:
        print bcolors.FAIL + 'The tag ' + options.tag + ' was not found.' + bcolors.ENDC
        exit()
   
    ############################ Give information about the setup ##################################
    print bcolors.OKBLUE + "Processing " + bcolors.OKGREEN + str(len(filesToProcess)) + bcolors.OKBLUE + " files in " + bcolors.OKGREEN + options.directory + bcolors.ENDC
    print ''
    outputFilesToProcess = []
    for f in filesToProcess:
        outputname = options.output + '/' + f[0:f.find('.root')] + '_' + options.skimName + '.root'
        outputFilesToProcess.append(outputname)
        print bcolors.OKGREEN + f + bcolors.WARNING + ' -> ' + bcolors.OKGREEN + outputname + bcolors.ENDC    	

    
    #################################### Preparing the launcher ##################################
    launcher = options.workdirectory + "/launcher.sh"
    print bcolors.OKBLUE + "Generating: " + bcolors.OKGREEN + launcher + bcolors.ENDC
    thelauncher = open(launcher, 'w')
    thelauncher.write('#/bin/bash\n')
    thelauncher.write('cd ' + options.cmsswpath + '\n')
    thelauncher.write('eval `scramv1 runtime -sh`\n')
    thelauncher.write('cd ' + options.workdirectory + '\n')
    thelauncher.write('./runSkim $1 $2\n')
    thelauncher.close()

    ############################### Preparing the run submitter #################################
    print ''
    nameOfSh = options.workdirectory + '/skimmer_' + options.tag + '_' + options.skimName + ".sh"
    print bcolors.OKBLUE + "Generating: " + bcolors.OKGREEN + nameOfSh + bcolors.ENDC
   
    thesh = open(nameOfSh, 'w')
    thesh.write('#/bin/bash\n')
    thesh.write('chmod +x launcher.sh\n')
    for forigin, foutput in zip(filesToProcess, outputFilesToProcess):
        origin = options.output + '/' + forigin
        destiny = foutput
        log = foutput + ".log"
        err = foutput + ".err"
        thesh.write('bsub -o ' + log + ' -e ' + err +  ' -q ' + options.queue + ' launcher.sh ' + options.directory + '/' + forigin + ' ' + destiny + '\n')
    thesh.close()
 
    print ''
    print bcolors.OKGREEN + "Program successful "  + bcolors.ENDC
    print ''

   

    





