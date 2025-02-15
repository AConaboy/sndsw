import ROOT
import os,sys
import re
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-f", "--inputFile", dest="inputFile", help="input file data and MC",default="",required=True)
parser.add_argument("--simMode", dest="simMode", required=True)
# parser.add_argument("--outFile", dest="outFile", help="outfile data and MC",default="./output.root",required=False)
options = parser.parse_args()
input_file = options.inputFile

# Make output file name by swapping in my eos path and dropping digi at end of filename
eospath=f'/eos/experiment/sndlhc/users/aconsnd/simulation/{options.simMode}/data/'
if not os.path.exists(eospath):
    print('No existing path:\n', eospath)
    os.exit()
outFile = f'{eospath}{options.inputFile.split("/")[-1].replace("-digCPP", "")}'

# print(outFile)
# sys.exit()

logger = ROOT.FairLogger.GetLogger()
logger.SetColoredLog(True)
logger.SetLogVerbosityLevel('low')
logger.SetLogScreenLevel('WARNING')
logger.SetLogToScreen(True)

run = ROOT.FairRunAna()
ioman = ROOT.FairRootManager.Instance()

f=ROOT.TFile.Open(input_file)
eventTree = f.cbmsim
runId = 'sim'

source = ROOT.FairFileSource(f)
run.SetSource(source)
sink = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)

#avoiding some error messages
xrdb = ROOT.FairRuntimeDb.instance()
xrdb.getContainer("FairBaseParSet").setStatic()
xrdb.getContainer("FairGeoParSet").setStatic()

run.Init()
OT = sink.GetOutTree()
eventTree = ioman.GetInTree()


# Write the TBranchList to the file
branches_to_copy = {"MCTrack":"ShipMCTrack","EmulsionDetPoint":"EmulsionDetPoints","ScifiPoint":"ScifiPoints","MuFilterPoint":"MuFilterPoints"}

#out_f = ROOT.TFile("SimonaWay.root","RECREATE")
#out_f.cd()
for b_name in branches_to_copy:
    b=ioman.GetObject(b_name)
    print(b_name,branches_to_copy[b_name])
    ioman.Register(b_name,branches_to_copy[b_name],b,ROOT.kTRUE)
#run.Run(0,100)
B = ROOT.TList()
B.SetName('BranchList')
for aBranch in branches_to_copy:
    B.Add(ROOT.TObjString(aBranch))
ioman.SetBranchNameList(B)
ioman.WriteFolder()
for n in range(ioman.GetInTree().GetEntries()):
    ioman.GetInTree().GetEvent(n)
#    ioman.GetOutTree().Fill()
    ioman.Fill()
#out_f.cd()
ioman.Write()    
#ioman.UpdateBranches()
#ioman.WriteFolder()
#out_f.cd()
#ioman.GetOutTree().Write()
#B.Add(ROOT.TObjString("ScifiPoints"))
#out_f.WriteObject(B,"BranchList",str(ROOT.TObject.kSingleKey))
#T = ROOT.TList()
#T.SetName('TimeBasedBranchList')
#out_f.WriteObject(T,'TimeBasedBranchList',str(ROOT.TObject.kSingleKey))

# Close the ROOT file
#out_f.Close()
