#include "MuFilterHit.h"
#include "MuFilter.h"
#include "ShipUnit.h"
#include "TROOT.h"
#include "FairRunSim.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include <TRandom.h>
#include <iomanip> 

using std::cout;
using std::endl;
using std::to_string;
using std::string;

// -----   Default constructor   -------------------------------------------
MuFilterHit::MuFilterHit()
  : SndlhcHit()
{
 flag = true;
 for (Int_t i=0;i<16;i++){fMasked[i]=kFALSE;}
}
// -----   Standard constructor   ------------------------------------------
MuFilterHit::MuFilterHit(Int_t detID)
  : SndlhcHit(detID)
{
 flag = true;
 for (Int_t i=0;i<16;i++){fMasked[i]=kFALSE;}
}
MuFilterHit::MuFilterHit(Int_t detID,Int_t nP,Int_t nS)
  : SndlhcHit(detID,nP,nS)
{
 flag = true;
 for (Int_t i=0;i<16;i++){
      fMasked[i]=kFALSE;
      signals[i]  = -999;
      times[i]    = -999;
      fDaqID[i]  = -1;
   }
}

// -----   constructor from MuFilterPoint   ------------------------------------------
MuFilterHit::MuFilterHit(Int_t detID, std::vector<MuFilterPoint*> V)
  : SndlhcHit()
{
     Int_t subsystem = floor(detID/10000);
     MuFilter* MuFilterDet = dynamic_cast<MuFilter*> (gROOT->GetListOfGlobals()->FindObject("MuFilter"));
     // get parameters from the MuFilter detector for simulating the digitized information
     nSiPMs = MuFilterDet->GetnSiPMs(detID);
     nSides = MuFilterDet->GetnSides(detID);
     
     Float_t attLength=0;
     Float_t SiPMcalibration=0;
     Float_t SiPMcalibrationS=0;
     Float_t timeresol_left, timeresol_right;
     Float_t signalspeed_left, signalspeed_right;

    //  Float_t DSpropspeed =0;
     if (subsystem==3) { 
              if (nSides==2){ attLength = MuFilterDet->GetConfParF("MuFilter/DsAttenuationLength"); }
              else { attLength = MuFilterDet->GetConfParF("MuFilter/DsTAttenuationLength"); }
              SiPMcalibration = MuFilterDet->GetConfParF("MuFilter/DsSiPMcalibration");
              
              signalspeed_left = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
              signalspeed_right = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");

              timeresol_left = MuFilterDet->GetConfParF("MuFilter/timeResol");
              timeresol_right = MuFilterDet->GetConfParF("MuFilter/timeResol");
     }
     else if (subsystem==2) {
              attLength = MuFilterDet->GetConfParF("MuFilter/VandUpAttenuationLength");
              SiPMcalibration = MuFilterDet->GetConfParF("MuFilter/VandUpSiPMcalibration");
              SiPMcalibrationS = MuFilterDet->GetConfParF("MuFilter/VandUpSiPMcalibrationS");

              timeresol_left = MuFilterDet->GetBarSideTimeResolution(detID, "left");
              timeresol_right = MuFilterDet->GetBarSideTimeResolution(detID, "right");

              signalspeed_left = MuFilterDet->GetBarSideSignalSpeed(detID, "left");
              signalspeed_right = MuFilterDet->GetBarSideSignalSpeed(detID, "right");

     }
     else { // if detID in veto
              if (subsystem==1 && nSides==1){ attLength = 2*MuFilterDet->GetConfParF("MuFilter/VandUpAttenuationLength"); }
              else { attLength = 2*MuFilterDet->GetConfParF("MuFilter/VandUpAttenuationLength"); }
              SiPMcalibration = MuFilterDet->GetConfParF("MuFilter/VandUpSiPMcalibration");

              signalspeed_left = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
              signalspeed_right = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
              
              timeresol_left = MuFilterDet->GetConfParF("MuFilter/timeResol");
              timeresol_right = MuFilterDet->GetConfParF("MuFilter/timeResol");
     }

     for (unsigned int j=0; j<16; ++j){
        signals[j] = -1;
        times[j] =-1;
     }
     LOG(DEBUG) << "detid "<<detID<< " size "<<nSiPMs<< "  side "<<nSides;

     fDetectorID  = detID;

     // Instance a load of floats used in next block
     Float_t signalLeft = 0, signalRight = 0, signal = 0;
     Float_t earliestToAL = 1E20, earliestToAR = 1E20;
     Float_t dxLphys=0, dxRphys=0, dxL=0, dxR=0;
     TVector3 vLeft,vRight;
     Float_t distance_Left=0, distance_Right=0;

     // for the timing, find earliest particle and smear with time resolution
     Float_t t_Left=0, t_Right=0;
     Float_t time;

     // Load bar positions
     MuFilterDet->GetPosition(fDetectorID,vLeft, vRight);
     auto x_ref = 0.5*( vLeft[0] + vRight[0] );

     for (auto p = std::begin(V); p!= std::end(V); ++p) {

        signal = (*p)->GetEnergyLoss();
        time = (*p)->GetTime();

        // Find distances from MCPoint centre to ends of bar 
        TVector3 impact((*p)->GetX(),(*p)->GetY() ,(*p)->GetZ() );
        distance_Left = (vLeft-impact).Mag();
        distance_Right = (vRight-impact).Mag();

        signalLeft+=signal/nSides*TMath::Exp(-distance_Left/attLength);
        signalRight+=signal/nSides*TMath::Exp(-distance_Right/attLength);

        if (subsystem==3) {
          signalspeed_left = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
          signalspeed_right = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
        }
        else if (subsystem==1) {
          signalspeed_left = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
          signalspeed_right = MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");
        }        
        else { // If detID in upstream
          signalspeed_left = MuFilterDet->GetBarSideSignalSpeed(detID, "left");
          signalspeed_right = MuFilterDet->GetBarSideSignalSpeed(detID, "right");
        }
        
        // Assume earliest arriving photons set the time on each side
        t_Left = time + distance_Left/signalspeed_left;
        t_Right = time + distance_Right/signalspeed_right;

        if ( t_Left < earliestToAL){
          earliestToAL = t_Left;
          dxLphys = impact[0]; // For aligned times in data, they correspond to distances from x=0 in the physics FoR
          dxL = -1*x_ref + dxLphys; // At the moment, the time alignment in data is done with respect to the bar centre. This quantity is wrt the bar centre 
          }
        
        if ( t_Right < earliestToAR){
          earliestToAR = t_Right;
          dxRphys = impact[0];
          dxR = x_ref - dxRphys; // At the moment, the time alignment in data is done with respect to the bar centre. This quantity is wrt the bar centre
          }
     }

     // In the SndlhcHit class the 'signals' array starts from 0.
     Float_t timeResol=0, aligned_time=0, SiPMcalibrationConstant=0; 
     TString side;
     for (unsigned int j=0; j<nSiPMs*nSides; ++j){
        
        // If small SiPM (j==2, 5, 10, 13) in HCAL
        if ( (subsystem==2) and (j%8==2 or j%8==5) ) { SiPMcalibrationConstant = SiPMcalibrationS;}
        else { SiPMcalibrationConstant = SiPMcalibration; }

        if ( (subsystem!=3 and j<8) || (subsystem==3 and j==0) ) { // If left-side or top channel
          timeResol=timeresol_left;
          aligned_time = dxL/signalspeed_left;
          signal=signalLeft;
        }
        else { // If right-side channel {
          timeResol=timeresol_right;
          aligned_time = dxR/signalspeed_right;          
          signal=signalRight;
        }    

        signals[j] = signal/float(nSiPMs) * SiPMcalibrationConstant;  // most simplest model, divide signal individually. Small SiPMs special
        times[j] = gRandom->Gaus(aligned_time, timeResol);
     }

     // Hard coding 0.720 MeV energy cut :sweat_smiling:
     // Assume energy deposition MPV of 1.8 MeV is 900 keV per side, assume 80% efficiency => 720 keV
     if (signalLeft < 0.4E-3 or signalRight < 0.4E-3) {flag=false;}
     else {flag = true;}
     
     for (Int_t i=0;i<16;i++){fMasked[i]=kFALSE;}
     LOG(DEBUG) << "signal created";
}

// -----   Destructor   ----------------------------------------------------
MuFilterHit::~MuFilterHit() { }
// -------------------------------------------------------------------------

// -----   Public method GetEnergy   -------------------------------------------
Float_t MuFilterHit::GetEnergy()
{
  // to be calculated from digis and calibration constants, missing!
  Float_t E = 0;
  for (unsigned int j=0; j<nSiPMs; ++j){
        E+=signals[j];
        if (nSides>1){ E+=signals[j+nSiPMs];}
  }
  return E;
}

bool MuFilterHit::isVertical(){
  if  ( (floor(fDetectorID/10000)==3&&fDetectorID%1000>59) ||
         (floor(fDetectorID/10000)==1&&int(fDetectorID/1000)%10==2) ) {  
      return kTRUE;
  }
  if  ( (floor(fDetectorID/10000)==3&&fDetectorID%1000>59) ||
         (floor(fDetectorID/10000)==1&&int(fDetectorID/1000)%10==2) ) {  
      return kTRUE;
  }
  else{return kFALSE;}
}

bool MuFilterHit::isShort(Int_t i){
  if (i%8==2 || i%8==5) {return kTRUE;}
  else{return kFALSE;}
}

// -----   Public method Get List of signals   -------------------------------------------
std::map<Int_t,Float_t> MuFilterHit::GetAllSignals(Bool_t mask,Bool_t positive)
{
          std::map<Int_t,Float_t> allSignals;
          for (unsigned int s=0; s<nSides; ++s){
              for (unsigned int j=0; j<nSiPMs; ++j){
               unsigned int channel = j+s*nSiPMs;
               if (signals[channel]<-900){continue;}
               if (signals[channel]> 0 || !positive){
                 if (!fMasked[channel] || !mask){
                    allSignals[channel] = signals[channel];
                    }
                }
              }
          }
          return allSignals;
}

// -----   Public method Get List of time measurements   -------------------------------------------
std::map<Int_t,Float_t> MuFilterHit::GetAllTimes(Bool_t mask, Bool_t apply_t_corr, Double_t SiPMDistance)
{
          OptForTimeCorrections(mask, apply_t_corr, SiPMDistance);
          std::map<Int_t,Float_t> allTimes;
          for (unsigned int s=0; s<nSides; ++s){
              for (unsigned int j=0; j<nSiPMs; ++j){
               unsigned int channel = j+s*nSiPMs;
               if (signals[channel]> 0){
                 if (!fMasked[channel] || !mask){
                    allTimes[channel] = fTimesHelper[channel];
                    }
                }
              }
          }
          return allTimes;
}

// -----   Public method Get time difference mean Left - mean Right   -----------------
Float_t MuFilterHit::GetDeltaT(Bool_t mask, Bool_t apply_t_corr, Double_t SiPMDistance)
// based on mean TDC measured on Left and Right
{
          OptForTimeCorrections(mask, apply_t_corr, SiPMDistance);
          Float_t mean[] = {0,0}; 
          Int_t count[] = {0,0}; 
          Float_t dT = -999.;
          for (unsigned int s=0; s<nSides; ++s){
              for (unsigned int j=0; j<nSiPMs; ++j){
               unsigned int channel = j+s*nSiPMs;
               if (signals[channel]> 0){
                 if (!fMasked[channel] || !mask){
                    mean[s] += fTimesHelper[channel];
                    count[s] += 1;
                    }
                }
              }
          }
          if (count[0]>0 && count[1]>0) {
                dT = mean[0]/count[0] - mean[1]/count[1];
          }
          return dT;
}
Float_t MuFilterHit::GetFastDeltaT(Bool_t mask, Bool_t apply_t_corr, Double_t SiPMDistance)
// based on fastest (earliest) TDC measured on Left and Right
{
          OptForTimeCorrections(mask, apply_t_corr, SiPMDistance);
          Float_t first[] = {1E20,1E20}; 
          Float_t dT = -999.;
          for (unsigned int s=0; s<nSides; ++s){
              for (unsigned int j=0; j<nSiPMs; ++j){
               unsigned int channel = j+s*nSiPMs;
               if (signals[channel]> 0){
                 if (!fMasked[channel] || !mask){
                    if  (times[channel]<first[s]) {first[s] = fTimesHelper[channel];}
                    }
                }
              }
          }
          if (first[0]<1E10 && first[1]<1E10) {
                dT = first[0] - first[1];
          }
          return dT;
}


// -----   Public method Get mean time  -----------------
Float_t MuFilterHit::GetImpactT(Bool_t mask, Bool_t apply_t_corr, Double_t SiPMDistance)
{
          OptForTimeCorrections(mask, apply_t_corr, SiPMDistance);
          Float_t mean[] = {0,0}; 
          Int_t count[] = {0,0}; 
          Float_t dT = -999.;
          Float_t dL;
          MuFilter* MuFilterDet = dynamic_cast<MuFilter*> (gROOT->GetListOfGlobals()->FindObject("MuFilter"));
          if (floor(fDetectorID/10000==3)) { 
             dL = MuFilterDet->GetConfParF("MuFilterDet/DownstreamBarX") / MuFilterDet->GetConfParF("MuFilter/DsPropSpeed");}
          else if (floor(fDetectorID/10000==2)) { 
             dL = MuFilterDet->GetConfParF("MuFilterDet/UpstreamBarX") / MuFilterDet->GetConfParF("MuFilter/VandUpPropSpeed");}
          else { 
             dL = MuFilterDet->GetConfParF("MuFilterDet/VetoBarX") / MuFilterDet->GetConfParF("MuFilter/VandUpPropSpeed");}

          for (unsigned int s=0; s<nSides; ++s){
              for (unsigned int j=0; j<nSiPMs; ++j){
               unsigned int channel = j+s*nSiPMs;
               if (signals[channel]> 0){
                 if (!fMasked[channel] || !mask){
                    mean[s] += fTimesHelper[channel];
                    count[s] += 1;
                    }
                }
              }
          }
          if (count[0]>0 && count[1]>0) {
                dT = (mean[0]/count[0] + mean[1]/count[1])/2.*6.25 -  dL/2.; // TDC to ns = 6.25
          }
          return dT;
}

std::map<TString,Float_t> MuFilterHit::SumOfSignals(Bool_t mask)
{   
/*    use cases, for Veto and DS small/large ignored
        sum of signals left large SiPM:    LL
        sum of signals right large SiPM: RL
        sum of signals left small SiPM:    LS
        sum of signals right small SiPM: RS
        sum of signals left and right:  
*/
          Float_t theSumL     = 0;
          Float_t theSumR    = 0;
          Float_t theSumLS   = 0;
          Float_t theSumRS  = 0;
          for (unsigned int s=0; s<nSides; ++s){
              for (unsigned int j=0; j<nSiPMs; ++j){
               unsigned int channel = j+s*nSiPMs;
               if (signals[channel]> 0){
                 if (!fMasked[channel] || !mask){
                    if (s==0 and !isShort(j)){theSumL+= signals[channel];}
                    if (s==0 and isShort(j)){theSumLS+= signals[channel];}
                    if (s==1 and !isShort(j)){theSumR+= signals[channel];}
                    if (s==1 and isShort(j)){theSumRS+= signals[channel];}
                    }
                }
              }
          }
         std::map<TString,Float_t> sumSignals;
         sumSignals["SumL"]=theSumL;
         sumSignals["SumR"]=theSumR;
         sumSignals["SumLS"]=theSumLS;
         sumSignals["SumRS"]=theSumRS;
         sumSignals["Sum"]=theSumL+theSumR;
         sumSignals["SumS"]=theSumLS+theSumRS;
         return sumSignals;
}

// -----   Public method Print   -------------------------------------------
void MuFilterHit::Print() const
{
  std::cout << "-I- MuFilterHit: MuFilter hit " << " in detector " << fDetectorID;

  if ( floor(fDetectorID/10000)==3&&fDetectorID%1000>59) {
     std::cout << " with vertical bars"<<std::endl;
     std::cout << "top digis:";
     for (unsigned int j=0; j<nSiPMs; ++j){
         std::cout << signals[j] <<" ";
     }
  }else{
     std::cout << " with horizontal bars"<<std::endl;
     for (unsigned int s=0; s<nSides; ++s){
       if (s==0) {std::cout << "left digis:";}
       else {std::cout << "right digis:";}
       for (unsigned int j=0; j<nSiPMs; ++j){
         std::cout << signals[j] <<" ";
      }
     }
 }
std::cout << std::endl;
}

// -----   Private method to opt for time corrections   -------------------------------------------
void MuFilterHit::OptForTimeCorrections(Bool_t mask, Bool_t apply, Double_t SiPMDistance)
{
          // simply copy the times if no time corrections is required
          if (apply == kFALSE) memcpy(fTimesHelper, times, sizeof(fTimesHelper));
          else { // apply the time corrections
            MuFilter* MuFilterDet = dynamic_cast<MuFilter*> (gROOT->GetListOfGlobals()->FindObject("MuFilter"));
            for (unsigned int s=0; s<nSides; ++s){
                for (unsigned int j=0; j<nSiPMs; ++j){
                 unsigned int channel = j+s*nSiPMs;
                 if (signals[channel]> 0){
                   if (!fMasked[channel] || !mask){
                     fTimesHelper[channel] = (MuFilterDet->GetCorrectedTime(fDetectorID,channel,
                                                                    times[channel]*6.25,
                                                                    SiPMDistance,
                                                                    signals[channel])*6.25);
                      }
                 }
                 else fTimesHelper[channel] = times[channel];
                }
            }
          }
}
// -------------------------------------------------------------------------

ClassImp(MuFilterHit)

