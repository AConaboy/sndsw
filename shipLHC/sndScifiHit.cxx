#include "sndScifiHit.h"
#include "Scifi.h"
#include "TROOT.h"
#include "FairRunSim.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include <TRandom.h>
#include <iomanip> 

namespace {
     //parameters for simulating the digitized information
     const Float_t ly_loss_params[4] = {20., 300.}; //x_0, lambda
     const Float_t npix_to_qdc_params[4] = {0.172, -1.31, 0.006, 0.33}; // A, B, sigma_A, sigma_B
}

// -----   Default constructor   -------------------------------------------
sndScifiHit::sndScifiHit()
  : SndlhcHit()
{
 flag = true;
}
// -----   Standard constructor   ------------------------------------------
sndScifiHit::sndScifiHit(Int_t detID)
  : SndlhcHit(detID)
{
 flag = true;
}
// -----   constructor from point class  ------------------------------------------
sndScifiHit::sndScifiHit (int SiPMChan, std::vector<ScifiPoint*> V, std::vector<Float_t> W)
{
     Scifi* ScifiDet = dynamic_cast<Scifi*> (gROOT->GetListOfGlobals()->FindObject("Scifi") );
     Float_t nphe_min = ScifiDet->GetConfParF("Scifi/nphe_min");
     Float_t nphe_max = ScifiDet->GetConfParF("Scifi/nphe_max");
     Float_t timeResol = ScifiDet->GetConfParF("Scifi/timeResol");
     Float_t signalSpeed = ScifiDet->GetConfParF("Scifi/signalSpeed");
     
     nSides   = 1;
     nSiPMs   = 1;
     for (unsigned int j=0; j<16; ++j){
        signals[j] = -1;
        times[j]    =-1;
     }

     fDetectorID  = SiPMChan;
     Float_t ly_total = 0;
     Float_t earliestToA   = 1E20;
     for( int i = 0; i <V.size();i++) {
        
        Double_t signal = V[i]->GetEnergyLoss()*W[i];

	// Find distances from MCPoint centre to ends of fibre
        TVector3 a, b;
        TVector3 impact(V[i]->GetX(),V[i]->GetY() ,V[i]->GetZ() );
        ScifiDet->GetSiPMPosition(V[i]->GetDetectorID(),a, b);
	     
	Double_t distance;

	bool verticalHit = int(fDetectorID/100000)%10 == 1;

	// Calculate distance from energy deposit to SiPM.
	// Vertical
	if (verticalHit) {
	  distance = (b - impact).Mag();
	  // Horizontal
	} else {
	  distance = (impact - a).Mag();
	}

	// convert energy deposit to light yield (here Np.e. == avg. N fired pixels)
        Float_t ly = signal*1E+6*0.16; //0.16 p.e per 1 keV
        // account for the light attenuation in the fibers
        ly*= ly_loss(distance);
	ly_total+= ly;
	
	// for the timing, find earliest light to arrive at SiPM and smear with time resolution
	Float_t arrival_time = V[i]->GetTime() + distance/signalSpeed;
	if (arrival_time < earliestToA){earliestToA = arrival_time;}
	
     }
     // smear the total light yield using Poisson distribution
     ly_total= gRandom->Poisson(ly_total);

     // account for limited SiPM dyn. range
     Float_t Npix = sipm_saturation(ly_total,nphe_max);
     // convert Npix to QDC
     signals[0] = npix_to_qdc(Npix);
     
     if (ly_total > nphe_min){   // nominal threshold at 3.5 p.e.
            flag=true;
      }else{
            flag=false;
      }
      times[0] = gRandom->Gaus(earliestToA, timeResol);

     LOG(DEBUG) << "signal created";
}

// -----   Destructor   ----------------------------------------------------
sndScifiHit::~sndScifiHit() { }
// -------------------------------------------------------------------------

// -----   Public method GetEnergy   -------------------------------------------
Float_t sndScifiHit::GetEnergy()
{
  // to be calculated from digis and calibration constants, missing!
  return signals[0];
}

Float_t sndScifiHit::ly_loss(Float_t distance){
//	It returns the light yield attenuation depending on the distance to SiPM
	return TMath::Exp(-(distance-ly_loss_params[0])/ly_loss_params[1]);
}

Float_t sndScifiHit::sipm_saturation(Float_t ly, Float_t nphe_max){
//	It returns the number of fired pixels per channel
        Float_t factor = 1 - TMath::Exp(-ly/nphe_max);
	return nphe_max*factor;
}

Float_t sndScifiHit::npix_to_qdc(Float_t npix){
//	It returns QDC per channel after Gaussian smearing of the parameters
        Float_t A = gRandom->Gaus(npix_to_qdc_params[0], npix_to_qdc_params[2]);
        Float_t B = gRandom->Gaus(npix_to_qdc_params[1], npix_to_qdc_params[3]);
        return A*npix + B;
}

// -----   Public method Print   -------------------------------------------
void sndScifiHit::Print()
{
  std::cout << "-I- sndScifiHit: Scifi hit " << " in station " << GetStation();
  if ( isVertical()) {     std::cout << " vertical plane ";}
  else {     std::cout << " horizontal plane ";}
  std::cout << "SiPM nr "<<GetSiPM()<< " SiPM channel "<<GetSiPMChan()<<std::endl;
}
// -------------------------------------------------------------------------

ClassImp(sndScifiHit)

