//
// Prints the events in the input GHEP event file
//
// Usage:
// shell% genie
// genie[0] .x ghep_dump.C("/path/to/my/genie_file.root");
//
// Costas Andreopoulos, Dec 13, 2006
//

void inv_mass(const char * filename, const char * outfilename)
{
 int pdg;
 // open the ROOT file and get the event TTree 
 TFile inpfile(filename,"READ");
 TTree * tree = dynamic_cast <TTree *> (inpfile.Get("gtree"));
 FILE *fp1;
 fp1 = fopen(outfilename, "w");
 // set the branch address
 genie::NtpMCEventRecord * mcrec = 0;
 tree->SetBranchAddress("gmcrec", &mcrec);

 // define some particle code

 const int kPdgNuE          =   12;
 const int kPdgAntiNuE      =  -12;
 const int kPdgNuMu         =   14;
 const int kPdgAntiNuMu     =  -14;
 const int kPdgNuTau        =   16;
 const int kPdgAntiNuTau    =  -16;
 const int kPdgGamma        =    22; // photon
 const int kPdgProton       =  2212;
 const int kPdgAntiProton   = -2212;
 const int kPdgNeutron      =  2112;
 const int kPdgAntiNeutron  = -2112;
 const int kPdgPiP          =   211; // pi+
 const int kPdgPiM          =  -211; // pi-
 const int kPdgPi0          =   111; // pi0
 const int kPdgEta          =   221; // eta
 const int kPdgKP           =   321; // K+
 const int kPdgKM           =  -321; // K-
 const int kPdgK0           =   311; // K0
 const int kPdgAntiK0       =  -311; // \bar{K0}
 const int kPdgLambda       =  3122; // Lambda
 const int kPdgAntiLambda   = -3122; // \bar{Lambda}
 //
double Px = 0;
double Py = 0;
double Pz = 0;
double E = 0;
double Minv = 0;
 // loop over event tree 
 for(int i = 0; i< tree->GetEntries(); i++) {
    tree->GetEntry(i);

    genie::NtpMCRecHeader rec_header = mcrec->hdr;
    genie::EventRecord &  event      = *(mcrec->event);
    genie::GHepParticle* p = 0;
    TObjArrayIter event_iter(&event);
   // print-out
    //cout << rec_header;
    if(!(event.HitNucleon())) continue;
    //if(event.Summary()->ProcInfo().IsWeakCC()) continue;
    //if(event.Summary()->ProcInfo().IsDeepInelastic())
    //{
    /*
    Px = 0;
    Py = 0;
    Pz = 0;
    E = 0;
    Minv = 0;
    */
	 genie::Interaction* in = event.Summary();
         const genie::ProcessInfo & proc  = in->ProcInfo();
        const genie::Kinematics  & kine = in->Kine();
        bool selected = true;
        double Q2s = kine.Q2(selected);
        double Ws = kine.W(selected);
        double x = kine.x(selected);
        double nu = (-193.42*193.42 + Q2s + Ws*Ws)/(2*193.42);
        //if(nu < 2) continue;
        //if(pdg == 211 || pdg == -211 || pdg == 111 || pdg == 221)
        //p->X4()->Print();
        //if(p->Status() != genie::kIStStableFinalState) continue;
            //E += p->E();
            //Px += p->Px();
            //Py += p->Py();
            //Pz += p->Pz();
        int hit_id = 0;
    while((p=dynamic_cast<genie::GHepParticle *>(event_iter.Next())))
{
	pdg = p->Pdg();
	//if(p->Status() != genie::kIStStableFinalState) continue;
	    //E += p->E();
	    //Px += p->Px();
	    //Py += p->Py();
	    //Pz += p->Pz();
	    fprintf(fp1,"%d  %d  %d  %d  %s  %6.5e  %6.5e  %6.5e %6.5e  %d  %6.5e  %6.5e  %6.5e \n", i, hit_id, (p->Pdg()), p->FirstMother(), p->Name().c_str(), p->Px(), p->Py(), p->Pz(), p->E(), p->Status() == genie::kIStStableFinalState, p->Vx(), p->Vy(), p->Vz());
	hit_id++;

}
    //fprintf(fp1, "%d \n", -88888);
    //cout << (event.Summary()) << "\n";
    //cout << event.E << 
    mcrec->Clear();
    // }
 }

 inpfile.Close();
}

