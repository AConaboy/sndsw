#include <iostream>
#include <fstream>
#include <string>
#include "TFile.h"
#include "TTree.h"

#ifdef __cplusplus
static_assert(__cplusplus >= 201103L, "C++11 or higher required");
#endif

// Helper function to extract the parent directory name from a file path
std::string ExtractParentDirectoryName(const std::string &filepath) {
    std::size_t lastslash = filepath.find_last_of("/\\");
    if (lastslash == std::string::npos)
        return ""; // No directory found, return empty string

    std::string directory = filepath.substr(0, lastslash);
    std::size_t second_lastslash = directory.find_last_of("/\\");
    if (second_lastslash == std::string::npos)
        return directory; // No parent directory found, return the current directory
    
    return directory.substr(second_lastslash + 1);
}

void ExtractStruckQuarkMomentum(const char* filename) {
    // Convert filename to std::string for easier manipulation
    std::string filepath(filename);

    // Extract directory and base filename
    std::string directory = ExtractParentDirectoryName(filepath);

    // Construct the output filename
    std::string output_filename = "/eos/experiment/sndlhc/users/aconsnd/simulation/neutrino/data/sndlhc_13TeV_down_volTarget_100fb-1_SNDG18_02a_01_000/nueFilter/filteredMC_"+directory+"_quark_momentum.txt";

    // Open the ROOT file containing the GENIE events
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Access the TTree containing the event data
    TTree *tree = (TTree*)file->Get("gst");
    if (!tree) {
        std::cerr << "Error: TTree 'gst' not found in file: " << filename << std::endl;
        file->Close();
        return;
    }

    // Variables to hold momentum components and PDG code
    double vtxx, vtxy, vtxz, px, py, pz;
    int hitqrk_pdg;

    // Set branch addresses
    tree->SetBranchAddress("pxi", &px);
    tree->SetBranchAddress("pyi", &py);
    tree->SetBranchAddress("pzi", &pz);
    tree->SetBranchAddress("vtxx", &vtxx);
    tree->SetBranchAddress("vtxy", &vtxy);
    tree->SetBranchAddress("vtxz", &vtxz);
    tree->SetBranchAddress("hitqrk", &hitqrk_pdg);

    // Open a text file to write the output
    std::ofstream outfile(output_filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening output file: " << output_filename << std::endl;
        file->Close();
        return;
    }

    // Loop over all entries in the tree
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        // Write the event number and momentum vector of the quark to the file
        outfile << i << " " << hitqrk_pdg << " " << vtxx << " " << vtxy << " " << vtxz << " " << px << " " << py << " " << pz << std::endl;
    }

    // Close the output file
    outfile.close();

    if (file) {
        // Clean up and close the ROOT file
        file->Close();
        delete file;
        file=nullptr;
    }

    std::cout << "Momentum vectors for the quark (from hitqrk branch) have been written to " << output_filename << std::endl;
}
