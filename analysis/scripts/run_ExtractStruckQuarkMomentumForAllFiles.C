#include <iostream>
#include <fstream>
#include <string>
#include "TFile.h"
#include "TTree.h"

#ifdef __cplusplus
static_assert(__cplusplus >= 201103L, "C++11 or higher required");
#endif

// Forward declaration of the existing macro
void ExtractStruckQuarkMomentum(const char* filename);

void run_ExtractStruckQuarkMomentumForAllFiles() {

    // Open the file containing the list of input file paths
    const char* list_filename = "/afs/cern.ch/user/a/aconsnd/Timing/condor_scripts/simlist.txt";
    std::ifstream infile(list_filename);
    std::cout << "File loaded as ifstream" << std::endl;
    if (!infile.is_open()) {
        std::cerr << "Error opening list file: " << list_filename << std::endl;
        return;
    }

    std::string input_file;
    // Read each line (input file path) from the list file
    while (std::getline(infile, input_file)) {
        // Call the existing macro for each input file
        std::cout << "Processing file: " << input_file << std::endl;
        ExtractStruckQuarkMomentum(input_file.c_str());
    }

    // Close the list file
    infile.close();

    std::cout << "Finished processing all files." << std::endl;
}
