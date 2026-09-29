// =============================================================================
// cells.cc - Main entry point for OpenFOAM-preCICE-BioDynaMo coupling
// Following ABM4bio pattern: minimal main, logic in headers
// =============================================================================

#include "cells.h"
#include "config/parameter_loader.h"

int main(int argc, const char** argv) {
    // Load runtime parameters from CSV (follows ABM4bio pattern)
    std::string params_file = "params.csv";
    if (argc >= 2) {
        params_file = argv[1];
    }
    bdm::LoadParametersFromCsv(params_file);

    // Run simulation (mode selected via BDM_SIMULATION_MODE env var)
    return bdm::Simulate(argc, argv);
}
