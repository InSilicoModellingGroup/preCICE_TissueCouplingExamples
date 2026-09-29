#ifndef BDM_CELLS_PARAMETER_LOADER_H_
#define BDM_CELLS_PARAMETER_LOADER_H_

#include <string>

#include "parameter_store.h"

namespace bdm {

//! Loads parameters from a CSV file (parameter_name,type,parameter_value).
//! Returns true on success. Missing files are tolerated (defaults remain).
bool LoadParametersFromCsv(const std::string& filepath);

//! Returns the global parameter store populated by LoadParametersFromCsv.
const ParameterStore& GetParameterStore();

//! Print parameters filtered by prefix (e.g., "fluid.", "simulation.")
void PrintParametersByPrefix(const std::string& prefix);

}  // namespace bdm

#endif  // BDM_CELLS_PARAMETER_LOADER_H_
