#include "parameter_loader.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "biodynamo.h"
#include "parameters.h"

namespace bdm {
namespace {

ParameterStore g_parameter_store;

std::string Trim(const std::string& value) {
  auto left = value.begin();
  while (left != value.end() && std::isspace(static_cast<unsigned char>(*left))) {
    ++left;
  }
  auto right = value.end();
  do {
    --right;
  } while (right >= left && std::isspace(static_cast<unsigned char>(*right)));
  if (right < left) {
    return "";
  }
  return std::string(left, right + 1);
}

bool ParseBool(const std::string& raw, bool* result) {
  std::string value = raw;
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  if (value == "true" || value == "1") {
    *result = true;
    return true;
  }
  if (value == "false" || value == "0") {
    *result = false;
    return true;
  }
  return false;
}

std::vector<std::string> SplitList(const std::string& raw) {
  std::vector<std::string> entries;
  std::stringstream ss(raw);
  std::string item;
  while (std::getline(ss, item, ';')) {
    auto trimmed = Trim(item);
    if (!trimmed.empty()) {
      entries.push_back(trimmed);
    }
  }
  return entries;
}

void ApplyPreciceOverrides() {
  auto params = GetPreciceAdapterParameters();
  PreciceAdapterParameters updated = params;
  const auto& store = g_parameter_store;

  auto maybe_string = [&](const std::string& name, std::string* target) {
    if (store.Has<std::string>(name)) {
      *target = store.Get<std::string>(name);
    }
  };

  auto maybe_bool = [&](const std::string& name, bool* target) {
    if (store.Has<bool>(name)) {
      *target = store.Get<bool>(name);
    }
  };

  auto maybe_double = [&](const std::string& name, double* target) {
    if (store.Has<double>(name)) {
      *target = store.Get<double>(name);
    }
  };

  auto maybe_string_list = [&](const std::string& name, std::vector<std::string>* target) {
    if (store.Has<std::vector<std::string>>(name)) {
      *target = store.Get<std::vector<std::string>>(name);
    }
  };

  // preCICE configuration
  maybe_string("precice.participant_name", &updated.precice.participant_name);
  maybe_string("precice.config_file", &updated.precice.config_file);
  maybe_bool("precice.enable_profiling", &updated.precice.enable_profiling);
  maybe_string("precice.profiling_variable", &updated.precice.profiling_variable);
  maybe_string("precice.profiling_value", &updated.precice.profiling_value);
  maybe_string("precice.profiling_mode_variable", &updated.precice.profiling_mode_variable);
  maybe_string("precice.profiling_mode_value", &updated.precice.profiling_mode_value);
  maybe_string("precice.velocity_data_name", &updated.precice.velocity_data_name);
  maybe_string("precice.pressure_data_name", &updated.precice.pressure_data_name);
  maybe_string("precice.pressure_gradient_data_name", &updated.precice.pressure_gradient_data_name);
  maybe_string("precice.temperature_data_name", &updated.precice.temperature_data_name);
  maybe_string("precice.position_data_name", &updated.precice.position_data_name);
  maybe_string("precice.wall_geometry_data_name", &updated.precice.wall_geometry_data_name);
  maybe_string("precice.wall_normals_data_name", &updated.precice.wall_normals_data_name);
  maybe_string("precice.wall_areas_data_name", &updated.precice.wall_areas_data_name);
  maybe_string("precice.wall_velocity_data_name", &updated.precice.wall_velocity_data_name);
  maybe_string("precice.fluid_density_data_name", &updated.precice.fluid_density_data_name);
  maybe_string("precice.fluid_dynamic_viscosity_data_name", &updated.precice.fluid_dynamic_viscosity_data_name);

  // Mesh configuration
  maybe_string("mesh.particle_mesh", &updated.mesh.particle_mesh);
  maybe_string("mesh.primary_volume_mesh", &updated.mesh.primary_volume_mesh);
  maybe_string("mesh.foam_case_directory", &updated.mesh.foam_case_directory);
  maybe_string("mesh.block_mesh_dict", &updated.mesh.block_mesh_dict);
  maybe_string_list("mesh.default_volume_meshes", &updated.mesh.default_volume_meshes);
  maybe_string_list("mesh.default_wall_meshes", &updated.mesh.default_wall_meshes);

  // Agent configuration
  maybe_double("agent.default_diameter", &updated.agent.default_diameter);
  maybe_double("agent.default_mass", &updated.agent.default_mass);
  maybe_string("agent.custom_positions_csv", &updated.agent.custom_positions_csv);
  maybe_bool("agent.custom_positions_required", &updated.agent.custom_positions_required);

  // Environment configuration
  maybe_string("environment.volume_meshes_variable", &updated.environment.volume_meshes_variable);
  maybe_string("environment.primary_volume_variable", &updated.environment.primary_volume_variable);
  maybe_string("environment.wall_meshes_variable", &updated.environment.wall_meshes_variable);

  SetPreciceAdapterParameters(updated);
}

// =============================================================================
// Optional parameter defaults
//
// Every entry here makes the corresponding CSV key *optional*.
// CSV values always win — AddDefault() is a no-op when the key is already set.
// All other parameters remain strictly required (GetParam<T>(key) will throw
// if they are absent from the CSV).
//
// HOW TO ADD A NEW OPTIONAL PARAMETER
//   1. Append one line below with the canonical default.
//   2. Call sites may keep using GetParam<T>(key) — no default needed there.
// =============================================================================
void RegisterOptionalDefaults(ParameterStore& store) {
  auto d = [&](const std::string& k, double v)      { store.AddDefault(k, v); };
  auto b = [&](const std::string& k, bool   v)      { store.AddDefault(k, v); };
  (void)b;  // suppress unused-lambda warning when no bool defaults are listed

  // --- Interaction / contact mechanics ---
  // Particle-particle mechanics and optional rebound coefficients. Geometric
  // wall confinement itself is mass-independent.
  d("interaction.youngs_modulus",            1.0e-5);
  d("interaction.poisson_ratio",             0.30);
  d("interaction.restitution_coefficient",   0.1);
  d("interaction.friction_coefficient",      0.0);
  d("interaction.cohesion_energy",           0.0);
  d("interaction.normal_damping",            0.02);
  d("interaction.tangential_damping",        0.02);
  d("interaction.damping_factor",            0.02);
  d("interaction.interaction_range",         2.0e-3);
  d("interaction.force_ramp_time",           0.0);
  d("interaction.wall_youngs_modulus",       1.0e-3);
  d("interaction.wall_poisson_ratio",        0.30);
  d("interaction.wall_friction_coefficient", 0.0);
  d("interaction.wall_restitution",          0.1);

  // --- Fluid properties ---
  // Required only when flow coupling (velocity/pressure/drag) is active.
  // Water at ~20 °C is used as a safe placeholder; override in params.csv for
  // blood (1060 / 0.0035) or any other fluid.
  d("fluid.density",           1000.0);
  d("fluid.dynamic_viscosity",    0.001);
}

}  // namespace

const ParameterStore& GetParameterStore() { return g_parameter_store; }

bool LoadParametersFromCsv(const std::string& filepath) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    Log::Warning("ParameterLoader", "Unable to open parameter file '", filepath, "'. Using defaults.");
    return false;
  }

  g_parameter_store.Clear();

  std::string header;
  std::getline(file, header);  // Skip header line

  std::string line;
  size_t line_no = 1;
  while (std::getline(file, line)) {
    ++line_no;
    auto trimmed = Trim(line);
    if (trimmed.empty() || trimmed[0] == '#') {
      continue;
    }

    std::vector<std::string> columns;
    std::stringstream ss(trimmed);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
      columns.push_back(Trim(cell));
    }
    if (columns.size() < 3) {
      Log::Warning("ParameterLoader", "Line ", line_no, ": expected 3 columns, found ", columns.size());
      continue;
    }

    const std::string& name = columns[0];
    const std::string& type = columns[1];
    const std::string& value = columns[2];

    if (name.empty()) {
      continue;
    }

    try {
      if (type == "string") {
        g_parameter_store.Set(name, value);
      } else if (type == "bool") {
        bool parsed = false;
        if (!ParseBool(value, &parsed)) {
          throw std::runtime_error("invalid boolean value: " + value);
        }
        g_parameter_store.Set(name, parsed);
      } else if (type == "int") {
        g_parameter_store.Set(name, std::stoi(value));
      } else if (type == "double" || type == "float") {
        g_parameter_store.Set(name, std::stod(value));
      } else if (type == "string_list") {
        g_parameter_store.Set(name, SplitList(value));
      } else {
        Log::Warning("ParameterLoader", "Line ", line_no, ": unsupported type '", type, "'");
      }
    } catch (const std::exception& e) {
      Log::Warning("ParameterLoader", "Line ", line_no, ": ", e.what());
    }
  }

  ApplyPreciceOverrides();
  RegisterOptionalDefaults(g_parameter_store);
  Log::Info("ParameterLoader", "Loaded ", g_parameter_store.Size(), " parameters from ", filepath);

  // Print all parameters if logging.print_parameters is true
  if (g_parameter_store.GetOr<bool>("logging.print_parameters", false)) {
    std::cout << "\n";
    g_parameter_store.PrintAll(std::cout);
    std::cout << std::endl;
  }

  return true;
}

void PrintParametersByPrefix(const std::string& prefix) {
  GetParameterStore().PrintByPrefix(prefix, std::cout);
}

}  // namespace bdm
