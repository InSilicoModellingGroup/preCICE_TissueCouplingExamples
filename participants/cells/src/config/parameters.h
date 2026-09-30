#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>
#include <vector>

#include "parameter_store.h"

namespace bdm {

struct PreciceAdapterParameters {
  struct PreciceConfig {
    std::string participant_name;
    std::string config_file;
    bool enable_profiling;
    std::string profiling_variable;
    std::string profiling_value;
    std::string profiling_mode_variable;
    std::string profiling_mode_value;
    std::string velocity_data_name;
    std::string pressure_data_name;
    std::string pressure_gradient_data_name;  // ∇p from OpenFOAM FV
    std::string temperature_data_name;
    std::string position_data_name;
    std::string wall_geometry_data_name;
    std::string wall_normals_data_name;
    std::string wall_areas_data_name;
    std::string wall_velocity_data_name;  // For moving walls
    std::string fluid_density_data_name;
    std::string fluid_dynamic_viscosity_data_name;
  } precice;

  struct MeshConfig {
    std::string particle_mesh;
    std::string primary_volume_mesh;
    std::vector<std::string> default_volume_meshes;
    std::vector<std::string> default_wall_meshes;
    std::string foam_case_directory;
    std::string block_mesh_dict;
  } mesh;

  struct AgentConfig {
    double default_diameter;
    double default_mass;
    std::string custom_positions_csv;
    bool custom_positions_required;
  } agent;

  struct EnvironmentConfig {
    std::string volume_meshes_variable;
    std::string primary_volume_variable;
    std::string wall_meshes_variable;
  } environment;
};

const PreciceAdapterParameters& GetPreciceAdapterParameters();

void SetPreciceAdapterParameters(const PreciceAdapterParameters& params);

} // namespace bdm

#endif
