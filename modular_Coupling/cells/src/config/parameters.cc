#include "parameters.h"
#include <limits>

namespace bdm {
namespace {

PreciceAdapterParameters::PreciceConfig DefaultPreciceConfig() {
  return {
      /*participant_name=*/"cells",
      /*config_file=*/"../precice-config.xml",
      /*enable_profiling=*/true,
      /*profiling_variable=*/"PRECICE_PROFILING",
      /*profiling_value=*/"1",
      /*profiling_mode_variable=*/"PRECICE_PROFILING_MODE",
      /*profiling_mode_value=*/"all",
      /*velocity_data_name=*/"U",
      /*pressure_data_name=*/"p",
      /*pressure_gradient_data_name=*/"gradP",
      /*temperature_data_name=*/"T",
      /*position_data_name=*/"Position",
      /*wall_geometry_data_name=*/"WallGeometry",
      /*wall_normals_data_name=*/"WallNormals",
      /*wall_areas_data_name=*/"WallAreas",
      /*wall_velocity_data_name=*/"WallVelocity",
      /*fluid_density_data_name=*/"rho",
      /*fluid_dynamic_viscosity_data_name=*/"mu"};
}

PreciceAdapterParameters::MeshConfig DefaultMeshConfig() {
  return {
      /*particle_mesh=*/"ABMmesh",
      /*primary_volume_mesh=*/"FVMmesh",
      /*default_volume_meshes=*/{"FVMmesh"},
      /*default_wall_meshes=*/{"WallSurfaceMesh"},
      /*foam_case_directory=*/"../participants/simple",
      /*block_mesh_dict=*/"system/blockMeshDict"};
}

PreciceAdapterParameters::AgentConfig DefaultAgentConfig() {
  return {
      /*default_diameter=*/0.004,
      /*default_mass=*/std::numeric_limits<double>::quiet_NaN(),
      /*custom_positions_csv=*/"agents_positions.csv",
      /*custom_positions_required=*/false};
}

PreciceAdapterParameters::EnvironmentConfig DefaultEnvironmentConfig() {
  return {
      /*volume_meshes_variable=*/"BDM_VOLUME_MESHES",
      /*primary_volume_variable=*/"BDM_PRIMARY_VOLUME_MESH",
      /*wall_meshes_variable=*/"BDM_WALL_MESHES"};
}

PreciceAdapterParameters g_precice_parameters{DefaultPreciceConfig(),
                                              DefaultMeshConfig(),
                                              DefaultAgentConfig(),
                                              DefaultEnvironmentConfig()};

}  // namespace

const PreciceAdapterParameters& GetPreciceAdapterParameters() {
  return g_precice_parameters;
}

void SetPreciceAdapterParameters(const PreciceAdapterParameters& params) {
  g_precice_parameters = params;
}

}  // namespace bdm
