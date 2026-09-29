#ifndef AGENT_TRAJECTORY_LOGGER_H_
#define AGENT_TRAJECTORY_LOGGER_H_

// System includes MUST come before biodynamo.h to avoid rootcling namespace
// aliasing issues (BioDynaMo creates bdm::std which confuses rootcling).
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cerrno>
#include <dirent.h>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <sys/stat.h>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

#include "biodynamo.h"
#include "../agents/coupled_cell.h"

namespace bdm {

class AgentTrajectoryLogger {
 public:
  struct AgentSample {
    uint64_t id;
    Real3 position;
  };

  AgentTrajectoryLogger() {
    const auto interval = ParseEnvDouble("BDM_AGENT_EXPORT_INTERVAL");
    if (!interval || *interval <= 0.0) {
      return;
    }

    export_interval_ = *interval;

    const auto start_time = ParseEnvDouble("BDM_COUPLING_START_TIME");
    if (start_time && *start_time >= 0.0) {
      start_time_ = *start_time;
    }

    output_dir_ = ResolveOutputDir();
    PrepareOutputDir();

    aggregate_path_ = JoinPath(output_dir_, "agents_all_times.csv");

    double first_snapshot_time = start_time_;
    if (export_interval_ > 0.0) {
      const double scaled = start_time_ / export_interval_;
      double multiple = std::ceil(scaled - 1e-12);
      if (multiple < 0.0) {
        multiple = 0.0;
      }
      first_snapshot_time = multiple * export_interval_;
    }

    next_write_time_ = first_snapshot_time;
    enabled_ = true;

    Log::Info("AgentTrajectoryLogger", "Enabled agent logging every ",
              export_interval_, " s starting at t=", first_snapshot_time,
              " (requested start t=", start_time_, ", output: ",
              output_dir_, ")");
  }

  bool Enabled() const { return enabled_; }
  double StartTime() const { return start_time_; }
  double ExportInterval() const { return export_interval_; }

  void MaybeWrite(double current_time, Simulation& simulation) {
    if (!enabled_) {
      return;
    }

    constexpr double kTolerance = 1e-12;
    while (current_time + kTolerance >= next_write_time_) {
      WriteSnapshot(next_write_time_, simulation);
      next_write_time_ += export_interval_;
    }
  }

  void Finalize() {
    if (aggregate_stream_.is_open()) {
      aggregate_stream_.flush();
      aggregate_stream_.close();
    }
  }

 private:
  bool enabled_ = false;
  double start_time_ = 0.0;
  double export_interval_ = 0.0;
  double next_write_time_ = 0.0;
  std::string output_dir_;
  std::string aggregate_path_;
  std::ofstream aggregate_stream_;
  bool aggregate_ready_ = false;
  size_t snapshots_written_ = 0;

  static std::string JoinPath(const std::string& lhs, const std::string& rhs) {
    if (lhs.empty()) {
      return rhs;
    }
    if (lhs.back() == '/') {
      return lhs + rhs;
    }
    return lhs + "/" + rhs;
  }

  static bool IsDirectory(const std::string& path) {
    struct stat st {};
    return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
  }

  static void EnsureDirectoryRecursive(const std::string& path) {
    if (path.empty()) {
      return;
    }

    std::string partial;
    partial.reserve(path.size());

    if (path[0] == '/') {
      partial = "/";
    }

    size_t start = (path[0] == '/') ? 1 : 0;
    while (start <= path.size()) {
      size_t slash = path.find('/', start);
      const size_t end = (slash == std::string::npos) ? path.size() : slash;
      const std::string component = path.substr(start, end - start);

      if (!component.empty()) {
        partial = partial.empty() ? component : JoinPath(partial, component);
        if (!IsDirectory(partial)) {
          if (mkdir(partial.c_str(), 0775) != 0 && errno != EEXIST) {
            throw std::runtime_error("Failed to create directory: " + partial);
          }
        }
      }

      if (slash == std::string::npos) {
        break;
      }
      start = slash + 1;
    }
  }

  static bool EndsWithCsv(const std::string& name) {
    return name.size() >= 4 && name.compare(name.size() - 4, 4, ".csv") == 0;
  }

  static std::optional<double> ParseEnvDouble(const char* name) {
    if (const char* value = std::getenv(name)) {
      try {
        return std::stod(value);
      } catch (...) {
        Log::Warning("AgentTrajectoryLogger", "Unable to parse ", name,
                     "='", value, "'.");
      }
    }
    return std::nullopt;
  }

  static std::string ResolveOutputDir() {
    if (const char* value = std::getenv("BDM_AGENT_OUTPUT_DIR")) {
      return std::string(value);
    }
    return "output/agent_history";
  }

  static std::string FormatTime(double value) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(6) << value;
    std::string token = oss.str();
    if (token.find('.') != std::string::npos) {
      while (!token.empty() && token.back() == '0') {
        token.pop_back();
      }
      if (!token.empty() && token.back() == '.') {
        token.pop_back();
      }
    }
    if (token.empty()) {
      token = "0";
    }
    return token;
  }

  void PrepareOutputDir() {
    EnsureDirectoryRecursive(output_dir_);

    DIR* dir = opendir(output_dir_.c_str());
    if (dir == nullptr) {
      throw std::runtime_error("Failed to open agent output directory: " + output_dir_);
    }

    while (dirent* entry = readdir(dir)) {
      const std::string name(entry->d_name);
      if (name == "." || name == ".." || !EndsWithCsv(name)) {
        continue;
      }

      const std::string full_path = JoinPath(output_dir_, name);
      struct stat st {};
      if (stat(full_path.c_str(), &st) == 0 && S_ISREG(st.st_mode)) {
        unlink(full_path.c_str());
      }
    }

    closedir(dir);
  }

  void EnsureAggregateReady() {
    if (aggregate_ready_) {
      return;
    }
    aggregate_stream_.open(aggregate_path_);
    if (!aggregate_stream_) {
      throw std::runtime_error("Unable to open " + aggregate_path_ +
                               " for writing.");
    }
    aggregate_stream_ << "ID,time,x,y,z\n";
    aggregate_ready_ = true;
  }

  void WriteSnapshot(double target_time, Simulation& simulation) {
    auto* rm = simulation.GetResourceManager();
    std::vector<AgentSample> samples;
    samples.reserve(rm->GetNumAgents());

    rm->ForEachAgent([&](Agent* agent) {
      if (auto* cell = dynamic_cast<MyCell*>(agent)) {
        samples.push_back(
            AgentSample{cell->GetUid().GetIndex(), cell->GetPosition()});
      }
      return true;
    });

    if (samples.empty()) {
      return;
    }

    const std::string time_token = FormatTime(target_time);
    const auto snapshot_path = JoinPath(output_dir_, time_token + ".csv");
    WriteCsv(snapshot_path, time_token, samples);
    AppendAggregateRows(time_token, samples);
    ++snapshots_written_;

    Log::Info("AgentTrajectoryLogger", "Wrote agent snapshot ", snapshots_written_,
              " at t=", time_token);
  }

  static void WriteCsv(const std::string& path,
                       const std::string& time_token,
                       const std::vector<AgentSample>& samples) {
    std::ofstream stream(path);
    if (!stream) {
      throw std::runtime_error("Failed to write agent snapshot " +
                               path);
    }

    stream << "ID,time,x,y,z\n";
    stream.setf(std::ios::fixed);
    stream << std::setprecision(9);
    for (const auto& sample : samples) {
      stream << sample.id << ',' << time_token << ',' << sample.position[0]
             << ',' << sample.position[1] << ',' << sample.position[2] << '\n';
    }
  }

  void AppendAggregateRows(const std::string& time_token,
                           const std::vector<AgentSample>& samples) {
    EnsureAggregateReady();
    aggregate_stream_.setf(std::ios::fixed);
    aggregate_stream_ << std::setprecision(9);
    for (const auto& sample : samples) {
      aggregate_stream_ << sample.id << ',' << time_token << ','
                        << sample.position[0] << ',' << sample.position[1]
                        << ',' << sample.position[2] << '\n';
    }
  }
};

}  // namespace bdm

#endif  // AGENT_TRAJECTORY_LOGGER_H_
