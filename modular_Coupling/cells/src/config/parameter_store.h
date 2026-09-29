#ifndef BDM_CELLS_PARAMETER_STORE_H_
#define BDM_CELLS_PARAMETER_STORE_H_

#include <algorithm>
#include <any>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <vector>

namespace bdm {

/**
 * @brief Enhanced parameter store with type info, printing, and cloning
 *
 * Following ABM4bio pattern for runtime parameter management.
 * Supports type-safe storage, retrieval, printing, and cloning.
 */
class ParameterStore {
 public:
  ParameterStore() = default;
  ParameterStore(const ParameterStore&) = default;
  ParameterStore& operator=(const ParameterStore&) = default;

  // ==========================================================================
  // Existence Checks
  // ==========================================================================

  bool Contains(const std::string& name) const {
    return values_.find(name) != values_.end();
  }

  template <typename T>
  bool Has(const std::string& name) const {
    auto it = values_.find(name);
    if (it == values_.end()) return false;
    return it->second.type() == typeid(T);
  }

  size_t Size() const { return values_.size(); }
  bool Empty() const { return values_.empty(); }

  // ==========================================================================
  // Setters
  // ==========================================================================

  // Set only when the key is absent (used by RegisterOptionalDefaults).
  // If the key was already provided by the CSV it is left unchanged.
  template <typename T>
  void AddDefault(const std::string& name, T value) {
    if (!Has<T>(name)) Set(name, std::move(value));
  }

  template <typename T>
  void Set(const std::string& name, T value) {
    values_[name] = std::move(value);
    type_names_[name] = GetTypeName<T>();
    // Register printer for this type
    if constexpr (std::is_same_v<T, std::vector<std::string>>) {
      printers_[name] = [](const std::any& v) -> std::string {
        const auto& vec = std::any_cast<const std::vector<std::string>&>(v);
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < vec.size(); ++i) {
          if (i > 0) oss << ", ";
          oss << vec[i];
        }
        oss << "]";
        return oss.str();
      };
    } else {
      printers_[name] = [](const std::any& v) -> std::string {
        std::ostringstream oss;
        oss << std::any_cast<const T&>(v);
        return oss.str();
      };
    }
  }

  // ==========================================================================
  // Getters
  // ==========================================================================

  template <typename T>
  const T& Get(const std::string& name) const {
    auto it = values_.find(name);
    if (it == values_.end()) {
      throw std::runtime_error("Missing parameter: " + name);
    }
    if (it->second.type() != typeid(T)) {
      throw std::runtime_error("Type mismatch for parameter '" + name +
                               "': expected " + GetTypeName<T>() +
                               ", got " + GetStoredTypeName(name));
    }
    return std::any_cast<const T&>(it->second);
  }

  template <typename T>
  T GetOr(const std::string& name, T default_value) const {
    auto it = values_.find(name);
    if (it == values_.end()) return default_value;
    if (it->second.type() != typeid(T)) {
      throw std::runtime_error("Type mismatch for parameter: " + name);
    }
    return std::any_cast<const T&>(it->second);
  }

  // ==========================================================================
  // Type Information
  // ==========================================================================

  std::string GetStoredTypeName(const std::string& name) const {
    auto it = type_names_.find(name);
    return (it != type_names_.end()) ? it->second : "unknown";
  }

  template <typename T>
  static std::string GetTypeName() {
    if constexpr (std::is_same_v<T, double>) return "double";
    else if constexpr (std::is_same_v<T, int>) return "int";
    else if constexpr (std::is_same_v<T, bool>) return "bool";
    else if constexpr (std::is_same_v<T, std::string>) return "string";
    else if constexpr (std::is_same_v<T, std::vector<std::string>>) return "string_list";
    else return typeid(T).name();
  }

  // ==========================================================================
  // Printing
  // ==========================================================================

  std::string ValueToString(const std::string& name) const {
    auto it = printers_.find(name);
    if (it == printers_.end()) return "<unprintable>";
    auto val_it = values_.find(name);
    if (val_it == values_.end()) return "<missing>";
    try {
      return it->second(val_it->second);
    } catch (...) {
      return "<error>";
    }
  }

  void PrintAll(std::ostream& os = std::cout) const {
    os << "ParameterStore (" << values_.size() << " parameters):\n";
    std::vector<std::string> keys;
    for (const auto& [k, v] : values_) keys.push_back(k);
    std::sort(keys.begin(), keys.end());

    for (const auto& key : keys) {
      os << "  " << key << " [" << GetStoredTypeName(key) << "] = "
         << ValueToString(key) << "\n";
    }
  }

  void PrintByPrefix(const std::string& prefix, std::ostream& os = std::cout) const {
    os << "Parameters with prefix '" << prefix << "':\n";
    for (const auto& [key, val] : values_) {
      if (key.rfind(prefix, 0) == 0) {
        os << "  " << key << " [" << GetStoredTypeName(key) << "] = "
           << ValueToString(key) << "\n";
      }
    }
  }

  // ==========================================================================
  // Cloning
  // ==========================================================================

  ParameterStore Clone() const {
    ParameterStore copy;
    copy.values_ = values_;
    copy.type_names_ = type_names_;
    copy.printers_ = printers_;
    return copy;
  }

  // ==========================================================================
  // Iteration
  // ==========================================================================

  std::vector<std::string> GetAllKeys() const {
    std::vector<std::string> keys;
    keys.reserve(values_.size());
    for (const auto& [k, v] : values_) keys.push_back(k);
    return keys;
  }

  std::vector<std::string> GetKeysByPrefix(const std::string& prefix) const {
    std::vector<std::string> keys;
    for (const auto& [k, v] : values_) {
      if (k.rfind(prefix, 0) == 0) keys.push_back(k);
    }
    return keys;
  }

  // ==========================================================================
  // Utilities
  // ==========================================================================

  void Clear() {
    values_.clear();
    type_names_.clear();
    printers_.clear();
  }

  void Remove(const std::string& name) {
    values_.erase(name);
    type_names_.erase(name);
    printers_.erase(name);
  }

 private:
  std::unordered_map<std::string, std::any> values_;
  std::unordered_map<std::string, std::string> type_names_;
  std::unordered_map<std::string, std::function<std::string(const std::any&)>> printers_;
};

}  // namespace bdm

#endif  // BDM_CELLS_PARAMETER_STORE_H_
