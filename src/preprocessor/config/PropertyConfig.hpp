#pragma once
#include <vector>
#include <string>

enum class VariableType
{
  flow,
  mechanics,
  service
};

struct DomainConfig
{
  int label;
  std::vector<std::string> expressions;
  std::vector<std::string> variables;

  std::vector< size_t > labels;
  std::unordered_map< size_t, std::vector<std::string> > tagged_vars;
  std::unordered_map< size_t, std::vector<double> > tagged_expr;

  VariableType type;
  bool coupled = true;  // whether to generate mapping
};


struct CellPropertyConfig
{
  std::vector<DomainConfig> domains;
  DomainConfig files;
  std::string x_kwd = "x";
  std::string y_kwd = "y";
  std::string z_kwd = "z";
  // are required for flow discretization
  std::string perm_kwd = "PERM";
  std::string permx_kwd = "PERMX";
  std::string permy_kwd = "PERMY";
  std::string permz_kwd = "PERMZ";
  std::string poro_kwd = "PORO";
  std::string vmult_kwd = "VFACTOR";
  std::vector<std::string> service_variables{x_kwd, y_kwd, z_kwd, vmult_kwd};
  std::vector<std::string> extra_variables{perm_kwd, permx_kwd, permy_kwd, permz_kwd, poro_kwd};
  static constexpr double default_volume_factor = 1.f;
};
