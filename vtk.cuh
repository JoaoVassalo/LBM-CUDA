#pragma once

#include <string>

__host__ void write_vti(int step, const std::string &out_dir, float *rho_host, float *ux_host, float *uy_host);
