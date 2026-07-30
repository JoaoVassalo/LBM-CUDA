#pragma once

#include "presets/config.h"

#include <string>

__host__ void write_vti(int step, const std::string &out_dir, varUnit *mom_host);
