#pragma once

#include <string>
#include <fstream>

namespace output
{
    std::string vtkPath = "./plot";
    std::ofstream tkePath("io/out/tke.csv");
}