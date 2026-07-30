#include "vtk.cuh"

#include "./presets/geometry.h"
#include "./functions/grid_plot.cuh"

#include "cuda_config/var.cuh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>

__host__ void write_vti(int step, const std::string &out_dir, varUnit *mom_host)
{
    namespace fs = std::filesystem;

    fs::path vtk_dir = fs::path(out_dir) / "vtk";
    fs::create_directories(vtk_dir);

    std::ostringstream filename;
    filename << "output_" << std::setw(6) << std::setfill('0') << step << ".vti";

    fs::path filepath = vtk_dir / filename.str();

    std::ofstream file(filepath.string());
    if (!file.is_open())
    {
        std::cerr << "Could not open VTI file for writing: " << filepath.string() << "\n";
        return;
    }

    const int nx = Nx;
    const int ny = Ny;

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <ImageData WholeExtent=\"0 " << (nx - 1)
         << " 0 " << (ny - 1)
         << " 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    file << "    <Piece Extent=\"0 " << (nx - 1)
         << " 0 " << (ny - 1)
         << " 0 0\">\n";

    file << "      <PointData>\n";

    // rho
    file << "        <DataArray type=\"Float32\" Name=\"rho\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y)
    {
        for (int x = 0; x < nx; ++x)
        {
            const size_t idx = grid_plot(x, y);
            const varUnit rho = mom_host[momIdx<MomentId::rho>(idx)]; // rho
            file << "          " << static_cast<varUnit>(rho) << "\n";
        }
    }
    file << "        </DataArray>\n";

    // ux
    file << "        <DataArray type=\"Float32\" Name=\"ux\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y)
    {
        for (int x = 0; x < nx; ++x)
        {
            const size_t idx = grid_plot(x, y);
            const varUnit ux = mom_host[momIdx<MomentId::ux>(idx)]; // ux
            file << "          " << static_cast<varUnit>(ux) << "\n";
        }
    }
    file << "        </DataArray>\n";

    // uy
    file << "        <DataArray type=\"Float32\" Name=\"uy\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y)
    {
        for (int x = 0; x < nx; ++x)
        {
            const size_t idx = grid_plot(x, y);
            const varUnit uy = mom_host[momIdx<MomentId::uy>(idx)]; // uy
            file << "          " << static_cast<varUnit>(uy) << "\n";
        }
    }
    file << "        </DataArray>\n";

    // velocity
    file << "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int y = 0; y < ny; ++y)
    {
        for (int x = 0; x < nx; ++x)
        {
            const size_t idx = grid_plot(x, y);
            const varUnit ux = mom_host[momIdx<MomentId::ux>(idx)]; // ux
            const varUnit uy = mom_host[momIdx<MomentId::uy>(idx)]; // uy

            file << "          "
                 << static_cast<varUnit>(ux) << " "
                 << static_cast<varUnit>(uy) << " "
                 << 0.0f << "\n";
        }
    }
    file << "        </DataArray>\n";

    file << "      </PointData>\n";
    file << "      <CellData>\n";
    file << "      </CellData>\n";
    file << "    </Piece>\n";
    file << "  </ImageData>\n";
    file << "</VTKFile>\n";

    file.close();
}
