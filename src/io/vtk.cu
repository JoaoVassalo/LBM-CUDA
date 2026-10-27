#include "vtk.cuh"

#include <cstdint>

namespace
{
void writeUInt64(std::ofstream &file, std::uint64_t value)
{
    unsigned char bytes[8];
    for (int i = 0; i < 8; ++i)
        bytes[i] = static_cast<unsigned char>((value >> (8 * i)) & 0xff);
    file.write(reinterpret_cast<const char *>(bytes), sizeof(bytes));
}

void writeArray(std::ofstream &file, const varType *data, std::uint64_t byteCount)
{
    writeUInt64(file, byteCount);
    file.write(reinterpret_cast<const char *>(data), static_cast<std::streamsize>(byteCount));
}
}

__host__ void writeVTI(size_t t, const std::string &outDir, varType *mom_host)
{
    namespace fs = std::filesystem;

    fs::path vtkDir = fs::path(outDir) / "vtk";
    fs::create_directories(vtkDir);

    std::ostringstream filename;
    filename << "output" << std::setw(6) << std::setfill('0') << t << ".vti";

    fs::path filepath = vtkDir / filename.str();
    std::ofstream file(filepath.string(), std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "Could not open VTI file for writing: " << filepath.string() << "\n";
        return;
    }

    const size_t pointCount = D2Q9::momSize;
    const std::uint64_t scalarBytes = pointCount * sizeof(varType);
    const std::uint64_t velocityBytes = 3 * scalarBytes;
    const std::uint64_t blockHeaderBytes = sizeof(std::uint64_t);
    const std::uint64_t uxOffset = blockHeaderBytes + scalarBytes;
    const std::uint64_t uyOffset = 2 * uxOffset;
    const std::uint64_t velocityOffset = 3 * uxOffset;
    const char *vtkType = sizeof(varType) == sizeof(double) ? "Float64" : "Float32";

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    file << "  <ImageData WholeExtent=\"0 " << (Geometry::NX - 1)
         << " 0 " << (Geometry::NY - 1)
         << " 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    file << "    <Piece Extent=\"0 " << (Geometry::NX - 1)
         << " 0 " << (Geometry::NY - 1)
         << " 0 0\">\n";
    file << "      <PointData>\n";
    file << "        <DataArray type=\"" << vtkType << "\" Name=\"rho\" format=\"appended\" offset=\"0\"/>\n";
    file << "        <DataArray type=\"" << vtkType << "\" Name=\"ux\" format=\"appended\" offset=\"" << uxOffset << "\"/>\n";
    file << "        <DataArray type=\"" << vtkType << "\" Name=\"uy\" format=\"appended\" offset=\"" << uyOffset << "\"/>\n";
    file << "        <DataArray type=\"" << vtkType << "\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << velocityOffset << "\"/>\n";
    file << "      </PointData>\n";
    file << "      <CellData/>\n";
    file << "    </Piece>\n";
    file << "  </ImageData>\n";
    file << "  <AppendedData encoding=\"raw\">_";

    writeArray(file, mom_host, scalarBytes);
    writeArray(file, mom_host + pointCount, scalarBytes);
    writeArray(file, mom_host + 2 * pointCount, scalarBytes);
    writeArray(file, mom_host + 3 * pointCount, velocityBytes);

    file << "</AppendedData>\n";
    file << "</VTKFile>\n";
}
