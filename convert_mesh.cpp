#include "read_mesh.hpp"
#include <fmt/format.h>

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        fmt::print(stderr, "Usage: convert_mesh [mesh order] [mesh file]\n");
        return 1;
    }
    const int order = atoi(argv[1]);
    if (!(order >= 1 && order <= 4))
    {
        fmt::print(stderr, "Mesh order should be between 1 and 4\n");
        return 2;
    }

    fmt::print("Converting mesh from file {:s}...\n", argv[2]);
    std::string output_file_name = fmt::format("{:s}.dat", argv[2]);
    FILE *outfile = fopen(output_file_name.c_str(), "wb");

    if (!outfile)
    {
        fmt::print(stderr, "Failed to open {:s} for writing\n", output_file_name);
        return 3;
    }
    auto mesh = Elasticity::read_mesh(argv[2], order);
    std::visit([=](const auto &m) { m.serialize(outfile); }, mesh);
    fclose(outfile);

    return 0;
}
