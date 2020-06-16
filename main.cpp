#include <fstream>
#include <string>

#include <boost/program_options.hpp>
#include <SurfaceMesh/SurfaceMesh.h>
#include <SurfaceMesh/IO.h>
#include <dbg.h>

#include "utils/FileIOUtils.h"
#include "utils/SurfaceMeshUtils.h"
#include "utils/EigenUtils.h"
#include "solver/LGSolver.h"
#include "arap.h"

namespace po = boost::program_options;

namespace LIM_para {
    struct argument {
        std::string base_mesh;
        float alpha{};
        float beta{};
        float mu_max{};
        float sigma_max{};
        float r{};
        float t{};
        float s_j{};
    };
}

int main(int argc, char *argv[]) {
    std::ifstream ifs(argv[1]);
    if (!ifs) return __LINE__;
    po::options_description desc("Available options");
    desc.add_options()
            ("help,h", "produce help message")
            ("base_mesh,b", "set the base mesh")
            ("alpha", "set alpha")
            ("beta", "set beta")
            ("mu_max", "set mu_max")
            ("sigma_max", "set sigma_max")
            ("r", "set r")
            ("t", "set t")
            ("s_j", "set s_j");
    po::variables_map vm;
    po::store(po::parse_config_file(ifs, desc), vm);
    ifs.close();

    po::notify(vm);
    if (vm.count("help")) {
        dbg(desc);
        return __LINE__;
    }

    LIM_para::argument args;
    {
        args.base_mesh = vm["base_mesh"].as<std::string>();
        args.alpha = std::stof(vm["alpha"].as<std::string>());
        args.beta = std::stof(vm["beta"].as<std::string>());
        args.mu_max = std::stof(vm["mu_max"].as<std::string>());
        args.sigma_max = std::stof(vm["sigma_max"].as<std::string>());
        args.s_j = std::stof(vm["s_j"].as<std::string>());
        args.r = std::stof(vm["r"].as<std::string>());
        args.t = std::stof(vm["t"].as<std::string>());
    }
    dbg(args.base_mesh);
    Surface_Mesh::SurfaceMesh mesh;
    Surface_Mesh::read_obj(mesh, args.base_mesh);
    dbg(args.base_mesh, mesh.n_vertices(), mesh.n_faces(), mesh.n_edges());
    std::vector<std::pair<size_t , Eigen::Vector2d>> pos_constrain = xry_mesh::readPosFile2D<double>(argv[2], 2);
    std::vector<std::pair<int , Eigen::VectorXf>> pos_constrain1 = xry_mesh::readPosFile<float>(argv[2], 2);
    Eigen::Matrix2Xd Vout;
    Eigen::Matrix3Xd V3d = xry_mesh::getGeometryMatrix<double>(mesh);
    auto V_ = V3d.block(0, 0, 2, V3d.cols());
    auto F_ = xry_mesh::getTopoMatrix(mesh);
    //jy_mesh::arap_deformation(V_, F_, pos_constrain, Vout, 100);
    LGSolver solver(mesh, pos_constrain1, 100);
    solver.init();
    auto x = solver.solve();
    auto v2d = Eigen::Map<Eigen::Matrix2Xf>(x.data(), 2, x.size() / 2);
    xry_mesh::rebuild2DMesh<float>(mesh, v2d);
    Surface_Mesh::write_obj(mesh, "result2D.obj");
    return 0;
}
