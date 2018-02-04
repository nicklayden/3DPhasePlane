#pragma once
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>

namespace {
    const int CMD_LINE_ERROR = 1;
    const int SUCCESS = 0;
    const int ERROR_EXCEPTION = 2;
}


namespace po = boost::program_options;

int cmdline_settings(int argc, char** argv, double* radius, double* thetamax,
                    double* a, double* b, int* N)
{

    std::ifstream config_file("config.ini");
    po::options_description init("Options");

    init.add_options()
        ("init.conds.radius", po::value<double>(radius), "Radius of circle to start initial condition")
        ("init.conds.thetamax",po::value<double>(thetamax), "Maximum angle through which to start initial conditions")
        ("init.conds.xcenter",po::value<double>(a), "Circle center in the x direction")
        ("init.conds.ycenter",po::value<double>(b), "Circle center in the y direction")
        ("init.conds.numinit",po::value<int>(N), "Number of initial conditions to create");


    po::variables_map vm;
    try 
    {
        po::store(po::parse_config_file(config_file,init), vm);

        po::notify(vm);


    }
    catch(po::error& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << init << std::endl;
        return CMD_LINE_ERROR;
    }

}