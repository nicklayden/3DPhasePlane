#include "cmdline.hpp"


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

int config_mapping(int argc, char** argv, boost::program_options::variables_map& dictionary)
{
    /*
        To map from config.ini into a dictionary of values (which is a map object)
        we have to cast them into the dictionary, then they lose type, then when
        extracting the values for actual use, we must cast back into their original
        type via:
            dictionary["init.conds.radius"].as<double>()
        for example.
    
    */

    std::ifstream config_file("config.ini");
    po::options_description init("Options");

    init.add_options()
        ("init.conds.radius", po::value<double>(), "Radius of circle to start initial condition")
        ("init.conds.dr", po::value<double>(), "Thickness of annulus if specified.")
        ("init.conds.thetamax",po::value<double>(), "Maximum angle through which to start initial conditions")
        ("init.conds.xcenter",po::value<double>(), "Circle center in the x direction")
        ("init.conds.ycenter",po::value<double>(), "Circle center in the y direction")
        ("init.conds.numinit",po::value<int>(), "Number of initial conditions to create")
        ("system.file", po::value<std::string>(), "Radius of circle to start initial condition")
        ("system.stepsize",po::value<double>(), "Maximum angle through which to start initial conditions")
        ("system.method",po::value<std::string>(), "Circle center in the x direction")
        ("system.timemin",po::value<double>(), "initial 'time' value for integration")
        ("system.timemax",po::value<double>(), "maximum time value")
        ("system.dt",po::value<double>(), "time step")
        ("gui.bkg_r", po::value<float>(), "background red colour (0,1)")
        ("gui.bkg_g",po::value<float>(), "background green colour (0,1)")
        ("gui.bkg_b",po::value<float>(), "background blue colour (0,1)")
        ("gui.bkg_alpha",po::value<float>(), "background alpha (0,1)")
        ("gui.fovy",po::value<float>(), "field of view in the vertical direction")
        ("gui.aspect",po::value<float>(), "aspect ratio")
        ("gui.znear",po::value<float>(), "Near cutting plane (nearest drawing distance)")
        ("gui.zfar",po::value<float>(), "Far cutting plane (farthest drawing distance)")
        ("gui.camera_dist",po::value<float>(), "distance of camera from origin");




    // po::variables_map vm;
    try 
    {
        po::store(po::parse_config_file(config_file,init), dictionary);

        po::notify(dictionary);
        std::cout << "Mapped " << dictionary.size() << " configurations." << std::endl;

    }
    catch(po::error& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << init << std::endl;
        return CMD_LINE_ERROR;
    }

}
