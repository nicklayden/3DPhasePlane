#include "cmdline.hpp"
#include <iterator>
#include <string>
#include <map>
#include <iostream>


void testmap(std::map<int, std::string>& m);

int main(int argc, char** argv){
    double r,t,a,b;
    int n;
    // cmdline_settings(argc, argv, &r, &t, &a, &b, &n);

    std::map<int, std::string> maptest;
    maptest.insert(std::make_pair(0,"fuck"));
    std::cout << "MAPPED: " << maptest[0] << std::endl;
    testmap(maptest);
    std::cout << "MAPPED: " << maptest[0] << std::endl;
    // std::cout << r << " " << t << " " << a << " " << b << " " << " " << n << std::endl;
    boost::program_options::variables_map vm;
    // vm.insert({"radius", 0});
    config_mapping(argc, argv, vm);
    std::cout << vm.size() << std::endl;
    std::cout << vm["init.conds.radius"].as<double>() << std::endl;
    std::cout << vm["init.conds.ycenter"].as<double>() << std::endl;


for (po::variables_map::iterator it = vm.begin(); it != vm.end(); it++) {
        std::cout << "> " << it->first;
        if (((boost::any)it->second.value()).empty()) {
            std::cout << "(empty)";
        }
        if (vm[it->first].defaulted() || it->second.defaulted()) {
            std::cout << "(default)";
        }
        std::cout << "=";

        bool is_char;
        try {
            boost::any_cast<const char *>(it->second.value());
            is_char = true;
        } catch (const boost::bad_any_cast &) {
            is_char = false;
        }
        bool is_str;
        try {
            boost::any_cast<std::string>(it->second.value());
            is_str = true;
        } catch (const boost::bad_any_cast &) {
            is_str = false;
        }

        if (((boost::any)it->second.value()).type() == typeid(int)) {
            std::cout << vm[it->first].as<int>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(bool)) {
            std::cout << vm[it->first].as<bool>() << std::endl;
        } else if (((boost::any)it->second.value()).type() == typeid(double)) {
            std::cout << vm[it->first].as<double>() << std::endl;
        } else if (is_char) {
            std::cout << vm[it->first].as<const char * >() << std::endl;
        } else if (is_str) {
            std::string temp = vm[it->first].as<std::string>();
            if (temp.size()) {
                std::cout << temp << std::endl;
            } else {
                std::cout << "true" << std::endl;
            }
        } else { // Assumes that the only remainder is vector<string>
            try {
                std::vector<std::string> vect = vm[it->first].as<std::vector<std::string> >();
                uint i = 0;
                for (std::vector<std::string>::iterator oit=vect.begin();
                     oit != vect.end(); oit++, ++i) {
                    std::cout << "\r> " << it->first << "[" << i << "]=" << (*oit) << std::endl;
                }
            } catch (const boost::bad_any_cast &) {
                std::cout << "UnknownType(" << ((boost::any)it->second.value()).type().name() << ")" << std::endl;
            }
        }
    }
}



void testmap(std::map<int, std::string>& testmap) {
    testmap[0] = "ffs.";
}