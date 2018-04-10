#pragma once
/**
 * 
 *  Class for openGL window handling for 3DPhasePlane program.
 * 
 *  
 * 
 * 
 * */
#include <iostream>
#include <SFML/Graphics.hpp>


template<typename Real = float, int DIM = 3 >
class Display
{
public:
    Display();

private:
    Real xlim;
    Real ylim;
    int xsize;
    int ysize;




    void creation();
    void set_default_context();
    template<typename... Args>  
    void print(Args&&... args);

    sf::ContextSettings settings;

};

template<typename Real, int DIM >
Display<Real, DIM>::Display()
{   
    // Instantiation message.
    creation();

    // Set default context settings
    set_default_context();

    // Initialize openGL context

}

template<typename Real, int DIM>
void Display<Real, DIM>::creation()
{
    // Display message on object creation, called in constructor.
    std::cout << "Created Display with  " << DIM << " dimensions.\n";
}

template<typename Real, int DIM>
void Display<Real, DIM>::set_default_context()
{
    settings.depthBits=24; settings.majorVersion=3; settings.minorVersion=0;

}

template<typename Real, int DIM>
template<typename ... Args>
void Display<Real, DIM>::print(Args&&... args)
{
    (std::cout << ... << args) << "\n";
}


/***
 * 
 *  Pre-defined types for usage without template calling.
 * 
 * 
 * 
*/

typedef Display<float, 3> Display3f;
typedef Display<float, 2> Display2f;
typedef Display<double, 3> Display3d;
typedef Display<double, 2> Display2d;
