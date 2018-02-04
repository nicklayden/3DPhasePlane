# 3DPhasePlane
Solve systems of ODEs and plot phase planes in up to 3 dimensions

## Dependencies
OpenGL
SFML 
Boost libraries

## Basic preparation for install
Check that you have SFML:
``` 
    sudo apt-get install libsfml-dev 

```
Check that you have boost:
```
    sudo apt-get install libboost-all-dev
```
or
```
    sudo apt-get install libboost-all
```
Check that you have OpenGL libraries:
```
    sudo apt-get install libgl1-mesa-dev libglu1-mesa-dev freeglut3-dev
```
Compile with:
```
    g++ minimal-cplane.cpp -lGL -lGLU -lGLEW -lglut -lsfml-graphics -lsfml-system -lsfml-window -lboost_program_options -O3 -o plane
```
