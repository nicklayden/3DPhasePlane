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

## Compiling and Running the Program
Compile with:
```
    make
```
Run the program:
```
    ./run
```
## Configuration File Options
This program comes with a config.ini file for altering some of the graphics and numerics parameters without having to touch the cpp file, or recompile for the changes. 

| Option      | Flag Type | Comment                                                                                                       |
|-------------|-----------|---------------------------------------------------------------------------------------------------------------|
| radius      | Real      | Sets the radius for initial conditions to be drawn in a circle                                                |
| thetamax    | Real      | Maximum angle theta through with we will populate a circle with initial conditions.                           |
| xcenter     | Real      | Center of the circles x coordinate.                                                                           |
| ycenter     | Real      | Center of the circles y coordinate.                                                                           |
| numinit     | Integer   | Number of randomly generated initial conditions to spawn inside above circle                                  |
| bkg_r       | Real      | Background R value for color, range is [0,1]                                                                  |
| bkg_g       | Real      | Background G value for color, range is [0,1]                                                                  |
| bkg_b       | Real      | Background B value for color, range is [0,1]                                                                  |
| bkg_alpha   | Real      | Background Alpha transparency value. range is [0,1]                                                           |
| fovy        | Real      | Field of view from bottom to top of window. Measured in degrees.                                              |
| aspect      | Real      | Aspect Ratio of window. 1. = equivalent height and widths. resizing the window will keep this value constant. |
| znear       | Real      | Near cutting plane of the graphics. Objects within this distance from the camera position will not be drawn.  |
| zfar        | Real      | Far cutting plane of the graphics. Objects past this distance will not be drawn.                              |
| camera_dist | Real      | Initial camera distance away from the origin in the z direction.                                              |
