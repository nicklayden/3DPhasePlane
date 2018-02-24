/*
    CPLANE Version 2.0.
    Author: Nick Layden, 2018
        Dependencies:
            OpenGL: for drawing and displaying. (freeglut, glut, glu, glew, gl)
            SFML  : for handling the main window instance (SFML)
            boost : for Integrating N-dimensional ODEs


    There are a considerable number of optimizations to be done, so just compile with -O3
    
*/
#include <GL/glut.h>
#include <SFML/Graphics.hpp>
#include <vector>
#include <math.h>
#include <random>
#include <iostream>
#include <map>

#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint.hpp>

#include "cmdline.hpp"

// compile: g++ minimal-cplane.cpp -o n -lGL -lGLU -lGLEW -lglut -lsfml-system -lsfml-window -lsfml-graphics -O3
// run    : ./n
 

#define PI 3.14159265

inline double x_start(double y) {
    return sqrt(1 - y*y);
}

inline double decel(double x, double y) {
    return 2*x*x - y*y;
}

inline double random(double a, double b) {
    std::random_device rd;
    std::mt19937 eng(rd());
    // how to seed:
    // std::mt19937 eng;
    // eng.seed(1993);
    std::uniform_real_distribution<double> distr(a,b);

    return distr(eng);
}

double generate_random_circle(std::vector<double>& xin, std::vector<double>& yin, double a, double b, double radius, double thetamax, int N) {
    double r, theta, xt, yt;
    for (int i=0; i < N; i++) {
        r = random(0.,radius*radius);
        theta = random(0.,thetamax);
        xt = sqrt(r)*cos(theta) + a;
        yt = sqrt(r)*sin(theta) + b;
        xin.push_back(xt);
        yin.push_back(yt);
    }

}

double generate_random_annulus(std::vector<double>& xin, std::vector<double>& yin, double a, double b, double radius, double thetamax, int N) {
    double r, theta, xt, yt;
    for (int i=0; i < N; i++) {
        r = random(0.,radius*radius);
        theta = random(0.,thetamax);
        xt = sqrt(r)*cos(theta) + a;
        yt = sqrt(r)*sin(theta) + b;
        xin.push_back(xt);
        yin.push_back(yt);
    }

}

double generate_random_sphere(std::vector<double>& xin, std::vector<double>& yin, double a, double b, double radius, double thetamax, int N) {
    return 0; // on windows computer implementation
}

struct push_back_state_and_time 
{
    std::vector<std::vector<double> >& m_states;
    std::vector<double>& m_times;

    push_back_state_and_time(std::vector<std::vector<double> >& states, std::vector<double>& times)
    :m_states(states), m_times(times) { }

    void operator()(const std::vector<double>& x, double t)
    {
        m_states.push_back(x);
        m_times.push_back(t);
    }


};

/****************************************************************************
 * 
 *  CHANGE ONLY THE VARIABLES IN THIS DYNAMICAL SYSTEM 
 * 
 * **************************************************************************/

class sean_problem
{
public:
    double m_par;
    sean_problem(double par): m_par(par) {}
    // exponential potential with constant epsilon background

    void operator() (const std::vector<double>& xin, std::vector<double>& dxdt, const double /* t */) {
        double x,y,z;
        double c = sqrt(6.)/2.;
        x = xin[0]; y = xin[1]; z = xin[2];
    
        // dxdt[0] = x*(-2 + 2*x*x - y*y - z*z) + c*m_par*z*z;
        // dxdt[1] = y*(1 + 2*x*x - y*y - z*z);
        // dxdt[2] = z*(-c*m_par*x + 2*x*x + 1 - y*y - z*z);

        // dxdt[0] =  -2*x -  (x*x + y*y);
        // dxdt[1] = -(1.-x)*y;
        // dxdt[2] = 0;

        // dxdt[0] = -3*x + 3*x*x - m_par*y*(1/(sqrt(x*x + y*y)));
        // dxdt[1] = m_par*x*(1/(sqrt(x*x + y*y)));
        // dxdt[2] = 0;

        // Harmonic potential equations.
        // dxdt[0] = -3*x - z*y - x*(y*y - 2*x*x - 1);
        // dxdt[1] = z*x - y*(y*y - 2*x*x - 1);
        // dxdt[2] = z*(2*x*x - y*x + 1);
        double chi; double m = m_par;
        chi = -(1./3.)*z*z - (2./3.)*x*x + (1./3.)*y*y;
        // pick m = 25 to do numerics.
        // Alans suggested equations for harmonic potential
        dxdt[0] = -z*x - sqrt(1 - z*z)*m*y - x*z*chi;
        dxdt[1] = m*x*sqrt(1 - z*z) - z*y*chi;
        dxdt[2] = (1 - z*z)*chi;



        //lorenz
        // dxdt[0] = 10.*(y - x);
        // dxdt[1] = x*(28. - z) - y;
        // dxdt[2] = x*y - (8./3.)*z;
    }

    public:
        // Define Equilibrium points! - do later...

};




int main(int argc, char** argv ){

    /* #####################################################################################
        LOCAL VARIABLE DECLARATIONS AND PROGRAM SETTINGS
       #####################################################################################
    */
    // Create SFML window instance
    sf::RenderWindow mainwin(sf::VideoMode(800,800), "NICE TITLE");
    sf::Clock Clock;

    // x       : state vector (holds x,y,z coordinates of system)
    // t       : stores time values of the solution. (for autonomous systems, isn't useful)
    // y       : matrix containing solutions for all N coordinates
    // stepper : numerical method used to integrate the ODE
    std::vector<double> x(3);
    std::vector<double> t, xc, yc, inflx, infly;
    std::vector<std::vector<double> > y;
    boost::numeric::odeint::runge_kutta_dopri5<std::vector<double> > stepper;
    // boost::numeric::odeint::adams_bashforth<4,std::vector<double> > stepper;
    // boost::numeric::odeint::euler<std::vector<double> > stepper;


/****************************************************************************
 * 
 *  DEFINE YOUR SPECIFIC INITAL CONDITIONS HERE, YOU ADD THEM TO THE VECTORS
 *  AS I DID BELOW
 * 
 * **************************************************************************/

    std::vector<double> xt,yt,zt;
    // xt.push_back(1/sqrt(2));
    // yt.push_back(1/sqrt(2));
    // zt.push_back(0.);

    // xt.push_back(-1/sqrt(2));
    // yt.push_back(1/sqrt(2));
    // zt.push_back(0.);
    // xt.push_back(1); // this looks like an eq point for the harmonic potential
    // yt.push_back(1);
    // zt.push_back(sqrt(2));


    xt.push_back(-sqrt(2)/4);
    yt.push_back(-sqrt(2)/4);
    zt.push_back(1./2.);

    xt.push_back(sqrt(2)/4);
    yt.push_back(-sqrt(2)/4);
    zt.push_back(1./2.);
    
    xt.push_back(0);
    yt.push_back(0);
    zt.push_back(1e-5);


    xt.push_back(sqrt(2)/8);
    yt.push_back(-sqrt(2)/8);
    zt.push_back(1./2.);


    // xt.push_back(0);
    // yt.push_back(1);
    // zt.push_back(1e-5);

    // xt.push_back(1);
    // yt.push_back(0);
    // zt.push_back(1e-5);
    

    // xt.push_back(-1);
    // yt.push_back(0);
    // zt.push_back(1e-5);
    


    // Camera rotation flags and properties.
    bool rotate = false;
    bool rotatex = false; 
    bool rotatey = false;
    bool rotatez = false;
    float angle = 0;
    float anglex = 0;
    float angley = 0;
    float anglez = 0;
    float rotspeed = 10;
    // mainwin.setFramerateLimit(30);
    // Scalar field parameter.
    double lambda;



    // Initial conditions, dynamical system settings, and gui settings for the window.
    int numinit;
    double r, tmax, xcenter, ycenter, sysstep;
    float bkg_alpha, bkg_r, bkg_g, bkg_b, \
          gui_aspect, gui_zfar, gui_znear, \
          gui_fovy, gui_camera_dist;
    std::string sysmethod, sysfile;


    /* ######################################################################################
        LOAD PROGRAM CONFIGURATION FROM CONFIG.INI!!!!
    
        need to capture variables from the config file, and assign them in this block of main.
        #####################################################################################
    */
    boost::program_options::variables_map vm;
    config_mapping(argc, argv, vm);

    r = vm["init.conds.radius"].as<double>();
    tmax = vm["init.conds.thetamax"].as<double>();
    xcenter = vm["init.conds.xcenter"].as<double>();
    ycenter = vm["init.conds.ycenter"].as<double>();
    numinit = vm["init.conds.numinit"].as<int>();
    sysfile = vm["system.file"].as<std::string>();
    sysstep = vm["system.stepsize"].as<double>();
    sysmethod = vm["system.method"].as<std::string>();
    bkg_r = vm["gui.bkg_r"].as<float>();
    bkg_g = vm["gui.bkg_g"].as<float>();
    bkg_b = vm["gui.bkg_b"].as<float>();
    bkg_alpha = vm["gui.bkg_alpha"].as<float>();
    gui_fovy = vm["gui.fovy"].as<float>();
    gui_aspect = vm["gui.aspect"].as<float>();
    gui_znear = vm["gui.znear"].as<float>();
    gui_zfar = vm["gui.zfar"].as<float>();
    gui_camera_dist = vm["gui.camera_dist"].as<float>();


    // Generate uniform random numbers in a circle centerd at a,b with radius r:
    generate_random_circle(xc,yc,xcenter,ycenter,r,tmax,numinit);







    // preparing opengl surface. defining default window color.
    glClearDepth(1.f);
    glClearColor(bkg_r, bkg_g, bkg_b, bkg_alpha);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glEnable(GL_BLEND);
    glHint(GL_LINE_SMOOTH, GL_NICEST);

    // setup perspective projection and the cameras position:
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(gui_fovy, gui_aspect, gui_znear, gui_zfar); // fovy, aspect ratio, znear, zfar // objects farther than zfar or closer than znear will NOT be rendered.
    glPointSize(5); // size of point sprites on the window. default is 1 pixel.
    // glRotatef(90.,1.,0.,0.);


    // Drawing circles around the axis for spatial reference.
    std::vector<float> xcirc,ycirc;
    float th;
    for (int i = 0; i < 80; i++) {
        th = 2.*PI*i/80.;
        xcirc.push_back(cos(th));
        ycirc.push_back(sin(th));
    }

    // glutInit(&argc, argv);
    bool running = true; // This is used so that the loop can end and opengl frees resources properly!!
    // SFML main window instance. Drawing handled in pure opengl context.
    while (running) {
        // for (int q = 0; q < 250; q++) {

            // Handle events through the SFML interface for the window. (keyboard press, closing, etc.)
            sf::Event event;
            while (mainwin.pollEvent(event)) {
                if (event.type == sf::Event::Closed) {
                    // mainwin.close();
                    running = false;
                }
                if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Right)){
				    rotspeed += 1;
			    }
                if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Left)){
				    rotspeed -= 1;
			    }
                
                // X ROTATIONS
                if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Numpad7)){
				    rotatex=!rotatex; // toggle continuous rotation in x
			    }
                if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::Numpad1))){
                    anglex += 0.1; // hold down positive rotation in x
			    }
                if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::Numpad4))){
                    anglex -= 0.1; // hold down positive rotation in x
			    }

                // Y ROTATIONS
                if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Numpad8)){
				    rotatey=!rotatey; // toggle continuous rotation in y
			    }
                if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::Numpad2))){
                    angley += 0.1; // hold down positive rotation in y
			    }
                if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::Numpad5))){
                    angley -= 0.1; // hold down positive rotation in y
			    }
                
                // Z ROTATIONS
                if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Numpad9)){
				    rotatez=!rotatez; // toggle continuous rotation in z
			    }
                if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::Numpad3))){
                    anglez += 0.1; // hold down positive rotation in z
			    }
                if ((event.type == sf::Event::KeyPressed) && (sf::Keyboard::isKeyPressed(sf::Keyboard::Numpad6))){
                    anglez -= 0.1; // hold down negative rotation in z
			    } 
                if (event.type == sf::Event::Resized)
                {
                    // adjust the viewport when the window is resized, adjusting
                    // the aspect ratio so we don't stretch the image.
                    float new_aspect = (float)mainwin.getSize().x/(float)mainwin.getSize().y;
                    // Adjusts the size of the viewport to reflect window resizing.
                    glViewport(0, 0, event.size.width, event.size.height);
                    // this is old opengl, switch to projection matrix mode,
                    // load identity so we dont stack projections by mistake,
                    // then change our perspective to the new aspect ratio
                    // then finalize by switching back to modelview matrix mode
                    // so that the drawings arent fucked up after.
                    glMatrixMode(GL_PROJECTION);
                    glLoadIdentity();
                    gluPerspective(gui_fovy, new_aspect, gui_znear, gui_zfar);
                    glMatrixMode(GL_MODELVIEW);
                }
                if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Up)){
				    gui_camera_dist -= 0.2; // toggle continuous rotation in z
			    }
                if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Down)){
				    gui_camera_dist += 0.2; // toggle continuous rotation in z
			    }


            } // event loop

            // clear window to prepare for drawing. fill window background with default color (defined above)
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            glTranslatef(0.f,0.f,-gui_camera_dist);
            glRotatef(-90.,1.,0.,0.); // this rotates to put XY plane into perspective

            // glColor4f(1,1,1,0.);
            // glutSolidSphere(1,50,50);
            // perspective transformations... rotate the coordinates if toggles are on.
            if (rotatex) {
                anglex += .05;
            }
            if (rotatey) {
                angley += .05;
            }
            if (rotatez) {
                anglez += .05;
            }
            // opengl rotation functions.
            glRotatef(anglex*rotspeed,1.0f,0.f,0.f);
            glRotatef(angley*rotspeed,0.0f,1.0f,0.0f);
            glRotatef(anglez*rotspeed,0.0f,0.0f,1.0f);

            // numerical solutions...
            double ystart = 0.01; // value for initial conditions. Will change to random seeded values later?
            double xp, yp, zp; // temporary variables for holding numbers. These are (probably) removed with -O3 optimization
            lambda = 25; // scalar field parameter. can be fixed or change through the q loop. if fixed, q loop will be optimized away.
    
            // Numerically solve and draw lines...
            for (int j = 0; j < xc.size(); j++) {
                
                // System to solve!
                sean_problem test(lambda); 



                //FORWARD TIME INTEGRATION.
                // INITIAL CONDITIONS SET FROM VECTORS DEFINED ABOVE
                // x[0] = xc[j]; x[1] = yc[j]; x[2] = 0.5; // This line defines randomly generated initial conditions 
                x[0] = xt[j]; x[1] = yt[j]; x[2] = zt[j];
                boost::numeric::odeint::integrate_const(stepper,test, x, 0.,100.,0.01,push_back_state_and_time(y,t));

                // DRAWS FORWARD SOLUTION
                glBegin(GL_LINE_STRIP);
                    for (int i = 0; i < y.size(); i++) {
                        xp = y[i][0]; yp = y[i][1]; zp = y[i][2]; //placeholders for readability
                        glColor3f(0.4,0,0.8);
                        glVertex3f(xp, yp, zp);
                    }
                glEnd();
                y.clear(); t.clear();


                // DRAWS BACKWARD SOLUTION
                boost::numeric::odeint::integrate_const(stepper,test, x, 0.,-100.,-0.01,push_back_state_and_time(y,t));
                
                // Draw solution curves in 3D phase space.
                glBegin(GL_LINE_STRIP);
                    for (int i = 0; i < y.size(); i++) {
                        xp = y[i][0]; yp = y[i][1]; zp = y[i][2]; //placeholders for readability
                        glColor3f(0.4,0,0.8);
                        glVertex3f(xp, yp, zp);
                    }
                glEnd();

                glBegin(GL_POINTS);
                glColor3f(1,0,0);
                glVertex3f(xt[j],yt[j],zt[j]);
                glEnd();

                // delete current solution curve to prepare for the next one...
                y.clear(); t.clear();

            }

            // These three line loops are drawing the circles around the 3 axes just for perspective
            glBegin(GL_LINE_LOOP);
                glColor3f(1,0,0);
                for (int i = 0; i < xcirc.size(); i++) {
                    glVertex3f(xcirc[i], ycirc[i], 0);
                }
            glEnd();
            glBegin(GL_LINE_LOOP);
                glColor3f(0,0.5,1);
                for (int i = 0; i < xcirc.size(); i++) {
                    glVertex3f(xcirc[i], 0, ycirc[i]);
                }
            glEnd();
            glBegin(GL_LINE_LOOP);
                glColor3f(1,0,1);
                for (int i = 0; i < xcirc.size(); i++) {
                    glVertex3f(0, xcirc[i], ycirc[i]);
                }
            glEnd();


            // DRAW EQUILIBRIUM POINTS
            // these floats are for the non trivial eq points.
            // float fepx = lambda/sqrt(6.);
            // float fepz = sqrt(1 - lambda*lambda/6.);
            // float hepx = sqrt(6.)/(3.*lambda);
            // float hepz = 2./(sqrt(3.)*lambda); 
            glBegin(GL_POINTS);
                // glColor3f(1,0,1);
                // glVertex3f(0,0,0);
                // glVertex3f(1,0,0);
                // glVertex3f(-1,0,0);
                // glVertex3f(0,1,0);
                // glVertex3f(0,-1,0);
                // glVertex3f(fepx,0,fepz);
                // glVertex3f(hepx,0,hepz);
            glEnd();

            // Display everything we've done to the SFML instance
            mainwin.display();
        // } // end q loop
    } // end display loop

} // end main
