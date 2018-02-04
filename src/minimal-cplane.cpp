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

#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/numeric/odeint.hpp>

#include "cmdline.hpp"

// compile: g++ minimal-cplane.cpp -o n -lGL -lGLU -lGLEW -lglut -lsfml-system -lsfml-window -lsfml-graphics -O3
// run    : ./n
 

#define PI 3.14159265

inline double x_start(double y) {
    return sqrt(1 - y*y);
}

double random(double a, double b) {
    std::random_device rd;
    std::mt19937 eng(rd());
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
    
        dxdt[0] = x*(-2 + 2*x*x - y*y - z*z) + c*m_par*z*z;
        dxdt[1] = y*(1 + 2*x*x - y*y - z*z);
        dxdt[2] = z*(-c*m_par*x + 2*x*x + 1 - y*y - z*z);

        //lorenz
        // dxdt[0] = 10.*(y - x);
        // dxdt[1] = x*(28. - z) - y;
        // dxdt[2] = x*y - (8./3.)*z;
    }
};




int main(int argc, char** argv ){
    // Create SFML window instance
    sf::RenderWindow mainwin(sf::VideoMode(800,800), "NICE TITLE");
    sf::Clock Clock;

    // x       : state vector (holds x,y,z coordinates of system)
    // t       : stores time values of the solution. (for autonomous systems, isn't useful)
    // y       : matrix containing solutions for all N coordinates
    // stepper : numerical method used to integrate the ODE
    std::vector<double> x(3);
    std::vector<double> t, xc, yc;
    std::vector<std::vector<double> > y;
    boost::numeric::odeint::runge_kutta_dopri5<std::vector<double> > stepper;
    // boost::numeric::odeint::adams_bashforth<4,std::vector<double> > stepper;
    // boost::numeric::odeint::euler<std::vector<double> > stepper;

    // Generate uniform random numbers in a circle centerd at a,b with radius r:
    double a,b,r,tmax;
    int n;
    // a = 1.5/sqrt(6);
    // b = sqrt(1 - (1.5*1.5/6));
    // a=0;b=0;
    // r = 1.;
    cmdline_settings(argc, argv, &r, &tmax, &a, &b, &n);
    generate_random_circle(xc,yc,a,b,r,tmax,n);


    // preparing opengl surface. defining default window color.
    // glClearDepth(1.f);
    glClearColor(0.f, 0.f, 0.f, 0.f);
    // glEnable(GL_DEPTH_TEST);
    // glDepthMask(GL_TRUE);

    // setup perspective projection and the cameras position:
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.f, 1., 1.f, 5.); // fovy, aspect ratio, znear, zfar // objects farther than zfar or closer than znear will NOT be rendered.
    glPointSize(5); // size of point sprites on the window. default is 1 pixel.

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

    // Scalar field parameter.
    double lambda;

    // Drawing circles around the axis for spatial reference.
    std::vector<float> xcirc,ycirc;
    float th;
    for (int i = 0; i < 80; i++) {
        th = 2.*PI*i/80.;
        xcirc.push_back(cos(th));
        ycirc.push_back(sin(th));
    }


    // SFML main window instance. Drawing handled in pure opengl context.
    while (mainwin.isOpen()) {
        for (int q = 0; q < 1; q++) {

            // Handle events through the SFML interface for the window. (keyboard press, closing, etc.)
            sf::Event event;
            while (mainwin.pollEvent(event)) {
                if (event.type == sf::Event::Closed) {
                    mainwin.close();
                }
                if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Up)){
				    rotspeed += 1;
			    }
                if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Down)){
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
            } // event loop

            // clear window to prepare for drawing. fill window background with default color (defined above)
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            glTranslatef(0.f,0.f,-2.f);
            
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
            lambda = 10; // scalar field parameter. can be fixed or change through the q loop. if fixed, q loop will be optimized away.
            
            // Numerically solve and draw lines...
            for (int j = 0; j < xc.size(); j++) {
                
                // System to solve!
                sean_problem test(lambda); 

                //forward solution
                // x[0] = 0.01*j; x[1] = 0.01*j; x[2] = 0.01*j;// generic initial conditions
                // x[0] = -x_start(ystart) + 2*x_start(ystart)*j/20; x[2] = ystart; x[1] = 0.01; // XZ plane solution initial conditions (same as 2D 'simple' solution)
                x[0] = xc[j]; x[1] = 0.01; x[2] = yc[j];
                boost::numeric::odeint::integrate_const(stepper,test, x, 0.,10.,0.01,push_back_state_and_time(y,t)); // forward solution
                glBegin(GL_LINE_STRIP);
                    for (int i = 0; i < y.size(); i++) {
                        glColor3f(0,1,0);
                        xp = y[i][0]; yp = y[i][1]; zp = y[i][2];
                        glVertex3f(xp, yp, zp);
                    }
                glEnd();
                y.clear(); t.clear();
                
                // backward solution
                // x[0] = 0.01*j; x[1] = 0.01*j; x[2] = 0.01*j;
                // x[0] = -x_start(ystart) + 2*x_start(ystart)*j/20; x[2] = ystart; x[1] = 0.01; // XZ plane solution conditions.
                x[0] = xc[j]; x[1] = 0.01; x[2] = yc[j];
                boost::numeric::odeint::integrate_const(stepper,test, x, 10.,0.,-0.01,push_back_state_and_time(y,t)); // forward solution

                // Draw solution curves in 3D phase space.
                glBegin(GL_LINE_STRIP);
                    for (int i = 0; i < y.size(); i++) {
                        glColor3f(0,1,0);
                        xp = y[i][0]; yp = y[i][1]; zp = y[i][2];
                        glVertex3f(xp, yp, zp);
                    }
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
            float fepx = -lambda/sqrt(6.);
            float fepz = sqrt(1 - lambda*lambda/6.);
            float hepx = sqrt(6.)/(3.*lambda);
            float hepz = 2./(sqrt(3.)*lambda); 
            glBegin(GL_POINTS);
                glColor3f(1,1,1);
                glVertex3f(0,0,0);
                glVertex3f(1,0,0);
                glVertex3f(-1,0,0);
                glVertex3f(0,1,0);
                glVertex3f(0,-1,0);
                glVertex3f(fepx,0,fepz);
                glVertex3f(hepx,0,hepz);
                for (int i =0; i < xc.size(); i++) {
                    glVertex3f(xc[i],0,yc[i]);
                }
            glEnd();

            // Display everything we've done to the SFML instance
            mainwin.display();
        } // end q loop
    } // end display loop

} // end main