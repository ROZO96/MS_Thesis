/*
 * This file was written by Ben Morton (bmorton@ed.ac.uk) and Zhenyu Wu (zhenyu.wu@ed.ac.uk).
 */

#include <vector>

#include "constants.h"
#include "vertex2D.h"

#ifdef ANALYTIC_GRAVITY
void plummer_gravity(std::vector<VERTEX> &MY_POINTS, double DT, int N_POINTS){
        // Fixed Plummer potential at (XC,YC)
        int i;
        double X, Y, GM, AX, AY, DU[4];
        double XC = 5.0, YC = 5.0, MPERT = 3.28E+05, EPS = 0.145;
        double MASS_DENSITY, DELTAX, DELTAY, RAD2;
        VERTEX MY_VERTEX;

        DU[0] = DU[3] = 0.0;

        for(i=0;i<N_POINTS;++i){
                MY_VERTEX = MY_POINTS[i];
                X = MY_VERTEX.get_x();
                Y = MY_VERTEX.get_y();
                MASS_DENSITY = MY_VERTEX.get_mass_density();
                DELTAX = X - XC, DELTAY = Y - YC;
                RAD2 = DELTAX*DELTAX + DELTAY*DELTAY;
                GM = GRAV * MPERT / (sqrt(RAD2 + EPS*EPS)*sqrt(RAD2 + EPS*EPS)*sqrt(RAD2 + EPS*EPS));
                AX = DELTAX * GM;
                AY = DELTAY * GM;
                DU[1] = -1.0*AX * DT * MASS_DENSITY;
                DU[2] = -1.0*AY * DT * MASS_DENSITY;
                MY_POINTS[i].update_du(DU);
        }
        return ;
}
#endif

#ifdef SELF_GRAVITY
void direct_gravity(std::vector<VERTEX> &MY_POINTS, double DT, int N_POINTS){
        double X0, X1, Y0, Y1, X_DIFF, Y_DIFF, AX, AY, MASS, MASS_DENSITY, DU[4];
        DU[0] = DU[3] = 0.0;
        for(int i=0; i<N_POINTS; ++i){
                X0 = MY_POINTS[i].get_x();
                Y0 = MY_POINTS[i].get_y();
                MASS = MY_POINTS[i].get_mass();
                MASS_DENSITY = MY_POINTS[i].get_mass_density();
                for(int j=0; j<N_POINTS; ++j){
                        if(j!=i){
                                X1 = MY_POINTS[j].get_x();
                                Y1 = MY_POINTS[j].get_y();
                                X_DIFF = X1 - X0;
                                Y_DIFF = Y1 - Y0;
                                AX = AX - GRAV * MASS * X_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF)),3);
                                AY = AY - GRAV * MASS * Y_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF)),3);
                                DU[1] = -1.0*AX*DT*MASS_DENSITY;
                                DU[2] = -1.0*AY*DT*MASS_DENSITY;
                                MY_POINTS[i].update_du(DU);
                        }
                }
        }
}
#endif

void sources(std::vector<VERTEX> &MY_POINTS, double DT, int N_POINTS){
#ifdef ANALYTIC_GRAVITY
        plummer_gravity(MY_POINTS, DT, N_POINTS);
#endif
#ifdef SELF_GRAVITY
        direct_gravity(MY_POINTS, DT, N_POINTS);
#endif
        return ;
}
