/*
 * This file was written by Ben Morton (bmorton@ed.ac.uk).
 */

#include <vector>

#include "constants3D.h"
#include "vertex3D.h"

#ifdef ANALYTIC_GRAVITY
void plummer_gravity(std::vector<VERTEX> &MY_POINTS, double DT, int N_POINTS){
        // Fixed Plummer potential at (XC,YC)
        int i;
        double X, Y, Z, GM, AX, AY, AZ, DU[5];
        double XC = 5.0, YC = 5.0, ZC = 5.0, MPERT = 3.28E+05, EPS = 0.145;
        double MASS_DENSITY, DELTAX, DELTAY, DELTAZ, RAD2;
        VERTEX MY_VERTEX;

        DU[0] = DU[4] = 0.0;

        for(i=0;i<N_POINTS;++i){
                MY_VERTEX = MY_POINTS[i];
                X = MY_VERTEX.get_x();
                Y = MY_VERTEX.get_y();
                Z = MY_VERTEX.get_z();
                MASS_DENSITY = MY_VERTEX.get_mass_density();
                DELTAX = X - XC, DELTAY = Y - YC, DELTAZ = Z - ZC;
                RAD2 = DELTAX*DELTAX + DELTAY*DELTAY + DELTAZ*DELTAZ;
                GM = GRAV * MPERT / (sqrt(RAD2 + EPS*EPS)*sqrt(RAD2 + EPS*EPS)*sqrt(RAD2 + EPS*EPS));
                AX = DELTAX * GM;
                AY = DELTAY * GM;
                AZ = DELTAZ * GM;
                DU[1] = -1.0*AX * DT * MASS_DENSITY;
                DU[2] = -1.0*AY * DT * MASS_DENSITY;
                DU[3] = -1.0*AZ * DT * MASS_DENSITY;
                MY_POINTS[i].update_du(DU);
        }
        return ;
}
#endif

#ifdef SELF_GRAVITY
void direct_gravity(std::vector<VERTEX> MY_POINTS, double DT, int N_POINTS){
        double X0, X1, Y0, Y1, Z0, Z1, X_DIFF, Y_DIFF, Z_DIFF, AX, AY, AZ, MASS, MASS_DENSITY, DU[5];
        DU[0] = DU[4] = 0.0;
        for(int i=0; i<N_POINTS; ++i){
                X0 = MY_POINTS[i].get_x();
                Y0 = MY_POINTS[i].get_y();
                Z0 = MY_POINTS[i].get_z();
                MASS = MY_POINTS[i].get_mass();
                for(int j=0; j<N_POINTS; ++j){
                        if(j!=i){
                                X1 = MY_POINTS[j].get_x();
                                Y1 = MY_POINTS[j].get_y();
                                Z1 = MY_POINTS[j].get_z();
                                X_DIFF = X1 - X0;
                                Y_DIFF = Y1 - Y0;
                                Z_DIFF = Z1 - Z0;
                                AX = AX - GRAV * MASS * X_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF + Z_DIFF*Z_DIFF)),3);
                                AY = AY - GRAV * MASS * Y_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF + Z_DIFF*Z_DIFF)),3);
                                AZ = AZ - GRAV * MASS * Z_DIFF/pow((sqrt(X_DIFF*X_DIFF + Y_DIFF*Y_DIFF + Z_DIFF*Z_DIFF)),3);
                                DU[1] = -1.0*AX*DT*MASS_DENSITY;
                                DU[2] = -1.0*AY*DT*MASS_DENSITY;
                                DU[3] = -1.0*AZ*DT*MASS_DENSITY;
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
}