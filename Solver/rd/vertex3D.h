/*
 * This file was written by Ben Morton (bmorton@ed.ac.uk).
 * class containing conserved and primative variables of fluid at the position of the vertex
        X = x position of vertex
        Y = y position of vertex
        DX = change in x between vertices (should no longer be used)
        DY = change in y between vertices
        DT_REQ = timestep required by vertex state
        DUAL = dual cell area
        LEN_VEL_SUM = sum of product of edge length and velocity for every edge (ised in dt calc)
        U_VARIABLES = vector of fluid variables
        DU = sum of change in fluid variables for first half timestep
        MASS_DENSITY = mass_density of material at vertex
        X_VELOCITY = x velocity of material at vertex
        Y_VELOCITY = y velocity of material at vertex
        PRESSURE = pressure at vertex
        SPECIFIC_ENERGY = specific energy density at vertex
        U_HALF = vector of fluid variables for intermediate state
        DU_HALF = sum of change in fluid variables for second half timestep
        MASS_DENSIT_HALF = mass_density of material for vertex at intermediate state
        X_VELOCITY_HALF = x velocity of material at vertex at intermediate state
        Y_VELOCITY_HALF = y velocity of material at vertex at intermediate state
        PRESSURE_HALF = pressure at vertex at intermediate state
        SPECIFIC_ENERGY_HALF = specific energy density at vertex at intermediate state
 */

#include <vector>

#include "constants3D.h"

#pragma once

class VERTEX{

private:

        int ID,TBIN_LOCAL;
        double X, Y, Z, DX, DY, DZ;
        double DT_REQ;
        double DUAL,LEN_VEL_SUM;
        double U_VARIABLES[5], DU[5];
        double MASS_DENSITY, X_VELOCITY, Y_VELOCITY, Z_VELOCITY;
        double PRESSURE, SPECIFIC_ENERGY;
        double U_HALF[5], DU_HALF[5];
        double MASS_DENSITY_HALF, X_VELOCITY_HALF, Y_VELOCITY_HALF, Z_VELOCITY_HALF;
        double PRESSURE_HALF, SPECIFIC_ENERGY_HALF;

        std::vector<int> ASSOC_TRIANG; // not used yet

public:

        // setter functions preventing varaibles being changed accidentally
        // (no setter functions for U and F(U) as these are set by the other variables)
        void set_id(int NEW_ID){ID = NEW_ID;};
        void set_tbin_local(int NEW_TBIN){TBIN_LOCAL = NEW_TBIN;}
        void set_x(   double NEW_X){X   = NEW_X;}
        void set_y(   double NEW_Y){Y   = NEW_Y;}
        void set_z(   double NEW_Z){Z   = NEW_Z;}
        void set_dx(  double NEW_DX){DX = NEW_DX;}
        void set_dy(  double NEW_DY){DY = NEW_DY;}
        void set_dz(  double NEW_DZ){DZ = NEW_DZ;}
        void set_dual(double NEW_DUAL){DUAL = NEW_DUAL;}
        void set_mass_density( double NEW_MASS_DENSITY){MASS_DENSITY  = NEW_MASS_DENSITY;}
        void set_x_velocity(   double NEW_X_VELOCITY){  X_VELOCITY    = NEW_X_VELOCITY;}
        void set_y_velocity(   double NEW_Y_VELOCITY){  Y_VELOCITY    = NEW_Y_VELOCITY;}
        void set_z_velocity(   double NEW_Z_VELOCITY){  Z_VELOCITY    = NEW_Z_VELOCITY;}
        void set_pressure(     double NEW_PRESSURE){    PRESSURE      = NEW_PRESSURE;}
        void set_pressure_half(double NEW_PRESSURE){    PRESSURE_HALF = NEW_PRESSURE;}

        void set_u0(double NEW){U_VARIABLES[0] = NEW;}
        void set_u1(double NEW){U_VARIABLES[1] = NEW;}
        void set_u2(double NEW){U_VARIABLES[2] = NEW;}
        void set_u3(double NEW){U_VARIABLES[3] = NEW;}
        void set_u4(double NEW){U_VARIABLES[4] = NEW;}

        // add triangle to list of those associated with this vertex

        void add_triang(int NEW_TRIANGLE){ASSOC_TRIANG.push_back(NEW_TRIANGLE);} // not used yet

        // getter functions for eXtracting values of variables
        int get_id(){return ID;}
        int get_tbin_local(){return TBIN_LOCAL;}
        double get_x(){      return X;}
        double get_y(){      return Y;}
        double get_z(){      return Z;}
        double get_dx(){     return DX;}
        double get_dy(){     return DY;}
        double get_dz(){     return DZ;}
        double get_dt_req(){ return DT_REQ;}
        double get_dual(){   return DUAL;}
        double get_mass(){   return DUAL*MASS_DENSITY;}

        double get_specific_energy(){return SPECIFIC_ENERGY;}
        double get_mass_density(){   return MASS_DENSITY;}
        double get_x_velocity(){     return X_VELOCITY;}
        double get_y_velocity(){     return Y_VELOCITY;}
        double get_z_velocity(){     return Z_VELOCITY;}
        double get_pressure(){       return PRESSURE;}
        double get_u0(){return U_VARIABLES[0];}
        double get_u1(){return U_VARIABLES[1];}
        double get_u2(){return U_VARIABLES[2];}
        double get_u3(){return U_VARIABLES[3];}
        double get_u4(){return U_VARIABLES[4];}

        double get_specific_energy_half(){return SPECIFIC_ENERGY_HALF;}
        double get_mass_density_half(){   return MASS_DENSITY_HALF;}
        double get_x_velocity_half(){     return X_VELOCITY_HALF;}
        double get_y_velocity_half(){     return Y_VELOCITY_HALF;}
        double get_z_velocity_half(){     return Z_VELOCITY_HALF;}
        double get_pressure_half(){       return PRESSURE_HALF;}
        double get_u0_half(){return U_HALF[0];}
        double get_u1_half(){return U_HALF[1];}
        double get_u2_half(){return U_HALF[2];}
        double get_u3_half(){return U_HALF[3];}
        double get_u4_half(){return U_HALF[4];}


        // set up the specific energy varaible, as well as u and f arrays
        void setup_specific_energy(){
                double VEL_SQ_SUM = X_VELOCITY*X_VELOCITY + Y_VELOCITY*Y_VELOCITY + Z_VELOCITY*Z_VELOCITY;
                SPECIFIC_ENERGY = PRESSURE/((GAMMA-1.0)*MASS_DENSITY) + VEL_SQ_SUM/2.0; // calculate specific energy
        }

        // recacluate pressure based on current primitive varaibles
        void recalculate_pressure(){
                double VEL_SQ_SUM = X_VELOCITY*X_VELOCITY + Y_VELOCITY*Y_VELOCITY + Z_VELOCITY*Z_VELOCITY;
                PRESSURE = (GAMMA-1.0) * MASS_DENSITY * (SPECIFIC_ENERGY - VEL_SQ_SUM/2.0);
                if(PRESSURE <= E_LIM){PRESSURE = E_LIM;}
        }

        void recalculate_pressure_half(){
                double VEL_SQ_SUM = X_VELOCITY_HALF*X_VELOCITY_HALF + Y_VELOCITY_HALF*Y_VELOCITY_HALF + Z_VELOCITY_HALF*Z_VELOCITY_HALF;
                PRESSURE_HALF = (GAMMA-1.0) * MASS_DENSITY_HALF * (SPECIFIC_ENERGY_HALF - VEL_SQ_SUM/2.0);
                if(PRESSURE_HALF <= E_LIM){PRESSURE_HALF = E_LIM;}
        }

        void calculate_dual(double CONTRIBUTION){DUAL = DUAL + CONTRIBUTION;}

        void prim_to_con(){
                U_VARIABLES[0] = MASS_DENSITY;                          // U0 = mass density
                U_VARIABLES[1] = MASS_DENSITY * X_VELOCITY;             // U1 = x momentum
                U_VARIABLES[2] = MASS_DENSITY * Y_VELOCITY;             // U2 = y momentum
                U_VARIABLES[3] = MASS_DENSITY * Z_VELOCITY;             // U3 = z momentum
                U_VARIABLES[4] = MASS_DENSITY * SPECIFIC_ENERGY;        // U4 = energy density
        }

        void prim_to_con_half(){
                U_HALF[0] = MASS_DENSITY_HALF;                          // U0 = mass density
                U_HALF[1] = MASS_DENSITY_HALF * X_VELOCITY_HALF;             // U1 = x momentum
                U_HALF[2] = MASS_DENSITY_HALF * Y_VELOCITY_HALF;             // U2 = y momentum
                U_HALF[3] = MASS_DENSITY_HALF * Z_VELOCITY_HALF;             // U3 = z momentum
                U_HALF[4] = MASS_DENSITY_HALF * SPECIFIC_ENERGY_HALF;        // U4 = energy density
        }

        // reset half state to new intial state
        void reset_u_half(){
                U_HALF[0] = U_VARIABLES[0];
                U_HALF[1] = U_VARIABLES[1];
                U_HALF[2] = U_VARIABLES[2];
                U_HALF[3] = U_VARIABLES[3];
                U_HALF[4] = U_VARIABLES[4];
        }

        // convert conserved variables to primitive variables
        void con_to_prim(){
                MASS_DENSITY    = U_VARIABLES[0];
                X_VELOCITY      = U_VARIABLES[1]/MASS_DENSITY;
                Y_VELOCITY      = U_VARIABLES[2]/MASS_DENSITY;
                Z_VELOCITY      = U_VARIABLES[3]/MASS_DENSITY;
                SPECIFIC_ENERGY = U_VARIABLES[4]/MASS_DENSITY;
                recalculate_pressure();
                // prim_to_con();
        }

        void con_to_prim_half(){
                MASS_DENSITY_HALF    = U_HALF[0];
                X_VELOCITY_HALF      = U_HALF[1]/MASS_DENSITY_HALF;
                Y_VELOCITY_HALF      = U_HALF[2]/MASS_DENSITY_HALF;
                Z_VELOCITY_HALF      = U_HALF[3]/MASS_DENSITY_HALF;
                SPECIFIC_ENERGY_HALF = U_HALF[4]/MASS_DENSITY_HALF;
                recalculate_pressure_half();
                // prim_to_con_half();
        }

        // reset the changes in primative variables
        void reset_du(){DU[0] = DU[1] = DU[2] = DU[3] = DU[4] = 0.0;}
        void reset_du_half(){DU_HALF[0] = DU_HALF[1] = DU_HALF[2] = DU_HALF[3] = DU_HALF[4] = 0.0;}

        // reset sum for timestep calculation
        void reset_len_vel_sum(){LEN_VEL_SUM = 0.0;}

        // update DU with value from face
        void update_du(double NEW_DU[5]){
                // #pragma omp atomic update
                DU[0] = DU[0] + NEW_DU[0];
                // #pragma omp atomic update
                DU[1] = DU[1] + NEW_DU[1];
                // #pragma omp atomic update
                DU[2] = DU[2] + NEW_DU[2];
                // #pragma omp atomic update
                DU[3] = DU[3] + NEW_DU[3];
                // #pragma omp atomic update
                DU[4] = DU[4] + NEW_DU[4];
        }

        void update_du_half(double NEW_DU[5]){
                // #pragma omp atomic update
                DU_HALF[0] = DU_HALF[0] + NEW_DU[0];
                // #pragma omp atomic update
                DU_HALF[1] = DU_HALF[1] + NEW_DU[1];
                // #pragma omp atomic update
                DU_HALF[2] = DU_HALF[2] + NEW_DU[2];
                // #pragma omp atomic update
                DU_HALF[3] = DU_HALF[3] + NEW_DU[3];
                // #pragma omp atomic update
                DU_HALF[4] = DU_HALF[4] + NEW_DU[4];
        }

        // update fluid varaiables based on sum of changes
        void update_u_variables(){
                // std::cout << "DU =\t" << DU[0] << "\t" << DU[1] << "\t" << DU[2] << "\t" << DU[3] << std::endl;
                U_VARIABLES[0] = U_HALF[0] - DU[0];
                U_VARIABLES[1] = U_HALF[1] - DU[1];
                U_VARIABLES[2] = U_HALF[2] - DU[2];
                U_VARIABLES[3] = U_HALF[3] - DU[3];
                U_VARIABLES[4] = U_HALF[4] - DU[4];
                // if(U_VARIABLES[0] <= MASS_LIM){U_VARIABLES[0] = MASS_LIM;}
                // if(U_VARIABLES[4] <= PRES_LIM){U_VARIABLES[4] = PRES_LIM;}
        }

        void update_u_half(){
                // std::cout << "DU_HALF  =\t" << DU_HALF[0] << "\t" << DU_HALF[1] << "\t" << DU_HALF[2] << "\t" << DU_HALF[3] << std::endl;
                U_HALF[0] = U_VARIABLES[0] - DU_HALF[0];
                U_HALF[1] = U_VARIABLES[1] - DU_HALF[1];
                U_HALF[2] = U_VARIABLES[2] - DU_HALF[2];
                U_HALF[3] = U_VARIABLES[3] - DU_HALF[3];
                U_HALF[4] = U_VARIABLES[4] - DU_HALF[4];
                // if(U_HALF[0] <= MASS_LIM){U_HALF[0] = MASS_LIM;}
                // if(U_HALF[4] <= PRES_LIM){U_HALF[4] = PRES_LIM;}
        }

        // calculate sum of length and velocity (used to calculate dt)
        void update_len_vel_sum(double CONTRIBUTION){
                LEN_VEL_SUM = LEN_VEL_SUM + CONTRIBUTION;
        }

        void check_values(){
                if (U_VARIABLES[0] < M_LIM){
                        // U_VARIABLES[3] = PRES_LIM;
                        // std::cout << "B WARNING: Exiting on negative density\t\t\t";
                        // std::cout << ID << "\tPosition =\t" << X << "\t" << Y << "\tMASS_DENSITY =\t" << U_VARIABLES[0] << std::endl;
                        U_VARIABLES[0] = M_LIM;
                        // std::cout << "B WARNING: Exiting on negative density\t\t\t";
                        // std::cout << ID << "\tPosition =\t" << X << "\t" << Y << "\tMASS_DENSITY =\t" << U_VARIABLES[0] << std::endl;
                        // exit(0);
                }

                if (U_VARIABLES[4] < E_LIM){
                        // U_VARIABLES[3] = PRES_LIM;
                        // std::cout << "B WARNING: Exiting on negative energy\t\t\t";
                        // std::cout << ID << "\tPosition =\t" << X << "\t" << Y << "\tSPECIFIC_ENERGY =\t" << U_VARIABLES[3] << std::endl;
                        U_VARIABLES[4] = E_LIM;
                        // std::cout << "B WARNING: Exiting on negative energy\t\t\t";
                        // std::cout << ID << "\tPosition =\t" << X << "\t" << Y << "\tSPECIFIC_ENERGY =\t" << U_VARIABLES[3] << std::endl;
                        // exit(0);
                }
                return ;
        }
        
        void check_values_half(){
                if (U_HALF[0] < M_LIM){
                        // U_HALF[3] = PRES_LIM;
                        // std::cout << "B WARNING: Exiting on negative half state density\t";
                        // std::cout << ID << "\tPosition =\t" << X << "\t" << Y << "\tMASS_DENSITY_HALF =\t" << U_HALF[0] << std::endl;
                        U_HALF[0] = M_LIM;
                        // std::cout << "B WARNING: Exiting on negative half state density\t";
                        // std::cout << ID << "\tPosition =\t" << X << "\t" << Y << "\tMASS_DENSITY_HALF =\t" << U_HALF[0] << std::endl;
                        // exit(0);
                }

                if (U_HALF[4] < E_LIM){
                        // U_VARIABLES[3] = PRES_LIM;
                        // std::cout << "B WARNING: Exiting on negative half state energy\t\t\t";
                        // std::cout << ID << "\tPosition =\t" << X << "\t" << Y << "\tSPECIFIC_ENERGY_HALF =\t" << U_HALF[3] << std::endl;
                        U_HALF[4] = E_LIM;
                        // std::cout << "B WARNING: Exiting on negative half state energy\t\t\t";
                        // std::cout << ID << "\tPosition =\t" << X << "\t" << Y << "\tSPECIFIC_ENERGY_HALF =\t" << U_HALF[3] << std::endl;
                        // exit(0);
                }
                return ;
        }

        void accelerate(double AX, double AY, double AZ, double DT){
                U_VARIABLES[1] = U_VARIABLES[1] - AX*DT*MASS_DENSITY;
                U_VARIABLES[2] = U_VARIABLES[2] - AY*DT*MASS_DENSITY;
                U_VARIABLES[3] = U_VARIABLES[3] - AZ*DT*MASS_DENSITY;
        }

#ifdef ANALYTIC_GRAVITY
        void calc_newtonian_gravity(double DT){
                // Fixed Plummer potential at (XC,YC)
                double GM,AX,AY,AZ;
                double MPERT = 3.28E+05;
                double XC = 5.0, YC = 5.0, ZC = 5.0;
                double DELTAX = X - XC, DELTAY = Y - YC, DELTAZ = Z - ZC;
                double EPS = 0.145;
                double RAD2 = DELTAX*DELTAX + DELTAY*DELTAY + DELTAZ*DELTAZ;
                GM = GRAV * MPERT / sqrt((RAD2 + EPS*EPS) * (RAD2 + EPS*EPS) * (RAD2 + EPS*EPS));
                AX = DELTAX * GM;
                AY = DELTAY * GM;
                AZ = DELTAZ * GM;
                accelerate(AX, AY, AZ, DT);
                return ;
        }
#endif

        // calculate min timestep this cell requires
        double calc_next_dt(){
                double NEXT_DT;
                NEXT_DT = CFL*3.0*DUAL/LEN_VEL_SUM;
                DT_REQ  = NEXT_DT;
                return NEXT_DT;
        }

        void reset_tbin_local(int INC_TBIN){
                if(INC_TBIN < TBIN_LOCAL){TBIN_LOCAL = INC_TBIN;}
        }

        // maximum value between A and B 
        double max_val(double A, double B){
                if(A>B){
                        return A;
                }else{
                        return B;
                }
        }

};
