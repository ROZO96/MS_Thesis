/*
 * This file was written by Ben Morton (bmorton@ed.ac.uk) and Zhenyu Wu (zhenyu.wu@ed.ac.uk).
 * class containing values and functions associated with triangles
        ID = ID number of triangle
        *VERTEX_0 => pointers to VERTEX 0 of triangle (counter clockwise order)
        *VERTEX_1 => pointers to VERTEX 1 of triangle
        *VERTEX_2 => pointers to VERTEX 2 of triangle
        BOUNDARY  => 0 or 1, denoted whether triangle crosses boundary
        AREA => area of triangle
        X,Y  => vertex coordinates for 0,1,2
        DUAL => area of dual cells corresponding to each vertex
        U_N  => fluid state at each vertex at start of timestep
        U_HALF => intermediate fluid state at each vertex
        FLUC_LDA => nodal residuals for each fluid variable based on initial state (LDA scheme)
        FLUC_N   => nodal residuals for each fluid variable based on initial state (N scheme)
        FLUC_B   => nodal residuals for each fluid variable based on initial state (B scheme)
        FLUC_HALF_LDA => nodal residuals for each fluid variable based on intermediate state (LDA scheme)
        FLUC_HALF_N   => nodal residuals for each fluid variable based on intermediate state (N scheme)
        FLUC_HALF_B   => nodal residuals for each fluid variable based on intermediate state (B scheme)
        PRESSURE = pressure at position of each vertex
        PRESSURE_HALF = pressure of iintermediate state at each vertex
        PHI => element residual
        BETA => distribution coefficient defined by chosen scheme
        MAG => length of normal to each edge
*/

#include <iostream>
#include <numeric>
#include <cmath>

#include "base.h"
#include "inverse.h"
#include "vertex2D.h"

using namespace std;
#pragma once

class TRIANGLE{

private:
        int ID;
        VERTEX *VERTEX_0,*VERTEX_1,*VERTEX_2;

        int BOUNDARY;
        int TBIN;

        double AREA;

        double X[3],Y[3],DUAL[3];
        double X_MOD[3],Y_MOD[3];
        double NORMAL[3][2];

        double U_N[4][3];
        double U_HALF[4][3];

        double FLUC_LDA[4][3],FLUC_N[4][3],FLUC_B[4][3],FLUC_Bx[4][3];
        double FLUC_HALF_LDA[4][3],FLUC_HALF_N[4][3],FLUC_HALF_B[4][3],FLUC_HALF_Bx[4][3];

        double PRESSURE[3], PRESSURE_HALF[3];

        double PHI[4];
        double BETA[4][4][3];

        double DU0[4],DU1[4],DU2[4];
        double DU0_HALF[4],DU1_HALF[4],DU2_HALF[4];

        double MAG[3];

        int PRINT;

public:

        void set_id(int NEW_ID){ID = NEW_ID;}

        void set_vertex_0(VERTEX* NEW_VERTEX){VERTEX_0 = NEW_VERTEX;}
        void set_vertex_1(VERTEX* NEW_VERTEX){VERTEX_1 = NEW_VERTEX;}
        void set_vertex_2(VERTEX* NEW_VERTEX){VERTEX_2 = NEW_VERTEX;}

        void set_boundary(int NEW_BOUNDARY){BOUNDARY = NEW_BOUNDARY;}
        void set_tbin(    int NEW_TBIN){TBIN = NEW_TBIN;}

        int get_id(){return ID;}

        VERTEX* get_vertex_0(){return VERTEX_0;}
        VERTEX* get_vertex_1(){return VERTEX_1;}
        VERTEX* get_vertex_2(){return VERTEX_2;}

        int get_boundary(){return BOUNDARY;}
        int get_tbin(){ return TBIN;}

        double get_un00(){
                U_N[0][0] = VERTEX_0->get_u0();
                return U_N[0][0];
        }
        double get_un01(){
                U_N[0][1] = VERTEX_1->get_u0();
                return U_N[0][1];
        }
        double get_un02(){
                U_N[0][2] = VERTEX_2->get_u0();
                return U_N[0][2];
        }

        void print_triangle_state(){
                std::cout << "U[0] =\t" << U_N[0][0] << "\t" << U_N[0][1] << "\t" << U_N[0][2] << std::endl;
                std::cout << "U[1] =\t" << U_N[1][0] << "\t" << U_N[1][1] << "\t" << U_N[1][2] << std::endl;
                std::cout << "U[2] =\t" << U_N[2][0] << "\t" << U_N[2][1] << "\t" << U_N[2][2] << std::endl;
                std::cout << "U[3] =\t" << U_N[3][0] << "\t" << U_N[3][1] << "\t" << U_N[3][2] << std::endl;

                std::cout << "U_HALF[0] =\t" << U_HALF[0][0] << "\t" << U_HALF[0][1] << "\t" << U_HALF[0][2] << std::endl;
                std::cout << "U_HALF[1] =\t" << U_HALF[1][0] << "\t" << U_HALF[1][1] << "\t" << U_HALF[1][2] << std::endl;
                std::cout << "U_HALF[2] =\t" << U_HALF[2][0] << "\t" << U_HALF[2][1] << "\t" << U_HALF[2][2] << std::endl;
                std::cout << "U_HALF[3] =\t" << U_HALF[3][0] << "\t" << U_HALF[3][1] << "\t" << U_HALF[3][2] << std::endl;
        }

        // import x and y for all vertices
        void setup_positions(){
                X[0] = VERTEX_0->get_x();
                X[1] = VERTEX_1->get_x();
                X[2] = VERTEX_2->get_x();

                Y[0] = VERTEX_0->get_y();
                Y[1] = VERTEX_1->get_y();
                Y[2] = VERTEX_2->get_y();
        }

        // import initial fluid state and pressure for all vertices
        void setup_initial_state(){
                U_N[0][0] = VERTEX_0->get_u0();
                U_N[0][1] = VERTEX_1->get_u0();
                U_N[0][2] = VERTEX_2->get_u0();

                U_N[1][0] = VERTEX_0->get_u1();
                U_N[1][1] = VERTEX_1->get_u1();
                U_N[1][2] = VERTEX_2->get_u1();

                U_N[2][0] = VERTEX_0->get_u2();
                U_N[2][1] = VERTEX_1->get_u2();
                U_N[2][2] = VERTEX_2->get_u2();

                U_N[3][0] = VERTEX_0->get_u3();
                U_N[3][1] = VERTEX_1->get_u3();
                U_N[3][2] = VERTEX_2->get_u3();

                PRESSURE[0] = VERTEX_0->get_pressure();
                PRESSURE[1] = VERTEX_1->get_pressure();
                PRESSURE[2] = VERTEX_2->get_pressure();
        }

        // import intermediate fluid state and pressure for all vertices
        void setup_half_state(){
                U_HALF[0][0] = VERTEX_0->get_u0_half();
                U_HALF[0][1] = VERTEX_1->get_u0_half();
                U_HALF[0][2] = VERTEX_2->get_u0_half();

                U_HALF[1][0] = VERTEX_0->get_u1_half();
                U_HALF[1][1] = VERTEX_1->get_u1_half();
                U_HALF[1][2] = VERTEX_2->get_u1_half();

                U_HALF[2][0] = VERTEX_0->get_u2_half();
                U_HALF[2][1] = VERTEX_1->get_u2_half();
                U_HALF[2][2] = VERTEX_2->get_u2_half();

                U_HALF[3][0] = VERTEX_0->get_u3_half();
                U_HALF[3][1] = VERTEX_1->get_u3_half();
                U_HALF[3][2] = VERTEX_2->get_u3_half();

                PRESSURE_HALF[0] = VERTEX_0->get_pressure_half();
                PRESSURE_HALF[1] = VERTEX_1->get_pressure_half();
                PRESSURE_HALF[2] = VERTEX_2->get_pressure_half();
        }

        //**********************************************************************************************************************

        // Calculate first half timestep change, passing change to vertice
        void calculate_first_half(double T, double DT,double P_DIFF =0){
                int i,j,m,p;
                double INFLOW[4][4][3][3];
                double C_SOUND[3];

                // Import conditions and positions of vertices

                setup_positions();
                setup_initial_state();

#ifdef CLOSED
                for(m=0;m<3;++m){
                        if(X[m] > 0.99*SIDE_LENGTH_X or X[m] < 0.01*SIDE_LENGTH_X or Y[m] > 0.99*SIDE_LENGTH_Y or Y[m] < 0.01*SIDE_LENGTH_Y){
                                return ;
                        }
                }
#endif

#ifdef DEBUG
                std::cout << "-- FIRST  -------------------------------------------------------" << std::endl;
                std::cout << "Time     =\t" << T << std::endl;
                std::cout << "0 =\t" << X[0] << "\t" << Y[0] << std::endl;
                std::cout << "1 =\t" << X[1] << "\t" << Y[1] << std::endl;
                std::cout << "2 =\t" << X[2] << "\t" << Y[2] << std::endl;
                std::cout << "State    =" << "\trho" << "\tx_mom" << "\ty_mom" << "\tenergy" << std::endl;
                for(i=0;i<3;i++){std::cout << i << " =\t" << U_N[0][i] << "\t" << U_N[1][i] << "\t" << U_N[2][i] << "\t" << U_N[3][i] << std::endl;}
                std::cout << "Pressure =\t" << PRESSURE[0] << "\t" << PRESSURE[1] << "\t" << PRESSURE[2] << std::endl;
#endif

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Calculate inflow parameters

                double H[3];
                double RHO,C,U,U_C,V,V_C,H_AVG,H_C,ALPHA,ALPHA_C,W;
                double Z[4][3];
                double VALUE1,VALUE2,VALUE3,VALUE4,VALUE12,VALUE123;
                double PRESSURE_AVG;
                double LAMBDA[4][3],LAMBDA_PLUS[4][3],LAMBDA_MINUS[4][3];
                double N_X[3],N_Y[3];

                double Z_BAR[4],W_HAT[4][3];

                // Construct Roe vector Z

                for(m=0;m<3;++m){
                        Z[0][m] = sqrt(U_N[0][m]);
                        Z[1][m] = U_N[1][m]/Z[0][m];
                        Z[2][m] = U_N[2][m]/Z[0][m];
                        Z[3][m] = (U_N[3][m] + PRESSURE[m])/Z[0][m];

                        N_X[m]  = NORMAL[m][0];
                        N_Y[m]  = NORMAL[m][1];

                        H[m] = (U_N[3][m] + PRESSURE[m])/U_N[0][m];
#ifdef DEBUG
                        std::cout << "Z =\t" << Z[0][m] << "\t" << Z[1][m] << "\t" << Z[2][m] << "\t" << Z[3][m] << std::endl;
#endif
                }

                for(i=0; i<4; ++i){Z_BAR[i] = (Z[i][0] + Z[i][1] + Z[i][2])/3.0;}

#ifdef DEBUG
                for(i=0; i<4; ++i){std::cout << "Z_BAR " << i << " =\t" << Z_BAR[i] << std::endl;}
                std::cout << std::endl;
#endif


                for(m=0; m<3; ++m){
                        W_HAT[0][m] =  2.0*Z_BAR[0]*Z[0][m];
                        W_HAT[1][m] =  Z_BAR[1]*Z[0][m] + Z_BAR[0]*Z[1][m];
                        W_HAT[2][m] =  Z_BAR[2]*Z[0][m] + Z_BAR[0]*Z[2][m];
                        W_HAT[3][m] = (Z_BAR[3]*Z[0][m] + GAMMA_1*Z_BAR[1]*Z[1][m] + GAMMA_1*Z_BAR[2]*Z[2][m] + Z_BAR[0]*Z[3][m])/GAMMA;

                        //std::cout << Z_BAR[1] << "\t" << Z[0][m] << "\t" << Z_BAR[0] << "\t" << Z[2][m] << std::endl;

#ifdef DEBUG
                        for(i=0;i<4;++i){std::cout << "W_HAT " << i << "\t" << m << " =\t" << W_HAT[i][m] << std::endl;}
                        std::cout << std::endl;
#endif
                }

                // Construct average state for element

                RHO   = pow((sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]))/3.0, 2);
                U     = (sqrt(U_N[0][0])*U_N[1][0]/U_N[0][0] + sqrt(U_N[0][1])*U_N[1][1]/U_N[0][1] + sqrt(U_N[0][2])*U_N[1][2]/U_N[0][2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));        // U now represents x velocity
                V     = (sqrt(U_N[0][0])*U_N[2][0]/U_N[0][0] + sqrt(U_N[0][1])*U_N[2][1]/U_N[0][1] + sqrt(U_N[0][2])*U_N[2][2]/U_N[0][2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));        // V represents y velocity
                H_AVG = (sqrt(U_N[0][0])*H[0] + sqrt(U_N[0][1])*H[1] + sqrt(U_N[0][2])*H[2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));
               
                // E     = (sqrt(U_N[0][0])*H[0]/U_N[0][0] + sqrt(U_N[0][1])*H[1]/U_N[0][1] + sqrt(U_N[0][2])*H[2]/U_N[0][2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));

                PRESSURE_AVG = (PRESSURE[0] + PRESSURE[1] + PRESSURE[2])/3.0;
                C = sqrt((GAMMA-1.0) * H_AVG - (GAMMA-1.0) * (U*U + V*V)/2.0);
                if(C <= C_LIM){C = C_LIM;}

#ifdef DEBUG
                std::cout << "PRESSURE_AVG =\t" << PRESSURE_AVG << std::endl;
                std::cout << "C_SOUND_AVG  =\t" << C << std::endl;
#endif

                // Reassign variables to local equivalents

                U_C = U/C;
                V_C = V/C;
                H_C = H_AVG/C;

                ALPHA   = GAMMA_1*(U*U + V*V)/2.0;
                ALPHA_C = ALPHA/C;

#ifdef DEBUG
                std::cout << "U =\t" << U << "\tU_C =\t" << U_C << std::endl;
                std::cout << "V =\t" << V << "\tV_C =\t" << V_C << std::endl;
                std::cout << "H =\t" << H_AVG << "\tH_C =\t" << H_C << std::endl;
                std::cout << "ALPHA =\t" << ALPHA << "\tALPHA_C =\t" << ALPHA_C << std::endl;
                std::cout << std::endl;
#endif

                // Calculate K+,K- and K matrices for each vertex i,j,k

                for(m=0;m<3;++m){

                        W = U*N_X[m] + V*N_Y[m];

#ifdef DEBUG
                        std::cout << "W =\t" << W << std::endl;
#endif

                        LAMBDA[0][m] = W + C;
                        LAMBDA[1][m] = W - C;
                        LAMBDA[2][m] = W;
                        LAMBDA[3][m] = W;

                        for(i=0;i<4;++i){
                                LAMBDA_PLUS[i][m]  = max_val(0.0,LAMBDA[i][m]);
                                LAMBDA_MINUS[i][m] = min_val(0.0,LAMBDA[i][m]);
                        } 

                        for(p=0;p<3;++p){
                                if(p==0){       // Identify and select positive eigenvalues
                                        VALUE1 = LAMBDA_PLUS[0][m];
                                        VALUE2 = LAMBDA_PLUS[1][m];
                                        VALUE3 = LAMBDA_PLUS[2][m];
                                        VALUE4 = LAMBDA_PLUS[3][m];
#ifdef DEBUG
                                        std::cout << "K+" << std::endl;
#endif
                                }else if(p==1){ // Identify and select negative eigenvalues
                                        VALUE1 = LAMBDA_MINUS[0][m];
                                        VALUE2 = LAMBDA_MINUS[1][m];
                                        VALUE3 = LAMBDA_MINUS[2][m];
                                        VALUE4 = LAMBDA_MINUS[3][m];
#ifdef DEBUG
                                        std::cout << "K-" << std::endl;
#endif
                                }else{          // Select all eigenvalues
                                        VALUE1 = LAMBDA[0][m];
                                        VALUE2 = LAMBDA[1][m];
                                        VALUE3 = LAMBDA[2][m];
                                        VALUE4 = LAMBDA[3][m];
#ifdef DEBUG
                                        std::cout << "K" << std::endl;
#endif
                                }

                                VALUE12  = (VALUE1 - VALUE2)/2.0;
                                VALUE123 = (VALUE1 + VALUE2 - 2.0*VALUE3)/2.0;

                                INFLOW[0][0][m][p] = 0.5*MAG[m]*(ALPHA_C*VALUE123/C - W*VALUE12/C + VALUE3);
                                INFLOW[0][1][m][p] = 0.5*MAG[m]*(-1.0*GAMMA_1*U_C*VALUE123/C + N_X[m]*VALUE12/C);
                                INFLOW[0][2][m][p] = 0.5*MAG[m]*(-1.0*GAMMA_1*V_C*VALUE123/C + N_Y[m]*VALUE12/C);
                                INFLOW[0][3][m][p] = 0.5*MAG[m]*(GAMMA_1*VALUE123/(C*C));

                                INFLOW[1][0][m][p] = 0.5*MAG[m]*((ALPHA_C*U_C - W*N_X[m])*VALUE123 + (ALPHA_C*N_X[m] - U_C*W)*VALUE12);
                                INFLOW[1][1][m][p] = 0.5*MAG[m]*((N_X[m]*N_X[m] - GAMMA_1*U_C*U_C)*VALUE123 - (GAMMA_2*U_C*N_X[m]*VALUE12) + VALUE3);
                                INFLOW[1][2][m][p] = 0.5*MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*U_C*V_C)*VALUE123 + (U_C*N_Y[m] - GAMMA_1*V_C*N_X[m])*VALUE12);
                                INFLOW[1][3][m][p] = 0.5*MAG[m]*(GAMMA_1*U_C*VALUE123/C + GAMMA_1*N_X[m]*VALUE12/C);

                                INFLOW[2][0][m][p] = 0.5*MAG[m]*((ALPHA_C*V_C - W*N_Y[m])*VALUE123 + (ALPHA_C*N_Y[m] - V_C*W)*VALUE12);
                                INFLOW[2][1][m][p] = 0.5*MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*U_C*V_C)*VALUE123 + (V_C*N_X[m] - GAMMA_1*U_C*N_Y[m])*VALUE12);
                                INFLOW[2][2][m][p] = 0.5*MAG[m]*((N_Y[m]*N_Y[m] - GAMMA_1*V_C*V_C)*VALUE123 - (GAMMA_2*V_C*N_Y[m]*VALUE12) + VALUE3);
                                INFLOW[2][3][m][p] = 0.5*MAG[m]*(GAMMA_1*V_C*VALUE123/C + GAMMA_1*N_Y[m]*VALUE12/C);

                                INFLOW[3][0][m][p] = 0.5*MAG[m]*((ALPHA_C*H_C - W*W)*VALUE123 + W*(ALPHA_C - H_C)*VALUE12);
                                INFLOW[3][1][m][p] = 0.5*MAG[m]*((W*N_X[m] - U - ALPHA_C*U_C)*VALUE123 + (H_C*N_X[m] - GAMMA_1*U_C*W)*VALUE12);
                                INFLOW[3][2][m][p] = 0.5*MAG[m]*((W*N_Y[m] - V - ALPHA_C*V_C)*VALUE123 + (H_C*N_Y[m] - GAMMA_1*V_C*W)*VALUE12);
                                INFLOW[3][3][m][p] = 0.5*MAG[m]*(GAMMA_1*H_C*VALUE123/C + GAMMA_1*W*VALUE12/C + VALUE3);
                        }
#ifdef DEBUG
                        for(i=0; i<4; ++i){
                                for (j=0; j<4; ++j){
                                        std::cout << INFLOW[i][j][m][0] << "\t";
                                }
                                std::cout << std::endl;
                        }
                        std::cout << std::endl;
#endif
                }


#ifdef DEBUG
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[0][0] << "\t" << LAMBDA_PLUS[0][1] << "\t" << LAMBDA_PLUS[0][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[1][0] << "\t" << LAMBDA_PLUS[1][1] << "\t" << LAMBDA_PLUS[1][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[2][0] << "\t" << LAMBDA_PLUS[2][1] << "\t" << LAMBDA_PLUS[2][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[3][0] << "\t" << LAMBDA_PLUS[3][1] << "\t" << LAMBDA_PLUS[3][2] << std::endl;

                std::cout << "Lambda - =\t" << LAMBDA_MINUS[0][0] << "\t" << LAMBDA_MINUS[0][1] << "\t" << LAMBDA_MINUS[0][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[1][0] << "\t" << LAMBDA_MINUS[1][1] << "\t" << LAMBDA_MINUS[1][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[2][0] << "\t" << LAMBDA_MINUS[2][1] << "\t" << LAMBDA_MINUS[2][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[3][0] << "\t" << LAMBDA_MINUS[3][1] << "\t" << LAMBDA_MINUS[3][2] << std::endl;

                std::cout << "Lambda   =\t" << LAMBDA[0][0] << "\t" << LAMBDA[0][1] << "\t" << LAMBDA[0][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[1][0] << "\t" << LAMBDA[1][1] << "\t" << LAMBDA[1][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[2][0] << "\t" << LAMBDA[2][1] << "\t" << LAMBDA[2][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[3][0] << "\t" << LAMBDA[3][1] << "\t" << LAMBDA[3][2] << std::endl;
#endif

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                for(i=0;i<4;++i){
                        PHI[i] = 0.0;
#ifdef DEBUG
                        std::cout << "Calculating Phi " << i << std::endl;
#endif
                        for(m=0;m<3;++m){
                                PHI[i] += INFLOW[i][0][m][2]*W_HAT[0][m] + INFLOW[i][1][m][2]*W_HAT[1][m] + INFLOW[i][2][m][2]*W_HAT[2][m] + INFLOW[i][3][m][2]*W_HAT[3][m];
                        }
                }

#ifdef DEBUG
                std::cout << "PHI =\t" << PHI[0] << "\t" << PHI[1] << "\t" << PHI[2] << "\t" << PHI[3] << std::endl;
#endif

                double INFLOW_MINUS_SUM[4][4];

                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                INFLOW_MINUS_SUM[i][j] = 0.0;
                                for(m=0;m<3;++m){
                                        INFLOW_MINUS_SUM[i][j] += INFLOW[i][j][m][1];
                                }
                        }
                }

                mat_inv(&INFLOW_MINUS_SUM[0][0],4,X[0],Y[0],ID,1);

                // std::cout << "Post-inversion =" << std::endl;

#ifdef DEBUG
                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                std::cout << INFLOW_MINUS_SUM[i][j] << "\t";
                        }
                        std::cout << std::endl;
                }
#endif

                // Calculate spatial splitting for first half timestep

#if defined(LDA_SCHEME) or defined(BLENDED) or defined(Bx_SCHEME)

#ifdef DEBUG
                std::cout << "BETA_0 =" << std::endl;
#endif

                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                for(m=0;m<3;++m){
                                        BETA[i][j][m] = -1.0*(INFLOW[i][0][m][0] * INFLOW_MINUS_SUM[0][j] + INFLOW[i][1][m][0] * INFLOW_MINUS_SUM[1][j] + INFLOW[i][2][m][0] * INFLOW_MINUS_SUM[2][j] + INFLOW[i][3][m][0] * INFLOW_MINUS_SUM[3][j]);
                                }
                        }
#ifdef DEBUG
                        std::cout << BETA[i][0][0] << "\t" << BETA[i][1][0] << "\t" << BETA[i][2][0] << "\t" << BETA[i][3][0] << std::endl;
#endif
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                FLUC_LDA[i][m] = BETA[i][0][m] * PHI[0] + BETA[i][1][m] * PHI[1] + BETA[i][2][m] * PHI[2] + BETA[i][3][m] * PHI[3];
                                // std::cout << FLUC_LDA[i][m] << std::endl;
                        }
                }
#endif


#if defined(N_SCHEME) or defined(BLENDED) or defined(Bx_SCHEME)
                double BRACKET[4][3];
                double KZ_SUM[4];

                for(i=0;i<4;++i){
                        KZ_SUM[i] = 0.0;
                        for(m=0;m<3;++m){
                                KZ_SUM[i] += INFLOW[i][0][m][1] * W_HAT[0][m] + INFLOW[i][1][m][1] * W_HAT[1][m] + INFLOW[i][2][m][1] * W_HAT[2][m] + INFLOW[i][3][m][1] * W_HAT[3][m];
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                BRACKET[i][m] = W_HAT[i][m] - (INFLOW_MINUS_SUM[i][0]*KZ_SUM[0] + INFLOW_MINUS_SUM[i][1]*KZ_SUM[1] + INFLOW_MINUS_SUM[i][2]*KZ_SUM[2] + INFLOW_MINUS_SUM[i][3]*KZ_SUM[3]);
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                FLUC_N[i][m] = INFLOW[i][0][m][0]*BRACKET[0][m] + INFLOW[i][1][m][0]*BRACKET[1][m] + INFLOW[i][2][m][0]*BRACKET[2][m] + INFLOW[i][3][m][0]*BRACKET[3][m];
                                // std::cout << FLUC_N[i][m] << std::endl;
                        }
                }
#endif

#ifdef BLENDED
                double THETA_E[4][4];
                double IDENTITY[4][4];
                double SUM_FLUC_N[4];

                THETA_E[0][0] = IDENTITY[0][0] = 1.0;
                THETA_E[0][1] = IDENTITY[0][1] = 0.0;
                THETA_E[0][2] = IDENTITY[0][2] = 0.0;
                THETA_E[0][3] = IDENTITY[0][3] = 0.0;

                THETA_E[1][0] = IDENTITY[1][0] = 0.0;
                THETA_E[1][1] = IDENTITY[1][1] = 1.0;
                THETA_E[1][2] = IDENTITY[1][2] = 0.0;
                THETA_E[1][3] = IDENTITY[1][3] = 0.0;

                THETA_E[2][0] = IDENTITY[2][0] = 0.0;
                THETA_E[2][1] = IDENTITY[2][1] = 0.0;
                THETA_E[2][2] = IDENTITY[2][2] = 1.0;
                THETA_E[2][3] = IDENTITY[2][3] = 0.0;

                THETA_E[3][0] = IDENTITY[3][0] = 0.0;
                THETA_E[3][1] = IDENTITY[3][1] = 0.0;
                THETA_E[3][2] = IDENTITY[3][2] = 0.0;
                THETA_E[3][3] = IDENTITY[3][3] = 1.0;

                for(i=0;i<4;i++){
                        SUM_FLUC_N[i] = abs(FLUC_N[i][0]) + abs(FLUC_N[i][1]) + abs(FLUC_N[i][2]);
                        THETA_E[i][i] = abs(PHI[i])/SUM_FLUC_N[i];
                        FLUC_B[i][0] = THETA_E[i][i]*FLUC_N[i][0] + (IDENTITY[i][i] - THETA_E[i][i])*FLUC_LDA[i][0];
                        FLUC_B[i][1] = THETA_E[i][i]*FLUC_N[i][1] + (IDENTITY[i][i] - THETA_E[i][i])*FLUC_LDA[i][1];
                        FLUC_B[i][2] = THETA_E[i][i]*FLUC_N[i][2] + (IDENTITY[i][i] - THETA_E[i][i])*FLUC_LDA[i][2];
                }
#endif

                DUAL[0] = VERTEX_0->get_dual();
                DUAL[1] = VERTEX_1->get_dual();
                DUAL[2] = VERTEX_2->get_dual();
                
#ifdef Bx_SCHEME              
             
     		double  Xc, Yc;
     		double interpol[3][5];
     		double grad_p[2];
     		
     		
     		Xc=std::accumulate(std::begin(X), std::end(X), 0);
     		Yc=std::accumulate(std::begin(Y), std::end(Y), 0);
     		Xc/=3;
     		Yc/=3;
     		
     		cell_interpolation_N(interpol);
     		
             	grad_p[0]=interpol[0][4];
             	grad_p[1]=interpol[1][4];
             	
             	double Vxc, Vyc,Vmag;
             	
             	Vxc=interpol[0][1]*Xc+interpol[1][1]*Yc+interpol[2][1];
           	Vyc=interpol[0][2]*Xc+interpol[1][2]*Yc+interpol[2][2];
             	Vmag=sqrt(Vxc*Vxc+Vyc*Vyc);
             	
             	double sc=0;
             	double h=sqrt(4*AREA/3.1415);
             	
             	if(Vmag==0 ||P_DIFF==0){sc=0;}
             	else{sc=(h/(Vmag*P_DIFF))*((grad_p[0]*Vxc)+(grad_p[1]*Vyc));}
             	
             	if (sc<0){sc=0;}
             	if (sc>1){sc=1;}
             
             	
             	double THETA=sc*sc*h;
             	
             	
             	
             	for(i=0;i<4;i++){
                        FLUC_Bx[i][0] = THETA*FLUC_N[i][0] + (1 - THETA)*FLUC_LDA[i][0];
                        FLUC_Bx[i][1] = THETA*FLUC_N[i][1] + (1 - THETA)*FLUC_LDA[i][1];
                        FLUC_Bx[i][2] = THETA*FLUC_N[i][2] + (1 - THETA)*FLUC_LDA[i][2];
                }
             	
#endif

#ifdef LDA_SCHEME
                for(i=0;i<4;i++){
                        DU0_HALF[i] = -1.0*DT*FLUC_LDA[i][0]/DUAL[0];
                        DU1_HALF[i] = -1.0*DT*FLUC_LDA[i][1]/DUAL[1];
                        DU2_HALF[i] = -1.0*DT*FLUC_LDA[i][2]/DUAL[2];
                }
#endif

#ifdef N_SCHEME
                for(i=0;i<4;i++){
                        DU0_HALF[i] = -1.0*DT*FLUC_N[i][0]/DUAL[0];
                        DU1_HALF[i] = -1.0*DT*FLUC_N[i][1]/DUAL[1];
                        DU2_HALF[i] = -1.0*DT*FLUC_N[i][2]/DUAL[2];
                }
#endif

#ifdef BLENDED 
                for(i=0;i<4;i++){
                        DU0_HALF[i] = -1.0*DT*FLUC_B[i][0]/DUAL[0];
                        DU1_HALF[i] = -1.0*DT*FLUC_B[i][1]/DUAL[1];
                        DU2_HALF[i] = -1.0*DT*FLUC_B[i][2]/DUAL[2];
                }

#endif

#ifdef Bx_SCHEME
                for(i=0;i<4;i++){
                        DU0_HALF[i] = -1.0*DT*FLUC_Bx[i][0]/DUAL[0];
                        DU1_HALF[i] = -1.0*DT*FLUC_Bx[i][1]/DUAL[1];
                        DU2_HALF[i] = -1.0*DT*FLUC_Bx[i][2]/DUAL[2];
                }

#endif



        }

        void pass_update_half(){
                VERTEX_0->update_du_half(DU0_HALF);
                VERTEX_1->update_du_half(DU1_HALF);
                VERTEX_2->update_du_half(DU2_HALF);
        }

        //**********************************************************************************************************************


        void calculate_second_half(double T, double DT, double P_DIFF =0){
                int i,j,m,p;
                double INFLOW[4][4][3][3];

                setup_half_state();

#ifdef FIRST_ORDER
                for(i=0;i<4;i++){
                        DU0[i] = 0.0;
                        DU1[i] = 0.0;
                        DU2[i] = 0.0;
                }

                VERTEX_0->update_du(DU0);
                VERTEX_1->update_du(DU1);
                VERTEX_2->update_du(DU2);

                return ;
#endif

#ifdef CLOSED
                for(m=0;m<3;++m){
                        if(X[m] > 0.99*SIDE_LENGTH_X or X[m] < 0.01*SIDE_LENGTH_X or Y[m] > 0.99*SIDE_LENGTH_Y or Y[m] < 0.01*SIDE_LENGTH_Y){
                                return ;
                        }
                }
#endif

#ifdef DEBUG
                std::cout << "-- SECOND -------------------------------------------------------" << std::endl;
                std::cout << "Time     =\t" << T << std::endl;
                std::cout << "0 =\t" << X[0] << "\t" << Y[0] << std::endl;
                std::cout << "1 =\t" << X[1] << "\t" << Y[1] << std::endl;
                std::cout << "2 =\t" << X[2] << "\t" << Y[2] << std::endl;
                std::cout << "Half State    =" << "\trho" << "\tx_mom" << "\ty_mom" << "\tenergy" << std::endl;
                for(i=0;i<3;i++){std::cout << i << " =\t" << U_HALF[0][i] << "\t" << U_HALF[1][i] << "\t" << U_HALF[2][i] << "\t" << U_HALF[3][i] << std::endl;}
                std::cout << "Pressure =\t" << PRESSURE_HALF[0] << "\t" << PRESSURE_HALF[1] << "\t" << PRESSURE_HALF[2] << std::endl;
#endif

                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Calculate inflow parameters

                double H[3];
                double RHO,C,U,U_C,V,V_C,H_AVG,H_C,ALPHA,ALPHA_C,W;
                double Z[4][3];
                double VALUE1,VALUE2,VALUE3,VALUE4,VALUE12,VALUE123;
                double PRESSURE_AVG,C_SOUND_AVG;
                double LAMBDA[4][3],LAMBDA_PLUS[4][3],LAMBDA_MINUS[4][3];
                double N_X[3],N_Y[3];

                double Z_BAR[4],W_HAT[4][3];

                // Construct Roe vector Z

                for(m=0;m<3;++m){
                        Z[0][m] = sqrt(U_HALF[0][m]);
                        Z[1][m] = U_HALF[1][m]/Z[0][m];
                        Z[2][m] = U_HALF[2][m]/Z[0][m];
                        Z[3][m] = (U_HALF[3][m] + PRESSURE_HALF[m])/Z[0][m];

                        N_X[m]  = NORMAL[m][0];
                        N_Y[m]  = NORMAL[m][1];

                        H[m] = (U_HALF[3][m] + PRESSURE_HALF[m])/U_HALF[0][m];
#ifdef DEBUG
                        std::cout << "Z =\t" << Z[0][m] << "\t" << Z[1][m] << "\t" << Z[2][m] << "\t" << Z[3][m] << std::endl;
#endif
                }

                for(i=0; i<4; ++i){Z_BAR[i] = (Z[i][0] + Z[i][1] + Z[i][2])/3.0;}

#ifdef DEBUG
                for(i=0; i<4; ++i){std::cout << "Z_BAR " << i << " =\t" << Z_BAR[i] << std::endl;}
                std::cout << std::endl;
#endif


                for(m=0; m<3; ++m){
                        W_HAT[0][m] =  2.0*Z_BAR[0]*Z[0][m];
                        W_HAT[1][m] =  Z_BAR[1]*Z[0][m] + Z_BAR[0]*Z[1][m];
                        W_HAT[2][m] =  Z_BAR[2]*Z[0][m] + Z_BAR[0]*Z[2][m];
                        W_HAT[3][m] = (Z_BAR[3]*Z[0][m] + GAMMA_1*Z_BAR[1]*Z[1][m] + GAMMA_1*Z_BAR[2]*Z[2][m] + Z_BAR[0]*Z[3][m])/GAMMA;

                        //std::cout << Z_BAR[1] << "\t" << Z[0][m] << "\t" << Z_BAR[0] << "\t" << Z[2][m] << std::endl;

#ifdef DEBUG
                        for(i=0;i<4;++i){std::cout << "W_HAT " << i << "\t" << m << " =\t" << W_HAT[i][m] << std::endl;}
                        std::cout << std::endl;
#endif
                }

                // Construct average state for element

                RHO   = pow((sqrt(U_HALF[0][0]) + sqrt(U_HALF[0][1]) + sqrt(U_HALF[0][2]))/3.0, 2);
                U     = (sqrt(U_HALF[0][0])*U_HALF[1][0]/U_HALF[0][0] + sqrt(U_HALF[0][1])*U_HALF[1][1]/U_HALF[0][1] + sqrt(U_HALF[0][2])*U_HALF[1][2]/U_HALF[0][2])/(sqrt(U_HALF[0][0]) + sqrt(U_HALF[0][1]) + sqrt(U_HALF[0][2]));        // U now represents x velocity
                V     = (sqrt(U_HALF[0][0])*U_HALF[2][0]/U_HALF[0][0] + sqrt(U_HALF[0][1])*U_HALF[2][1]/U_HALF[0][1] + sqrt(U_HALF[0][2])*U_HALF[2][2]/U_HALF[0][2])/(sqrt(U_HALF[0][0]) + sqrt(U_HALF[0][1]) + sqrt(U_HALF[0][2]));        // V represents y velocity
                H_AVG = (sqrt(U_HALF[0][0])*H[0] + sqrt(U_HALF[0][1])*H[1] + sqrt(U_HALF[0][2])*H[2])/(sqrt(U_HALF[0][0]) + sqrt(U_HALF[0][1]) + sqrt(U_HALF[0][2]));
               
                // E     = (sqrt(U_N[0][0])*H[0]/U_N[0][0] + sqrt(U_N[0][1])*H[1]/U_N[0][1] + sqrt(U_N[0][2])*H[2]/U_N[0][2])/(sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]));

                PRESSURE_AVG = (PRESSURE_HALF[0] + PRESSURE_HALF[1] + PRESSURE_HALF[2])/3.0;
                C = sqrt((GAMMA-1.0) * H_AVG - (GAMMA-1.0) * (U*U + V*V)/2.0);
                if(C <= C_LIM){C = C_LIM;}

#ifdef DEBUG
                std::cout << "PRESSURE_AVG =\t" << PRESSURE_AVG << std::endl;
                std::cout << "C_SOUND_AVG  =\t" << C_SOUND_AVG  << std::endl;
#endif

                // Reassign variables to local equivalents

                U_C = U/C;
                V_C = V/C;
                H_C = H_AVG/C;

                ALPHA   = GAMMA_1*(U*U + V*V)/2.0;
                ALPHA_C = ALPHA/C;

#ifdef DEBUG
                std::cout << "U =\t" << U << "\tU_C =\t" << U_C << std::endl;
                std::cout << "V =\t" << V << "\tV_C =\t" << V_C << std::endl;
                std::cout << "H =\t" << H_AVG << "\tH_C =\t" << H_C << std::endl;
                std::cout << "ALPHA =\t" << ALPHA << "\tALPHA_C =\t" << ALPHA_C << std::endl;
                std::cout << std::endl;
#endif

                // Calculate K+,K- and K matrices for each vertex i,j,k

                for(m=0;m<3;++m){

                        W = U*N_X[m] + V*N_Y[m];

#ifdef DEBUG
                        std::cout << "W =\t" << W << std::endl;
#endif

                        LAMBDA[0][m] = W + C;
                        LAMBDA[1][m] = W - C;
                        LAMBDA[2][m] = W;
                        LAMBDA[3][m] = W;

                        for(i=0;i<4;++i){
                                LAMBDA_PLUS[i][m]  = max_val(0.0,LAMBDA[i][m]);
                                LAMBDA_MINUS[i][m] = min_val(0.0,LAMBDA[i][m]);
                        } 

                        for(p=0;p<3;++p){
                                if(p==0){       // Identify and select positive eigenvalues
                                        VALUE1 = LAMBDA_PLUS[0][m];
                                        VALUE2 = LAMBDA_PLUS[1][m];
                                        VALUE3 = LAMBDA_PLUS[2][m];
                                        VALUE4 = LAMBDA_PLUS[3][m];
#ifdef DEBUG
                                        std::cout << "K+" << std::endl;
#endif
                                }else if(p==1){ // Identify and select negative eigenvalues
                                        VALUE1 = LAMBDA_MINUS[0][m];
                                        VALUE2 = LAMBDA_MINUS[1][m];
                                        VALUE3 = LAMBDA_MINUS[2][m];
                                        VALUE4 = LAMBDA_MINUS[3][m];
#ifdef DEBUG
                                        std::cout << "K-" << std::endl;
#endif
                                }else{          // Select all eigenvalues
                                        VALUE1 = LAMBDA[0][m];
                                        VALUE2 = LAMBDA[1][m];
                                        VALUE3 = LAMBDA[2][m];
                                        VALUE4 = LAMBDA[3][m];
#ifdef DEBUG
                                        std::cout << "K" << std::endl;
#endif
                                }


                                VALUE12  = (VALUE1 - VALUE2)/2.0;
                                VALUE123 = (VALUE1 + VALUE2 - 2.0*VALUE3)/2.0;

                                INFLOW[0][0][m][p] = 0.5*MAG[m]*(ALPHA_C*VALUE123/C - W*VALUE12/C + VALUE3);
                                INFLOW[0][1][m][p] = 0.5*MAG[m]*(-1.0*GAMMA_1*U_C*VALUE123/C + N_X[m]*VALUE12/C);
                                INFLOW[0][2][m][p] = 0.5*MAG[m]*(-1.0*GAMMA_1*V_C*VALUE123/C + N_Y[m]*VALUE12/C);
                                INFLOW[0][3][m][p] = 0.5*MAG[m]*(GAMMA_1*VALUE123/(C*C));

                                INFLOW[1][0][m][p] = 0.5*MAG[m]*((ALPHA_C*U_C - W*N_X[m])*VALUE123 + (ALPHA_C*N_X[m] - U_C*W)*VALUE12);
                                INFLOW[1][1][m][p] = 0.5*MAG[m]*((N_X[m]*N_X[m] - GAMMA_1*U_C*U_C)*VALUE123 - (GAMMA_2*U_C*N_X[m]*VALUE12) + VALUE3);
                                INFLOW[1][2][m][p] = 0.5*MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*U_C*V_C)*VALUE123 + (U_C*N_Y[m] - GAMMA_1*V_C*N_X[m])*VALUE12);
                                INFLOW[1][3][m][p] = 0.5*MAG[m]*(GAMMA_1*U_C*VALUE123/C + GAMMA_1*N_X[m]*VALUE12/C);

                                INFLOW[2][0][m][p] = 0.5*MAG[m]*((ALPHA_C*V_C - W*N_Y[m])*VALUE123 + (ALPHA_C*N_Y[m] - V_C*W)*VALUE12);
                                INFLOW[2][1][m][p] = 0.5*MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*U_C*V_C)*VALUE123 + (V_C*N_X[m] - GAMMA_1*U_C*N_Y[m])*VALUE12);
                                INFLOW[2][2][m][p] = 0.5*MAG[m]*((N_Y[m]*N_Y[m] - GAMMA_1*V_C*V_C)*VALUE123 - (GAMMA_2*V_C*N_Y[m]*VALUE12) + VALUE3);
                                INFLOW[2][3][m][p] = 0.5*MAG[m]*(GAMMA_1*V_C*VALUE123/C + GAMMA_1*N_Y[m]*VALUE12/C);

                                INFLOW[3][0][m][p] = 0.5*MAG[m]*((ALPHA_C*H_C - W*W)*VALUE123 + W*(ALPHA_C - H_C)*VALUE12);
                                INFLOW[3][1][m][p] = 0.5*MAG[m]*((W*N_X[m] - U - ALPHA_C*U_C)*VALUE123 + (H_C*N_X[m] - GAMMA_1*U_C*W)*VALUE12);
                                INFLOW[3][2][m][p] = 0.5*MAG[m]*((W*N_Y[m] - V - ALPHA_C*V_C)*VALUE123 + (H_C*N_Y[m] - GAMMA_1*V_C*W)*VALUE12);
                                INFLOW[3][3][m][p] = 0.5*MAG[m]*(GAMMA_1*H_C*VALUE123/C + GAMMA_1*W*VALUE12/C + VALUE3);

#ifdef DEBUG
                                for(i=0; i<4; ++i){
                                        for (j=0; j<4; ++j){
                                                std::cout << INFLOW[i][j][m][p] << "\t";
                                        }
                                        std::cout << std::endl;
                                }
                                std::cout << std::endl;
#endif
                        }
                }


#ifdef DEBUG
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[0][0] << "\t" << LAMBDA_PLUS[0][1] << "\t" << LAMBDA_PLUS[0][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[1][0] << "\t" << LAMBDA_PLUS[1][1] << "\t" << LAMBDA_PLUS[1][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[2][0] << "\t" << LAMBDA_PLUS[2][1] << "\t" << LAMBDA_PLUS[2][2] << std::endl;
                std::cout << "Lambda + =\t" << LAMBDA_PLUS[3][0] << "\t" << LAMBDA_PLUS[3][1] << "\t" << LAMBDA_PLUS[3][2] << std::endl;

                std::cout << "Lambda - =\t" << LAMBDA_MINUS[0][0] << "\t" << LAMBDA_MINUS[0][1] << "\t" << LAMBDA_MINUS[0][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[1][0] << "\t" << LAMBDA_MINUS[1][1] << "\t" << LAMBDA_MINUS[1][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[2][0] << "\t" << LAMBDA_MINUS[2][1] << "\t" << LAMBDA_MINUS[2][2] << std::endl;
                std::cout << "Lambda - =\t" << LAMBDA_MINUS[3][0] << "\t" << LAMBDA_MINUS[3][1] << "\t" << LAMBDA_MINUS[3][2] << std::endl;

                std::cout << "Lambda   =\t" << LAMBDA[0][0] << "\t" << LAMBDA[0][1] << "\t" << LAMBDA[0][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[1][0] << "\t" << LAMBDA[1][1] << "\t" << LAMBDA[1][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[2][0] << "\t" << LAMBDA[2][1] << "\t" << LAMBDA[2][2] << std::endl;
                std::cout << "Lambda   =\t" << LAMBDA[3][0] << "\t" << LAMBDA[3][1] << "\t" << LAMBDA[3][2] << std::endl;
#endif

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                // std::cout << SECOND_FLUC_N[0][0] << "\t" << SECOND_FLUC_N[0][1] << "\t" << SECOND_FLUC_N[0][2] << std::endl;

#if defined(LDA_SCHEME) or defined(BLENDED) or defined(Bx_SCHEME)

                double PHI_HALF[4];
                double SECOND_FLUC_LDA[4][3];

                for(i=0;i<4;++i){
                        PHI_HALF[i] = 0.0;
                        for(m=0;m<3;++m){
                                PHI_HALF[i] += INFLOW[i][0][m][2]*W_HAT[0][m] + INFLOW[i][1][m][2]*W_HAT[1][m] + INFLOW[i][2][m][2]*W_HAT[2][m] + INFLOW[i][3][m][2]*W_HAT[3][m];
                        }
                }

                // Calculate spatial splitting for first half timestep

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                FLUC_HALF_LDA[i][m] = BETA[i][0][m] * (PHI_HALF[0]) + BETA[i][1][m] * (PHI_HALF[1]) + BETA[i][2][m] * (PHI_HALF[2]) + BETA[i][3][m] * (PHI_HALF[3]);
                        }
                }

#ifdef DEBUG
                for(m=0;m<3;++m){
                        for(i=0;i<4;++i){
                                std::cout << "BETA =\t" << m << "\t" <<  BETA[i][0][m] << "\t" << BETA[i][1][m] << "\t" << BETA[i][2][m] << "\t" << BETA[i][3][m] << std::endl;
                        }
                }
#endif
                double DIFF[4][3];
                double MASS[4][4][3];
                double MASS_DIFF[4][3];
                double SUM_MASS[4];

                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                for(m=0;m<3;++m){MASS[i][j][m] = AREA * BETA[i][j][m]/3.0;}
                        }
                }

#ifdef DEBUG
                for(m=0;m<3;++m){
                        for(i=0;i<4;++i){
                                std::cout << "MASS =\t" << i << "\t" <<  MASS[i][0][m] << "\t" << MASS[i][1][m] << "\t" << MASS[i][2][m] << "\t" << MASS[i][3][m] << std::endl;
                        }
                }
#endif

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                DIFF[i][m] = U_HALF[i][m] - U_N[i][m];
                        }
                }

#ifdef DEBUG
                for(i=0;i<4;++i){
                        std::cout << "DIFF =\t" << i << "\t" <<  DIFF[i][0] << "\t" << DIFF[i][1] << "\t" << DIFF[i][2] << std::endl;
                }
#endif

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                MASS_DIFF[i][m] = MASS[i][0][m] * DIFF[0][m] + MASS[i][1][m] * DIFF[1][m] + MASS[i][2][m] * DIFF[2][m] + MASS[i][3][m] * DIFF[3][m];
                        }
                }

#ifdef DEBUG
                for(i=0;i<4;++i){
                        std::cout << "MASS_DIFF =\t" << i << "\t" <<  MASS_DIFF[i][0] << "\t" << MASS_DIFF[i][1] << "\t" << MASS_DIFF[i][2] << std::endl;
                }
#endif

                for(i=0;i<4;i++){
                        SUM_MASS[i] = 0.0;
                        for(m=0;m<3;++m){
                                if(DT==0.0){
                                        SUM_MASS[i] = 0.0;
                                }else{
                                        SUM_MASS[i] += MASS_DIFF[i][m]/DT;
                                }
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                SECOND_FLUC_LDA[i][m] = SUM_MASS[i] + 0.5*(FLUC_LDA[i][m] + FLUC_HALF_LDA[i][m]);
                        }
                }
                // std::cout << "2nd half fluctiation (LDA) =\t" << SECOND_FLUC_LDA[0][0] << "\t" << SECOND_FLUC_LDA[0][1] << "\t" << SECOND_FLUC_LDA[0][2] << std::endl;

#endif

#if defined(N_SCHEME) or defined(BLENDED) or defined(Bx_SCHEME)

                double INFLOW_MINUS_SUM[4][4];
                double SECOND_FLUC_N[4][3];

                for(i=0;i<4;++i){
                        for(j=0;j<4;++j){
                                INFLOW_MINUS_SUM[i][j] = 0.0;
                                for(m=0;m<3;++m){
                                        INFLOW_MINUS_SUM[i][j] += INFLOW[i][j][m][1];
                                }
                        }
                }

                mat_inv(&INFLOW_MINUS_SUM[0][0],4,X[0],Y[0],ID,2);

                double AREA_DIFF[4][3];
                double BRACKET[4][3];
                double KZ_SUM[4];

                for(i=0;i<4;++i){
                        KZ_SUM[i] = 0.0;
                        for(m=0;m<3;++m){
                                KZ_SUM[i] += INFLOW[i][0][m][1] * W_HAT[0][m] + INFLOW[i][1][m][1] * W_HAT[1][m] + INFLOW[i][2][m][1] * W_HAT[2][m] + INFLOW[i][3][m][1] * W_HAT[3][m];
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                BRACKET[i][m] = W_HAT[i][m] - (INFLOW_MINUS_SUM[i][0]*KZ_SUM[0] + INFLOW_MINUS_SUM[i][1]*KZ_SUM[1] + INFLOW_MINUS_SUM[i][2]*KZ_SUM[2] + INFLOW_MINUS_SUM[i][3]*KZ_SUM[3]);
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                FLUC_HALF_N[i][m] = INFLOW[i][0][m][0]*BRACKET[0][m] + INFLOW[i][1][m][0]*BRACKET[1][m] + INFLOW[i][2][m][0]*BRACKET[2][m] + INFLOW[i][3][m][0]*BRACKET[3][m];
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                AREA_DIFF[i][m] = AREA*(U_HALF[i][m] - U_N[i][m])/3.0;
                        }
                }

                for(i=0;i<4;++i){
                        for(m=0;m<3;++m){
                                if(DT == 0.0){
                                        SECOND_FLUC_N[i][m] = 0.0;
                                }else{
                                        SECOND_FLUC_N[i][m] = AREA_DIFF[i][m]/DT + 0.5*(FLUC_N[i][m] + FLUC_HALF_N[i][m]);
                                }
                        }
                }

#endif
                DUAL[0] = VERTEX_0->get_dual();
                DUAL[1] = VERTEX_1->get_dual();
                DUAL[2] = VERTEX_2->get_dual();

                // std::cout << SECOND_FLUC_LDA[0][0] << "\t" << SECOND_FLUC_LDA[0][1] << "\t" << SECOND_FLUC_LDA[0][2] << std::endl;

#ifdef LDA_SCHEME
                for(i=0;i<4;i++){
                        DU0[i] = -1.0*(DT/DUAL[0])*SECOND_FLUC_LDA[i][0];
                        DU1[i] = -1.0*(DT/DUAL[1])*SECOND_FLUC_LDA[i][1];
                        DU2[i] = -1.0*(DT/DUAL[2])*SECOND_FLUC_LDA[i][2];
                }
#endif

                // std::cout << SECOND_FLUC_N[0][0] << "\t" << SECOND_FLUC_N[0][1] << "\t" << SECOND_FLUC_N[0][2] << std::endl;

#ifdef N_SCHEME
                for(i=0;i<4;i++){
                        DU0[i] = -1.0*(DT/DUAL[0])*SECOND_FLUC_N[i][0];
                        DU1[i] = -1.0*(DT/DUAL[1])*SECOND_FLUC_N[i][1];
                        DU2[i] = -1.0*(DT/DUAL[2])*SECOND_FLUC_N[i][2];
                }
#endif

#ifdef BLENDED
                double THETA_E[4][4];
                double IDENTITY[4][4];
                double SUM_FLUC_N[4];

                THETA_E[0][0] = IDENTITY[0][0] = 1.0;
                THETA_E[0][1] = IDENTITY[0][1] = 0.0;
                THETA_E[0][2] = IDENTITY[0][2] = 0.0;
                THETA_E[0][3] = IDENTITY[0][3] = 0.0;

                THETA_E[1][0] = IDENTITY[1][0] = 0.0;
                THETA_E[1][1] = IDENTITY[1][1] = 1.0;
                THETA_E[1][2] = IDENTITY[1][2] = 0.0;
                THETA_E[1][3] = IDENTITY[1][3] = 0.0;

                THETA_E[2][0] = IDENTITY[2][0] = 0.0;
                THETA_E[2][1] = IDENTITY[2][1] = 0.0;
                THETA_E[2][2] = IDENTITY[2][2] = 1.0;
                THETA_E[2][3] = IDENTITY[2][3] = 0.0;

                THETA_E[3][0] = IDENTITY[3][0] = 0.0;
                THETA_E[3][1] = IDENTITY[3][1] = 0.0;
                THETA_E[3][2] = IDENTITY[3][2] = 0.0;
                THETA_E[3][3] = IDENTITY[3][3] = 1.0;

                for(i=0;i<4;i++){
                        SUM_FLUC_N[i] = abs(SECOND_FLUC_N[i][0]) + abs(SECOND_FLUC_N[i][1]) + abs(SECOND_FLUC_N[i][2]);

                        THETA_E[i][i] = abs(PHI[i])/SUM_FLUC_N[i];

                        FLUC_B[i][0] = THETA_E[i][i]*SECOND_FLUC_N[i][0] + (IDENTITY[i][i] - THETA_E[i][i])*SECOND_FLUC_LDA[i][0];
                        FLUC_B[i][1] = THETA_E[i][i]*SECOND_FLUC_N[i][1] + (IDENTITY[i][i] - THETA_E[i][i])*SECOND_FLUC_LDA[i][1];
                        FLUC_B[i][2] = THETA_E[i][i]*SECOND_FLUC_N[i][2] + (IDENTITY[i][i] - THETA_E[i][i])*SECOND_FLUC_LDA[i][2];

                        DU0[i] = -1.0*DT*FLUC_B[i][0]/DUAL[0];
                        DU1[i] = -1.0*DT*FLUC_B[i][1]/DUAL[1];
                        DU2[i] = -1.0*DT*FLUC_B[i][2]/DUAL[2];
                }
#endif


#ifdef Bx_SCHEME              
             
     		double  Xc, Yc;
     		double interpol[3][5];
     		double interpol_new[3][5];
     		double grad_p[2];
     		
     		
     		Xc=std::accumulate(std::begin(X), std::end(X), 0);
     		Yc=std::accumulate(std::begin(Y), std::end(Y), 0);
     		Xc/=3;
     		Yc/=3;
     		
     		cell_interpolation_N(interpol);
     		cell_interpolation_HALF(interpol_new);
     		
     		
             	grad_p[0]=interpol_new[0][4];
             	grad_p[1]=interpol_new[1][4];
             	
             	double Vxc, Vyc, Vmag;
             	
             	Vxc=interpol_new[0][1]*Xc+interpol_new[1][1]*Yc+interpol_new[2][1];
           	Vyc=interpol_new[0][2]*Xc+interpol_new[1][2]*Yc+interpol_new[2][2];
             	Vmag=sqrt(Vxc*Vxc+Vyc*Vyc);
             	double P_old=interpol[0][4]*Xc+interpol[1][4]*Yc+interpol[2][4];
             	double P_half=interpol_new[0][4]*Xc+interpol_new[1][4]*Yc+interpol_new[2][4];
             	
             	double sc;
             	double h=sqrt(4*AREA/3.1415);
             	
             	if(Vmag==0 ||P_DIFF==0){sc=0;}
             	else{sc=(h/(Vmag*P_DIFF))*(((P_half-P_old)/DT)+(grad_p[0]*Vxc)+(grad_p[1]*Vyc));}
             	//printf("%f, %f,%f ",sc,P_half,P_old);
             	if (sc<0){sc=0;}
             	if (sc>1){sc=1;}
             	
             	
             	
             	double THETA=sc*sc*h;
             	
             	
             	
             	for(i=0;i<4;i++){
                        FLUC_Bx[i][0] = THETA*SECOND_FLUC_N[i][0] + (1 - THETA)*SECOND_FLUC_LDA[i][0];
                        FLUC_Bx[i][1] = THETA*SECOND_FLUC_N[i][1] + (1 - THETA)*SECOND_FLUC_LDA[i][1];
                        FLUC_Bx[i][2] = THETA*SECOND_FLUC_N[i][2] + (1 - THETA)*SECOND_FLUC_LDA[i][2];
                        
                        DU0[i] = -1.0*DT*FLUC_Bx[i][0]/DUAL[0];
                        DU1[i] = -1.0*DT*FLUC_Bx[i][1]/DUAL[1];
                        DU2[i] = -1.0*DT*FLUC_Bx[i][2]/DUAL[2];
                }
             	
#endif

        }

        void pass_update(){
                VERTEX_0->update_du(DU0);
                VERTEX_1->update_du(DU1);
                VERTEX_2->update_du(DU2);
        }

        // Returns Roe average of left and right states
        double roe_avg(double L1, double L2, double R1, double R2){
                double AVG;
                AVG = (sqrt(L1)*L2+sqrt(R1)*R2)/(sqrt(L1)+sqrt(R1));
                return AVG;
        }




	void cell_interpolation_N(double (&coeff)[3][5]){
	
		double X1,X2,X3,Y1,Y2,Y3;
		double U[3];
		double A_inv[3][3];
		double det_A;
		//double coeff[3][5];
	
	
	
		X1 = X[0];
                X2 = X[1];
                X3 = X[2];

                Y1 = Y[0];
                Y2 = Y[1];
                Y3 = Y[2];
                
                det_A=X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2);
                
		A_inv[0][0]=(Y2-Y3);
		A_inv[0][1]=-(X2-X3);
		A_inv[0][2]=(X2*Y3-X3*Y2);
		
		A_inv[1][0]=-(Y3-Y1);
		A_inv[1][1]=(X1-X3);
		A_inv[1][2]=-(X1*Y3-X3*Y1);
		
		A_inv[0][0]=(Y1-Y2);
		A_inv[0][1]=-(X1-X2);
		A_inv[0][2]=(X1*Y2-X2*Y1);
		            
		for (int i=0;i<5;i++){
		
			if (i==0){
				for (int j=0;j<3;j++){
                			U[j]=U_N[i][j];
                		}
			}
			else if (i<=3){
                		for (int j=0;j<3;j++){
                			U[j]=U_N[i][j]/U_N[0][j];
                		}
                	}
                
                	else if (i==4){
                		for (int j=0;j<3;j++){
                			U[j]=PRESSURE[j];
                		}
                	}
                
		
		
			for (int j=0;j<3;j++){
				coeff[j][i]=(1/det_A)*(A_inv[j][0]*U[0]+A_inv[j][1]*U[1]+A_inv[j][2]*U[2]);
			
                	}
		
		}	
                

        	//return coeff;        
	}
	
	
	
	void cell_interpolation_HALF(double (&coeff)[3][5]){
	
		double X1,X2,X3,Y1,Y2,Y3;
		double U[3];
		double A_inv[3][3];
		double det_A;

	
	
		X1 = X[0];
                X2 = X[1];
                X3 = X[2];

                Y1 = Y[0];
                Y2 = Y[1];
                Y3 = Y[2];
                
                det_A=X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2);
                
		A_inv[0][0]=(Y2-Y3);
		A_inv[0][1]=-(X2-X3);
		A_inv[0][2]=(X2*Y3-X3*Y2);
		
		A_inv[1][0]=-(Y3-Y1);
		A_inv[1][1]=(X1-X3);
		A_inv[1][2]=-(X1*Y3-X3*Y1);
		
		A_inv[0][0]=(Y1-Y2);
		A_inv[0][1]=-(X1-X2);
		A_inv[0][2]=(X1*Y2-X2*Y1);
		            
		for (int i=0;i<5;i++){
		
			if (i==0){
				for (int j=0;j<3;j++){
                			U[j]=U_HALF[i][j];
                		}
			}
			else if (i<=3){
                		for (int j=0;j<3;j++){
                			U[j]=U_HALF[i][j]/U_HALF[0][j];
                		}
                	}
                
                	else if (i==4){
                		for (int j=0;j<3;j++){
                			U[j]=PRESSURE_HALF[j];
                		}
                	}
                
		
		
			for (int j=0;j<3;j++){
				coeff[j][i]=(1/det_A)*(A_inv[j][0]*U[0]+A_inv[j][1]*U[1]+A_inv[j][2]*U[2]);
			
                	}
		
		}	       
	}
	
	
        void setup_normals(){
                // Calculate normals (just in first timestep for static grid)

                setup_positions();

                for(int m=0; m<3; ++m){X_MOD[m] = X[m];Y_MOD[m] = Y[m];}

#ifdef PERIODIC_BOUNDARY
                if(BOUNDARY == 1){
                        for(int i=0; i<3; ++i){
                                for(int j=0; j<3; ++j){
                                        if(X[j] - X[i] > 0.5*SIDE_LENGTH_X){
                                                X_MOD[i] = X[i] + SIDE_LENGTH_X;
                                        }
                                        if(Y[j] - Y[i] > 0.5*SIDE_LENGTH_Y){
                                                Y_MOD[i] = Y[i] + SIDE_LENGTH_Y;
                                        }
                                }
                        }
                }
#endif

                // check vertices are ordered counter-clockwise

                double X0,X1,X2,Y0,Y1,Y2;
                double L1X,L1Y,L2X,L2Y,CROSS;

                X0 = X_MOD[0];
                X1 = X_MOD[1];
                X2 = X_MOD[2];

                Y0 = Y_MOD[0];
                Y1 = Y_MOD[1];
                Y2 = Y_MOD[2];

                L1X = X1 - X0;
                L1Y = Y1 - Y0;

                L2X = X2 - X0;
                L2Y = Y2 - Y0;

                CROSS = L1X*L2Y - L1Y*L2X;

                if(CROSS < 0.0){
                        reorder_vertices();
                }
                
                calculate_normals(X_MOD,Y_MOD);


        }

        void calculate_normals(double X_input[3],double Y_input[3]){
                int i;
                double PERP[3][2];

                PERP[0][0] = (Y_input[1] - Y_input[2]);
                PERP[0][1] = (X_input[2] - X_input[1]);

                PERP[1][0] = (Y_input[2] - Y_input[0]);
                PERP[1][1] = (X_input[0] - X_input[2]);

                PERP[2][0] = (Y_input[0] - Y_input[1]);
                PERP[2][1] = (X_input[1] - X_input[0]);

                // calculate area of triangle and pass 1/3 to each vertex for dual

                double THETA0 = atan(PERP[0][1]/PERP[0][0]);
                double THETA1 = atan(PERP[1][1]/PERP[1][0]);

                double THETA = std::abs(THETA0 - THETA1);

                if(THETA > 3.14159/2.0){THETA = 3.14159 - THETA;}

                AREA = 0.5*(sqrt(PERP[0][0]*PERP[0][0] + PERP[0][1]*PERP[0][1])*sqrt(PERP[1][0]*PERP[1][0] + PERP[1][1]*PERP[1][1]))*sin(THETA);

                VERTEX_0->calculate_dual(AREA/3.0);
                VERTEX_1->calculate_dual(AREA/3.0);
                VERTEX_2->calculate_dual(AREA/3.0);

                for(i=0;i<3;i++){
                        MAG[i] = sqrt(PERP[i][0]*PERP[i][0]+PERP[i][1]*PERP[i][1]);
                        NORMAL[i][0] = PERP[i][0]/MAG[i];
                        NORMAL[i][1] = PERP[i][1]/MAG[i];
                }

        }

        void calculate_len_vel_contribution(){
                int m;
                double L02,L12,L22;
                double H,U,V,VEL[3];
                double C_SOUND[3];
                double LMAX,VMAX,CONT;

                setup_initial_state();

                L02 = (X_MOD[1] - X_MOD[0])*(X_MOD[1] - X_MOD[0]) + (Y_MOD[1] - Y_MOD[0])*(Y_MOD[1] - Y_MOD[0]);
                L12 = (X_MOD[2] - X_MOD[0])*(X_MOD[2] - X_MOD[0]) + (Y_MOD[2] - Y_MOD[0])*(Y_MOD[2] - Y_MOD[0]);
                L22 = (X_MOD[2] - X_MOD[1])*(X_MOD[2] - X_MOD[1]) + (Y_MOD[2] - Y_MOD[1])*(Y_MOD[2] - Y_MOD[1]);

                LMAX = max_val(L02,L12);
                LMAX = max_val(LMAX,L22);

                LMAX = sqrt(LMAX);

                // std::cout << LMAX << std::endl;

                for(m=0;m<3;++m){
                        H = (U_N[3][m] + PRESSURE[m])/U_N[0][m];                       
                        U = U_N[1][m]/U_N[0][m];
                        V = U_N[2][m]/U_N[0][m];
                        VEL[m] = sqrt(U*U + V*V);
                        C_SOUND[m] = sqrt((GAMMA-1.0) * H - (GAMMA-1.0) * (U*U + V*V)/2.0);
                }

                VMAX = max_val((VEL[0] + C_SOUND[0]),(VEL[1] + C_SOUND[1]));
                VMAX = max_val(VMAX,(VEL[2] + C_SOUND[2]));

                CONT = LMAX * VMAX;
               	/*printf("%f\n",VEL[0]);
               	printf("%f\n",VEL[1]);
               	printf("%f\n",VEL[2]);*/

                VERTEX_0->update_len_vel_sum(CONT);
                VERTEX_1->update_len_vel_sum(CONT);
                VERTEX_2->update_len_vel_sum(CONT);

        }

#ifdef DRIFT_SHELL
        void send_tbin_limit(){
                VERTEX_0->reset_tbin_local(2*TBIN);
                VERTEX_1->reset_tbin_local(2*TBIN);
                VERTEX_2->reset_tbin_local(2*TBIN);
        }

        void check_tbin(){
                int TBIN0,TBIN1,TBIN2;
                TBIN0 = VERTEX_0->get_tbin_local();
                TBIN1 = VERTEX_1->get_tbin_local();
                TBIN2 = VERTEX_2->get_tbin_local();
                if(TBIN0<TBIN){TBIN=TBIN0;}
                if(TBIN1<TBIN){TBIN=TBIN1;}
                if(TBIN2<TBIN){TBIN=TBIN2;}
        }
#endif

        void reorder_vertices(){

//              std::cout<<"reordering !!!!!!"<<std::endl;
//              Warning: X_MOD, Y_MOD haven't been reordered!
                VERTEX *TEMP_VERTEX;

                TEMP_VERTEX = VERTEX_1;
                VERTEX_1 = VERTEX_2;
                VERTEX_2 = TEMP_VERTEX;

        }     
       
};
