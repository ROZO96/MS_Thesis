/*
 * This file was written by Ben Morton (bmorton@ed.ac.uk).
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

#include "base.h"
#include "inverse.h"
#include "vertex2D.h"

using namespace std;
#pragma once

class TRIANGLE{

private:
        int ID;
        VERTEX *VERTEX_0,*VERTEX_1,*VERTEX_2,*VERTEX_3;

        int BOUNDARY;
        int TBIN;

        double VOLUME;

        double X[4],Y[4],Z[4],DUAL[4];
        double X_MOD[4],Y_MOD[4],Z_MOD[4];
        double NORMAL[4][3];

        double U_N[5][4];
        double U_HALF[5][4];

        double FLUC_LDA[5][4],FLUC_N[5][4],FLUC_B[5][4];
        double FLUC_HALF_LDA[5][4],FLUC_HALF_N[5][4],FLUC_HALF_B[5][4];

        double PRESSURE[4], PRESSURE_HALF[4];

        double PHI[5];
        double BETA[5][5][4];

        double DU0[5],DU1[5],DU2[5],DU3[5];
        double DU0_HALF[5],DU1_HALF[5],DU2_HALF[5],DU3_HALF[5];

        double MAG[4];

        int PRINT;

public:

        void set_id(int NEW_ID){ID = NEW_ID;}

        void set_vertex_0(VERTEX* NEW_VERTEX){VERTEX_0 = NEW_VERTEX;}
        void set_vertex_1(VERTEX* NEW_VERTEX){VERTEX_1 = NEW_VERTEX;}
        void set_vertex_2(VERTEX* NEW_VERTEX){VERTEX_2 = NEW_VERTEX;}
        void set_vertex_3(VERTEX* NEW_VERTEX){VERTEX_3 = NEW_VERTEX;}

        void set_boundary(int NEW_BOUNDARY){BOUNDARY = NEW_BOUNDARY;}
        void set_tbin(    int NEW_TBIN){TBIN = NEW_TBIN;}

        int get_id(){return ID;}

        VERTEX* get_vertex_0(){return VERTEX_0;}
        VERTEX* get_vertex_1(){return VERTEX_1;}
        VERTEX* get_vertex_2(){return VERTEX_2;}
        VERTEX* get_vertex_3(){return VERTEX_3;}

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
        double get_un03(){
                U_N[0][3] = VERTEX_3->get_u0();
                return U_N[0][3];
        }

        void print_triangle_state(){
                std::cout << ID << "\tU[0] =\t" << U_N[0][0] << "\t" << U_N[0][1] << "\t" << U_N[0][2] << "\t" << U_N[0][3] << std::endl;
                std::cout << ID << "\tU[1] =\t" << U_N[1][0] << "\t" << U_N[1][1] << "\t" << U_N[1][2] << "\t" << U_N[1][3] << std::endl;
                std::cout << ID << "\tU[2] =\t" << U_N[2][0] << "\t" << U_N[2][1] << "\t" << U_N[2][2] << "\t" << U_N[2][3]<< std::endl;
                std::cout << ID << "\tU[3] =\t" << U_N[3][0] << "\t" << U_N[3][1] << "\t" << U_N[3][2] << "\t" << U_N[3][3]<< std::endl;
                std::cout << ID << "\tU[4] =\t" << U_N[4][0] << "\t" << U_N[4][1] << "\t" << U_N[4][2] << "\t" << U_N[4][3]<< std::endl;

                std::cout << ID << "\tU_HALF[0] =\t" << U_HALF[0][0] << "\t" << U_HALF[0][1] << "\t" << U_HALF[0][2] << "\t" << U_HALF[0][3]<< std::endl;
                std::cout << ID << "\tU_HALF[1] =\t" << U_HALF[1][0] << "\t" << U_HALF[1][1] << "\t" << U_HALF[1][2] << "\t" << U_HALF[1][3]<< std::endl;
                std::cout << ID << "\tU_HALF[2] =\t" << U_HALF[2][0] << "\t" << U_HALF[2][1] << "\t" << U_HALF[2][2] << "\t" << U_HALF[2][3]<< std::endl;
                std::cout << ID << "\tU_HALF[3] =\t" << U_HALF[3][0] << "\t" << U_HALF[3][1] << "\t" << U_HALF[3][2] << "\t" << U_HALF[3][3]<< std::endl;
                std::cout << ID << "\tU_HALF[4] =\t" << U_HALF[4][0] << "\t" << U_HALF[4][1] << "\t" << U_HALF[4][2] << "\t" << U_HALF[4][3]<< std::endl;
        }

        // import x and y for all vertices
        void setup_positions(){
                X[0] = VERTEX_0->get_x();
                X[1] = VERTEX_1->get_x();
                X[2] = VERTEX_2->get_x();
                X[3] = VERTEX_3->get_x();

                Y[0] = VERTEX_0->get_y();
                Y[1] = VERTEX_1->get_y();
                Y[2] = VERTEX_2->get_y();
                Y[3] = VERTEX_3->get_y();

                Z[0] = VERTEX_0->get_z();
                Z[1] = VERTEX_1->get_z();
                Z[2] = VERTEX_2->get_z();
                Z[3] = VERTEX_3->get_z();
        }

        void setup_dual(){
                DUAL[0] = VERTEX_0->get_dual();
                DUAL[1] = VERTEX_1->get_dual();
                DUAL[2] = VERTEX_2->get_dual();
                DUAL[3] = VERTEX_3->get_dual();
        }

        // import initial fluid state and pressure for all vertices
        void setup_initial_state(){
                U_N[0][0] = VERTEX_0->get_u0();
                U_N[0][1] = VERTEX_1->get_u0();
                U_N[0][2] = VERTEX_2->get_u0();
                U_N[0][3] = VERTEX_3->get_u0();

                U_N[1][0] = VERTEX_0->get_u1();
                U_N[1][1] = VERTEX_1->get_u1();
                U_N[1][2] = VERTEX_2->get_u1();
                U_N[1][3] = VERTEX_3->get_u1();

                U_N[2][0] = VERTEX_0->get_u2();
                U_N[2][1] = VERTEX_1->get_u2();
                U_N[2][2] = VERTEX_2->get_u2();
                U_N[2][3] = VERTEX_3->get_u2();

                U_N[3][0] = VERTEX_0->get_u3();
                U_N[3][1] = VERTEX_1->get_u3();
                U_N[3][2] = VERTEX_2->get_u3();
                U_N[3][3] = VERTEX_3->get_u3();

                U_N[4][0] = VERTEX_0->get_u4();
                U_N[4][1] = VERTEX_1->get_u4();
                U_N[4][2] = VERTEX_2->get_u4();
                U_N[4][3] = VERTEX_3->get_u4();

                PRESSURE[0] = VERTEX_0->get_pressure();
                PRESSURE[1] = VERTEX_1->get_pressure();
                PRESSURE[2] = VERTEX_2->get_pressure();
                PRESSURE[3] = VERTEX_3->get_pressure();
        }

        // import intermediate fluid state and pressure for all vertices
        void setup_half_state(){
                U_HALF[0][0] = VERTEX_0->get_u0_half();
                U_HALF[0][1] = VERTEX_1->get_u0_half();
                U_HALF[0][2] = VERTEX_2->get_u0_half();
                U_HALF[0][3] = VERTEX_3->get_u0_half();

                U_HALF[1][0] = VERTEX_0->get_u1_half();
                U_HALF[1][1] = VERTEX_1->get_u1_half();
                U_HALF[1][2] = VERTEX_2->get_u1_half();
                U_HALF[1][3] = VERTEX_3->get_u1_half();

                U_HALF[2][0] = VERTEX_0->get_u2_half();
                U_HALF[2][1] = VERTEX_1->get_u2_half();
                U_HALF[2][2] = VERTEX_2->get_u2_half();
                U_HALF[2][3] = VERTEX_3->get_u2_half();

                U_HALF[3][0] = VERTEX_0->get_u3_half();
                U_HALF[3][1] = VERTEX_1->get_u3_half();
                U_HALF[3][2] = VERTEX_2->get_u3_half();
                U_HALF[3][3] = VERTEX_3->get_u3_half();

                U_HALF[4][0] = VERTEX_0->get_u4_half();
                U_HALF[4][1] = VERTEX_1->get_u4_half();
                U_HALF[4][2] = VERTEX_2->get_u4_half();
                U_HALF[4][3] = VERTEX_3->get_u4_half();

                PRESSURE_HALF[0] = VERTEX_0->get_pressure_half();
                PRESSURE_HALF[1] = VERTEX_1->get_pressure_half();
                PRESSURE_HALF[2] = VERTEX_2->get_pressure_half();
                PRESSURE_HALF[3] = VERTEX_3->get_pressure_half();
        }

        //**********************************************************************************************************************

        // Calculate first half timestep change, passing change to vertex
        void calculate_first_half(double T, double DT){
                int i,j,m,p;
                double INFLOW[5][5][4][3];        // K+, K-, K matrices for each vertex (m index for vertices, p index for +,-,0)

                // Import conditions and positions of vertices

                setup_positions();
                setup_dual();
                setup_initial_state();

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Calculate inflow parameters

                double H[4];
                double RHO,C,VX,VX_C,VY,VY_C,VZ,VZ_C,H_AVG,H_C,ALPHA,ALPHA_C,W,NORM=1.0/3.0;
                double Z_ROE[5][4];
                double VALUE1,VALUE2,VALUE3,VALUE12,VALUE123;
                double PRESSURE_AVG;
                double LAMBDA[5][4],LAMBDA_PLUS[5][4],LAMBDA_MINUS[5][4];
                double N_X[4],N_Y[4],N_Z[4];

                double Z_BAR[5],W_HAT[5][4];

                // Construct Roe vector Z

                for(m=0;m<4;++m){
                        Z_ROE[0][m] = sqrt(U_N[0][m]);
                        Z_ROE[1][m] = U_N[1][m]/Z_ROE[0][m];
                        Z_ROE[2][m] = U_N[2][m]/Z_ROE[0][m];
                        Z_ROE[3][m] = U_N[3][m]/Z_ROE[0][m];
                        Z_ROE[4][m] = (U_N[4][m] + PRESSURE[m])/Z_ROE[0][m];

                        N_X[m]  = NORMAL[m][0];
                        N_Y[m]  = NORMAL[m][1];
                        N_Z[m]  = NORMAL[m][2];

                        H[m] = (U_N[4][m] + PRESSURE[m])/U_N[0][m];
                }

                for(i=0; i<5; ++i){Z_BAR[i] = (Z_ROE[i][0] + Z_ROE[i][1] + Z_ROE[i][2] + Z_ROE[i][3])/4.0;}

                for(m=0; m<4; ++m){
                        W_HAT[0][m] =  2.0*Z_BAR[0]*Z_ROE[0][m];
                        W_HAT[1][m] =  Z_BAR[1]*Z_ROE[0][m] + Z_BAR[0]*Z_ROE[1][m];
                        W_HAT[2][m] =  Z_BAR[2]*Z_ROE[0][m] + Z_BAR[0]*Z_ROE[2][m];
                        W_HAT[3][m] =  Z_BAR[3]*Z_ROE[0][m] + Z_BAR[0]*Z_ROE[3][m];
                        W_HAT[4][m] = (Z_BAR[4]*Z_ROE[0][m] + GAMMA_1*Z_BAR[1]*Z_ROE[1][m] + GAMMA_1*Z_BAR[2]*Z_ROE[2][m] + GAMMA_1*Z_BAR[3]*Z_ROE[3][m] + Z_BAR[0]*Z_ROE[4][m])/GAMMA;
                }

                // Construct average state for element

                RHO   = pow((sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]) + sqrt(U_N[0][3]))/4.0, 2);
                VX    = (sqrt(U_N[0][0])*U_N[1][0]/U_N[0][0] + sqrt(U_N[0][1])*U_N[1][1]/U_N[0][1] + sqrt(U_N[0][2])*U_N[1][2]/U_N[0][2] + sqrt(U_N[0][3])*U_N[1][3]/U_N[0][3]) / (sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]) + sqrt(U_N[0][3]));
                VY    = (sqrt(U_N[0][0])*U_N[2][0]/U_N[0][0] + sqrt(U_N[0][1])*U_N[2][1]/U_N[0][1] + sqrt(U_N[0][2])*U_N[2][2]/U_N[0][2] + sqrt(U_N[0][3])*U_N[2][3]/U_N[0][3]) / (sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]) + sqrt(U_N[0][3]));
                VZ    = (sqrt(U_N[0][0])*U_N[3][0]/U_N[0][0] + sqrt(U_N[0][1])*U_N[3][1]/U_N[0][1] + sqrt(U_N[0][2])*U_N[3][2]/U_N[0][2] + sqrt(U_N[0][3])*U_N[3][3]/U_N[0][3]) / (sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]) + sqrt(U_N[0][3]));
                H_AVG = (sqrt(U_N[0][0])*H[0] + sqrt(U_N[0][1])*H[1] + sqrt(U_N[0][2])*H[2] + sqrt(U_N[0][3])*H[3]) / (sqrt(U_N[0][0]) + sqrt(U_N[0][1]) + sqrt(U_N[0][2]) + sqrt(U_N[0][3]));

                PRESSURE_AVG = (PRESSURE[0] + PRESSURE[1] + PRESSURE[2] + PRESSURE[3])/4.0;
                C = sqrt((GAMMA-1.0) * H_AVG - (GAMMA-1.0) * (VX*VX + VY*VY + VZ*VZ)/2.0);
                // if(C <= C_LIM){C = C_LIM;}

                VX_C = VX/C;
                VY_C = VY/C;
                VZ_C = VZ/C;

                H_C = H_AVG/C;

                ALPHA   = GAMMA_1*(VX*VX + VY*VY + VZ*VZ)/2.0;
                ALPHA_C = ALPHA/C;

                // Calculate K+,K- and K matrices for each vertex i,j,k

                for(m=0;m<4;++m){

                        W = VX*N_X[m] + VY*N_Y[m] + VZ*N_Z[m];

                        LAMBDA[0][m] = W + C;
                        LAMBDA[1][m] = W - C;
                        LAMBDA[2][m] = W;
                        LAMBDA[3][m] = W;
                        LAMBDA[4][m] = W;

                        for(i=0;i<5;++i){
                                LAMBDA_PLUS[i][m]  = max_val(0.0,LAMBDA[i][m]);
                                LAMBDA_MINUS[i][m] = min_val(0.0,LAMBDA[i][m]);
                        }

                        for(p=0;p<3;++p){
                                if(p==0){       // Identify and select positive eigenvalues
                                        VALUE1 = LAMBDA_PLUS[0][m];
                                        VALUE2 = LAMBDA_PLUS[1][m];
                                        VALUE3 = LAMBDA_PLUS[2][m];
#ifdef DEBUG
                                        std::cout << "K+" << std::endl;
#endif
                                }else if(p==1){ // Identify and select negative eigenvalues
                                        VALUE1 = LAMBDA_MINUS[0][m];
                                        VALUE2 = LAMBDA_MINUS[1][m];
                                        VALUE3 = LAMBDA_MINUS[2][m];
#ifdef DEBUG
                                        std::cout << "K-" << std::endl;
#endif
                                }else{          // Select all eigenvalues
                                        VALUE1 = LAMBDA[0][m];
                                        VALUE2 = LAMBDA[1][m];
                                        VALUE3 = LAMBDA[2][m];
#ifdef DEBUG
                                        std::cout << "K" << std::endl;
#endif
                                }

                                VALUE12  = (VALUE1 - VALUE2)/2.0;
                                VALUE123 = (VALUE1 + VALUE2 - 2.0*VALUE3)/2.0;

                                INFLOW[0][0][m][p] = NORM * MAG[m]*(ALPHA_C*VALUE123/C - W*VALUE12/C + VALUE3);
                                INFLOW[0][1][m][p] = NORM * MAG[m]*(-1.0*GAMMA_1*VX_C*VALUE123/C + N_X[m]*VALUE12/C);
                                INFLOW[0][2][m][p] = NORM * MAG[m]*(-1.0*GAMMA_1*VY_C*VALUE123/C + N_Y[m]*VALUE12/C);
                                INFLOW[0][3][m][p] = NORM * MAG[m]*(-1.0*GAMMA_1*VZ_C*VALUE123/C + N_Z[m]*VALUE12/C);
                                INFLOW[0][4][m][p] = NORM * MAG[m]*(GAMMA_1*VALUE123/(C*C));

                                INFLOW[1][0][m][p] = NORM * MAG[m]*((ALPHA_C*VX_C - W*N_X[m])*VALUE123 + (ALPHA_C*N_X[m] - VX_C*W)*VALUE12);
                                INFLOW[1][1][m][p] = NORM * MAG[m]*((N_X[m]*N_X[m] - GAMMA_1*VX_C*VX_C)*VALUE123 - (GAMMA_2*VX_C*N_X[m]*VALUE12) + VALUE3);
                                INFLOW[1][2][m][p] = NORM * MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*VX_C*VY_C)*VALUE123 + (VX_C*N_Y[m] - GAMMA_1*VY_C*N_X[m])*VALUE12);
                                INFLOW[1][3][m][p] = NORM * MAG[m]*((N_X[m]*N_Z[m] - GAMMA_1*VX_C*VZ_C)*VALUE123 + (VX_C*N_Z[m] - GAMMA_1*VZ_C*N_X[m])*VALUE12);
                                INFLOW[1][4][m][p] = NORM * MAG[m]*(GAMMA_1*VX_C*VALUE123/C + GAMMA_1*N_X[m]*VALUE12/C);

                                INFLOW[2][0][m][p] = NORM * MAG[m]*((ALPHA_C*VY_C - W*N_Y[m])*VALUE123 + (ALPHA_C*N_Y[m] - VY_C*W)*VALUE12);
                                INFLOW[2][1][m][p] = NORM * MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*VX_C*VY_C)*VALUE123 + (VY_C*N_X[m] - GAMMA_1*VX_C*N_Y[m])*VALUE12);
                                INFLOW[2][2][m][p] = NORM * MAG[m]*((N_Y[m]*N_Y[m] - GAMMA_1*VY_C*VY_C)*VALUE123 - (GAMMA_2*VY_C*N_Y[m]*VALUE12) + VALUE3);
                                INFLOW[2][3][m][p] = NORM * MAG[m]*((N_Z[m]*N_Y[m] - GAMMA_1*VZ_C*VY_C)*VALUE123 + (VY_C*N_Z[m] - GAMMA_1*VZ_C*N_Y[m])*VALUE12);
                                INFLOW[2][4][m][p] = NORM * MAG[m]*(GAMMA_1*VY_C*VALUE123/C + GAMMA_1*N_Y[m]*VALUE12/C);

                                INFLOW[3][0][m][p] = NORM * MAG[m]*((ALPHA_C*VZ_C - W*N_Z[m])*VALUE123 + (ALPHA_C*N_Z[m] - VZ_C*W)*VALUE12);
                                INFLOW[3][1][m][p] = NORM * MAG[m]*((N_X[m]*N_Z[m] - GAMMA_1*VX_C*VZ_C)*VALUE123 + (VZ_C*N_X[m] - GAMMA_1*VX_C*N_Z[m])*VALUE12);
                                INFLOW[3][2][m][p] = NORM * MAG[m]*((N_Y[m]*N_Z[m] - GAMMA_1*VY_C*VZ_C)*VALUE123 + (VZ_C*N_Y[m] - GAMMA_1*VY_C*N_Z[m])*VALUE12);
                                INFLOW[3][3][m][p] = NORM * MAG[m]*((N_Z[m]*N_Z[m] - GAMMA_1*VZ_C*VZ_C)*VALUE123 - (GAMMA_2*VZ_C*N_Z[m]*VALUE12) + VALUE3);
                                INFLOW[3][4][m][p] = NORM * MAG[m]*(GAMMA_1*VZ_C*VALUE123/C + GAMMA_1*N_Z[m]*VALUE12/C);

                                INFLOW[4][0][m][p] = NORM * MAG[m]*((ALPHA_C*H_C - W*W)*VALUE123 + W*(ALPHA_C - H_C)*VALUE12);
                                INFLOW[4][1][m][p] = NORM * MAG[m]*((W*N_X[m] - VX - ALPHA_C*VX_C)*VALUE123 + (H_C*N_X[m] - GAMMA_1*VX_C*W)*VALUE12);
                                INFLOW[4][2][m][p] = NORM * MAG[m]*((W*N_Y[m] - VY - ALPHA_C*VY_C)*VALUE123 + (H_C*N_Y[m] - GAMMA_1*VY_C*W)*VALUE12);
                                INFLOW[4][3][m][p] = NORM * MAG[m]*((W*N_Z[m] - VZ - ALPHA_C*VZ_C)*VALUE123 + (H_C*N_Z[m] - GAMMA_1*VZ_C*W)*VALUE12);
                                INFLOW[4][4][m][p] = NORM * MAG[m]*(GAMMA_1*H_C*VALUE123/C + GAMMA_1*W*VALUE12/C + VALUE3);
                        }
                }


                for(i=0;i<5;++i){
                        PHI[i] = 0.0;
                        for(m=0;m<4;++m){
                                PHI[i] += INFLOW[i][0][m][2]*W_HAT[0][m] + INFLOW[i][1][m][2]*W_HAT[1][m] + INFLOW[i][2][m][2]*W_HAT[2][m] + INFLOW[i][3][m][2]*W_HAT[3][m] + INFLOW[i][4][m][2]*W_HAT[4][m];
                        }
                }

                double INFLOW_MINUS_SUM[5][5];

                // if(ID == 402){std::cout << INFLOW[0][0][0][1] << "\t" << INFLOW[0][0][1][1] << "\t" << INFLOW[0][0][2][1] << "\t" << INFLOW[0][0][3][1] << "\t" << std::endl;}

                for(i=0;i<5;++i){
                        for(j=0;j<5;++j){
                                INFLOW_MINUS_SUM[i][j] = 0.0;
                                for(m=0;m<4;++m){
                                        INFLOW_MINUS_SUM[i][j] += INFLOW[i][j][m][1];
                                }
                        }
                }

                mat_inv(&INFLOW_MINUS_SUM[0][0],5,X[0],Y[0],ID,1);

                // Calculate spatial splitting for first half timestep

#if defined(LDA_SCHEME) or defined(BLENDED)

                for(i=0;i<5;++i){
                        for(j=0;j<5;++j){
                                for(m=0;m<4;++m){
                                        BETA[i][j][m] = -1.0*(INFLOW[i][0][m][0]*INFLOW_MINUS_SUM[0][j] + INFLOW[i][1][m][0]*INFLOW_MINUS_SUM[1][j] + INFLOW[i][2][m][0]*INFLOW_MINUS_SUM[2][j] + INFLOW[i][3][m][0]*INFLOW_MINUS_SUM[3][j] + INFLOW[i][4][m][0]*INFLOW_MINUS_SUM[4][j]);
                                }
                        }
                }

                for(i=0;i<5;++i){
                        for(m=0;m<4;++m){
                                FLUC_LDA[i][m] = 0.5*(BETA[i][0][m]*PHI[0] + BETA[i][1][m]*PHI[1] + BETA[i][2][m]*PHI[2] + BETA[i][3][m]*PHI[3] + BETA[i][4][m]*PHI[4]);
                        }
                }
#endif


#if defined(N_SCHEME) or defined(BLENDED)
                double BRACKET[5][4];
                double KZ_SUM[5];

                for(i=0;i<5;++i){
                        KZ_SUM[i] = 0.0;
                        for(m=0;m<4;++m){
                                KZ_SUM[i] += INFLOW[i][0][m][1] * W_HAT[0][m] + INFLOW[i][1][m][1] * W_HAT[1][m] + INFLOW[i][2][m][1] * W_HAT[2][m] + INFLOW[i][3][m][1] * W_HAT[3][m] + INFLOW[i][4][m][1] * W_HAT[4][m];
                        }
                }

                for(i=0;i<5;++i){
                        for(m=0;m<4;++m){
                                BRACKET[i][m] = W_HAT[i][m] - (INFLOW_MINUS_SUM[i][0]*KZ_SUM[0] + INFLOW_MINUS_SUM[i][1]*KZ_SUM[1] + INFLOW_MINUS_SUM[i][2]*KZ_SUM[2] + INFLOW_MINUS_SUM[i][3]*KZ_SUM[3] + INFLOW_MINUS_SUM[i][4]*KZ_SUM[4]);
                        }
                }

                for(i=0;i<5;++i){
                        for(m=0;m<4;++m){
                                FLUC_N[i][m] = 0.5*(INFLOW[i][0][m][0]*BRACKET[0][m] + INFLOW[i][1][m][0]*BRACKET[1][m] + INFLOW[i][2][m][0]*BRACKET[2][m] + INFLOW[i][3][m][0]*BRACKET[3][m] + INFLOW[i][4][m][0]*BRACKET[4][m]);
                        }
                }
#endif

#ifdef BLENDED
                double THETA_E[5][5];
                double IDENTITY[5][5];
                double SUM_FLUC_N[5];

                THETA_E[0][0] = IDENTITY[0][0] = 1.0;
                THETA_E[0][1] = IDENTITY[0][1] = 0.0;
                THETA_E[0][2] = IDENTITY[0][2] = 0.0;
                THETA_E[0][3] = IDENTITY[0][3] = 0.0;
                THETA_E[0][4] = IDENTITY[0][4] = 0.0;

                THETA_E[1][0] = IDENTITY[1][0] = 0.0;
                THETA_E[1][1] = IDENTITY[1][1] = 1.0;
                THETA_E[1][2] = IDENTITY[1][2] = 0.0;
                THETA_E[1][3] = IDENTITY[1][3] = 0.0;
                THETA_E[1][4] = IDENTITY[1][4] = 0.0;

                THETA_E[2][0] = IDENTITY[2][0] = 0.0;
                THETA_E[2][1] = IDENTITY[2][1] = 0.0;
                THETA_E[2][2] = IDENTITY[2][2] = 1.0;
                THETA_E[2][3] = IDENTITY[2][3] = 0.0;
                THETA_E[2][4] = IDENTITY[2][4] = 0.0;

                THETA_E[3][0] = IDENTITY[3][0] = 0.0;
                THETA_E[3][1] = IDENTITY[3][1] = 0.0;
                THETA_E[3][2] = IDENTITY[3][2] = 0.0;
                THETA_E[3][3] = IDENTITY[3][3] = 1.0;
                THETA_E[3][4] = IDENTITY[3][4] = 0.0;

                THETA_E[4][0] = IDENTITY[4][0] = 0.0;
                THETA_E[4][1] = IDENTITY[4][1] = 0.0;
                THETA_E[4][2] = IDENTITY[4][2] = 0.0;
                THETA_E[4][3] = IDENTITY[4][3] = 0.0;
                THETA_E[4][4] = IDENTITY[4][4] = 1.0;

                for(i=0;i<5;i++){
                        SUM_FLUC_N[i] = abs(FLUC_N[i][0]) + abs(FLUC_N[i][1]) + abs(FLUC_N[i][2]) + abs(FLUC_N[i][3]);

                        THETA_E[i][i] = abs(PHI[i])/SUM_FLUC_N[i];

                        FLUC_B[i][0] = THETA_E[i][i]*FLUC_N[i][0] + (IDENTITY[i][i] - THETA_E[i][i])*FLUC_LDA[i][0];
                        FLUC_B[i][1] = THETA_E[i][i]*FLUC_N[i][1] + (IDENTITY[i][i] - THETA_E[i][i])*FLUC_LDA[i][1];
                        FLUC_B[i][2] = THETA_E[i][i]*FLUC_N[i][2] + (IDENTITY[i][i] - THETA_E[i][i])*FLUC_LDA[i][2];
                        FLUC_B[i][3] = THETA_E[i][i]*FLUC_N[i][3] + (IDENTITY[i][i] - THETA_E[i][i])*FLUC_LDA[i][3];
                }
#endif

#ifdef LDA_SCHEME
                for(i=0;i<5;i++){
                        DU0_HALF[i] = DT*FLUC_LDA[i][0] / DUAL[0];
                        DU1_HALF[i] = DT*FLUC_LDA[i][1] / DUAL[1];
                        DU2_HALF[i] = DT*FLUC_LDA[i][2] / DUAL[2];
                        DU3_HALF[i] = DT*FLUC_LDA[i][3] / DUAL[3];
                        // if(PRINT == 1){std::cout << ID << "\tDU i =\t" << i << "\t" << DU0[i] << "\t" << DU1[i] << "\t" << DU2[i] << "\t" << DU3[i] << std::endl;}
                }
#endif

#ifdef N_SCHEME
                for(i=0;i<5;i++){
                        DU0_HALF[i] = DT*FLUC_N[i][0]/DUAL[0];
                        DU1_HALF[i] = DT*FLUC_N[i][1]/DUAL[1];
                        DU2_HALF[i] = DT*FLUC_N[i][2]/DUAL[2];
                        DU3_HALF[i] = DT*FLUC_N[i][3]/DUAL[3];
                }
#endif

#ifdef BLENDED 
                for(i=0;i<5;i++){
                        DU0_HALF[i] = DT*FLUC_B[i][0]/DUAL[0];
                        DU1_HALF[i] = DT*FLUC_B[i][1]/DUAL[1];
                        DU2_HALF[i] = DT*FLUC_B[i][2]/DUAL[2];
                        DU3_HALF[i] = DT*FLUC_B[i][3]/DUAL[3];
                }
#endif
                return ;
        }

        void pass_update_half(){
                VERTEX_0->update_du_half(DU0_HALF);
                VERTEX_1->update_du_half(DU1_HALF);
                VERTEX_2->update_du_half(DU2_HALF);
                VERTEX_3->update_du_half(DU3_HALF);
                return ;
        }

        //**********************************************************************************************************************

        void calculate_second_half(double T, double DT){
                int i,j,m,p;
                double INFLOW[5][5][4][3];
                double DU0[5],DU1[5],DU2[5],DU3[5];

                // double DT = DT_TOT;

                setup_half_state();

#ifdef FIRST_ORDER
                for(i=0;i<5;i++){
                        DU0[i] = 0.0;
                        DU1[i] = 0.0;
                        DU2[i] = 0.0;
                        DU3[i] = 0.0;
                }

                VERTEX_0->update_du(DU0);
                VERTEX_1->update_du(DU1);
                VERTEX_2->update_du(DU2);
                VERTEX_3->update_du(DU3);

                return ;
#endif

                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // Calculate inflow parameters

                double H[4];
                double RHO,C,VX,VX_C,VY,VY_C,VZ,VZ_C,H_AVG,H_C,ALPHA,ALPHA_C,W,NORM=1.0/3.0;
                double Z[5][4];
                double VALUE1,VALUE2,VALUE3,VALUE12,VALUE123;
                double PRESSURE_AVG;
                double LAMBDA[5][4],LAMBDA_PLUS[5][4],LAMBDA_MINUS[5][4];
                double N_X[4],N_Y[4],N_Z[4];

                double Z_BAR[5],W_HAT[5][4];

                // Construct Roe vector Z

                for(m=0;m<4;++m){
                        Z[0][m] = sqrt(U_HALF[0][m]);
                        Z[1][m] = U_HALF[1][m]/Z[0][m];
                        Z[2][m] = U_HALF[2][m]/Z[0][m];
                        Z[3][m] = U_HALF[3][m]/Z[0][m];
                        Z[4][m] = (U_HALF[4][m] + PRESSURE_HALF[m])/Z[0][m];

                        N_X[m]  = NORMAL[m][0];
                        N_Y[m]  = NORMAL[m][1];
                        N_Z[m]  = NORMAL[m][2];

                        H[m] = (U_HALF[4][m] + PRESSURE_HALF[m])/U_HALF[0][m];
                }

                for(i=0; i<5; ++i){Z_BAR[i] = (Z[i][0] + Z[i][1] + Z[i][2] + Z[i][3])/4.0;}


                for(m=0; m<4; ++m){
                        W_HAT[0][m] =  2.0*Z_BAR[0]*Z[0][m];
                        W_HAT[1][m] =  Z_BAR[1]*Z[0][m] + Z_BAR[0]*Z[1][m];
                        W_HAT[2][m] =  Z_BAR[2]*Z[0][m] + Z_BAR[0]*Z[2][m];
                        W_HAT[3][m] =  Z_BAR[3]*Z[0][m] + Z_BAR[0]*Z[3][m];
                        W_HAT[4][m] = (Z_BAR[4]*Z[0][m] + GAMMA_1*Z_BAR[1]*Z[1][m] + GAMMA_1*Z_BAR[2]*Z[2][m] + GAMMA_1*Z_BAR[3]*Z[3][m] + Z_BAR[0]*Z[4][m])/GAMMA;
                        //std::cout << Z_BAR[1] << "\t" << Z[0][m] << "\t" << Z_BAR[0] << "\t" << Z[2][m] << std::endl;
                }

                // Construct average state for element

                RHO   = pow((sqrt(U_HALF[0][0]) + sqrt(U_HALF[0][1]) + sqrt(U_HALF[0][2]) + sqrt(U_HALF[0][3]))/4.0, 2);
                VX    = (sqrt(U_HALF[0][0])*U_HALF[1][0]/U_HALF[0][0] + sqrt(U_HALF[0][1])*U_HALF[1][1]/U_HALF[0][1] + sqrt(U_HALF[0][2])*U_HALF[1][2]/U_HALF[0][2] + sqrt(U_HALF[0][3])*U_HALF[1][3]/U_HALF[0][3])/(sqrt(U_HALF[0][0]) + sqrt(U_HALF[0][1]) + sqrt(U_HALF[0][2]) + sqrt(U_HALF[0][3]));
                VY    = (sqrt(U_HALF[0][0])*U_HALF[2][0]/U_HALF[0][0] + sqrt(U_HALF[0][1])*U_HALF[2][1]/U_HALF[0][1] + sqrt(U_HALF[0][2])*U_HALF[2][2]/U_HALF[0][2] + sqrt(U_HALF[0][3])*U_HALF[2][3]/U_HALF[0][3])/(sqrt(U_HALF[0][0]) + sqrt(U_HALF[0][1]) + sqrt(U_HALF[0][2]) + sqrt(U_HALF[0][3]));
                VZ    = (sqrt(U_HALF[0][0])*U_HALF[3][0]/U_HALF[0][0] + sqrt(U_HALF[0][1])*U_HALF[3][1]/U_HALF[0][1] + sqrt(U_HALF[0][2])*U_HALF[3][2]/U_HALF[0][2] + sqrt(U_HALF[0][3])*U_HALF[3][3]/U_HALF[0][3])/(sqrt(U_HALF[0][0]) + sqrt(U_HALF[0][1]) + sqrt(U_HALF[0][2]) + sqrt(U_HALF[0][3]));
                H_AVG = (sqrt(U_HALF[0][0])*H[0] + sqrt(U_HALF[0][1])*H[1] + sqrt(U_HALF[0][2])*H[2] + sqrt(U_HALF[0][3])*H[3])/(sqrt(U_HALF[0][0]) + sqrt(U_HALF[0][1]) + sqrt(U_HALF[0][2]) + sqrt(U_HALF[0][3]));

                PRESSURE_AVG = (PRESSURE_HALF[0] + PRESSURE_HALF[1] + PRESSURE_HALF[2])/4.0;
                C = sqrt((GAMMA-1.0) * H_AVG - (GAMMA-1.0) * (VX*VX + VY*VY + VZ*VZ)/2.0);
                // if(C <= C_LIM){C = C_LIM;}

                // Reassign variables to local equivalents

                VX_C = VX/C;
                VY_C = VY/C;
                VZ_C = VZ/C;

                H_C = H_AVG/C;

                ALPHA   = GAMMA_1*(VX*VX + VY*VY + VZ*VZ)/2.0;
                ALPHA_C = ALPHA/C;

                // Calculate K+,K- and K matrices for each vertex i,j,k

                for(m=0;m<4;++m){

                        W = VX*N_X[m] + VY*N_Y[m] + VZ*N_Z[m];

                        LAMBDA[0][m] = W + C;
                        LAMBDA[1][m] = W - C;
                        LAMBDA[2][m] = W;
                        LAMBDA[3][m] = W;
                        LAMBDA[4][m] = W;

                        for(i=0;i<5;++i){
                                LAMBDA_PLUS[i][m]  = max_val(0.0,LAMBDA[i][m]);
                                LAMBDA_MINUS[i][m] = min_val(0.0,LAMBDA[i][m]);
                        } 

                        for(p=0;p<3;++p){
                                if(p==0){       // Identify and select positive eigenvalues
                                        VALUE1 = LAMBDA_PLUS[0][m];
                                        VALUE2 = LAMBDA_PLUS[1][m];
                                        VALUE3 = LAMBDA_PLUS[2][m];
#ifdef DEBUG
                                        std::cout << "K+" << std::endl;
#endif
                                }else if(p==1){ // Identify and select negative eigenvalues
                                        VALUE1 = LAMBDA_MINUS[0][m];
                                        VALUE2 = LAMBDA_MINUS[1][m];
                                        VALUE3 = LAMBDA_MINUS[2][m];
#ifdef DEBUG
                                        std::cout << "K-" << std::endl;
#endif
                                }else{          // Select all eigenvalues
                                        VALUE1 = LAMBDA[0][m];
                                        VALUE2 = LAMBDA[1][m];
                                        VALUE3 = LAMBDA[2][m];
#ifdef DEBUG
                                        std::cout << "K" << std::endl;
#endif
                                }

                                VALUE12  = (VALUE1 - VALUE2)/2.0;
                                VALUE123 = (VALUE1 + VALUE2 - 2.0*VALUE3)/2.0;

                                INFLOW[0][0][m][p] = NORM * MAG[m]*(ALPHA_C*VALUE123/C - W*VALUE12/C + VALUE3);
                                INFLOW[0][1][m][p] = NORM * MAG[m]*(-1.0*GAMMA_1*VX_C*VALUE123/C + N_X[m]*VALUE12/C);
                                INFLOW[0][2][m][p] = NORM * MAG[m]*(-1.0*GAMMA_1*VY_C*VALUE123/C + N_Y[m]*VALUE12/C);
                                INFLOW[0][3][m][p] = NORM * MAG[m]*(-1.0*GAMMA_1*VZ_C*VALUE123/C + N_Z[m]*VALUE12/C);
                                INFLOW[0][4][m][p] = NORM * MAG[m]*(GAMMA_1*VALUE123/(C*C));

                                INFLOW[1][0][m][p] = NORM * MAG[m]*((ALPHA_C*VX_C - W*N_X[m])*VALUE123 + (ALPHA_C*N_X[m] - VX_C*W)*VALUE12);
                                INFLOW[1][1][m][p] = NORM * MAG[m]*((N_X[m]*N_X[m] - GAMMA_1*VX_C*VX_C)*VALUE123 - (GAMMA_2*VX_C*N_X[m]*VALUE12) + VALUE3);
                                INFLOW[1][2][m][p] = NORM * MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*VX_C*VY_C)*VALUE123 + (VX_C*N_Y[m] - GAMMA_1*VY_C*N_X[m])*VALUE12);
                                INFLOW[1][3][m][p] = NORM * MAG[m]*((N_X[m]*N_Z[m] - GAMMA_1*VX_C*VZ_C)*VALUE123 + (VX_C*N_Z[m] - GAMMA_1*VZ_C*N_X[m])*VALUE12);
                                INFLOW[1][4][m][p] = NORM * MAG[m]*(GAMMA_1*VX_C*VALUE123/C + GAMMA_1*N_X[m]*VALUE12/C);

                                INFLOW[2][0][m][p] = NORM * MAG[m]*((ALPHA_C*VY_C - W*N_Y[m])*VALUE123 + (ALPHA_C*N_Y[m] - VY_C*W)*VALUE12);
                                INFLOW[2][1][m][p] = NORM * MAG[m]*((N_X[m]*N_Y[m] - GAMMA_1*VX_C*VY_C)*VALUE123 + (VY_C*N_X[m] - GAMMA_1*VX_C*N_Y[m])*VALUE12);
                                INFLOW[2][2][m][p] = NORM * MAG[m]*((N_Y[m]*N_Y[m] - GAMMA_1*VY_C*VY_C)*VALUE123 - (GAMMA_2*VY_C*N_Y[m]*VALUE12) + VALUE3);
                                INFLOW[2][3][m][p] = NORM * MAG[m]*((N_Z[m]*N_Y[m] - GAMMA_1*VZ_C*VY_C)*VALUE123 + (VY_C*N_Z[m] - GAMMA_1*VZ_C*N_Y[m])*VALUE12);
                                INFLOW[2][4][m][p] = NORM * MAG[m]*(GAMMA_1*VY_C*VALUE123/C + GAMMA_1*N_Y[m]*VALUE12/C);

                                INFLOW[3][0][m][p] = NORM * MAG[m]*((ALPHA_C*VZ_C - W*N_Z[m])*VALUE123 + (ALPHA_C*N_Z[m] - VZ_C*W)*VALUE12);
                                INFLOW[3][1][m][p] = NORM * MAG[m]*((N_X[m]*N_Z[m] - GAMMA_1*VX_C*VZ_C)*VALUE123 + (VZ_C*N_X[m] - GAMMA_1*VX_C*N_Z[m])*VALUE12);
                                INFLOW[3][2][m][p] = NORM * MAG[m]*((N_Y[m]*N_Z[m] - GAMMA_1*VY_C*VZ_C)*VALUE123 + (VZ_C*N_Y[m] - GAMMA_1*VY_C*N_Z[m])*VALUE12);
                                INFLOW[3][3][m][p] = NORM * MAG[m]*((N_Z[m]*N_Z[m] - GAMMA_1*VZ_C*VZ_C)*VALUE123 - (GAMMA_2*VZ_C*N_Z[m]*VALUE12) + VALUE3);
                                INFLOW[3][4][m][p] = NORM * MAG[m]*(GAMMA_1*VZ_C*VALUE123/C + GAMMA_1*N_Z[m]*VALUE12/C);

                                INFLOW[4][0][m][p] = NORM * MAG[m]*((ALPHA_C*H_C - W*W)*VALUE123 + W*(ALPHA_C - H_C)*VALUE12);
                                INFLOW[4][1][m][p] = NORM * MAG[m]*((W*N_X[m] - VX - ALPHA_C*VX_C)*VALUE123 + (H_C*N_X[m] - GAMMA_1*VX_C*W)*VALUE12);
                                INFLOW[4][2][m][p] = NORM * MAG[m]*((W*N_Y[m] - VY - ALPHA_C*VY_C)*VALUE123 + (H_C*N_Y[m] - GAMMA_1*VY_C*W)*VALUE12);
                                INFLOW[4][3][m][p] = NORM * MAG[m]*((W*N_Z[m] - VZ - ALPHA_C*VZ_C)*VALUE123 + (H_C*N_Z[m] - GAMMA_1*VZ_C*W)*VALUE12);
                                INFLOW[4][4][m][p] = NORM * MAG[m]*(GAMMA_1*H_C*VALUE123/C + GAMMA_1*W*VALUE12/C + VALUE3);
                        }
                }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#if defined(LDA_SCHEME) or defined(BLENDED)

                double PHI_HALF[5];
                double SECOND_FLUC_LDA[5][4];

                for(i=0;i<5;++i){
                        PHI_HALF[i] = 0.0;
                        for(m=0;m<4;++m){
                                PHI_HALF[i] += INFLOW[i][0][m][2]*W_HAT[0][m] + INFLOW[i][1][m][2]*W_HAT[1][m] + INFLOW[i][2][m][2]*W_HAT[2][m] + INFLOW[i][3][m][2]*W_HAT[3][m] + INFLOW[i][4][m][2]*W_HAT[4][m];
                        }
                }

                // Calculate spatial splitting for first half timestep

                for(i=0;i<5;++i){
                        for(m=0;m<4;++m){
                                FLUC_HALF_LDA[i][m] = 0.5*(BETA[i][0][m] * PHI_HALF[0] + BETA[i][1][m] * PHI_HALF[1] + BETA[i][2][m] * PHI_HALF[2] + BETA[i][3][m] * PHI_HALF[3] + BETA[i][4][m] * PHI_HALF[4]);
                        }
                }

                double DIFF[5][4];
                double MASS[5][5][4];
                double MASS_DIFF[5][4];
                double SUM_MASS[5];

                for(i=0;i<5;++i){
                        for(j=0;j<5;++j){
                                for(m=0;m<4;++m){MASS[i][j][m] = VOLUME * BETA[i][j][m]/4.0;}
                        }
                }

                for(i=0;i<5;++i){
                        for(m=0;m<4;++m){
                                DIFF[i][m] = U_HALF[i][m] - U_N[i][m];
                        }
                }

                for(i=0;i<5;++i){
                        for(m=0;m<4;++m){
                                MASS_DIFF[i][m] = MASS[i][0][m] * DIFF[0][m] + MASS[i][1][m] * DIFF[1][m] + MASS[i][2][m] * DIFF[2][m] + MASS[i][3][m] * DIFF[3][m] + MASS[i][4][m] * DIFF[4][m];
                        }
                }

                for(i=0;i<5;i++){
                        SUM_MASS[i] = 0.0;
                        for(m=0;m<4;++m){
                                if(DT==0.0){
                                        SUM_MASS[i] = 0.0;
                                }else{
                                        SUM_MASS[i] += MASS_DIFF[i][m]/DT;
                                }
                        }
                }

                for(i=0;i<5;++i){
                        for(m=0;m<4;++m){
                                SECOND_FLUC_LDA[i][m] = SUM_MASS[i] + 0.5*(FLUC_LDA[i][m] + FLUC_HALF_LDA[i][m]);
                        }
                }
#endif

#if defined(N_SCHEME) or defined(BLENDED)

                double INFLOW_MINUS_SUM[5][5];
                double SECOND_FLUC_N[5][4];

                for(i=0;i<5;++i){
                        for(j=0;j<5;++j){
                                INFLOW_MINUS_SUM[i][j] = 0.0;
                                for(m=0;m<4;++m){
                                        INFLOW_MINUS_SUM[i][j] += INFLOW[i][j][m][1];
                                }
                        }
                }

                mat_inv(&INFLOW_MINUS_SUM[0][0],5,X[0],Y[0],ID,2);

                double AREA_DIFF[5][4];
                double BRACKET[5][4];
                double KZ_SUM[5];

                for(i=0;i<5;++i){
                        KZ_SUM[i] = 0.0;
                        for(m=0;m<4;++m){
                                KZ_SUM[i] += INFLOW[i][0][m][1]*W_HAT[0][m] + INFLOW[i][1][m][1]*W_HAT[1][m] + INFLOW[i][2][m][1]*W_HAT[2][m] + INFLOW[i][3][m][1]*W_HAT[3][m] + INFLOW[i][4][m][1]*W_HAT[4][m];
                        }
                }

                for(i=0;i<5;++i){
                        for(m=0;m<4;++m){
                                BRACKET[i][m] = W_HAT[i][m] - (INFLOW_MINUS_SUM[i][0]*KZ_SUM[0] + INFLOW_MINUS_SUM[i][1]*KZ_SUM[1] + INFLOW_MINUS_SUM[i][2]*KZ_SUM[2] + INFLOW_MINUS_SUM[i][3]*KZ_SUM[3] + INFLOW_MINUS_SUM[i][4]*KZ_SUM[4]);
                        }
                }

                for(i=0;i<5;++i){
                        for(m=0;m<4;++m){
                                FLUC_HALF_N[i][m] = 0.5*(INFLOW[i][0][m][0]*BRACKET[0][m] + INFLOW[i][1][m][0]*BRACKET[1][m] + INFLOW[i][2][m][0]*BRACKET[2][m] + INFLOW[i][3][m][0]*BRACKET[3][m] + INFLOW[i][4][m][0]*BRACKET[4][m]);
                        }
                }

                for(i=0;i<5;++i){
                        for(m=0;m<4;++m){
                                AREA_DIFF[i][m] = VOLUME*(U_HALF[i][m] - U_N[i][m])/4.0;
                        }
                }

                for(i=0;i<5;++i){
                        for(m=0;m<4;++m){
                                if(DT == 0.0){
                                        SECOND_FLUC_N[i][m] = 0.0;
                                }else{
                                        SECOND_FLUC_N[i][m] = AREA_DIFF[i][m]/DT + 0.5*(FLUC_N[i][m] + FLUC_HALF_N[i][m]);
                                }
                        }
                }
#endif

#ifdef LDA_SCHEME
                for(i=0;i<5;i++){
                        DU0[i] = (DT/DUAL[0])*SECOND_FLUC_LDA[i][0];
                        DU1[i] = (DT/DUAL[1])*SECOND_FLUC_LDA[i][1];
                        DU2[i] = (DT/DUAL[2])*SECOND_FLUC_LDA[i][2];
                        DU3[i] = (DT/DUAL[3])*SECOND_FLUC_LDA[i][3];
                }
#endif

#ifdef N_SCHEME
                for(i=0;i<5;i++){
                        DU0[i] = (DT/DUAL[0])*SECOND_FLUC_N[i][0];
                        DU1[i] = (DT/DUAL[1])*SECOND_FLUC_N[i][1];
                        DU2[i] = (DT/DUAL[2])*SECOND_FLUC_N[i][2];
                        DU3[i] = (DT/DUAL[3])*SECOND_FLUC_N[i][3];
                }
#endif

#ifdef BLENDED
                double THETA_E[5][5];
                double IDENTITY[5][5];
                double SUM_FLUC_N[5];

                THETA_E[0][0] = IDENTITY[0][0] = 1.0;
                THETA_E[0][1] = IDENTITY[0][1] = 0.0;
                THETA_E[0][2] = IDENTITY[0][2] = 0.0;
                THETA_E[0][3] = IDENTITY[0][3] = 0.0;
                THETA_E[0][4] = IDENTITY[0][4] = 0.0;

                THETA_E[1][0] = IDENTITY[1][0] = 0.0;
                THETA_E[1][1] = IDENTITY[1][1] = 1.0;
                THETA_E[1][2] = IDENTITY[1][2] = 0.0;
                THETA_E[1][3] = IDENTITY[1][3] = 0.0;
                THETA_E[1][4] = IDENTITY[1][4] = 0.0;

                THETA_E[2][0] = IDENTITY[2][0] = 0.0;
                THETA_E[2][1] = IDENTITY[2][1] = 0.0;
                THETA_E[2][2] = IDENTITY[2][2] = 1.0;
                THETA_E[2][3] = IDENTITY[2][3] = 0.0;
                THETA_E[2][4] = IDENTITY[2][4] = 0.0;

                THETA_E[3][0] = IDENTITY[3][0] = 0.0;
                THETA_E[3][1] = IDENTITY[3][1] = 0.0;
                THETA_E[3][2] = IDENTITY[3][2] = 0.0;
                THETA_E[3][3] = IDENTITY[3][3] = 1.0;
                THETA_E[3][4] = IDENTITY[3][4] = 0.0;

                THETA_E[4][0] = IDENTITY[4][0] = 0.0;
                THETA_E[4][1] = IDENTITY[4][1] = 0.0;
                THETA_E[4][2] = IDENTITY[4][2] = 0.0;
                THETA_E[4][3] = IDENTITY[4][3] = 0.0;
                THETA_E[4][4] = IDENTITY[4][4] = 1.0;

                for(i=0;i<5;i++){
                        SUM_FLUC_N[i] = abs(SECOND_FLUC_N[i][0]) + abs(SECOND_FLUC_N[i][1]) + abs(SECOND_FLUC_N[i][2]) + abs(SECOND_FLUC_N[i][3]);

                        THETA_E[i][i] = abs(PHI[i])/SUM_FLUC_N[i];

                        FLUC_B[i][0] = THETA_E[i][i]*SECOND_FLUC_N[i][0] + (IDENTITY[i][i] - THETA_E[i][i])*SECOND_FLUC_LDA[i][0];
                        FLUC_B[i][1] = THETA_E[i][i]*SECOND_FLUC_N[i][1] + (IDENTITY[i][i] - THETA_E[i][i])*SECOND_FLUC_LDA[i][1];
                        FLUC_B[i][2] = THETA_E[i][i]*SECOND_FLUC_N[i][2] + (IDENTITY[i][i] - THETA_E[i][i])*SECOND_FLUC_LDA[i][2];
                        FLUC_B[i][3] = THETA_E[i][i]*SECOND_FLUC_N[i][3] + (IDENTITY[i][i] - THETA_E[i][i])*SECOND_FLUC_LDA[i][3];

                        DU0[i] = DT*FLUC_B[i][0]/DUAL[0];
                        DU1[i] = DT*FLUC_B[i][1]/DUAL[1];
                        DU2[i] = DT*FLUC_B[i][2]/DUAL[2];
                        DU3[i] = DT*FLUC_B[i][3]/DUAL[3];
                }
#endif

                return ;
        }

        void pass_update(){
                VERTEX_0->update_du(DU0);
                VERTEX_1->update_du(DU1);
                VERTEX_2->update_du(DU2);
                VERTEX_3->update_du(DU3);
        }

        //**********************************************************************************************************************

        // Returns Roe average of left and right states
        // double roe_avg(double L1, double L2, double R1, double R2){
        //         double AVG;
        //         AVG = (sqrt(L1)*L2+sqrt(R1)*R2)/(sqrt(L1)+sqrt(R1));
        //         return AVG;
        // }

        void check_boundary(){
#ifdef PERIODIC_BOUNDARY
                if(BOUNDARY == 1){
                        for(int i=0; i<4; ++i){
                                for(int j=0; j<4; ++j){
                                        if(X[j] - X[i] > BND_TOL*SIDE_LENGTH_X){
                                                X_MOD[i] = X[i] + SIDE_LENGTH_X;
                                        }
                                        if(Y[j] - Y[i] > BND_TOL*SIDE_LENGTH_Y){
                                                Y_MOD[i] = Y[i] + SIDE_LENGTH_Y;
                                        }
                                        if(Z[j] - Z[i] > BND_TOL*SIDE_LENGTH_Z){
                                                Z_MOD[i] = Z[i] + SIDE_LENGTH_Z;
                                        }
                                }
                        }
                }
#endif
                return;
        }

        void check_theta(double THETA){
                if (THETA > 3.148/2.0 and THETA < 3.15/2.0){
                        std::cout << "B ERROR: TRIANGLE EXTREMLY ELONGATED: THETA =\t" << THETA << "\tID =\t" << ID << std::endl;
                        // exit(0);
                }
        }

        void setup_normals(){
                // Calculate normals (just in first timestep for static grid)
                int m;

                setup_positions();

                for(int m=0; m<4; ++m){X_MOD[m] = X[m]; Y_MOD[m] = Y[m]; Z_MOD[m] = Z[m];}

                check_boundary();

                calculate_normals(X_MOD,Y_MOD,Z_MOD);

                return ;

        }

        void calculate_normals(double X[4],double Y[4],double Z[4]){
                int m;
                double ADX,ADY,ADZ;
                double BDX,BDY,BDZ;
                double CDX,CDY,CDZ;
                double CROSSX,CROSSY,CROSSZ;
                double V0[3],V1[3],V2[3],V3[3];
                double V21[3],V31[3],V20[3],V30[3],V10[3],V01[3];
                double PERP[4][3];
                double RPLUS[4],RMINUS[4];

                // calculate volume of tetrahedron and pass 1/4 to each vertex for dual

                V0[0] = X[0];                                // vertices
                V0[1] = Y[0];
                V0[2] = Z[0];

                V1[0] = X[1];
                V1[1] = Y[1];
                V1[2] = Z[1];

                V2[0] = X[2];
                V2[1] = Y[2];
                V2[2] = Z[2];

                V3[0] = X[3];
                V3[1] = Y[3];
                V3[2] = Z[3];

                V21[0] = V2[0] - V1[0];                        // edge vectors
                V21[1] = V2[1] - V1[1];
                V21[2] = V2[2] - V1[2];

                V31[0] = V3[0] - V1[0];
                V31[1] = V3[1] - V1[1];
                V31[2] = V3[2] - V1[2];

                V20[0] = V2[0] - V0[0];
                V20[1] = V2[1] - V0[1];
                V20[2] = V2[2] - V0[2];

                V30[0] = V3[0] - V0[0];
                V30[1] = V3[1] - V0[1];
                V30[2] = V3[2] - V0[2];

                V10[0] = V1[0] - V0[0];
                V10[1] = V1[1] - V0[1];
                V10[2] = V1[2] - V0[2];

                V01[0] = V0[0] - V1[0];
                V01[1] = V0[1] - V1[1];
                V01[2] = V0[2] - V1[2];

                PERP[0][0] = cross_product_x(V21,V31);        // normal opposite vertex 0
                PERP[0][1] = cross_product_y(V21,V31);
                PERP[0][2] = cross_product_z(V21,V31);

                PERP[1][0] = cross_product_x(V20,V30);        // normal opposite vertex 1
                PERP[1][1] = cross_product_y(V20,V30);
                PERP[1][2] = cross_product_z(V20,V30);

                PERP[2][0] = cross_product_x(V10,V30);        // normal opposite vertex 2
                PERP[2][1] = cross_product_y(V10,V30);
                PERP[2][2] = cross_product_z(V10,V30);

                PERP[3][0] = cross_product_x(V20,V10);        // normal opposite vertex 3
                PERP[3][1] = cross_product_y(V20,V10);
                PERP[3][2] = cross_product_z(V20,V10);

                // Check if normals point inwards or outwards

                double DOT001,DOT110,DOT220,DOT330;
                double THETA0,THETA1,THETA2,THETA3;

                DOT001 = PERP[0][0]*V01[0] + PERP[0][1]*V01[1] + PERP[0][2]*V01[2];
                DOT110 = PERP[1][0]*V10[0] + PERP[1][1]*V10[1] + PERP[1][2]*V10[2];
                DOT220 = PERP[2][0]*V20[0] + PERP[2][1]*V20[1] + PERP[2][2]*V20[2];
                DOT330 = PERP[3][0]*V30[0] + PERP[3][1]*V30[1] + PERP[3][2]*V30[2];

                THETA0 = acos(DOT001/(sqrt(PERP[0][0]*PERP[0][0] + PERP[0][1]*PERP[0][1] + PERP[0][2]*PERP[0][2])*sqrt(V01[0]*V01[0] + V01[1]*V01[1] + V01[2]*V01[2])));
                THETA1 = acos(DOT110/(sqrt(PERP[1][0]*PERP[1][0] + PERP[1][1]*PERP[1][1] + PERP[1][2]*PERP[1][2])*sqrt(V10[0]*V10[0] + V10[1]*V10[1] + V10[2]*V10[2])));
                THETA2 = acos(DOT220/(sqrt(PERP[2][0]*PERP[2][0] + PERP[2][1]*PERP[2][1] + PERP[2][2]*PERP[2][2])*sqrt(V20[0]*V20[0] + V20[1]*V20[1] + V20[2]*V20[2])));
                THETA3 = acos(DOT330/(sqrt(PERP[3][0]*PERP[3][0] + PERP[3][1]*PERP[3][1] + PERP[3][2]*PERP[3][2])*sqrt(V30[0]*V30[0] + V30[1]*V30[1] + V30[2]*V30[2])));

                // std::cout << "n0\t" << PERP[0][0] << "\t" << PERP[0][1] << "\t" << PERP[0][2] << std::endl;
                // std::cout << "n1\t" << PERP[1][0] << "\t" << PERP[1][1] << "\t" << PERP[1][2] << std::endl;
                // std::cout << "n2\t" << PERP[2][0] << "\t" << PERP[2][1] << "\t" << PERP[2][2] << std::endl;
                // std::cout << "n3\t" << PERP[3][0] << "\t" << PERP[3][1] << "\t" << PERP[3][2] << std::endl;
                // std::cout << std::endl;

                // check_theta(THETA0);
                // check_theta(THETA1);
                // check_theta(THETA2);
                // check_theta(THETA3);

                if(THETA0 > 3.1416/2.0){
                        // std::cout << "Flipping 0th Normal" << std::endl;
                        PERP[0][0] = -1.0*PERP[0][0];
                        PERP[0][1] = -1.0*PERP[0][1];
                        PERP[0][2] = -1.0*PERP[0][2];
                }
                if(THETA1 > 3.1416/2.0){
                        // std::cout << "Flipping 1st Normal" << std::endl;
                        PERP[1][0] = -1.0*PERP[1][0];
                        PERP[1][1] = -1.0*PERP[1][1];
                        PERP[1][2] = -1.0*PERP[1][2];
                }
                if(THETA2 > 3.1416/2.0){
                        // std::cout << "Flipping 2nd Normal" << std::endl;
                        PERP[2][0] = -1.0*PERP[2][0];
                        PERP[2][1] = -1.0*PERP[2][1];
                        PERP[2][2] = -1.0*PERP[2][2];
                }
                if(THETA3 > 3.1416/2.0){
                        // std::cout << "Flipping 3rd Normal" << std::endl;
                        PERP[3][0] = -1.0*PERP[3][0];
                        PERP[3][1] = -1.0*PERP[3][1];
                        PERP[3][2] = -1.0*PERP[3][2];
                }

                // std::cout << "V0\t" << V0[0] << "\t" << V0[1] << "\t" << V0[2] << std::endl;
                // std::cout << "V1\t" << V1[0] << "\t" << V1[1] << "\t" << V1[2] << std::endl;
                // std::cout << "V2\t" << V2[0] << "\t" << V2[1] << "\t" << V2[2] << std::endl;
                // std::cout << "V3\t" << V3[0] << "\t" << V3[1] << "\t" << V3[2] << std::endl;
                // std::cout << std::endl;

                // std::cout << "n0\t" << PERP[0][0] << "\t" << PERP[0][1] << "\t" << PERP[0][2] << std::endl;
                // std::cout << "n1\t" << PERP[1][0] << "\t" << PERP[1][1] << "\t" << PERP[1][2] << std::endl;
                // std::cout << "n2\t" << PERP[2][0] << "\t" << PERP[2][1] << "\t" << PERP[2][2] << std::endl;
                // std::cout << "n3\t" << PERP[3][0] << "\t" << PERP[3][1] << "\t" << PERP[3][2] << std::endl;
                // std::cout << std::endl;

                ADX = X[0] - X[3];
                ADY = Y[0] - Y[3];
                ADZ = Z[0] - Z[3];

                BDX = X[1] - X[3];
                BDY = Y[1] - Y[3];
                BDZ = Z[1] - Z[3];

                CDX = X[2] - X[3];
                CDY = Y[2] - Y[3];
                CDZ = Z[2] - Z[3];

                CROSSX = BDY*CDZ - BDZ*CDY;
                CROSSY = -1.0*(BDX*CDZ - BDZ*CDX);
                CROSSZ = BDX*CDY - BDY*CDX;

                VOLUME = std::abs(ADX*CROSSX + ADY*CROSSY + ADZ*CROSSZ)/6.0;

                // std::cout << VOLUME << std::endl;
                // exit(0);

                VERTEX_0->calculate_dual(VOLUME/4.0);
                VERTEX_1->calculate_dual(VOLUME/4.0);
                VERTEX_2->calculate_dual(VOLUME/4.0);
                VERTEX_3->calculate_dual(VOLUME/4.0);

                for(m=0;m<4;++m){
                        MAG[m] = sqrt(PERP[m][0]*PERP[m][0] + PERP[m][1]*PERP[m][1] + PERP[m][2]*PERP[m][2]);
                        NORMAL[m][0] = PERP[m][0]/MAG[m];
                        NORMAL[m][1] = PERP[m][1]/MAG[m];
                        NORMAL[m][2] = PERP[m][2]/MAG[m];
                        // std::cout << NORMAL[i][0] << "\t" << NORMAL[i][1] << "\t" << NORMAL[i][2] << "\t" << MAG[i] << std::endl;
                }

                return ;
        }

        double area_triangle(double V0[3], double V1[3], double V2[3]){
                double AREA;

                double L01 = sqrt((V1[0] - V0[0])*(V1[0] - V0[0]) + (V1[1] - V0[1])*(V1[1] - V0[1]) + (V1[2] - V0[2])*(V1[2] - V0[2]));
                double L02 = sqrt((V2[0] - V0[0])*(V2[0] - V0[0]) + (V2[1] - V0[1])*(V2[1] - V0[1]) + (V2[2] - V0[2])*(V2[2] - V0[2]));

                double THETA = acos(((V1[0] - V0[0])*(V2[0] - V0[0]) + (V1[1] - V0[1])*(V2[1] - V0[1]) + (V1[2] - V0[2])*(V2[2] - V0[2]))/(L01*L02));

                if(THETA > 3.14159/2.0){THETA = 3.14159 - THETA;}

                AREA = 0.5 * (sqrt((V1[0] - V0[0])*(V1[0] - V0[0]) + (V1[1] - V0[1])*(V1[1] - V0[1]) + (V1[2] - V0[2])*(V1[2] - V0[2])) * sqrt((V2[0] - V0[0])*(V2[0] - V0[0]) + (V2[1] - V0[1])*(V2[1] - V0[1]) + (V2[2] - V0[2])*(V2[2] - V0[2]))) * sin(THETA);

                // std::cout << AREA << std::endl;

                return AREA;
        }

        void calculate_len_vel_contribution(){
                int m;
                double A0,A1,A2,A3;
                double H,VX,VY,VZ,VEL[4];
                double C_SOUND[4];
                double AMAX,VMAX,CONT;

                setup_initial_state();

                double V0[3],V1[3],V2[3],V3[3];

                V0[0] = X_MOD[0];
                V0[1] = Y_MOD[0];
                V0[2] = Z_MOD[0];

                V1[0] = X_MOD[1];
                V1[1] = Y_MOD[1];
                V1[2] = Z_MOD[1];

                V2[0] = X_MOD[2];
                V2[1] = Y_MOD[2];
                V2[2] = Z_MOD[2];

                V3[0] = X_MOD[3];
                V3[1] = Y_MOD[3];
                V3[2] = Z_MOD[3];

                A0 = area_triangle(V1,V2,V3);
                A1 = area_triangle(V0,V2,V3);
                A2 = area_triangle(V0,V1,V3);
                A3 = area_triangle(V0,V1,V2);

                AMAX = max_val(A0,A1);
                AMAX = max_val(AMAX,A2);
                AMAX = max_val(AMAX,A3);

                // std::cout << LMAX << std::endl;

                for(m=0;m<4;++m){
                        H = (U_N[4][m] + PRESSURE[m])/U_N[0][m];
                        VX = U_N[1][m]/U_N[0][m];
                        VY = U_N[2][m]/U_N[0][m];
                        VZ = U_N[3][m]/U_N[0][m];
                        VEL[m] = sqrt(VX*VX + VY*VY + VZ*VZ);
                        C_SOUND[m] = sqrt((GAMMA-1.0) * H - (GAMMA-1.0) * (VX*VX + VY*VY + VZ*VZ)/2.0);
                }

                VMAX = max_val((VEL[0] + C_SOUND[0]),(VEL[1] + C_SOUND[1]));
                VMAX = max_val(VMAX,(VEL[2] + C_SOUND[2]));
                VMAX = max_val(VMAX,(VEL[3] + C_SOUND[3]));

                CONT = AMAX * VMAX;

                VERTEX_0->update_len_vel_sum(CONT);
                VERTEX_1->update_len_vel_sum(CONT);
                VERTEX_2->update_len_vel_sum(CONT);
                VERTEX_3->update_len_vel_sum(CONT);

                return ;
        }

        double cross_product_x(double A[3], double B[3]){
                // std::cout << A[0] << "\t" << A[1] << "\t" << A[2] << "\t" << B[0] << "\t" << B[1] << "\t" << B[2] << "\t" <<  A[1]*B[2] - A[2]*B[1] << std::endl;
                return A[1]*B[2] - A[2]*B[1];
        }
        double cross_product_y(double A[3], double B[3]){
                return -1.0*(A[0]*B[2] - A[2]*B[0]);
        }
        double cross_product_z(double A[3], double B[3]){
                return A[0]*B[1] - A[1]*B[0];
        }

        void reorder_vertices(){
                VERTEX *TEMP_VERTEX;

                TEMP_VERTEX = VERTEX_1;
                VERTEX_1 = VERTEX_2;
                VERTEX_2 = TEMP_VERTEX;

                return;
        }

};
