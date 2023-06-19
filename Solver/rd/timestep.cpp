/*
 * This file was written by Ben Morton (bmorton@ed.ac.uk) and Zhenyu Wu (zhenyu.wu@ed.ac.uk).
 */

#include <vector>
#include "triangle2D.h"
#include "base.h"

void drift_update_half(int TBIN_CURRENT, int N_TRIANG, double T, double DT, std::vector<TRIANGLE> &RAND_MESH , double P_DIFF, double Vmean, double V_DIFF){
        int TBIN;
        for(int j=0;j<N_TRIANG;++j){                                                                         // loop over all triangles in MESH
                // printf("%i\n",j);
                TBIN = RAND_MESH[j].get_tbin();
                // printf("%i\t%i\t%i\n",TBIN_CURRENT,TBIN);
                if(TBIN == 0){
                        RAND_MESH[j].calculate_first_half(T,DT,P_DIFF,Vmean,V_DIFF);
                }else if(TBIN_CURRENT % TBIN == 0){
                        RAND_MESH[j].calculate_first_half(T,DT,P_DIFF,Vmean,V_DIFF);
                }
                RAND_MESH[j].pass_update_half();
        }
}

#ifdef JUMP
void jump_update_half(int TBIN_CURRENT, int N_TRIANG, double T, double DT, std::vector<TRIANGLE> &RAND_MESH){
        int TBIN;
        for(int j=0;j<N_TRIANG;++j){                                                                         // loop over all triangles in MESH
                        TBIN = RAND_MESH[j].get_tbin();
                        if(TBIN_CURRENT == 0 or (TBIN_CURRENT == 1 and  TBIN == 1)\
                                             or (TBIN_CURRENT == 2 and (TBIN == 2 or TBIN == 1))\
                                             or (TBIN_CURRENT == 3 and  TBIN == 1)\
                                             or (TBIN_CURRENT == 4 and (TBIN == 4 or TBIN == 2 or TBIN == 1))\
                                             or (TBIN_CURRENT == 5 and  TBIN == 1)\
                                             or (TBIN_CURRENT == 6 and (TBIN == 2 or TBIN == 1))\
                                             or (TBIN_CURRENT == 7 and  TBIN == 1)\
                                             ){
                                RAND_MESH[j].calculate_first_half(T);
                        }
                        if(MAX_TBIN == 1){
                                RAND_MESH[j].pass_update_half(DT);
                        }else if(MAX_TBIN == 2){
                                if(TBIN_CURRENT == 0 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 1 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 1 and (TBIN == 2 or TBIN == 4 or TBIN == 8)){RAND_MESH[j].pass_update_half(2.0*DT);}
                        }else if(MAX_TBIN == 4){
                                if(TBIN_CURRENT == 0 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 1 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 1 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}

                                if(TBIN_CURRENT == 3 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 3 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}
                                if(TBIN_CURRENT == 3 and (TBIN == 4 or TBIN == 8)){RAND_MESH[j].pass_update_half(4.0*DT);}
                        }else if(MAX_TBIN == 8){
                                if(TBIN_CURRENT == 0 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 1 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 1 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}

                                if(TBIN_CURRENT == 2 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 3 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 3 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}
                                if(TBIN_CURRENT == 3 and TBIN == 4){RAND_MESH[j].pass_update_half(4.0*DT);}

                                if(TBIN_CURRENT == 4 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 5 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 5 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}

                                if(TBIN_CURRENT == 6 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 7 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 7 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}
                                if(TBIN_CURRENT == 7 and TBIN == 4){RAND_MESH[j].pass_update_half(4.0*DT);}
                                if(TBIN_CURRENT == 7 and TBIN == 8){RAND_MESH[j].pass_update_half(8.0*DT);}
                        }
                }
}
#endif

void drift_update(int TBIN_CURRENT, int N_TRIANG, double T, double DT, std::vector<TRIANGLE> &RAND_MESH, double P_DIFF, double Vmean, double V_DIFF){
        int TBIN;
        for(int j=0;j<N_TRIANG;++j){                                                                         // loop over all triangles in MESH
                TBIN = RAND_MESH[j].get_tbin();
                if(TBIN == 0){
                        RAND_MESH[j].calculate_second_half(T,DT,P_DIFF,Vmean,V_DIFF);
                }else if(TBIN_CURRENT % TBIN == 0){
                        RAND_MESH[j].calculate_second_half(T,DT,P_DIFF,Vmean,V_DIFF);
                }
                RAND_MESH[j].pass_update();
        }
}

void reset_tbins(double T, double DT, int N_TRIANG, int N_POINTS, double &NEXT_DT, std::vector<TRIANGLE> &RAND_MESH, std::vector<VERTEX> &RAND_POINTS){
        double POSSIBLE_DT, MIN_DT;
       	
        for(int j=0;j<N_TRIANG;++j){                                       // loop over all triangles in MESH
                RAND_MESH[j].calculate_len_vel_contribution();             // calculate contribution from each edge TRIANGLE
        }
        NEXT_DT = T_TOT - (T + DT);                                        // set next timestep to max possible value (time remaining to end)å
        for(int i=0;i<N_POINTS;++i){                                       // loop over all vertices
                POSSIBLE_DT = RAND_POINTS[i].calc_next_dt();               // calculate next timestep based on new state
                //printf("%f\n",POSSIBLE_DT );
                if(POSSIBLE_DT < NEXT_DT){NEXT_DT = POSSIBLE_DT;}              
                RAND_POINTS[i].reset_len_vel_sum();
                RAND_POINTS[i].set_tbin_local(MAX_TBIN);
        }
        for(int j=0;j<N_TRIANG;++j){                                        // bin triangles by minimum timestep of vertices
                MIN_DT = RAND_MESH[j].get_vertex_0()->get_dt_req();
                if(RAND_MESH[j].get_vertex_1()->get_dt_req() < MIN_DT){MIN_DT = RAND_MESH[j].get_vertex_1()->get_dt_req();}
                if(RAND_MESH[j].get_vertex_2()->get_dt_req() < MIN_DT){MIN_DT = RAND_MESH[j].get_vertex_2()->get_dt_req();}
#ifdef THREE_D
                if(RAND_MESH[j].get_vertex_3()->get_dt_req() < MIN_DT){MIN_DT = RAND_MESH[j].get_vertex_3()->get_dt_req();}
#endif
                RAND_MESH[j].set_tbin( min_val( MAX_TBIN,int(pow(2.0,int(log2(MIN_DT/NEXT_DT)) ) ) ) );
#ifdef DRIFT_SHELL
                RAND_MESH[j].send_tbin_limit();
#endif
        }
#ifdef DRIFT_SHELL
        for(int j=0;j<N_TRIANG;++j){
                RAND_MESH[j].check_tbin();
        }
#endif
}








