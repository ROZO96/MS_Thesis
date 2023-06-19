/*
 * This file was written by Ben Morton (bmorton@ed.ac.uk) and Zhenyu Wu (zhenyu.wu@ed.ac.uk).
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <omp.h> 

#include "constants.h"
#include "all_functions.h"

#ifdef TWO_D
#include "vertex2D.h"
#include "triangle2D.h"
#include "source2D.cpp"
#include "timestep.cpp"
#include "io2D.cpp"
#include "setup2D.cpp"

#endif

#ifdef SEDOV
    int POINT_CHECK = 0;
#endif

int main(int ARGC, char *ARGV[]){
        /*
        Setup and run simulation from input file constants.h, using precalculated triangulation
        */
        int i, j, l = 0, m;                                        // ******* declare variables and vectors ******
        int SNAP_ID = 0;
        double DT, T = 0.0;                                        //
        double NEXT_TIME = 0.0;                                    // NEXT_TIME    = time of next snapshot
        double NEXT_DT = T_TOT, POSSIBLE_DT = T_TOT;               // NEXT_DT.     = timestep for upcoming time iteration
        double MIN_DT;
        VERTEX                               NEW_VERTEX;           // NEW_VERTEX   = dummy variable for setting up vertices
        TRIANGLE                             NEW_TRIANGLE;         // NEW_TRIABLE  = dummy variable for setting up triangles
        std::vector<VERTEX>                  RAND_POINTS;          // X_POINTS     = vector of x vertices
        std::vector<TRIANGLE>                RAND_MESH;            // X_MESH       = vector of x triangles

        // Initialise seed for random number generator (rand)
        std::srand(68315);

        // read_parameter_file(ARGC, ARGV);

        /****** Setup simulation options ******/

        printf("*********************************************************\n");

        printf("LAIRDS 2D\n");

#ifdef LDA_SCHEME
        printf("Using LDA Scheme\n");
#endif

#ifdef N_SCHEME
        printf("Using N Scheme\n");
#endif

#ifdef BLENDED
        printf("Using B Scheme\n");
#endif

#ifdef FIRST_ORDER
        printf("Using 1st order\n");
#else
        printf("Using 2nd order\n");
#endif

        printf("Building vertices and mesh\n");

        std::ofstream LOGFILE;
        LOGFILE << std::setprecision(12);
        LOGFILE.open(LOG_DIR);

#ifdef READ_IC
#ifdef QHULL_IC
        int N_POINTS, N_TRIANG;
        std::string   POSITIONS_FILE_NAME, TRIANGLES_FILE_NAME;
        std::ifstream POSITIONS_FILE, TRIANGLES_FILE;

        /****** Setup vertices ******/

        printf("Reading QHULL vertex positions ...\n");

        POSITIONS_FILE_NAME = "triangulation/points.txt";
        TRIANGLES_FILE_NAME = "triangulation/ordered_triangles.txt";

        POSITIONS_FILE.open(POSITIONS_FILE_NAME);
        TRIANGLES_FILE.open(TRIANGLES_FILE_NAME);

        N_POINTS = qhull_read_positions_header(POSITIONS_FILE);
        N_TRIANG = qhull_read_triangles_header(TRIANGLES_FILE);

        printf("Number of vertices = %d\n", N_POINTS);

        for(i=0; i<N_POINTS; ++i){
                NEW_VERTEX = qhull_read_positions_line(POSITIONS_FILE);
                NEW_VERTEX.reset_len_vel_sum();
                RAND_POINTS.push_back(NEW_VERTEX);
        }


        /****** Setup mesh ******/

        printf("Reading QHULL triangles ...");

        printf("Number of triangles = %d\n", N_TRIANG);

        for(j=0; j<N_TRIANG; ++j){
                NEW_TRIANGLE = qhull_read_triangles_line(TRIANGLES_FILE,RAND_POINTS);
                NEW_TRIANGLE.set_tbin(1);
                RAND_MESH.push_back(NEW_TRIANGLE);
        }

        POSITIONS_FILE.close();
        TRIANGLES_FILE.close();
#endif
#ifdef CGAL_IC
        int N_POINTS, N_TRIANG;
        std::string   CGAL_FILE_NAME;
        std::ifstream CGAL_FILE;

        /****** Setup vertices ******/


        printf("Reading CGAL vertex positions ...");

        CGAL_FILE_NAME = "Delaunay2D.txt";
        CGAL_FILE.open(CGAL_FILE_NAME);
        
        N_POINTS = cgal_read_positions_header(CGAL_FILE);
	printf("Delunay archive read \n");
	
        printf("Number of vertices = %d\n", N_POINTS);

        for(i=0; i<N_POINTS; ++i){
                NEW_VERTEX = cgal_read_positions_line(CGAL_FILE);
                NEW_VERTEX.reset_len_vel_sum();
                NEW_VERTEX.set_id(i);
                RAND_POINTS.push_back(NEW_VERTEX);
        }

        /****** Setup mesh ******/

        printf("Reading CGAL triangles ...");

        N_TRIANG = cgal_read_triangles_header(CGAL_FILE);

        printf("Number of triangles = %d\n", N_TRIANG);

        for(j=0; j<N_TRIANG; ++j){
                NEW_TRIANGLE = cgal_read_triangles_line(CGAL_FILE,RAND_POINTS,j);
                NEW_TRIANGLE.set_tbin(1);
                RAND_MESH.push_back(NEW_TRIANGLE);
        }

#endif
#endif

#ifdef SEDOV
        /****** Inject pressure for Sedov test  ******/

        double AREA_CHECK  = 0.0;
        double ETOT = 0.0,ETOT_AIM = 300000.0,PRESSURE_AIM;
        for(i=0; i<N_POINTS; ++i){
                if((RAND_POINTS[i].get_x()-5.0)*(RAND_POINTS[i].get_x()-5.0) + (RAND_POINTS[i].get_y()-5.0)*(RAND_POINTS[i].get_y()-5.0) < R_BLAST*R_BLAST){
                        AREA_CHECK = AREA_CHECK + RAND_POINTS[i].get_dual();
                }
        }
        for(i=0; i<N_POINTS; ++i){
                if((RAND_POINTS[i].get_x()-5.0)*(RAND_POINTS[i].get_x()-5.0) + (RAND_POINTS[i].get_y()-5.0)*(RAND_POINTS[i].get_y()-5.0) < R_BLAST*R_BLAST){
                        PRESSURE_AIM = (ETOT_AIM * GAMMA_1 / RAND_POINTS[i].get_dual()) * (RAND_POINTS[i].get_dual() / (AREA_CHECK));
                        RAND_POINTS[i].set_pressure(PRESSURE_AIM);
                        ETOT = ETOT + RAND_POINTS[i].get_pressure()*RAND_POINTS[i].get_dual()/GAMMA_1;
                        printf("%f\t%f\t%f\n", PRESSURE_AIM, RAND_POINTS[i].get_pressure(), ETOT);
                        RAND_POINTS[i].setup_specific_energy();
                        RAND_POINTS[i].prim_to_con();
                }
        }
#endif

        /****** Set initial timestep  ******/

        printf("Finding initial timestep ...\n");

        for(j=0;j<N_TRIANG;++j){                                           // loop over all triangles in MESH
                RAND_MESH[j].calculate_len_vel_contribution();             // calculate flux through TRIANGLE
        }

	double P_MAX,P_MIN,POSSIBLE_P, P_DIFF, Vmean_x, Vmean_y, Vmean,V_MIN,V_MAX,V_DIFF,V_MAG;
	P_MIN=10000000;
	P_MAX=-100000000;
	Vmean=0;
	Vmean_x=0;
	Vmean_y=0;
	V_MIN=10000000;
	V_MAX=-100000000;
        for(i=0; i<N_POINTS; ++i){
        	
                NEXT_DT = RAND_POINTS[i].calc_next_dt();
                POSSIBLE_P= RAND_POINTS[i].get_pressure();
                //Vmean_x+=RAND_POINTS[i].get_x_velocity();
                //Vmean_y+=RAND_POINTS[i].get_y_velocity();
                Vmean_x=RAND_POINTS[i].get_x_velocity();
                Vmean_y=RAND_POINTS[i].get_y_velocity();
                V_MAG=sqrt(Vmean_x*Vmean_x+Vmean_y*Vmean_y);
                
                //printf("%f\n",NEXT_DT  );      // check dt is min required by CFL
                if(POSSIBLE_DT < NEXT_DT){NEXT_DT=POSSIBLE_DT;}
                if(POSSIBLE_P>P_MAX){P_MAX=POSSIBLE_P;}
                if(POSSIBLE_P<P_MIN){P_MIN=POSSIBLE_P;}
                if(V_MAG>V_MAX){V_MAX=V_MAG;}
                if(V_MAG<V_MIN){V_MIN=V_MAG;}
                
                RAND_POINTS[i].reset_len_vel_sum();
                
                }
        P_DIFF=P_MAX-P_MIN;
        V_DIFF=V_MAX-V_MIN;
        Vmean=sqrt(((Vmean_x/N_POINTS)*(Vmean_x/N_POINTS)) +((Vmean_y/N_POINTS)*(Vmean_y/N_POINTS)));
        printf("V_DIFF=%f\n",V_DIFF);
	//printf("P_DIFF=%f\n",P_DIFF);
	//printf("Vmean=%f,denominator=%f\n",Vmean,Vmean*P_DIFF);
	
        printf("Checking mesh size ...");
        printf("Mesh Size = %d\n",int(RAND_MESH.size()));
        printf("Evolving fluid ...\n");

        int TBIN, TBIN_CURRENT = 0;
        NEXT_DT = 0.000001;                                                            // set first timestep to zero

        /****** Loop over time until total time T_TOT is reached *****************************************************************************************************/
        while(T<T_TOT){

        /****** Update time step to new value ******/
                DT = NEXT_DT;   
                                                                // set timestep based on caclulation from previous timestep

#ifdef FIXED_DT
        /****** Reset time step if fixed ******/
                DT = DT_FIX;
#endif
                printf("STEP =\t%d\tTIME =\t%f\tTIMESTEP =\t%f\t%f/100\n", l, T, DT, 100.0*T/T_TOT);

            /****** Write snapshot *****************************************************************************************************/
                if(T >= NEXT_TIME){                                       // write out densities at given interval
                        write_snap(RAND_POINTS,T,DT,N_POINTS,SNAP_ID,LOGFILE);
                        write_active(RAND_MESH, N_TRIANG, SNAP_ID, TBIN_CURRENT);
                        NEXT_TIME = NEXT_TIME + T_TOT/float(N_SNAP);
                        if(NEXT_TIME > T_TOT){NEXT_TIME = T_TOT;}
                        SNAP_ID ++;
                }

        /****** 1st order update ***************************************************************************************************/

#ifdef DRIFT
                /****** Update residual for active bins (Drift method) ******/
                drift_update_half(TBIN_CURRENT, N_TRIANG, T, DT, RAND_MESH,P_DIFF,Vmean,V_DIFF);
#endif
#ifdef JUMP
                /****** Update residual for active bins (Jump method) ******/
                jump_update_half(TBIN_CURRENT, N_TRIANG, T, DT, RAND_MESH);
#endif




#if !defined(DRIFT) && !defined(JUMP)
#ifdef PARA_RES
                #pragma omp parallel for
#endif
                /****** Update residual for all bins (No adaptive method) ******/
                for(j=0;j<N_TRIANG;++j){                                                                         // loop over all triangles in MESH
                        RAND_MESH[j].calculate_first_half(T,DT,P_DIFF,Vmean,V_DIFF);                                        	// calculate flux through TRIANGLE
                        RAND_MESH[j].pass_update_half();
                }
#endif

#ifdef PARA_UP
                #pragma omp parallel for
#endif

		P_MIN=10000000;
		P_MAX=-100000000;
		V_MIN=10000000;
		V_MAX=-100000000;
		
		Vmean_x=0;
		Vmean_y=0;
                for(i=0;i<N_POINTS;++i){                                       // loop over all vertices
                        RAND_POINTS[i].update_u_half();                        // update the half time state
                        RAND_POINTS[i].reset_du_half();                        // reset du value to zero for next timestep
                        RAND_POINTS[i].check_values_half();
                        RAND_POINTS[i].con_to_prim_half();
                        //Vmean_x+=RAND_POINTS[i].get_x_velocity_half();
                	//Vmean_y+=RAND_POINTS[i].get_y_velocity_half();
                	
                	Vmean_x=RAND_POINTS[i].get_x_velocity_half();
                	Vmean_y=RAND_POINTS[i].get_y_velocity_half();
                	V_MAG=sqrt(Vmean_x*Vmean_x+Vmean_y*Vmean_y);
                        POSSIBLE_P= RAND_POINTS[i].get_pressure_half();
                       	
                        if(POSSIBLE_P>P_MAX){P_MAX=POSSIBLE_P;}
                	if(POSSIBLE_P<P_MIN){P_MIN=POSSIBLE_P;}
                	if(V_MAG>V_MAX){V_MAX=V_MAG;}
                	if(V_MAG<V_MIN){V_MIN=V_MAG;}
                }
                P_DIFF=P_MAX-P_MIN;
                V_DIFF=V_MAX-V_MIN;
                Vmean=sqrt(((Vmean_x/N_POINTS)*(Vmean_x/N_POINTS)) +((Vmean_y/N_POINTS)*(Vmean_y/N_POINTS)));
		//printf("P_DIFF=%f\n",P_DIFF);
		printf("V_DIFF=%f\n",V_DIFF);
		//printf("Vmean=%f,denominator=%f\n",Vmean,Vmean*P_DIFF);
        /****** 2nd order update ***************************************************************************************************/

#ifdef DRIFT
                /****** Update residual for active bins (Drift method) ******/
                drift_update(TBIN_CURRENT, N_TRIANG, T, DT, RAND_MESH,P_DIFF,Vmean,V_DIFF);
#endif

#if !defined(DRIFT) && !defined(JUMP)
#ifdef PARA_RES
                #pragma omp parallel for
#endif
                /****** Update residual for all bins (No adaptive method) ******/
                for(j=0;j<N_TRIANG;++j){                                       // loop over all triangles in MESH
                        RAND_MESH[j].calculate_second_half(T,DT,P_DIFF,Vmean,V_DIFF);             // calculate flux through TRIANGLE
                        RAND_MESH[j].pass_update();
                }
#endif

                sources(RAND_POINTS, DT, N_POINTS);

#ifdef PARA_UP
                #pragma omp parallel for
#endif
		P_MIN=10000000;
		P_MAX=-100000000;
		V_MIN=10000000;
		V_MAX=-100000000;
		Vmean_x=0;
		Vmean_y=0;

                for(i=0;i<N_POINTS;++i){                                       // loop over all vertices
                        RAND_POINTS[i].update_u_variables();                   // update the fluid state at vertex
                        RAND_POINTS[i].reset_du();                             // reset du value to zero for next timestep
                        RAND_POINTS[i].check_values();
                        RAND_POINTS[i].con_to_prim();                          // convert these to their corresponding conserved
                        //Vmean_x+=RAND_POINTS[i].get_x_velocity();
                	//Vmean_y+=RAND_POINTS[i].get_y_velocity();
                	Vmean_x=RAND_POINTS[i].get_x_velocity();
                	Vmean_y=RAND_POINTS[i].get_y_velocity();
                	V_MAG=sqrt(Vmean_x*Vmean_x+Vmean_y*Vmean_y);
                        POSSIBLE_P= RAND_POINTS[i].get_pressure();
                        if(POSSIBLE_P>P_MAX){P_MAX=POSSIBLE_P;}
                	if(POSSIBLE_P<P_MIN){P_MIN=POSSIBLE_P;}
                	if(V_MAG>V_MAX){V_MAX=V_MAG;}
                	if(V_MAG<V_MIN){V_MIN=V_MAG;}
                }
		P_DIFF=P_MAX-P_MIN;
		Vmean=sqrt(((Vmean_x/N_POINTS)*(Vmean_x/N_POINTS)) +((Vmean_y/N_POINTS)*(Vmean_y/N_POINTS)));
		V_DIFF=V_MAX-V_MIN;
		//printf("P_DIFF=%f\n",P_DIFF);
		printf("V_DIFF=%f\n",V_DIFF);
		//printf("Vmean=%f,denominator=%f\n",Vmean,Vmean*P_DIFF);
                if(TBIN_CURRENT == 0){
                        reset_tbins(T, DT, N_TRIANG, N_POINTS, NEXT_DT, RAND_MESH, RAND_POINTS);
                }

// #if defined(FIXED_BOUNDARY) && (defined(NOH) || defined(DF))
//                 for(j=0;j<N_TRIANG;++j){                                         // loop over all triangles in MESH
//                         RAND_MESH[j].check_boundary();                           // calculate flux through TRIANGLE
//                 }
// #endif
                TBIN_CURRENT = (TBIN_CURRENT + 1) % MAX_TBIN;                     // increment time step bin
                T += DT;                                                         // increment time
                l += 1;    // increment step number

        }

        write_snap(RAND_POINTS,T,DT,N_POINTS,SNAP_ID,LOGFILE);


        return 0;
}
