/*
 * This file was written by Ben Morton (bmorton@ed.ac.uk).
 * IO routines to write simple ASCII output for python plotting
        open_files => opens ouptut files to write positions, dens, pressure, vel maps, and column of values for 1D plot
 */

#include <iostream>
#include <fstream>
#include <vertex3D.h>
#include <triangle3D.h>
#include "all_functions3D.h"

void open_snap(std::ofstream &SNAPFILE, int i){
        SNAPFILE.open(OUT_DIR+"snapshot3D_"+std::to_string(i)+".txt");
        return;
}

void open_active(std::ofstream &SNAPFILE, int i){
        SNAPFILE.open("output/active_"+std::to_string(i)+".txt");
        return;
}

void write_snap(std::vector<VERTEX> POINTS, double T, double DT, int N_POINTS, int SNAP_ID, std::ofstream &LOGFILE){
        std::ofstream SNAPFILE;
        open_snap(SNAPFILE,SNAP_ID);
        double TOTAL_DENSITY = 0.0;
        double TOTAL_ENERGY = 0.0;
        SNAPFILE << N_POINTS << "\t" << T << std::endl;
        for(int i=0;i<N_POINTS;++i){
                // write         X                           Y                                 rho                                   v_x                                    v_y                                 p                                      e                                         |S|
                SNAPFILE << POINTS[i].get_x() << "\t" << POINTS[i].get_y() << "\t" << POINTS[i].get_z() << "\t" << POINTS[i].get_mass_density() << "\t" << POINTS[i].get_x_velocity() << "\t" << POINTS[i].get_y_velocity() << "\t" << POINTS[i].get_z_velocity() << "\t" << POINTS[i].get_pressure() << "\t" << POINTS[i].get_specific_energy() << "\t" << POINTS[i].get_dual() << std::endl;
                TOTAL_DENSITY += POINTS[i].get_mass_density()*POINTS[i].get_dual();
                TOTAL_ENERGY += POINTS[i].get_specific_energy()*POINTS[i].get_dual() * POINTS[i].get_mass_density();
        }
        std::cout << "*************************************************************************************************" << std::endl;            // right out time and total density to terminal
        std::cout << "time\t" << T << " \t-> total mass =\t" << TOTAL_DENSITY << " \t-> total energy =\t" << TOTAL_ENERGY << "\ttime step = \t" << DT << std::endl;
        LOGFILE << T << "\t" << TOTAL_DENSITY << "\t" << TOTAL_ENERGY << "\t" << DT << "\t" << N_TBINS << "\t" << N_POINTS << std::endl;
        SNAPFILE.close();
        return;
}

void write_active(std::vector<TRIANGLE> MESH, int N_TRIANG, int SNAP_ID, int TBIN_CURRENT){
        std::ofstream SNAPFILE;
        SNAPFILE << std::setprecision(12);
        double X0,X1,X2,Y0,Y1,Y2;
        open_active(SNAPFILE,SNAP_ID);
        SNAPFILE << N_TRIANG << "\t" << std::endl;
        for(int j=0;j<N_TRIANG;++j){
                if(MESH[j].get_boundary() == 0){
                        X0 = MESH[j].get_vertex_0()->get_x();
                        X1 = MESH[j].get_vertex_1()->get_x();
                        X2 = MESH[j].get_vertex_2()->get_x();
                        Y0 = MESH[j].get_vertex_0()->get_y();
                        Y1 = MESH[j].get_vertex_1()->get_y();
                        Y2 = MESH[j].get_vertex_2()->get_y();
                        // write         X        Y          TBIN
                        SNAPFILE << X0 << "\t" << Y0 << "\t"  << X1 << "\t" << Y1 << "\t"  << X2 << "\t" << Y2 << "\t" << MESH[j].get_tbin() << "\t" << MESH[j].get_un00() << "\t" << MESH[j].get_un01() << "\t" << MESH[j].get_un02() << std::endl;
                }
        }
}

void read_parameter_file(int ARGC, char *ARGV[]){
        printf("Parameter file = %s\n", ARGV[1]);
}

// if using CGAL triangulation file, read vertex header info
#ifdef CGAL_IC
int cgal_read_positions_header(std::ifstream &CGAL_FILE){
        int N_POINTS, XSHEETS, YSHEETS, ZSHEETS;
        float XLOW, YLOW, ZLOW, XHIGH, YHIGH, ZHIGH;

        CGAL_FILE >> XLOW >> YLOW >> ZLOW >> XHIGH >> YHIGH >> ZHIGH;

        // check boundaries match
        if(XLOW < 0.0 or YLOW < 0.0 or int(XHIGH) != SIDE_LENGTH_X or int(YHIGH) != SIDE_LENGTH_Y){
                std::cout << "BWARNING: Exiting on mismatching boundaries. Have (X): " << XLOW  << "\t" << XHIGH << "\tNeed (X): " << 0.0 << "\t" << SIDE_LENGTH_X << std::endl;
                exit(0);
        }

        CGAL_FILE >> XSHEETS >> YSHEETS >> ZSHEETS;

        // check data in 1-sheet format
        if(XSHEETS !=1 or YSHEETS!=1 or ZSHEETS!=1){
                std::cout << "BWARNING: Exiting on CGAL file not in 1-sheet format." << std::endl;
                exit(0);
        }

        CGAL_FILE >> N_POINTS;

        return N_POINTS;
}

// read CGAL triangle header info
int cgal_read_triangles_header(std::ifstream &CGAL_FILE){
        int N_TRIANG;

        CGAL_FILE >> N_TRIANG;

        return N_TRIANG;
}

// read postion of one CGAL vertex
VERTEX cgal_read_positions_line(std::ifstream &CGAL_FILE){
        double X,Y,Z;
        VERTEX NEW_VERTEX;

        CGAL_FILE >> X >> Y >> Z;

        // std::cout << X << "\t" << Y << "\t" << Z << std::endl;

        NEW_VERTEX = setup_vertex(X,Y,Z);

        // std::cout << NEW_VERTEX.get_x() << "\t" << NEW_VERTEX.get_y() << "\t" << NEW_VERTEX.get_z() << std::endl;

        return NEW_VERTEX;
}

// read indices of vertices for one CGAL triangle
TRIANGLE cgal_read_triangles_line(std::ifstream &CGAL_FILE, std::vector<VERTEX> &POINTS, int ID){
        int VERT0,VERT1,VERT2,VERT3;
        TRIANGLE NEW_TRIANGLE;

        CGAL_FILE >> VERT0 >> VERT1 >> VERT2 >> VERT3;

        if(VERT0 == VERT1 or VERT0 == VERT2 or VERT0 == VERT3 or VERT1 == VERT2 or VERT1 == VERT3 or VERT2 == VERT3){
                std::cout << "Repeated vertex in triangle" << VERT0 << "\t" << VERT1 << "\t" << VERT2 << "\t" << VERT3 << std::endl;
                exit(0);
        }

        NEW_TRIANGLE.set_vertex_0(&POINTS[VERT0]);
        NEW_TRIANGLE.set_vertex_1(&POINTS[VERT1]);
        NEW_TRIANGLE.set_vertex_2(&POINTS[VERT2]);
        NEW_TRIANGLE.set_vertex_3(&POINTS[VERT3]);

        NEW_TRIANGLE.set_id(ID);

        // std::cout << POINTS[VERT0].get_x() << "\t" << POINTS[VERT1].get_x() << "\t" << POINTS[VERT2].get_x() << std::endl;

        double X0,X1,X2,X3,Y0,Y1,Y2,Y3,Z0,Z1,Z2,Z3;

        // check if boundary triangle

        X0 = POINTS[VERT0].get_x();
        X1 = POINTS[VERT1].get_x();
        X2 = POINTS[VERT2].get_x();
        X3 = POINTS[VERT3].get_x();

        Y0 = POINTS[VERT0].get_y();
        Y1 = POINTS[VERT1].get_y();
        Y2 = POINTS[VERT2].get_y();
        Y3 = POINTS[VERT3].get_y();

        Z0 = POINTS[VERT0].get_z();
        Z1 = POINTS[VERT1].get_z();
        Z2 = POINTS[VERT2].get_z();
        Z3 = POINTS[VERT3].get_z();

        if(abs(X0 - X1) > BND_TOL*SIDE_LENGTH_X or abs(X0 - X2) > BND_TOL*SIDE_LENGTH_X or abs(X0 - X3) > BND_TOL*SIDE_LENGTH_X or abs(X1 - X2) > BND_TOL*SIDE_LENGTH_X or abs(X1 - X3) > BND_TOL*SIDE_LENGTH_X or abs(X2 - X3) > BND_TOL*SIDE_LENGTH_X or
           abs(Y0 - Y1) > BND_TOL*SIDE_LENGTH_Y or abs(Y0 - Y2) > BND_TOL*SIDE_LENGTH_Y or abs(Y0 - Y3) > BND_TOL*SIDE_LENGTH_Y or abs(Y1 - Y2) > BND_TOL*SIDE_LENGTH_Y or abs(Y1 - Y3) > BND_TOL*SIDE_LENGTH_Y or abs(Y2 - Y3) > BND_TOL*SIDE_LENGTH_Y or
           abs(Z0 - Z1) > BND_TOL*SIDE_LENGTH_Z or abs(Z0 - Z2) > BND_TOL*SIDE_LENGTH_Z or abs(Z0 - Z3) > BND_TOL*SIDE_LENGTH_Z or abs(Z1 - Z2) > BND_TOL*SIDE_LENGTH_Z or abs(Z1 - Z3) > BND_TOL*SIDE_LENGTH_Z or abs(Z2 - Z3) > BND_TOL*SIDE_LENGTH_Z){
                NEW_TRIANGLE.set_boundary(1);
                // std::cout << "Boundary " << std::endl;
        }else{
                NEW_TRIANGLE.set_boundary(0);
        } 
        
        NEW_TRIANGLE.setup_normals();

        return NEW_TRIANGLE;
}
#endif
