/*                                                                           
 * This file was written by Ben Morton (bmorton@ed.ac.uk) and Zhenyu Wu (zhenyu.wu@ed.ac.uk).
 * IO routines to write simple ASCII output for python plotting
        open_files => opens ouptut files to write positions, dens, pressure, vel maps, and column of values for 1D plot
*/

#include <iostream>
#include <fstream>
#include "vertex2D.h"
#include "triangle2D.h"
//#include "all_functions.h"

void open_snap(std::ofstream &SNAPFILE, int i){
        SNAPFILE.open(OUT_DIR+"snapshot_"+std::to_string(i)+".txt");
        return;
}

void open_active(std::ofstream &SNAPFILE, int i){
        SNAPFILE.open(OUT_DIR+"active_"+std::to_string(i)+".txt");
        return;
}

void open_sc_file(std::ofstream &SNAPFILE, int i){
        SNAPFILE.open(OUT_DIR+"sc_file_"+std::to_string(i)+".txt");
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
                SNAPFILE << POINTS[i].get_x() << "\t" << POINTS[i].get_y() << "\t" << POINTS[i].get_mass_density() << "\t" << POINTS[i].get_x_velocity() << "\t" << POINTS[i].get_y_velocity() << "\t" << POINTS[i].get_pressure() << "\t" << POINTS[i].get_specific_energy() << "\t" << POINTS[i].get_dual() << std::endl;
                TOTAL_DENSITY += POINTS[i].get_mass_density()*POINTS[i].get_dual();
                TOTAL_ENERGY += POINTS[i].get_specific_energy()*POINTS[i].get_dual() * POINTS[i].get_mass_density();
        }
        std::cout << "*************************************************************************************************" << std::endl;            // right out time and total density to terminal
        std::cout << "time\t" << T << " \t-> total mass =\t" << TOTAL_DENSITY << " \t-> total energy =\t" << TOTAL_ENERGY << "\ttime step = \t" << DT << std::endl;
        LOGFILE << T << "\t" << TOTAL_DENSITY << "\t" << TOTAL_ENERGY << "\t" << DT << "\t" << MAX_TBIN << "\t" << N_POINTS << std::endl;
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


void write_sc(std::vector<TRIANGLE> MESH, int N_TRIANG, int SNAP_ID, double T){
	std::ofstream SNAPFILE;
        SNAPFILE << std::setprecision(12);
        double X0,X1,X2,Y0,Y1,Y2,SC;
        open_sc_file(SNAPFILE, SNAP_ID);
        SNAPFILE << N_TRIANG << "\t" << T << "\t" <<  std::endl;
         for(int j=0;j<N_TRIANG;++j){
                if(MESH[j].get_boundary() == 0){
                        X0 = MESH[j].get_vertex_0()->get_x();
                        X1 = MESH[j].get_vertex_1()->get_x();
                        X2 = MESH[j].get_vertex_2()->get_x();
                        Y0 = MESH[j].get_vertex_0()->get_y();
                        Y1 = MESH[j].get_vertex_1()->get_y();
                        Y2 = MESH[j].get_vertex_2()->get_y();
                        SC = MESH[j].get_sc();
                        // write         X        Y          TBIN
                        SNAPFILE << X0 << "\t" << Y0 << "\t"  << X1 << "\t" << Y1 << "\t"  << X2 << "\t" << Y2 << "\t" << SC << "\t" << std::endl;
		}
	}
}


void read_parameter_file(int ARGC, char *ARGV[]){
        printf("Parameter file = %s\n", ARGV[1]);
}

// if using qhull triangulation (closed boundaries only), read vertex header info on triangulation
#ifdef QHULL_IC
int qhull_read_positions_header(std::ifstream &POSITIONS_FILE){
        std::string INFO;
        int N_POINTS;

        std::getline(POSITIONS_FILE,INFO);

        std::cout << "Details of point creation =\t" << INFO << std::endl;

        POSITIONS_FILE >> N_POINTS;

        return N_POINTS;
}

// read qhull triangles header info
int qhull_read_triangles_header(std::ifstream &TRIANGLES_FILE){
        int N_TRIANG;

        TRIANGLES_FILE >> N_TRIANG;

        return N_TRIANG;
}

// read one qhull vertex position
VERTEX qhull_read_positions_line(std::ifstream &POSITIONS_FILE){
        double X,Y;
        VERTEX NEW_VERTEX;

        POSITIONS_FILE >> X >> Y;

        X = SIDE_LENGTH_X*(X + 0.5);
        Y = SIDE_LENGTH_Y*(Y + 0.5);

        NEW_VERTEX = setup_vertex(X,Y);

        return NEW_VERTEX;
}

// read one qhull triangle (indices of vertices)
TRIANGLE qhull_read_triangles_line(std::ifstream &TRIANGLES_FILE, std::vector<VERTEX> &POINTS){
        int N_VERT,VERT0,VERT1,VERT2;
        TRIANGLE NEW_TRIANGLE;

        TRIANGLES_FILE >> N_VERT >> VERT0 >> VERT1 >> VERT2;

        // std::cout << VERT0 << "\t" << VERT1 << "\t" << VERT2 << std::endl;

        NEW_TRIANGLE.set_vertex_0(&POINTS[VERT0]);
        NEW_TRIANGLE.set_vertex_1(&POINTS[VERT1]);
        NEW_TRIANGLE.set_vertex_2(&POINTS[VERT2]);

        NEW_TRIANGLE.setup_normals();

        return NEW_TRIANGLE;
}
#endif

// if using CGAL triangulation file, read vertex header info
#define CGAL_IC
#ifdef CGAL_IC
int cgal_read_positions_header(std::ifstream &CGAL_FILE){
        int N_POINTS, XSHEETS, YSHEETS;
        float XLOW, YLOW, XHIGH, YHIGH;

        CGAL_FILE >> XLOW >> YLOW >> XHIGH >> YHIGH;

        // check boundaries match
        if(XLOW < 0.0 or YLOW < 0.0 or int(XHIGH) != SIDE_LENGTH_X or int(YHIGH) != SIDE_LENGTH_Y){
                std::cout << "BWARNING: Exiting on mismatching boundaries. Have (X): " << XLOW  << "\t" << XHIGH << "\tNeed (X): " << 0.0 << "\t" << SIDE_LENGTH_X << std::endl;
                exit(0);
        }

        CGAL_FILE >> XSHEETS >> YSHEETS;

        // check data in 1-sheet format
        if(XSHEETS !=1 or YSHEETS!=1){
                std::cout << "BWARNING: Exiting on CGAL file not in 1-sheet format." << std::endl;
                exit(0);
        }

        CGAL_FILE >> N_POINTS;

        return N_POINTS;
}

// read CGAL triangle header info
int cgal_read_triangles_header(std::ifstream &CGAL_FILE){
        int N_TRIANG;
        int EMPTY;

        CGAL_FILE >> N_TRIANG;

        return N_TRIANG;
}

// read postion of one CGAL vertex
VERTEX cgal_read_positions_line(std::ifstream &CGAL_FILE){
        double X,Y;
        VERTEX NEW_VERTEX;

        CGAL_FILE >> X >> Y;

        NEW_VERTEX = setup_vertex(X,Y);

        return NEW_VERTEX;
}

// read indices of vertices for one CGAL triangle
TRIANGLE cgal_read_triangles_line(std::ifstream &CGAL_FILE, std::vector<VERTEX> &POINTS, int ID){
        int VERT0,VERT1,VERT2;
        TRIANGLE NEW_TRIANGLE;

        CGAL_FILE >> VERT0 >> VERT1 >> VERT2;

        NEW_TRIANGLE.set_vertex_0(&POINTS[VERT0]);
        NEW_TRIANGLE.set_vertex_1(&POINTS[VERT1]);
        NEW_TRIANGLE.set_vertex_2(&POINTS[VERT2]);

        NEW_TRIANGLE.set_id(ID);

        // std::cout << POINTS[VERT0].get_x() << "\t" << POINTS[VERT1].get_x() << "\t" << POINTS[VERT2].get_x() << std::endl;

        double X0,X1,X2,Y0,Y1,Y2;

        // check if boundary triangle

        X0 = POINTS[VERT0].get_x();
        X1 = POINTS[VERT1].get_x();
        X2 = POINTS[VERT2].get_x();

        Y0 = POINTS[VERT0].get_y();
        Y1 = POINTS[VERT1].get_y();
        Y2 = POINTS[VERT2].get_y();

        if(abs(X0 - X1) > 0.5*SIDE_LENGTH_X or abs(X0 - X2) > 0.5*SIDE_LENGTH_X or abs(X1 - X2) > 0.5*SIDE_LENGTH_X or
           abs(Y0 - Y1) > 0.5*SIDE_LENGTH_Y or abs(Y0 - Y2) > 0.5*SIDE_LENGTH_Y or abs(Y1 - Y2) > 0.5*SIDE_LENGTH_Y){
                NEW_TRIANGLE.set_boundary(1);
        }else{
                NEW_TRIANGLE.set_boundary(0);
        } 
        
        NEW_TRIANGLE.setup_normals();

        return NEW_TRIANGLE;
}



#endif



