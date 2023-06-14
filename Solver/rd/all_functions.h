
#include <vector>
#include <iostream>
#include <fstream>
#include "vertex2D.h"
#include "triangle2D.h"
//source2D.cpp
void sources(std::vector<VERTEX>&, double, int);
void plummer_gravity(std::vector<VERTEX>&, double, int);
void direct_gravity(std::vector<VERTEX>&, double, int);

//timestep.cpp
void drift_update_half(int, int, double, double, std::vector<TRIANGLE> &);
void jump_update_half(int, int, double, double, std::vector<TRIANGLE> &);
void drift_update(int, int, double, double, std::vector<TRIANGLE> &);
void reset_tbins(double, double, int, int, double &, std::vector<TRIANGLE>&, std::vector<VERTEX> &);

//io2D.cpp
void open_snap(std::ofstream &, int);
void open_active(std::ofstream &, int);
void write_snap(std::vector<VERTEX>, double, double, int, int, std::ofstream &);
void write_active(std::vector<TRIANGLE>, int, int, int);
void read_parameter_file(int ARGC, char *ARGV[]);
int cgal_read_positions_header(std::ifstream &);
int cgal_read_triangles_header(std::ifstream &);
VERTEX cgal_read_positions_line(std::ifstream &);
TRIANGLE cgal_read_triangles_line(std::ifstream &, std::vector<VERTEX> &, int);

//setup2D.cpp
double F(double);
double G(double);
VERTEX setup_vertex(double, double);


