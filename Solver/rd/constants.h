#ifndef RD_constants_H
#define RD_constants_H

//-----------------------------------------
/* choose hydro test */
//----------------------------------------

//#define SODX
// #define SODY
// #define SINEX
//#define SEDOV
// #define GAUSSX
// #define GAUSSY
// #define UNIFORM
// #define NOH
// #define KHX
#define RTY
// #define KHY
// #define KHXSMOOTH
// #define KHYSMOOTH
// #define BLOB
// #define DF
//#define GRESHOVORTEX
//#define YEEVORTEX

//-----------------------------------------
/* set dimensionality */
//-----------------------------------------
#define TWO_D
// #define THREE_D

//-----------------------------------------
/* set umber of snapshots */
//-----------------------------------------
#define N_SNAP 80

//-----------------------------------------
/* debug flag for debug output */
//-----------------------------------------
//#define DEBUG

//-----------------------------------------
/* define flag to read positions and triangles from file (only option currently) */
//-----------------------------------------
#define READ_IC           // doesn't work yet


//-----------------------------------------
/* define flag for using either QHULL or CGAL triangulation */
//-----------------------------------------
// #define QHULL_IC
#define CGAL_IC

//-----------------------------------------
/* define boundar  y conditions (none for periodic) */
//-----------------------------------------
//#define CLOSED_BOUNDARY
#define PERIODIC_BOUNDARY
// #define FIXED_BOUNDARY
// #define REFLECTIVE_BOUNDARY        // doesn't work yet

//-----------------------------------------
/* define flag for timesteps */
//-----------------------------------------
//#define FIXED_DT
#define DRIFT
// #define DRIFT_SHELL
// #define JUMP
#include<cmath>
#include<string>

const double N_TBINS = 4; // set maximum time bin (must be power of 2)
const int MAX_TBIN = pow(2,N_TBINS);

//-----------------------------------------
/* define type of grid (none for square grid of vertices) */
//-----------------------------------------
// #define OFFSET_GRID
// #define EQUILATERAL_GRID

//-----------------------------------------
/* define distribution scheme */
//-----------------------------------------


//#define SU_SCHEME
#define LDA_SCHEME
//#define N_SCHEME
//#define BLENDED
//#define Bx_SCHEME
//-----------------------------------------
/* set order of scheme (none for 2nd order) */
//-----------------------------------------
//#define FIRST_ORDER

// #define SELF_GRAVITY // !!! NOT PERIODIC !!!
// #define ANALYTIC_GRAVITY
// #define PARA_RES
#define PARA_UP

//const double GRAV = 6.67e-11;
const double MSOLAR = 1.989e+30;

const double M_LIM = 0.0001;    // change for different tests
const double E_LIM = 0.0001;
const double C_LIM = 0.0001;

const std::string OUT_DIR = "Mesh_case/Case_1/LDA_Scheme/";
const std::string LOG_DIR = OUT_DIR + "log.txt";

// Sod Shock Tube (Varied in X)
#ifdef SODX
const double CFL = 0.5;
const double T_TOT = 0.2;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 2.0;
const double SIDE_LENGTH_Y = 2.0;
#endif

// Sod Shock Tube (Varied in Y)
#ifdef SODY
const double CFL = 0.1;
const double T_TOT = 0.1;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 2.0;
const double SIDE_LENGTH_Y = 2.0;
#endif

// Sine Wave Tube
#ifdef SINEX
const double CFL = 0.1;
const double T_TOT = 5.0;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 2.0;
const double SIDE_LENGTH_Y = 2.0;
#endif

// Sedov Blast Wave
#ifdef SEDOV
const double CFL = 0.2;
const double T_TOT = 0.01;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 2.0; // if altered, change setup.cpp as well
const double SIDE_LENGTH_Y = 2.0;

const double BLAST_E_TOT = 0.0;
const double R_BLAST     = 0.25;
#endif

// Gaussian pulse advection (x-direction)
#ifdef GAUSSX
const double CFL = 0.4;
const double T_TOT = 0.5;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 1.0;
const double SIDE_LENGTH_Y = 1.0;
#endif

// Gaussian pulse advection (y-direction)
#ifdef GAUSSY
const double CFL = 0.1;
const double T_TOT = 0.1;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 1.0;
const double SIDE_LENGTH_Y = 1.0;
#endif

// Uniform flow
#ifdef UNIFORM
const double CFL = 0.1;
const double T_TOT = 1.0;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 1.0;
const double SIDE_LENGTH_Y = 1.0;
#endif

// 2D Noh problem
#ifdef NOH
const double CFL = 0.1;
const double T_TOT = 1.0;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 1.0;
const double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability (x flow)
#ifdef KHX
const double CFL = 0.4;
const double T_TOT = 2.0;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 1.0;
const double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability (y flow)
#ifdef KHY
const double CFL = 0.4;
const double T_TOT = 4.0;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 1.0;
const double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability - smoothed (x flow)
#ifdef KHXSMOOTH
const double CFL = 0.4;
const double T_TOT = 4.0;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 1.0;
const double SIDE_LENGTH_Y = 1.0;
#endif

// KH instability - smoothed (y flow)
#ifdef KHYSMOOTH
const double CFL = 0.4;
const double T_TOT = 4.0;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 1.0;
const double SIDE_LENGTH_Y = 1.0;
#endif

// Blob test
#ifdef BLOB
const double CFL = 0.4;
const double T_TOT = 4.0;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 10.0;
const double SIDE_LENGTH_Y = 10.0;
#endif

// Grav Test
#ifdef DF
const double CFL = 0.4;
const double T_TOT = 120.0;
const double GAMMA = 5.0/3.0;
const double SIDE_LENGTH_X = 10.0;
const double SIDE_LENGTH_Y = 10.0;
const double MACH = 0.0;
#endif



#ifdef GRESHOVORTEX
const double CFL=0.9;
const double T_TOT=3;
const double GAMMA=5.0/3.0;
const double SIDE_LENGTH_X = 1.0;
const double SIDE_LENGTH_Y = 1.0;
const double R_INNER_VORTEX=0.2;
const double R_OUTER_VORTEX=0.4;


#endif 

#ifdef YEEVORTEX
const double CFL=0.9;
const double T_TOT=10;
const double GAMMA=5.0/3.0;
const double SIDE_LENGTH_X = 10.0;
const double SIDE_LENGTH_Y = 10.0;
#endif


#ifdef RTY
const double CFL=0.9;
const double T_TOT=2;
const double GAMMA=5.0/3.0;
const double SIDE_LENGTH_X = 0.25;
const double SIDE_LENGTH_Y = 2.0;
const double Y_HALF=0.5*SIDE_LENGTH_Y;
const double GRAV=1;
#endif


const double GAMMA_1 = GAMMA - 1.0;
const double GAMMA_2 = GAMMA - 2.0;

#ifdef FIXED_DT
const double DT_FIX = 0.000001;
#endif


#endif //RD_constants_H
