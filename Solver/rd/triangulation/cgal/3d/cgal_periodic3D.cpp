#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/periodic_3_triangulation_3_io.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>       Gt;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<Gt>             P3DT3;
typedef P3DT3::Point             Point;
typedef P3DT3::Iso_cuboid        Iso_cuboid;
// typedef P3DT3::Vertex_handle     Vertex_handle;
// typedef P3DT3::Cell_handle       Cell_handle;
// typedef P3DT3::Locate_type       Locate_type;


#define RANDOMIC
// #define UNIFORMIC
// #define UNIFORMOFFSETIC
// #define PERTUNIFORMOFFSETIC

int main(){
        float xmax=10.0,ymax=10.0,zmax=10.0;
        Iso_cuboid domain(0, 0, 0, xmax, ymax, zmax); // The cube for the periodic domain

        // construction from a list of points :
        std::list<Point> L;

        int i,j,k;
        int nx=256, ny=256, nz=256, count=nx*ny*nz;
        float x,y,z,xmove,ymove,zmove;

#ifdef RANDOMIC
        for (int i = 0; i < count; ++i){
	  if(int(100.0 * float(i) / float(count)) % 2 == 0){printf("%f\n",(float(i) / float(count)));}
                x = xmax*(rand() % 10000)/10000.0;
                y = ymax*(rand() % 10000)/10000.0;
                z = zmax*(rand() % 10000)/10000.0;
                L.push_back(Point(x,y,z));
        }
#endif

#ifdef UNIFORMIC
        for(i=0; i < (nx); ++i){
                x = (xmax) * (float(i)+0.5) / float(nx);
                for(j=0; j < (ny); ++j){
                        y = (ymax) * (float(j)+0.5) / float(ny);
                        for(k=0; k < (nz); ++k){
                                z = (zmax) * (float(k)+0.5) / float(nz);
                                L.push_back(Point(x,y,z));
                        }
                }
        }
#endif

#ifdef UNIFORMOFFSETIC        // check it works !!!!
        for(i=0; i < (nx); ++i){
                x = (xmax) * (float(i)+0.5) / float(nx);
                for(j=0; j < (ny); ++j){
                        y = (ymax) * (float(j)+0.5) / float(ny);
                        for(k=0; k < (nz); ++k){
                                if(k % 2 == 0){
                                        z = (zmax) * (float(k)+0.5) / float(nz);
                                }else{
                                        z = (zmax) * (float(k)+1.0) / float(nz);
                                }
                                L.push_back(Point(x,y,z));
                        }
                }
        }
#endif

        P3DT3 T(L.begin(), L.end(), domain); // Put the domain with the constructor

        assert( T.is_valid() ); // checking validity of T

        // std::cout << T.number_of_sheets()[0] << "\t" << T.number_of_sheets()[1] << "\t" << T.number_of_sheets()[2] << std::endl;

        T.convert_to_1_sheeted_covering();

        // std::cout << T.number_of_sheets()[0] << "\t" << T.number_of_sheets()[1] << "\t" << T.number_of_sheets()[2] << std::endl;


        std::ofstream oFileT("../../../Delaunay3D.txt", std::ios::out);
        // writing file output;
        oFileT << T;

        return 0;
}
