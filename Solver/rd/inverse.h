//
// Created by Wu Zhenyu on 04/01/2022.
//

#ifndef RD_INVERSE_H
#define RD_INVERSE_H

#include <cblas.h>
#include <lapacke.h>
//lapack_int mat_inv(double*, unsigned, double, double, int, int);

lapack_int mat_inv(double *A, unsigned n, double X, double Y, int ID, int POINT)
{
    int ipiv[n+1];
    lapack_int ret;
    /*
	for(int i=0;i<n*n;++i){
                        std::cout << A[i] << "\t";
                        if((i+1)%n == 0){std::cout << std::endl;}
                }
	*/
    ret =  LAPACKE_dgetrf(LAPACK_COL_MAJOR,
                          n,
                          n,
                          A,
                          n,
                          ipiv);

    if(ret !=0){              
         return ret;
        }
        


    ret = LAPACKE_dgetri(LAPACK_COL_MAJOR,
                       n,
                       A,
                       n,
                       ipiv);
    return ret;
}

lapack_int matFac(double*, unsigned);

#endif //RD_INVERSE_H
