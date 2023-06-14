/*
 * This file was written by Ben Morton (bmorton@ed.ac.uk) and Zhenyu Wu (zhenyu.wu@ed.ac.uk).
 */
 
#include <cblas.h>
#include <lapacke.h>
#include <iostream>

lapack_int mat_inv(double *A, unsigned n, double X, double Y, int ID, int POINT){
        int ipiv[n+1];
        lapack_int ret;

        ret =  LAPACKE_dgetrf(LAPACK_COL_MAJOR,n,n,A,n,ipiv);

#ifdef DEBUG
        std::cout << "ret =\t" << ret << "\t(0 = done, <0 = illegal arguement, >0 = singular)" << std::endl;
#endif

        if(ret !=0){
                std::cout << "B WARNING: MATRIX CANNOT BE INVERTED\t" << ret << "\t(0 = done, <0 = illegal arguement, >0 = singular) at " << X << "\t" << Y << "\tPoint =\t" << POINT << "\t" << ID << std::endl;
                for(int i=0;i<n*n;++i){
                        std::cout << A[i] << "\t";
                        if((i+1)%n == 0){std::cout << std::endl;}
                }
                exit(0);
        }
        ret = LAPACKE_dgetri(LAPACK_COL_MAJOR,n,A,n,ipiv);

        return ret;
}

lapack_int matFac(double *A, unsigned n){
        int ipiv[n+1];
        lapack_int ret;

        ret =  LAPACKE_dgetrf(LAPACK_COL_MAJOR,n,n,A,n,ipiv);
#ifdef DEBUG
        std::cout << "ret =\t" << ret << "\t(0 = done, <0 = illegal arguement, >0 = singular)" << std::endl;
#endif
        return ret;
}
