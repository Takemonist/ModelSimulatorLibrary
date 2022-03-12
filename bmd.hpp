
#ifndef NMPC_BMD
#define NMPC_BMD
#include "NMPC.cpp"

class BMD:public Model{
public:
BMD():Model(1,2,1){}
virtual matrix f(matrix x,matrix u,double t = 0.0);
public:
virtual matrix L(matrix x,matrix u,double t = 0.0);
matrix Huu(matrix x,matrix u,double t = 0.0);
virtual matrix phi(matrix x,double t = 0.0);
virtual matrix C(matrix x,matrix u,double t = 0.0);
/*
matrix Hu(matrix x,matrix u,matrix lmd,double t = 0.0){
    matrix dHdu(1,dim_u+dim_c*2);
    dHdu[0] = 1*u[0] + lmd[1]*x[1] + u[2]*(2*u[0] - 0.5);
    dHdu[1] = 2*u[1]*u[2] - 0.01;
    dHdu[2] = C(x,u,t)[0];
    return dHdu;
}

matrix Hx(matrix x,matrix u,matrix lmd,double t = 0.0){
    matrix dHdx(1,dim_x);
    dHdx[0] = -1.0*lmd[1] + 1.0*x[0];
    dHdx[1] = -1.0*lmd[1]*u[0] + 10.0*x[1]+lmd[0];
    return dHdx;
}

matrix phix(matrix x,double t = 0.0){
    matrix dphidx(1,dim_x);
    dphidx[0] = 1.0*x[0];
    dphidx[1] = 10.0*x[1];
    return dphidx;
}*/

};

matrix BMD::f(matrix x,matrix u,double t){
    matrix f(dim_x,1);
	f[0] = x[1];
	f[1] = -1.0*x[0]-1.0*x[1]*u[0];
    return f;
}
matrix BMD::L(matrix x,matrix u,double t){
    matrix L(1,1);
    L = x.T()*diag({1.0,10.0})*x+u.T()*diag({0.01,0.00})*u-0.02*u[1];
    return L/2.0;
}
matrix BMD::C(matrix x,matrix u,double t){
    matrix min(1,1);
    matrix max(1,1);
    matrix val(1,1);
    min(0) = -1.0;
    max(0) = 1.0;
    val(0) = u(0);
    return Constraints(min,max,val,u.block(0,0,dim_u,dim_u+dim_c-1));
}
matrix BMD::phi(matrix x,double t){
    return (x.T()*diag({1.0,10.0})*x)/2.0;
}
matrix BMD::Huu(matrix x,matrix u,double t){
    matrix Huu(1,1);
    return Huu;
}

#endif