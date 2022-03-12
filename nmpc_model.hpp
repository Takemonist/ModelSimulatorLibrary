
#ifndef NMPC_TWOWHEEL
#define NMPC_TWOWHEEL
#include "NMPC.cpp"

class TwoWheel:public Model{
public:
TwoWheel():Model(2,3,1){}
void addObj(double x,double y);
virtual matrix f(matrix x,matrix u,double t = 0.0);
private:
vector<matrix> obj; 
virtual matrix L(matrix x,matrix u,double t = 0.0);
virtual matrix phi(matrix x,double t = 0.0);
virtual matrix C(matrix x,matrix u,matrix dummy,double t = 0.0);

/*
matrix Hu(matrix x,matrix u,matrix lmd,double t = 0.0){
    matrix dHdu(1,dim_u+dim_mu);
    dHdu[0] = lmd[0]*cos(x[2])+lmd[1]*sin(x[2])+0.1*u[0];
    dHdu[1] = lmd[2]+2.0*u[1];
    return dHdu;
}

matrix Hx(matrix x,matrix u,matrix lmd,double t = 0.0){
    matrix dHdx(1,dim_x);
    dHdx[0] = 10.0*x[0];
    dHdx[1] = 10.0*x[1];
    dHdx[2] = -lmd[0]*u[0]*sin(x[2])+lmd[1]*u[0]*cos(x[2]);
    return dHdx;
}

matrix phix(matrix x,double t = 0.0){
    matrix dphidx(1,dim_x);
    dphidx[0] = 0.0*x[0];
    dphidx[1] = 0.0*x[1];
    dphidx[2] = 0.0;
    return dphidx;
}*/

};

void TwoWheel::addObj(double x,double y){
    matrix p(3,1);
    p[0] = x;
    p[1] = y;
    obj.push_back(p);
}
matrix TwoWheel::f(matrix x,matrix u,double t){
	x[0] = u[0]*cos(x[2]);
	x[1] = u[0]*sin(x[2]);
	x[2] = u[1];
    return x;
}
double lim(double x){
    if(x < 0.0) return 0;
    else return -log(x);
}
matrix TwoWheel::L(matrix x,matrix u,double t){
    matrix L(1,1);
    L = x.T()*diag({10.0,10.0,0.0})*x+u.T()*diag({0.0,1.0})*u;
    return L;
}
matrix TwoWheel::C(matrix x,matrix u,matrix dummy,double t){
    matrix min(1,1);
    matrix max(1,1);
    matrix val(1,1);
    min(0) = 10.0;
    max(0) = 1000.0;
    val(0) = u(0);
    return Constraints(min,max,val,dummy);
}
matrix TwoWheel::phi(matrix x,double t){
    return x.T()*diag({14.0,14.0,0.0})*x;
}

#endif