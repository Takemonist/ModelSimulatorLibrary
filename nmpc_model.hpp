
#ifndef NMPC_TWOWHEEL
#define NMPC_TWOWHEEL
#include "NMPC.cpp"

class TwoWheel:public Model{
public:
TwoWheel():Model(2,3){}
void addObj(double x,double y);
virtual matrix f(matrix x,matrix u,double t = 0.0);
private:
vector<matrix> obj; 
virtual matrix L(matrix x,matrix u,double t = 0.0);
virtual matrix phi(matrix x,double t = 0.0);
/*
virtual matrix Hu(matrix x,matrix u,matrix lmd,double t = 0.0);
virtual matrix Hx(matrix x,matrix u,matrix lmd,double t = 0.0);
virtual matrix phix(matrix x,double t = 0.0);*/
};

void TwoWheel::addObj(double x,double y){
    matrix p(3,1);
    p[0] = x;
    p[1] = y;
    obj.push_back(p);
}
matrix TwoWheel::f(matrix x,matrix u,double t){
	x.at(0,0) = u.at(0,0)*cos(x.at(2,0));
	x.at(1,0) = u.at(0,0)*sin(x.at(2,0));
	x.at(2,0) = u.at(1,0);
    return x;
}
double lim(double x){
    if(x < 0.0) return 0;
    else return -log(x);
}
matrix TwoWheel::L(matrix x,matrix u,double t){
    matrix L(1,1);
    L = x.T()*diag({2.0,2.0,0.0})*x+u.T()*diag({0.1,0.5})*u;
    //for(int i = 0;i < obj.size();i++) L += lim(0.2-fabs(obj[i][0]-x[0])) + lim(0.2-fabs(obj[i][1]-x[1]));
    return L;
}
matrix TwoWheel::phi(matrix x,double t){
    return x.T()*diag({4.0,4.0,0.0})*x;
}
/*
matrix TwoWheel::Hu(matrix x,matrix u,matrix lmd,double t){
    matrix dHdu(1,dim_u+dim_mu);
    dHdu[0] = lmd[0]*cos(x[2])+lmd[1]*sin(x[2])+0.2*u[0];
    dHdu[1] = lmd[2]-1.0/(u[1]-0.2)+1.0*u[1];
    return dHdu;
}

matrix TwoWheel::Hx(matrix x,matrix u,matrix lmd,double t){
    matrix dHdx(1,dim_x);
    dHdx[0] = 4.0*x[0];
    dHdx[1] = 4.0*x[1];
    dHdx[2] = -lmd[0]*u[0]*sin(x[2])+lmd[1]*u[0]*cos(x[2]);
    return dHdx;
}

matrix TwoWheel::phix(matrix x,double t){
    matrix dphidx(1,dim_x);
    dphidx[0] = 8.0*x[0];
    dphidx[1] = 8.0*x[1];
    dphidx[2] = 0.0;
    return dphidx;
}*/
#endif