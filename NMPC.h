#ifndef PTL_NMPC_MODEL
#define PTL_NMPC_MODEL

#include "matrix.h"

/*
:--------------------------:
f(x,u,t) 状態方程式
L(x,u,t) ステージコスト
phi(x,t) 終端コスト
phix(x,t) dphi/dx
Hx(x,u,lmd,t) dH/dx
Hu(x,u,lmd,t) dH/du
:--------------------------:
偏微分いやなら
f,L,phiの実装
偏微分頑張るなら
f,phix,Hx,Huの実装
:--------------------------:
拘束条件がある場合 dim_u += dim_mu
f    row:dim_x col:1
u    row:dim_u col:1
x    row:dim_x col:1
lmd  row:1     col:dim_x
phi  row:1     col:1 (scalar)
L    row:1     col:1 (scalar)
phix row:1     col:dim_x
Hx   row:1     col:dim_x
Hu   row:1     col:dim_u
U    row:N     col:dim_u
X    row:N+1   col:dim_x
LMD  row:N     col:dim_x
F    row:N     col:dim_u
:--------------------------:
*/


using matrix = Matrix<double>;
class Model{
protected:
    uint dim_u,dim_x,dim_mu;
    const double h = 10e-6;
public:
    Model(uint input_dim,uint state_dim,uint constraint_dim = 0):dim_u(input_dim),dim_x(state_dim),dim_mu(constraint_dim){}
    uint dim_input() {return dim_u;}
    uint dim_state() {return dim_x;}
    uint dim_constraint() {return dim_mu;}
    virtual matrix f(matrix x,matrix u,double t = 0.0);
    virtual matrix L(matrix x,matrix u,double t = 0.0);
    virtual matrix phi(matrix x,double t = 0.0);
    virtual matrix phix(matrix x,double t = 0.0);
    virtual matrix Hx(matrix x,matrix u,matrix lmd,double t = 0.0);
    virtual matrix Hu(matrix x,matrix u,matrix lmd,double t = 0.0);
};

class NMPC{
private:
    Model* M;
    uint N;
    double dt;
    double zeta;
    uint dim_u,dim_x,dim_mu;
public:
    NMPC(Model* model,uint step,double dt,double zeta):M(model),N(step),dt(dt),zeta(zeta),dim_u(model->dim_input()),dim_x(model->dim_state()),dim_mu(model->dim_constraint()){}
    matrix X(matrix x,matrix U,double t = 0.0);
    matrix LMD(matrix X,matrix U,double t = 0.0);
    matrix F(matrix x,matrix U,double t = 0.0);
    matrix GMRES(matrix x,matrix U,double t = 0.0);
    matrix GRADIENT(matrix x,matrix U,double t = 0.0);
    matrix BFGS(matrix x,matrix U,double t = 0.0);
};

#endif /*PTL_NMPC_MODEL*/