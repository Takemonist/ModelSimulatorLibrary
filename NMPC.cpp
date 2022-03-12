#ifndef PTL_NMPC_CPP
#define PTL_NMPC_CPP

#include "NMPC.h"

double circle(double x,double y,double min,double max){
    double r = (max-min)/2.0;
    double d = (max+min)/2.0;
    return (x-d)*(x-d) + y*y - r*r;
}

matrix Model::Constraints(matrix min,matrix max,matrix val,matrix dummy){
    matrix C(dim_c,1);
    for(int i = 0;i < dim_c;i++){
        C(i) = circle(val(i),dummy(i),min(i),max(i));
    }
    return C;
}

matrix Model::f(matrix x,matrix u,double t){
    return x;
}
matrix Model::L(matrix x,matrix u,double t){
    return x.T()*x;
}
matrix Model::C(matrix x,matrix u,double t){
    matrix o(0,0);
    return o;
}
double Model::H(matrix x,matrix u,matrix lmd,double t){
    matrix mu; 
    mu = u.block(0,0,dim_u+dim_c,-1);
    u = u.block(0,0,0,dim_u+dim_c-1).T();
    x = x.T();
    double H = 0.0;
    if(dim_c > 0) H = ( L(x,u,t) + lmd*f(x,u,t) + mu*C(x,u,t) ).at(0);
    else H = ( L(x,u,t) + lmd*f(x,u,t) ).at(0);
    return H;
}
matrix Model::phi(matrix x,double t){
    return x.T()*x;
}
matrix Model::phix(matrix x,double t){
    matrix dphidx(1,dim_x);
    double phi0 = phi(x,t).at(0);
    for(uint i = 0;i < dim_x;i++){
        x.at(i) += h;
        dphidx.at(i) = phi(x,t).at(0);
        x.at(i) -= h;
    }
    dphidx -= phi0;
    dphidx /= h;
    return dphidx;
}
matrix Model::Hx(matrix x,matrix u,matrix lmd,double t){
    matrix dHdx(1,dim_x);
    double H_ = H(x,u,lmd,t);
    for(uint i = 0;i < dim_x;i++){
        x.at(i) += h;
        dHdx.at(i) = H(x,u,lmd,t);
        x.at(i) -= h;
    }
    dHdx -= H_;
    dHdx /= h;
    return dHdx;
}
matrix Model::Hu(matrix x,matrix u,matrix lmd,double t){
    matrix dHdu(1,dim_u+dim_c*2);
    double H_ = H(x,u,lmd,t);
    for(uint i = 0;i < dim_u+dim_c;i++){
        u.at(i) += h;
        dHdu.at(i) = H(x,u,lmd,t) - H_;
        dHdu.at(i) /= h;
        u.at(i) -= h;
    }
    matrix C_ = C(x,u,t);
    for(uint i = 0;i < dim_c;i++){
        dHdu.at(dim_u+dim_c+i) = C_(i);
    }
    return dHdu;
}

matrix NMPC::X(matrix x,matrix U,double t){
    matrix X(N+1,dim_x);
    matrix u;
    X.row(0) = x.T();
    for(unsigned i = 0;i < N;i++){
        u = U.block(i,i,0,dim_u+dim_c-1).T();
        x += M->f(x,u,t)*dt;
        X.row(i+1) = x.T();
    }
    return X;
}
matrix NMPC::LMD(matrix X,matrix U,double t){
    matrix LMD(N,dim_x);
    matrix lmd(1,dim_x);
    matrix u(1,dim_U);
    matrix x(1,dim_x);
    x = X.row(N).T();
    lmd = M->phix(x,t);
    LMD.row(N-1) = lmd;
    for(uint i = N-1;i > 0;i--){
        u = U.row(i);
        x = X.row(i);
        lmd += M->Hx(x,u,lmd,t)*dt;
        LMD.row(i-1) = lmd;
    }
    return LMD;
}
matrix NMPC::F(matrix x,matrix U,double t){
    matrix X(N+1,dim_x);
    matrix LMD(N,dim_x);
    matrix F(N,dim_U);
    matrix dHdu(1,dim_U);
    X = this->X(x,U,t);
    LMD = this->LMD(X,U,t);
    for(uint i = 0;i < N;i++){
        dHdu = M->Hu(X.row(i),U.row(i),LMD.row(i),t);
        F.row(i) = dHdu;
    }
    return F.reshape(N*dim_U,1);
}
matrix NMPC::FU(matrix x,matrix U,double t){
    matrix dFdU(N*(dim_u+dim_c),N*(dim_u+dim_c));
    matrix f = F(x,U,t);
    for(int i = 0;i < N*(dim_u+dim_c);i++){
        U[i] += 10e7;
        dFdU.col(i) = (F(x,U,t) - f)/10e7;
        U[i] -= 10e7;
    }
    return dFdU;
}
matrix NMPC::GMRES(matrix x,matrix U,double t){
    int m = (dim_U)*N;
    matrix dx = M->f(x,U.block(0,0,0,dim_u-1).T(),t)*dt;
    matrix Fxt = F(x+dx,U,t+dt);
    matrix F = this->F(x,U,t);
    matrix dU = U;
    matrix Fuxt = this->F(x+dx,U+dU*dt,t+dt);
    matrix Fut = this->F(x,U+dU*dt,t+dt);
    matrix right = -zeta*F - ((Fxt - F) / dt);
    matrix left = ((Fuxt - Fxt) / dt);
    
    matrix r0 = right - left;
    double r0_norm = norm(r0);
	matrix V(r0.rows(),m);
	matrix g(m+1,1);
	matrix H(m+1,m);
    V.col(0) = r0/r0_norm;
    g(0) = r0_norm;
    matrix v;
    matrix AVj;
    for(int j = 0;j < m;j++){
        v = V.col(j);
        v.reshape(N,dim_U);
        AVj = ((this->F(x+dx,U+v*dt,t+dt) - Fxt) / dt);
        for(int i = 0;i <= j;i++) H(i,j) = dot(AVj,V.col(i));
        for(int i = 0;i <= j;i++) AVj -= H(i,j)*V.col(i);
        H(j+1,j) = norm(AVj);
        if(j < m-1) V.col(j+1) = AVj/H(j+1,j);
    }
    matrix y = H.pinv()*g;
    dU += (V.block(0,-1,0,m-1)*y).reshape(N,dim_U);
    U += dU*dt;
	printf("r0:%lf F:%lf\n",norm(g-H*y),norm(this->F(x,U,t)));

    return U;
}
matrix NMPC::GRADIENT(matrix x,matrix U,double t){
    double a = 0.01;
    for(int i = 0;norm(F(x,U)) > 0.1 and i < 1000;i++){
        U += -a*(F(x,U)).reshape(N,(dim_U));
    }
    printf("F:%lf\n",norm(F(x,U,t)));
    return U;
}

matrix NMPC::BFGS(matrix x,matrix U,double t){
    return U;
}

#endif