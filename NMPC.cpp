#ifndef PTL_NMPC_CPP
#define PTL_NMPC_CPP

#include "NMPC.h"

matrix Model::f(matrix x,matrix u,double t){
    return x;
}
matrix Model::L(matrix x,matrix u,double t){
    return x.T()*x;
}
matrix Model::phi(matrix x,double t){
    return x.T()*x;
}
matrix Model::phix(matrix x,double t){
    matrix dphidx(1,dim_x);
    double phi0 = phi(x,t).at(0);
    for(uint i = 0;i < dim_x;i++){
        x.at(i) += h;
        dphidx.at(i) = phi(x).at(0);
        x.at(i) -= h;
    }
    dphidx -= phi0;
    dphidx /= h;
    return dphidx;
}
matrix Model::Hx(matrix x,matrix u,matrix lmd,double t){
    matrix dHdx(1,dim_x);
    double H = (L(x,u,t)+lmd*f(x,u,t)).at(0);
    for(uint i = 0;i < dim_x;i++){
        x.at(i) += h;
        dHdx.at(i) = (L(x,u,t)+lmd*f(x,u,t)).at(0);
        x.at(i) -= h;
    }
    dHdx -= H;
    dHdx /= h;
    return dHdx;
}
matrix Model::Hu(matrix x,matrix u,matrix lmd,double t){
    matrix dHdu(1,dim_u+dim_mu);
    double H = (L(x,u,t)+lmd*f(x,u,t)).at(0);
    for(uint i = 0;i < dim_u+dim_mu;i++){
        u.at(i) += h;
        dHdu.at(i) = (L(x,u,t)+lmd*f(x,u,t)).at(0);
        u.at(i) -= h;
    }
    dHdu -= H;
    dHdu /= h;
    return dHdu;
}

matrix NMPC::X(matrix x,matrix U,double t){
    matrix X(N+1,dim_x);
    matrix u(dim_u+dim_mu,1);
    X.row(0) = x.T();
    for(unsigned i = 0;i < N;i++){
        u = U.row(i).T();
        x += M->f(x,u,t)*dt;
        X.row(i+1) = x.T();
    }
    return X;
}
matrix NMPC::LMD(matrix X,matrix U,double t){
    matrix LMD(N,dim_x);
    matrix lmd(dim_x,1);
    matrix u(dim_u+dim_mu,1);
    matrix x(dim_x,1);
    x = X.row(N).T();
    lmd = M->phix(x,t);
    LMD.row(N-1) = lmd;
    for(uint i = N-1;i > 0;i--){
        u = U.row(i).T();
        x = X.row(i).T();
        lmd += M->Hx(x,u,lmd,t)*dt;
        LMD.row(i-1) = lmd;
    }
    return LMD;
}
matrix NMPC::F(matrix x,matrix U,double t){
    matrix X(N+1,dim_x);
    matrix LMD(N,dim_x);
    matrix F(N,dim_u+dim_mu);
    matrix dHdu(dim_u+dim_mu,1);
    X = this->X(x,U,t);
    LMD = this->LMD(X,U,t);
    for(uint i = 0;i < N;i++){
        dHdu = M->Hu(X.row(i).T(),U.row(i).T(),LMD.row(i));
        F.row(i) = dHdu;
    }
    return F;
}
matrix NMPC::GMRES(matrix x,matrix U,double t){
    uint K = (dim_u+dim_mu)*N;
    matrix dx = M->f(x,U.row(0).T(),t)*dt;
    matrix Fxt = F(x+dx,U,t+dt);
    matrix F = this->F(x,U,t);
    matrix Fuxt;
    matrix dU = U;
    matrix right = -zeta*F - ((Fxt - F) / dt);
    matrix left;
    
    Fuxt = this->F(x+dx,U+dU*dt,t+dt);
    left = ((Fuxt - Fxt) / dt);
    right.reshape(K,1);
    left.reshape(K,1);
    matrix r0 = right - left;
    double r0_norm = norm(r0);
    matrix e1(K+1,1);e1.at(0) = 1.0; 
    matrix V(K,K);
    V.col(0) = r0/(r0_norm+10e-7);
    matrix H(K+1,K);
    matrix Av;
    matrix y;
    matrix Hm;
    for(int i = 0;i < K;i++){
        matrix v = V.col(i)*dt;
        v.reshape(N,dim_u+dim_mu);
        Fuxt = this->F(x+dx,U+v);
        Av = ((Fuxt - Fxt) / dt);
        Av.reshape(K,1);
        for(int j = 0;j <= i;j++) H(j,i) = dot(Av,V.col(j));
        matrix sum_Av(K,1);
        for(int j = 0;j <= i;j++) sum_Av += V.col(j)*H(j,i);
        matrix v_hat = Av - sum_Av;
        H(i+1,i) = norm(v_hat);
        V.col(i+1) = v_hat/(H(i+1,i)+10e-7);
        Hm = H.block(0,i+1,0,i);
        y = Hm.pinv()*(e1.block(0,i+1,0,0)*r0_norm);
        if(norm(e1.block(0,i+1,0,0)*r0_norm-Hm*y) < 0.001 || i == K-1){
            matrix update_value = V.block(0,-1,0,i)*y.block(0,i,0,0);
            update_value.reshape(N,dim_u+dim_mu);
            dU += update_value;
            U += dU*dt;
			printf("m:%d F:%lf\n",i+1,norm(this->F(x,U,t).reshape(K,1)));
            break;
        }
    }
    return U;
}
matrix NMPC::GRADIENT(matrix x,matrix U,double t){
    for(int i = 0;i < 10;i++){
        U += -0.05*F(x,U);
    }
    //printf("F:%lf\n",norm(F(x,U,t).reshape(N*(dim_u+dim_mu),1)));
    return U;
}

#endif