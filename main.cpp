#include <unistd.h>
//#include "matrix.h"
#include "nmpc_model.hpp"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
typedef Matrix<double> mat;

mat state_func(mat x,mat u){
	x.at(0,0) = u.at(0,0)*cos(x.at(2,0));
	x.at(1,0) = u.at(0,0)*sin(x.at(2,0));
	x.at(2,0) = u.at(1,0);
	return x;
}
mat state_cost(mat x,mat u){
	return x.T()*diag({0.01,0.01,0.0})*x;
}
mat input_cost(mat x,mat u){
	return u.T()*diag({0.2,0.0})*u;
}
mat terminal_cost(mat x){
	return x.T()*diag({14.0,14.0,0.0})*x;
}


void GMRES(mat A,mat b,uint N){
	int K = N;
	mat x(N,1);
	for(int k = 0;k < 1;k++){
		mat r0 = b - A*x;
		double r0_norm = norm(r0);
		mat e1(K+1,1);e1.at(0) = 1.0; 
		mat V(N,K);
		V.col(0) = r0/r0_norm;
		mat H(K+1,K);
		mat Av;
		mat y;
		mat h;
		for(int i = 0;i < K;i++){
			Av = A*V.col(i);
			for(int j = 0;j <= i;j++) H(j,i) = dot(Av,V.col(j));
			mat sum_Av(K,1);
			for(int j = 0;j <= i;j++) sum_Av += V.col(j)*H(j,i);
			mat v = Av - sum_Av;
			H(i+1,i) = norm(v);
			V.col(i+1) = v/H(i+1,i);
			h = H.block(0,i+1,0,i);
			y = h.pinv()*(e1.block(0,i+1,0,0)*r0_norm);
			if(norm(e1.block(0,i+1,0,0)*r0_norm-h*y) < 0.01 || i == K-1){
				x += V.block(0,-1,0,i)*y.block(0,i,0,0);
				printf("m:%d\n",i);
				x.show();
				break;
			}
		}
	}
}/*
int main(){
	mat A(N,N),d(N,1),b(N,1);
	for(int i = 0;i < N*N;i++) A[i] = rand()%10;
	for(int i = 0;i < N;i++) b[i] = rand()%10;
	LU<double> lu(A);
	d = lu.solve(b);
	//d.show();
	printf("A\n");
	A.show();
	printf("b\n");
	b.show();
	printf("x = A⁻¹b\n");
	(A.pinv()*b).show();
	printf("x (LU分解)\n");
	d.show();
	printf("x (GMRES)\n");
	GMRES(A,b,N);*/
	/*
	Model m;
	m.state_func = state_func;
	NMPC nmpc(&m,step,dt,0.2,2,3);
	nmpc.setCostFunc(state_cost,input_cost,terminal_cost);
	mat U(step,2);
	mat X(3,step);
	mat x(3,1);
	x.set({1,1,-M_PI/2.0});
	U.set({2,0});
	U = nmpc.newton(x,U);
	U = nmpc.newton(x,U);
	U = nmpc.newton(x,U);
	vector<double> xx(step),yy(step),dx,dy;
	for(int i = 0;i < 1000;i++){
		//U.show();
		U = nmpc.GMRES(x,U);
		x = x + m.state_func(x,U.row(0).T())*dt;
		printf("F:%lf\n",norm(nmpc.F_(x,U)));
		X = nmpc.X_(x,U);
		for(int k = 0;k < step;k++){
			xx[k] = X.row(k).at(0);
			yy[k] = X.row(k).at(1);
		}
		dx.push_back(x.at(0));
		dy.push_back(x.at(1));
		plt::cla();
		plt::plot({0},{0},"-gD");
		plt::plot(dx,dy,"-b");
		plt::plot(xx,yy,"-r");
		plt::pause(0.01);
	}
	plt::show();
*/



int N = 10;
double step = 10;
double dt = 0.01;
int main(){
	TwoWheel tw;
	NMPC nmpc(&tw,step,dt,0.01);
	mat U(step,2);
	mat X(3,step);
	mat x(3,1);
	x.set({-10,1.0,-M_PI/2.0});
	vector<double> dx,dy;
	//tw.addObj(-8,0.3);
	//tw.addObj(-6,0.1);
	for(int i = 0;i < 2;i++) U = nmpc.GRADIENT(x,U);
	for(int i = 0;i < 1000;i++){
		//U.show();
		//nmpc.F(x,U).show();
		U = nmpc.GRADIENT(x,U);
		X = nmpc.X(x,U);
		x += tw.f(x,U.row(0).T())*dt;
		U.row(0).show();
		dx.push_back(x.at(0));
		dy.push_back(x.at(1));
		plt::cla();
		plt::xlim(-13,13);
		plt::ylim(-13,13);
		plt::plot({0.0},{0.0},"gD");
		plt::plot({X[0]},{X[1]},"bD");
		//plt::plot(dx,dy,"-b");
		plt::plot(X.col(0).getelem(),X.col(1).getelem(),"-r");
		plt::pause(0.01);
		if(fabs(x[0]) < 0.5 && fabs(x[1]) < 0.5) {
			x[0] = rand()%20-10;
			x[1] = rand()%20-10;
		}
	}
	plt::show();
	return 0;
}