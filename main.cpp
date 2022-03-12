#include <unistd.h>
#include "matrix.h"
#include "bmd.hpp"

#include <chrono>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
typedef Matrix<double> matrix;

matrix givensGMRES(matrix A,matrix b,matrix x,uint m){
	for(int k = 0;k < 30;k++){
		matrix r0 = b - A*x;
		double r0_norm = norm(r0);
		matrix V(r0.rows(),m);
		matrix g(m+1,1);
		matrix H(m+1,m);
		V.col(0) = r0/r0_norm;
		g(0) = r0_norm;
		for(int j = 0;j < m;j++){
			matrix AVj = A*V.col(j);
			for(int i = 0;i <= j;i++) H(i,j) = dot(AVj,V.col(i));
			for(int i = 0;i <= j;i++) AVj -= H(i,j)*V.col(i);
			H(j+1,j) = norm(AVj);
			if(j < m-1) V.col(j+1) = AVj/H(j+1,j);
		}
		matrix y = H.pinv()*g;
		x += V.block(0,-1,0,m-1)*y;
	}
	return x;
}

int N = 10;
double step = 5;
double dt = 0.01;
int main(){
	BMD bmd;
	NMPC nmpc(&bmd,step,dt,100);
	matrix U(step,3,0.1);
	matrix x(2,1);

	x.set({2.0,0});
	double t = 0;
	vector<double> x1,x2,u1,u2,T;
	U.col(0) += 1.;
	matrix u(1,3);
	for(int i = 0;i < 1;i++) (bmd.Hu(x.T(),u,bmd.phix(x,t),t)).show();
	for(int i = 0;i < 2000;i++){
		printf("-------%d-------\n",i);
		int start = clock();

		U = nmpc.GMRES(x,U);
		int end = clock();
		u = U.row(0).T();
		u = u.block(0,0,0,1);
		x += bmd.f(x,u)*dt;
		t += dt;
		printf("x1:%lf,x2:%lf\n",x[0],x[1]);
		printf("time:%lf\n",double(end-start)/CLOCKS_PER_SEC);
		printf("-------%d-------\n",i);
		x1.push_back(x(0));
		x2.push_back(x(1));
		u1.push_back(u(0));
		u2.push_back(u(1));
		T.push_back(t);
		//plt::subplot2grid(5,1,0,0,4,1);
		//plt::subplot2grid(5,1,0,0,4,2);
	}
	plt::grid(true);
	plt::plot(T,x1,"-r");
	plt::plot(T,x2,"-g");
	plt::plot(T,u1,"-y");
	plt::plot(T,u2,"-b");
	plt::show();
	return 0;
	
}