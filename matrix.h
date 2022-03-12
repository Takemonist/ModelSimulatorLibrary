#ifndef PTL_MATRIX
#define PTL_MATRIX
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>

using std::vector;
using std::string;

template <class type> class Matrix;
template <class type> class BlockMatrix;
template <class type> class LU;
Matrix<double> E(uint N);
Matrix<double> diag(uint N,double val);
Matrix<double> diag(vector<double> vals);
double norm(Matrix<double> A);
double dot(Matrix<double> A,Matrix<double> B);

template <class type> class LU{
public:
    LU(Matrix<type>& A){
        if(A.rows() == A.cols()){
            N = A.rows();
            L = Matrix<type>(N,N);
            U = E(N);
            for(uint i = 0; i < N; ++i){
                for(uint j = 0; j <= i; ++j){
                    double lu = A(i,j);
                    for(uint k = 0; k < j; ++k){
                        lu -= L(i,k)*U(k,j);
                    }
                    L(i,j) = lu;
                }
                // u_ijの計算(i < j)
                for(uint j = i+1; j < N; ++j){
                    double lu = A(i,j);
                    for(uint k = 0; k < i; ++k){
                        lu -= L(i,k)*U(k,j);
                    }
                    U(i,j) = lu/L(i,i);
                }
            }
        }
    }
    Matrix<type> solve(Matrix<type> b){
        Matrix<type> Lb(N,N+1);
        Lb.block(0,-1,0,-2) = L;
        Lb.col(N) = b;
        for(uint i = 0;i < N;i++){
            Lb.row(i) /= Lb.at(i,i);
            for(uint j = i+1;j < N;j++){
                Lb.row(j) -= Lb.row(i)*Lb.at(j,i);
            }
        }
        Matrix<type> Uy(N,N+1);
        Uy.block(0,-1,0,-2) = U;
        Uy.col(N) = Lb.col(N);
        for(int i = N-1;i >= 0;i--){
            for(int j = i-1;j >= 0;j--){
                Uy.row(j) -= Uy.row(i)*Uy.at(j,i);
            }
        }
        return Uy.col(N);
    }
    uint N;
    Matrix<type> L;
    Matrix<type> U;
};

template <class type> class Matrix
{
private:
    uint R,C;
    vector<type> elements;
public:
void show(string str = "")
    {
        std::cout << str << std::endl;
        printf("[");
        for(uint i = 0; i < R; i++) {
            if(i == 0) printf("[");
            else printf(" [");
            for(uint j = 0; j < C; j++) {
                std::cout << std::fixed << std::setw(10) << std::right << std::setprecision(3) << at(i,j);
                if(j != C-1) printf(" ");
            }
            printf("]");
            if(1 != R && i != R-1) printf("\n");
        }
        printf("]\n\n");
    }
    Matrix():R(0),C(0),elements(0,0){}
    Matrix(uint R,uint C,type init_val = 0):R(R),C(C),elements(R*C,init_val){}
    Matrix(BlockMatrix<type> mat):R(mat.rows()),C(mat.cols()),elements(R*C){
        for(uint i = 0; i < R*C; i++) at(i) = mat.at(i);
    }
    virtual type& at(uint r,uint c) {return elements[C*r+c];}
    virtual type& at(uint n) {return elements[n];}
    type& operator ()(uint r,uint c) {return at(r,c);}
    BlockMatrix<type> operator ()(int rb,int re,int cb,int ce) {return block(rb,re,cb,ce);}
    type& operator ()(uint n) {return at(n);}
    type& operator [](uint n) {return at(n);}
    Matrix<type>& reshape(uint r,uint c) {
        R = r;
        C = c;
        elements.resize(R*C);
        return *this;
    }
    vector<type> getelem(){
        vector<type> ret(R*C);
        for(uint i = 0;i < R*C;i++) ret.at(i) = at(i);
        return ret;
    }
    virtual uint rows() {return R;}
    virtual uint cols() {return C;}
    BlockMatrix<type> row(uint r){
        BlockMatrix<type> ret(1,C);
        for(uint i = 0;i < C;i++) ret.ptr_at(0,i) = &at(r,i);
        return ret;
    }
    BlockMatrix<type> col(uint c){
        BlockMatrix<type> ret(R,1);
        for(uint i = 0;i < R;i++) ret.ptr_at(i,0) = &at(i,c);
        return ret;
    }
    BlockMatrix<type> block(int rb,int re,int cb,int ce) {
        if(re < 0) re += int(R);
        if(ce < 0) ce += int(C);
        BlockMatrix<type> ret(re-rb+1,ce-cb+1);
        for(uint i = rb;i <= re;i++){
            for(uint j = cb;j <= ce;j++){
                ret.ptr_at(i-rb,j-cb) = &at(i,j);
            }
        }
        return ret;
    }
    void set(vector<type> elem) {for(uint i = 0; i < R*cols() && i < elem.size(); i++) at(i) = elem.at(i);}
    Matrix<type>& operator = (Matrix<type> rhs) {
        reshape(rhs.R,rhs.C);
        elements = rhs.elements;
        return *this;
    }
    Matrix<type> operator + (Matrix<type> rhs) {
        Matrix<type> ret(R,C);
        if(R != rhs.R || C != rhs.C) std::cout << "rows or cols does not match." << std::endl;
        for(uint i = 0; i < R*cols(); i++) ret.at(i) = at(i) + rhs.at(i);
        return ret;
    }
    Matrix<type> operator - (Matrix<type> rhs) {
        Matrix<type> ret(R,C);
        if(R != rhs.R || C != rhs.C) std::cout << "rows or cols does not match." << std::endl;
        for(uint i = 0; i < R*cols(); i++) ret.at(i) = at(i) - rhs.at(i);
        return ret;
    }
    Matrix<type> operator * (Matrix<type> rhs) {
        Matrix<type> ret(R,rhs.C);
        if(C != rhs.R) std::cout << "rows or cols does not match." << std::endl;
        for(uint i = 0; i < R; i++) {
            for(uint j = 0; j < rhs.C; j++) {
                for(uint k = 0; k < C; k++) {
                    ret.at(i,j) += at(i,k) * rhs.at(k,j);
                }
            }
        }
        return ret;
    }
    Matrix<type>& operator += (Matrix<type> rhs) {
        if(R != rhs.R || C != rhs.C) std::cout << "rows or cols does not match." << std::endl;
        else for(uint i = 0; i < R*cols(); i++) at(i) += rhs.at(i);
        return *this;
    }
    Matrix<type>& operator -= (Matrix<type> rhs) {
        if(R != rhs.R || C != rhs.C) std::cout << "rows or cols does not match." << std::endl;
        else for(uint i = 0; i < R*cols(); i++) at(i) -= rhs.at(i);
        return *this;
    }
    Matrix<type>& operator = (BlockMatrix<type> rhs) {
        reshape(rhs.rows(),rhs.cols());
        for(uint i = 0; i < R*C; i++) at(i) = rhs.at(i);
        return *this;
    }
    Matrix<type> operator + (BlockMatrix<type> rhs) {
        Matrix<type> ret(R,C);
        if(R != rhs.rows() || C != rhs.cols()) std::cout << "rows or cols does not match." << std::endl;
        for(uint i = 0; i < R*C; i++) ret.at(i) = at(i) + rhs.at(i);
        return ret;
    }
    Matrix<type> operator - (BlockMatrix<type> rhs) {
        Matrix<type> ret(R,C);
        if(R != rhs.rows() || C != rhs.cols()) std::cout << "rows or cols does not match." << std::endl;
        for(uint i = 0; i < R*C; i++) ret.at(i) = at(i) - rhs.at(i);
        return ret;
    }
    Matrix<type> operator * (BlockMatrix<type> rhs) {
        Matrix<type> ret(R,rhs.cols());
        if(C != rhs.rows()) std::cout << "rows or cols does not match." << std::endl;
        for(uint i = 0; i < R; i++) {
            for(uint j = 0; j < rhs.cols(); j++) {
                for(uint k = 0; k < C; k++) {
                    ret.at(i,j) += at(i,k) * rhs.at(k,j);
                }
            }
        }
        return ret;
    }
    Matrix<type>& operator += (BlockMatrix<type> rhs) {
        if(R != rhs.rows() || C != rhs.cols()) std::cout << "rows or cols does not match." << std::endl;
        else for(uint i = 0; i < R*C; i++) at(i) += rhs.at(i);
        return *this;
    }
    Matrix<type>& operator -= (BlockMatrix<type> rhs) {
        if(R != rhs.rows() || C != rhs.cols()) std::cout << "rows or cols does not match." << std::endl;
        else for(uint i = 0; i < R*C; i++) at(i) -= rhs.at(i);
        return *this;
    }

    Matrix<type> operator + (type scalar) {
        Matrix<type> ret(R,C);
        for(uint i = 0; i < R*C; i++) ret.at(i) = at(i) + scalar;
        return ret;
    }
    Matrix<type> operator - (type scalar) {
        Matrix<type> ret(R,C);
        for(uint i = 0; i < R*C; i++) ret.at(i) = at(i) - scalar;
        return ret;
    }
    Matrix<type> operator * (type scalar) {
        Matrix<type> ret(R,C);
        for(uint i = 0; i < R*C; i++) ret.at(i) = at(i) * scalar;
        return ret;
    }
    Matrix<type> operator / (type scalar) {
        Matrix<type> ret(R,C);
        for(uint i = 0; i < R*C; i++) ret.at(i) = at(i) / scalar;
        return ret;
    }
    Matrix<type> operator % (type scalar) {
        Matrix<type> ret(R,C);
        for(uint i = 0; i < R*C; i++) ret.at(i) = std::remainder(at(i),scalar);
        return ret;
    }
    Matrix<type> operator - () {
        return *this*-1;
    }
    Matrix<type>& operator += (type scalar) {
        for(uint i = 0; i < R*C; i++) at(i) += scalar;
        return *this;
    }
    Matrix<type>& operator -= (type scalar) {
        for(uint i = 0; i < R*C; i++) at(i) -= scalar;
        return *this;
    }
    Matrix<type>& operator *= (type scalar) {
        for(uint i = 0; i < R*C; i++) at(i) *= scalar;
        return *this;
    }
    Matrix<type>& operator /= (type scalar) {
        for(uint i = 0; i < R*C; i++) at(i) /= scalar;
        return *this;
    }
    Matrix<type> T()
    {
        Matrix<type> ret(C,R);
        for(uint i = 0; i < C; i++) {
            for(uint j = 0; j < R; j++) {
                ret.at(i,j) = at(j,i);
            }
        }
        return ret;
    }
    double det(){
        if(R != C) return 0;
        else{
            uint N = R;
            double ret = 1.0;
            Matrix<type> A = *this;
            for(uint i = 0;i < N; i++){
                for(uint j = 0;j < i; j++){
                    if(A(i,j) != 0.0){
                        if(A(j,j) == 0.0) A.row(j) += A.row(i);
                        A.row(i) -= A.row(j)*A(i,j)/A(j,j);
                    }
                }
            }
            for(uint i = 0;i < N;i++) ret *= A(i,i);
            return ret;
        }
    }
    Matrix<type> inv(){
        if(R != C) return *this;
        else{ 
            LU<type> lu(*this);
            Matrix<type> ret = E(R);
            for(uint i = 0;i < R;i++){
                ret.col(i) = lu.solve(ret.col(i)); 
            }
            return ret;
        }
    }
    Matrix<type> pinv(){
        Matrix<type> A = *this;
        return (A.T()*A).inv()*A.T();
    }
    Matrix<type>& cleanup(){
        for(uint i = 0;i < R*C;i++) if(fabs(at(i)) < 10e-8) at(i) = 0;
        return *this;
    }
};

template <class type> Matrix<type> operator * (type scalar,Matrix<type> rhs) {
        return rhs*scalar;
}
template <class type> Matrix<type> operator * (type scalar,BlockMatrix<type> rhs) {
        return rhs*scalar;
}
template <class type> class BlockMatrix:public Matrix<type>
{
private:
    uint R,C;
    vector<type*> elements;
public:
    BlockMatrix(uint R,uint C):R(R),C(C),elements(R*C,nullptr){this->reshape(R,C);}
    type& at(uint r,uint c) {return *elements[C*r+c];}
    type& at(uint n) {return *elements[n];}
    type& operator ()(uint r,uint c) {return at(r,c);}
    type& operator ()(uint n) {return at(n);}
    type& operator [](uint n) {return at(n);}
    type*& ptr_at(uint r,uint c) {return elements[C*r+c];}
    uint rows() {return R;}
    uint cols() {return C;}
    void set(vector<type*> elem) {for(uint i = 0; i < R*C && i < elem.size(); i++) elements[i] = elem[i];}
    BlockMatrix<type>& operator = (Matrix<type> rhs){
        if(R != rhs.rows() || C != rhs.cols()) std::cout << "rows or cols does not match." << std::endl;
        else for(uint i = 0; i < R*C; i++) at(i) = rhs.at(i);
        return *this;
    }
    BlockMatrix<type>& operator = (BlockMatrix<type> rhs){
        if(R != rhs.rows() || C != rhs.cols()) std::cout << "rows or cols does not match." << std::endl;
        else for(uint i = 0; i < R*C; i++) at(i) = rhs.at(i);
        return *this;
    }
    operator Matrix<type>(){
        Matrix<type> ret(R,C);
        for(uint i = 0; i < R*C; i++) ret.at(i) = at(i);
        return ret;
    }
};

Matrix<double> E(uint N){
    Matrix<double> ret(N,N);
    for(uint i = 0;i < N;i++ ) ret(i,i) = 1.0;
    return ret;
}
Matrix<double> diag(uint N,double val){
    Matrix<double> ret(N,N);
    for(uint i = 0;i < N;i++ ) ret(i,i) = val;
    return ret;
}
Matrix<double> diag(vector<double> vals){
    uint N = vals.size();
    Matrix<double> ret(N,N);
    for(uint i = 0;i < N;i++ ) ret(i,i) = vals.at(i);
    return ret;
}
double norm(Matrix<double> A){
    return sqrt((A.T()*A).at(0));
}
double dot(Matrix<double> A,Matrix<double> B){
    return (A.T()*B).at(0);
}
#endif