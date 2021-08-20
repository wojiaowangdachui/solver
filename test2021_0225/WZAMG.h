#pragma once
#ifndef WZAMG_H
#define WZAMG_H
#include<math.h>
#include<vector>
#include<algorithm>
#include<fstream>
#include <numeric>
#include<Eigen\Sparse>//包含稀疏矩阵求解;
#include<Eigen\Dense>
#include<unordered_map>
#include"matrix.h"
#include"ST.h"

//从包含角度上讲是只要包含此文件就行，不用再包含"matrix.h";
typedef Eigen::SparseMatrix<double> SparseMatrixType;
typedef Eigen::Triplet<double> T;
typedef Eigen::SimplicialCholesky<SparseMatrixType> Solve;



class AMG {
public:
	AMG() = delete;
	AMG(const TernArray& temp,double* tempb, int tempcycle, int tempcoarsest, double tempSW_BOUND);
	~AMG();
	int UseEigen(int level,OneDimensionalArray* tempx, OneDimensionalArray* tempb);
	double Smooth(int level,int max_iter, OneDimensionalArray* b, OneDimensionalArray* u_in, OneDimensionalArray* u_apx, double tempErr=1.0e-10);
//	int	SolveEquationsJPCG(TernArray* tempA, OneDimensionalArray* temp_b, int maxcycle, OneDimensionalArray& pX);
//	int	JPCG(TernArray* tempA, OneDimensionalArray* temp_b, int maxcycle, OneDimensionalArray& pX);
	void my_strong(int level, double local_SW_BOUND);
	void amg_setup_wz(int level);
	void my_AMG(int level);
	double amg_cycle(int level, OneDimensionalArray* tempb, OneDimensionalArray* u_in, OneDimensionalArray* u_out);
	int MainReSult(OneDimensionalArray* x, int v_max);
	int SolveEquationsCG(TernArray* tempA, int m_nFiledMaxItTime, int nMatrixNum, double* pX, double* pB);
	int IsFieldConverge(double* temp, int nmax);
	double NormErr(TernArray* tempA, OneDimensionalArray* x, OneDimensionalArray* b, OneDimensionalArray* tempu_out);
	double proCGmethode(double* R, double* X);
	void amg_get_Gnode(int& nc, int* Gnode,TernArray* A);
	void amg_get_W(int level);

private:

	TernArray* Amatrix = NULL;
	double* b = NULL;
	int rowmax;
	int colmax;
	CompleteMatrix* Cmatrix = NULL;
	SparseMatrixType* A = NULL;
	int COARSEST;
	int CYCLES;
	double SW_BOUND;
	TernArray** A_matrix = NULL;
	TernArray** Wweight=NULL;
	TernArray** TrangeWweight=NULL;

	/////////
//	int* WCOL = NULL;
	////////
};

#endif