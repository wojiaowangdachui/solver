#pragma once
#ifndef MATRIX_H
#define MATRIX_H

#include<vector>
#include<iostream>
#include<math.h>
#include<Eigen\Sparse>//包含稀疏矩阵求解;
#include<Eigen\Dense>
#include<numeric>
#include<iterator>
#include<utility>
#include<omp.h>
#include"OneArray.h"

typedef Eigen::SparseMatrix<double> SparseMatrixType;
typedef Eigen::Triplet<double> T;
typedef Eigen::SimplicialCholesky<SparseMatrixType> Solve;

class AMG;

struct node
{
	int col;
	double num;
};

template<class T>
struct csrnode
{
	int row;
	int col;
	T num;
};

class WZmatrix
{
public:
	WZmatrix() = delete;
	explicit WZmatrix(int nmax) {
		Node = new node[nmax];
		counter = 0;
		max = nmax;
	}
	~WZmatrix() {
		if (Node != NULL) {
			delete Node;
			Node = NULL;
		}
	}

	void OrderAdd(int i, double j);
//	void Sort(int ColPos);
	int Return(int* Col, double* Num, int count);
	int ReturnAllnum();

//	std::vector<int>Pos;
	node* Node;
private:
	int counter;
	int max;
};


class CompleteMatrix
{
public:
	CompleteMatrix() = delete;
	CompleteMatrix(int nmax,int col):max(nmax),colmax(col)
	{
		matrix = new double* [max];
		for (int i = 0; i < max; i++) {
			matrix[i] = new double[colmax];
			std::memset(matrix[i], 0, sizeof(double) * colmax);
		}
	}

	~CompleteMatrix()
	{
		if (matrix != NULL) {
			for (int i = 0; i < max; i++) {
				delete[] matrix[i];
				matrix[i] = NULL;
			}
			delete [] matrix;
			matrix=NULL;
		}
	}

	void Assignment(int i, int j, double num)
	{
		matrix[i][j] = num;
	}

	void Printf()
	{
		for (int i = 0; i < max; i++) {
			for (int j = 0; j < colmax; j++) {
				std::cout << matrix[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}

private:
	//CompleteMatrix() {}
	double** matrix;
	int max;
	int colmax;
};



class TernArray
{
public:
	friend class AMG;
	friend class matrix;
	//vs2010不支持C++11
	TernArray() = delete;
//	explicit TernArray(int max);
	TernArray(int* temprow,int* tempcol,double* tempnum,int temprowmax,int tempcolmax,int nummax);
//	TernArray(int max, int nummax);
	TernArray(int rowmax, int colmax, int nummax);
	TernArray(const TernArray& other);
//	TernArray(const TernArray& other1, const TernArray& other2);
	~TernArray();
	TernArray& operator=(const TernArray& temp);
	void add(TernArray* result, TernArray* add2);
	void reduce(TernArray* result, const TernArray* minute);
	//void putTranspose(const TernArray& TempArray);
	void transpose( TernArray* TempArray);
	double GetNum(int i, int j);
	void PrintfTernArray();
	void PrintfMatrix();
	void triu();
	int GetAllnum();
	int GetNmax() { return this->nmax; }
	int GetColMax() { return this->colmax; }
	void GetRowArray(const int k,node* TempK);
	void FindRowArray(int tempi,std::vector<int>& temp);
	int GetRowNoZero(const int k);
	void OrderAssignment(int i, int j, double num);
	void Multiply(TernArray* MultiplyOther, TernArray* ResultOther);
	double Addallnum();
	void MultiplyVector(OneDimensionalArray* mutiply, OneDimensionalArray* result);
	void MutiplyAllNum(double tempb);
	void PowerAllNum(double tempb,TernArray * tempArray);
	void AddMemcpy(int* start1, int* last1, double* start2, double* last2, int length);
	void ChangeRow(int temprow);
	double MutiplyIJ(int tempi, int j,std::vector<std::pair<int, double>> Tij);
	int SolveSPDMLU(double* pX, double* pB);

private:
	TernArray* Dialog = NULL;
	TernArray* LowerTrangle = NULL;
	TernArray* UpperTrangle = NULL;
	//TernArray* Transpose = NULL;
	int* row=NULL;
	int* col=NULL;
	double* num=NULL;
	int allnum = 0;
	int nmax;
	int colmax;
	CompleteMatrix* Cmatrix=NULL;
	SparseMatrixType* A=NULL;

};




#endif