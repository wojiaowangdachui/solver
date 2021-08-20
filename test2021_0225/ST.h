#pragma once
#ifndef ST_h
#define ST_h
#include<vector>
#include<algorithm>
#include <numeric>
#include<set>
#include"matrix.h"

class STMatrix 
{
public:
	STMatrix()=delete;
	STMatrix(const STMatrix& temp)=delete;
	STMatrix& operator=(const STMatrix& temp) = delete;
	STMatrix(int temprowmax,int tempcolmax,int tempallnum);
	STMatrix(int** tempST,int *len,int nmax, int tempall);
	void OrderInitial(int* i, int length,int row);
	~STMatrix();
	void GetAMGNode(TernArray** Weight, int level);
	int RowMax(int* STRow,int length);

private:
	int* STRowArray = NULL;
	int** ST=NULL ;
	int rowmax;
	int colmax;
	int allnum=0;

};


#endif // 

