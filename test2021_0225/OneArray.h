#pragma once
#ifndef OneArray
#define OneArray
#include<stdio.h>
class OneDimensionalArray
{
public:
	OneDimensionalArray() = delete;
	explicit OneDimensionalArray(int nmax);
	OneDimensionalArray(const OneDimensionalArray& Temp);
	~OneDimensionalArray();
	//void Trange(OneDimensionalArray* TempB);
	//void Change(int length,double* Temp);
	OneDimensionalArray& operator=(const OneDimensionalArray& Temp);
	void Initial(int length, double* b);
	void copy(const OneDimensionalArray& temp);
	void Assignment(int i, double num);
	void FindMax(int& imax, double& vmax);
	void FindMin(int& imax, double& vmax);
	double ReturnNum(int i);
	void add(OneDimensionalArray* add1, OneDimensionalArray* result=NULL);
	void reduce(OneDimensionalArray* reduce1, OneDimensionalArray* result);
	void reduce(OneDimensionalArray* reduce1, bool chose);
	void reduce(double* b1);
	double Multiply(const OneDimensionalArray& multiply);
	double MultiplyOneSelf();
	void MultiplyNum(double temp,OneDimensionalArray* result=NULL);
	void Printf();
	double norm2();
	double ReturnMax() { return this->max; }
	void ReverseB(double** b);
	void InitialB(double* b);
	double Fabs();
	void Free(); 
	void quickSort(int left, int right);
private:
	int max;
	double* b;
};


#endif // !OneArray
