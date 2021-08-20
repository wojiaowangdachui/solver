#include"OneArray.h"
#include<iostream>

void OneDimensionalArray::Printf()
{
	for (int i = 0; i < this->max; i++) {
		std::cout << this->b[i] << std::endl;
	}
}

void OneDimensionalArray::copy(const OneDimensionalArray& temp)
{
	this->Initial(temp.max, temp.b);
}

void OneDimensionalArray::reduce(double* b1)
{
	for (int i = 0; i < max; i++) {
		this->b[i] = this->b[i] - b1[i];
	}
}

void OneDimensionalArray::ReverseB(double** temp)
{
	if (*temp != NULL) {
		delete[] *temp;
	}
	*temp = new double[this->max];
	std::memcpy(*temp, this->b, sizeof(double) * max);
}

void OneDimensionalArray::InitialB(double* b)
{ 
	memcpy(b, this->b, sizeof(double) * max); 
}

double OneDimensionalArray::Fabs()
{
	double temp = 0;
	for (int i = 0; i < this->max; i++) {
		temp = temp + (this->b[i]);
		//if (temp == NAN) {
		//	system("pause");
		//}
	 //  if (temp ==nan) {
		//	system("pause");
		//}
	}
	return temp;
}

void OneDimensionalArray::Free()
{
	delete this->b;
	this->b = new double[this->max];
	memset(this->b, 0, sizeof(double) * this->max);
}


OneDimensionalArray& OneDimensionalArray::operator=(const OneDimensionalArray& Temp)
{
	if (this->b != NULL) {
		OneDimensionalArray temp(Temp);
		this->max = temp.max;
		double* tempb = this->b;
		this->b = temp.b;
		temp.b = tempb;
		return *this;
	}
	this->max = Temp.max;
	this->b = new double[this->max];
	std::memcpy(this->b, Temp.b, sizeof(double) * max);
	return *this;
}

void OneDimensionalArray::MultiplyNum(double temp,OneDimensionalArray* result)
{
	if (result != NULL) {
		for (int i = 0; i < this->max; i++) {
			result->b[i] = this->b[i] * temp;
		}
		return;
	}
	for (int i = 0; i < this->max; i++) {
		this->b[i] = this->b[i] * temp;
	}
}

double OneDimensionalArray::Multiply(const OneDimensionalArray& multiply)
{
	double temp = 0;
	for (int i = 0; i < this->max; i++) {
		temp = temp + this->b[i] * multiply.b[i];
	}
	return temp;
}

double OneDimensionalArray::MultiplyOneSelf()
{
	double temp = 0;
	for (int i = 0; i < this->max; i++) {
		temp = temp + b[i] * b[i];
	}
	return temp;
}

OneDimensionalArray::OneDimensionalArray(const OneDimensionalArray& Temp)
{
	this->max = Temp.max;
	this->b = new double[max];
	std::memcpy(this->b, Temp.b, sizeof(double) * max);
}


//void OneDimensionalArray::Change(int length,double* Temp)
//{
//	this->max = length;
//	delete[] this->b;
//	this->b = new double[this->max];
//	memcpy(this->b, Temp, sizeof(double) * max);
//}

//void OneDimensionalArray::Change(OneDimensionalArray& Temp)
//{
//	if (this->max != 0) {
//		this->max = Temp.max;
//	}
//	delete this->b;
//	this->b = new double[max];
//	memcpy(this->b, Temp.b, sizeof(double) * max);
//}

double OneDimensionalArray::norm2()
{
	double temp = 0;
	for (int i = 0; i < this->max; i++) {
		if (fabs(this->b[i]) > 1e-10) {
			temp = this->b[i] * this->b[i] + temp;
		}
	}
	double  temp1 = pow(temp, 0.5);
	return temp1;
}

void OneDimensionalArray::add(OneDimensionalArray* add1, OneDimensionalArray* result)
{
	if (result != NULL) {
		for (int i = 0; i < this->max; i++) {
			result->b[i] = this->b[i] + add1->b[i];
		}
	}
	else {
		for (int i = 0; i < this->max; i++) {
			this->b[i] = this->b[i] + add1->b[i];
		}
	}
}

void OneDimensionalArray::reduce(OneDimensionalArray* reduce1, OneDimensionalArray* result)
{
	for (int i = 0; i < this->max; i++) {
		result->b[i] =  this->b[i]- reduce1->b[i] ;
	}
}

void OneDimensionalArray::reduce(OneDimensionalArray* reduce1,bool chose)
{
	if (chose) {
		for (int i = 0; i < this->max; i++) {
			this->b[i] = this->b[i] - reduce1->b[i];
		}
	}
	else {
		for (int i = 0; i < this->max; i++) {
			reduce1->b[i] = this->b[i] - reduce1->b[i];
		}
	}
}

double OneDimensionalArray::ReturnNum(int i)
{
	return this->b[i];
}

OneDimensionalArray::OneDimensionalArray(int nmax)
{
	max = nmax;
	b = new double[max];
	std::memset(b, 0.0, sizeof(double) * max);
}

OneDimensionalArray::~OneDimensionalArray()
{
	if (this->b != NULL) {
		delete[] this->b;
		this->b = NULL;
	}
}

//void OneDimensionalArray::Trange(OneDimensionalArray* TempB)
//{
//	trange = 1;
//}

void OneDimensionalArray::Initial(int length, double* b)
{
	this->max = length;
	std::memcpy(this->b, b, sizeof(double) * max);
}


void OneDimensionalArray::quickSort(int left, int right)
{
	if (left >= right)
		return;
	int i, j, base, temp;
	i = left, j = right;
	base = this->b[left];  //取最左边的数为基准数
	while (i < j)
	{
		while (this->b[j] >= base && i < j)
			j--;	
		while (this->b[i] <= base && i < j)
			i++;
		if (i < j)
		{
			temp = this->b[i];
			this->b[i] = this->b[j];
			this->b[j] = temp;
		}
	}
	//基准数归位
	this->b[left] = this->b[i];
	this->b[i] = base;
	quickSort(left, i - 1);//递归左边
	quickSort(i + 1, right);//递归右边
}

void  OneDimensionalArray::Assignment(int i, double num)
{
	this->b[i] = num;
}

void OneDimensionalArray::FindMax(int& imax, double& vmax)
{
	imax = 0;
	vmax = this->b[0];
	for (int i = 0; i < this->max; i++) {
		if (this->b[i] > vmax) {
			imax = i;
			vmax = this->b[i];
		}
	}
}

void OneDimensionalArray::FindMin(int& imax, double& vmax)
{
	imax = 0;
	vmax = this->b[0];
	for (int i = 0; i < this->max; i++) {
		if (this->b[i] < vmax) {
			imax = i;
			vmax = this->b[i];
		}
	}
}