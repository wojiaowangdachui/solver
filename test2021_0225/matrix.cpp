#include<omp.h>
#include<algorithm>
#include<memory>
#include<map>
#include"matrix.h"

/*******************ֻ���˳�����*****************/
void WZmatrix::OrderAdd(int i, double j)
{
	int num = omp_get_thread_num();
	//no sure
	if (counter == max) {
		{
			WZmatrix TempNode(2 * max);
			std::memcpy(TempNode.Node, this->Node, sizeof(node) * counter);
			delete this->Node;
			int temp = 2 * max;
			this->Node = new node[temp];
			std::memcpy(this->Node, TempNode.Node, sizeof(node) * counter);
			max = 2 * max;
		}
	}
	Node[counter].col = i;
	Node[counter].num = j;
	counter++;
}

/// <summary>
/// //////////////////////////////////////////////////////////////////////////////////////////////////����
/// </summary>
/// <returns></returns>

/*******************���ؾ���һ�е�ֵ*****************/
int WZmatrix::ReturnAllnum()
{
	return counter;
}

/*******************����ת�þ���*****************/
int WZmatrix::Return(int* Col, double* Num, int count)
{
	//no sure
	if (count != counter) return 0;
	for (int i = 0; i < counter; i++) {
		Col[i] = Node[i].col;
		Num[i] = Node[i].num;
	}
	return 1;
}

/*******************ѹ������Ĺ���*****************/
TernArray::TernArray(int* temprow, int* tempcol, double* tempnum, int temprowmax, int tempcolmax, int nummax) :nmax(temprowmax), colmax(tempcolmax), allnum(nummax)
{
	int temp1 = nmax + 1;
	row = new int[temp1];
	col = new int[nummax];
	num = new double[nummax];
	std::memcpy(row, temprow, sizeof(int) * temp1);
	std::memcpy(col, tempcol, sizeof(int) * allnum);
	std::memcpy(num, tempnum, sizeof(double) * allnum);

}



void TernArray::AddMemcpy(int* start1, int* last1, double* start2, double* last2, int length)
{
	if (length) {
		std::memcpy(start1, last1, sizeof(int) * length);
		std::memcpy(start2, last2, sizeof(double) * length);
	}
}

/*******************���������*****************/
void TernArray::MultiplyVector(OneDimensionalArray* mutiply1, OneDimensionalArray* result)
{
	if (this->colmax != mutiply1->ReturnMax()|| this->nmax != result->ReturnMax()) {
		std::cout << "������������������������ȣ��޷���" << std::endl;
		return;
	}
	for (int i = 0; i < nmax; i++) {
		double temp = 0;

		for (int j = this->row[i]; j < row[i + 1]; j++) {
			temp = this->num[j] * mutiply1->ReturnNum(this->col[j]) + temp;
		}
		result->Assignment(i, temp);
	}
}

/*******************��������ֵ��tempbָ����***********************/
void TernArray::PowerAllNum(double tempb, TernArray* tempArray = NULL)
{
	for (int i = 0; i < this->row[this->nmax]; i++) {
		this->num[i] = pow(this->num[i], tempb);
	}
}

/*******************���ؾ������Ԫ����*****************/
int TernArray::GetAllnum()
{
	int all = this->row[this->nmax];
	return all;
}
/*******************��������ֵtempb��*****************/
void TernArray::MutiplyAllNum(double tempb)
{
	///////////////////////////////////////////////////////////////////////////�ؼ��ö����allnum�Ƿ���ʵ
	for (int i = 0; i < this->row[this->nmax]; i++) {
		this->num[i] = this->num[i] * tempb;
	}
}

/*******************����*****************/
TernArray::TernArray(int max, int colnmax, int nummax)
{
	int temp1 = max + 1;
	row = new int[temp1];
	col = new int[nummax];
	num = new double[nummax];
	nmax = max;
	colmax = colnmax;
	allnum = nummax;
	std::memset(row, 0, sizeof(int) * temp1);
	std::memset(col, 0, sizeof(int) * allnum);
	std::memset(num, 0, sizeof(double) * allnum);
}

/******************�ı���*****************/
void TernArray::ChangeRow(int temprow)
{
	int tempmax = this->nmax;
	int tempa= temprow+1;
	this->nmax = temprow;
	int* rowarray = new int[tempmax];
	std::memcpy(rowarray, this->row, sizeof(int)*tempa);
	delete this->row;
	this->row = rowarray;
}

/******************ȡ����i�з���Ԫ��ֵ*****************/
void TernArray::FindRowArray(int tempi, std::vector<int>& temp)
{
	for (int i = row[tempi]; i < row[tempi + 1]; i++) {
		temp.push_back(this->col[i]);
	}
}

/*******************���0����ľ������˳��ֵ*****************/
void TernArray::OrderAssignment(int i, int j, double num)
{
	if (this->row[this->nmax] == this->allnum) {

		int tem = 2 * allnum;
		int* newcol = new int[tem];
		double* newnum = new double[tem];
		memcpy(newcol, this->col, sizeof(int) * this->allnum);
		memcpy(newnum, this->num, sizeof(double) * this->allnum);
		delete[] this->col;
		this->col = newcol;
		delete[] this->num;
		this->num = newnum;
		allnum = 2 * allnum;
	}

	for (int k = i + 1; k < (nmax + 1); k++) {
		this->row[k] = this->row[k] + 1;
	}
	this->col[this->row[i + 1]-1] = j;
	this->num[this->row[i + 1]-1] = num;

}

/*******************����i�еķ���ֵ����*****************/
int TernArray::GetRowNoZero(const int k)
{
	int temp = this->row[k + 1] - this->row[k];
	return temp;
}


/*******************����k�еģ�node��ʽ�ķ���Ԫ����ʽ*****************/
void TernArray::GetRowArray(const int k, node* TempK)
{
	int temp = 0;
	for (int i = this->row[k]; i < this->row[k + 1]; i++) {
		TempK[temp].col = this->col[i];
		TempK[temp].num = this->num[i];
	}

}

//void TernArray::triu(TernArray* UpperTriangle, TernArray* LowerTriangle, TernArray* Diagonal)

/*******************ѹ�����������Ǿ���*****************/
void TernArray::triu()
{
	if (this->nmax != colmax) {
		std::cout << "�޷�������" << std::endl;
		return;
	}
	if (this->UpperTrangle == NULL) {
		this->UpperTrangle = new TernArray(this->nmax, this->colmax, this->allnum);
		this->LowerTrangle = new TernArray(this->nmax, this->colmax, this->allnum);
		this->Dialog = new TernArray(this->nmax, this->colmax, this->nmax + 1);

		for (int i = 0; i < nmax; i++) {
			//		int TempRow = this->row[i + 1] - this->row[i];
			this->LowerTrangle->row[i + 1] = this->LowerTrangle->row[i];
			this->UpperTrangle->row[i + 1] = this->UpperTrangle->row[i];
			this->Dialog->row[i + 1] = this->Dialog->row[i];
			for (int j = this->row[i]; j < this->row[i + 1]; j++) {
				if (col[j] < i) {
					//	LowerTriangle->row[i + 1]++;
					this->LowerTrangle->col[this->LowerTrangle->row[i + 1]] = this->col[j];
					this->LowerTrangle->num[this->LowerTrangle->row[i + 1]] = this->num[j];
					this->LowerTrangle->row[i + 1]++;
				}
				else if (col[j] > i) {
					//	UpperTriangle->row[i + 1]++;
					this->UpperTrangle->col[this->UpperTrangle->row[i + 1]] = this->col[j];
					this->UpperTrangle->num[this->UpperTrangle->row[i + 1]] = this->num[j];
					this->UpperTrangle->row[i + 1]++;
				}
				else {
					//	Diagonal->row[i + 1]++;
					this->Dialog->col[this->Dialog->row[i + 1]] = this->col[j];
					this->Dialog->num[this->Dialog->row[i + 1]] = this->num[j];
					this->Dialog->row[i + 1]++;
				}
			}
		}
	}
}

/******************���ѹ������*****************/
void TernArray::PrintfMatrix()
{
	Cmatrix = new CompleteMatrix(this->nmax, this->colmax);

	for (int i = 0; i < nmax; i++) {
		for (int j = this->row[i]; j < row[i + 1]; j++) {
			Cmatrix->Assignment(i, col[j], num[j]);
		}
	}

	std::cout << std::endl;
	Cmatrix->Printf();
	std::cout << std::endl;
	//delete Cmatrix;
	//Cmatrix = NULL;
}


/*******************���ѹ������*****************/
void TernArray::PrintfTernArray()
{
	std::cout << "row =";
	for (int i = 0; i <= nmax; i++) {
		std::cout << row[i] << "  ";
	}
	std::cout << std::endl << "col =";
	for (int i = 0; i < nmax; i++) {
		for (int j = row[i]; j<row[i+1]; j++) {
			std::cout << col[j] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << "num =";
	for (int i = 0; i < allnum; i++) {
		std::cout << num[i] << "  ";
	}
	std::cout << std::endl;
}


TernArray::~TernArray()
{
	delete row;
	row = NULL;
	delete col;
	col = NULL;
	delete num;
	num = NULL;

	if (A != NULL) {
		delete A;
		A = NULL;
	}
	if (this->LowerTrangle != NULL) {
		delete this->LowerTrangle;
		this->LowerTrangle = NULL;
	}

	if (this->UpperTrangle != NULL) {
		delete this->UpperTrangle;
		this->UpperTrangle = NULL;
	}

	if (this->Dialog != NULL) {
		delete this->Dialog;
		this->Dialog = NULL;
	}

	if (this->Cmatrix != NULL) {
		delete Cmatrix;
		Cmatrix = NULL;
	}

	//if (this->Transpose != NULL) {
	//	delete Transpose;
	//	Transpose = NULL;
	//}

}

/*******************��ֵ�����*****************/
TernArray& TernArray::operator=(const TernArray& temp) {
	if (this != &temp) {
		TernArray just(temp);
		this->allnum = just.allnum;
		this->nmax = just.allnum;
		this->colmax = just.colmax;

		if (just.Dialog) {
			TernArray* tempDialog = just.Dialog;
			TernArray* tempLowerTrangle = just.LowerTrangle;
			TernArray* tempUpperTrangle = just.UpperTrangle;

			just.Dialog = this->Dialog;
			just.LowerTrangle = this->LowerTrangle;
			just.UpperTrangle = this->UpperTrangle;

			this->Dialog=tempDialog;
			this->LowerTrangle=tempLowerTrangle;
			this->UpperTrangle =tempUpperTrangle;
		}
		//if (just.Transpose) {
		//	TernArray* tempTrange = just.Transpose;

		//	just.Transpose = this->Transpose;

		//	this->Transpose = tempTrange;
		//}

		int* temprow = just.row;
		int* tempcol = just.col;
		double* tempnum = just.num;
	    
		just.row=this->row;
	    just.col=this->col;
        just.num=this->num;

		this->row=temprow;
		this->col=tempcol;
		this->num=tempnum;
	}
	return *this;
}

/*******************���ƹ���*****************/
TernArray::TernArray(const TernArray& other)
{
	nmax = other.nmax;
	colmax = other.colmax;
	allnum = other.allnum;
	int temp1 = nmax + 1;
	row = new int[temp1];
	col = new int[allnum];
	num = new double[allnum];
	std::memcpy(row, other.row, sizeof(int) * temp1);
	std::memcpy(col, other.col, sizeof(int) * allnum);
	std::memcpy(num, other.num, sizeof(double) * allnum);
	if (other.Dialog != NULL) {
		this->Dialog = new TernArray(*other.Dialog);
		this->LowerTrangle = new TernArray(*other.LowerTrangle);
		this->UpperTrangle = new TernArray(*other.UpperTrangle);
	}
	//if (other.Transpose !=NULL) {
	//	this->Transpose = new TernArray(*other.Transpose);
	//}
	

}

//TernArray::TernArray(int max, int nummax)
//{
//	int temp1 = max + 1;
//	row = new int[temp1];
//	col = new int[nummax];
//	num = new double[nummax];
//	b = new double[max];
//	nmax = max;
//	colmax = max;
//	allnum = nummax;
//	std::memset(row, 0, sizeof(int) * temp1);
//	std::memset(col, 0, sizeof(int) * allnum);
//	std::memset(num, 0, sizeof(double) * allnum);
//	std::memset(b, 0, sizeof(double) * nmax);
//}

/*******************�������*****************/
void TernArray::reduce(TernArray* result, const TernArray* minute)
{
	if (this->nmax != minute->nmax || this->colmax != minute->colmax) {
		std::cout << "�޷����" << std::endl;
	}
	for (int i = 0; i < this->nmax; i++)
	{
		int max = minute->row[i + 1];
		int start = minute->row[i];
		for (int j = this->row[i]; j < this->row[i + 1]; j++)
		{
			if (start == max){
				result->OrderAssignment(i, this->col[j], this->num[j]);
			}
			if (this->col[j] < minute->col[start]) {
				result->OrderAssignment(i, this->col[j], this->num[j]);
			}
			else if (this->col[j] == minute->col[start]) {
				result->OrderAssignment(i, this->col[j], this->num[j] - minute->num[start]);
				start++;
			}
			else {
				result->OrderAssignment(i, minute->col[start], -1*minute->num[start]);
				start++;
			}
		}
		while (start != max) {
			result->OrderAssignment(i, minute->col[start], -1*minute->num[start]);
			start++;
		}
	}
}

double TernArray::Addallnum()
{
	double ret = 0;
	for (int i = 0; i < this->nmax; i++) {
		for (int j = this->row[i]; j < this->row[i + 1]; j++) {
			ret += fabs(this->num[j]);
		}
	}
	return ret;
}
/*******************�������*****************/
void TernArray::add(TernArray* result, TernArray* add2)
{
	if ((this->nmax != add2->nmax) || (this->colmax != add2->colmax)) {
		std::cout << "�����޷����";
		return;
	}
	//if (*result == NULL) {
	//	*result = new TernArray(this->nmax, this->allnum + add2.allnum);
	//}

	int TempNum = 0;

	for (int i = 0; i < nmax; i++)
	{
		//int RowNoZero1 = this->row[i + 1] - this->row[i];
		//int RowNoZero2 = add2->row[i + 1] - add2->row[i];
		int TempNum1 = this->row[i];
		int TempNum2 = add2->row[i];

		//if (!RowNoZero1) {
		//	//std::memcpy((*result)->col + TempNum, add2.col + TempNum2, sizeof(int) * (RowNoZero2));
		//	//std::memcpy((*result)->num + TempNum, add2.num + TempNum2, sizeof(double) * (RowNoZero2));
		//	this->AddMemcpy((*result)->col + TempNum, add2.col + TempNum2, (*result)->num + TempNum, add2.num + TempNum2, RowNoZero2);
		//	TempNum = TempNum + RowNoZero2;
		//	(*result)->row[i + 1] = TempNum;
		//	continue;
		//}
		//else if (!RowNoZero2) {
		//	//std::memcpy((*result)->col + TempNum, this->col + TempNum1, sizeof(int) * RowNoZero1);
		//	//std::memcpy((*result)->num + TempNum, this->num + TempNum1, sizeof(double) * RowNoZero1);
		//	this->AddMemcpy((*result)->col + TempNum, this->col + TempNum1, (*result)->num + TempNum, this->num + TempNum1, RowNoZero1);
		//	TempNum = TempNum + RowNoZero1;
		//	(*result)->row[i + 1] = TempNum;
		//	continue;
		//}

		while (TempNum1 < this->row[i + 1] && TempNum2 < add2->row[i + 1])
		{
			if (this->col[TempNum1] < add2->col[TempNum2])
			{
				result->col[TempNum] = this->col[TempNum1];
				result->num[TempNum] = this->num[TempNum1];
				TempNum1++;
			}
			else if (this->col[TempNum1] > add2->col[TempNum2])
			{
				result->col[TempNum] = add2->col[TempNum2 ];
				result->num[TempNum] = add2->num[TempNum2 ];
				TempNum2++;
			}
			else
			{
				result->col[TempNum] = this->col[TempNum1 ];
				result->num[TempNum] = this->num[TempNum1 ] + add2->num[TempNum2 ];
				if (!result->num[TempNum]) {
					TempNum = TempNum - 1;
				}
				TempNum1++;
				TempNum2++;
			}
			TempNum++;
		}

		if ((TempNum1 == this->row[i + 1])&&(TempNum2 < add2->row[i + 1])) {
			for (int j = TempNum2; j < add2->row[i + 1]; j++) {
				result->col[TempNum] = add2->col[j];
				result->num[TempNum] = add2->num[j];
				TempNum++;
			}
		}
		else if ((TempNum1 < this->row[i + 1]) && (TempNum2 == add2->row[i + 1])) {
			for (int j = TempNum1;j < this->row[i + 1]; j++) {
				result->col[TempNum] = this->col[j];
				result->num[TempNum] = this->num[j];
				TempNum++;
			}
		}
		//TempNum2 = TempNum2 + k;
		//RowNoZero2 = RowNoZero2 - k;
		//TempNum1 = TempNum1 + j;
		//RowNoZero1 = RowNoZero1 - j;
		////����ֵ
		//if (this->col[this->row[i + 1] - 1] < add2.col[add2.row[i + 1] - 1])
		//{
		//	this->AddMemcpy((*result)->col + TempNum, add2.col + TempNum2, (*result)->num + TempNum, add2.num + TempNum2, RowNoZero2);
		//	TempNum = TempNum + RowNoZero2;
		//}
		//else if (this->col[this->row[i + 1] - 1] > add2.col[add2.row[i + 1] - 1])
		//{
		//	//std::memcpy((*result)->col + TempNum, this->col + TempNum1, sizeof(int) * RowNoZero1);
		//	//std::memcpy((*result)->num + TempNum, this->num + TempNum1 , sizeof(double) * RowNoZero1);
		//	this->AddMemcpy((*result)->col + TempNum, this->col + TempNum1, (*result)->num + TempNum, this->num + TempNum1, RowNoZero1);
		//	TempNum = TempNum + RowNoZero1;
		//}
		////
		result->row[i + 1] = TempNum;
	}
	result->allnum = TempNum;
	//�Ƿ���Ŀ���ʵ��,��Ҫ��֤��


}

/*******************���ؾ���i,j���ֵ*****************/
double TernArray::GetNum(int i, int j)
{
	int TempNowRow = row[i + 1] - row[i];
	int TempPreRow = row[i];
	if (TempNowRow < 300)
	{
		for (int k = TempPreRow; k < (TempPreRow + TempNowRow); k++)
		{
			if (j == col[k])
				return num[k];
		}
		return 0;
	}
	else
	{
		int RightInterval = TempNowRow - 1;
		int LeftInterval = 0;
		while (1)
		{
			int temp = floor((RightInterval - LeftInterval) / 2) + LeftInterval;
			if (j < col[TempPreRow + temp])
			{
				RightInterval = temp;
			}
			else if (j > col[TempPreRow + temp])
			{
				LeftInterval = temp;
			}
			else
			{
				return num[TempPreRow + temp];
			}

			if (RightInterval - LeftInterval == 1)
			{
				if (j == col[TempPreRow + LeftInterval])
					return num[TempPreRow + LeftInterval];
				else if (j == col[TempPreRow + RightInterval])
					return num[TempPreRow + LeftInterval];
				else
					return 0;
			}

		}
	}
}

/*******************����ת��*****************/
void TernArray::transpose(TernArray* TempArray)
{

	/* ����ʵ��*/
	if ((this->nmax != TempArray->colmax) || (this->colmax != TempArray->nmax)) {
		std::cout << "������������޷�ת��" << std::endl;
		return;
	}


	TempArray->row[0] = 0;

	WZmatrix** TempMatrix = new WZmatrix * [this->colmax];
	int temp = 0;

	//�������� ��θպ÷�����ֵ�����ڴ棬��ô��Ҫ�����ֵʱ�趨һ������Ч��,�����
	for (; temp < this->colmax; temp++) {
		TempMatrix[temp] = new WZmatrix(100);
	}

//#pragma omp parallel for num_threads(4)
	for (int i = 0; i < nmax; i++) 
	{
		for (int j = row[i]; j < row[i + 1]; j++) 
		{
			TempMatrix[this->col[j]]->OrderAdd(i, this->num[j]);
			TempArray->row[this->col[j] + 1]++;
		}
	}

//#pragma omp parallel for 
	for (int i = 0; i < this->colmax; i++) 
	{
		int err = TempMatrix[i]->Return(TempArray->col + TempArray->row[i], TempArray->num + TempArray->row[i], TempArray->row[i + 1]);
		if (!err) {
			std::cout << "ת�ô��� i=" << i << std::endl;
		}
		TempArray->row[i + 1] = TempArray->row[i + 1] + TempArray->row[i];
	}

	for (int i = 0; i < this->colmax; i++) {
		delete TempMatrix[i];
		TempMatrix[i] = NULL;
	}
		delete[] TempMatrix;
		TempMatrix = NULL;

	TempArray->colmax = this->nmax;
	TempArray->nmax = this->colmax;

	//this->Transpose = new TernArray(*TempArray);
	//TempArray->Transpose = new TernArray(*this);
}

/*******************�������*****************/
void TernArray::Multiply(TernArray* MultiplyOther, TernArray* ResultOther)
{

	if (this->colmax != MultiplyOther->nmax) {
		std::cout << "�޷��������" << std::endl;
		return;
	}

	if (this->nmax != ResultOther->nmax || MultiplyOther->colmax != ResultOther->colmax) {
		std::cout << "�������޷��洢" << std::endl;
		return;
	}

    /*csr*/
    std::map<int,double>* allnode = new std::map<int, double>  [this->nmax];
	
#pragma omp parallel for shared(allnode) 
	for (int i = 0; i < this->nmax; i++) 
	{
		for (int j = this->row[i]; j < this->row[i + 1]; j++)
		{
			for (int otherRow = MultiplyOther->row[this->col[j]]; \
				otherRow < MultiplyOther->row[this->col[j] + 1]; otherRow++)
			{
				int tempcol = MultiplyOther->col[otherRow];
				double tempnum= MultiplyOther->num[otherRow] * this->num[j];
                 #pragma omp atomic
				allnode[i][tempcol] += tempnum;
			}
		}
	}

	for (int i = 0; i < this->nmax; i++)
	{
		for (auto it : allnode[i])
		{
			ResultOther->OrderAssignment(i, it.first, it.second);
		}
	}

	delete[] allnode;
}



/*�������*/
double TernArray::MutiplyIJ(int tempi,int tempj, std::vector<std::pair<int, double>> Tij)
{
	int j = this->row[tempj];
	double re=0;
	int i = 0;

	while (i < Tij.size()&&j<this->row[tempj+1]) {
		if (Tij[i].first >= tempi || this->col[j] >= tempj)break;
		if (Tij[i].first == this->col[j]) {
			re = Tij[i].second * this->num[j] + re;
			i++;j++;
			continue;
		}
		if (Tij[i].first < this->col[j]) {
			i++;
		}
		else {
			j++;
		}
	}
	return re;
}



/*LU�ֽ������*/
int TernArray::SolveSPDMLU(double* pX, double* pB)
{

	//double** m_ppInitialwzMatrix = new double* [nmax];
	//or_mat(nmax, m_ppInitialwzMatrix);
	//creatLU(nmax, m_ppInitialwzMatrix);

	this->triu();
	/*���A�����ǶԳ���������*/
	TernArray tempL(nmax, colmax, nmax * colmax / 2);
	TernArray tempD(*(this->Dialog));

	for (int i = 1; i < nmax; i++) {
		std::vector<std::pair<int, double>> Tij;
		for (int j = 0; j < i; j++) {
			double temp = this->GetNum(i, j);
			temp = temp - tempL.MutiplyIJ(i, j, Tij);
			if (fabs(temp) > 1e-7) {
				std::pair<int, double> temp1(j, temp);
				Tij.push_back(temp1);
			}
		}

		double Aii = tempD.num[i];
		for (std::vector<std::pair<int, double>>::iterator it = Tij.begin(); it < Tij.end(); it++) {
			if (fabs(tempD.num[it->first]) > 1e-15) {
				double temp = it->second / tempD.num[it->first];
				tempL.OrderAssignment(i, it->first, temp);
				Aii = Aii - temp * it->second;
			}
			else return 0;
		}
		tempD.num[i] = Aii;
	}

    

	TernArray TransposeL(nmax, colmax, tempL.GetAllnum());
	tempL.transpose(&TransposeL);

	/*����*/
	//tempL.PrintfMatrix();
	//TransposeL.PrintfMatrix();
	//tempD.PrintfMatrix();



	memset(pX, 0, sizeof(double) * nmax);
	double* mat_y;
	mat_y = new double[nmax];
	memcpy(mat_y, pB, sizeof(double) * nmax);
	for (int i = 1; i < nmax; i++)
	{
		double tempYi=0.0;
		for (int j = tempL.row[i]; j < tempL.row[i + 1]; j++) {
			tempYi = tempYi + tempL.num[j] * mat_y[tempL.col[j]];
		}
	}
	pX[nmax - 1] = mat_y[nmax - 1] / tempD.num[nmax - 1];
	for (int i = nmax - 2; i >= 0; i--)
	{
//wz20210620
		double tempXi = 0.0;
		for (int j = TransposeL.row[i]; j < TransposeL.row[i + 1]; j++) {
			tempXi = tempXi + TransposeL.num[j] * mat_y[TransposeL.col[j]];
		}
		pX[i] = mat_y[i] / tempD.num[i] - tempXi;
	}

	delete[] mat_y;
	mat_y = NULL;
	return 1;
}