#include <memory>
#include<algorithm>
#include"WZAMG.h"

/*******************构造*****************/
AMG::AMG(const TernArray& temp, double* tempb,int tempcycle,int tempcoarsest,double tempSW_BOUND):CYCLES(tempcycle),COARSEST(tempcoarsest),SW_BOUND(tempSW_BOUND)
{
	Amatrix = new TernArray(temp);
	rowmax = Amatrix->GetNmax();
	colmax = Amatrix->GetColMax();
	b = new double[Amatrix->GetColMax()];
	A_matrix = new TernArray * [COARSEST];
	for (int i = 0; i < COARSEST; i++) {
		A_matrix[i] = NULL;
	}
	Wweight = new TernArray * [COARSEST];
	for (int i = 0; i < COARSEST; i++) {
		Wweight[i] = NULL;
	}
	TrangeWweight = new TernArray * [COARSEST];
	for (int i = 0; i < COARSEST; i++) {
		TrangeWweight[i] = NULL;
	}

	memcpy(b, tempb, sizeof(double) * colmax);

	A_matrix[0] = new TernArray(*Amatrix);

}

AMG::~AMG()
{
	delete Amatrix;
	Amatrix = NULL;

	delete[] b;
	b = NULL;

	for (int i = 0; i < COARSEST; i++) {
		if(A_matrix[i]!=NULL)
		delete A_matrix[i];
		A_matrix[i] = NULL;
	}
	A_matrix = NULL;

	for (int i = 0; i < COARSEST; i++) {
		if (Wweight[i] != NULL)
			delete Wweight[i];
		Wweight[i] = NULL;
	}
	Wweight = NULL;

	for (int i = 0; i < COARSEST; i++) {
		if (TrangeWweight[i] != NULL)
			delete TrangeWweight[i];
		TrangeWweight[i] = NULL;
	}
    TrangeWweight = NULL;

}

/*******************求解方程组*****************/
int AMG::UseEigen(int level,OneDimensionalArray* tempx, OneDimensionalArray* tempb)
{
	if (A != NULL) {
		delete A;
	}
	int row_A, col_A, row_b;
	int length = A_matrix[level]->GetNmax();
	col_A = length; row_A = length;
	row_b = row_A;
	A = new SparseMatrixType(length, length);
	Eigen::VectorXd x;
	Eigen::VectorXd b1;
	std::vector<T> tripletlist;

	//给向量b赋值;
	b1.resize(row_b);
	for (int i = 0; i < row_b; i++) {
		b1(i) = tempb->ReturnNum(i);
	}

	//给稀疏矩阵A赋值;
	//for (int i = 0; i < row_A; i++) {
	//	int tempLength = Amatrix->GetRowNoZero(i);
	//	if (!tempLength)continue;
	//	node* temp = new node[tempLength];
	//	for (int j = 0; j < tempLength; j++) {
	//		tripletlist.push_back(T(i, temp[j].col, temp[j].num));
	//	}
	//}
		for (int i = 0; i < row_A; i++){
		for (int j = A_matrix[level]->row[i]; j < A_matrix[level]->row[i+1]; j++){
			tripletlist.push_back(T(i, A_matrix[level]->col[j], A_matrix[level]->num[j]));
		}
	}
	A->setFromTriplets(tripletlist.begin(), tripletlist.end());
	A->makeCompressed();

	////test
	//SparseMatrixType B(*A);
	//SparseMatrixType c;
	//SparseMatrixType d;
	//d=B* c;

	//求解;
	Solve* p_A = new Solve(*A);
	x = p_A->solve(b1);
	for (int i = 0; i < row_b; i++) {
		tempx->Assignment(i, x(i));
	}

	return 1;
}

/*******************JACOBI方法光滑误差*****************/
double AMG::Smooth(int level, int max_iter, OneDimensionalArray* b, OneDimensionalArray* u_in, OneDimensionalArray* u_apx,double tempErr)
{

	if (A_matrix[level]->LowerTrangle == NULL) {
		A_matrix[level]->triu();
	}
	int length = A_matrix[level]->nmax;
	TernArray B(length,length,2*(A_matrix[level]->GetAllnum()));
	A_matrix[level]->UpperTrangle->add(&B, (A_matrix[level]->LowerTrangle));
	for (int i = 0; i < A_matrix[level]->nmax; i++) {
		for (int j = B.row[i]; j < B.row[i + 1]; j++) {
			if (A_matrix[level]->Dialog->num[i] > 1e-5) {
				B.num[j] = (double)-1 * (B.num[j] / (A_matrix[level]->Dialog->num[i]));
			}
		}
	}
	//   B.PrintfMatrix();
	OneDimensionalArray Multiplyb(length);
	for (int i = 0; i < length; i++) {
		if (A_matrix[level]->Dialog->num[i] > 1e-5) {
			double tempb = (b->ReturnNum(i) / (A_matrix[level]->Dialog->num[i]));
			Multiplyb.Assignment(i, tempb);
		}
	}

	//OneDimensionalArray TempX(*u_in);
	//vector<double>err;
	
	int i = 0;
	double err = 1;

	u_apx->copy(*u_in);
	while (i < max_iter&&err> tempErr) {
		OneDimensionalArray TempMultiplyB(length);
		OneDimensionalArray TempPlusX(length);
		B.MultiplyVector(u_apx, &TempMultiplyB);
		TempMultiplyB.add(&Multiplyb, &TempPlusX);
		u_apx->copy(TempPlusX);
		TempMultiplyB.~OneDimensionalArray();
		TempPlusX.~OneDimensionalArray();

		OneDimensionalArray TempErr01(length);
		OneDimensionalArray TempErr02(length);
		A_matrix[level]->MultiplyVector(u_apx, &TempErr01);
		TempErr01.reduce(b, &TempErr02);
		err = TempErr02.norm2();

		i++;
	}

	return err;

}



///*******************共轭梯度法求解方程组*****************/
//int	AMG::JPCG(TernArray* tempA, OneDimensionalArray* temp_b, int maxcycle, OneDimensionalArray& pX)
//{
//	//此处针对的是对称正定矩阵，所以A的转置=A
////最大循环次数maxcycle
//
//
//	int length = tempA->nmax;
//	OneDimensionalArray Xk(pX);
//	OneDimensionalArray Rk(length);
//	OneDimensionalArray Pk(length);
//	tempA->MultiplyVector(&Xk, &Rk);
//	temp_b->reduce(&Rk,&Pk);
//	Rk = Pk;
////	OneDimensionalArray Pk(Rk);
//	int k = 0;
//	double tempRR = Rk.norm2();
//
//
//	while (k<maxcycle) {
//		
//		if (tempRR < 1e-6) {
//			std::cout << "成功 k=" << k << std::endl;
//			break;
//		}
//		OneDimensionalArray tempApk(length);
//		OneDimensionalArray Xk2(length);
//		OneDimensionalArray Rk2(length);
//		OneDimensionalArray Pk2(length);
//		/*   1   */
//
//		tempA->MultiplyVector(&Pk, &tempApk);
//
//		double temp = Pk.Multiply(tempApk);
//
//		double ak = tempRR / temp;
//		/*   2   */
//		OneDimensionalArray a_p(Pk);
//	
//		a_p.MultiplyNum(ak);
//
//		Xk.add(&a_p, &Xk2);
//		/*   3   */
//		tempApk.MultiplyNum(ak);
//
//		Rk.reduce(&tempApk, &Rk2);
//		/*   4   */
//		double tempR2R2 = Rk2.MultiplyOneSelf();
//
//		double bk = tempR2R2 / tempRR;
//
//		/*    5    */
//		OneDimensionalArray b_p(Pk);
//
//		b_p.MultiplyNum(bk);
//
//		Rk2.add(&b_p, &Pk2);
//		
//		/*   6    */
//		Rk = Rk2;
//		Xk = Xk2;
//		Pk = Pk2;
//		tempRR = tempR2R2;
//
//		k++;
//	}
//
//	if (k == maxcycle) {
//		std::cout << "最后的共轭梯度计算不达标" << std::endl;
//		return 0;
//	}
//
//	pX = Xk;
//	return 1;
//
//}




//*******************共轭梯度法求解方程组*****************/

/*原始存档*/
int AMG::SolveEquationsCG(TernArray* tempA, int m_nFiledMaxItTime, int nMatrixNum, double* pX, double* pB)
{
	int i = 0;
	int nColumnNum = 0;
	int nColumn = 0;
	int j = 0;
	int nItTime = 0;					//迭代次数
	double dResult = 0.0;
	double	dBeta = 0.0;
	double	dAlpha = 0.0;
	double* pTemp = NULL;				//保存中间过程
	double* pTemp1 = NULL;
	double* pMargin = NULL;			//保存余量
	double* pVectorP = NULL;			//搜索向量
	double	dProduct[2] = { 1.0, 1.0 };	//保存前两次余量的内积
	pTemp = new double[nMatrixNum]();
	pTemp1 = new double[nMatrixNum]();
	pMargin = new double[nMatrixNum]();
	pVectorP = new double[nMatrixNum]();


	for (i = 0; i < nMatrixNum; i++)
	{
		nColumnNum = tempA->row[i + 1] - tempA->row[i];

		for (j = 0; j < nColumnNum; j++)
		{
			nColumn = tempA->row[i] + j;
			pTemp[i] += tempA->num[nColumn] * pX[tempA->col[nColumn]];
			if (tempA->col[nColumn] != i)
			{
				pTemp[tempA->col[nColumn]] += tempA->num[nColumn] * pX[i]; //对称矩阵
			}
		}

		pMargin[i] = pB[i] - pTemp[i];

		dResult += pMargin[i] * pMargin[i];
	}

	dProduct[0] = dResult;

	nItTime = 0;
	double err = 0;
	while (!IsFieldConverge(pMargin, nMatrixNum))
		//	while (err>1|| nItTime==0)
	{

		nItTime++;
		//		PostMessage(g_hMTDEWnd, WM_EOS_PROG_NUM_DOWN, nItTime, (LONG)g_hSelfWnd);
		if (nItTime == 1)
		{
			/*test 添加Qz=r p0=z0*/

			//////////////////////
			for (i = 0; i < nMatrixNum; i++)
			{
				pVectorP[i] = pMargin[i];
			}
		}
		else if (nItTime > m_nFiledMaxItTime)
		{
			//SetEosMessage(0x2020020A);
			//SendOutputText(g_hWnd, pDoc, "！！！达到场求解的最大迭代次数，结果可能不收敛，请增大最大迭代次数", TEXT_OUTPUT_OK);
			std::cout << "达到场求解的最大迭代次数" << std::endl;
			return 0;
		}
		else
		{
			dBeta = dProduct[0] / dProduct[1];
			for (i = 0; i < nMatrixNum; i++)
			{
				pTemp[i] = dBeta * pVectorP[i];
				pVectorP[i] = pMargin[i] + pTemp[i];
			}
		}

		dProduct[1] = dProduct[0];
		double dResult_2 = 0.0;

		memset(pTemp, 0, sizeof(double) * nMatrixNum);//使pTemp内的元素重置，因为此处原是调用矩阵乘向量函数，此函数重置了pTemp.
		//以下循环代替矩阵乘向量函数。
		for (i = 0; i < nMatrixNum; i++)
		{
			nColumnNum = tempA->row[i + 1] - tempA->row[i];

			for (j = 0; j < nColumnNum; j++)
			{
				nColumn = tempA->row[i] + j;
				pTemp[i] += tempA->num[nColumn] * pVectorP[tempA->col[nColumn]];
				if (tempA->col[nColumn] != i)
				{
					pTemp[tempA->col[nColumn]] += tempA->num[nColumn] * pVectorP[i]; //对称矩阵
				}
			}
			int iii = 0;

		}

		double temp_num = 0;
		for (int k = 0; k < nMatrixNum; k++) {
			double PV = pVectorP[k];
			double PT = pTemp[k];
			temp_num = temp_num + PV * PT;
		}

		dAlpha = dProduct[0] / temp_num;

		for (i = 0; i < nMatrixNum; i++)
		{
			//这个循环代替数乘向量函数调用。
			//pTemp[i] = dAlpha * pTemp[i];

			//ComputeVectorSubVector(pMargin, pTemp, m_nMatrixNum, pMargin);
			//这个循环代替向量减向量函数的调用。
			pMargin[i] = pMargin[i] - dAlpha * pTemp[i];

			//原来就是这个循环，没有函数调用。
			//pTemp[i] = pVectorP[i];

			//ComputeNumberMulVector(dAlpha, pTemp, m_nMatrixNum); 
			//下面循环代替数乘向量函数调用。
			//pTemp[i] = dAlpha * pVectorP[i];

			//ComputeVectorAddVector(pX, pTemp, m_nMatrixNum, pX);
			//下面循环代替向量加向量函数调用。
			pX[i] = pX[i] + dAlpha * pVectorP[i];

			//dProduct[0] = ComputeVectorMulVector(pMargin, pMargin, m_nMatrixNum);
			//下面循环代替向量乘向量函数调用。


			dResult_2 += pMargin[i] * pMargin[i];
		}
		dProduct[0] = dResult_2;


	}

	delete[]pTemp;
	pTemp = NULL;
	delete[]pMargin;
	pMargin = NULL;
	delete[]pVectorP;
	pVectorP = NULL;

	return nItTime;
}
/*************/
//
//int AMG::SolveEquationsCG(TernArray* tempA,int m_nFiledMaxItTime, int nMatrixNum, double* pX, double* pB)
//{
//	int i = 0;
//	int nColumnNum = 0;
//	int nColumn = 0;
//	int j = 0;
//	int nItTime = 0;					//迭代次数
//	double dResult = 0.0;
//	double	dBeta = 0.0;
//	double	dAlpha = 0.0;
//	double* pTemp = NULL;				//保存中间过程
//	double* pTemp1 = NULL;
//	double* pMargin = NULL;			//保存余量
//	double* pVectorP = NULL;			//搜索向量
//	double* pZ = NULL;                  //R预处理后的Z
//	double	dProduct[2] = { 1.0, 1.0 };	//保存前两次余量的内积
//	pTemp = new double[nMatrixNum]();
//	pTemp1 = new double[nMatrixNum]();
//	pMargin = new double[nMatrixNum]();
//	pVectorP = new double[nMatrixNum]();
//	pZ= new double[nMatrixNum]();
//
//	for (i = 0; i < nMatrixNum; i++)
//	{
//		nColumnNum = tempA->row[i + 1] - tempA->row[i];
//
//		for (j = 0; j < nColumnNum; j++)
//		{
//			nColumn = tempA->row[i] + j;
//			pTemp[i] += tempA->num[nColumn] * pX[tempA->col[nColumn]];
//			if (tempA->col[nColumn] != i)
//			{
//				pTemp[tempA->col[nColumn]] += tempA->num[nColumn] * pX[i]; //对称矩阵
//			}
//		}
//
//		pMargin[i] = pB[i] - pTemp[i];
//
//	}
//	 
//	dResult=this->proCGmethode(pMargin, pZ);
//
//	dProduct[0] = dResult;
//
//	nItTime = 0;
//	double err=0;
//	while (!IsFieldConverge(pMargin, nMatrixNum))
//	{
//
//		nItTime++;
//		if (nItTime == 1)
//		{
//			for (i = 0; i < nMatrixNum; i++)
//			{
////				pVectorP[i] = pMargin[i];
//				pVectorP[i] = pZ[i];
//			}
//		}
//		else if (nItTime > m_nFiledMaxItTime)
//		{
//			std::cout << "达到场求解的最大迭代次数" << std::endl;
//			return 0;
//		}
//		else
//		{
//			dBeta = dProduct[0] / dProduct[1];
//			for (i = 0; i < nMatrixNum; i++)
//			{
//				pTemp[i] = dBeta * pVectorP[i];
//				//pVectorP[i] = pMargin[i] + pTemp[i];
//				pVectorP[i] = pZ[i] + pTemp[i];
//			}
//		}
//
//		dProduct[1] = dProduct[0];
//		double dResult_2 = 0.0;
//
//		memset(pTemp, 0, sizeof(double) * nMatrixNum);//使pTemp内的元素重置，因为此处原是调用矩阵乘向量函数，此函数重置了pTemp.
//		//以下循环代替矩阵乘向量函数。
//		for (i = 0; i < nMatrixNum; i++)
//		{
//			nColumnNum = tempA->row[i + 1] - tempA->row[i];
//
//			for (j = 0; j < nColumnNum; j++)
//			{
//				nColumn = tempA->row[i] + j;
//				pTemp[i] += tempA->num[nColumn] * pVectorP[tempA->col[nColumn]];
//				if (tempA->col[nColumn] != i)
//				{
//					pTemp[tempA->col[nColumn]] += tempA->num[nColumn] * pVectorP[i]; //对称矩阵
//				}
//			}
//			int iii = 0;
//
//		}
//		
//		double temp_num = 0;
//		for (int k = 0; k < nMatrixNum; k++) {
//			double PV = pVectorP[k];
//			double PT= pTemp[k];
//			temp_num = temp_num +PV*PT;
//		}
//
//		dAlpha = dProduct[0] / temp_num;
//
//		for (i = 0; i < nMatrixNum; i++)
//		{
//			pMargin[i] = pMargin[i] - dAlpha * pTemp[i];
//			pX[i] = pX[i] + dAlpha * pVectorP[i];
//		}
//		dResult_2= this->proCGmethode(pMargin, pZ);
//		dProduct[0] = dResult_2;
//		
//	}
//
//	delete[]pTemp;
//	pTemp = NULL;
//	delete[]pMargin;
//	pMargin = NULL;
//	delete[]pVectorP;
//	pVectorP = NULL;
//	delete[] pZ;
//	pZ = NULL;
//    	return nItTime;
//}

/*CG法中的预处理方式*/
double AMG::proCGmethode(double* R,double *X)
{
	OneDimensionalArray tempb(this->rowmax);
	OneDimensionalArray uin(this->rowmax);
	OneDimensionalArray uout(this->rowmax);
	tempb.Initial(this->rowmax, R);
	this->amg_cycle(0,&tempb, &uin, &uout);
	uout.InitialB(X);
	double ret=uout.Multiply(tempb);
	return ret;
}
/*************************************/

int AMG::IsFieldConverge(double* temp, int nmax)
{
	for (int i = 0; i < nmax; i++) {
		if (fabs(temp[i]) > 1.e-15) {
			return false;
		}
	}
	return true;
}

/*******************求解网格的强连接点*****************/
void AMG::my_strong(int level, double local_SW_BOUND)
{
	int length = A_matrix[level]->GetNmax();

	/*生成稀疏矩阵 列项量数组和行向量数组、一个存每行对应值的列向量一个存每行对应值的非零元个数 */
	int** tempST = new int* [A_matrix[level]->GetNmax()];

	for (int i = 0; i < length; i++) {
		int templeng = A_matrix[level]->GetRowNoZero(i) + 1;
		tempST[i] = new int[templeng];
		memset(tempST[i], 0, sizeof(int) * templeng);
	}
	int* rowarray = new int[length];

	/*对应值非零元总数*/
	int tempallnum = 0;
	Wweight[level] = new TernArray(length, length, 2 * length);
	/*test2*/
    Wweight[level] = new TernArray(length, length, 2 * length);
    std::vector<std::pair<int,int>>templambda(length);
	int maxTrans = 100;
	int** trans = new int*[length];
	for (int i =0 ; i < length; i++) {
		trans[i] = new int[maxTrans];
		memset(trans[i], 0, sizeof(int) * maxTrans);
	}
	int* lambdaPos = new int[length];
	int* lambda = new int[length];
	int* v_set = new int[length];
	memset(lambda, 0, sizeof(int) * length);
	memset(v_set, 0, sizeof(int)* length);
	for (int i = 0; i < length; i++) {
		int rownum = A_matrix[level]->row[i];
		double mini = A_matrix[level]->num[rownum];
		if (A_matrix[level]->col[rownum] == i) {
			mini = A_matrix[level]->num[rownum + 1];
		}
		for (int j = rownum; j < A_matrix[level]->row[i + 1]; j++) {
			if (A_matrix[level]->col[j] != i) {
				mini = A_matrix[level]->num[j] > mini ? mini : A_matrix[level]->num[j];
			}
		}
		int time = 0;
		for (int j = rownum; j < A_matrix[level]->row[i + 1]; j++) {
			if (A_matrix[level]->col[j] != i) {
				double compare = A_matrix[level]->num[j];
				int tempcol = A_matrix[level]->col[j];
				if (compare < mini * local_SW_BOUND) {
					lambda[tempcol]++;
					tempST[i][time] = tempcol;
					time++;
				}
			}
		}
		rowarray[i] = time;
		tempallnum = tempallnum + time;
	}
	/*生成粗网格*/
	std::pair<int, int>temppair;
	for (int i = 0; i < length; i++) {
		temppair.first = lambda[i];
		temppair.second = i;
		templambda[i]=temppair;
	}
	std::sort(templambda.begin(), templambda.end(), [](std::pair<int, int>a, std::pair<int, int>b) {return a.first < b.first; });
	int tempn = 0;
	for (std::vector<std::pair<int, int>>::iterator it = templambda.begin(); it != templambda.end(); it++)
	{
		lambda[tempn] = it->second;
		lambdaPos[it->second] = tempn;
		tempn++;
	}
	int nodenumber = 0;
	int nG = 0;
	//如果set花费时间很长
	//int* inNode = new int[length];
	//inNode[0] = 0;
	for (int i = 0; i < length; i++) 
	{
		std::set<int> setNode;
		if (nodenumber == length)break;
		if (lambda[i]== -1)continue;
		for (int j = 0; j < rowarray[lambda[i]]; j++)
		{
			if (v_set[tempST[lambda[i]][j]] != -1) {
				setNode.insert(tempST[lambda[i]][j]);
			}
		}	
		if (v_set[lambda[i]] != -1) {
			setNode.insert(lambda[i]);
		}
		for (std::set<int>::iterator it = setNode.begin(); it != setNode.end(); it++) {
			Wweight[level]->OrderAssignment(nG, *it, 1);
			v_set[*it] = -1;
			lambda[lambdaPos[*it]] = -1;
			nodenumber++;
			trans[*it][0]++;
			if (trans[*it][0]==maxTrans) {
				int* temp = new int[2*maxTrans];
				memcpy(temp, trans[*it], sizeof(int) * maxTrans);
				maxTrans = 2 * maxTrans;
			}
			trans[*it][trans[*it][0]] = nG;
		}
		if (setNode.size() != 0) {
			nG++;
		}
	}
	Wweight[level]->ChangeRow(nG);

	TrangeWweight[level] = new TernArray(length, nG, Wweight[level]->GetAllnum());
	for (int i = 0; i < length; i++) 
	{
		for (int j = 1; j <= trans[i][0]; j++) 
		{
			TrangeWweight[level]->OrderAssignment(i,trans[i][j],1);
		}
	}

	for (int i = 0; i < length; i++) {
		delete[] trans[i];
		trans[i]=NULL;
	}
	delete[] trans;
	trans = NULL;
	delete[] lambda;
	lambda = NULL;
	delete[] v_set;
	v_set = NULL;
	delete[] lambdaPos;
	lambdaPos = NULL;
//	/*test*/
//	int length = A_matrix[level]->GetNmax();
//	std::vector<int> minS(A_matrix[level]->GetNmax());
//	std::vector<int>* tempS = new std::vector<int>[A_matrix[level]->GetNmax()];
//	Wweight[level] = new TernArray(length, length, 2 * length);
//
//	/*取出每行的最小值min，判断该行其他值<min*local_SW_BOUND,满足条件的点稀疏矩阵tempST相应位置=1   */
//	for (int i = 0; i < A_matrix[level]->nmax; i++) {
//		double mini = 0;
//		for (int j = A_matrix[level]->row[i]; j < A_matrix[level]->row[i + 1]; j++) {
//			if (A_matrix[level]->col[j]!=i) {
//				mini = A_matrix[level]->num[j] > mini ? mini : A_matrix[level]->num[j];
//			}
//		}
////		mini = -1 * fabs(mini);
//		for (int j = A_matrix[level]->row[i]; j < A_matrix[level]->row[i + 1]; j++) {
//			double compare = A_matrix[level]->num[j];
//			int tempcol = A_matrix[level]->col[j];
//			if (compare < mini * local_SW_BOUND) {
//				tempS[i].push_back(tempcol);
//				minS[tempcol]++;
//			}
//		}
//	}
//
////	std::vector<std::set<int>>trans(length);
//	int* G1node = new int[2 * length];
//	memset(G1node, 0, sizeof(int) * 2 * length);
//	std::vector<int>U(A_matrix[level]->GetNmax(),1);
//	int alltime = A_matrix[level]->GetNmax();
//	
//	int nc = 0;
//
//	while (alltime > 0) {
//	  //求i点
//		int minValue = std::distance(minS.begin(), min_element(minS.begin(), minS.end()));
//		if (minS[minValue] == 0) {
//			U[minValue] = -1;
//			minS[minValue] = 2 * A_matrix[level]->GetNmax();
//		}
//		int j = -1;
//		double minj = 1.0e5;
//		for (int i = A_matrix[level]->row[minValue]; i < A_matrix[level]->row[minValue + 1]; i++) {
//			U[A_matrix[level]->col[i]];
//			if (U[A_matrix[level]->col[i]] != -1&& A_matrix[level]->col[i]!=minValue) {
//				minj = minj < A_matrix[level]->num[i] ? minj : A_matrix[level]->num[i];
//				j = A_matrix[level]->col[i];
//			}
//		}
//		std::vector<int>::iterator it=std::find(tempS[minValue].begin(), tempS[minValue].end(), j);
//		if (it != tempS[minValue].end()) {
//			if (j < minValue) {
//				/////////algorithm1
//				Wweight[level]->OrderAssignment(nc, j, 1);
//				Wweight[level]->OrderAssignment(nc, minValue, 1);
//			}
//			else {
//				//////////algorithm1
//				Wweight[level]->OrderAssignment(nc, minValue, 1);
//				Wweight[level]->OrderAssignment(nc, j, 1);
//			}
//			U[j] = -1;
//			minS[j] =2*A_matrix[level]->GetNmax();
//			alltime--;
//			for (int l = 0; l < tempS[j].size(); l++) {
//				minS[tempS[j][l]]--;
//			}
//		}
//		//////////algorithm1
//		else {
//			Wweight[level]->OrderAssignment(nc, minValue, 1);
//		}
//		alltime--;
//		U[minValue] = -1;
//		minS[minValue] = A_matrix[level]->GetNmax();
//		for (int l = 0; l < tempS[minValue].size(); l++) {
//			minS[tempS[minValue][l]]--;
//		}
//		nc++;
//	}
//	////////algorithm1
//	Wweight[level]->ChangeRow(nc);
//
//	//////////////////////

}

/*algorithm1*/
void AMG::amg_get_Gnode(int& nc, int* Gnode,TernArray* A)
{
	/*test*/
	int local_SW_BOUND = 0.25;
	int length = A->GetNmax();
	std::vector<int> minS(length);
	std::vector<int>* tempS = new std::vector<int>[length];

	/*取出每行的最小值min，判断该行其他值<min*local_SW_BOUND,满足条件的点稀疏矩阵tempST相应位置=1   */
	for (int i = 0; i < A->nmax; i++) {
		double mini = 0;
		for (int j = A->row[i]; j < A->row[i + 1]; j++) {
			if (A->col[j] != i) {
				mini = A->num[j] > mini ? mini : A->num[j];
			}
		}
		for (int j = A->row[i]; j < A->row[i + 1]; j++) {
			double compare = A->num[j];
			int tempcol = A->col[j];
			if (compare < mini * local_SW_BOUND) {
				tempS[i].push_back(tempcol);
				minS[tempcol]++;
			}
		}
	}

	std::vector<int>U(length, 1);
	int alltime = length;

	nc = 0;

	while (alltime > 0) {
		//求i点
		int minValue = std::distance(minS.begin(), min_element(minS.begin(), minS.end()));
		if (minS[minValue] == 0) {
			U[minValue] = -1;
			minS[minValue] = 2 * length;
		}
		int j = -1;
		double minj = 1.0e5;
		for (int i = A->row[minValue]; i < A->row[minValue + 1]; i++) {
			U[A->col[i]];
			if (U[A->col[i]] != -1 && A->col[i] != minValue) {
				minj = minj < A->num[i] ? minj : A->num[i];
				j = A->col[i];
			}
		}
		std::vector<int>::iterator it = std::find(tempS[minValue].begin(), tempS[minValue].end(), j);
		if (it != tempS[minValue].end()) {
			if (j < minValue) {
				Gnode[2 * nc] = j;
				Gnode[2 * nc + 1] = minValue;
			}
			else {
				Gnode[2 * nc] = minValue;
				Gnode[2 * nc + 1] = j;
			}
			U[j] = -1;
			minS[j] = 2 * length;
			alltime--;
			for (int l = 0; l < tempS[j].size(); l++) {
				minS[tempS[j][l]]--;
			}
		}
		else {
			Gnode[2 * nc] = minValue;
			Gnode[2 * nc + 1] = -1;
		}
		alltime--;
		U[minValue] = -1;
		minS[minValue] = length;
		for (int l = 0; l < tempS[minValue].size(); l++) {
			minS[tempS[minValue][l]]--;
		}
		nc++;
	}

	//////////////////////
}



/*******************algorithm2*****************/
void AMG::amg_get_W(int level)
{
	int length=A_matrix[level]->GetNmax();
	TernArray* A = A_matrix[level];
	int* G1node = new int[2 * length]; 
	memset(G1node, 0, sizeof(int) * 2 * length);
	int nc = 0;
	this->amg_get_Gnode(nc, G1node, A);
	TernArray A1(nc, nc, 100 * nc);
	std::vector<std::pair<int,double>>* A1node = new std::vector<std::pair<int,double>>[nc];
#pragma omp  parallel for 
	for (int i = 0; i < nc; i++) 
	{
		for (int j = 0; j < nc; j++) 
		{
			double tempnum = A->GetNum(G1node[2*i],G1node[2*j]);
			if (G1node[2 * i + 1] != -1) {
			   tempnum+= A->GetNum(G1node[2 * i+1], G1node[2 * j]);
			}
			if (G1node[2 * j + 1] != -1) {
				tempnum += A->GetNum(G1node[2 * i], G1node[2 * j+1]);
				if (G1node[2 * i + 1] != -1) {
					tempnum += A->GetNum(G1node[2 * i + 1], G1node[2 * j + 1]);
				}
			}

			if (!tempnum) {
				std::pair<int, double>temp(j, tempnum);
				A1node[i].push_back(temp);
			}
		}
	}

	for (int i = 0; i < nc; i++) 
	{
		std::sort(A1node[i].begin(), A1node[i].end(), [](std::pair<int, double>a, std::pair<int, double>b) {return a.first < b.first; });
	}

	for (int i = 0; i < nc; i++) 
	{
		for (std::vector<std::pair<int, double>>::iterator it = A1node[i].begin();it != A1node[i].end(); it++)
		{ 
			A1.OrderAssignment(i, it->first, it->second);
		}
	}

	int nc1 = 0;
	int* G2node = new int[2 * nc];
	this->amg_get_Gnode(nc1, G2node, &A1);
	Wweight[level] = new TernArray(nc1, nc1, 4 * nc1);
	int** wnode = new int* [nc1];
	for (int i = 0; i < 5; i++) 
	{
		wnode[i] = new int[5];
		memset(wnode[i], 0, sizeof(int) * 5);
	}

#pragma omp parallel for 
	for (int i = 0; i < nc1; i++) 
	{
		std::set<int>tempnode;
		int tempj1 = G2node[2 * i];
		int tempj2 = G2node[2 * i];
		tempnode.insert(G1node[2 * tempj1]);
		if (G1node[2 * tempj1 + 1] != -1) {
			tempnode.insert(G1node[2 * tempj1 + 1]);
		}
		if (tempj2 != -1) {
			tempnode.insert(G1node[2 * tempj2]);
			if (G1node[2 * tempj2 + 1] != -1) {
				tempnode.insert(G1node[2 * tempj2 + 1]);
			}
		}

		int temp = 0;
		for (std::set<int>::iterator it = tempnode.begin(); it != tempnode.end(); it++) {
			temp++;
			wnode[i][temp] = *it;
		}
		wnode[i][0] = temp;
	}

	for (int i = 0; i < nc1; i++)
	{
		for (int j = 0; j < wnode[i][0]; j++)
		{
			Wweight[level]->OrderAssignment(i, wnode[i][j + 1], 1);
		}
	}

	for (int i = 0; i < nc1; i++)
	{
		delete wnode[i];
		wnode[i] = NULL;
	}
	delete wnode;
	wnode = NULL;
	delete A1node;
	A1node = NULL;
	delete[] G1node;
	G1node = NULL;
}

/*******************循环的到足够的粗网格*****************/
void AMG::amg_setup_wz(int level)
{
	/*递归的构造粗网格*/
	if (level != (COARSEST-1))
	{
		my_AMG(level);
		amg_setup_wz(level+1);
	}
}


/*******************求解粗网格 插值算子 光滑算子*****************/
void AMG::my_AMG(int level)
{
	/*构造粗网格*/
	this->my_strong(level, SW_BOUND);

	/*  求A(LEVEL+1)=W'*A(LEVEL)*W  */
	TernArray TempTernArray(Wweight[level]->GetNmax(), Wweight[level]->GetColMax(), Wweight[level]->GetAllnum());
	Wweight[level]->Multiply(A_matrix[level], &TempTernArray);
	A_matrix[level + 1] = new TernArray(Wweight[level]->GetNmax(), Wweight[level]->GetNmax(), TempTernArray.GetAllnum());
	TempTernArray.Multiply(TrangeWweight[level], A_matrix[level+1]);


}

/*******************在粗网格上求解方程组 并返回直到最细的网格上*****************/
double AMG::amg_cycle(int level, OneDimensionalArray* tempb, OneDimensionalArray* u_in, OneDimensionalArray* u_out)
{
	
	double Err;
	if (level == (COARSEST-1)) {

		this->UseEigen(level,u_out,tempb);
		
	}
	else {
		/*u_apx是 x=Bx+f在40次JACOBI迭代后的返回值*/

		OneDimensionalArray u_apx(u_in->ReturnMax());
		double Ierr = this->Smooth(level,1, tempb, u_in, &u_apx,1e-15);

		/*resid=b-Ax*/
		OneDimensionalArray resid(tempb->ReturnMax());
		//A_matrix[level]->MultiplyVector(&u_apx, &tempr);
		//tempb->reduce(&tempr, &resid);
		double RH1= this->NormErr(A_matrix[level], &u_apx, tempb, &resid);
//		double z = resid.norm2();

		/*apx_coarse=W*u_apx  b_coarse=W*resid */
		OneDimensionalArray apx_coarse(Wweight[level]->GetNmax());
		Wweight[level]->MultiplyVector(&u_apx, &apx_coarse);

		OneDimensionalArray b_coarse(Wweight[level]->GetNmax());
		Wweight[level]->MultiplyVector(&resid, &b_coarse);

		/*递归*/
		OneDimensionalArray u_coarse(Wweight[level]->GetNmax());
		amg_cycle( level + 1, &b_coarse, &apx_coarse, &u_coarse);

		/*correct=TransposeW*u_coarse*/
		OneDimensionalArray correct(TrangeWweight[level]->GetNmax());
		TrangeWweight[level]->MultiplyVector(&u_coarse, &correct); 

		/*u_out=u_apx+correct*/
		u_apx.add(&correct);

//		Err = (this->NormErr(A_matrix[level], &u_apx, tempb, u_out));
		Err = this->Smooth(level,1, tempb, &u_apx, u_out, 1e-15);

		return Err;

		/*u_apx是 x=Bx+f在40次JACOBI迭代后的返回值*/

		//OneDimensionalArray u_apx01(u_in->ReturnMax());
		//A_matrix[level]->MultiplyVector(u_in, &u_apx01);
		//tempb->reduce(&u_apx01, false);

		//OneDimensionalArray u_apz(u_in->ReturnMax());
		//OneDimensionalArray u_apx02(u_in->ReturnMax());
		//double Ierr = this->Smooth(level, 100, &u_apx01, &u_apx02, &u_apz, 1e-10);

		//
		//OneDimensionalArray Rk(u_apx01);
		//A_matrix[level]->MultiplyVector(&u_apz, &u_apx02);
		//Rk.reduce(&u_apx02,true);

		//OneDimensionalArray Rkother(Wweight[level]->GetNmax());
		//Wweight[level]->MultiplyVector(&Rk, &Rkother);

		///*递归*/
		//OneDimensionalArray u_coarse(Wweight[level]->GetNmax());
		//OneDimensionalArray u_apxother(Wweight[level]->GetNmax());
		//amg_cycle(level + 1, &Rkother, &u_apxother, &u_coarse);

		///*correct=TransposeW*u_coarse*/
		//OneDimensionalArray u_zother(TrangeWweight[level]->GetNmax());
		//OneDimensionalArray u_zother3(TrangeWweight[level]->GetNmax());
		//OneDimensionalArray temp(TrangeWweight[level]->GetNmax());
		//TrangeWweight[level]->MultiplyVector(&u_coarse, &u_zother);
		//A_matrix[level]->MultiplyVector(&u_zother, &temp);
		//Rkother.reduce(&temp, false);	
		//Rkother.Free();

		//this->Smooth(level, 100, &temp,&Rkother, &u_zother3);

		///*u_out=u_apx+correct*/
		//u_apx.add(&correct, u_out);


		//OneDimensionalArray u_temp(u_in->ReturnMax());
		//Err = (this->NormErr(A_matrix[level], u_out, tempb, &u_temp));
	}
}

//容器vector中元素的去重
//std::vector<int> AMG::unique_element_in_vector(std::vector<int> v)
//{
//	std::vector<int>::iterator vector_iterator;
//	sort(v.begin(), v.end());
//	vector_iterator = unique(v.begin(), v.end());
//	if (vector_iterator != v.end()) {
//		v.erase(vector_iterator, v.end());
//	}
//	return v;
//}

////两个std::vector求交集
//void AMG::vectors_intersection(std::vector<int>& v1, std::vector<int>& v2, std::vector<int>& v)
//{
//	sort(v1.begin(), v1.end());
//	sort(v2.begin(), v2.end());
//	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v));//求交集 
//}
//
////两个std::vector求并集
//void AMG::vectors_set_union(std::vector<int>& v1, std::vector<int>& v2, std::vector<int>& v)
//{
//	sort(v1.begin(), v1.end());
//	sort(v2.begin(), v2.end());
//	set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v));//求交集 
//}
//
////两个std::vector求差集
//void AMG::vectors_set_different(std::vector<int>& v1, std::vector<int>& v2)
//{
//	sort(v1.begin(), v1.end());
//	sort(v2.begin(), v2.end());
//	std::vector<int> ivec(v1.size());
//	auto iter = set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), ivec.begin());	//ivec为：2,3,5,6,7
//	ivec.resize(iter - ivec.begin());//重新确定ivec大小
//	v1.swap(ivec);
//}
//
////判断std::vector的某一元素是否存在
//bool AMG::is_element_in_vector(std::vector<int>& v, int element)
//{
//	std::vector<int>::iterator it;
//	it = find(v.begin(), v.end(), element);
//	if (it != v.end()) {
//		return true;
//	}
//	else {
//		return false;
//	}
//}

/*******************输出误差***************************/
double AMG::NormErr(TernArray* tempA, OneDimensionalArray* x, OneDimensionalArray* b, OneDimensionalArray* tempu_out=NULL)
{
	if (tempu_out == NULL) {
		tempu_out = new OneDimensionalArray(tempA->GetNmax());
		tempA->MultiplyVector(x, tempu_out);
		tempu_out->reduce(b,true);
		double err001 = tempu_out->norm2();
		delete tempu_out;
		tempu_out = NULL;
		return err001;
	}
	tempA->MultiplyVector(x, tempu_out);
	tempu_out->reduce(b,true);
	double err001 = tempu_out->norm2();
	return err001;
}

/*******************AMG的主程序*****************/
int AMG::MainReSult(OneDimensionalArray* x,int v_max)
{
	/*0阶A矩阵*/
	int level = 0;

	/*构造需要的x、b和中间变量*/
	OneDimensionalArray FirstX(colmax);
	OneDimensionalArray LastX(colmax);
	OneDimensionalArray Lastb(rowmax);
	Lastb.Initial(this->rowmax,this->b);

	/*构造粗网格矩阵和插值矩阵、光滑算子*/
	this->amg_setup_wz(level);

	/*由上一部生成的矩阵*/
	int time=0;
	double Err = 1;
	v_max = 4;
//	Err = this->Smooth(0, 1000, &Lastb, &FirstX, &LastX, 1.0e-15);
	while (time < v_max&&Err>1.0e-15) {
	    Err=this->amg_cycle(level, &Lastb, &FirstX, &LastX);
		FirstX = LastX;
		LastX.Free();
		time++;
	}

	if (Err > 1.0e-15) {
		std::vector<double>err;
		err.push_back(this->NormErr(Amatrix, &LastX, &Lastb));

		/*共轭梯度法*/
		double* tempx = NULL;
		LastX.ReverseB(&tempx);
		int chose = this->SolveEquationsCG(Amatrix, 1e6, rowmax, tempx, this->b);
		OneDimensionalArray tempX(rowmax);
		tempX.Initial(rowmax, tempx);
		err.push_back(this->NormErr(Amatrix, &tempX, &Lastb));

		/*判断是否成功*/
		if (!chose) {

			std::cout << "共轭梯度法失败" << std::endl;
			return 0;
		}
		*x = tempX;
		return chose;
	}
	return 0;
}




