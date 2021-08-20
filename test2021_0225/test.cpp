#include"WZAMG.h"
#include<numeric>
#include<algorithm>
#include<fstream>
#include<string>
#include<time.h>
#include<windows.h>
#include <psapi.h>
#pragma comment(lib,"psapi.lib")
#include"WZAMG.h"



void getAmatrix(int max, int CYCLES,int COARSEST,double SW_BOUND)
{
	clock_t start, end;
	start = clock();
	/*s三元法记录压缩矩阵- -先定义指针*/
	int nmax = max;
	int colmax = max;
	int temp1 = max + 1;
	int temp2 = max * 20;
	int* row;
	int* col;
	double* num;
	double* b;
	/*四个文件对应四个输入流的文件名*/
	const char* filename1=NULL;
	const char* filename2=NULL;
    const char* filename3=NULL;
	const char* filename4=NULL;

	double* WZ_Voltage = new double[max];
	memset(WZ_Voltage, 0, sizeof(double) * max);
	
    /*根据矩阵阶数选择对应的文件地址*/
	switch (max)
	{
	case 28109: {
		filename1 = "C:\\Users\\wangzhi\\Desktop\\matrix\\newmatrix\\matrix\\row.txt";
		filename2 = "C:\\Users\\wangzhi\\Desktop\\matrix\\newmatrix\\matrix\\col.txt";
		filename3 = "C:\\Users\\wangzhi\\Desktop\\matrix\\newmatrix\\matrix\\num.txt";
		filename4 = "C:\\Users\\wangzhi\\Desktop\\matrix\\newmatrix\\matrix\\b.txt";
		break; 
	}
	case 3: {
		filename1 = "C:\\Users\\wangzhi\\Desktop\\matrix\\lu2\\row.txt";
		filename2 = "C:\\Users\\wangzhi\\Desktop\\matrix\\lu2\\col.txt";
		filename3 = "C:\\Users\\wangzhi\\Desktop\\matrix\\lu2\\num.txt";
		filename4 = "C:\\Users\\wangzhi\\Desktop\\matrix\\lu2\\b.txt";
		break;
	}
	case 4: {
		filename1 = "C:\\Users\\wangzhi\\Desktop\\matrix\\lu\\row.txt";
		filename2 = "C:\\Users\\wangzhi\\Desktop\\matrix\\lu\\col.txt";
		filename3 = "C:\\Users\\wangzhi\\Desktop\\matrix\\lu\\num.txt";
		filename4 = "C:\\Users\\wangzhi\\Desktop\\matrix\\lu\\b.txt";
		break;
	}
	case 9: {
		filename1 = "C:\\Users\\wangzhi\\Desktop\\matrix\\matrix2\\row.txt";
		filename2 = "C:\\Users\\wangzhi\\Desktop\\matrix\\matrix2\\col.txt";
		filename3 = "C:\\Users\\wangzhi\\Desktop\\matrix\\matrix2\\num.txt";
		filename4 = "C:\\Users\\wangzhi\\Desktop\\matrix\\matrix2\\b.txt";
		break;
	}
	case 36: {
		filename1 = "C:\\Users\\wangzhi\\Desktop\\matrix\\new2\\row.txt";
		filename2 = "C:\\Users\\wangzhi\\Desktop\\matrix\\new2\\col.txt";
		filename3 = "C:\\Users\\wangzhi\\Desktop\\matrix\\new2\\num.txt";
		filename4 = "C:\\Users\\wangzhi\\Desktop\\matrix\\new2\\b.txt";
		break;
	}
	case 14133: {
		filename1 = "C:\\Users\\wangzhi\\Desktop\\matrix\\new3\\row.txt";
		filename2 = "C:\\Users\\wangzhi\\Desktop\\matrix\\new3\\col.txt";
		filename3 = "C:\\Users\\wangzhi\\Desktop\\matrix\\new3\\num.txt";
		filename4 = "C:\\Users\\wangzhi\\Desktop\\matrix\\new3\\b.txt";
		//const char* filename5 = NULL;
		//filename5 = "C:\\Users\\wangzhi\\Desktop\\AW-COL.txt";
		//std::ifstream readwCOL;
		//readwCOL.open(filename5);

		//if (!readwCOL.good())
		//{
		//	std::cout << "num读取数据出错";
		//}
		//for (int i=0; (!readwCOL.eof()); i++)
		//{
		//	readwCOL >> wCOL[i];
		//}
		break;
	}

	case 117171: {
		filename1 = "C:\\Users\\wangzhi\\Desktop\\matrix\\mesh_0_05\\Colptr.txt";
		filename2 = "C:\\Users\\wangzhi\\Desktop\\matrix\\mesh_0_05\\Rowind.txt";
		filename3 = "C:\\Users\\wangzhi\\Desktop\\matrix\\mesh_0_05\\Crs_val.txt";
		filename4 = "C:\\Users\\wangzhi\\Desktop\\matrix\\mesh_0_05\\pB.txt";
		const char* filename5 = NULL;
		filename5 = "C:\\Users\\wangzhi\\Desktop\\matrix\\mesh_0_05\\WZ_Voltage.txt";
		std::ifstream readwCOL;
		readwCOL.open(filename5);

		if (!readwCOL.good())
		{
			std::cout << "num读取数据出错";
		}
		for (int i=0; (!readwCOL.eof()); i++)
		{
			readwCOL >> WZ_Voltage[i];
		}
		break;
	}
	default:
		break;
	}
	/*存储文件数组的指针数组*/
	row = new int[temp1];
	col = new int[temp2];
	num = new double[temp2];
	b = new double[temp1]; 

	std::ifstream readrow, readcol, readnum, readb;
	readrow.open(filename1);
	readcol.open(filename2);
	readnum.open(filename3);
	readb.open(filename4);

	/*记录读取数字的个数*/
	int allnum = 0;
	int tempnmax = 0;
	int tempcolmax = 0;

	/*读取文件*/
	if (!readrow.good())
	{
		std::cout << "row读取数据出错";
	}
	for (; (!readrow.eof()); tempnmax++)
	{
		readrow >> row[tempnmax];
	}
	if (!readcol.good())
	{
		std::cout << "col读取数据出错";
	}
	for (; (!readcol.eof()); tempcolmax++)
	{
		readcol >> col[tempcolmax];
	}
	if (!readnum.good())
	{
		std::cout << "num读取数据出错";
	}
	for (; (!readnum.eof()); allnum++)
	{
		readnum >> num[allnum];
	}

	if (!readb.good())
	{
		std::cout << "b读取数据出错";
	}
	for (int i = 0; (!readb.eof()); i++)
	{
		readb >> b[i];
	}

	if (tempcolmax != allnum)
		std::cout << "输入的列数和非零元个数不匹配"<<std::endl;
	if (row[tempnmax-1] != allnum)
		std::cout << "输入的非零元个数和行记录不匹配"<<std::endl;
	 
	/*构造压缩矩阵*/
    TernArray temp(row, col, num, nmax, colmax, allnum);

	/*************************/
	//temp.PrintfMatrix();
	//std::cout << std::endl;
	//double* pX = new double[nmax];
	//if (!temp.SolveSPDMLU(pX, b)) {
	//	std::cout << "LU失败" << std::endl;
	//}
	/************************/
    
	/*构造AMG对象*/
	AMG TestAmg(temp, b, CYCLES, COARSEST, SW_BOUND);

	/*构造初始变量*/
	OneDimensionalArray lastx(colmax);
	
	/*开始执行AMG算法*/
	int v_max = 1;
	int times=TestAmg.MainReSult(&lastx,v_max);
	
	///*测试CG*/
	//double* tempX = new double[nmax];
	//memset(tempX, 0, sizeof(double)* nmax);
	//int times = TestAmg.SolveEquationsCG(&temp, 1e6, nmax,tempX,b);
	//lastx.Initial(nmax, tempX);
	//delete tempX;
	//tempX = NULL;
	///////////////////

	///*测试UseEigen*/
	//OneDimensionalArray tempB(nmax);
	//tempB.Initial(nmax, b);
	//TestAmg.UseEigen(0,&lastx,&tempB);
	///////////////////


	const char* filenamepX = "C:\\Users\\wangzhi\\Desktop\\matrix\\mesh_0_05\\pX.txt";
	std::ofstream writepX;
	writepX.open(filenamepX);
	if (!writepX.good())
	{
		std::cout << "写数据出错" << std::endl;
	}
	for (int i = 0; i < (tempnmax - 1); i++)
	{
		writepX << lastx.ReturnNum(i) << " ";
	}

	end = clock();
	end = end - start;
	const char* filename_res = "C:\\Users\\wangzhi\\Desktop\\matrix\\mesh_0_05\\result.txt";
	std::ofstream writefile;
	writefile.open(filename_res);
	if (!writefile.good())
	{
		std::cout << "写数据出错" << std::endl;
	}
	writefile << "A矩阵的阶数: " << max << std::endl;
	writefile << "A矩阵的非零元个数: " << allnum << std::endl;
	writefile << "运行时间: " << end << "ms" << std::endl;
	writefile << "网格层数: " << COARSEST << std::endl;
	OneDimensionalArray tempb(nmax);
    tempb.Initial(nmax, b);
    OneDimensionalArray tempu_in(colmax);
	temp.MultiplyVector(&lastx, &tempu_in);
	tempu_in.reduce(&tempb,true);
	writefile << "绝对误差: " << tempu_in.norm2() << std::endl;
	writefile << "相对误差: " << tempu_in.norm2() / tempb.norm2() << std::endl;
//	lastx.reduce(WZ_Voltage);
	writefile << "与MTSS的误差: " <<lastx.norm2()<< std::endl;
	writefile << "V循环次数: " << v_max << std::endl;
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	K32GetProcessMemoryInfo(handle, &pmc, sizeof(pmc));
	unsigned long long vmPeakMem = pmc.PeakWorkingSetSize / (1024 * 1024);
	writefile << "最大峰值内存: " << vmPeakMem << std::endl;
	writefile << "共轭梯度法迭代次数: " << times << std::endl;



	return;
}




int main()
{


	int CYCLES = 20;
	/*网格层数*/
	int COARSEST = 2;
	/*判断强关联点的系数*/
	double SW_BOUND = 0.25;

	/*选择矩阵--通过每个矩阵对应的阶数来选择需要的矩阵*/
	int chose = 4;
	int num = 0;
	switch (chose)
	{
	case 0:
		num = 28109;
		break;
	case 1:
		num = 4;
		break;
	case 2:
		num = 9;
		break;
	case 3:
		num = 36;
		break;
	case 4:
		num = 14133;
		break;
	case 5:
		num = 117171;
		break;
	default:
		break;
	}
	/*转到生成矩阵的函数*/
	getAmatrix(num, CYCLES, COARSEST, SW_BOUND);	
	

	return 0;
}