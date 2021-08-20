#include"ST.h"

STMatrix::STMatrix(int temprowmax, int tempcolmax, int tempallnum) :rowmax(temprowmax),colmax(tempcolmax)
{
	 this->STRowArray = new int[rowmax];
	 this->ST = new int* [rowmax];
	 
}

STMatrix::STMatrix(int** tempST, int* len, int nmax,int tempall) :rowmax(nmax), colmax(nmax),allnum(tempall)
{
	//�����õ�����һ�����ڴ��Ǹ�����
	ST = tempST;
	STRowArray = len;
	
}

STMatrix::~STMatrix()
{
	delete[] STRowArray;
	STRowArray = NULL;
	for (int i = 0; i < rowmax; i++) {
		delete[] ST[i];
		ST[i] = NULL;
	}
	delete[] ST;
	ST = NULL;

}

void STMatrix::OrderInitial(int* i,int length,int row)
{
	ST[row] = i;
	STRowArray[row]=length;
	this->allnum = this->allnum + length;
}


void STMatrix::GetAMGNode(TernArray** Weight, int level)
{
	/*���з���Ԫ*/
	std::vector<std::vector<node>>array(this->colmax);
	   
	for (int i = 0; i < rowmax; i++) {
		for (int tempj = 0; tempj < STRowArray[i]; tempj++) {
			node temp;
			temp.col = i;
			temp.num = tempj;
			array[ST[i][tempj]].push_back(temp);
		}
	}

	int* STRow = new int[rowmax];
	memcpy(STRow, STRowArray, sizeof(int) * rowmax);
	int* STColArray=new int[colmax];
	for (int i = 0; i < colmax; i++) {
		STColArray[i] = array[i].size();
	}

	/*�����ֵ����*/
	int nG = 0;
	Weight[level] = new TernArray(this->colmax, this->rowmax, this->allnum);
	std::set<int>* allnode = new std::set<int>[this->rowmax];
	/*�����з���Ԫ���ĵ�i��Ȼ��i�����з���Ԫ����ֵj��ʹWeight[level](nG,j)=1��Ȼ�󽫷���Ԫ�Ķ�Ӧ��һ�е�ֵ���㣬�����Щ����������е�k��ȫΪ0�ͽ��õ��Weight[level](nG,k)=1*/
	std::vector<int> v_set(colmax,1);
	int RowMax = this->RowMax(STColArray, colmax);
	while (STColArray[RowMax]!=-1) {
	/////////////////////////////
		std::vector<int>temptime;
		for (int j = 0; j < STRowArray[RowMax]; j++) {
				if (v_set[ST[RowMax][j]] != 0) {
					temptime.push_back(ST[RowMax][j]);
					v_set[ST[RowMax][j]] = 0;
					STColArray[ST[RowMax][j]] = -1;
				}
		}
		if (v_set[RowMax] != 0) {
			temptime.push_back(RowMax);
			STColArray[RowMax] = -1;
			v_set[RowMax] = 0;
		}


	    /// �Ƿ���ȷδ֪
		for (int it = 0; it < temptime.size(); it++) {
			for (std::vector<node>::iterator it01 = array[temptime[it]].begin(); it01 < array[temptime[it]].end(); it01++) {
				if (v_set[(*it01).col] != 0) {
					STColArray[(*it01).col] = STColArray[(*it01).col] - 1;
					if (STColArray[(*it01).col] == 0) {
						// indI.push_back((*it01).col);
						temptime.push_back((*it01).col);
						STColArray[(*it01).col] = -1;
						v_set[(*it01).col] = 0;
					}
				}
			}
		}
	    /// </summary>
	    /// <param name="Weight"></param>
	    /// <param name="level"></param>
	    sort(temptime.begin(), temptime.end());
	    for (int it = 0; it < temptime.size(); it++) {
				Weight[level]->OrderAssignment(nG, temptime[it], 1);
				allnode[temptime[it]].insert(nG);
	    }
		nG++;
		RowMax = this->RowMax(STColArray, rowmax);
	}
	/*������ֵ*/
	Weight[level]->ChangeRow(nG);

	TernArray* transW=new TernArray(this->rowmax, nG, this->allnum);
	for (int i = 0; i < this->rowmax; i++) {
		for (auto it : allnode[i]) {
			transW->OrderAssignment(i, it, 1);
		}
	}
	//ȱ��ת��transW��ֵ
}

int STMatrix::RowMax(int* STRow,int length)
{
	int b = 0;
	for (int i = 0; i < length; i++) {
		b = (STRow[b] < STRow[i]) ? i : b;
	}
	return b;
}