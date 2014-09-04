// MatrixFamily.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Matrix.h"
#include "DenseMatrix.h"
#include "DiagonalMatrix.h"
#include "SparseMatrix.h"
using namespace std;

//矩阵输出
template <class U> inline 
ostream& operator << (ostream& ostrm, const const Matrix<U>& m)
{
	for (size_t i=0; i < m.RowNo(); i++)
	{
		for (size_t j=0; j < m.ColNo(); j++)
		{
			U x = m(i,j);
			ostrm << x << '\t';
		}
		ostrm << endl;
	} 
	return ostrm;
}

//矩阵输入
template <class U> inline  
istream& operator >> (istream& istrm, Matrix<U>& m)
{
   for (size_t i=0; i < m.RowNo(); i++)
      for (size_t j=0; j < m.ColNo(); j++)
      {
         U x;
         istrm >> x;
         m(i,j) = x;
      }
   return istrm;
}


int _tmain(int argc, _TCHAR* argv[])
{
	
	cout<<"cm开头的为稠密矩阵，dj开头的为对角矩阵，xs开头的为稀疏矩阵"<<endl;
	cout<<"\n第一段为稠密矩阵的测试"<<endl;
    system("pause");
	//------------稠密矩阵测试------------
	float** cmData = new float*[4];
	for(int i=0; i<4; i++)
	{
		cmData[i] = new float[4];
	}
	for(int i=0; i<4; i++)
	{
		for (int j=0; j<4; j++)
		{
			cmData[i][j] = i*j;
		}
	}
	DenseMatrix<float> cm1(4,4, cmData);
	DenseMatrix<float> cm2(4,4, cmData);
	cout<<"cm1=\n"<<cm1<<endl;
	cout<<"cm2=\n"<<cm2<<endl;

	cm1+=cm2;
	cout<<"cm1+=cm2,cm1=\n"<<cm1<<endl;
	cm1*=cm2;
	cout<<"cm1*=cm2,cm1=\n"<<cm1<<endl;
    cm2(0,2)=5;  
	cout<<"cm1(0,2)=5,cm1=\n"<<cm2<<endl;
	cout<<"cm2转置=\n"<<cm2.Transpose()<<endl;
	
	//-------------稠密矩阵测试完毕-------------
	
	cout<<"下一段为对角矩阵的测试"<<endl;
    system("pause");

	//-------------对角矩阵测试-----------------
	float* djData1 = new float[5];
	float* djData2 = new float[5];
	for(int i=0; i<3; i++)
	{
			djData1[i] = i+1;
			djData2[i] = (i+1)*(i+1);
	}
	DiagonalMatrix<float> dj1(3,djData1);
	DiagonalMatrix<float> dj2(3,djData2);
	//dj1(1,2)=5.0;抛出异常
	cout<<"dj1=\n"<<dj1<<endl;
	cout<<"dj2=\n"<<dj2<<endl;
	cout<<"dj1的转置=\n"<<dj1.Transpose()<<endl;
	dj1+=dj2;
	cout<<"dj1+=dj2,dj1=\n"<<dj1<<endl;
	//-----------对角矩阵测试完毕-----------

	cout<<"下一段为稀疏矩阵的测试"<<endl;
    system("pause");
	//--------------稀疏矩阵测试----------------------
	       //-------乘法测试---------
	SparseData<float>* sdata1 = new SparseData<float>[7];
	for (size_t i=0;i<7;i++)
	{
		sdata1[i].c=6-i;
		sdata1[i].r=i;
		sdata1[i].Val=i;
	}
	SparseMatrix<float> xs1(7,7,7,sdata1);
	cout<<"xs1=\n"<<xs1<<endl;
	cout<<"xs1的转置=\n"<<xs1.Transpose()<<endl;
	
	SparseData<float>* sdata2 = new SparseData<float>[8];
	SparseData<float>* sdata3 = new SparseData<float>[5];
	sdata2[0].r=0;sdata2[0].c=0;sdata2[0].Val=1;
	sdata2[1].r=0;sdata2[1].c=2;sdata2[1].Val=2;
	sdata2[2].r=1;sdata2[2].c=1;sdata2[2].Val=1;
	sdata2[3].r=1;sdata2[3].c=2;sdata2[3].Val=-1;
	sdata2[4].r=2;sdata2[4].c=0;sdata2[4].Val=-1;
	sdata2[5].r=2;sdata2[5].c=1;sdata2[5].Val=2;
	sdata2[6].r=3;sdata2[6].c=0;sdata2[6].Val=2;
	sdata2[7].r=3;sdata2[7].c=1;sdata2[7].Val=1;

	sdata3[0].r=0;sdata3[0].c=0;sdata3[0].Val=1;
	sdata3[1].r=0;sdata3[1].c=1;sdata3[1].Val=2;
	sdata3[2].r=1;sdata3[2].c=0;sdata3[2].Val=2;
	sdata3[3].r=1;sdata3[3].c=1;sdata3[3].Val=1;
	sdata3[4].r=2;sdata3[4].c=1;sdata3[4].Val=3;
    SparseMatrix<float> xs2(4,3,8,sdata2);
    SparseMatrix<float> xs3(3,2,5,sdata3);
	cout<<"xs2=\n"<<xs2<<endl;
	cout<<"xs3=\n"<<xs3<<endl;
	xs2*=xs3;
	cout<<"xs2*=xs3,xs2=\n"<<xs2<<endl;

	//-------------加减法运算------------
	SparseData<float>* sdata4 = new SparseData<float>[2];
	SparseData<float>* sdata5 = new SparseData<float>[5];
	sdata4[0].r=0;sdata4[0].c=1;sdata4[0].Val=1;
	sdata4[1].r=2;sdata4[1].c=2;sdata4[1].Val=1;

	sdata5[0].r=0;sdata5[0].c=0;sdata5[0].Val=1;
	sdata5[1].r=0;sdata5[1].c=1;sdata5[1].Val=1;
	sdata5[2].r=2;sdata5[2].c=2;sdata5[2].Val=3;
	
    SparseMatrix<float> xs4(3,3,2,sdata4);
    SparseMatrix<float> xs5(3,3,3,sdata5);
	cout<<"xs4=\n"<<xs4<<endl;
	cout<<"xs5=\n"<<xs5<<endl;
	xs4+=xs5;
	cout<<"xs4+=xs5,xs5=\n"<<xs4<<endl;
	xs4-=xs5;
	cout<<"xs4-=xs5,xs4=\n"<<xs4<<endl;
    //-----------稀疏矩阵测试完毕
	//------------------------------------
	cout<<"输入稠密矩阵，测试与其他类型矩阵运算!"<<endl;
	system("pause");

	//由于对角矩阵和稀疏矩阵数据结构的原因，无法随意取矩阵元素的引用
	//否则会抛出逻辑错误，所以不能通过cin输入对角矩阵和稀疏矩阵
	DenseMatrix<float> cm3(3,3);
	cout<<"输出3X3稠密矩阵cm3:\n";
	cin>>cm3;
	cout<<"输入完毕!\n  克隆cm3到cm3_1  ;  xs4到xs4_1"<<endl;
	DenseMatrix<float> cm3_1=cm3;
	SparseMatrix<float> xs4_1=xs4;
	cout<<"cm3=\n"<<cm3<<endl;
	cout<<"xs4=\n"<<xs4<<endl;
	cout<<"dj2=\n"<<dj2<<endl;
	cm3*=xs4;
	xs4*=dj2;
	xs4_1*=cm3_1;
	cout<<"cm3*=xs4,cm3=\n"<<cm3<<endl;
	cout<<"xs4*=dj2,xs4=\n"<<xs4<<endl;
	cout<<"xs4_1*=cm3_1,xs4_1=\n"<<xs4_1<<endl;
    //由于数据结构限制，对角矩阵无法直接左乘稀疏或稠密矩阵（会出现逻辑错误）
	//需要将它转换到稠密或稀疏再乘
	DenseMatrix<float> dj2cm=dj2.ToDenseMatrix();
	dj2cm*=cm3_1;
	cout<<"将dj2转换成稠密矩阵,dj2cm*=cm3_1,dj2cm="<<endl;
	cout<<dj2cm<<endl;
    cout<<"测试完毕!"<<endl;
	system("pause");
	
}

