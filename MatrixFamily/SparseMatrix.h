#pragma once
#include "matrix.h"
#include <math.h>





//稀疏矩阵--大部分元素为0，一般用三元组存贮数据
template <class U>
struct SparseData
{
	size_t r,c;
	U Val;
	SparseData(){ }
	template <class U>
	SparseData(size_t x,size_t y,U v)
	{
		r=x;
		c=y;
		Val=v;
	}	
};

template <class U>
class SparseMatrix :public Matrix<U>
{

public:
	virtual ~SparseMatrix(void);
	//拷贝构造函数
	SparseMatrix<U> (const SparseMatrix& m);
	//构造函数
	SparseMatrix<U> (const size_t row,const size_t col,const size_t NozeroCount,SparseData<U>* val=0);
	//赋值构造函数
	SparseMatrix<U>& operator = (const SparseMatrix& m);

	virtual size_t RowNo () const { return _sm->_row; }
	virtual size_t ColNo () const { return _sm->_col; }
	size_t GetNozeroCount() const { return _sm->_NozeroCount; }

	//取指定位置的值（引用）
	virtual U& operator () (size_t row, size_t col);
	//返回值
	virtual U  operator () (size_t row, size_t col) const;


	/*矩阵操作符重载
	(原矩阵被修改)*/
	virtual Matrix<U>& operator += (const Matrix& m);
	virtual Matrix<U>& operator -= (const Matrix& m);
	virtual Matrix<U>& operator *= (const Matrix& m);
	virtual Matrix<U>& operator *= (const U& c);
	virtual Matrix<U>& operator /= (const U& c);

	////转置，返回一个转置后的新矩阵
	virtual SparseMatrix<U> Transpose() const;
	virtual void Clone();
	
private:
	template <class U>
	struct Spar_mat//必须把数据成员和引用计数封装到一起
	{
		SparseData<U>* Spar_Data;
		size_t _row, _col,_NozeroCount;
		int Refcnt;

		Spar_mat (size_t row, size_t col,size_t NozeroCount, SparseData<U>* val)
		{
			_row = row; 
			_col = col; 
			Refcnt = 1;
			_NozeroCount = NozeroCount;
			size_t rowlen = NozeroCount * sizeof(SparseData<U>);
			Spar_Data = new SparseData<U>[rowlen];
			if(val)
			{
				memcpy(Spar_Data, val, rowlen);
			}
		}
		~Spar_mat ()
		{

			delete [] Spar_Data;
		}
	};
	Spar_mat<U> *_sm;
};

template <class U> inline
SparseMatrix<U>::SparseMatrix(const size_t row,const size_t col,const size_t NozeroCount,SparseData<U>* val=0)
{
	_sm = new Spar_mat<U>( row, col,NozeroCount ,val);
}

template <class U> inline
SparseMatrix<U>::~SparseMatrix(void)
{
	if (--_sm->Refcnt == 0) 
	{
		delete _sm;
	}
}

template <class U> inline
SparseMatrix<U>::SparseMatrix (const SparseMatrix<U>& m)
{
	_sm = m._sm;
	_sm->Refcnt++;
}

template <class U> inline
SparseMatrix<U>& SparseMatrix<U>::operator = (const SparseMatrix<U>& m)
{
	m._sm->Refcnt++;
	if (--_sm->Refcnt == 0) delete _sm;
	_sm = m._sm;
	return *this;
}

template <class U> inline
U& SparseMatrix<U>::operator () (size_t row, size_t col)
{
	if (row>=_sm->_row||col>=_sm->_col)
	{
		throw logic_error("矩阵行数或列数越界!");
	}
	if (_sm->Refcnt > 1) Clone();
	for (size_t i=0;i<_sm->_NozeroCount;i++)
	{
		if(row==_sm->Spar_Data[i].r&&col==_sm->Spar_Data[i].c)
		{
			return _sm->Spar_Data[i].Val;
		}
	}
	throw logic_error("稀疏矩阵中0元素没有引用!");
}//存疑

template <class U> inline
U SparseMatrix<U>::operator () (size_t row, size_t col) const
{
	if (row>=_sm->_row||col>=_sm->_col)
	{
		throw logic_error("矩阵行数或列数越界!");
	}
	for (size_t i=0;i<_sm->_NozeroCount;i++)
	{
		if(row==_sm->Spar_Data[i].r&&col==_sm->Spar_Data[i].c)
		{
			return _sm->Spar_Data[i].Val;
		}
	}
	return U(0);
}

template <class U> inline
Matrix<U>& SparseMatrix<U>::operator += (const Matrix<U>& m)
{
	if(_sm->_row!=m.RowNo()||_sm->_col!=m.ColNo())
	{
		throw logic_error("矩阵行数或列数不匹配!");
	}
	vector<SparseData<U>> tvector;
	int flag=0;
	for (size_t i=0; i < _sm->_NozeroCount ; i++)
	{
		tvector.push_back(_sm->Spar_Data[i]);
	}
	for (size_t i=0;i<m.RowNo();i++)
	{
		for (size_t j=0;j<m.ColNo();j++)
		{
			for (vector<SparseData<U>>::iterator it=tvector.begin();it!=tvector.end();it++)
			{
				if (it->r==i&&it->c==j)
				{
					it->Val+=m(i,j);
					flag=1;
					break;
				}
			}
			if (0==flag)
			{
				SparseData<U> tem;
				tem.r=i;tem.c=j;tem.Val=m(i,j);
				tvector.push_back(tem);
			}
		}
	}
	size_t sdata_len = tvector.size();
	SparseData<U>* sdata=new SparseData<U>[sdata_len];
	for (size_t i=0;i<sdata_len;i++)
	{
		sdata[i] = tvector.at(i);
	}
	SparseMatrix<U> temp(_sm->_col,_sm->_row,sdata_len,sdata);
	*this =temp;
	return *this;

}

template <class U> inline
Matrix<U>& SparseMatrix<U>::operator -= (const Matrix<U>& m)
{
	if(_sm->_row!=m.RowNo()||_sm->_col!=m.ColNo())
	{
		throw logic_error("矩阵行数或列数不匹配!");
	}
	vector<SparseData<U>> tvector;
	int flag=0;
	for (size_t i=0; i < _sm->_NozeroCount ; i++)
	{
		tvector.push_back(_sm->Spar_Data[i]);
	}
	for (size_t i=0;i<m.RowNo();i++)
	{
		for (size_t j=0;j<m.ColNo();j++)
		{
			for (vector<SparseData<U>>::iterator it=tvector.begin();it!=tvector.end();it++)
			{
				if (it->r==i&&it->c==j)
				{
					it->Val-=m(i,j);
					flag=1;
					break;
				}
			}
			if (0==flag)
			{
				SparseData<U> tem;
				tem.r=i;tem.c=j;tem.Val=(-1)*m(i,j);
				tvector.push_back(tem);
			}
		}
	}
	size_t sdata_len = tvector.size();
	SparseData<U>* sdata=new SparseData<U>[sdata_len];
	for (size_t i=0;i<sdata_len;i++)
	{
		sdata[i] = tvector.at(i);
	}
	SparseMatrix<U> temp(_sm->_col,_sm->_row,sdata_len,sdata);
	*this =temp;
	return *this;
}


template <class U> inline
Matrix<U>& SparseMatrix<U>::operator *= (const Matrix<U>& m)
{
	vector<SparseData<U>> tvector;
	
	for (size_t i=0; i < _sm->_NozeroCount; i++)
	{
		int flag=0;
		SparseData<U> tem;//结果矩阵的元素
		size_t r = _sm->Spar_Data[i].r;
		size_t c = _sm->Spar_Data[i].c;
		U v = _sm->Spar_Data[i].Val;
		for (size_t j=0;j<m.ColNo();j++)
		{
			tem.r=r;
			tem.c=j;
			tem.Val=v*m(c,j);
			for (vector<SparseData<U>>::iterator it=tvector.begin();it!=tvector.end();it++)
			{
				if (it->r==tem.r&&it->c==tem.c)
				{
					it->Val+=tem.Val;
					flag=1;
					break;
				}
			}
			if(0==flag){
				tvector.push_back(tem);
			}
		}
	 }
	size_t sdata_len = tvector.size();
	SparseData<U>* sdata=new SparseData<U>[sdata_len];
	for (size_t i=0;i<sdata_len;i++)
	{
		sdata[i] = tvector.at(i);
	}
	SparseMatrix<U> temp(_sm->_row,m.ColNo(),sdata_len,sdata);
	*this =temp;
	return *this;

}

template <class U> inline
Matrix<U>& SparseMatrix<U>::operator *= (const U& c)
{
	if (_sm->Refcnt > 1) Clone();
	for (size_t i=0; i < _sm->_NozeroCount; i++)
	{
		_sm->Spar_Data[i].Val *= c;
	}
	return *this;
}

template <class U> inline
Matrix<U>& SparseMatrix<U>::operator /= (const U& c)
{
	if (_sm->Refcnt > 1) Clone();
	for (size_t i=0; i < _sm->_NozeroCount; i++)
	{
		_sm->Spar_Data[i].Val /= c;
	}
	return *this;
}

template <class U> inline
SparseMatrix<U> SparseMatrix<U>::Transpose() const
{
	//SparseData<U>* Spar_Data;
	SparseMatrix<U> temp(_sm->_col,_sm->_row,_sm->_NozeroCount,_sm->Spar_Data);
	if (_sm->_NozeroCount>0)
	{
		for(size_t i=0;i<_sm->_NozeroCount;i++)
		{
			size_t r=temp._sm->Spar_Data[i].r;
			size_t c=temp._sm->Spar_Data[i].c;
			temp._sm->Spar_Data[i].r=c;
			temp._sm->Spar_Data[i].c=r;
			temp._sm->Spar_Data[i].Val=_sm->Spar_Data[i].Val;
		}
	}
	return temp;

}


template <class U> inline
void SparseMatrix<U>::Clone()
{
	_sm->Refcnt--;
	_sm = new Spar_mat<U>( _sm->_row, _sm->_col, _sm->_NozeroCount,_sm->Spar_Data);
}