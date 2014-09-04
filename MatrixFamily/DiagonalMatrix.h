#pragma once
#include "matrix.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include <math.h>
#include <iostream>

//对角矩阵--为行数等于列数，仅对角线值不全为0的矩阵
template <class U>
class DiagonalMatrix :public Matrix<U>
{
public:
	virtual ~DiagonalMatrix(void);
	//拷贝构造函数
	DiagonalMatrix (const DiagonalMatrix<U>& m);
	//构造函数
	DiagonalMatrix (const size_t Diag_size, U* val=0);
	//赋值构造函数
	DiagonalMatrix<U>& operator = (const DiagonalMatrix<U>& m);

	//矩阵操作符重载
	//(原矩阵被修改)
	virtual Matrix<U>& operator += (const Matrix<U>& m);//与其他类型矩阵运算就不再是对角矩阵了，全部返回稠密矩阵??
	virtual Matrix<U>& operator -= (const Matrix<U>& m);
	virtual Matrix<U>& operator *= (const Matrix<U>& m);//不管与哪种矩阵进行乘或除运算都还是对角矩阵
	virtual Matrix<U>& operator *= (const U& c);
	virtual Matrix<U>& operator /= (const U& c);
    //仅实现自身类型加减??
	virtual DiagonalMatrix<U>& operator += (const DiagonalMatrix<U>& m);
	virtual DiagonalMatrix<U>& operator -= (const DiagonalMatrix<U>& m);

	////(原矩阵不被修改)
	//virtual DiagonalMatrix<U> operator + (const DiagonalMatrix<U>& right);
	//virtual DiagonalMatrix<U> operator - (const DiagonalMatrix<U>& right);
	//virtual DiagonalMatrix<U> operator * (const DiagonalMatrix<U>& right);
	//virtual DiagonalMatrix<U> operator / (const U& c);

	//---实现抽象类的接口---

	virtual size_t RowNo () const { return _dm->_Diag_size; }
	virtual size_t ColNo () const { return _dm->_Diag_size; }
	//取指定位置的值（引用）
	virtual U& operator () (size_t row, size_t col);
	//返回值
	virtual U  operator () (size_t row, size_t col) const;
	//转置，返回一个转置后的新矩阵
	virtual DiagonalMatrix<U> Transpose() const;
	virtual void Clone();
	virtual DenseMatrix<U> ToDenseMatrix();
	virtual SparseMatrix<U> ToSparseMatrix();
	//--------------------------
private:
	struct Diag_mat//必须把数据成员和引用计数封装到一起
	{
		U* Diag_data;
		size_t _Diag_size;
		int Refcnt;

		Diag_mat (size_t Diag_size, U* v)
		{
			_Diag_size = Diag_size; 
			Refcnt = 1;

			Diag_data = new U [Diag_size];
			size_t rowlen = _Diag_size * sizeof(U);
			for(size_t i=0;i<_Diag_size;i++)
			{
				Diag_data[i]=0;
			}
			if (v) memcpy(Diag_data, v, rowlen);
		}
		~Diag_mat ()
		{
			delete [] Diag_data;
		}
	};
	Diag_mat *_dm;
};

template <class U> inline
DiagonalMatrix<U>::DiagonalMatrix(const size_t Diag_size, U* val=0)
{
	_dm = new Diag_mat(Diag_size, val);

}

template <class U> inline
DiagonalMatrix<U>::~DiagonalMatrix(void)
{	
	if (--_dm->Refcnt == 0) 
	{
		delete _dm;
	}
}

template <class U> inline
DiagonalMatrix<U>::DiagonalMatrix (const DiagonalMatrix<U>& m)
{
	_dm = m._dm;
	_dm->Refcnt++;
}

template <class U> inline
DiagonalMatrix<U>& DiagonalMatrix<U>:: operator = (const DiagonalMatrix<U>& m)
{
	m._dm->Refcnt++;
	if (--_dm->Refcnt == 0) delete _dm;
	_dm = m._dm;
	return *this;
}

template <class U> inline
U& DiagonalMatrix<U>::operator () (size_t row, size_t col)
{
	if(row >= _dm->_Diag_size||col >= _dm->_Diag_size)
	{
		throw logic_error("矩阵行数或列数越界!");
	}else
	{
		if (row==col)
		{
			if (_dm->Refcnt > 1) Clone();
			return _dm->Diag_data[row];
		}else
		{
			throw logic_error("0值元素不能更改!");
		}
	}

}

template <class U> inline
U  DiagonalMatrix<U>::operator () (size_t row, size_t col) const
{
	if(row >= _dm->_Diag_size||col >= _dm->_Diag_size)
	{
		throw logic_error("矩阵行数或列数越界!");
	}else
	{
		if (row==col)
		{
			return _dm->Diag_data[row];
		}else
		{
			return U(0);
		}
	}
}

template <class U> inline
Matrix<U>& DiagonalMatrix<U>::operator += (const Matrix<U>& m)
{
	
	if(_dm->_Diag_size!=m.RowNo()||_dm->_Diag_size!=m.ColNo())
	{
		throw logic_error("矩阵行数或列数不匹配!");
	}
	//如果不是对角矩阵，则抛出异常
	for(size_t i=0;i<m.RowNo();i++)
	{
		for(size_t j=0;j<m.ColNo();j++)
		{
			if(i!=j&&fabs((float)m(i,j))>0.00000001)
			{
				throw logic_error("对角矩阵不能与其他类型矩阵相加!");
			}
		}
	}
	if (_dm->Refcnt > 1) Clone();
	for (size_t i=0;i<_dm->_Diag_size;i++)
	{
		_dm->Diag_data[i]+=m(i,i);
	}	
	return *this;
}

template <class U> inline
Matrix<U>& DiagonalMatrix<U>::operator -= (const Matrix<U>& m)
{
	if(_dm->_Diag_size!=m.RowNo()||_dm->_Diag_size!=m.ColNo())
	{
		throw logic_error("矩阵行数或列数不匹配!");
	}
	//如果不是对角矩阵，则抛出异常
	for(size_t i=0;i<m.RowNo();i++)
	{
		for(size_t j=0;j<m.ColNo();j++)
		{
			if(i!=j&&fabs((float)m(i,j))>0.00000001)
			{
				throw logic_error("对角矩阵不能与其他类型矩阵相减!");
			}
		}
	}
	if (_dm->Refcnt > 1) Clone();
	for (size_t i=0;i<_dm->_Diag_size;i++)
	{
		_dm->Diag_data[i]-=m(i,i);
	}	
	return *this;

}

template <class U> inline
Matrix<U>& DiagonalMatrix<U>::operator *= (const Matrix<U>& m)
{
	if(_dm->_Diag_size!=m.RowNo())
	{
		throw logic_error("矩阵行数或列数不匹配!");
	}
	//如果不是对角矩阵，则抛出异常
	for(size_t i=0;i<m.RowNo();i++)
	{
		for(size_t j=0;j<m.ColNo();j++)
		{
			if(i!=j&&fabs((float)m(i,j))>0.00000001)
			{
				throw logic_error("对角矩阵不能与其他类型矩阵相乘!");
			}
		}
	}
	if (_dm->Refcnt > 1) Clone();
	for (size_t i=0;i<_dm->_Diag_size;i++)
	{
		_dm->Diag_data[i]*=m(i,i);
	}	
	return *this;

}

template <class U> inline
Matrix<U>& DiagonalMatrix<U>::operator *= (const U& c)
{
	if (_dm->Refcnt > 1) Clone();
	for (size_t i=0;i<_dm->_Diag_size;i++)
	{
		_dm->Diag_data[i]*=c;
	}	
	return *this;
}

template <class U> inline
Matrix<U>& DiagonalMatrix<U>::operator /= (const U& c)
{
	if (_dm->Refcnt > 1) Clone();
	for (size_t i=0;i<_dm->_Diag_size;i++)
	{
		_dm->Diag_data[i]/=c;
	}	
	return *this;
}

//转置，即原矩阵
template <class U> inline
DiagonalMatrix<U> DiagonalMatrix<U>::Transpose() const
{
	return *this;
}


template <class U> inline
DiagonalMatrix<U>& DiagonalMatrix<U>::operator += (const DiagonalMatrix<U>& m)
{
	if(_dm->_Diag_size!=m.RowNo())
	{
		throw logic_error("矩阵行数或列数不匹配!");
	}
	if (_dm->Refcnt > 1) Clone();
	for (size_t i=0;i<_dm->_Diag_size;i++)
	{
		_dm->Diag_data[i]+=m(i,i);
	}	

	return *this;
}

template <class U> inline
DiagonalMatrix<U>& DiagonalMatrix<U>::operator -= (const DiagonalMatrix<U>& m)
{
	if(_dm->_Diag_size!=m.RowNo())
	{
		throw logic_error("矩阵行数或列数不匹配!");
	}
	if (_dm->Refcnt > 1) Clone();
	for (size_t i=0;i<_dm->_Diag_size;i++)
	{
		_dm->Diag_data[i]-=m(i,i);
	}	

	return *this;

}

template <class U> inline
DenseMatrix<U> DiagonalMatrix<U>::ToDenseMatrix()
{
	DenseMatrix<U> temp(_dm->_Diag_size,_dm->_Diag_size);
	for (size_t i=0;i<_dm->_Diag_size;i++)
	{
		temp(i,i)+=_dm->Diag_data[i];
	}
	return temp;
	
}

template <class U> inline
SparseMatrix<U> DiagonalMatrix<U>::ToSparseMatrix()
{
	SparseMatrix<U> temp(_dm->_Diag_size,_dm->_Diag_size,_dm->_Diag_size);
	for (size_t i=0;i<_dm->_Diag_size;i++)
	{
		temp(i,i)+=_dm->Diag_data[i];
	}
	return temp;

}


template <class U> inline
void DiagonalMatrix<U>::Clone()
{
	_dm->Refcnt--;
	_dm = new Diag_mat(_dm->_Diag_size,_dm->Diag_data);
}