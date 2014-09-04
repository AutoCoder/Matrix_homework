#pragma once
#include "matrix.h"
#include <vector>


//稠密矩阵
template <class U>

class DenseMatrix :public Matrix<U>
{
	
public:	
	virtual ~DenseMatrix(void);
	//拷贝构造函数
	DenseMatrix (const DenseMatrix<U>& m);
	//构造函数
	DenseMatrix (const size_t row,const size_t col, U** val=0);
	//赋值构造函数
	DenseMatrix<U>& operator = (const DenseMatrix<U>& m);

	virtual size_t RowNo () const { return this->_m->_row; }
	virtual size_t ColNo () const { return this->_m->_col; }

	////转换成稀疏矩阵
	//SparseMatrix<U> ToSparseMatrix() const; 

	//取指定位置的值（引用）
	virtual U& operator () (size_t row, size_t col);
	//返回值
	virtual U  operator () (size_t row, size_t col) const;

	//矩阵操作符重载
	//(原矩阵被修改)
	virtual Matrix<U>& operator += (const Matrix<U>& m);
	virtual Matrix<U>& operator -= (const Matrix<U>& m);
	virtual Matrix<U>& operator *= (const Matrix<U>& m);
	virtual Matrix<U>& operator *= (const U& c);
	virtual Matrix<U>& operator /= (const U& c);
	
	////转置，返回一个转置后的新矩阵
	virtual DenseMatrix<U> Transpose() const;
	virtual void Clone();

private:
	//为了实现多个引用操作同一块内存而又不会被不经意释放，必须实现引用计数
	struct base_mat//必须把数据成员和引用计数封装到一起
	{
		U **Data;
		size_t _row, _col;
		int Refcnt;

		base_mat (size_t row, size_t col, U** v)
		{
			_row = row; 
			_col = col; 
			Refcnt = 1;

			Data = new U* [row];
			size_t rowlen = _col * sizeof(U);

			for (size_t i=0; i < _row; i++)
			{
				Data[i] = new U [_col];
				for(size_t j=0;j<_col;j++)
				{
					Data[i][j]=0;
				}
				if (v) memcpy( Data[i], v[i], rowlen);
			}
		}
		~base_mat ()
		{
			for (size_t i=0; i < _row; i++)
			{	
				delete [] Data[i];
			}
			delete [] Data;
		}
	};
	base_mat *_m;
};

template <class U> inline
DenseMatrix<U>::DenseMatrix(size_t row, size_t col, U** val=0)
{
	_m = new base_mat( row, col, val);
}

template <class U> inline
DenseMatrix<U>::~DenseMatrix(void)
{
	if (--_m->Refcnt == 0) 
	{
		delete _m;
	}
}

template <class U> inline
DenseMatrix<U>::DenseMatrix (const DenseMatrix<U>& m)
{
	_m = m._m;
	_m->Refcnt++;
}

template <class U> inline
DenseMatrix<U>& DenseMatrix<U>::operator = (const DenseMatrix<U>& m)
{
	m._m->Refcnt++;
	if (--_m->Refcnt == 0) delete _m;
	_m = m._m;
	return *this;
}

//[][]难以重载，重载(i,j)去矩阵值
template <class U> inline
 U& DenseMatrix<U>::operator () (size_t row, size_t col)
{
	if(row >= _m->_row||col >= _m->_col)
	{
		throw logic_error("矩阵行数或列数越界!");
	}
	if (_m->Refcnt > 1) Clone();
	return _m->Data[row][col];
}

//矩阵常量返回值类型
template <class U> inline
U DenseMatrix<U>::operator () (size_t row, size_t col) const
{
	if(row >= _m->_row||col >= _m->_col)
	{
		throw logic_error("矩阵行数或列数越界!");
	}
	return _m->Data[row][col];
}

template <class U> inline
Matrix<U>& DenseMatrix<U>::operator += (const Matrix<U>& m)
{
	if(_m->_row!=m.RowNo()||_m->_col!=m.ColNo())
	{
		throw logic_error("矩阵行数或列数不匹配!");
	}
	if (_m->Refcnt > 1) Clone();
	for (size_t i=0;i<_m->_row;i++)
	{
		for (size_t j=0;j<_m->_col;j++)
		{
			_m->Data[i][j]+=m(i,j);
		}
	}
	return *this;

}

template <class U> inline
Matrix<U>& DenseMatrix<U>::operator -= (const Matrix<U>& m)
{

	if(_m->_row!=m.RowNo()||_m->_col!=m.ColNo())
	{
		throw logic_error("矩阵行数或列数不匹配!");
	}
	if (_m->Refcnt > 1) Clone();
	for (size_t i=0;i<_m->_row;i++)
	{
		for (size_t j=0;j<_m->_col;j++)
		{
			_m->Data[i][j] -= m(i,j);
		}
	}
	return *this;

}

//乘法：（原矩阵）在左，（乘矩阵）在右
template <class U> inline
Matrix<U>& DenseMatrix<U>::operator *= (const Matrix<U>& m)
{
	if(_m->_col!=m.RowNo())
	{
		throw logic_error("相乘矩阵行数或列数不匹配!");
	}
	
	DenseMatrix<U> temp(_m->_row,m.ColNo());

	for (size_t i=0; i < _m->_row; i++)
	{	
		for (size_t j=0; j < m.ColNo(); j++)
		{
			temp(i,j) = U(0);
			for (size_t k=0; k < _m->_col; k++)
			{
				temp(i,j) += _m->Data[i][k] * m(k,j);
			}
		}
	}
	*this = temp;
	return *this;

}

//乘法：矩阵乘以系数
template <class U> inline
Matrix<U>& DenseMatrix<U>::operator *= (const U& c)
{
	if (_m->Refcnt > 1) Clone();
	for (size_t i=0; i < _m->_row; i++)
	{
		for (size_t j=0; j < _m->_col; j++)
		{
			_m->Data[i][j] *= c;
		}
	}
	return *this;
}

//矩阵除系数
template <class U> inline
Matrix<U>& DenseMatrix<U>::operator /= (const U& c)
{
	if (_m->Refcnt > 1) Clone();
	for (size_t i=0; i < _m->_row; i++)
	{	
		for (size_t j=0; j < _m->_col; j++)	
		{
			_m->Data[i][j] /= c;
		}
	}
	return *this;
}

//转置，不修改原矩阵
template <class U> inline
DenseMatrix<U> DenseMatrix<U>::Transpose() const
{
	DenseMatrix<U> temp(_m->_col,_m->_row);

	for (size_t i=0; i < _m->_row; i++)
		for (size_t j=0; j < _m->_col; j++)
		{
			U x = _m->Data[i][j];
			temp(j,i) = x;
		}
		return temp;
}

template <class U> inline
void DenseMatrix<U>::Clone()
{
	_m->Refcnt--;
	_m = new base_mat( _m->_row, _m->_col, _m->Data);

}


