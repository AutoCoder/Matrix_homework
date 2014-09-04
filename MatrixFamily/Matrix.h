#pragma once
#include "stdafx.h"
#include <iostream>

//基类
template <class U>
class Matrix
{
public:

	virtual ~Matrix(void);
public:
	//////构造函数
	Matrix ();
	

	//返回行数和列数
	virtual size_t RowNo () const=0;
	virtual size_t ColNo () const=0; 
	virtual void Clone()=0;

	//取值（引用）
	virtual U& operator () (size_t row, size_t col)=0;
	//返回值
	virtual U  operator () (size_t row, size_t col) const=0;

	//矩阵操作符重载
	//(原矩阵被修改)
	virtual Matrix<U>& operator += (const Matrix<U>& m)=0;
	virtual Matrix<U>& operator -= (const Matrix<U>& m)=0;
	virtual Matrix<U>& operator *= (const Matrix<U>& m)=0;
	virtual Matrix<U>& operator *= (const U& c)=0;
	virtual Matrix<U>& operator /= (const U& c)=0;
};

template <class U>
Matrix<U>::Matrix ()
{

}

template <class U>
Matrix<U>::~Matrix ()
{

}