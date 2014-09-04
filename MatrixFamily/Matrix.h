#pragma once
#include "stdafx.h"
#include <iostream>

//����
template <class U>
class Matrix
{
public:

	virtual ~Matrix(void);
public:
	//////���캯��
	Matrix ();
	

	//��������������
	virtual size_t RowNo () const=0;
	virtual size_t ColNo () const=0; 
	virtual void Clone()=0;

	//ȡֵ�����ã�
	virtual U& operator () (size_t row, size_t col)=0;
	//����ֵ
	virtual U  operator () (size_t row, size_t col) const=0;

	//�������������
	//(ԭ�����޸�)
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