#pragma once
#include "matrix.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include <math.h>
#include <iostream>

//�ԽǾ���--Ϊ�����������������Խ���ֵ��ȫΪ0�ľ���
template <class U>
class DiagonalMatrix :public Matrix<U>
{
public:
	virtual ~DiagonalMatrix(void);
	//�������캯��
	DiagonalMatrix (const DiagonalMatrix<U>& m);
	//���캯��
	DiagonalMatrix (const size_t Diag_size, U* val=0);
	//��ֵ���캯��
	DiagonalMatrix<U>& operator = (const DiagonalMatrix<U>& m);

	//�������������
	//(ԭ�����޸�)
	virtual Matrix<U>& operator += (const Matrix<U>& m);//���������;�������Ͳ����ǶԽǾ����ˣ�ȫ�����س��ܾ���??
	virtual Matrix<U>& operator -= (const Matrix<U>& m);
	virtual Matrix<U>& operator *= (const Matrix<U>& m);//���������־�����г˻�����㶼���ǶԽǾ���
	virtual Matrix<U>& operator *= (const U& c);
	virtual Matrix<U>& operator /= (const U& c);
    //��ʵ���������ͼӼ�??
	virtual DiagonalMatrix<U>& operator += (const DiagonalMatrix<U>& m);
	virtual DiagonalMatrix<U>& operator -= (const DiagonalMatrix<U>& m);

	////(ԭ���󲻱��޸�)
	//virtual DiagonalMatrix<U> operator + (const DiagonalMatrix<U>& right);
	//virtual DiagonalMatrix<U> operator - (const DiagonalMatrix<U>& right);
	//virtual DiagonalMatrix<U> operator * (const DiagonalMatrix<U>& right);
	//virtual DiagonalMatrix<U> operator / (const U& c);

	//---ʵ�ֳ�����Ľӿ�---

	virtual size_t RowNo () const { return _dm->_Diag_size; }
	virtual size_t ColNo () const { return _dm->_Diag_size; }
	//ȡָ��λ�õ�ֵ�����ã�
	virtual U& operator () (size_t row, size_t col);
	//����ֵ
	virtual U  operator () (size_t row, size_t col) const;
	//ת�ã�����һ��ת�ú���¾���
	virtual DiagonalMatrix<U> Transpose() const;
	virtual void Clone();
	virtual DenseMatrix<U> ToDenseMatrix();
	virtual SparseMatrix<U> ToSparseMatrix();
	//--------------------------
private:
	struct Diag_mat//��������ݳ�Ա�����ü�����װ��һ��
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
		throw logic_error("��������������Խ��!");
	}else
	{
		if (row==col)
		{
			if (_dm->Refcnt > 1) Clone();
			return _dm->Diag_data[row];
		}else
		{
			throw logic_error("0ֵԪ�ز��ܸ���!");
		}
	}

}

template <class U> inline
U  DiagonalMatrix<U>::operator () (size_t row, size_t col) const
{
	if(row >= _dm->_Diag_size||col >= _dm->_Diag_size)
	{
		throw logic_error("��������������Խ��!");
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
		throw logic_error("����������������ƥ��!");
	}
	//������ǶԽǾ������׳��쳣
	for(size_t i=0;i<m.RowNo();i++)
	{
		for(size_t j=0;j<m.ColNo();j++)
		{
			if(i!=j&&fabs((float)m(i,j))>0.00000001)
			{
				throw logic_error("�ԽǾ��������������;������!");
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
		throw logic_error("����������������ƥ��!");
	}
	//������ǶԽǾ������׳��쳣
	for(size_t i=0;i<m.RowNo();i++)
	{
		for(size_t j=0;j<m.ColNo();j++)
		{
			if(i!=j&&fabs((float)m(i,j))>0.00000001)
			{
				throw logic_error("�ԽǾ��������������;������!");
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
		throw logic_error("����������������ƥ��!");
	}
	//������ǶԽǾ������׳��쳣
	for(size_t i=0;i<m.RowNo();i++)
	{
		for(size_t j=0;j<m.ColNo();j++)
		{
			if(i!=j&&fabs((float)m(i,j))>0.00000001)
			{
				throw logic_error("�ԽǾ��������������;������!");
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

//ת�ã���ԭ����
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
		throw logic_error("����������������ƥ��!");
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
		throw logic_error("����������������ƥ��!");
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