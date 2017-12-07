//==========================================
//Authors: yin le 
//emila:   yinle@neolix.cn
//time:    2017/11/10
//address: Beijing
//==========================================

#ifndef MATRIX_H_
#define  MATRIX_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

#ifndef _MSC_VER
	#include <stdint.h>
#else
typedef __int8			 int8_ta;
typedef __int16			 int16_t;
typedef __int32			 int32_t;
typedef __int64			 int64_t;
typedef unsigned __int8  uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
	typedef double FLOAT;
#endif // !_MSC_VER

#define SWAP(a,b) {temp=a;a=b;b=temp;}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
	static FLOAT sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
	static FLOAT maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
	static int32_t iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))


class CMatrix
{
public:
	//��ʼ�����캯���Լ��ڲ�struct��Ա
	CMatrix(){ this->SMatrix.row = 0; this->SMatrix.col = 0; this->SMatrix.VAL = 0; };
	CMatrix(const int32_t m_, const int32_t n_){
		if (m_<0||n_<0)
		{
			std::cout << "input must up zero!" << std::endl;
			exit(0);
		}
		this->SMatrix.row = m_;
		this->SMatrix.col = n_; 
		allocateMemory(m_, n_);
	};
	CMatrix(const int32_t m_, const int32_t n_, const FLOAT* VAL_)
	{ 
		this->SMatrix.row = m_;
		this->SMatrix.col = n_;
		allocateMemory(m_, n_);
		int32_t k = 0;
		for (int32_t i = 0; i < m_; i++)
			for (int32_t j = 0; j < n_; j++)
				this->SMatrix.VAL[i][j] = VAL_[k++];
	}
	CMatrix(const CMatrix &M){ 
		this->SMatrix.row = M.SMatrix.row;
		this->SMatrix.col = M.SMatrix.col; 
		allocateMemory(M.SMatrix.row, M.SMatrix.col);
	for (int32_t i = 0; i < M.SMatrix.row; i++)
		memcpy(this->SMatrix.VAL[i], M.SMatrix.VAL[i], M.SMatrix.col*sizeof(FLOAT));
	};
	struct SMatrix
	{
	public:
		int32_t row;
		int32_t col;
		FLOAT **VAL;
	}SMatrix;
	~CMatrix();
	
	
public:
	//�з������Ͳ�������
	CMatrix getMat(const int32_t row1, const int32_t row2, const int32_t col1, const int32_t col2);
	CMatrix getrow(const int32_t row);
	CMatrix getcol(const int32_t col);
	CMatrix LU();
	CMatrix& t();
	CMatrix& invMat();//�����
	
public:
	//�޷������Ͳ�������
	void setMat(const CMatrix &M, const int32_t x, const int32_t y);
	void clearMat();
	void eyes();
	void ones();
	void printfMat();
	FLOAT sqatMat();
	FLOAT meanMat();
	FLOAT det();//����ʽ
	FLOAT sumMat();
	FLOAT pythag(FLOAT a, FLOAT b);
	//�ж����߶��Ƿ��н��㣬����з��ص��������ã����������ʽx y��x y ...
	bool JudgelineSegmentIntersection(const CMatrix& A,const CMatrix& B, CMatrix& out_);
	void svd(CMatrix &U2, CMatrix &W, CMatrix &V);
public:
	//���ڲ���������
	CMatrix& operator= (const CMatrix& M);
	CMatrix& operator- (const CMatrix& M);
	CMatrix& operator- (const FLOAT M);
	CMatrix& operator+ (const CMatrix& M);
	CMatrix& operator+ (const FLOAT M);
	CMatrix operator* (const CMatrix& M);
	CMatrix& operator* (const FLOAT M);
	CMatrix& operator/ (const CMatrix& M);
	CMatrix& operator/ (const FLOAT M);
public:
	//���־�̬��Ա����������ֱ���๹�����
	static CMatrix eyes(const int32_t row, const int32_t col);
	static CMatrix ones(const int32_t row, const int32_t col);
	static CMatrix RotMatX(const FLOAT angle);
	static CMatrix RotMatY(const FLOAT angle);
	static CMatrix RotMatZ(const FLOAT angle);
	static CMatrix pow(CMatrix& M,int32_t n);
	static CMatrix solve(CMatrix& A,const CMatrix& C);
private:
	void LU(CMatrix& L, CMatrix& U, CMatrix& P);
	CMatrix LU(int32_t& changeNum);
	void allocateMemory(const int32_t m_, const int32_t n_);
	void releaseMemory();
	void swap(FLOAT &a, FLOAT& b);
};


#endif //Matrix.h
