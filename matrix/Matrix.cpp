#include "Matrix.h"
#include <assert.h>
#include <memory>
#include <algorithm>
//#define NDEBUG
CMatrix::~CMatrix()
{
	releaseMemory();
}

void CMatrix::allocateMemory(const int32_t m_, const int32_t n_)
{
	if (abs(m_)==0 || abs(n_)==0)
	{
		this->SMatrix.VAL = 0;
		return;
	}
	this->SMatrix.VAL = (FLOAT**)malloc(m_*sizeof(FLOAT*));
	*this->SMatrix.VAL = (FLOAT *)calloc(m_*n_, sizeof(FLOAT));
   for (int32_t i = 1; i < m_; i++)
   {
	   this->SMatrix.VAL[i] = this->SMatrix.VAL[i - 1] + n_;
   }
   (*this).clearMat();
}

void CMatrix::releaseMemory(){

	if (this->SMatrix.VAL != 0)
	{
		free(this->SMatrix.VAL[0]);
		free(this->SMatrix.VAL);
	}
}

void CMatrix::swap(FLOAT &a, FLOAT& b)
{
	if (a != b)
	{
		FLOAT temp = a;
		a = b;
		b = temp;
	}
}

void CMatrix::setMat(const CMatrix &M, const int32_t x, const int32_t y)
{
	if (x<0||y<0||(x+M.SMatrix.col>this->SMatrix.col)||(y+M.SMatrix.row>this->SMatrix.row))
	{
		std::cout << "setMat Beyond the original matrix size or the starting point is wrong!" << std::endl;
		std::exit(0);
	}
	for (int32_t i = 0; i < M.SMatrix.row; i++)
	{
		for (int32_t j = 0; j < M.SMatrix.col; j++)
		{
			this->SMatrix.VAL[y + i][x + j] = M.SMatrix.VAL[i][j];
		}
	}
}

void CMatrix::clearMat()
{
	for (int32_t i = 0; i < this->SMatrix.row; i++)
	{
		for (int32_t j = 0; j < this->SMatrix.col; j++)
			this->SMatrix.VAL[i][j] = 0;
	}
}

void CMatrix::eyes()
{
	clearMat();
	for (int32_t i = 0; i < this->SMatrix.row; i++)
		for (int32_t j = 0; j < this->SMatrix.col; j++)
			if (i==j)
				this->SMatrix.VAL[i][j] = 1;

}

void CMatrix::printfMat()
{
	for (int32_t i = 0; i < this->SMatrix.row; i++)
	{
		for (int32_t j = 0; j < this->SMatrix.col; j++)
		{
			std::cout << this->SMatrix.VAL[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

void CMatrix::ones()
{
	clearMat();
	for (int32_t i = 0; i < this->SMatrix.row; i++)
	{
		for (int32_t j = 0; j < this->SMatrix.col;j++)
		{
			this->SMatrix.VAL[i][j] = 1;
		}
	}
}

void CMatrix::LU(CMatrix& L, CMatrix& U, CMatrix& P)
{
	//注解说明部分：LU分解
	CMatrix tempMat(*this);

	CMatrix temprow(1, tempMat.SMatrix.col);
	CMatrix temprowP(1, tempMat.SMatrix.col);
	L.eyes();


	for (int32_t i = 0; i < tempMat.SMatrix.row; i++)
	{
		//交换列住元
		FLOAT big = tempMat.SMatrix.VAL[i][i];
		int32_t sub = i;
		for (int32_t k = i; k < tempMat.SMatrix.col; k++)
			if (big < tempMat.SMatrix.VAL[k][i])
			{
				big = tempMat.SMatrix.VAL[k][i];
				sub = k;
			}
		if (sub != i)
		{
			//tempMat.printfMat();
			temprow = tempMat.getrow(i);
			temprowP = P.getrow(i);
			//temprow.printfMat();
			tempMat.setMat(tempMat.getrow(sub), 0, i);
			P.setMat(P.getrow(sub), 0, i);

			tempMat.setMat(temprow, 0, sub);
			P.setMat(temprowP, 0, sub);

		}
		//tempMat.printfMat();
	}
	for (int32_t i = 0; i < tempMat.SMatrix.row; i++)
	{

		for (int32_t j = 0; j < tempMat.SMatrix.col; j++)
		{


			if (i == 0)
			{
				U.SMatrix.VAL[i][j] = tempMat.SMatrix.VAL[i][j];
				L.SMatrix.VAL[j][i] = tempMat.SMatrix.VAL[j][i] / U.SMatrix.VAL[0][0];
				//std::cout << U.SMatrix.VAL[i][j]<<std::endl;
				//std::cout << L.SMatrix.VAL[j][i]<<std::endl;
			}
			else
			{
				if (j >= i)
				{
					//L.getrow(i).printfMat();
					//U.getcol(j).printfMat();
					//(L.getrow(i)*U.getcol(j)).printfMat();
					U.SMatrix.VAL[i][j] = tempMat.SMatrix.VAL[i][j] - (L.getrow(i)*U.getcol(j)).sumMat();
					//std::cout << "U=" << U.SMatrix.VAL[i][j] << std::endl;
					if (j != tempMat.SMatrix.col - 1)
					{
						L.SMatrix.VAL[j + 1][i] = (tempMat.SMatrix.VAL[j + 1][i] - (L.getrow(j + 1)
							*U.getcol(i)).SMatrix.VAL[0][0]) / U.SMatrix.VAL[i][i];
						//std::cout << "L=" << L.SMatrix.VAL[j+1][i] << std::endl;
					}
				}
			}
			if (j >= i)
				this->SMatrix.VAL[i][j] = U.SMatrix.VAL[i][j];
			else
				this->SMatrix.VAL[i][j] = L.SMatrix.VAL[i][j];
		}
	}
	//std::cout << "===============" << std::endl;
	//U.printfMat();
	//L.printfMat();
	//(L*U).printfMat();

}

bool CMatrix::JudgelineSegmentIntersection(const CMatrix& A, const CMatrix& B, CMatrix& out_)
{
	//JudgelineSegmentIntersection注解见readme
	//input x y;x y
	assert(A.SMatrix.row == A.SMatrix.col);
	assert(B.SMatrix.row == B.SMatrix.col);
	assert(out_.SMatrix.row == 1 && out_.SMatrix.col == 2);

	CMatrix tempA(2, 2);
	CMatrix tempB(2, 1);
	tempA.SMatrix.VAL[0][0] = A.SMatrix.VAL[1][0] - A.SMatrix.VAL[0][0];
	tempA.SMatrix.VAL[0][1] = -(B.SMatrix.VAL[1][0] - B.SMatrix.VAL[0][0]);
	tempA.SMatrix.VAL[1][0] = A.SMatrix.VAL[1][1] - A.SMatrix.VAL[0][1];
	tempA.SMatrix.VAL[1][1] = -(B.SMatrix.VAL[1][1] - B.SMatrix.VAL[0][1]);
	tempB.SMatrix.VAL[0][0] = B.SMatrix.VAL[0][0] - A.SMatrix.VAL[0][0];
	tempB.SMatrix.VAL[1][0] = B.SMatrix.VAL[0][1] - A.SMatrix.VAL[0][1];
	if (tempA.SMatrix.VAL[1][1] * tempA.SMatrix.VAL[0][0] == 
		tempA.SMatrix.VAL[0][1]*tempA.SMatrix.VAL[1][0])
	{
		std::cout << "Don't have any solve!These lines can't have cross" << std::endl;
		return false;
		exit(0);
	}
	CMatrix tempout_ = solve(tempA, tempB);
	if (tempout_.SMatrix.VAL[0][0] >= 0 && tempout_.SMatrix.VAL[0][0] <= 1
		&& tempout_.SMatrix.VAL[1][0] >= 0 && tempout_.SMatrix.VAL[1][0] <= 0)
	{
		std::cout << "find any cross on the Line extension cord !" << std::endl;
		return false;
		exit(0);
	}
	out_.SMatrix.VAL[0][0] = A.SMatrix.VAL[0][0] + tempout_.SMatrix.VAL[0][0] * (A.SMatrix.VAL[1][0] - A.SMatrix.VAL[0][0]);
	out_.SMatrix.VAL[0][1] = A.SMatrix.VAL[0][1] + tempout_.SMatrix.VAL[0][0] * (A.SMatrix.VAL[1][1] - A.SMatrix.VAL[0][1]);
	return true;
}

CMatrix& CMatrix::t()
{
	assert(this->SMatrix.row == this->SMatrix.col);

	for (int32_t i =0; i < this->SMatrix.row; i++)
		for (int32_t j = 0; j < this->SMatrix.col;j++)
			if (j>i)
				swap(this->SMatrix.VAL[i][j], this->SMatrix.VAL[j][i]);
	return *this;
}

FLOAT CMatrix::det()
{
	assert(this->SMatrix.row == this->SMatrix.col);
	FLOAT detData = 1.0;
	int32_t changeNum = 0;
	this->LU(changeNum);
	for (int32_t i = 0; i < this->SMatrix.row;i++)
		detData *= this->SMatrix.VAL[i][i];
	if (changeNum % 2 == 0)
		return detData;
	else
		return (-detData);

}

FLOAT CMatrix::sumMat()
{
	FLOAT sum = 0.0;
	for (int32_t i = 0;i< this->SMatrix.row; i++)
	{
		for (int32_t j = 0; j < this->SMatrix.col; j++)
		{
			sum += this->SMatrix.VAL[i][j];
		}
	}
	return sum;
}

FLOAT CMatrix::sqatMat()
{
	FLOAT data = 0.0;
	for (int32_t row = 0; row < this->SMatrix.row; row++)
		for (int32_t col = 0; col < this->SMatrix.col; col++)
			data += std::pow(this->SMatrix.VAL[row][col],2);
	
	return std::sqrt(data);

}

FLOAT CMatrix::meanMat()
{
	FLOAT data = 0.0;
	int32_t n = 0;
	for (int32_t row = 0; row < this->SMatrix.row; row++)
		for (int32_t col = 0; col < this->SMatrix.col; col++)
		{
			data += this->SMatrix.VAL[row][col];
			n++;
		}

	return data/n;
}

CMatrix& CMatrix::invMat()
{
	CMatrix L(this->SMatrix.row, this->SMatrix.col);
	CMatrix U(this->SMatrix.row, this->SMatrix.col);
	CMatrix Linv(this->SMatrix.row, this->SMatrix.row);
	CMatrix Uinv(this->SMatrix.row, this->SMatrix.row);
	CMatrix E(this->SMatrix.row, this->SMatrix.row);
	CMatrix P(this->SMatrix.row, this->SMatrix.row);
	E.eyes();
	P.eyes();
	this->LU(L, U, P);
	////===============debug
	//P.printfMat();
	//L.printfMat();
	//U.printfMat();
	int32_t row = this->SMatrix.row;
	int32_t col = this->SMatrix.col;
	for (int32_t i = 0; i < row; i++)
		for (int32_t j = 0; j < col; j++)
		{//E.ij - getrow(A)*getcol(A - 1) / A.ii;
			Linv.SMatrix.VAL[i][j] = (E.SMatrix.VAL[i][j] - (L.getrow(i)*Linv.getcol(j)).sumMat()) / L.SMatrix.VAL[i][i];
			Uinv.SMatrix.VAL[row - i - 1][col - j - 1] = (E.SMatrix.VAL[row - i - 1][col - j - 1]
				- (U.getrow(row - i - 1)*Uinv.getcol(col - j - 1)).sumMat()) / U.SMatrix.VAL[row - i - 1][row - i - 1];
	
			//std::cout << "=========================" << std::endl;
			//Linv.printfMat();
			//Uinv.printfMat();
			//std::cout << "=========================" << std::endl;
		}

	return (*this = Uinv*Linv*P);
}

CMatrix CMatrix::getMat(const int32_t row1, const int32_t row2, const int32_t col1, const int32_t col2)
{
	if (row1 > row2 || col1 > col2 || row2 == 0 || col2 == 0)
	{
		std::cout << "first input must higher second!" << std::endl;
		std::exit(0);
	}
	CMatrix M(row2 - row1 + 1, col2 - col1 + 1);
	for (int32_t row = 0; row < M.SMatrix.row; row++)
	{
		for (int32_t col = 0; col < M.SMatrix.col; col++)
		{
			M.SMatrix.VAL[row][col] = this->SMatrix.VAL[row1 + row][col1 + col];
		}
	}
	return M;
}

CMatrix CMatrix::ones(const int32_t row, const int32_t col)
{
	assert(row >= 0 && col >= 0);

	CMatrix M(row, col);
	for (int32_t i = 0; i < M.SMatrix.row; i++)
		for (int32_t j = 0; j < M.SMatrix.col; j++)
			M.SMatrix.VAL[i][j] = 1;

	return M;
}

CMatrix CMatrix::eyes(const int32_t row, const int32_t col)
{
	assert(row >= 0 && col >= 0);

	CMatrix M(row, col);
	M.eyes();
	return M;
}

CMatrix CMatrix::RotMatX(const FLOAT angle)
{
	FLOAT Angle = angle*3.141592635 / 180;
	CMatrix R(3, 3);
	FLOAT C = std::cos(Angle);
	FLOAT S = std::sin(Angle);
	R.SMatrix.VAL[0][0] = 1;
	R.SMatrix.VAL[1][1] = C;
	R.SMatrix.VAL[1][2] = -S;
	R.SMatrix.VAL[2][1] = S;
	R.SMatrix.VAL[2][2] = C;
	return R;
}

CMatrix CMatrix::RotMatY(const FLOAT angle)
{
	FLOAT Angle = angle*3.141592635 / 180;
	CMatrix R(3, 3);
	FLOAT C = std::cos(Angle);
	FLOAT S = std::sin(Angle);
	R.SMatrix.VAL[0][0] = C;
	R.SMatrix.VAL[0][2] = S;
	R.SMatrix.VAL[1][1] = 1;
	R.SMatrix.VAL[2][0] = -S;
	R.SMatrix.VAL[2][2] = C;
	return R;
}

CMatrix CMatrix::RotMatZ(const FLOAT angle)
{
	FLOAT Angle = angle*3.141592635 / 180;
	CMatrix R(3, 3);
	FLOAT C = std::cos(Angle);
	FLOAT S = std::sin(Angle);
	R.SMatrix.VAL[0][0] = C;
	R.SMatrix.VAL[0][1] = -S;
	R.SMatrix.VAL[1][0] = S;
	R.SMatrix.VAL[1][1] = C;
	R.SMatrix.VAL[2][2] = 1;
	return R;
}

CMatrix CMatrix::getrow(const int32_t row)
{
	assert(row <= this->SMatrix.row);

	CMatrix M(1, this->SMatrix.col);
	for (int32_t j = 0; j < this->SMatrix.col; j++)
	{
		M.SMatrix.VAL[0][j] = this->SMatrix.VAL[row][j];
	}
	return M;
}

CMatrix CMatrix::getcol(const int32_t col)
{
	assert(col <= this->SMatrix.col);

	CMatrix M(this->SMatrix.col, 1);
	for (int32_t i = 0; i < this->SMatrix.row; i++)
	{
		M.SMatrix.VAL[i][0] = this->SMatrix.VAL[i][col];
	}
	return M;
}

CMatrix CMatrix::LU()
{
	//注解说明部分：LU分解
	CMatrix tempMat(*this);
	CMatrix L(tempMat.SMatrix.row, tempMat.SMatrix.col);
	CMatrix U(tempMat.SMatrix.row, tempMat.SMatrix.col);
	CMatrix temprow(1, tempMat.SMatrix.col);
	L.eyes();
	

	for (int32_t i = 0; i < tempMat.SMatrix.row;i++)
	{
		//交换列住元
		FLOAT big = tempMat.SMatrix.VAL[i][i];
		int32_t sub = i;
		for (int32_t k = i; k < tempMat.SMatrix.col; k++)
			if (big < tempMat.SMatrix.VAL[k][i])
			{
				big = tempMat.SMatrix.VAL[k][i];
				sub = k;
			}
		if (sub != i)
		{
			//tempMat.printfMat();
			temprow = tempMat.getrow(i);
			//temprow.printfMat();
			tempMat.setMat(tempMat.getrow(sub), 0, i);
			tempMat.setMat(temprow, 0, sub);
		}
		//tempMat.printfMat();
	}
	for (int32_t i = 0; i < tempMat.SMatrix.row;i++)
	{

		for (int32_t j = 0; j < tempMat.SMatrix.col;j++)
		{
			

			if (i == 0)
			{
				U.SMatrix.VAL[i][j] = tempMat.SMatrix.VAL[i][j];
				L.SMatrix.VAL[j][i] = tempMat.SMatrix.VAL[j][i] / U.SMatrix.VAL[0][0];
				//std::cout << U.SMatrix.VAL[i][j]<<std::endl;
				//std::cout << L.SMatrix.VAL[j][i]<<std::endl;
			}
			else
			{
				if (j >= i)
				{
					//L.getrow(i).printfMat();
					//U.getcol(j).printfMat();
					//(L.getrow(i)*U.getcol(j)).printfMat();
					U.SMatrix.VAL[i][j] = tempMat.SMatrix.VAL[i][j] - (L.getrow(i)*U.getcol(j)).sumMat();
					//std::cout << "U=" << U.SMatrix.VAL[i][j] << std::endl;
					if (j != tempMat.SMatrix.col - 1)
					{
						L.SMatrix.VAL[j + 1][i] = (tempMat.SMatrix.VAL[j + 1][i] - (L.getrow(j + 1)
							*U.getcol(i)).SMatrix.VAL[0][0]) / U.SMatrix.VAL[i][i];
						//std::cout << "L=" << L.SMatrix.VAL[j+1][i] << std::endl;
					}
				}
			}
			if (j >= i)
				this->SMatrix.VAL[i][j] = U.SMatrix.VAL[i][j];
			else
				this->SMatrix.VAL[i][j] = L.SMatrix.VAL[i][j];
		}
	}
	//U.printfMat();
	//L.printfMat();
	//(L*U).printfMat();
	return *this;
}

CMatrix CMatrix::LU(int32_t& changeNum)
{
	//注解说明部分：LU分解
	CMatrix tempMat(*this);
	CMatrix L(tempMat.SMatrix.row, tempMat.SMatrix.col);
	CMatrix U(tempMat.SMatrix.row, tempMat.SMatrix.col);
	CMatrix temprow(1, tempMat.SMatrix.col);
	L.eyes();


	for (int32_t i = 0; i < tempMat.SMatrix.row; i++)
	{
		//交换列住元
		FLOAT big = tempMat.SMatrix.VAL[i][i];
		int32_t sub = i;
		for (int32_t k = i; k < tempMat.SMatrix.col; k++)
			if (big < tempMat.SMatrix.VAL[k][i])
			{
				big = tempMat.SMatrix.VAL[k][i];
				sub = k;
			}
		if (sub != i)
		{
			//tempMat.printfMat();
			temprow = tempMat.getrow(i);
			//temprow.printfMat();
			tempMat.setMat(tempMat.getrow(sub), 0, i);
			tempMat.setMat(temprow, 0, sub);
			changeNum++;
		}
		//tempMat.printfMat();
	}
	for (int32_t i = 0; i < tempMat.SMatrix.row; i++)
	{

		for (int32_t j = 0; j < tempMat.SMatrix.col; j++)
		{


			if (i == 0)
			{
				U.SMatrix.VAL[i][j] = tempMat.SMatrix.VAL[i][j];
				L.SMatrix.VAL[j][i] = tempMat.SMatrix.VAL[j][i] / U.SMatrix.VAL[0][0];
				//std::cout << U.SMatrix.VAL[i][j]<<std::endl;
				//std::cout << L.SMatrix.VAL[j][i]<<std::endl;
			}
			else
			{
				if (j >= i)
				{
					//L.getrow(i).printfMat();
					//U.getcol(j).printfMat();
					//(L.getrow(i)*U.getcol(j)).printfMat();
					U.SMatrix.VAL[i][j] = tempMat.SMatrix.VAL[i][j] - (L.getrow(i)*U.getcol(j)).sumMat();
					//std::cout << "U=" << U.SMatrix.VAL[i][j] << std::endl;
					if (j != tempMat.SMatrix.col - 1)
					{
						L.SMatrix.VAL[j + 1][i] = (tempMat.SMatrix.VAL[j + 1][i] - (L.getrow(j + 1)
							*U.getcol(i)).SMatrix.VAL[0][0]) / U.SMatrix.VAL[i][i];
						//std::cout << "L=" << L.SMatrix.VAL[j+1][i] << std::endl;
					}
				}
			}
			if (j >= i)
				this->SMatrix.VAL[i][j] = U.SMatrix.VAL[i][j];
			else
				this->SMatrix.VAL[i][j] = L.SMatrix.VAL[i][j];
		}
	}
	std::cout << "===============" << std::endl;
	//U.printfMat();
	//L.printfMat();
	(L*U).printfMat();
	return *this;
}

CMatrix CMatrix::pow(CMatrix& M, int32_t n)
{
	if (n==0)
	{
		return CMatrix::eyes(M.SMatrix.row,M.SMatrix.col);
	}
	return M * CMatrix::pow(M, --n);
}

CMatrix CMatrix::solve(CMatrix& A, const CMatrix& C)
{
	if (A.SMatrix.row!=A.SMatrix.col)
	{
		std::cout << "this solve can't solve overdetermined equation " << std::endl;
		exit(0);
	}
	assert(A.SMatrix.col == C.SMatrix.row);

	//采用智能指针，自动释放内存，省去在main函数调用时候，再次free
	//同时函数不可以采用带引用的返回，因为返回前share_prt已经将指针释放，引用为未知。
	std::shared_ptr<CMatrix> B(new CMatrix(C.SMatrix.row, C.SMatrix.col));
	A.invMat();
	*B = A*C;
	return *B;
}

FLOAT CMatrix::pythag(FLOAT a, FLOAT b) {
	FLOAT absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb)
		return absa*sqrt(1.0 + SQR(absb / absa));
	else
		return (absb == 0.0 ? 0.0 : absb*sqrt(1.0 + SQR(absa / absb)));
}

void CMatrix::svd(CMatrix &U2, CMatrix &W, CMatrix &V) {

	CMatrix U = CMatrix(*this);
	U2 = CMatrix(this->SMatrix.row, this->SMatrix.row);
	V = CMatrix(this->SMatrix.col, this->SMatrix.col);
	int32_t m = this->SMatrix.row;
	int32_t n = this->SMatrix.col;
	FLOAT* w = (FLOAT*)malloc(this->SMatrix.col*sizeof(FLOAT));
	FLOAT* rv1 = (FLOAT*)malloc(this->SMatrix.col*sizeof(FLOAT));

	int32_t flag, i, its, j, jj, k, l, nm;
	FLOAT   anorm, c, f, g, h, s, scale, x, y, z;

	g = scale = anorm = 0.0; // Householder reduction to bidiagonal form.
	for (i = 0; i < n; i++) {
		l = i + 1;
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		if (i < m) {
			for (k = i; k < m; k++) scale += fabs(U.SMatrix.VAL[k][i]);
			if (scale) {
				for (k = i; k < m; k++) {
					U.SMatrix.VAL[k][i] /= scale;
					s += U.SMatrix.VAL[k][i] * U.SMatrix.VAL[k][i];
				}
				f = U.SMatrix.VAL[i][i];
				g = -SIGN(sqrt(s), f);
				h = f*g - s;
				U.SMatrix.VAL[i][i] = f - g;
				for (j = l; j < n; j++) {
					for (s = 0.0, k = i; k < m; k++) s += U.SMatrix.VAL[k][i] * U.SMatrix.VAL[k][j];
					f = s / h;
					for (k = i; k < m; k++) U.SMatrix.VAL[k][j] += f*U.SMatrix.VAL[k][i];
				}
				for (k = i; k < m; k++) U.SMatrix.VAL[k][i] *= scale;
			}
		}
		w[i] = scale*g;
		g = s = scale = 0.0;
		if (i < m && i != n - 1) {
			for (k = l; k < n; k++) scale += fabs(U.SMatrix.VAL[i][k]);
			if (scale) {
				for (k = l; k < n; k++) {
					U.SMatrix.VAL[i][k] /= scale;
					s += U.SMatrix.VAL[i][k] * U.SMatrix.VAL[i][k];
				}
				f = U.SMatrix.VAL[i][l];
				g = -SIGN(sqrt(s), f);
				h = f*g - s;
				U.SMatrix.VAL[i][l] = f - g;
				for (k = l; k < n; k++) rv1[k] = U.SMatrix.VAL[i][k] / h;
				for (j = l; j < m; j++) {
					for (s = 0.0, k = l; k < n; k++) s += U.SMatrix.VAL[j][k] * U.SMatrix.VAL[i][k];
					for (k = l; k < n; k++) U.SMatrix.VAL[j][k] += s*rv1[k];
				}
				for (k = l; k < n; k++) U.SMatrix.VAL[i][k] *= scale;
			}
		}
		anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}
	for (i = n - 1; i >= 0; i--) { // Accumulation of right-hand transformations.
		if (i < n - 1) {
			if (g) {
				for (j = l; j < n; j++) // Double division to avoid possible underflow.
					V.SMatrix.VAL[j][i] = (U.SMatrix.VAL[i][j] / U.SMatrix.VAL[i][l]) / g;
				for (j = l; j < n; j++) {
					for (s = 0.0, k = l; k < n; k++) s += U.SMatrix.VAL[i][k] * V.SMatrix.VAL[k][j];
					for (k = l; k < n; k++) V.SMatrix.VAL[k][j] += s*V.SMatrix.VAL[k][i];
				}
			}
			for (j = l; j < n; j++) V.SMatrix.VAL[i][j] = V.SMatrix.VAL[j][i] = 0.0;
		}
		V.SMatrix.VAL[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i = IMIN(m, n) - 1; i >= 0; i--) { // Accumulation of left-hand transformations.
		l = i + 1;
		g = w[i];
		for (j = l; j < n; j++) U.SMatrix.VAL[i][j] = 0.0;
		if (g) {
			g = 1.0 / g;
			for (j = l; j < n; j++) {
				for (s = 0.0, k = l; k < m; k++) s += U.SMatrix.VAL[k][i] * U.SMatrix.VAL[k][j];
				f = (s / U.SMatrix.VAL[i][i])*g;
				for (k = i; k < m; k++) U.SMatrix.VAL[k][j] += f*U.SMatrix.VAL[k][i];
			}
			for (j = i; j < m; j++) U.SMatrix.VAL[j][i] *= g;
		}
		else for (j = i; j < m; j++) U.SMatrix.VAL[j][i] = 0.0;
		++U.SMatrix.VAL[i][i];
	}
	for (k = n - 1; k >= 0; k--) { // Diagonalization of the bidiagonal form: Loop over singular values,
		for (its = 0; its < 30; its++) { // and over allowed iterations.
			flag = 1;
			for (l = k; l >= 0; l--) { // Test for splitting.
				nm = l - 1;
				if ((FLOAT)(fabs(rv1[l]) + anorm) == anorm) { flag = 0; break; }
				if ((FLOAT)(fabs(w[nm]) + anorm) == anorm) { break; }
			}
			if (flag) {
				c = 0.0; // Cancellation of rv1[l], if l > 1.
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((FLOAT)(fabs(f) + anorm) == anorm) break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g*h;
					s = -f*h;
					for (j = 0; j < m; j++) {
						y = U.SMatrix.VAL[j][nm];
						z = U.SMatrix.VAL[j][i];
						U.SMatrix.VAL[j][nm] = y*c + z*s;
						U.SMatrix.VAL[j][i] = z*c - y*s;
					}
				}
			}
			z = w[k];
			if (l == k) { // Convergence.
				if (z < 0.0) { // Singular value is made nonnegative.
					w[k] = -z;
					for (j = 0; j < n; j++) V.SMatrix.VAL[j][k] = -V.SMatrix.VAL[j][k];
				}
				break;
			}
			if (its == 29)
				std::cout << "ERROR in SVD: No convergence in 30 iterations" << std::endl;
			x = w[l]; // Shift from bottom 2-by-2 minor.
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0*h*y);
			g = pythag(f, 1.0);
			f = ((x - z)*(x + z) + h*((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0; // Next QR transformation:
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x*c + g*s;
				g = g*c - x*s;
				h = y*s;
				y *= c;
				for (jj = 0; jj < n; jj++) {
					x = V.SMatrix.VAL[jj][j];
					z = V.SMatrix.VAL[jj][i];
					V.SMatrix.VAL[jj][j] = x*c + z*s;
					V.SMatrix.VAL[jj][i] = z*c - x*s;
				}
				z = pythag(f, h);
				w[j] = z; // Rotation can be arbitrary if z = 0.
				if (z) {
					z = 1.0 / z;
					c = f*z;
					s = h*z;
				}
				f = c*g + s*y;
				x = c*y - s*g;
				for (jj = 0; jj < m; jj++) {
					y = U.SMatrix.VAL[jj][j];
					z = U.SMatrix.VAL[jj][i];
					U.SMatrix.VAL[jj][j] = y*c + z*s;
					U.SMatrix.VAL[jj][i] = z*c - y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}

	// sort singular values and corresponding columns of u and v
	// by decreasing magnitude. Also, signs of corresponding columns are
	// flipped so as to maximize the number of positive elements.
	int32_t s2, inc = 1;
	FLOAT   sw;
	FLOAT* su = (FLOAT*)malloc(m*sizeof(FLOAT));
	FLOAT* sv = (FLOAT*)malloc(n*sizeof(FLOAT));
	do { inc *= 3; inc++; } while (inc <= n);
	do {
		inc /= 3;
		for (i = inc; i < n; i++) {
			sw = w[i];
			for (k = 0; k < m; k++) su[k] = U.SMatrix.VAL[k][i];
			for (k = 0; k < n; k++) sv[k] = V.SMatrix.VAL[k][i];
			j = i;
			while (w[j - inc] < sw) {
				w[j] = w[j - inc];
				for (k = 0; k < m; k++) U.SMatrix.VAL[k][j] = U.SMatrix.VAL[k][j - inc];
				for (k = 0; k < n; k++) V.SMatrix.VAL[k][j] = V.SMatrix.VAL[k][j - inc];
				j -= inc;
				if (j < inc) break;
			}
			w[j] = sw;
			for (k = 0; k < m; k++) U.SMatrix.VAL[k][j] = su[k];
			for (k = 0; k < n; k++) V.SMatrix.VAL[k][j] = sv[k];
		}
	} while (inc > 1);
	for (k = 0; k < n; k++) { // flip signs
		s2 = 0;
		for (i = 0; i < m; i++) if (U.SMatrix.VAL[i][k] < 0.0) s2++;
		for (j = 0; j < n; j++) if (V.SMatrix.VAL[j][k] < 0.0) s2++;
		if (s2 > (m + n) / 2) {
			for (i = 0; i < m; i++) U.SMatrix.VAL[i][k] = -U.SMatrix.VAL[i][k];
			for (j = 0; j < n; j++) V.SMatrix.VAL[j][k] = -V.SMatrix.VAL[j][k];
		}
	}

	// create vector and copy singular values
	W = CMatrix(std::min(m, n), 1, w);

	// extract mxm submatrix U
	U2.setMat(U.getMat(0, 0, m - 1, std::min(m - 1, n - 1)), 0, 0);

	// release temporary memory
	free(w);
	free(rv1);
	free(su);
	free(sv);
}

CMatrix& CMatrix::operator=(const CMatrix& M)
	{
		if (this != &M)//判断本地的CMatrix是否等于传入的M
		{
			if (M.SMatrix.row != this->SMatrix.row || M.SMatrix.col != this->SMatrix.col)
			{
				releaseMemory();
				allocateMemory(M.SMatrix.row, M.SMatrix.col);
				this->setMat(M, 0, 0);
			}
			else
			{
				this->setMat(M, 0, 0);
			}
		}

		return *this;
	}

CMatrix& CMatrix::operator+(const FLOAT M)
{
	for (int32_t i = 0; i < this->SMatrix.row; i++)
	{
		for (int32_t j=0; j < this->SMatrix.col; j++)
		{
			this->SMatrix.VAL[i][j] += M;
		}
	}
	return *this;
}

CMatrix& CMatrix::operator+(const CMatrix& M)
{
	assert(this->SMatrix.row == M.SMatrix.row&&this->SMatrix.col == M.SMatrix.col);

	for (int32_t i = 0; i < this->SMatrix.row; i++)
	{
		for (int32_t j = 0; j < this->SMatrix.col; j++)
		{
			this->SMatrix.VAL[i][j] += M.SMatrix.VAL[i][j];
		}
	}
	return *this;
}

CMatrix& CMatrix::operator-(const FLOAT M)
{
	for (int32_t i = 0; i < this->SMatrix.row; i++)
	{
		for (int32_t j=0; j < this->SMatrix.col; j++)
		{
			this->SMatrix.VAL[i][j] -= M;
		}
	}
	return *this;
}

CMatrix& CMatrix::operator-(const CMatrix& M)
{
	assert(this->SMatrix.row == M.SMatrix.row&&this->SMatrix.col == M.SMatrix.col);

	for (int32_t i = 0; i < this->SMatrix.row; i++)
	{
		for (int32_t j = 0; j < this->SMatrix.col; j++)
		{
			this->SMatrix.VAL[i][j] -= M.SMatrix.VAL[i][j];
		}
	}
	return *this;
}

CMatrix CMatrix::operator*(const CMatrix& M)
{
	assert(this->SMatrix.col == M.SMatrix.row);

	const CMatrix &A = *this;
	const CMatrix &B = M;

	CMatrix C(A.SMatrix.row, B.SMatrix.col);
	for (int32_t i = 0; i < A.SMatrix.row; i++)
		for (int32_t j = 0; j < B.SMatrix.col; j++)
			for (int32_t k = 0; k < A.SMatrix.col; k++)
				C.SMatrix.VAL[i][j] += A.SMatrix.VAL[i][k] * B.SMatrix.VAL[k][j];
	return C;
}

CMatrix& CMatrix::operator*(const FLOAT M)
{

	for (int32_t i = 0; i < this->SMatrix.row; i++)
	{
		for (int32_t j = 0; j < this->SMatrix.col; j++)
		{
			this->SMatrix.VAL[i][j] *= M;
		}
	}
	return *this;
}

CMatrix& CMatrix::operator/(const FLOAT M)
{
	assert(M != 0);

	for (int32_t i = 0; i < this->SMatrix.row; i++)
	{
		for (int32_t j = 0; j < this->SMatrix.col; j++)
		{
			this->SMatrix.VAL[i][j] /= M;
		}
	}
	return *this;
}

CMatrix& CMatrix::operator/(const CMatrix& M)
{
	assert(this->SMatrix.col == M.SMatrix.row);
	return *this;
	
}