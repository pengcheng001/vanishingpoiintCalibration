#include "Matrix.h"
int main()
{
	//while (1)
	//{
	//	//*一般就用引用类型去接收，或者就使用了引用的作用，如果用非引用类型接受，
	//	//*就等于将函数返回的引用的数据值，复制给了该接收对象，和函数返回非引用类型是一样的效果。
	//	double val[] = { 4, 2, 1, 5, 8, 7, 2, 10, 4, 8, 3, 6, 6, 8, 4, 9 };
	//	CMatrix CM(4, 4, val);
	//	CMatrix CM1(4, 4, val);
	//	CM1.t().printfMat();

	//	double A[] = { 1, 2, 2, 2 };
	//	double C[] = { 1, 2 };
	//	CMatrix A1(2, 2, A);
	//	CMatrix C1(2, 1, C);
	//	CMatrix B = CMatrix::solve(A1, C1);
	//	B.printfMat();
	//	std::cout << "============" << std::endl;
	//}
	CMatrix A(2,2);
	A.SMatrix.VAL[0][0] = 0;


}