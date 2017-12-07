#ifndef _LEVENBERG_GAUSS_NEWTON_H
#define  _LEVENBERG_GAUSS_NEWTON_H

#include "headGather.h"

class CLevenberGaussNewton
{
public:
	CLevenberGaussNewton(cv::Mat IterativeParameter,const int Max_iter)
	{
		this->IterativeParameter = IterativeParameter;
		this->DERIV_STEP = 1e-5;  
		this->Max_iter = Max_iter;
		this->dx = 0.009;
		this->dy = 0.009;
	};
	~CLevenberGaussNewton(){};

public:
	friend class CVirVainshingPoint;
	bool LMNewtonIterative(cv::Mat &IterativeParameter,std::vector<std::vector<cv::Mat>> inputData,cv::Mat& outputData);
	bool NewtonIterative(cv::Mat &IterativeParameter,std::vector<std::vector<cv::Mat>> &inputData,cv::Mat& outputData);
	void LMNewtonIterativeRunExam(std::vector<std::vector<cv::Mat>> &inputData,cv::Mat &outputData);
	void NewtonIterativeRunExam(std::vector<std::vector<cv::Mat>> &inputData,cv::Mat &outputData);
private:
	cv::Mat IterativeParameter;
	cv::Mat outputData;  
	std::vector<std::vector<cv::Mat>> inputData;  
private:
	double dx,dy;
	int Max_iter;
	double DERIV_STEP;
	double DerivativeSolving(cv::Mat IterativeParameter,const std::vector<cv::Mat> FunCinputData,const int n);
	//==========根据实际需要拟合的函数修改========
	double NewtonFunC(const cv::Mat IterativeParameter,std::vector<cv::Mat> FunCinputData);
	//==========根据实际需要拟合的函数修改========

};



#endif   //LevenbergGaussNewton.h



