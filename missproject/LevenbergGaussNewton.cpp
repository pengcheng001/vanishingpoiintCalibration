#include "headGather.h"


double CLevenberGaussNewton::NewtonFunC(const cv::Mat IterativeParameter, std::vector<cv::Mat> FunCinputData)
{
	//IterativeParameter: u0 ,v0 , fx ,fy;
	//根据实际需要拟合的函数修改



			double U[2] ={0},V[2]={0},L[2];
			double U_0 = IterativeParameter.at<double>(0,0);
			double V_0 = IterativeParameter.at<double>(1,0);
			double F_X = IterativeParameter.at<double>(2,0);
			double F_Y = IterativeParameter.at<double>(3,0);
			double K = IterativeParameter.at<double>(4,0);
			int PointNum = FunCinputData.size();
			std::vector<cv::Mat> temp(PointNum);

			for (int j=0;j<PointNum;j++)
			{
				double BiasU = FunCinputData[j].at<double>(0,0) ;
				double BiasV = FunCinputData[j].at<double>(1,0) ;
				double Media =  K*(std::pow((BiasU-U_0),2)+std::pow((BiasV-V_0),2));
				temp[j] = (cv::Mat_<double>(2,1)<<BiasU +Media*dx*dx*(BiasU-U_0),BiasV +Media*dy*dy*(BiasV-V_0));
				U[j] = temp[j].at<double>(0,0) ;
				V[j] = temp[j].at<double>(1,0) ;
				double tempX = (U[j] - U_0)/ F_X;
				double tempY = (V[j] - V_0)/ F_Y;
				double tempZ = 1.0;
				L[j] = std::sqrt(std::pow(tempX,2)+std::pow(tempY,2)+1);
			}

			double LConsequence = L[0] * L[1];
			double HConsequence = 1 + (U[0] - U_0)*(U[1] - U_0)/std::pow(F_X,2) + (V[0] - V_0)*(V[1] - V_0)/std::pow(F_Y,2);
 			return HConsequence/LConsequence;
		
	//根据实际需要拟合的函数修改
}


double CLevenberGaussNewton::DerivativeSolving(cv::Mat IterativeParameter,const std::vector<cv::Mat> FunCinputData,const int n)
{
	cv::Mat TempParameter1 = IterativeParameter.clone();
	cv::Mat TempParameter2 = IterativeParameter.clone();

	TempParameter1.at<double>(n,0) = TempParameter1.at<double>(n,0) + DERIV_STEP;
	TempParameter2.at<double>(n,0) = TempParameter2.at<double>(n,0) - DERIV_STEP;
	double d = (NewtonFunC(TempParameter1,FunCinputData) - NewtonFunC(TempParameter2,FunCinputData))/(2*DERIV_STEP);
	return d;
}

bool CLevenberGaussNewton::LMNewtonIterative(cv::Mat &IterativeParameter,std::vector<std::vector<cv::Mat>> inputData,cv::Mat& outputData)
 {
	   double bestErr=0.0;
	int inputDataPos = inputData.size();

	int ParameterNum = IterativeParameter.rows;
	cv::Mat J(cv::Size(ParameterNum,inputDataPos),CV_64F);
	cv::Mat R(cv::Size(1,inputDataPos),CV_64F);    //残差
	cv::Mat Rtemp(cv::Size(1,inputDataPos),CV_64F);
	cv::Mat I =cv::Mat::ones(ParameterNum,ParameterNum,CV_64F);
	double  u=1,v=2;
	for (int iter =0 ; iter < Max_iter ;iter++)
	{
		
		double TotalErr=0.0,TotalErrTemp=0.0;
		for (int i=0;i<inputDataPos;i++)
		{
			std::vector<cv::Mat> tempInputData = inputData[i];
			//std::cout << inputData[0][0] <<std::endl;
			//std::cout << inputData[0][1] <<std::endl;
			double tempOutputData = 0;
			R.at<double>(i,0) = tempOutputData - NewtonFunC(IterativeParameter,tempInputData);
			//std::cout << inputData[0][0] <<std::endl;
			//std::cout << inputData[0][1] <<std::endl;
			TotalErr+=std::pow(R.at<double>(i,0),2);
			for (int j=0;j<ParameterNum;j++)
			{
				J.at<double>(i,j) = DerivativeSolving(IterativeParameter ,tempInputData , j);
			}
		}
		TotalErr/=inputDataPos;//r^2就是目标函数Fx（0）
		//std::cout<<IterativeParameter<<std::endl;
		IterativeParameter = IterativeParameter + ((J.t()*J + u*I)).inv()*J.t()*R;
		//std::cout<<IterativeParameter<<std::endl;
		cv::Mat hlm = ((J.t()*J + u*I)).inv()*J.t()*R; //搜索步长
		cv::Mat tempParameter = IterativeParameter.clone();

		for (int i=0;i<inputDataPos;i++)
		{
			std::vector<cv::Mat> tempInputData = inputData[i];
			double tempOutputData = 0;
			Rtemp.at<double>(i,0) = tempOutputData - NewtonFunC(tempParameter,tempInputData);
			TotalErrTemp += std::pow(Rtemp.at<double>(i,0),2);
		}
		TotalErrTemp /=inputDataPos; 
		cv::Mat Q(cv::Size(1,1),CV_64F);  
		Q =(TotalErr - TotalErrTemp) / (0.5*hlm.t()*(u*hlm - J.t()*R));
		double q = Q.at<double>(0,0);
		if (q>0)
		{
			double s = 1/3;
			double t = 1 - std::pow((2*q-1),3);
			if (t > s)
			{	u = u*t;
			v = 2;	} 
			else
			{	u = u*s;
			v = 2;   }  
		}
		else
		{
			u = u*v;
			v = 2*v;
		}
		if (abs(bestErr-TotalErr)<1e-9)	
		{
			std::cout << "TotalErr = " <<TotalErr<<std::endl;
			std::cout << "IterativeParameter(1) = " <<IterativeParameter.at<double>(0,0)<<std::endl;
			std::cout << "IterativeParameter(2) = " <<IterativeParameter.at<double>(1,0)<<std::endl;
			std::cout << "IterativeParameter(3) = " <<IterativeParameter.at<double>(2,0)<<std::endl;
			std::cout << "IterativeParameter(4) = " <<IterativeParameter.at<double>(3,0)<<std::endl;
			std::cout << "IterativeParameter(5) = " <<IterativeParameter.at<double>(4,0)<<std::endl;
			return true;
		}

		std::cout << "TotalErr = " <<TotalErr<<std::endl;
		bestErr = TotalErr;
	}
	return false;
}


bool CLevenberGaussNewton::NewtonIterative(cv::Mat &IterativeParameter,std::vector<std::vector<cv::Mat>> &inputData,cv::Mat& outputData)
{
	double bestErr=0.0;
	int inputDataPos = inputData.size();
	int ParameterNum = IterativeParameter.rows;
	cv::Mat J(cv::Size(ParameterNum,inputDataPos),CV_64F);
	cv::Mat R(cv::Size(1,inputDataPos),CV_64F);
	for (int iter =0 ; iter < Max_iter ;iter++)
	{
		double TotalErr=0.0;
		for (int i=0;i<inputDataPos;i++)
		{
			std::vector<cv::Mat> tempInputData = inputData[i];
			double tempOutputData = 0;
			R.at<double>(i,0) = tempOutputData - NewtonFunC(IterativeParameter,tempInputData);
			TotalErr+=std::pow(R.at<double>(i,0),2);
			for (int j=0;j<ParameterNum;j++)
			{
				J.at<double>(i,j) = DerivativeSolving(IterativeParameter ,tempInputData , j);
			}
		}
		TotalErr/=inputDataPos;
		if (abs(bestErr-TotalErr)<1e-5)		
		{
			std::cout << "TotalErr = " <<TotalErr<<std::endl;
			std::cout << "IterativeParameter(1) = " <<IterativeParameter.at<double>(0,0)<<std::endl;
			std::cout << "IterativeParameter(2) = " <<IterativeParameter.at<double>(1,0)<<std::endl;
			std::cout << "IterativeParameter(3) = " <<IterativeParameter.at<double>(2,0)<<std::endl;
			std::cout << "IterativeParameter(4) = " <<IterativeParameter.at<double>(3,0)<<std::endl;
			std::cout << "IterativeParameter(5) = " <<IterativeParameter.at<double>(4,0)<<std::endl;
			return true;
		}

		IterativeParameter = IterativeParameter + ((J.t()*J)).inv()*J.t()*R;
		std::cout << "TotalErr = " <<TotalErr<<std::endl;
		bestErr = TotalErr;
	}
	return false;
}

void CLevenberGaussNewton::LMNewtonIterativeRunExam(std::vector<std::vector<cv::Mat>> &inputData,cv::Mat&outputData)
{

	int DataVecNum = inputData.size();
	for (int i=0;i<DataVecNum;i++)
	{
		(this->inputData[i]).assign(inputData[i].begin(),inputData[i].end());
	}
	this->outputData = outputData.clone();
	//=========================数据输入范例y = Aexp(Bx)
		//int num_params = 2;  
		//int total_data = 8;  
		//cv::Mat IterativeParameter = (cv::Mat_<double>(2,1)<<6,0.3);
		//cv::Mat inputData(total_data, 1, CV_64F);  
		//cv::Mat outputData(total_data, 1, CV_64F); 
		//for(int i=0; i < total_data; i++) {  
		//	inputData.at<double>(i,0) = i+1;  //load year  
		//}  
		//outputData.at<double>(0,0)= 8.3;  
		//outputData.at<double>(1,0)= 11.0;  
		//outputData.at<double>(2,0)= 14.7;  
		//outputData.at<double>(3,0)= 19.7;  
		//outputData.at<double>(4,0)= 26.7;  
		//outputData.at<double>(5,0)= 35.2;  
		//outputData.at<double>(6,0)= 44.4;  
		//outputData.at<double>(7,0)= 55.9;  
	//=========================数据输入范例y = Aexp(Bx)
		LMNewtonIterative(IterativeParameter,inputData,outputData);
}
void CLevenberGaussNewton::NewtonIterativeRunExam(cv::vector<cv::vector<cv::Mat>> &inputData,cv::Mat&outputData)
{
	int DataVecNum = inputData.size();
	for (int i=0;i<DataVecNum;i++)
	{
		(this->inputData[i]).assign(inputData[i].begin(),inputData[i].end());
	}
	this->outputData = outputData.clone();
	//=========================数据输入范例y = Aexp(Bx)
	//int num_params = 2;  
	//int total_data = 8;  
	//cv::Mat IterativeParameter = (cv::Mat_<double>(2,1)<<6,0.3);
	//cv::Mat inputData(total_data, 1, CV_64F);  
	//cv::Mat outputData(total_data, 1, CV_64F); 
	//for(int i=0; i < total_data; i++) {  
	//	inputData.at<double>(i,0) = i+1;  //load year  
	//}  
	//outputData.at<double>(0,0)= 8.3;  
	//outputData.at<double>(1,0)= 11.0;  
	//outputData.at<double>(2,0)= 14.7;  
	//outputData.at<double>(3,0)= 19.7;  
	//outputData.at<double>(4,0)= 26.7;  
	//outputData.at<double>(5,0)= 35.2;  
	//outputData.at<double>(6,0)= 44.4;  
	//outputData.at<double>(7,0)= 55.9;  
	//=========================数据输入范例y = Aexp(Bx)
	NewtonIterative(IterativeParameter,inputData,outputData);
}