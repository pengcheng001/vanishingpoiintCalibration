//#include <opencv2/opencv.hpp>
//const double DERIV_STEP = 1e-5;  
//double FunC(const cv::Mat &IterativeParameter,const double inputData)
//{
//	return IterativeParameter.at<double>(0,0)*exp(IterativeParameter.at<double>(1,0)*inputData);
//}
//double DerivativeSolving(cv::Mat IterativeParameter,const double inputData,const int n)
//{
//	cv::Mat TempParameter1 = IterativeParameter.clone();
//	cv::Mat TempParameter2 = IterativeParameter.clone();
//
//    TempParameter1.at<double>(n,0) = TempParameter1.at<double>(n,0) + DERIV_STEP;
//    TempParameter2.at<double>(n,0) = TempParameter2.at<double>(n,0) - DERIV_STEP;
//	double d = (FunC(TempParameter1,inputData) - FunC(TempParameter2,inputData))/(2*DERIV_STEP);
//	return ((FunC(TempParameter1,inputData) - FunC(TempParameter2,inputData))/(2*DERIV_STEP));
//}
//
//void NewtonIterative(cv::Mat &IterativeParameter,cv::Mat &inputData,cv::Mat& outputData)
//{
//	double bestErr=0.0;
//	int inputDataRowNum = inputData.rows;
//	int inputDataColNum = inputData.cols;
//	int ParameterNum = IterativeParameter.rows;
//	cv::Mat J(cv::Size(ParameterNum,inputDataRowNum),CV_64F);
//	cv::Mat R(cv::Size(1,inputDataRowNum),CV_64F);
//	for (int iter =0 ; iter < 100 ;iter++)
//	{
//		double TotalErr=0.0;
//	for (int i=0;i<inputDataRowNum;i++)
//	{
//		double tempInputData = inputData.at<double>(i,0);
//		double tempOutputData = outputData.at<double>(i,0);
//		R.at<double>(i,0) = tempOutputData - FunC(IterativeParameter,tempInputData);
//		TotalErr+=std::pow(R.at<double>(i,0),2);
//		for (int j=0;j<ParameterNum;j++)
//		{
//			J.at<double>(i,j) = DerivativeSolving(IterativeParameter ,tempInputData , j);
//		}
//	}
//	TotalErr/=inputDataRowNum;
//	if (abs(bestErr-TotalErr)<1e-5)		
//	{
//		std::cout << "TotalErr = " <<TotalErr<<std::endl;
//		std::cout << "IterativeParameter(1) = " <<IterativeParameter.at<double>(0,0)<<std::endl;
//		std::cout << "IterativeParameter(2) = " <<IterativeParameter.at<double>(1,0)<<std::endl;
//		break;
//	}
//
//	IterativeParameter = IterativeParameter + ((J.t()*J)).inv()*J.t()*R;
//	std::cout << "TotalErr = " <<TotalErr<<std::endl;
//	bestErr = TotalErr;
//	}
//}
//int main0()
//{
//	//参数数量，数据大小分配空间
//	int num_params = 2;  
//	int total_data = 8;  
//	cv::Mat IterativeParameter = (cv::Mat_<double>(2,1)<<6,0.3);
//	
//	
//	cv::Mat inputData(total_data, 1, CV_64F);  
//	cv::Mat outputData(total_data, 1, CV_64F);  
//
//	//load observation data  
//	for(int i=0; i < total_data; i++) {  
//		inputData.at<double>(i,0) = i+1;  //load year  
//	}  
//	//load America population  
//	outputData.at<double>(0,0)= 8.3;  
//	outputData.at<double>(1,0)= 11.0;  
//	outputData.at<double>(2,0)= 14.7;  
//	outputData.at<double>(3,0)= 19.7;  
//	outputData.at<double>(4,0)= 26.7;  
//	outputData.at<double>(5,0)= 35.2;  
//	outputData.at<double>(6,0)= 44.4;  
//	outputData.at<double>(7,0)= 55.9;  
//	NewtonIterative(IterativeParameter,inputData,outputData);
//}
//
