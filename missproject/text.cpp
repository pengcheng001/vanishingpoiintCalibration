#include "headGather.h"
#include <time.h>
#include <fstream>
using namespace cv;
using namespace std;
int main()
{
	//=============================================
	//直接无调试运行函数：观察输入输出点，比较像素误差
	//dx=0.009 dy=0.009
	//采用LM-Gauss求解局部最优
	//消失点理论求解内参
	//结论：算法效率较快，有畸变情况下，精度在3个像素内点集噪声，误差最大6个像素点。
	//无畸变情况下，精度在三个像素内，精度最大3个像素点。
	//=============================================


	//std::vector<cv::Point2f> image_points_buf;  /* 缓存每幅图像上检测到的角点 */  
	//std::vector<std::vector<cv::Point2f>> image_points_seq; /* 保存检测到的所有角点 */  
	// ifstream fin("calibdata.txt"); /* 标定所用图像文件的路径 */  
	// Size board_size = Size(4,6);    /* 标定板上每行、列的角点数 */
	// int  image_count = 0;
	// std::string filename;
	//while (getline(fin,filename))  
	//{  
	//	image_count++;        
	//	// 用于观察检验输出  
	//	cout<<"image_count = "<<image_count<<endl;          

	//	Mat imageInput=imread(filename);  
	//	if (image_count == 1)  //读入第一张图片时获取图像宽高信息  
	//	{    
	//		cout<<"image_size.width = "<<imageInput.cols<<endl;  
	//		cout<<"image_size.height = "<<imageInput.rows<<endl;  
	//	}  

	//	/* 提取角点 */  
	//	if (0 == findChessboardCorners(imageInput,board_size,image_points_buf))  
	//	{             
	//		cout<<"can not find chessboard corners!\n"; //找不到角点  
	//		exit(1);  
	//	}   
	//	else   
	//	{  
	//		int temp1 = (board_size.height/2-1)*board_size.width + board_size.width/2-1;
	//		Mat view_gray;  
	//		cvtColor(imageInput,view_gray,CV_RGB2GRAY);  
	//		/* 亚像素精确化 */  
	//		find4QuadCornerSubpix(view_gray,image_points_buf,Size(5,5)); //对粗提取的角点进行精确化  
	//		//cornerSubPix(view_gray,image_points_buf,Size(5,5),Size(-1,-1),TermCriteria(CV_TERMCRIT_EPS+CV_TERMCRIT_ITER,30,0.1));  
	//		image_points_seq.push_back(image_points_buf);  //保存亚像素角点  image_points_buf  std::vector<std::vector<cv::Point2f>>
	//		/* 在图像上显示角点位置 */  
	//		//drawChessboardCorners(view_gray,board_size,image_points_buf,false); //用于在图片中标记角点  
	//		//cv::circle(view_gray,image_points_buf[temp1+board_size.width],10,cv::Scalar(255,255,255));
	//		//cv::circle(view_gray,image_points_buf[temp1+board_size.width+1],15,cv::Scalar(255,255,255));
	//		//cv::circle(view_gray,image_points_buf[temp1+1],20,cv::Scalar(255,255,255));
	//		//cv::circle(view_gray,image_points_buf[temp1],25,cv::Scalar(255,255,255));
	//		//imshow("Camera Calibration",view_gray);//显示图片  
	//		//cv::waitKey(2000);//暂停0.5S         
	//	}  
	//}  
	//std::vector<std::vector<cv::Mat>> inputPointData;
	// //保存亚像素角点  image_points_buf  std::vector<std::vector<cv::Point2f>>
	//for (int i=0;i<image_points_seq.size();i++)
	//{
	//		std::vector<cv::Mat> tempVec;
	//		int temp1 = (board_size.height/2-1)*board_size.width + board_size.width/2-1;
	//		tempVec.push_back ((cv::Mat_<double>(2,1)<<image_points_seq[i][temp1+board_size.width].x,image_points_seq[i][temp1+board_size.width].y));
	//		tempVec.push_back ((cv::Mat_<double>(2,1)<<image_points_seq[i][temp1+board_size.width+1].x,image_points_seq[i][temp1+board_size.width+1].y));	
	//		tempVec.push_back ( (cv::Mat_<double>(2,1)<<image_points_seq[i][temp1+1].x,image_points_seq[i][temp1+1].y));
	//		tempVec.push_back ( (cv::Mat_<double>(2,1)<<image_points_seq[i][temp1].x,image_points_seq[i][temp1].y));
	//		inputPointData.push_back(tempVec);
	//		if (inputPointData.size()==1)
	//		{
	//			std::cout <<tempVec[0]<<std::endl;
	//			std::cout <<tempVec[1]<<std::endl;
	//			std::cout <<tempVec[2]<<std::endl;
	//			std::cout <<tempVec[3]<<std::endl;
	//		}
	//}
	//CVirVainshingPoint A;
	//std::vector<std::vector<cv::Mat>> CrossPoint;
	//A.GetCrossPoint(inputPointData,CrossPoint);
	//cv::Mat CalculateCoa_matrix;
	//A.CalculateCoamatrix(CrossPoint,CalculateCoa_matrix);
	//std::cout << "畸变情况下消失点计算的参数:" <<std::endl;
	//std::cout << CalculateCoa_matrix <<std::endl;

	//std::cout <<"===============开始进行迭代求解最佳逼近参数================"<<std::endl;
	//
	//cv::Mat outParameter = (cv::Mat_<double>(5,1)<<CalculateCoa_matrix.at<double>(0,2),CalculateCoa_matrix.at<double>(1,2),
	//	CalculateCoa_matrix.at<double>(0,0),CalculateCoa_matrix.at<double>(1,1),0);
	//std::vector<std::vector<cv::Mat>> BiasTonoBiasCrossPoint;
	//A.DistortionToUvnoBias(CrossPoint,BiasTonoBiasCrossPoint,outParameter);

	CVirVainshingPoint A;
	A.CrunVirVanishingbiasPoint();

	return 0;
}