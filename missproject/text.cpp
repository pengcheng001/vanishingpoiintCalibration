#include "headGather.h"
#include <time.h>
#include <fstream>
using namespace cv;
using namespace std;
int main()
{
	//=============================================
	//ֱ���޵������к������۲���������㣬�Ƚ��������
	//dx=0.009 dy=0.009
	//����LM-Gauss���ֲ�����
	//��ʧ����������ڲ�
	//���ۣ��㷨Ч�ʽϿ죬�л�������£�������3�������ڵ㼯������������6�����ص㡣
	//�޻�������£����������������ڣ��������3�����ص㡣
	//=============================================


	//std::vector<cv::Point2f> image_points_buf;  /* ����ÿ��ͼ���ϼ�⵽�Ľǵ� */  
	//std::vector<std::vector<cv::Point2f>> image_points_seq; /* �����⵽�����нǵ� */  
	// ifstream fin("calibdata.txt"); /* �궨����ͼ���ļ���·�� */  
	// Size board_size = Size(4,6);    /* �궨����ÿ�С��еĽǵ��� */
	// int  image_count = 0;
	// std::string filename;
	//while (getline(fin,filename))  
	//{  
	//	image_count++;        
	//	// ���ڹ۲�������  
	//	cout<<"image_count = "<<image_count<<endl;          

	//	Mat imageInput=imread(filename);  
	//	if (image_count == 1)  //�����һ��ͼƬʱ��ȡͼ������Ϣ  
	//	{    
	//		cout<<"image_size.width = "<<imageInput.cols<<endl;  
	//		cout<<"image_size.height = "<<imageInput.rows<<endl;  
	//	}  

	//	/* ��ȡ�ǵ� */  
	//	if (0 == findChessboardCorners(imageInput,board_size,image_points_buf))  
	//	{             
	//		cout<<"can not find chessboard corners!\n"; //�Ҳ����ǵ�  
	//		exit(1);  
	//	}   
	//	else   
	//	{  
	//		int temp1 = (board_size.height/2-1)*board_size.width + board_size.width/2-1;
	//		Mat view_gray;  
	//		cvtColor(imageInput,view_gray,CV_RGB2GRAY);  
	//		/* �����ؾ�ȷ�� */  
	//		find4QuadCornerSubpix(view_gray,image_points_buf,Size(5,5)); //�Դ���ȡ�Ľǵ���о�ȷ��  
	//		//cornerSubPix(view_gray,image_points_buf,Size(5,5),Size(-1,-1),TermCriteria(CV_TERMCRIT_EPS+CV_TERMCRIT_ITER,30,0.1));  
	//		image_points_seq.push_back(image_points_buf);  //���������ؽǵ�  image_points_buf  std::vector<std::vector<cv::Point2f>>
	//		/* ��ͼ������ʾ�ǵ�λ�� */  
	//		//drawChessboardCorners(view_gray,board_size,image_points_buf,false); //������ͼƬ�б�ǽǵ�  
	//		//cv::circle(view_gray,image_points_buf[temp1+board_size.width],10,cv::Scalar(255,255,255));
	//		//cv::circle(view_gray,image_points_buf[temp1+board_size.width+1],15,cv::Scalar(255,255,255));
	//		//cv::circle(view_gray,image_points_buf[temp1+1],20,cv::Scalar(255,255,255));
	//		//cv::circle(view_gray,image_points_buf[temp1],25,cv::Scalar(255,255,255));
	//		//imshow("Camera Calibration",view_gray);//��ʾͼƬ  
	//		//cv::waitKey(2000);//��ͣ0.5S         
	//	}  
	//}  
	//std::vector<std::vector<cv::Mat>> inputPointData;
	// //���������ؽǵ�  image_points_buf  std::vector<std::vector<cv::Point2f>>
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
	//std::cout << "�����������ʧ�����Ĳ���:" <<std::endl;
	//std::cout << CalculateCoa_matrix <<std::endl;

	//std::cout <<"===============��ʼ���е��������ѱƽ�����================"<<std::endl;
	//
	//cv::Mat outParameter = (cv::Mat_<double>(5,1)<<CalculateCoa_matrix.at<double>(0,2),CalculateCoa_matrix.at<double>(1,2),
	//	CalculateCoa_matrix.at<double>(0,0),CalculateCoa_matrix.at<double>(1,1),0);
	//std::vector<std::vector<cv::Mat>> BiasTonoBiasCrossPoint;
	//A.DistortionToUvnoBias(CrossPoint,BiasTonoBiasCrossPoint,outParameter);

	CVirVainshingPoint A;
	A.CrunVirVanishingbiasPoint();

	return 0;
}