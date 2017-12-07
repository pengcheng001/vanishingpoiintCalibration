#include "headGather.h"
#include <math.h>

//设置虚拟相机参数，以及空间坐标点
void CVirVainshingPoint::SetRotateTranslationParameter(double alpha ,  double beta ,  double gamma ,  double T[][1]  )
{
	//初始化平移矩阵
	cv::Mat Rotate,Translation;
	cv::Mat TT(cv::Size(1,3),CV_64F,T);
	/*std::cout << TT <<std::endl;*/
	TT.copyTo(Translation);
	Pos_Translation.push_back(Translation);
	alpha = alpha * 3.1415926535 / 180 ;
	beta = beta * 3.1415926535 / 180 ; 
	gamma = gamma * 3.1415926535 /180 ; 
	double r1[3][3] = {{cos(gamma),-sin(gamma),0},{sin(gamma),cos(gamma),0},{0,0,1}};
	double r2[3][3] = {{cos(beta),0,sin(beta)},{0,1,0},{-sin(beta) , 0 , cos(beta)}};
	double r3[3][3] = {{1,0,0},{0,cos(alpha),-sin(alpha)},{0,sin(alpha),cos(alpha)}};
	cv::Mat r11(3,3,CV_64F,r1);
	cv::Mat r22(3,3,CV_64F,r2);
	cv::Mat r33(3,3,CV_64F,r3);
	Rotate = r11 * r22 * r33 ;
	/*std::cout << Rotate <<std::endl;*/
	Pos_Rotate.push_back(Rotate);
	
}
void CVirVainshingPoint::SetVirtualCameraParameter(const double f_set , const double u0_set ,const double  v0_set)
{
	
	assert((dx!=0)||(dy!=0));
	this->f_set = f_set;this->u0_set = u0_set;
	this->v0_set = v0_set;this ->fx_set = f_set/dx;
	this->fy_set = f_set/dy;
	CoA_matrix =( cv::Mat_<double>(3,3)<<fx_set , 0 , u0_set , 0 , fy_set , v0_set , 0 , 0 , 1);
}
void CVirVainshingPoint::WorldCoordinateToUV(std::vector<cv::Mat> &Point,std::vector<std::vector<cv::Mat>> &noBiasPoint_Pos)
{
	int RotateSize = Pos_Rotate.size();
	int PointSize = Point.size();

	for (int j=0;j<RotateSize;j++)
	{
	std::vector<cv::Mat> tempVector;
		for (int i=0;i<PointSize;i++)
	{
		
		//std::cout << CoA_matrix <<std::endl;
		//std::cout << Pos_Rotate[j] <<std::endl;
		//std::cout <<Point[i] <<std::endl;
		//std::cout <<Pos_Translation[j]<<std::endl;

		cv::Mat tempPoint = (CoA_matrix*Pos_Rotate[j]*Point[i] + Pos_Translation[j]);
		cv::Mat tempMat =(cv::Mat_<double>(2,1)<<(tempPoint.at<double>(0,0)/tempPoint.at<double>(2,0)),(tempPoint.at<double>(1,0)/tempPoint.at<double>(2,0)));
		std::cout <<"原来的点："<<tempMat<<std::endl;
		 tempVector.push_back(tempMat);
	}
	noBiasPoint_Pos.push_back(tempVector);
}
}
//获取畸变坐标值
void CVirVainshingPoint::UvnoBiasToDistortion(const std::vector<std::vector<cv::Mat>> noBiasPoint_Pos,std::vector<std::vector<cv::Mat>> &BiasPoint_Pos)
{
	if (noBiasPoint_Pos.size()==0)
	{
		std::cout << "没有输入任何角点，请检查后重试"<<std::endl;
		exit(1);
	}
	int PosNum = noBiasPoint_Pos.size();

	for (int i=0;i<PosNum;i++)
	{
		int PointNum = noBiasPoint_Pos[i].size();
		std::vector<cv::Mat> tempVector;
		for (int j=0;j<PointNum;j++)
		{
			double D_x = 0.0,D_x1=1.0, D_y = 0.0;
				//迭代法求解畸变偏值
			while (abs(D_x-D_x1)>0.0001)
			{
				D_x1 = D_x ;
				D_y = std::pow((dy/dx),2) *  (noBiasPoint_Pos[i][j].at<double>(1,0)-v0_set) * D_x1 /(noBiasPoint_Pos[i][j].at<double>(0,0)-u0_set);
				D_x = -distortion_parameter*std::pow(dx,2)*(std::pow( (noBiasPoint_Pos[i][j].at<double>(0,0)+D_x1-u0_set),2) +
					std::pow( (noBiasPoint_Pos[i][j].at<double>(1,0)+D_y-v0_set),2)) * (noBiasPoint_Pos[i][j].at<double>(0,0)-u0_set+D_x1);
				//std::cout <<"D_x"<<D_x <<std::endl;
				//std::cout <<"D_x1 "<<D_x1 <<std::endl;
			}
			D_y = std::pow((dy/dx),2) *  (noBiasPoint_Pos[i][j].at<double>(1,0)-v0_set) * D_x1 /(noBiasPoint_Pos[i][j].at<double>(0,0)-u0_set);
			cv::Mat tempPos = (cv::Mat_<double>(2,1)<<noBiasPoint_Pos[i][j].at<double>(0,0) + D_x,noBiasPoint_Pos[i][j].at<double>(1,0) + D_y);
			
			//std::cout << tempPos <<std::endl;
			tempVector.push_back(tempPos);
		}
		BiasPoint_Pos.push_back(tempVector);
	}
}


//求取图像坐标系下的交点。
void CVirVainshingPoint::GetCrossPoint(const std::vector<std::vector<cv::Mat>>Point_Pos , std::vector<std::vector<cv::Mat>>& CrossPoint )
{
	int Point_PosNum = Point_Pos.size();
	for (int i=0;i<Point_PosNum;i++)
	{
		std::vector<cv::Mat> tempMat;

		
		//std::cout <<  Point_Pos[i][1].at<double>(1,0) << "  "<< Point_Pos[i][0].at<double>(1,0) << Point_Pos[i][1].at<double>(0,0) <<"  "<<Point_Pos[i][0].at<double>(0,0)<<std::endl;
		double Line12_co[2] = {Point_Pos[i][1].at<double>(1,0) - Point_Pos[i][0].at<double>(1,0) ,( -Point_Pos[i][1].at<double>(0,0) + Point_Pos[i][0].at<double>(0,0))} ;
		double Line12_b = Point_Pos[i][0].at<double>(0,0) * ( Point_Pos[i][1].at<double>(1,0) - Point_Pos[i][0].at<double>(1,0)) - 
			Point_Pos[i][0].at<double>(1,0) * (Point_Pos[i][1].at<double>(0,0)-Point_Pos[i][0].at<double>(0,0)) ;

		double Line34_co[2] ={ Point_Pos[i][3].at<double>(1,0) - Point_Pos[i][2].at<double>(1,0), - Point_Pos[i][3].at<double>(0,0) + Point_Pos[i][2].at<double>(0,0)} ;
		double Line34_b = Point_Pos[i][2].at<double>(0,0) * ( Point_Pos[i][3].at<double>(1,0) - Point_Pos[i][2].at<double>(1,0)) - 
			Point_Pos[i][2].at<double>(1,0) * (Point_Pos[i][3].at<double>(0,0)-Point_Pos[i][2].at<double>(0,0)) ;
	
		double Line14_co[2] ={ Point_Pos[i][3].at<double>(1,0) - Point_Pos[i][0].at<double>(1,0) ,- Point_Pos[i][3].at<double>(0,0) + Point_Pos[i][0].at<double>(0,0)} ;
		double Line14_b = Point_Pos[i][0].at<double>(0,0) * ( Point_Pos[i][3].at<double>(1,0) - Point_Pos[i][0].at<double>(1,0)) - 
			Point_Pos[i][0].at<double>(1,0) * (Point_Pos[i][3].at<double>(0,0)-Point_Pos[i][0].at<double>(0,0)) ;

		double Line23_co[2] ={ Point_Pos[i][2].at<double>(1,0) - Point_Pos[i][1].at<double>(1,0), - Point_Pos[i][2].at<double>(0,0) + Point_Pos[i][1].at<double>(0,0) };
		double Line23_b = Point_Pos[i][1].at<double>(0,0) * ( Point_Pos[i][2].at<double>(1,0) - Point_Pos[i][1].at<double>(1,0)) - 
			Point_Pos[i][1].at<double>(1,0) * (Point_Pos[i][2].at<double>(0,0)-Point_Pos[i][1].at<double>(0,0)) ;
		
		cv::Mat Across,Bcross;
		cv::Mat A_co = (cv::Mat_<double>(2,2)<<Line12_co[0],Line12_co[1],Line34_co[0],Line34_co[1]);
		cv::Mat A_b = (cv::Mat_<double>(2,1)<<Line12_b,Line34_b);
		cv::Mat B_co = (cv::Mat_<double>(2,2)<<Line14_co[0],Line14_co[1],Line23_co[0],Line23_co[1]);
		cv::Mat B_b = (cv::Mat_<double>(2,1)<<Line14_b,Line23_b);
		Across = A_co.inv() * A_b; 
		//std::cout << Across<<std::endl;
		Bcross = B_co.inv() * B_b;
		//std::cout << Bcross<<std::endl;
		tempMat.push_back(Across);
		tempMat.push_back(Bcross);

	CrossPoint.push_back(tempMat);
}
}
//由畸变点反推回相机坐标系下的坐标
void CVirVainshingPoint::DistortionToUvnoBias(const std::vector<std::vector<cv::Mat>> CrossPoint,
	std::vector<std::vector<cv::Mat>> &BiasTonoBiasCrossPoint,cv::Mat& outParameter)
{
	if (CrossPoint.size()==0)
	{
		std::cout << "输入原始相交点为0，请检查后重试" <<std::endl;
	}
	BiasTonoBiasCrossPoint =CrossPoint;
	cv::Mat tempOutput(cv::Size(1,10),CV_64F);
	//std::vector<std::vector<cv::Mat>> tempCrossPoint = CrossPoint ;
	CLevenberGaussNewton *B = new CLevenberGaussNewton(outParameter,100);
	double tempEnd = 0.0;
	//do 
	//{
	//	tempEnd = B->IterativeParameter.at<double>(0,0);
	//	BiasTonoBiasCrossPoint.clear();
	//	int PosNum = CrossPoint.size();
	//	for (int i=0;i<PosNum;i++)
	//	{
	//		int PointNum = CrossPoint[i].size();
	//		std::vector<cv::Mat> tempVec;
	//		for (int j=0;j<PointNum;j++)
	//		{
	//			double BiasU = tempCrossPoint[i][j].at<double>(0,0) ;
	//			double BiasV = tempCrossPoint[i][j].at<double>(1,0) ;
	//			double k = B->IterativeParameter.at<double>(4,0);
	//			double U_0 = B->IterativeParameter.at<double>(0,0);
	//			double V_0 = B->IterativeParameter.at<double>(1,0);
	//			tempCrossPoint[i][j].at<double>(0,0) = BiasU + k*(std::pow((BiasU-U_0),2)+std::pow((BiasV-V_0),2))*dx*dx*(BiasU-U_0);
	//			tempCrossPoint[i][j].at<double>(1,0) = BiasV + k*(std::pow((BiasU-U_0),2)+std::pow((BiasV-V_0),2))*dy*dy*(BiasV-V_0);
	//			double tempX = (tempCrossPoint[i][j].at<double>(0,0) - U_0)/ B->IterativeParameter.at<double>(2,0);
	//			double tempY = (tempCrossPoint[i][j].at<double>(1,0) - V_0)/ B->IterativeParameter.at<double>(3,0);
	//			double tempZ = 1.0;
	//			cv::Mat tempMat = (cv::Mat_<double>(3,1)<<tempX,tempY,tempZ);
	//			tempVec.push_back(tempMat);
	//		}
	//		BiasTonoBiasCrossPoint.push_back(tempVec);

	//	}
	//}while ((B->LMNewtonIterative(B->IterativeParameter,BiasTonoBiasCrossPoint,tempOutput))&&abs(B->IterativeParameter.at<double>(0,0)-tempEnd)<1e-5);
	while (!((B->LMNewtonIterative(B->IterativeParameter,BiasTonoBiasCrossPoint,tempOutput))&&abs(B->IterativeParameter.at<double>(0,0)-tempEnd)<1e-5))
	{
		tempEnd = B->IterativeParameter.at<double>(0,0);
	}
	delete B;
}
void CVirVainshingPoint::CalculateCoamatrix(std::vector<std::vector<cv::Mat>> CrossPoint ,cv::Mat &CalculateCoa_matrix)
{
	int CrossPoint_Pos = CrossPoint.size();
	if (CrossPoint_Pos==0)
	{
		std::cout << "调整拍摄位置，无相交点" << std::endl;
		exit(1);
	}
	std::vector<cv::Mat> AnsTemp;
	double CoMatrixU[3]={0},CoMatrixV[3]={0},CoMatrixVV[3]={0},CoMatrixUV[3]={0};

	for (int g=0;g<CrossPoint_Pos-3;g++)
	{
		int k=g;
		for (int i=0;i<3;i++)
		{
		int j = k + 1 ;
		//                                          F  A                        U 
		CoMatrixU[i] = CrossPoint[j][0].at<double>(0,0)+CrossPoint[j][1].at<double>(0,0) - CrossPoint[k][0].at<double>(0,0) - CrossPoint[k][1].at<double>(0,0);
		CoMatrixV[i] = CrossPoint[j][0].at<double>(1,0)+CrossPoint[j][1].at<double>(1,0) - CrossPoint[k][0].at<double>(1,0) - CrossPoint[k][1].at<double>(1,0);
		CoMatrixVV[i] = CrossPoint[k][0].at<double>(1,0) * CrossPoint[k][1].at<double>(1,0) - CrossPoint[j][0].at<double>(1,0) * CrossPoint[j][1].at<double>(1,0);
		CoMatrixUV[i] =  CrossPoint[j][0].at<double>(0,0) * CrossPoint[j][1].at<double>(0,0) - CrossPoint[k][0].at<double>(0,0) * CrossPoint[k][1].at<double>(0,0);
		k++;
				}
	cv::Mat CoMatrix = (cv::Mat_<double>(3,3)<< CoMatrixU[0],CoMatrixV[0],CoMatrixVV[0],
																	CoMatrixU[1],CoMatrixV[1],CoMatrixVV[1],
																	CoMatrixU[2],CoMatrixV[2],CoMatrixVV[2]);
	cv::Mat  bMatrix = (cv::Mat_<double>(3,1)<<CoMatrixUV[0],CoMatrixUV[1],CoMatrixUV[2]);
	cv::Mat Ans = CoMatrix.inv() * bMatrix ;//里面存放这xyz需要进一步变换到内参
	AnsTemp.push_back(Ans);
	}
	//========================debug
	
	double SumAns_0 = 0.0,SumAns_1 = 0.0,SumAns_2 = 0.0;
	for (int i=0;i<CrossPoint_Pos-3;i++)
	{
		SumAns_0 += AnsTemp[i].at<double>(0,0);
		SumAns_1 += AnsTemp[i].at<double>(1,0);
		SumAns_2 += AnsTemp[i].at<double>(2,0);
	}
	SumAns_1 = SumAns_1/(CrossPoint_Pos-3);
	SumAns_2 = SumAns_2/(CrossPoint_Pos-3);
	double u0_calculate=SumAns_0/(CrossPoint_Pos-3);
	double v0_calculate = SumAns_1 / SumAns_2 ;
	double z =  SumAns_2;
	double fx[50]={0.0},fy[50]={0.0};
	double fx_calculate = 0.0,fy_calculate = 0.0;
	for (int i=0;i<CrossPoint_Pos;i++)
	{
		//std::cout <<CrossPoint[i][0].at<double>(0,0)<<"  "<< CrossPoint[i][1].at<double>(0,0);
		fx[i] = -(u0_calculate - CrossPoint[i][0].at<double>(0,0)) * (u0_calculate - CrossPoint[i][1].at<double>(0,0)) -
			z * (v0_calculate - CrossPoint[i][0].at<double>(1,0)) *(v0_calculate - CrossPoint[i][1].at<double>(1,0));
		fx[i] = std::sqrt(fx[i]);
		fx_calculate += fx[i];
	}
	fx_calculate = fx_calculate/CrossPoint_Pos;
	fy_calculate = fx_calculate / z;
	CalculateCoa_matrix = (cv::Mat_<double>(3,3)<<fx_calculate,0,u0_calculate,0,fy_calculate,v0_calculate,0,0,1);
}
void CVirVainshingPoint::CrunVirVanishingNobiasPoint()
{
	Set_dxdyParameter(0.009,0.009,0.0002);
	//设置相机内参矩阵
	SetVirtualCameraParameter(30,340,270);
	//设置四种不同的外参数矩阵
	double T1[3][1]={3000,2000,1000};
	double T2[3][1]={2000,1000,1000};
	double T3[3][1]={1000,-1500,4000};
	double T4[3][1]={1000,1500,3000};
	//
	//SetRotateTranslationParameter(10,10,20,T1);
	//SetRotateTranslationParameter(15,-30,10,T2);
	SetRotateTranslationParameter(40,20,30,T3);
	SetRotateTranslationParameter(10,20,10,T4);
	SetRotateTranslationParameter(15,-30,18,T2);
	//SetRotateTranslationParameter(40,20,35,T1);
	SetRotateTranslationParameter(12,20,10,T1);
	SetRotateTranslationParameter(15,-35,10,T2);
	SetRotateTranslationParameter(40,20,35,T3);
	SetRotateTranslationParameter(16,20,10,T4);
	//虚拟使用的角点，世界坐标系下坐标
	const cv::Mat PointA = (cv::Mat_<double>(3,1)<<0 , 0.1 , 0);
	const cv::Mat PointB = (cv::Mat_<double>(3,1)<<0.1 , 30.1 , 0);
	const cv::Mat PointC = (cv::Mat_<double>(3,1)<<40 , 30.1 , 0);
	const cv::Mat PointD = (cv::Mat_<double>(3,1)<<40 , 0.1 , 0);
	PointTogether.push_back(PointA);
	PointTogether.push_back(PointB);
	PointTogether.push_back(PointC);
	PointTogether.push_back(PointD);
	//坐标系坐标转换
	WorldCoordinateToUV(PointTogether , noBiasPoint_Pos);//std::vector<std::vector<cv::Mat>> &noBiasPoint_Pos,noBiasPoint_Pos放着四组图四组点
	//求个幅图的消失点
	std::vector<std::vector<cv::Mat>> CrossPoint;
	GetCrossPoint(noBiasPoint_Pos , CrossPoint);

	cv::Mat CalculateCoa_matrix;
	CalculateCoamatrix(CrossPoint,CalculateCoa_matrix);
	std::cout << "虚拟设置的参数:" <<std::endl;
	std::cout << CoA_matrix<<std::endl;
	std::cout << "无畸变下消失点计算的参数:" <<std::endl;
	std::cout << CalculateCoa_matrix <<std::endl;
///////////////////////////////////debug
	//int RotateSize = Pos_Rotate.size();
	//int PointSize = PointTogether.size();

	//for (int j=0;j<RotateSize;j++)
	//{
	//	for (int i=0;i<PointSize;i++)
	//	{


	//		cv::Mat tempPoint = (CalculateCoa_matrix*Pos_Rotate[j]*PointTogether[i] + Pos_Translation[j]);
	//		cv::Mat tempMat =(cv::Mat_<double>(2,1)<<(tempPoint.at<double>(0,0)/tempPoint.at<double>(2,0)),(tempPoint.at<double>(1,0)/tempPoint.at<double>(2,0)));
	//		std::cout <<"重计算的点："<<tempMat<<std::endl;
	//	}
	//}

}

void CVirVainshingPoint::CrunVirVanishingbiasPoint()
{
	Set_dxdyParameter(0.009,0.009,0.0002);
	//设置相机内参矩阵
	SetVirtualCameraParameter(30,340,270);
	//设置四种不同的外参数矩阵
	double T1[3][1]={3000,2000,1000};
	double T2[3][1]={2000,1000,1000};
	double T3[3][1]={1000,-1500,4000};
	double T4[3][1]={1000,1500,3000};

	SetRotateTranslationParameter(10,10,20,T1);
	SetRotateTranslationParameter(15,-30,10,T2);
	SetRotateTranslationParameter(40,20,30,T3);
	SetRotateTranslationParameter(10,20,10,T4);

	//虚拟使用的角点，世界坐标系下坐标
	const cv::Mat PointA = (cv::Mat_<double>(3,1)<<40 , 0 , 0);
	const cv::Mat PointB = (cv::Mat_<double>(3,1)<<40 , 30.001 , 0);
	const cv::Mat PointC = (cv::Mat_<double>(3,1)<<80 , 30 , 0);
	const cv::Mat PointD = (cv::Mat_<double>(3,1)<<80 , 0 , 0);
	PointTogether.push_back(PointA);
	PointTogether.push_back(PointB);
	PointTogether.push_back(PointC);
	PointTogether.push_back(PointD);
	//坐标系坐标转换
	WorldCoordinateToUV(PointTogether , noBiasPoint_Pos);//std::vector<std::vector<cv::Mat>> &noBiasPoint_Pos,noBiasPoint_Pos放着四组图四组点
	//畸变点转换
	UvnoBiasToDistortion(noBiasPoint_Pos,BiasPoint_Pos);
	
	std::vector<std::vector<cv::Mat>> CrossPoint;
	GetCrossPoint(BiasPoint_Pos , CrossPoint);
	CalculateCoamatrix(CrossPoint,CalculateCoa_matrix);
	std::cout << "虚拟设置的参数:" <<std::endl;
	std::cout << CoA_matrix<<std::endl;
	std::cout << "畸变情况下消失点计算的参数:" <<std::endl;
	std::cout << CalculateCoa_matrix <<std::endl;
	std::cout <<"===============开始进行迭代求解最佳逼近参数================"<<std::endl;
	
	cv::Mat outParameter = (cv::Mat_<double>(5,1)<<CalculateCoa_matrix.at<double>(0,2),CalculateCoa_matrix.at<double>(1,2),
		CalculateCoa_matrix.at<double>(0,0),CalculateCoa_matrix.at<double>(1,1),0);
	DistortionToUvnoBias(CrossPoint,BiasTonoBiasCrossPoint,outParameter);

	cv::Mat DebugParameter = (cv::Mat_<double>(3,3)<<outParameter.at<double>(2,0),0,outParameter.at<double>(0,0),0,outParameter.at<double>(3,0),outParameter.at<double>(1,0),0,0,1);
	int RotateSize = Pos_Rotate.size();
	int PointSize = PointTogether.size();

	for (int j=0;j<RotateSize;j++)
	{
		for (int i=0;i<PointSize;i++)
		{


			cv::Mat tempPoint = (DebugParameter*Pos_Rotate[j]*PointTogether[i] + Pos_Translation[j]);
			cv::Mat tempMat =(cv::Mat_<double>(2,1)<<(tempPoint.at<double>(0,0)/tempPoint.at<double>(2,0)),(tempPoint.at<double>(1,0)/tempPoint.at<double>(2,0)));
			std::cout <<"重计算的点："<<tempMat<<std::endl;
		}
	}
}