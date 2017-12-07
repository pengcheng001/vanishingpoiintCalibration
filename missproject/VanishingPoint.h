#ifndef _VANISHINGPOINT_H
#define _VANISHINGPOINT_H

#include "headGather.h"

class CVirVainshingPoint
{
public:
	CVirVainshingPoint(){};
	~CVirVainshingPoint(){};

public:
	void CrunVirVanishingNobiasPoint();
	void CrunVirVanishingbiasPoint();

	void SetVirtualCameraParameter(const double f_set ,const double u0_set ,const double  v0_set);
	void SetRotateTranslationParameter( double alpha ,  double beta ,  double gamma ,  double T[][1] );
	void Set_dxdyParameter(const double dx ,const double dy,const double distortion_parameter){
		this->dx = dx ;
		this->dy = dy ;
		this->distortion_parameter = distortion_parameter;
	};
	void	WorldCoordinateToUV(std::vector<cv::Mat> &Point,std::vector<std::vector<cv::Mat>> &noBiasPoint_Pos);
	void UvnoBiasToDistortion(const std::vector<std::vector<cv::Mat>> noBiasPoint_Pos,std::vector<std::vector<cv::Mat>> &BiasPoint_Pos);
	void DistortionToUvnoBias(const std::vector<std::vector<cv::Mat>> CrossPoint,
		std::vector<std::vector<cv::Mat>> &BiasTonoBiasCrossPoint,cv::Mat& outParameter);
	void GetCrossPoint(const std::vector<std::vector<cv::Mat>>Point_Pos ,  std::vector<std::vector<cv::Mat>>& CrossPoint );
	void CalculateCoamatrix(std::vector<std::vector<cv::Mat>> CrossPoint ,cv::Mat &CalculateCoa_matrix);
	//==========�����������==========
private:
	double f_set , u0_set , v0_set;
	double fx_set , fy_set;
	double dx,dy;
	double distortion_parameter;
	cv::Mat CoA_matrix,CalculateCoa_matrix,outParameter;
	std::vector<cv::Mat> Pos_Translation;
	std::vector<cv::Mat> Pos_Rotate;
	//=========��ʼ�㼯�ϣ����޻���======
private:
	std::vector<cv::Mat> PointTogether;
	std::vector<std::vector<cv::Mat>> noBiasPoint_Pos; //�����޻����ͼ��
	std::vector<std::vector<cv::Mat>> BiasTonoBiasCrossPoint;//��ʧ�㰴�������ڲ����·����޻����
	std::vector<std::vector<cv::Mat>> BiasPoint_Pos;//�Ƴ������ͼ��
};



#endif //VanishingPoint.h
