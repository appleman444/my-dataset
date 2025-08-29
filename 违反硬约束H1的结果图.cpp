//#include <vtkAutoInit.h>
//VTK_MODULE_INIT(vtkRenderingOpenGL2);
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingFreeType);
//#include "vtkDICOMImageReader.h"//DCM医学文件读取类
//#include "vtkPolyDataWriter.h"//保存为.vtk图像类
//#include "vtkPolyDataReader.h"//读取.vtk图像类
//#include<vtkOBBTree.h>
//#include<vtkBooleanOperationPolyDataFilter.h>
//#include<vtkTubeFilter.h>
//#include<vtkLandmarkTransform.h>
//#include <vtkSmartPointer.h>
//#include<vtkTransform.h>
//#include<vtkMatrix4x4.h>
//#include <vtkJPEGReader.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkInteractorStyleImage.h>
//#include <vtkRenderer.h>
//#include <vtkStringArray.h>
//#include<vtkPlaneSource.h>
//#include <vtkRenderWindow.h>
//#include <vtkContourFilter.h>
//#include <vtkPolyDataNormals.h>
//#include <vtkStripper.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include<vtkVertexGlyphFilter.h>
//#include <vtkTextActor.h>
//#include<vtkTriangleFilter.h>
//#include <vtkTextProperty.h>
//#include<vtkAxesActor.h>
//#include<vtkIterativeClosestPointTransform.h>
//#include <vtkCamera.h>
//#include <vtkSphereSource.h>
//#include <vtkSTLWriter.h>
//#include<vtkRegularPolygonSource.h>
//#include <vtkSTLReader.h>
//#include <vtkSmoothPolyDataFilter.h>
//#include <vtkActor.h>
//#include <vtkMassProperties.h>
//#include <vtkFloatArray.h>
//#include <vtkLookupTable.h>
//#include <vtkPointData.h>
//#include<vtkTransformPolyDataFilter.h>
//#include <vtkPolyData.h>
//#include<vtkProperty.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderer.h>
//#include<vtkSelectEnclosedPoints.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkSmartPointer.h>
//#include <vtkSTLReader.h>
//#include<vtkPlane.h>
//#include<vtkNew.h>
//#include<vtkPoints.h>
//#include<vtkMassProperties.h>
//#include<vtkLine.h>
//#include<vtkMath.h>
//#include<vector>
//#include<vtkLineSource.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include<vtkParametricEllipsoid.h>
//#include<vtkParametricFunctionSource.h>
//#include<random>
//#include<vector>
//#include<iostream>
//#include<cmath>
//#include<ctime>
//#include <cstdlib>
//#include <bitset>
//#include<iomanip>
//#include <fstream>  // 用于文件操作
//#include <vtkBoundingBox.h>
//#include <vtkCubeSource.h>
//#include <vtkSmartPointer.h>
//#include <vtkPoints.h>
//#include <vtkCellArray.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkLinearTransform.h>
//#include <vtkTransformPolyDataFilter.h>
//#include <vtkCylinderSource.h>
//#include <vtkSmartPointer.h>
//#include <vtkCylinderSource.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkTransform.h>
//#include <vtkTransformPolyDataFilter.h>
//#include <vtkLine.h>
//#include <vtkMath.h>
//using namespace std;
//const double PI = 3.141592653589793;//定义一个不可改变的常量值PI
//const int Po_Size = 50;//种群规模
//const int Ev_Algebra = 18000;//进化代数
//const double Ov_Probability = 0.850; //交叉概率,交叉概率用于判断两两个体是否需要交叉
//const double Va_Probability = 0.10;//变异概率,变异概率用于判断任一个体是否需要变异
//vector<double>line_to_center;
//double tube_radius = 100;
//double Ellipsoid_long = 20;
//double Ellipsoid_short = 15;
//
////point 的actor 的构建的形式
//vtkSmartPointer<vtkActor> get_point_actor(vtkSmartPointer<vtkPoints> points) {
//	vtkSmartPointer<vtkPolyData> pointsPolydata =
//		vtkSmartPointer<vtkPolyData>::New();
//	pointsPolydata->SetPoints(points);
//	vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	vertexGlyphFilter->AddInputData(pointsPolydata);
//	vertexGlyphFilter->Update();
//	vtkSmartPointer<vtkPolyDataMapper> pointsMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	pointsMapper->SetInputConnection(vertexGlyphFilter->GetOutputPort());
//	vtkSmartPointer<vtkActor> pointsActor =
//		vtkSmartPointer<vtkActor>::New();
//	pointsActor->SetMapper(pointsMapper);
//	pointsActor->GetProperty()->SetPointSize(0.1);//定义点的尺寸大小，这样点才能在画布上显示出来
//	return pointsActor;
//}
//
//
////创建读取STL文件的函数
//vtkSmartPointer<vtkPolyData> Reader_STL(const char* path_name, const char* file_name) {
//	if (path_name == NULL || file_name == NULL) {
//		cout << "The path error!!!!" << endl;
//		return NULL;
//	}
//	char total_path[100] = { 0 };
//	strcat_s(total_path, path_name);
//	strcat_s(total_path, file_name);
//	vtkSmartPointer<vtkSTLReader> stl_Reader = vtkSmartPointer<vtkSTLReader>::New();
//	stl_Reader->SetFileName(total_path);
//	stl_Reader->Update();
//	vtkSmartPointer<vtkPolyData>Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//	Poly_Data->DeepCopy(stl_Reader->GetOutput());
//	return Poly_Data;
//}
////创建读取VTK 文件的函数
//vtkSmartPointer<vtkPolyData> Reader_VTK(const char* path_name, const char* file_name) {
//	if (path_name == NULL || file_name == NULL) {
//		cout << "The path error!!!!" << endl;
//		return NULL;
//	}
//	char total_path[100] = { 0 };
//	strcat_s(total_path, path_name);
//	strcat_s(total_path, file_name);
//	vtkSmartPointer<vtkPolyDataReader> vtk_Reader = vtkSmartPointer<vtkPolyDataReader>::New();
//	vtk_Reader->SetFileName(total_path);
//	vtk_Reader->Update();
//	vtkSmartPointer<vtkPolyData>Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//	Poly_Data->DeepCopy(vtk_Reader->GetOutput());
//	return Poly_Data;
//}
////创建对应的VTK 文件的actor 的形式
//vtkSmartPointer<vtkActor> Reader_VTK_Actor(const char* path_name, const char* file_name) {
//	if (path_name == NULL || file_name == NULL) {
//		cout << "The path error!!!!" << endl;
//		return NULL;
//	}
//	char total_path[100] = { 0 };
//	strcat_s(total_path, path_name);
//	strcat_s(total_path, file_name);
//	vtkSmartPointer<vtkPolyDataReader> vtk_Reader = vtkSmartPointer<vtkPolyDataReader>::New();
//	vtk_Reader->SetFileName(total_path);
//	vtk_Reader->Update();
//	vtkSmartPointer<vtkPolyData>Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//	Poly_Data->DeepCopy(vtk_Reader->GetOutput());
//	vtkSmartPointer<vtkPolyDataMapper>Poly_Mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	Poly_Mapper->SetInputData(Poly_Data);
//	Poly_Mapper->ScalarVisibilityOff();
//	vtkSmartPointer<vtkActor> Poly_actor =
//		vtkSmartPointer<vtkActor>::New();
//	Poly_actor->SetMapper(Poly_Mapper);
//	Poly_actor->GetProperty()->SetOpacity(0.5);
//	return Poly_actor;
//}
//void SavePointsToCSV(vtkSmartPointer<vtkPoints> points, const std::string& fileName) {
//	// 打开文件，若文件不存在则创建
//	std::ofstream outFile(fileName);
//
//	if (!outFile.is_open()) {
//		std::cerr << "无法打开文件进行写入: " << fileName << std::endl;
//		return;
//	}
//
//	// 写入CSV表头
//	outFile << "X,Y,Z\n";
//
//	// 遍历点，写入每个点的坐标
//	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
//		double* point = points->GetPoint(i);
//		outFile << point[0] << "," << point[1] << "," << point[2] << "\n";
//	}
//
//	// 关闭文件
//	outFile.close();
//	std::cout << "点已保存到文件: " << fileName << std::endl;
//}
//// 计算两点间的欧氏距离
//double calculateDistance(double* point1, double* point2) {
//	return sqrt(pow(point2[0] - point1[0], 2) + pow(point2[1] - point1[1], 2) + pow(point2[2] - point1[2], 2));
//}
//// 筛选距离小于指定值的皮肤点
//vtkSmartPointer<vtkPoints> filterPointsByDistance(vtkSmartPointer<vtkPoints> skinPoints, double* tumorCenter, double maxDistance) {
//	vtkSmartPointer<vtkPoints> filteredPoints = vtkSmartPointer<vtkPoints>::New();
//
//	// 遍历所有皮肤点
//	for (vtkIdType i = 0; i < skinPoints->GetNumberOfPoints(); i++) {
//		double* point = skinPoints->GetPoint(i);
//		double distance = calculateDistance(point, tumorCenter); // 计算当前点到肿瘤质心的距离
//		if (distance < maxDistance) {  // 如果距离小于指定最大值，则保留这个点
//			filteredPoints->InsertNextPoint(point);
//		}
//	}
//
//	return filteredPoints;  // 返回筛选后的点
//}
//// 计算两个向量之间的角度
//double calculateAngle(double pt1[3], double pt2[3]) {
//	double dotProduct = pt1[0] * pt2[0] + pt1[1] * pt2[1] + pt1[2] * pt2[2];
//	double magnitude1 = std::sqrt(pt1[0] * pt1[0] + pt1[1] * pt1[1] + pt1[2] * pt1[2]);
//	double magnitude2 = std::sqrt(pt2[0] * pt2[0] + pt2[1] * pt2[1] + pt2[2] * pt2[2]);
//
//	double cosTheta = dotProduct / (magnitude1 * magnitude2);
//	cosTheta = std::min(std::max(cosTheta, -1.0), 1.0);  // 防止浮动误差
//
//	return std::acos(cosTheta);  // 返回角度值（弧度）
//}
//
//// 筛选出满足角度条件的皮肤点
//vtkSmartPointer<vtkPoints> filterPointsByAngle(vtkSmartPointer<vtkPoints> skinPoints, double* tumorCenter, double angleThreshold) {
//	vtkSmartPointer<vtkPoints> filteredPoints = vtkSmartPointer<vtkPoints>::New();
//
//	// 遍历所有皮肤点
//	for (vtkIdType i = 0; i < skinPoints->GetNumberOfPoints(); ++i) {
//		double point[3];
//		skinPoints->GetPoint(i, point);  // 获取当前皮肤点
//
//		// 计算肿瘤质心和当前点之间的角度
//		double angle = calculateAngle(tumorCenter, point);
//
//		// 如果角度小于给定的阈值，保留该点
//		if (angle <= angleThreshold) {
//			filteredPoints->InsertNextPoint(point);  // 将该点添加到筛选结果中
//		}
//	}
//
//	return filteredPoints;  // 返回筛选后的vtkPoints对象
//}
//int main() {
//	//法向量的设定的形式要保证是单位向量
//
//	double angle_0 = -30;
//	//double angle_1 = 30;
//	double normal[3] = { cos(angle_0 * PI / 180.0) ,sin(angle_0 * PI / 180.0),0 };
//	//转化为的单位向量的形式
//	double total_distance = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
//
//	normal[0] = normal[0] / total_distance;
//	normal[1] = normal[1] / total_distance;
//	normal[2] = normal[2] / total_distance;
//	cout << "The distance :" << sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]) << endl;
//	cout << "normal[0]:" << normal[0] << '\t' << "normal[1]:" << normal[1] << '\t' << "normal[2]:" << normal[2] << endl;
//	//椭球体的生成
//
//
//
//	vtkSmartPointer<vtkRenderer> aRenderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renWin =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renWin->AddRenderer(aRenderer);
//	vtkSmartPointer<vtkRenderWindowInteractor> iren =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();//设置绘图窗口交互
//	iren->SetRenderWindow(renWin);
//	//vtk形式的数据的读取
//	//读取骨骼的stl格式的文件的形式（骨骼的计算）
//	vtkSmartPointer<vtkActor>poly_artery = Reader_VTK_Actor("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\MESHES_VTK\\", "artery.vtk");
//	vtkSmartPointer<vtkPolyData>poly_bone = Reader_VTK("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\MESHES_VTK\\", "bone.vtk");
//	vtkSmartPointer<vtkActor>poly_gallbladder = Reader_VTK_Actor("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\MESHES_VTK\\", "leftlung.vtk");
//	vtkSmartPointer<vtkActor>poly_liver = Reader_VTK_Actor("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\MESHES_VTK\\", "liver.vtk");
//	vtkSmartPointer<vtkActor>poly_liverkyste = Reader_VTK_Actor("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\MESHES_VTK\\", "liverkyst.vtk");
//	vtkSmartPointer<vtkPolyData>poly_tumor = Reader_VTK("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\MESHES_VTK\\", "livertumor01.vtk");
//	vtkSmartPointer<vtkActor>poly_portalvein = Reader_VTK_Actor("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\MESHES_VTK\\", "portalvein.vtk");
//	vtkSmartPointer<vtkPolyData>poly_skin = Reader_STL("D:\\\Data Disk\\LW\\Data\\3Dircadb1.1\\STL\\", "skin.stl");
//	vtkSmartPointer<vtkActor>poly_venacava = Reader_VTK_Actor("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\MESHES_VTK\\", "venoussystem.vtk");
//
//	//肿瘤的质心的形式的转化
//	double* center = poly_tumor->GetCenter();
//	//肿瘤的包围盒
//	double* bound = poly_tumor->GetBounds();
//	cout << "The center：" << center[0] << '\t' << center[1] << '\t' << center[2] << endl;
//	cout << "bound[0]:" << bound[0] << '\t' << "bound[1]：" << bound[1] << '\t' << "bound[2]:" << bound[2] << '\t' << "bound[3]:" << bound[3] << '\t' << "bound[4]:" << bound[4] << '\t' << "bound[5]:" << bound[5] << endl;
//	cout << "length_x:" << bound[1] - bound[0] << '\t' << "length_y :" << bound[3] - bound[2] << "length_z:" << bound[5] - bound[4] << endl;
//
//	int tumor_point_size = poly_tumor->GetNumberOfPoints();
//
//
//	//皮肤的裁剪的操作
//	vtkSmartPointer<vtkPoints>skin_Points = vtkSmartPointer<vtkPoints>::New();
//	vtkSmartPointer<vtkPoints>newskin_Points = vtkSmartPointer<vtkPoints>::New();
//	vtkSmartPointer<vtkPoints>bonenewskin_Points = vtkSmartPointer<vtkPoints>::New();
//	vtkSmartPointer<vtkPoints>newskin_Points2 = vtkSmartPointer<vtkPoints>::New();
//
//	vtkSmartPointer<vtkPoints>finnalnewskin_Points150 = vtkSmartPointer<vtkPoints>::New();
//
//	vtkSmartPointer<vtkPoints>original_points = vtkSmartPointer<vtkPoints>::New();
//	cout << "The size :" << poly_skin->GetNumberOfPoints() << endl;
//	for (int i = 0; i < poly_skin->GetNumberOfPoints(); i++) {
//		double* point = poly_skin->GetPoint(i);
//		//point[2]需要与肿瘤的中心相联系
//		if (point[2] >= (bound[4] - 50) && point[2] <= (bound[5] + 25) && point[0] <= (bound[1] + 50) && point[1] <= (bound[3] + 50))
//			skin_Points->InsertNextPoint(point);
//	}
//	for (int i = 0; i < poly_skin->GetNumberOfPoints(); i++) {
//		double* point = poly_skin->GetPoint(i);
//		//point[2]需要与肿瘤的中心相联系
//		if (point[2] >= (bound[4] - 50) && point[2] <= (bound[5] + 50))
//			original_points->InsertNextPoint(point);
//	}
//
//	vtkSmartPointer<vtkPolyData> skin_pointsPolydata =
//		vtkSmartPointer<vtkPolyData>::New();
//	skin_pointsPolydata->SetPoints(skin_Points);
//	vtkSmartPointer<vtkPolyData> newskin_pointsPolydata =
//		vtkSmartPointer<vtkPolyData>::New();
//	newskin_pointsPolydata->SetPoints(newskin_Points);
//	vtkSmartPointer<vtkPolyData> newskin_pointsPolydata2 =
//		vtkSmartPointer<vtkPolyData>::New();
//	newskin_pointsPolydata2->SetPoints(newskin_Points2);
//
//	vtkSmartPointer<vtkPolyData> original_skin_pointsPolydata =
//		vtkSmartPointer<vtkPolyData>::New();
//	original_skin_pointsPolydata->SetPoints(original_points);
//	//查看点皮肤点的个数
//
//	cout << "skin_points:" << skin_Points->GetNumberOfPoints() << endl;
//	cout << "original_points:" << original_points->GetNumberOfPoints() << endl;
//
//	vtkSmartPointer<vtkPolyDataReader> vtkReaderBone = vtkSmartPointer<vtkPolyDataReader>::New();
//	vtkSmartPointer<vtkPolyDataReader> vtkReaderartery = vtkSmartPointer<vtkPolyDataReader>::New();
//	vtkReaderBone->SetFileName("D://Data Disk//LW//Data//3Dircadb1.1//MyDo_VTK//bone.vtk");//读取VTK格式文件
//	vtkReaderBone->Update();//穿刺画线前的一个数据操作
//	vtkReaderartery->SetFileName("D://Data Disk//LW//Data//3Dircadb1.1//MyDo_VTK//artery.vtk");//读取VTK格式文件
//	vtkReaderartery->Update();//穿刺画线前的一个数据操作
//	// 3. 创建存储筛选点的容器
//	vtkNew<vtkOBBTree> tree;//骨骼与线的交点检测
//	vtkNew<vtkPoints> intersectPoints;
//	vtkSmartPointer<vtkPoints> newFinalSkinPoints = vtkSmartPointer<vtkPoints>::New();
//	for (int i = 0; i < skin_Points->GetNumberOfPoints(); i++) {
//		double* point = skin_Points->GetPoint(i);
//		vtkNew<vtkLineSource> line;
//		line->SetPoint1(center);
//		line->SetPoint2(point);
//
//		tree->SetDataSet(vtkReaderBone->GetOutput());
//		tree->BuildLocator();
//		//tree->SetDataSet(vtkReaderartery->GetOutput());
//	   // tree->BuildLocator();
//		bool pkx = tree->IntersectWithLine(center, point, intersectPoints, NULL);
//		// bool ppp=tree->IntersectWithLine(center, point, intersectPoints, NULL);
//		if (pkx == false) {
//			std::cout << "交点的坐标x，y，z : " << point[0] << " " << point[1] << " " << point[2] << std::endl;
//			newskin_Points->InsertNextPoint(point);
//		}
//		else {
//			std::cout << "碰撞骨头交点的坐标x，y，z : " << point[0] << " " << point[1] << " " << point[2] << std::endl;
//			bonenewskin_Points->InsertNextPoint(point);
//
//
//		}
//
//
//	}
//
//	
//
//
//	vtkSmartPointer<vtkPolyData> bonenewskin_pointsPolydata2 =
//		vtkSmartPointer<vtkPolyData>::New();
//	bonenewskin_pointsPolydata2->SetPoints(bonenewskin_Points);
//	vtkSmartPointer<vtkVertexGlyphFilter>boneskin =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	boneskin->SetInputData(bonenewskin_pointsPolydata2); //可视化点就在这里修改
//	boneskin->Update();
//	vtkSmartPointer<vtkPolyDataMapper>boneMapper1 =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	boneMapper1->SetInputData(boneskin->GetOutput());
//	vtkSmartPointer<vtkActor> boneskin_actor1 =
//		vtkSmartPointer<vtkActor>::New();
//	boneskin_actor1->SetMapper(boneMapper1);
//	boneskin_actor1->GetProperty()->SetColor(0, 0, 0);
//	boneskin_actor1->GetProperty()->SetPointSize(5);
//
//	cout << "newskin_points:" << newskin_Points->GetNumberOfPoints() << endl;
//
//	for (int i = 0; i < newskin_Points->GetNumberOfPoints(); i++) {
//		double* point = newskin_Points->GetPoint(i);
//		vtkNew<vtkLineSource> line;
//		line->SetPoint1(center);
//		line->SetPoint2(point);
//
//		//tree->SetDataSet(vtkReaderBone->GetOutput());
//		//tree->BuildLocator();
//		tree->SetDataSet(vtkReaderartery->GetOutput());
//		tree->BuildLocator();
//		bool ppp = tree->IntersectWithLine(center, point, intersectPoints, NULL);
//		// bool ppp=tree->IntersectWithLine(center, point, intersectPoints, NULL);
//		if (ppp == false) {
//			std::cout << "交点的坐标x，y，z : " << point[0] << " " << point[1] << " " << point[2] << std::endl;
//			newskin_Points2->InsertNextPoint(point);
//		}
//
//	}
//
//	cout << "newskin_points2:" << newskin_Points2->GetNumberOfPoints() << endl;
//
//	// 调用过滤函数，保留距离小于150的皮肤点
//	finnalnewskin_Points150 = filterPointsByDistance(newskin_Points2, center, 100.0);
//	cout << "finnalnewskin_points2:" << finnalnewskin_Points150->GetNumberOfPoints() << endl;
//
//	vtkSmartPointer<vtkPolyData> finnalnewskin_pointsPolydata150 =
//		vtkSmartPointer<vtkPolyData>::New();
//	finnalnewskin_pointsPolydata150->SetPoints(finnalnewskin_Points150);
//	//
//	////	设置角度阈值（以弧度为单位）
//	//		double angleThreshold = 20.0 * (PI / 180.0);  // 30度转为弧度
//	//
//	//		// 筛选满足角度条件的点
//	//		vtkSmartPointer<vtkPoints> filteredPoints = filterPointsByAngle(finnalnewskin_Points150, center, angleThreshold);
//	//
//	//	vtkSmartPointer<vtkPolyData> jiaodu =
//	//		vtkSmartPointer<vtkPolyData>::New();
//	//	jiaodu->SetPoints(filteredPoints);
//
//		// Step 1: 画线
//	vtkSmartPointer<vtkPoints> linepoints = vtkSmartPointer<vtkPoints>::New();
//	linepoints->InsertNextPoint(47.2498 ,74.0783, 118.938);  // 点1 (起点)
//	linepoints->InsertNextPoint(center);  // 点2 (终点)
//	 // Step 2: 创建一个单一的线段
//	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
//	lines->InsertNextCell(2);  // 线段需要两个点
//	lines->InsertCellPoint(0);
//	lines->InsertCellPoint(1);
//
//	// Step 3: 创建 PolyData，保存点和线信息
//	vtkSmartPointer<vtkPolyData> lineData = vtkSmartPointer<vtkPolyData>::New();
//	lineData->SetPoints(linepoints);
//	lineData->SetLines(lines);
//
//	// Step 4: 创建映射器
//	vtkSmartPointer<vtkPolyDataMapper> lineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	lineMapper->SetInputData(lineData);
//
//	// Step 5: 创建演员
//	vtkSmartPointer<vtkActor> lineActor = vtkSmartPointer<vtkActor>::New();
//	lineActor->SetMapper(lineMapper);
//
//
//	vtkSmartPointer<vtkLineSource>line = vtkSmartPointer<vtkLineSource>::New();
//	line->SetPoint1(47.2498, 74.0783, 118.938);
//	line->SetPoint2(center);
//	vtkSmartPointer<vtkTubeFilter>tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
//	tubeFilter->SetInputConnection(line->GetOutputPort());
//	tubeFilter->SetRadius(0.5);
//	tubeFilter->SetNumberOfSides(100);
//	tubeFilter->CappingOn();
//	tubeFilter->Update();
//	//圆柱体的三角面片的形式转化
//	vtkSmartPointer<vtkTriangleFilter> triFilter_source =
//		vtkSmartPointer<vtkTriangleFilter>::New();
//	triFilter_source->SetInputData(tubeFilter->GetOutput());
//	triFilter_source->Update();
//
//	//圆柱体的三角面片生成后的体积的计算的形式的展示
//	vtkSmartPointer<vtkMassProperties>massProp_source =
//		vtkSmartPointer<vtkMassProperties>::New();
//	massProp_source->SetInputData(triFilter_source->GetOutput());
//	float vol_source = massProp_source->GetVolume();
//	cout << "The original vol:" << vol_source << endl;
//
//	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper->SetInputData(tubeFilter->GetOutput());
//	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//	actor->GetProperty()->SetColor(0, 0, 0);
//	vtkSmartPointer<vtkPoints>cu1 = vtkSmartPointer<vtkPoints>::New();
//	// 设置聚类 1 的中心坐标
//	double clusterCenter1[3] = { 85.26,80.43, 117.59 };
//
//	// 将中心坐标添加到 vtkPoints 中
//	cu1->InsertNextPoint(clusterCenter1);
//
//	vtkSmartPointer<vtkPoints>cu2 = vtkSmartPointer<vtkPoints>::New();
//	double clusterCenter2[3] = { 105.73,  73.04, 117.30 };
//
//	// 将中心坐标添加到 vtkPoints 中
//	cu2->InsertNextPoint(clusterCenter2);
//	vtkSmartPointer<vtkPoints>cu3 = vtkSmartPointer<vtkPoints>::New();
//	double clusterCenter3[3] = { 96.14, 91.82,  122.79 };
//
//	// 将中心坐标添加到 vtkPoints 中
//	cu3->InsertNextPoint(clusterCenter3);
//	vtkSmartPointer<vtkPoints>cu4 = vtkSmartPointer<vtkPoints>::New();
//	double clusterCenter4[3] = { 111.84,  87.48,  102.92 };
//
//	// 将中心坐标添加到 vtkPoints 中
//	cu4->InsertNextPoint(clusterCenter4);
//	vtkSmartPointer<vtkPoints>cu5 = vtkSmartPointer<vtkPoints>::New();
//	double clusterCenter5[3] = { 112.77, 71.49, 109.66 };
//
//	// 将中心坐标添加到 vtkPoints 中
//	cu5->InsertNextPoint(clusterCenter5);
//
//	vtkSmartPointer<vtkPolyData> cuPolydata1 =
//		vtkSmartPointer<vtkPolyData>::New();
//	cuPolydata1->SetPoints(cu1);
//	vtkSmartPointer<vtkPolyData> cuPolydata2 =
//		vtkSmartPointer<vtkPolyData>::New();
//	cuPolydata2->SetPoints(cu2);
//	vtkSmartPointer<vtkPolyData> cuPolydata3 =
//		vtkSmartPointer<vtkPolyData>::New();
//	cuPolydata3->SetPoints(cu3);
//	vtkSmartPointer<vtkPolyData> cuPolydata4 =
//		vtkSmartPointer<vtkPolyData>::New();
//	cuPolydata4->SetPoints(cu4);
//	vtkSmartPointer<vtkPolyData> cuPolydata5 =
//		vtkSmartPointer<vtkPolyData>::New();
//	cuPolydata5->SetPoints(cu5);
//
//	vtkSmartPointer<vtkVertexGlyphFilter>cufilter1 =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	cufilter1->SetInputData(cuPolydata1); //可视化点就在这里修改
//	cufilter1->Update();
//	vtkSmartPointer<vtkVertexGlyphFilter>cufilter2 =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	cufilter2->SetInputData(cuPolydata2); //可视化点就在这里修改
//	cufilter2->Update();
//	vtkSmartPointer<vtkVertexGlyphFilter>cufilter3 =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	cufilter3->SetInputData(cuPolydata3); //可视化点就在这里修改
//	cufilter3->Update();
//	vtkSmartPointer<vtkVertexGlyphFilter>cufilter4 =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	cufilter4->SetInputData(cuPolydata4); //可视化点就在这里修改
//	cufilter4->Update();
//	vtkSmartPointer<vtkVertexGlyphFilter>cufilter5 =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	cufilter5->SetInputData(cuPolydata5); //可视化点就在这里修改
//	cufilter5->Update();
//	vtkSmartPointer<vtkPolyDataMapper>cuMapper1 =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	cuMapper1->SetInputData(cufilter1->GetOutput());
//	vtkSmartPointer<vtkPolyDataMapper>cuMapper2 =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	cuMapper2->SetInputData(cufilter2->GetOutput());
//	vtkSmartPointer<vtkPolyDataMapper>cuMapper3 =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	cuMapper3->SetInputData(cufilter3->GetOutput());
//	vtkSmartPointer<vtkPolyDataMapper>cuMapper4 =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	cuMapper4->SetInputData(cufilter4->GetOutput());
//	vtkSmartPointer<vtkPolyDataMapper>cuMapper5 =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	cuMapper5->SetInputData(cufilter5->GetOutput());
//
//	vtkSmartPointer<vtkActor> cu_actor1 =
//		vtkSmartPointer<vtkActor>::New();
//	cu_actor1->SetMapper(cuMapper1);
//	cu_actor1->GetProperty()->SetColor(0, 1, 0);
//	cu_actor1->GetProperty()->SetPointSize(5);
//	vtkSmartPointer<vtkActor> cu_actor2 =
//		vtkSmartPointer<vtkActor>::New();
//	cu_actor2->SetMapper(cuMapper2);
//	cu_actor2->GetProperty()->SetColor(1.0, 1.0, 0.0);
//	cu_actor2->GetProperty()->SetPointSize(5);
//	vtkSmartPointer<vtkActor> cu_actor3 =
//		vtkSmartPointer<vtkActor>::New();
//	cu_actor3->SetMapper(cuMapper3);
//	cu_actor3->GetProperty()->SetColor(0, 0, 1);
//	cu_actor3->GetProperty()->SetPointSize(5);
//	vtkSmartPointer<vtkActor> cu_actor4 =
//		vtkSmartPointer<vtkActor>::New();
//	cu_actor4->SetMapper(cuMapper4);
//	cu_actor4->GetProperty()->SetColor(1, 0, 0);
//	cu_actor4->GetProperty()->SetPointSize(5);
//	vtkSmartPointer<vtkActor> cu_actor5 =
//		vtkSmartPointer<vtkActor>::New();
//	cu_actor5->SetMapper(cuMapper5);
//	cu_actor5->GetProperty()->SetColor(0.5, 0, 0.5);
//	cu_actor5->GetProperty()->SetPointSize(5);
//
//	vtkSmartPointer<vtkVertexGlyphFilter>vert_filter =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	vert_filter->SetInputData(finnalnewskin_pointsPolydata150); //可视化点就在这里修改
//	vert_filter->Update();
//	vtkSmartPointer<vtkPolyDataMapper>pointsMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	pointsMapper->SetInputData(vert_filter->GetOutput());
//	pointsMapper->ScalarVisibilityOff();
//	vtkSmartPointer<vtkActor> skined_actor =
//		vtkSmartPointer<vtkActor>::New();
//	skined_actor->SetMapper(pointsMapper);
//	skined_actor->GetProperty()->SetColor(0, 1.0, 0);
//	skined_actor->GetProperty()->SetPointSize(5);
//	//points_actor->GetProperty()->SetOpacity(0.5);
//	//原始皮肤节点的绘制
//	/*vtkSmartPointer<vtkVertexGlyphFilter>vert_filter =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	vert_filter->SetInputData(original_skin_pointsPolydata);
//	vert_filter->Update();
//	vtkSmartPointer<vtkPolyDataMapper>pointsMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	pointsMapper->SetInputData(vert_filter->GetOutput());
//	pointsMapper->ScalarVisibilityOff();
//	vtkSmartPointer<vtkActor> skined_actor =
//		vtkSmartPointer<vtkActor>::New();
//	skined_actor->SetMapper(pointsMapper);
//	skined_actor->GetProperty()->SetColor(0.5, 0.5, 0);
//	skined_actor->GetProperty()->SetPointSize(2.5);*/
//	//骨头的绘制的方式
//	vtkSmartPointer<vtkPolyDataMapper>bone_Mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	bone_Mapper->SetInputData(poly_bone);
//	bone_Mapper->ScalarVisibilityOff();
//	vtkSmartPointer<vtkActor> bone_actor =
//		vtkSmartPointer<vtkActor>::New();
//	bone_actor->SetMapper(bone_Mapper);
//	bone_actor->GetProperty()->SetColor(1.0, 1.0, 1.0);
//	//bone_actor->GetProperty()->SetOpacity(0.5);
//
//	//肿瘤的展示的形式 poly_tumor
//	vtkSmartPointer<vtkPolyDataMapper>tumor_Mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	tumor_Mapper->SetInputData(poly_tumor);
//	tumor_Mapper->ScalarVisibilityOff();
//	vtkSmartPointer<vtkActor> tumor_actor =
//		vtkSmartPointer<vtkActor>::New();
//	tumor_actor->SetMapper(tumor_Mapper);
//	tumor_actor->GetProperty()->SetColor(0.5, 0, 1);
//	tumor_actor->GetProperty()->SetOpacity(0.7);
//
//	vtkSmartPointer<vtkCamera> aCamera =
//		vtkSmartPointer<vtkCamera>::New();
//	aCamera->SetViewUp(0, 0, -1);//视图
//	aCamera->SetPosition(0, 1, 0);//位置
//	aCamera->SetFocalPoint(0, 0, 0);//焦点
//	aCamera->ComputeViewPlaneNormal();//视图平面法线;
//	aCamera->Azimuth(30.0);
//	aCamera->Elevation(30.0);
//	aCamera->Dolly(1.5);
//	aRenderer->AddActor(bone_actor);
//	//aRenderer->AddActor(boneskin_actor1);
//	aRenderer->AddActor(actor);
//	//aRenderer->AddActor(lineActor);
//	//aRenderer->AddActor(skined_actor);
//	aRenderer->AddActor(tumor_actor);
//	//aRenderer->AddActor(cu_actor1);
//	//aRenderer->AddActor(cu_actor2);
//	//aRenderer->AddActor(cu_actor3);
//	//aRenderer->AddActor(cu_actor4);
//	//aRenderer->AddActor(cu_actor5);
//
//	poly_artery->GetProperty()->SetColor(1, 0.5, 0);
//	aRenderer->AddActor(poly_artery);
//
//	poly_liverkyste->GetProperty()->SetColor(0.4, 0.7, 0.2);
//	aRenderer->AddActor(poly_liverkyste);
//	//liver肝脏（保留）
//	poly_liver->GetProperty()->SetColor(1.0, 1.0, 1.0);
//	aRenderer->AddActor(poly_liver);
//	//portalvein门静脉（保留）
//	poly_portalvein->GetProperty()->SetColor(1, 0.2, 0);
//	aRenderer->AddActor(poly_portalvein);
//	//spleen脾脏（保留）
//	//venoussystem静脉系统(保留)
//	//liver肝脏（保留）
//	poly_liver->GetProperty()->SetColor(1.0, 1.0, 1.0);
//	aRenderer->AddActor(poly_liver);
//
//	//portalvein门静脉（保留）
//	poly_portalvein->GetProperty()->SetColor(1, 0.2, 0);
//	aRenderer->AddActor(poly_portalvein);
//	//venoussystem静脉系统(保留)
//	poly_venacava->GetProperty()->SetColor(1, 0, 0.2);
//	aRenderer->AddActor(poly_venacava);
//	poly_gallbladder->GetProperty()->SetColor(1, 0.4, 0.2);
//	aRenderer->AddActor(poly_gallbladder);
//	//aRenderer->AddActor(selected_skin_points_actor);
//
//	aCamera->SetPosition(0, 0, 1);  // 设置相机位置
//	aCamera->SetFocalPoint(0, 0, 0); // 设置相机焦点（目标位置）
//	aCamera->SetViewUp(0, 1, 0);    // 设置“上”方向，通常是 (0, 1, 0) 表示“Y轴”向上
//
//
//	aRenderer->SetActiveCamera(aCamera);//设置相机
//	aRenderer->ResetCamera();//重置相机
//	aRenderer->SetBackground(1, 1, 1);//设置背景颜色，double类型
//	aRenderer->ResetCameraClippingRange();
//	renWin->Render();
//	iren->Initialize();
//	iren->Start();
//	return 0;
//}
