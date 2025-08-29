//#include <vtkAutoInit.h>
//VTK_MODULE_INIT(vtkRenderingOpenGL2);
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingFreeType);
//#include "vtkDICOMImageReader.h"//DCMҽѧ�ļ���ȡ��
//#include "vtkPolyDataWriter.h"//����Ϊ.vtkͼ����
//#include "vtkPolyDataReader.h"//��ȡ.vtkͼ����
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
//using namespace std;
//const double PI = 3.141592653589793;//����һ�����ɸı�ĳ���ֵPI
//const int Po_Size = 50;//��Ⱥ��ģ
//const int Ev_Algebra = 18000;//��������
//const double Ov_Probability = 0.850; //�������,������������ж����������Ƿ���Ҫ����
//const double Va_Probability = 0.10;//�������,������������ж���һ�����Ƿ���Ҫ����
//vector<double>line_to_center;
//double tube_radius = 100;
//double Ellipsoid_long = 20;
//double Ellipsoid_short = 15;
//
////point ��actor �Ĺ�������ʽ
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
//	pointsActor->GetProperty()->SetPointSize(0.1);//�����ĳߴ��С������������ڻ�������ʾ����
//	return pointsActor;
//}
//
//
////������ȡSTL�ļ��ĺ���
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
////������ȡVTK �ļ��ĺ���
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
////������Ӧ��VTK �ļ���actor ����ʽ
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
//
//int main() {
//	//���������趨����ʽҪ��֤�ǵ�λ����
//
//	double angle_0 = -30;
//	//double angle_1 = 30;
//	double normal[3] = { cos(angle_0 * PI / 180.0) ,sin(angle_0 * PI / 180.0),0 };
//	//ת��Ϊ�ĵ�λ��������ʽ
//	double total_distance = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
//
//	normal[0] = normal[0] / total_distance;
//	normal[1] = normal[1] / total_distance;
//	normal[2] = normal[2] / total_distance;
//	cout << "The distance :" << sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]) << endl;
//	cout << "normal[0]:" << normal[0] << '\t' << "normal[1]:" << normal[1] << '\t' << "normal[2]:" << normal[2] << endl;
//	//�����������
//
//
//
//	vtkSmartPointer<vtkRenderer> aRenderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renWin =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renWin->AddRenderer(aRenderer);
//	vtkSmartPointer<vtkRenderWindowInteractor> iren =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();//���û�ͼ���ڽ���
//	iren->SetRenderWindow(renWin);
//	//vtk��ʽ�����ݵĶ�ȡ
//	//��ȡ������stl��ʽ���ļ�����ʽ�������ļ��㣩
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
//	//���������ĵ���ʽ��ת��
//	double* center = poly_tumor->GetCenter();
//	//�����İ�Χ��
//	double* bound = poly_tumor->GetBounds();
//	cout << "The center��" << center[0] << '\t' << center[1] << '\t' << center[2] << endl;
//	cout << "bound[0]:" << bound[0] << '\t' << "bound[1]��" << bound[1] << '\t'<<"bound[2]:" << bound[2]<< '\t'<< "bound[3]:" << bound[3] << '\t'<<"bound[4]:"<<bound[4]<< '\t'<<"bound[5]:"<<bound[5] << endl;
//	cout << "length_x:" << bound[1] - bound[0] << '\t' << "length_y :" << bound[3] - bound[2] << "length_z:" << bound[5] - bound[4] << endl;
//
//	int tumor_point_size = poly_tumor->GetNumberOfPoints();
//	//Ƥ���Ĳü��Ĳ���
//	vtkSmartPointer<vtkPoints>skin_Points = vtkSmartPointer<vtkPoints>::New();
//	vtkSmartPointer<vtkPoints>original_points = vtkSmartPointer<vtkPoints>::New();
//	cout << "The size :" << poly_skin->GetNumberOfPoints() << endl;
//	for (int i = 0; i < poly_skin->GetNumberOfPoints(); i++) {
//		double* point = poly_skin->GetPoint(i);
//		//point[2]��Ҫ����������������ϵ
//		if (point[2] >= (bound[4] - 30) && point[2] <= (bound[5] + 40) && point[0] <= (bound[1] + 30) && point[1] <= (bound[3] + 60))
//			skin_Points->InsertNextPoint(point);
//	}
//	for (int i = 0; i < poly_skin->GetNumberOfPoints(); i++) {
//		double* point = poly_skin->GetPoint(i);
//		//point[2]��Ҫ����������������ϵ
//		if (point[2] >= (bound[4] - 30) && point[2] <= (bound[5] + 120))
//			original_points->InsertNextPoint(point);
//	}
//
//	vtkSmartPointer<vtkPolyData> skin_pointsPolydata =
//		vtkSmartPointer<vtkPolyData>::New();
//	skin_pointsPolydata->SetPoints(skin_Points);
//	vtkSmartPointer<vtkPolyData> original_skin_pointsPolydata =
//		vtkSmartPointer<vtkPolyData>::New();
//	original_skin_pointsPolydata->SetPoints(original_points);
//	//�鿴��Ƥ����ĸ���
//
//	cout << "skin_points:" << skin_Points->GetNumberOfPoints() << endl;
//	cout << "original_points:" << original_points->GetNumberOfPoints() << endl;
//
//	//ԭʼƤ���ڵ�Ļ���
//	vtkSmartPointer<vtkVertexGlyphFilter>vert_filter =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	vert_filter->SetInputData(skin_pointsPolydata);
//	vert_filter->Update();
//	vtkSmartPointer<vtkPolyDataMapper>pointsMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	pointsMapper->SetInputData(vert_filter->GetOutput());
//	pointsMapper->ScalarVisibilityOff();
//	vtkSmartPointer<vtkActor> skined_actor =
//		vtkSmartPointer<vtkActor>::New();
//	skined_actor->SetMapper(pointsMapper);
//	skined_actor->GetProperty()->SetColor(0.5, 0.5, 0);
//	skined_actor->GetProperty()->SetPointSize(2.5);
//	//points_actor->GetProperty()->SetOpacity(0.5);
//	//ԭʼƤ���ڵ�Ļ���
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
//	//��ͷ�Ļ��Ƶķ�ʽ
//	vtkSmartPointer<vtkPolyDataMapper>bone_Mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	bone_Mapper->SetInputData(poly_bone);
//	bone_Mapper->ScalarVisibilityOff();
//	vtkSmartPointer<vtkActor> bone_actor =
//		vtkSmartPointer<vtkActor>::New();
//	bone_actor->SetMapper(bone_Mapper);
//	bone_actor->GetProperty()->SetColor(0.50, 0.73, 1.0);
//	//bone_actor->GetProperty()->SetOpacity(0.5);
//
//	//������չʾ����ʽ poly_tumor
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
//	aCamera->SetViewUp(0, 0, -1);//��ͼ
//	aCamera->SetPosition(0, 1, 0);//λ��
//	aCamera->SetFocalPoint(0, 0, 0);//����
//	aCamera->ComputeViewPlaneNormal();//��ͼƽ�淨��;
//	aCamera->Azimuth(30.0);
//	aCamera->Elevation(30.0);
//	aCamera->Dolly(1.5);
//	aRenderer->AddActor(bone_actor);
//	aRenderer->AddActor(skined_actor);
//	aRenderer->AddActor(tumor_actor);
//
//	poly_artery->GetProperty()->SetColor(1, 0.5, 0);
//	aRenderer->AddActor(poly_artery);
//
//	poly_liverkyste->GetProperty()->SetColor(0.4, 0.7, 0.2);
//	aRenderer->AddActor(poly_liverkyste);
//	//liver���ࣨ������
//	poly_liver->GetProperty()->SetColor(0, 0.5, 0);
//	aRenderer->AddActor(poly_liver);
//	//portalvein�ž�����������
//	poly_portalvein->GetProperty()->SetColor(1, 0.2, 0);
//	aRenderer->AddActor(poly_portalvein);
//	//spleenƢ�ࣨ������
//	//venoussystem����ϵͳ(����)
//	//liver���ࣨ������
//	poly_liver->GetProperty()->SetColor(0, 0.5, 0);
//	aRenderer->AddActor(poly_liver);
//
//	//portalvein�ž�����������
//	poly_portalvein->GetProperty()->SetColor(1, 0.2, 0);
//	aRenderer->AddActor(poly_portalvein);
//	//venoussystem����ϵͳ(����)
//	poly_venacava->GetProperty()->SetColor(1, 0, 0.2);
//	aRenderer->AddActor(poly_venacava);
//	poly_gallbladder->GetProperty()->SetColor(1, 0.4, 0.2);
//	aRenderer->AddActor(poly_gallbladder);
//	//aRenderer->AddActor(selected_skin_points_actor);
//	aRenderer->SetActiveCamera(aCamera);//�������
//	aRenderer->ResetCamera();//�������
//	aRenderer->SetBackground(.2, .3, .4);//���ñ�����ɫ��double����
//	aRenderer->ResetCameraClippingRange();
//	renWin->Render();
//	iren->Initialize();
//	iren->Start();
//	return 0;
//}
