//#include <vtkAutoInit.h>
//VTK_MODULE_INIT(vtkRenderingOpenGL2);
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingFreeType);
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
//#include<vtkLineSource.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include<vtkParametricEllipsoid.h>
//#include<vtkParametricFunctionSource.h>
//using namespace std;
//double tube_radius = 15;
//double Ellipsoid_long = 20;
//double Ellipsoid_short = 15;
////��õ�λ����
//void GetNormVector(double* Direction, double* p1, double* p2)
//{
//	//����������ľ���ó��Ľ���� �����ƽ������ʽ
//	double length = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
//	cout << "The length:" << length << endl;
//	Direction[0] = (p2[0] - p1[0]) / length;
//	Direction[1] = (p2[1] - p1[1]) / length;
//	Direction[2] = (p2[2] - p1[2]) / length;
//}
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
////�������ת��
//vtkSmartPointer<vtkPolyData> trans_polydata_Ellipsoid(vtkSmartPointer<vtkTransformPolyDataFilter> transformed_filter) {
//	vtkSmartPointer<vtkPolyData>Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//	Poly_Data->DeepCopy(transformed_filter->GetOutput());
//	return Poly_Data;
//}
////����������е�ת��
//vector<vtkSmartPointer<vtkPolyData>> trans_polydata_Ellipsoid_vector(vector<vtkSmartPointer<vtkTransformPolyDataFilter>> transformed_filters) {
//	vector<vtkSmartPointer<vtkPolyData>> result;
//	for (vtkSmartPointer<vtkTransformPolyDataFilter>transformed_filter : transformed_filters) {
//		vtkSmartPointer<vtkPolyData>Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//		Poly_Data->DeepCopy(transformed_filter->GetOutput());
//		result.push_back(Poly_Data);
//	}
//	return result;
//}
////Բ�����ת��
//vtkSmartPointer<vtkPolyData> trans_polydata_tube(vtkSmartPointer<vtkTubeFilter>tubeFilter) {
//	vtkSmartPointer<vtkPolyData>Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//	Poly_Data->DeepCopy(tubeFilter->GetOutput());
//	return Poly_Data;
//}
////����actor �ĺ���(ͨ��)
//vtkSmartPointer<vtkActor> Build_Actor(vtkSmartPointer<vtkPolyData> PolyData) {
//	vtkSmartPointer<vtkPolyDataMapper> Mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	Mapper->SetInputData(PolyData);//�����������ӣ�����˿ڻ�ȡ����vtkReader������	
//	Mapper->ScalarVisibilityOff();    //�����������ɫ
//	vtkSmartPointer<vtkActor> Actor =
//		vtkSmartPointer<vtkActor>::New();
//	Actor->SetMapper(Mapper);//����ӳ����
//	return Actor;
//}
////����actor ���еĺ���(ͨ�õİ汾�������������Լ� Բ��������)
//vector<vtkSmartPointer<vtkActor>> Build_Actor_vector(vector<vtkSmartPointer<vtkPolyData>> PolyDatas) {
//	vector<vtkSmartPointer<vtkActor>> Actors;
//	for (vtkSmartPointer<vtkPolyData > PolyData : PolyDatas) {
//		vtkSmartPointer<vtkPolyDataMapper> Mapper =
//			vtkSmartPointer<vtkPolyDataMapper>::New();
//		Mapper->SetInputData(PolyData);//�����������ӣ�����˿ڻ�ȡ����vtkReader������	
//		Mapper->ScalarVisibilityOff();    //�����������ɫ
//		vtkSmartPointer<vtkActor> Actor =
//			vtkSmartPointer<vtkActor>::New();
//		Actor->SetMapper(Mapper);//����ӳ����
//		Actor->GetProperty()->SetColor(1, 0, 0);
//		Actor->GetProperty()->SetOpacity(0.5);
//		Actors.push_back(Actor);
//	}
//	return Actors;
//}
////Ѱ�� STL ģ�͵����ĵ�(������ʵ��)
//void Get_center(vtkSmartPointer<vtkPolyData> PolyData, double* center) {
//	double* bounds = PolyData->GetBounds();
//	double x_len = bounds[1] - bounds[0];
//	double y_len = bounds[3] - bounds[2];
//	double z_len = bounds[5] - bounds[4];
//	center[0] = bounds[0] + x_len / 2;
//	center[1] = bounds[2] + y_len / 2;
//	center[2] = bounds[4] + z_len / 2;
//	cout << "The bound_x:" << bounds[0] << '\t' << bounds[1] << endl;
//	cout << "The bound_y:" << bounds[2] << '\t' << bounds[3] << endl;
//	cout << "The bound_z:" << bounds[4] << '\t' << bounds[5] << endl;
//};
////һ��ֱ�ߺ�һ�����ݼ��Ľ�����������
//vtkSmartPointer<vtkPoints> Intersect_points(vtkSmartPointer<vtkPolyData>poly_data, double* start_point, double* end_point) {
//	vtkSmartPointer<vtkOBBTree> obbtree = vtkSmartPointer<vtkOBBTree>::New();
//	obbtree->SetDataSet(poly_data);
//	obbtree->BuildLocator();
//	vtkSmartPointer<vtkPoints>intersect_point = vtkSmartPointer<vtkPoints>::New();
//	obbtree->IntersectWithLine(start_point, end_point, intersect_point, NULL);
//	return intersect_point;
//}
////3����λ�����ݵ�����ݸ�ʽ��ת������
//vtkSmartPointer<vtkVertexGlyphFilter> Insert_point(double* left_point, double* mid_point, double* right_point) {
//	vtkSmartPointer<vtkPoints>Points = vtkSmartPointer<vtkPoints>::New();
//	Points->InsertNextPoint(left_point);
//	Points->InsertNextPoint(mid_point);
//	Points->InsertNextPoint(right_point);
//	vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
//	pointsPolydata->SetPoints(Points);
//	vtkSmartPointer<vtkVertexGlyphFilter>vert_filter =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//	vert_filter->SetInputData(pointsPolydata);
//	vert_filter->Update();
//	return vert_filter;
//}
////���任����Ȼ��ת���������λ�����
//vtkSmartPointer<vtkTransformPolyDataFilter> get_transfilter(vtkSmartPointer<vtkParametricFunctionSource> source, vtkSmartPointer<vtkVertexGlyphFilter>source_filter, vtkSmartPointer<vtkVertexGlyphFilter>target_filter) {
//	vtkSmartPointer<vtkIterativeClosestPointTransform>icp =
//		vtkSmartPointer<vtkIterativeClosestPointTransform>::New();
//	icp->SetSource(source_filter->GetOutput());
//	icp->SetTarget(target_filter->GetOutput());
//	icp->GetLandmarkTransform()->SetModeToRigidBody();
//	icp->SetMaximumNumberOfIterations(50);
//	icp->StartByMatchingCentroidsOn();
//	icp->Modified();
//	icp->Update();
//	vtkSmartPointer<vtkMatrix4x4> m = icp->GetMatrix();
//	//��������һ������Ĳ���
//	cout << "The result :" << *m << endl;
//	vtkSmartPointer<vtkTransform> trans = vtkSmartPointer<vtkTransform>::New();
//	trans->SetMatrix(m);
//	//�������ת���Ĺ�ϵ�����
//	vtkSmartPointer<vtkTransformPolyDataFilter> transformfilter =
//		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//	transformfilter->SetInputConnection(source->GetOutputPort());
//	transformfilter->SetTransform(trans);
//	transformfilter->Update();
//	return transformfilter;
//}
////����ֱ������һϵ�е�������
//vector<vtkSmartPointer<vtkTransformPolyDataFilter>> get_transfilter_vector(vtkSmartPointer<vtkParametricFunctionSource> source, vtkSmartPointer<vtkVertexGlyphFilter>source_filter, vtkSmartPointer<vtkVertexGlyphFilter>target_filter, double* left, double* center, double* right, double* normal) {
//	vector<vtkSmartPointer<vtkTransformPolyDataFilter>> result;
//	//�õ��任��ϵ����ʽ��չʾ
//	vtkSmartPointer<vtkIterativeClosestPointTransform>icp =
//		vtkSmartPointer<vtkIterativeClosestPointTransform>::New();
//	icp->SetSource(source_filter->GetOutput());
//	icp->SetTarget(target_filter->GetOutput());
//	icp->GetLandmarkTransform()->SetModeToRigidBody();
//	icp->SetMaximumNumberOfIterations(50);
//	icp->StartByMatchingCentroidsOn();
//	icp->Modified();
//	icp->Update();
//	vtkSmartPointer<vtkMatrix4x4> m = icp->GetMatrix();
//	//��������һ������Ĳ���
//	cout << "The result :" << *m << endl;
//	vtkSmartPointer<vtkTransform> trans = vtkSmartPointer<vtkTransform>::New();
//	trans->SetMatrix(m);
//	//�������ת���Ĺ�ϵ�����
//	vtkSmartPointer<vtkTransformPolyDataFilter> transformfilter =
//		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//	transformfilter->SetInputConnection(source->GetOutputPort());
//	transformfilter->SetTransform(trans);
//	transformfilter->Update();
//	result.push_back(transformfilter);
//	//int left_length = vtkMath::Distance2BetweenPoints(left, center);
//	int i = 1;
//	int count = Ellipsoid_long - 10;
//	cout << "The center[0]:" << center[0] << endl;
//	cout << "The left[0]:" << left[0] << endl;
//	cout << "The normal[0]:" << normal[0] << endl;
//	while ((center[0] - i * count * normal[0]) >= left[0]) {
//		vtkSmartPointer<vtkTransform>transform = vtkSmartPointer<vtkTransform>::New();
//		transform->Translate(-i * count * normal[0], -i * count * normal[1], -i * count * normal[2]);
//		vtkSmartPointer<vtkTransformPolyDataFilter> transformfilter_temp =
//			vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//		transformfilter_temp->SetInputData(transformfilter->GetOutput());
//		transformfilter_temp->SetTransform(transform);
//		transformfilter_temp->Update();
//		result.push_back(transformfilter_temp);
//		i++;
//	}
//	i = 1;
//	cout << "The center[0]:" << center[0] << endl;
//	cout << "The right[0]:" << right[0] << endl;
//	while ((center[0] + i * count * normal[0]) <= right[0]) {
//		vtkSmartPointer<vtkTransform>transform = vtkSmartPointer<vtkTransform>::New();
//		transform->Translate(i * count * normal[0], i * count * normal[1], i * count * normal[2]);
//		vtkSmartPointer<vtkTransformPolyDataFilter> transformfilter_temp =
//			vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//		transformfilter_temp->SetInputData(transformfilter->GetOutput());
//		transformfilter_temp->SetTransform(transform);
//		transformfilter_temp->Update();
//		result.push_back(transformfilter_temp);
//		i++;
//	}
//	return result;
//}
////Բ��������ɵĲ���
//vtkSmartPointer<vtkTubeFilter> Create_Tube(double r, double* start_point, double* end_point) {
//	vtkSmartPointer<vtkLineSource>line = vtkSmartPointer<vtkLineSource>::New();
//	line->SetPoint1(start_point);
//	line->SetPoint2(end_point);
//	vtkSmartPointer<vtkTubeFilter>tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
//	tubeFilter->SetInputConnection(line->GetOutputPort());
//	tubeFilter->SetRadius(tube_radius);
//	tubeFilter->SetNumberOfSides(100);
//	tubeFilter->CappingOn();
//	tubeFilter->Update();
//	return tubeFilter;
//}
////������PolyData ת��Ϊ������Ƭ����ʽ�ĺ���
//vtkSmartPointer<vtkTriangleFilter> Polydata_to_Triangle(vtkSmartPointer<vtkPolyData> poly_data) {
//	vtkSmartPointer<vtkTriangleFilter> triFilter_source =
//		vtkSmartPointer<vtkTriangleFilter>::New();
//	triFilter_source->SetInputData(poly_data);
//	triFilter_source->Update();
//	return triFilter_source;
//}
////��һϵ�е����������ת��������Ƭ�Ĳ�������ʽ
//vector<vtkSmartPointer<vtkTriangleFilter>> Polydata_to_Triangle_vector(vector<vtkSmartPointer<vtkPolyData>> poly_data_vector) {
//	vector<vtkSmartPointer<vtkTriangleFilter>>result;
//	for (vtkSmartPointer<vtkPolyData> temp : poly_data_vector) {
//		vtkSmartPointer<vtkTriangleFilter> triFilter_source =
//			vtkSmartPointer<vtkTriangleFilter>::New();
//		triFilter_source->SetInputData(temp);
//		triFilter_source->Update();
//		result.push_back(triFilter_source);
//	}
//	return result;
//}
////���������ĺ�����չ�ֵ���ʽ
//vtkSmartPointer<vtkBooleanOperationPolyDataFilter> Triangle_bool_operation(vtkSmartPointer<vtkTriangleFilter>source_1, vtkSmartPointer<vtkTriangleFilter>source_2) {
//	vtkSmartPointer<vtkBooleanOperationPolyDataFilter> shapes = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
//	//���Խ��ж����������Ĳ����Ĳ�������ʽ
//	//shapes->SetInputConnection(0, source_1->GetOutputPort());
//	shapes->SetInputData(0, source_1->GetOutput());
//	shapes->SetInputData(1, source_2->GetOutput());
//	//shapes->SetInputConnection(1, source_2->GetOutputPort());
//	//�����ļ������ʽ
//	shapes->SetOperationToIntersection();
//	//��ļ������ʽ
//	//shapes->SetOperationToDifference();
//	// �����ļ������ʽ
//	//two_shape->SetOperationToUnion();
//	shapes->Update();
//	return shapes;
//}
////��һϵ�е���������в�����������ʽ
//vtkSmartPointer<vtkBooleanOperationPolyDataFilter> Triangle_bool_operation_vector(vector<vtkSmartPointer<vtkTriangleFilter>>vector_source) {
//	vtkSmartPointer<vtkBooleanOperationPolyDataFilter> shapes = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
//	//���Խ��ж����������Ĳ����Ĳ�������ʽ
//	shapes->SetInputConnection(0, vector_source[0]->GetOutputPort());
//	//for (int i = 0; i < vector_source.size(); i++)
//		shapes->SetInputConnection(1, vector_source[1]->GetOutputPort());
//	//shapes->SetInputConnection(1, source_2->GetOutputPort());
//	//�����ļ������ʽ
//	shapes->SetOperationToIntersection();
//	//��ļ������ʽ
//	//shapes->SetOperationToDifference();
//	// �����ļ������ʽ
//	//two_shape->SetOperationToUnion();
//	shapes->Update();
//	return shapes;
//}
////�������������������BooleanOperationPolyDataFilter��������ļ������ʽ
//
//double get_vol_bool(vtkSmartPointer<vtkBooleanOperationPolyDataFilter>bool_polydata) {
//	vtkSmartPointer<vtkMassProperties> massProp_two_shape =
//		vtkSmartPointer<vtkMassProperties>::New();
//	massProp_two_shape->SetInputData(bool_polydata->GetOutput());
//	double vol_shape = massProp_two_shape->GetVolume();
//
//	//float area_shape = massProp_two_shape->GetSurfaceArea();
//	//float maxArea_shape = massProp_two_shape->GetMaxCellArea();
//	//float minArea_shape = massProp_two_shape->GetMinCellArea();
//	//��һ����ƫ�������ģ��뾶Խ����Խ��
//	/*cout << "vol:" << vol_shape << endl;
//	cout << "area:" << area_shape << endl;
//	cout << "maxArea:" << maxArea_shape << endl;
//	cout << "minArea:" << minArea_shape << endl;*/
//	return vol_shape;
//}
////��������Ƭ��ʽ���������ͽ���������㺯��
//double get_volume(vtkSmartPointer<vtkTriangleFilter> triFilter_source) {
//	vtkSmartPointer<vtkMassProperties>massProp_source =
//		vtkSmartPointer<vtkMassProperties>::New();
//	massProp_source->SetInputData(triFilter_source->GetOutput());
//	double vol = massProp_source->GetVolume();
//	return vol;
//}
//
//
//int main()
//{
//	//��һ���ӽ�
//	vtkSmartPointer<vtkRenderer> first_renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	//�ڶ������ӽ�
//	vtkSmartPointer<vtkRenderer> second_renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renWin =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renWin->AddRenderer(first_renderer);
//	renWin->AddRenderer(second_renderer);
//	renWin->SetSize(2000, 1000);
//	vtkSmartPointer<vtkRenderWindowInteractor> iren =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();//���û�ͼ���ڽ���
//	iren->SetRenderWindow(renWin);
//
//	//����������ɵĲ���
//	vtkSmartPointer<vtkParametricEllipsoid> ellipsoid = vtkSmartPointer<vtkParametricEllipsoid>::New();
//	ellipsoid->SetXRadius(Ellipsoid_long);
//	ellipsoid->SetYRadius(Ellipsoid_short);
//	ellipsoid->SetZRadius(Ellipsoid_short);
//	double left_ellipsoid[3] = { -Ellipsoid_long,0,0 };
//	double mid_ellipsoid[3] = { 0,0,0 };
//	double right_ellipsoid[3] = { Ellipsoid_long,0,0 };
//	vtkSmartPointer<vtkParametricFunctionSource> ellipsoid_source = vtkSmartPointer<vtkParametricFunctionSource>::New();
//	ellipsoid_source->SetParametricFunction(ellipsoid);
//	ellipsoid_source->Update();
//	vtkSmartPointer<vtkVertexGlyphFilter> ellipsoid_filter =
//		Insert_point(left_ellipsoid, mid_ellipsoid, right_ellipsoid);
//
//	//������STL �������ʽ
//	vtkSmartPointer<vtkPolyData>poly_tumor = Reader_STL("D:\\\Data Disk\\LW\\Data\\3Dircadb1.1\\STL\\", "livertumor01.stl");
//	//���������Ķ���ĸ�������ʽ
//	vtkIdType nums = poly_tumor->GetNumberOfPoints();
//	cout << "The nums of point :" << nums << endl;
//	vtkSmartPointer<vtkActor> tumor = Build_Actor(poly_tumor);
//	double center[3] = { 0 };
//	Get_center(poly_tumor, center);
//	cout << "The center[0]:" << center[0] << '\t' << center[1] << '\t' << center[2] << endl;
//	//ͳһ����Ϊǰ��60
//	double center_right[3] = { center[0] + 60,center[1]  ,center[2] };
//	double center_left[3] = { center[0] - 60,center[1]  ,center[2] };
//	tumor->GetProperty()->SetColor(0, 1, 1);
//	tumor->GetProperty()->SetOpacity(0.5);
//
//
//	//һ��ֱ��(�����ĺ����ĵ�)�����ݼ��Ľ�������
//	vtkSmartPointer<vtkPoints>left_points = Intersect_points(poly_tumor, center_left, center);
//	cout << "The number:" << left_points->GetNumberOfPoints() << endl;
//	double* left_point = left_points->GetPoint(0);
//
//	//һ��ֱ��(�����ĺ��Ҳ�ĵ�)�����ݼ��Ľ�������
//	vtkSmartPointer<vtkPoints>right_points = Intersect_points(poly_tumor, center_right, center);
//	cout << "The number:" << right_points->GetNumberOfPoints() << endl;
//	double* right_point = right_points->GetPoint(0);
//
//	//��ȡ��λ����
//	double normal_vector[3] = { 0,0,0 };
//	GetNormVector(normal_vector, left_point, right_point);  //�õ���λ����
//	cout << "left:" << left_point[0] << '\t' << left_point[1] << '\t' << left_point[2] << endl;
//	cout << "right:" << right_point[0] << '\t' << right_point[1] << '\t' << right_point[2] << endl;
//	cout << "normal:" << normal_vector[0] << '\t' << normal_vector[1] << '\t' << normal_vector[2] << endl;
//	//����ת�����left mid right ��λ�õ���Ϣ���
//	double trans_left[3] = { center[0] - Ellipsoid_long * normal_vector[0],center[1] - Ellipsoid_long * normal_vector[1],center[2] - Ellipsoid_long * normal_vector[2] };
//
//	double trans_right[3] = { center[0] + Ellipsoid_long * normal_vector[0],center[1] + Ellipsoid_long * normal_vector[1],center[2] + Ellipsoid_long * normal_vector[2] };
//
//	vtkSmartPointer<vtkVertexGlyphFilter>trans_filter =
//		Insert_point(trans_left, center, trans_right);
//
//	/*vtkSmartPointer<vtkTransformPolyDataFilter> transformed_filter =
//		get_transfilter(ellipsoid_source, ellipsoid_filter, trans_filter);*/
//	vector<vtkSmartPointer<vtkTransformPolyDataFilter>> result = get_transfilter_vector(ellipsoid_source, ellipsoid_filter, trans_filter, left_point, center, right_point, normal_vector);
//	//vtkSmartPointer<vtkPolyData> trans_polydata = trans_polydata_Ellipsoid(result[100]);
//	vector<vtkSmartPointer<vtkPolyData>> trans_polydatas = trans_polydata_Ellipsoid_vector(result);
//	//vtkSmartPointer<vtkActor> transformedActor = Build_Actor(trans_polydatas[0]);
//	vector<vtkSmartPointer<vtkTriangleFilter>> vector_Triangle = Polydata_to_Triangle_vector(trans_polydatas);
//	vtkSmartPointer<vtkBooleanOperationPolyDataFilter> vector_BooleanOperation = Triangle_bool_operation_vector(vector_Triangle);
//	vector<vtkSmartPointer<vtkActor>> transformedActors = Build_Actor_vector(trans_polydatas);
//	//transformedActor->GetProperty()->SetColor(1,0, 0);
//	//transformedActor->GetProperty()->SetOpacity(0.5);
//
//	//Բ��������ɵķ�ʽ
//	int num = 5;
//	double start_tube[3] = { left_point[0] - num * normal_vector[0],left_point[1] - num * normal_vector[1],left_point[2] - num * normal_vector[2] };
//	double end_tube[3] = { right_point[0] + num * normal_vector[0],right_point[1] + num * normal_vector[1],right_point[2] + num * normal_vector[2] };
//	vtkSmartPointer<vtkTubeFilter>tubeFilter = Create_Tube(50.0, start_tube, end_tube);
//	//Բ���������ļ������ʽ
//	vtkSmartPointer<vtkPolyData> tube_polydata = trans_polydata_tube(tubeFilter);
//	vtkSmartPointer<vtkTriangleFilter> tube_triangle = Polydata_to_Triangle(tube_polydata);
//	double vol_tube = get_volume(tube_triangle);
//	cout << "The original vol:" << vol_tube << endl;
//	vtkSmartPointer<vtkActor> tube_actor = Build_Actor(tube_polydata);
//	tube_actor->GetProperty()->SetColor(0.5, 0.5, 1);
//	tube_actor->GetProperty()->SetOpacity(0.5);
//
//	//���������������
//	vtkSmartPointer<vtkTriangleFilter> tumor_triangle = Polydata_to_Triangle(poly_tumor);
//	double vol_tumor = get_volume(tumor_triangle);
//	cout << "The tumor_vol :" << vol_tumor << endl;
//
//	//���Բ���������Ľ��������
//	vtkSmartPointer<vtkBooleanOperationPolyDataFilter> tumor_and_tube = Triangle_bool_operation(tube_triangle, tumor_triangle);
//
//	//double tumor_and_tube_vol = get_vol_bool(tumor_and_tube);
//	//cout << "The tumor_and_tube_vol :" << tumor_and_tube_vol << endl;
//
//	//Boolpolydata����ת��Ϊ polydata�ĸ�ʽ��
//	vtkSmartPointer<vtkPolyData>Bool_Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//	Bool_Poly_Data->DeepCopy(tumor_and_tube->GetOutput());
//	vtkSmartPointer<vtkActor>Bool_actor = Build_Actor(Bool_Poly_Data);
//	Bool_actor->GetProperty()->SetColor(0, 0, 1);
//	Bool_actor->GetProperty()->SetOpacity(0.5);
//	//���ò�ͬ����ͼ����ʽ
//	double first_viewpoint[4] = { 0.0,0.0,0.5,1.0 };
//	double second_viewpoint[4] = { 0.5,0.0,1.0,1.0 };
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
//	//��һ����ͼ�Ļ���
//	first_renderer->AddActor(tumor);
//	first_renderer->AddActor(tube_actor);
//	for (vtkSmartPointer<vtkActor> transformedActor : transformedActors)
//		first_renderer->AddActor(transformedActor);
//	first_renderer->SetActiveCamera(aCamera);//�������
//	first_renderer->ResetCamera();//�������
//	first_renderer->SetBackground(0, 1, 0);
//	first_renderer->ResetCameraClippingRange();
//	first_renderer->SetViewport(first_viewpoint);
//	//�ڶ�����ͼ�Ļ���
//	second_renderer->AddActor(Bool_actor);//���������ͶӰ��
//	second_renderer->ResetCamera();//�������
//	second_renderer->SetBackground(1, 0, 0);//���ñ�����ɫ��double����
//	second_renderer->ResetCameraClippingRange();
//	second_renderer->SetViewport(second_viewpoint);
//	renWin->Render();
//	iren->Initialize();
//	iren->Start();
//	return 0;
//}
//
