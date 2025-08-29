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
////获得单位向量
//void GetNormVector(double* Direction, double* p1, double* p2)
//{
//	//计算两个点的距离得出的结果是 距离的平方的形式
//	double length = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
//	cout << "The length:" << length << endl;
//	Direction[0] = (p2[0] - p1[0]) / length;
//	Direction[1] = (p2[1] - p1[1]) / length;
//	Direction[2] = (p2[2] - p1[2]) / length;
//}
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
////椭球体的转化
//vtkSmartPointer<vtkPolyData> trans_polydata_Ellipsoid(vtkSmartPointer<vtkTransformPolyDataFilter> transformed_filter) {
//	vtkSmartPointer<vtkPolyData>Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//	Poly_Data->DeepCopy(transformed_filter->GetOutput());
//	return Poly_Data;
//}
////椭球体的序列的转化
//vector<vtkSmartPointer<vtkPolyData>> trans_polydata_Ellipsoid_vector(vector<vtkSmartPointer<vtkTransformPolyDataFilter>> transformed_filters) {
//	vector<vtkSmartPointer<vtkPolyData>> result;
//	for (vtkSmartPointer<vtkTransformPolyDataFilter>transformed_filter : transformed_filters) {
//		vtkSmartPointer<vtkPolyData>Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//		Poly_Data->DeepCopy(transformed_filter->GetOutput());
//		result.push_back(Poly_Data);
//	}
//	return result;
//}
////圆柱体的转化
//vtkSmartPointer<vtkPolyData> trans_polydata_tube(vtkSmartPointer<vtkTubeFilter>tubeFilter) {
//	vtkSmartPointer<vtkPolyData>Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//	Poly_Data->DeepCopy(tubeFilter->GetOutput());
//	return Poly_Data;
//}
////创建actor 的函数(通用)
//vtkSmartPointer<vtkActor> Build_Actor(vtkSmartPointer<vtkPolyData> PolyData) {
//	vtkSmartPointer<vtkPolyDataMapper> Mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	Mapper->SetInputData(PolyData);//设置输入链接，输出端口获取来自vtkReader的数据	
//	Mapper->ScalarVisibilityOff();    //这样不会带颜色
//	vtkSmartPointer<vtkActor> Actor =
//		vtkSmartPointer<vtkActor>::New();
//	Actor->SetMapper(Mapper);//设置映射器
//	return Actor;
//}
////创建actor 序列的函数(通用的版本适用于椭球体以及 圆柱体的情况)
//vector<vtkSmartPointer<vtkActor>> Build_Actor_vector(vector<vtkSmartPointer<vtkPolyData>> PolyDatas) {
//	vector<vtkSmartPointer<vtkActor>> Actors;
//	for (vtkSmartPointer<vtkPolyData > PolyData : PolyDatas) {
//		vtkSmartPointer<vtkPolyDataMapper> Mapper =
//			vtkSmartPointer<vtkPolyDataMapper>::New();
//		Mapper->SetInputData(PolyData);//设置输入链接，输出端口获取来自vtkReader的数据	
//		Mapper->ScalarVisibilityOff();    //这样不会带颜色
//		vtkSmartPointer<vtkActor> Actor =
//			vtkSmartPointer<vtkActor>::New();
//		Actor->SetMapper(Mapper);//设置映射器
//		Actor->GetProperty()->SetColor(1, 0, 0);
//		Actor->GetProperty()->SetOpacity(0.5);
//		Actors.push_back(Actor);
//	}
//	return Actors;
//}
////寻找 STL 模型的中心的(肿瘤的实现)
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
////一条直线和一个数据集的交点的求解的情况
//vtkSmartPointer<vtkPoints> Intersect_points(vtkSmartPointer<vtkPolyData>poly_data, double* start_point, double* end_point) {
//	vtkSmartPointer<vtkOBBTree> obbtree = vtkSmartPointer<vtkOBBTree>::New();
//	obbtree->SetDataSet(poly_data);
//	obbtree->BuildLocator();
//	vtkSmartPointer<vtkPoints>intersect_point = vtkSmartPointer<vtkPoints>::New();
//	obbtree->IntersectWithLine(start_point, end_point, intersect_point, NULL);
//	return intersect_point;
//}
////3个定位的数据点的数据格式的转化函数
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
////求解变换矩阵，然后转化椭球体的位置情况
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
//	//这里面有一个输出的操作
//	cout << "The result :" << *m << endl;
//	vtkSmartPointer<vtkTransform> trans = vtkSmartPointer<vtkTransform>::New();
//	trans->SetMatrix(m);
//	//椭球体的转化的关系的求解
//	vtkSmartPointer<vtkTransformPolyDataFilter> transformfilter =
//		vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//	transformfilter->SetInputConnection(source->GetOutputPort());
//	transformfilter->SetTransform(trans);
//	transformfilter->Update();
//	return transformfilter;
//}
////沿着直线生成一系列的椭球体
//vector<vtkSmartPointer<vtkTransformPolyDataFilter>> get_transfilter_vector(vtkSmartPointer<vtkParametricFunctionSource> source, vtkSmartPointer<vtkVertexGlyphFilter>source_filter, vtkSmartPointer<vtkVertexGlyphFilter>target_filter, double* left, double* center, double* right, double* normal) {
//	vector<vtkSmartPointer<vtkTransformPolyDataFilter>> result;
//	//得到变换关系的形式的展示
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
//	//这里面有一个输出的操作
//	cout << "The result :" << *m << endl;
//	vtkSmartPointer<vtkTransform> trans = vtkSmartPointer<vtkTransform>::New();
//	trans->SetMatrix(m);
//	//椭球体的转化的关系的求解
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
////圆柱体的生成的操作
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
////将数据PolyData 转换为三角面片的形式的函数
//vtkSmartPointer<vtkTriangleFilter> Polydata_to_Triangle(vtkSmartPointer<vtkPolyData> poly_data) {
//	vtkSmartPointer<vtkTriangleFilter> triFilter_source =
//		vtkSmartPointer<vtkTriangleFilter>::New();
//	triFilter_source->SetInputData(poly_data);
//	triFilter_source->Update();
//	return triFilter_source;
//}
////对一系列的椭球体进行转换三角面片的操作的形式
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
////布尔操作的函数的展现的形式
//vtkSmartPointer<vtkBooleanOperationPolyDataFilter> Triangle_bool_operation(vtkSmartPointer<vtkTriangleFilter>source_1, vtkSmartPointer<vtkTriangleFilter>source_2) {
//	vtkSmartPointer<vtkBooleanOperationPolyDataFilter> shapes = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
//	//可以进行多个操作对象的布尔的操作的形式
//	//shapes->SetInputConnection(0, source_1->GetOutputPort());
//	shapes->SetInputData(0, source_1->GetOutput());
//	shapes->SetInputData(1, source_2->GetOutput());
//	//shapes->SetInputConnection(1, source_2->GetOutputPort());
//	//交集的计算的形式
//	shapes->SetOperationToIntersection();
//	//差集的计算的形式
//	//shapes->SetOperationToDifference();
//	// 并集的计算的形式
//	//two_shape->SetOperationToUnion();
//	shapes->Update();
//	return shapes;
//}
////对一系列的椭球体进行布尔操作的形式
//vtkSmartPointer<vtkBooleanOperationPolyDataFilter> Triangle_bool_operation_vector(vector<vtkSmartPointer<vtkTriangleFilter>>vector_source) {
//	vtkSmartPointer<vtkBooleanOperationPolyDataFilter> shapes = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
//	//可以进行多个操作对象的布尔的操作的形式
//	shapes->SetInputConnection(0, vector_source[0]->GetOutputPort());
//	//for (int i = 0; i < vector_source.size(); i++)
//		shapes->SetInputConnection(1, vector_source[1]->GetOutputPort());
//	//shapes->SetInputConnection(1, source_2->GetOutputPort());
//	//交集的计算的形式
//	shapes->SetOperationToIntersection();
//	//差集的计算的形式
//	//shapes->SetOperationToDifference();
//	// 并集的计算的形式
//	//two_shape->SetOperationToUnion();
//	shapes->Update();
//	return shapes;
//}
////布尔操作后的数据类型BooleanOperationPolyDataFilter进行体积的计算的形式
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
//	//有一定的偏差的情况的，半径越大差距越大
//	/*cout << "vol:" << vol_shape << endl;
//	cout << "area:" << area_shape << endl;
//	cout << "maxArea:" << maxArea_shape << endl;
//	cout << "minArea:" << minArea_shape << endl;*/
//	return vol_shape;
//}
////对三角面片形式的数据类型进行体积计算函数
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
//	//第一个视角
//	vtkSmartPointer<vtkRenderer> first_renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	//第二部分视角
//	vtkSmartPointer<vtkRenderer> second_renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renWin =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renWin->AddRenderer(first_renderer);
//	renWin->AddRenderer(second_renderer);
//	renWin->SetSize(2000, 1000);
//	vtkSmartPointer<vtkRenderWindowInteractor> iren =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();//设置绘图窗口交互
//	iren->SetRenderWindow(renWin);
//
//	//椭球体的生成的操作
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
//	//肿瘤的STL 导入的形式
//	vtkSmartPointer<vtkPolyData>poly_tumor = Reader_STL("D:\\\Data Disk\\LW\\Data\\3Dircadb1.1\\STL\\", "livertumor01.stl");
//	//测试肿瘤的顶点的个数的形式
//	vtkIdType nums = poly_tumor->GetNumberOfPoints();
//	cout << "The nums of point :" << nums << endl;
//	vtkSmartPointer<vtkActor> tumor = Build_Actor(poly_tumor);
//	double center[3] = { 0 };
//	Get_center(poly_tumor, center);
//	cout << "The center[0]:" << center[0] << '\t' << center[1] << '\t' << center[2] << endl;
//	//统一设置为前后60
//	double center_right[3] = { center[0] + 60,center[1]  ,center[2] };
//	double center_left[3] = { center[0] - 60,center[1]  ,center[2] };
//	tumor->GetProperty()->SetColor(0, 1, 1);
//	tumor->GetProperty()->SetOpacity(0.5);
//
//
//	//一根直线(已质心和左侧的点)和数据集的交点的情况
//	vtkSmartPointer<vtkPoints>left_points = Intersect_points(poly_tumor, center_left, center);
//	cout << "The number:" << left_points->GetNumberOfPoints() << endl;
//	double* left_point = left_points->GetPoint(0);
//
//	//一根直线(已质心和右侧的点)和数据集的交点的情况
//	vtkSmartPointer<vtkPoints>right_points = Intersect_points(poly_tumor, center_right, center);
//	cout << "The number:" << right_points->GetNumberOfPoints() << endl;
//	double* right_point = right_points->GetPoint(0);
//
//	//获取单位向量
//	double normal_vector[3] = { 0,0,0 };
//	GetNormVector(normal_vector, left_point, right_point);  //得到单位向量
//	cout << "left:" << left_point[0] << '\t' << left_point[1] << '\t' << left_point[2] << endl;
//	cout << "right:" << right_point[0] << '\t' << right_point[1] << '\t' << right_point[2] << endl;
//	cout << "normal:" << normal_vector[0] << '\t' << normal_vector[1] << '\t' << normal_vector[2] << endl;
//	//设置转化后的left mid right 的位置的信息情况
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
//	//圆柱体的生成的方式
//	int num = 5;
//	double start_tube[3] = { left_point[0] - num * normal_vector[0],left_point[1] - num * normal_vector[1],left_point[2] - num * normal_vector[2] };
//	double end_tube[3] = { right_point[0] + num * normal_vector[0],right_point[1] + num * normal_vector[1],right_point[2] + num * normal_vector[2] };
//	vtkSmartPointer<vtkTubeFilter>tubeFilter = Create_Tube(50.0, start_tube, end_tube);
//	//圆柱体的体积的计算的形式
//	vtkSmartPointer<vtkPolyData> tube_polydata = trans_polydata_tube(tubeFilter);
//	vtkSmartPointer<vtkTriangleFilter> tube_triangle = Polydata_to_Triangle(tube_polydata);
//	double vol_tube = get_volume(tube_triangle);
//	cout << "The original vol:" << vol_tube << endl;
//	vtkSmartPointer<vtkActor> tube_actor = Build_Actor(tube_polydata);
//	tube_actor->GetProperty()->SetColor(0.5, 0.5, 1);
//	tube_actor->GetProperty()->SetOpacity(0.5);
//
//	//求解的是肿瘤的体积
//	vtkSmartPointer<vtkTriangleFilter> tumor_triangle = Polydata_to_Triangle(poly_tumor);
//	double vol_tumor = get_volume(tumor_triangle);
//	cout << "The tumor_vol :" << vol_tumor << endl;
//
//	//求解圆柱与肿瘤的交集的情况
//	vtkSmartPointer<vtkBooleanOperationPolyDataFilter> tumor_and_tube = Triangle_bool_operation(tube_triangle, tumor_triangle);
//
//	//double tumor_and_tube_vol = get_vol_bool(tumor_and_tube);
//	//cout << "The tumor_and_tube_vol :" << tumor_and_tube_vol << endl;
//
//	//Boolpolydata类型转换为 polydata的格式的
//	vtkSmartPointer<vtkPolyData>Bool_Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//	Bool_Poly_Data->DeepCopy(tumor_and_tube->GetOutput());
//	vtkSmartPointer<vtkActor>Bool_actor = Build_Actor(Bool_Poly_Data);
//	Bool_actor->GetProperty()->SetColor(0, 0, 1);
//	Bool_actor->GetProperty()->SetOpacity(0.5);
//	//设置不同的视图的形式
//	double first_viewpoint[4] = { 0.0,0.0,0.5,1.0 };
//	double second_viewpoint[4] = { 0.5,0.0,1.0,1.0 };
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
//	//第一个视图的绘制
//	first_renderer->AddActor(tumor);
//	first_renderer->AddActor(tube_actor);
//	for (vtkSmartPointer<vtkActor> transformedActor : transformedActors)
//		first_renderer->AddActor(transformedActor);
//	first_renderer->SetActiveCamera(aCamera);//设置相机
//	first_renderer->ResetCamera();//重置相机
//	first_renderer->SetBackground(0, 1, 0);
//	first_renderer->ResetCameraClippingRange();
//	first_renderer->SetViewport(first_viewpoint);
//	//第二个视图的绘制
//	second_renderer->AddActor(Bool_actor);//添加肿瘤的投影点
//	second_renderer->ResetCamera();//重置相机
//	second_renderer->SetBackground(1, 0, 0);//设置背景颜色，double类型
//	second_renderer->ResetCameraClippingRange();
//	second_renderer->SetViewport(second_viewpoint);
//	renWin->Render();
//	iren->Initialize();
//	iren->Start();
//	return 0;
//}
//
