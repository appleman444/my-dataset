//#include <vtkAutoInit.h>
//VTK_MODULE_INIT(vtkRenderingOpenGL2);
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingFreeType);
//
//#include "vtkDICOMImageReader.h"
//#include "vtkPolyDataWriter.h"
//#include "vtkPolyDataReader.h"
//#include <vtkOBBTree.h>
//#include <vtkBooleanOperationPolyDataFilter.h>
//#include <vtkTubeFilter.h>
//#include <vtkLandmarkTransform.h>
//#include <vtkSmartPointer.h>
//#include <vtkTransform.h>
//#include <vtkMatrix4x4.h>
//#include <vtkJPEGReader.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkInteractorStyleImage.h>
//#include <vtkRenderer.h>
//#include <vtkStringArray.h>
//#include <vtkPlaneSource.h>
//#include <vtkRenderWindow.h>
//#include <vtkContourFilter.h>
//#include <vtkPolyDataNormals.h>
//#include <vtkStripper.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkVertexGlyphFilter.h>
//#include <vtkTextActor.h>
//#include <vtkTriangleFilter.h>
//#include <vtkTextProperty.h>
//#include <vtkAxesActor.h>
//#include <vtkIterativeClosestPointTransform.h>
//#include <vtkCamera.h>
//#include <vtkSphereSource.h>
//#include <vtkSTLWriter.h>
//#include <vtkRegularPolygonSource.h>
//#include <vtkSTLReader.h>
//#include <vtkSmoothPolyDataFilter.h>
//#include <vtkMassProperties.h>
//#include <vtkFloatArray.h>
//#include <vtkLookupTable.h>
//#include <vtkPointData.h>
//#include <vtkTransformPolyDataFilter.h>
//#include <vtkPolyData.h>
//#include <vtkProperty.h>
//#include <vtkSelectEnclosedPoints.h>
//#include <vtkPlane.h>
//#include <vtkNew.h>
//#include <vtkPoints.h>
//#include <vtkLine.h>
//#include <vtkMath.h>
//#include <vector>
//#include <vtkLineSource.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkParametricEllipsoid.h>
//#include <vtkParametricFunctionSource.h>
//#include <random>
//#include <iostream>
//#include <cmath>
//#include <ctime>
//#include <cstdlib>
//#include <bitset>
//#include <iomanip>
//#include <fstream>
//#include <vtkBoundingBox.h>
//#include <vtkCubeSource.h>
//#include <vtkDelaunay3D.h>
//#include <vtkCellArray.h>
//#include <vtkDataSetSurfaceFilter.h>
//#include <vtkSmartPointer.h>
//#include <vtkPoints.h>
//#include <vtkPolyData.h>
//#include <vtkVertexGlyphFilter.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <array>
//#include <string>
//#include <cstdio> // for fgets
//#include <vtkTable.h>
//#include <vtk-9.4/vtkDelimitedTextReader.h>
//#include <vtkSmartPointer.h>
//#include <vtkPoints.h>
//#include <vtkPolyData.h>
//#include <vtkSphereSource.h>
//#include <vtkGlyph3D.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkProperty.h>
//
//using namespace std;
//const double PI = 3.141592653589793;
//
//
//// VTK 文件读取
//vtkSmartPointer<vtkPolyData> Reader_VTK(const char* path_name, const char* file_name) {
//    if (!path_name || !file_name) {
//        cerr << "Error: Invalid file path!" << endl;
//        return vtkSmartPointer<vtkPolyData>::New();
//    }
//    string full_path = string(path_name) + string(file_name);
//    vtkSmartPointer<vtkPolyDataReader> vtk_Reader = vtkSmartPointer<vtkPolyDataReader>::New();
//    vtk_Reader->SetFileName(full_path.c_str());
//    vtk_Reader->Update();
//    vtkSmartPointer<vtkPolyData> Poly_Data = vtkSmartPointer<vtkPolyData>::New();
//    Poly_Data->DeepCopy(vtk_Reader->GetOutput());
//    return Poly_Data;
//}
//
//// 计算三角形法线
//void ComputeTriangleNormal(vtkSmartPointer<vtkPoints> points, const vtkIdType* ptIds, double* normal) {
//    double p1[3], p2[3], p3[3];
//    points->GetPoint(ptIds[0], p1);
//    points->GetPoint(ptIds[1], p2);
//    points->GetPoint(ptIds[2], p3);
//
//    double v1[3] = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
//    double v2[3] = { p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2] };
//
//    vtkMath::Cross(v1, v2, normal);
//    vtkMath::Normalize(normal);
//}
//
//// 计算两个向量的夹角
//double ComputeAngleBetweenVectors(double* v1, double* v2) {
//    double dotProduct = vtkMath::Dot(v1, v2);
//    double magnitudeV1 = vtkMath::Norm(v1);
//    double magnitudeV2 = vtkMath::Norm(v2);
//
//    return acos(dotProduct / (magnitudeV1 * magnitudeV2)) * 180.0 / PI;  // 返回角度，单位为度
//}
//
//// 计算肿瘤表面点与肝包膜的最小距离
//bool CheckMinDistanceToSkin(vtkSmartPointer<vtkPolyData> skin, double* tumorCenter, vtkSmartPointer<vtkPoints> points, const vtkIdType* ptIds) {
//    double minDistance = 1e6;  // 设置一个极大的初始最小距离
//
//    // 获取三角形的三个顶点
//    double p1[3], p2[3], p3[3];
//    points->GetPoint(ptIds[0], p1);
//    points->GetPoint(ptIds[1], p2);
//    points->GetPoint(ptIds[2], p3);
//
//    // 计算三角形的法向量
//    double normal[3] = { 0.0, 0.0, 0.0 };
//    double v1[3] = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
//    double v2[3] = { p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2] };
//    vtkMath::Cross(v1, v2, normal);
//    vtkMath::Normalize(normal);
//
//    // 计算肿瘤质心到三角形平面的距离
//    double vectorToPoint[3] = { tumorCenter[0] - p1[0], tumorCenter[1] - p1[1], tumorCenter[2] - p1[2] };
//    double distance = vtkMath::Dot(vectorToPoint, normal);
//    distance = fabs(distance); // 取绝对值，得到与平面的垂直距离
//
//    if (distance < minDistance) {
//        minDistance = distance;  // 更新最小距离
//    }
//
//    // 如果最小距离小于预设的安全距离阈值，则返回false
//    if (minDistance < 5.0) { // 假设安全距离小于5mm
//        return false;
//    }
//
//    return true;
//}
//
//// 判断消融针路径是否满足约束
//bool CheckConstraints(vtkSmartPointer<vtkPolyData> skin, double* tumorCenter, vtkSmartPointer<vtkOBBTree> obbTree, const vtkIdType* ptIds, vtkSmartPointer<vtkPoints> points) {
//    double avgDistance = 0.0;
//    vtkSmartPointer<vtkPoints> intersectPoints = vtkSmartPointer<vtkPoints>::New();
//
//    // 1. 安全距离约束
//    for (int i = 0; i < 3; i++) {
//        double point[3];
//        points->GetPoint(ptIds[i], point);
//
//        // 计算肿瘤表面点与肝包膜的最小距离
//        if (obbTree->IntersectWithLine(tumorCenter, point, intersectPoints, nullptr)) {
//            return false; // 路径与肝包膜有交点，违反了安全距离
//        }
//
//        double distance = sqrt(pow(point[0] - tumorCenter[0], 2) +
//            pow(point[1] - tumorCenter[1], 2) +
//            pow(point[2] - tumorCenter[2], 2));
//        avgDistance += distance / 3.0;
//    }
//
//    // 安全距离约束: 检查肿瘤表面到路径的平均距离是否小于5mm或大于60mm
//    if (avgDistance < 5.0 || avgDistance > 60.0) {
//        return false;
//    }
//
//    // 2. 切线角度约束
//    // 计算交点附近三角形的法线
//    double normal[3] = { 0.0, 0.0, 0.0 };
//    for (int i = 0; i < 3; i++) {
//        double point[3];
//        points->GetPoint(ptIds[i], point);
//        // 获取交点附近三角形的法线
//        ComputeTriangleNormal(skin->GetPoints(), ptIds, normal);
//
//        // 计算消融针路径与肝包膜表面法线之间的夹角
//        double direction[3] = { point[0] - tumorCenter[0], point[1] - tumorCenter[1], point[2] - tumorCenter[2] };
//        double angle = ComputeAngleBetweenVectors(direction, normal);
//
//        if (angle < 20.0) { // 切线角度小于20°，不满足约束
//            return false;
//        }
//    }
//
//    return true;
//}
//
//
//vtkSmartPointer<vtkPolyData> FilterValidSkinTriangles(vtkSmartPointer<vtkPolyData> skin,
//    vtkSmartPointer<vtkPolyData> bone,
//    vtkSmartPointer<vtkPolyData> artery,
//    vtkSmartPointer<vtkPolyData> venoussystem,
//    double* tumorCenter) {
//    // 将 skin 的数据转换为三角形
//    vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
//    triangleFilter->SetInputData(skin);
//    triangleFilter->Update();
//
//    // 建立骨骼、动脉和静脉的 OBBTree，供后续碰撞检测使用
//    vtkSmartPointer<vtkOBBTree> obbTree = vtkSmartPointer<vtkOBBTree>::New();
//    obbTree->SetDataSet(bone);
//    obbTree->BuildLocator();
//
//    vtkSmartPointer<vtkOBBTree> obbTree1 = vtkSmartPointer<vtkOBBTree>::New();
//    obbTree1->SetDataSet(artery);
//    obbTree1->BuildLocator();
//
//    vtkSmartPointer<vtkOBBTree> obbTree2 = vtkSmartPointer<vtkOBBTree>::New();
//    obbTree2->SetDataSet(venoussystem);
//    obbTree2->BuildLocator();
//
//    // 存储筛选后的三角形
//    vtkSmartPointer<vtkPolyData> filteredTriangles = vtkSmartPointer<vtkPolyData>::New();
//    vtkSmartPointer<vtkCellArray> newCells = vtkSmartPointer<vtkCellArray>::New();
//    vtkSmartPointer<vtkPoints> points = triangleFilter->GetOutput()->GetPoints();
//    vtkSmartPointer<vtkCellArray> polys = triangleFilter->GetOutput()->GetPolys();
//
//    vtkIdType npts;
//    const vtkIdType* ptIds;
//    polys->InitTraversal();
//
//    while (polys->GetNextCell(npts, ptIds)) {
//        if (npts == 3) {
//            // 获取三角形三个顶点的坐标
//            double p0[3], p1[3], p2[3];
//            points->GetPoint(ptIds[0], p0);
//            points->GetPoint(ptIds[1], p1);
//            points->GetPoint(ptIds[2], p2);
//
//            // 计算三角形质心
//            double centroid[3] = { (p0[0] + p1[0] + p2[0]) / 3.0,
//                                   (p0[1] + p1[1] + p2[1]) / 3.0,
//                                   (p0[2] + p1[2] + p2[2]) / 3.0 };
//
//            // 计算肿瘤中心到三角形质心的距离
//            double dist2 = vtkMath::Distance2BetweenPoints(tumorCenter, centroid);
//            // 仅保留距离小于150的三角形
//            if (dist2 < 45.0 * 45.0) {
//                // 检查其他约束条件
//                bool valid = CheckConstraints(skin, tumorCenter, obbTree, ptIds, points);
//                if (valid) {
//                    newCells->InsertNextCell(3, ptIds);
//                }
//            }
//        }
//    }
//    filteredTriangles->SetPoints(points);
//    filteredTriangles->SetPolys(newCells);
//    return filteredTriangles;
//}
//// 假设 Whale 类存储了优化过程中的每个解的位置和适应度
//class Whale {
//public:
//    std::vector<double> position; // 存储鲸鱼位置的向量 (x, y, z)
//    std::vector<double> fitness;  // 存储鲸鱼的适应度值
//};
//
//// 计算消融路径的长度（路径是由两点表示的直线）
//double CalculatePathLength(const double* point1, const double* point2) {
//    return vtkMath::Distance2BetweenPoints(point1, point2);  // 距离的平方
//}
//
//// 计算路径与血管/骨头的最小距离
//double CalculatePathDistanceToStructure(const double* pathStart, const double* pathEnd, vtkSmartPointer<vtkOBBTree> obbTree) {
//    // 计算路径与骨头或血管的最小距离
//    // 使用 OBBTree 查找路径和血管/骨头之间的距离（可以优化路径点到骨头、血管的距离）
//    double minDistance = std::numeric_limits<double>::max();
//    // 你可以通过遍历或查找来获取路径和结构之间的最小距离，假设这里计算最小距离为模拟
//    return minDistance;
//}
//
//// 计算路径与肿瘤质心的平行度
//double CalculatePathAlignmentWithTumorCenter(const double* pathStart, const double* pathEnd, const double* tumorCenter) {
//    // 计算路径的方向向量
//    double pathVector[3] = { pathEnd[0] - pathStart[0], pathEnd[1] - pathStart[1], pathEnd[2] - pathStart[2] };
//
//    // 计算肿瘤质心的方向向量
//    double tumorVector[3] = { tumorCenter[0] - pathStart[0], tumorCenter[1] - pathStart[1], tumorCenter[2] - pathStart[2] };
//
//    // 计算两个向量的点积，点积越大，表示越平行
//    double dotProduct = vtkMath::Dot(pathVector, tumorVector);
//    double pathLength = vtkMath::Norm(pathVector);
//    double tumorLength = vtkMath::Norm(tumorVector);
//
//    // 计算平行度，点积越大，方向越相似
//    double alignment = dotProduct / (pathLength * tumorLength);
//    return alignment;
//}
//
//// 从三角形中采样点（这里选择三角形重心）
//void SamplePointsFromTriangles(vtkSmartPointer<vtkPolyData> filteredSkin, std::vector<std::vector<double>>& sampledPoints) {
//    vtkSmartPointer<vtkPoints> points = filteredSkin->GetPoints();
//    vtkSmartPointer<vtkCellArray> polys = filteredSkin->GetPolys();
//
//    vtkIdType npts;
//    const vtkIdType* ptIds;
//    polys->InitTraversal();
//
//    while (polys->GetNextCell(npts, ptIds)) {
//        if (npts == 3) { // 如果是三角形
//            // 获取三角形三个顶点的坐标
//            double p0[3], p1[3], p2[3];
//            points->GetPoint(ptIds[0], p0);
//            points->GetPoint(ptIds[1], p1);
//            points->GetPoint(ptIds[2], p2);
//
//            // 计算三角形重心
//            double centroid[3] = { (p0[0] + p1[0] + p2[0]) / 3.0,
//                                   (p0[1] + p1[1] + p2[1]) / 3.0,
//                                   (p0[2] + p1[2] + p2[2]) / 3.0 };
//
//            // 存储重心
//            sampledPoints.push_back({ centroid[0], centroid[1], centroid[2] });
//        }
//    }
//}
//
////// 多目标鲸鱼优化算法
////std::vector<std::vector<double>> WhaleOptimizationAlgorithm(vtkSmartPointer<vtkPolyData> filteredSkin, vtkSmartPointer<vtkOBBTree> obbTree, double* tumorCenter) {
////    // 假设优化的种群数量为 10，每个鲸鱼有 3 个位置参数 (x, y, z)
////    int populationSize = 1;
////    std::vector<Whale> population(populationSize);
////    std::vector<std::vector<double>> optimizedPoints;  // 存储优化后的点（起始点）
////
////    // 从筛选后的三角面片中采样点
////    std::vector<std::vector<double>> sampledPoints;
////    SamplePointsFromTriangles(filteredSkin, sampledPoints);
////
////    // 初始化种群位置 (即路径的起始点和终止点)
////    for (int i = 0; i < populationSize; i++) {
////        // 随机选择采样的点作为路径的起始点
////        int startIndex = rand() % sampledPoints.size();
////        int endIndex = rand() % sampledPoints.size();
////
////        // 存储路径的起始点和终止点
////        population[i].position = { sampledPoints[startIndex][0], sampledPoints[startIndex][1], sampledPoints[startIndex][2] };
////        population[i].position.push_back(sampledPoints[endIndex][0]);
////        population[i].position.push_back(sampledPoints[endIndex][1]);
////        population[i].position.push_back(sampledPoints[endIndex][2]);
////    }
////
////    // 对每个鲸鱼进行优化（计算适应度）
////    for (int i = 0; i < populationSize; i++) {
////        double pathStart[3] = { population[i].position[0], population[i].position[1], population[i].position[2] };
////        double pathEnd[3] = { population[i].position[3], population[i].position[4], population[i].position[5] };
////
////        // 计算目标1：路径长度最小化
////        double pathLength = CalculatePathLength(pathStart, pathEnd);
////        population[i].fitness.push_back(pathLength);  // 目标1：路径长度最小化
////
////        // 计算目标2：路径与骨头/血管的距离最大化
////        double pathDistanceToStructure = CalculatePathDistanceToStructure(pathStart, pathEnd, obbTree);
////        population[i].fitness.push_back(-pathDistanceToStructure);  // 目标2：最大化路径与血管骨头的距离
////
////        // 计算目标3：路径平行度最大化
////        double pathAlignment = CalculatePathAlignmentWithTumorCenter(pathStart, pathEnd, tumorCenter);
////        population[i].fitness.push_back(pathAlignment);  // 目标3：最大化路径平行度
////    }
////
////    // 提取路径的起始点到结果数组
////    for (int i = 0; i < populationSize; i++) {
////        double pathStart[3] = { population[i].position[0], population[i].position[1], population[i].position[2] };
////        optimizedPoints.push_back({ pathStart[0], pathStart[1], pathStart[2] });
////    }
////
////    return optimizedPoints;  // 返回优化后的起点
////}
//// 计算两个点之间的距离
//double distance(double x1, double y1, double z1, double x2, double y2, double z2) {
//    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
//}
//
//// 计算单位方向向量
//void normalizeDirection(double& dx, double& dy, double& dz) {
//    double len = sqrt(dx * dx + dy * dy + dz * dz);
//    dx /= len;
//    dy /= len;
//    dz /= len;
//}
//vtkSmartPointer<vtkPolyData> ExtractPointsFromTrianglesAndProcess(vtkSmartPointer<vtkPolyData> skin, double* tumorCenter) {
//    // 将 skin 的数据转换为三角形
//    vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
//    triangleFilter->SetInputData(skin);
//    triangleFilter->Update();
//
//    // 存储筛选后的点
//    vtkSmartPointer<vtkPolyData> filteredPointsPolyData = vtkSmartPointer<vtkPolyData>::New();
//    vtkSmartPointer<vtkPoints> filteredPoints = vtkSmartPointer<vtkPoints>::New();
//
//    vtkSmartPointer<vtkCellArray> polys = triangleFilter->GetOutput()->GetPolys();
//    vtkSmartPointer<vtkPoints> points = triangleFilter->GetOutput()->GetPoints();
//
//    vtkIdType npts;
//    const vtkIdType* ptIds;
//    polys->InitTraversal();
//
//    while (polys->GetNextCell(npts, ptIds)) {
//        if (npts == 3) {
//            // 获取三角形三个顶点的坐标
//            double p0[3], p1[3], p2[3];
//            points->GetPoint(ptIds[0], p0);
//            points->GetPoint(ptIds[1], p1);
//            points->GetPoint(ptIds[2], p2);
//
//            // 计算三角形质心
//            double centroid[3] = { (p0[0] + p1[0] + p2[0]) / 3.0,
//                                   (p0[1] + p1[1] + p2[1]) / 3.0,
//                                   (p0[2] + p1[2] + p2[2]) / 3.0 };
//
//            // 计算肿瘤中心到三角形质心的距离
//            double dist2 = vtkMath::Distance2BetweenPoints(tumorCenter, centroid);
//            // 仅保留距离小于150的三角形
//            if (dist2 < 45.0 * 45.0) {
//                // 将三角形的顶点添加到筛选后的点集
//                filteredPoints->InsertNextPoint(p0);
//                filteredPoints->InsertNextPoint(p1);
//                filteredPoints->InsertNextPoint(p2);
//            }
//        }
//    }
//
//    // 设置点数据
//    filteredPointsPolyData->SetPoints(filteredPoints);
//
//    // 获取并遍历点集
//    vtkSmartPointer<vtkPoints> finalPoints = filteredPointsPolyData->GetPoints();
//
//    // 遍历点集
//    for (vtkIdType i = 0; i < finalPoints->GetNumberOfPoints(); i++) {
//        double pathStart[3];
//        finalPoints->GetPoint(i, pathStart);  // 获取每个点的坐标
//
//        // 使用 pathStart 进行后续操作
//        // 例如，可以在这里使用 pathStart 进行路径规划或其他操作
//        std::cout << "Point " << i << ": (" << pathStart[0] << ", " << pathStart[1] << ", " << pathStart[2] << ")" << std::endl;
//    }
//
//    return filteredPointsPolyData;
//}
//
//int main() {
//    vtkSmartPointer<vtkPolyData> bone = Reader_VTK("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\", "bone.vtk");
//    vtkSmartPointer<vtkPolyData> tumor = Reader_VTK("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\", "livertumor01.vtk");
//    vtkSmartPointer<vtkPolyData> skin = Reader_VTK("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\", "skin.vtk");
//    vtkSmartPointer<vtkPolyData> liver = Reader_VTK("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\", "liver.vtk");
//    vtkSmartPointer<vtkPolyData> artery = Reader_VTK("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\", "artery.vtk");
//    vtkSmartPointer<vtkPolyData> venoussystem = Reader_VTK("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\", "venoussystem.vtk");
//
//    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//    renderWindow->AddRenderer(renderer);
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    renderWindowInteractor->SetRenderWindow(renderWindow);
//
//    double* tumorCenter = tumor->GetCenter();
//    vtkSmartPointer<vtkPolyData> filteredSkin = FilterValidSkinTriangles(skin, bone, artery, venoussystem, tumorCenter);
//    double tumorCenter1[] = { 103.99, 74.89, 116.64 };
//    double tumorCenter2[] = { 102.06, 100.99, 107.93 };
//    double tumorCenter3[] = { 85.23,  87.69,  118.29 };
//    double tumorCenter4[] = { 85.37,  91.77,  92.54 };
//    double tumorCenter5[] = { 70.69, 101.64,  102.01 };
//    double tumorCenter6[] = { 106.74,  78.04,  98.72 };
//    vtkSmartPointer<vtkPolyData> filteredSkin1 = FilterValidSkinTriangles(skin, bone, artery, venoussystem, tumorCenter1);
//    vtkSmartPointer<vtkPolyData> filteredSkin2 = FilterValidSkinTriangles(skin, bone, artery, venoussystem, tumorCenter2);
//    vtkSmartPointer<vtkPolyData> filteredSkin3 = FilterValidSkinTriangles(skin, bone, artery, venoussystem, tumorCenter3);
//    vtkSmartPointer<vtkPolyData> filteredSkin4 = FilterValidSkinTriangles(skin, bone, artery, venoussystem, tumorCenter4);
//    vtkSmartPointer<vtkPolyData> filteredSkin5 = FilterValidSkinTriangles(skin, bone, artery, venoussystem, tumorCenter5);
//    vtkSmartPointer<vtkPolyData> filteredSkin6 = FilterValidSkinTriangles(skin, bone, artery, venoussystem, tumorCenter6);
//    std::string fileName = "D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\CSV\\tumorinpoints\\252020p1\\cluster_1.csv";
//    std::string fileName1 = "D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\CSV\\tumorinpoints\\252020p1\\cluster_2.csv";
//    std::string fileName2 = "D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\CSV\\tumorinpoints\\252020p1\\cluster_3.csv";
//    std::string fileName3 = "D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\CSV\\tumorinpoints\\252020p1\\cluster_4.csv";
//    std::string fileName4 = "D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\CSV\\tumorinpoints\\252020p1\\cluster_5.csv";
//    std::string fileName5 = "D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\CSV\\tumorinpoints\\252020p1\\cluster_6.csv";
//
//    // 创建 vtkDelimitedTextReader 读取 CSV 文件
//    vtkSmartPointer<vtkDelimitedTextReader> reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader->SetFileName(fileName.c_str());
//    reader->SetHaveHeaders(true); // 如果 CSV 文件有头部信息
//    reader->SetFieldDelimiterCharacters(","); // 使用逗号作为分隔符
//    reader->Update();
//    vtkSmartPointer<vtkDelimitedTextReader> reader1 = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader1->SetFileName(fileName1.c_str());
//    reader1->SetHaveHeaders(true); // 如果 CSV 文件有头部信息
//    reader1->SetFieldDelimiterCharacters(","); // 使用逗号作为分隔符
//    reader1->Update();
//    vtkSmartPointer<vtkDelimitedTextReader> reader2 = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader2->SetFileName(fileName2.c_str());
//    reader2->SetHaveHeaders(true); // 如果 CSV 文件有头部信息
//    reader2->SetFieldDelimiterCharacters(","); // 使用逗号作为分隔符
//    reader2->Update();
//    vtkSmartPointer<vtkDelimitedTextReader> reader3 = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader3->SetFileName(fileName3.c_str());
//    reader3->SetHaveHeaders(true); // 如果 CSV 文件有头部信息
//    reader3->SetFieldDelimiterCharacters(","); // 使用逗号作为分隔符
//    reader3->Update();
//    vtkSmartPointer<vtkDelimitedTextReader> reader4 = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader4->SetFileName(fileName4.c_str());
//    reader4->SetHaveHeaders(true); // 如果 CSV 文件有头部信息
//    reader4->SetFieldDelimiterCharacters(","); // 使用逗号作为分隔符
//    reader4->Update();
//    vtkSmartPointer<vtkDelimitedTextReader> reader5 = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader5->SetFileName(fileName5.c_str());
//    reader5->SetHaveHeaders(true); // 如果 CSV 文件有头部信息
//    reader5->SetFieldDelimiterCharacters(","); // 使用逗号作为分隔符
//    reader5->Update();
//
//    // 读取的数据将存储在 vtkTable 中
//    vtkSmartPointer<vtkTable> table = reader->GetOutput();
//    vtkSmartPointer<vtkTable> table1 = reader1->GetOutput();
//    vtkSmartPointer<vtkTable> table2 = reader2->GetOutput();
//    vtkSmartPointer<vtkTable> table3 = reader3->GetOutput();
//    vtkSmartPointer<vtkTable> table4 = reader4->GetOutput();
//    vtkSmartPointer<vtkTable> table5 = reader5->GetOutput();
//
//
//    // 创建 vtkPoints 来存储 CSV 文件中的点数据
//    vtkSmartPointer<vtkPoints> cupoints1 = vtkSmartPointer<vtkPoints>::New();
//    // 假设 CSV 文件的列是 x, y, z
//    for (vtkIdType i = 0; i < table->GetNumberOfRows(); ++i) {
//        double x = table->GetValue(i, 0).ToDouble(); // 获取 x 坐标
//        double y = table->GetValue(i, 1).ToDouble(); // 获取 y 坐标
//        double z = table->GetValue(i, 2).ToDouble(); // 获取 z 坐标
//
//        // 将点添加到 vtkPoints 中
//        cupoints1->InsertNextPoint(x, y, z);
//    }
//    // 创建 vtkPoints 来存储 CSV 文件中的点数据
//    vtkSmartPointer<vtkPoints> cupoints2 = vtkSmartPointer<vtkPoints>::New();
//    // 假设 CSV 文件的列是 x, y, z
//    for (vtkIdType i = 0; i < table1->GetNumberOfRows(); ++i) {
//        double x = table1->GetValue(i, 0).ToDouble(); // 获取 x 坐标
//        double y = table1->GetValue(i, 1).ToDouble(); // 获取 y 坐标
//        double z = table1->GetValue(i, 2).ToDouble(); // 获取 z 坐标
//
//        // 将点添加到 vtkPoints 中
//        cupoints2->InsertNextPoint(x, y, z);
//    }
//    // 创建 vtkPoints 来存储 CSV 文件中的点数据
//    vtkSmartPointer<vtkPoints> cupoints3 = vtkSmartPointer<vtkPoints>::New();
//    // 假设 CSV 文件的列是 x, y, z
//    for (vtkIdType i = 0; i < table2->GetNumberOfRows(); ++i) {
//        double x = table2->GetValue(i, 0).ToDouble(); // 获取 x 坐标
//        double y = table2->GetValue(i, 1).ToDouble(); // 获取 y 坐标
//        double z = table2->GetValue(i, 2).ToDouble(); // 获取 z 坐标
//
//        // 将点添加到 vtkPoints 中
//        cupoints3->InsertNextPoint(x, y, z);
//    }
//    // 创建 vtkPoints 来存储 CSV 文件中的点数据
//    vtkSmartPointer<vtkPoints> cupoints4 = vtkSmartPointer<vtkPoints>::New();
//    // 假设 CSV 文件的列是 x, y, z
//    for (vtkIdType i = 0; i < table3->GetNumberOfRows(); ++i) {
//        double x = table3->GetValue(i, 0).ToDouble(); // 获取 x 坐标
//        double y = table3->GetValue(i, 1).ToDouble(); // 获取 y 坐标
//        double z = table3->GetValue(i, 2).ToDouble(); // 获取 z 坐标
//
//        // 将点添加到 vtkPoints 中
//        cupoints4->InsertNextPoint(x, y, z);
//    }
//    // 创建 vtkPoints 来存储 CSV 文件中的点数据
//    vtkSmartPointer<vtkPoints> cupoints5 = vtkSmartPointer<vtkPoints>::New();
//    // 假设 CSV 文件的列是 x, y, z
//    for (vtkIdType i = 0; i < table4->GetNumberOfRows(); ++i) {
//        double x = table4->GetValue(i, 0).ToDouble(); // 获取 x 坐标
//        double y = table4->GetValue(i, 1).ToDouble(); // 获取 y 坐标
//        double z = table4->GetValue(i, 2).ToDouble(); // 获取 z 坐标
//
//        // 将点添加到 vtkPoints 中
//        cupoints5->InsertNextPoint(x, y, z);
//    }
//    // 创建 vtkPoints 来存储 CSV 文件中的点数据
//    vtkSmartPointer<vtkPoints> cupoints6= vtkSmartPointer<vtkPoints>::New();
//    // 假设 CSV 文件的列是 x, y, z
//    for (vtkIdType i = 0; i < table5->GetNumberOfRows(); ++i) {
//        double x = table5->GetValue(i, 0).ToDouble(); // 获取 x 坐标
//        double y = table5->GetValue(i, 1).ToDouble(); // 获取 y 坐标
//        double z = table5->GetValue(i, 2).ToDouble(); // 获取 z 坐标
//
//        // 将点添加到 vtkPoints 中
//        cupoints6->InsertNextPoint(x, y, z);
//    }
//    // 3. 使用 vtkSphereSource 创建球体
//    vtkSmartPointer<vtkSphereSource> cusphereSource1 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource1->SetRadius(0.2);  // 设置球体的半径
//    cusphereSource1->SetPhiResolution(10);  // 设置纵向分辨率
//    cusphereSource1->SetThetaResolution(10);  // 设置横向分辨率
//    cusphereSource1->Update();  // 更新球体源
//    vtkSmartPointer<vtkSphereSource> cusphereSource2 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource2->SetRadius(0.2);  // 设置球体的半径
//    cusphereSource2->SetPhiResolution(10);  // 设置纵向分辨率
//    cusphereSource2->SetThetaResolution(10);  // 设置横向分辨率
//    cusphereSource2->Update();  // 更新球体源
//    vtkSmartPointer<vtkSphereSource> cusphereSource3 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource3->SetRadius(0.2);  // 设置球体的半径
//    cusphereSource3->SetPhiResolution(10);  // 设置纵向分辨率
//    cusphereSource3->SetThetaResolution(10);  // 设置横向分辨率
//    cusphereSource3->Update();  // 更新球体源
//    vtkSmartPointer<vtkSphereSource> cusphereSource4 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource4->SetRadius(0.2);  // 设置球体的半径
//    cusphereSource4->SetPhiResolution(10);  // 设置纵向分辨率
//    cusphereSource4->SetThetaResolution(10);  // 设置横向分辨率
//    cusphereSource4->Update();  // 更新球体源
//    vtkSmartPointer<vtkSphereSource> cusphereSource5 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource5->SetRadius(0.2);  // 设置球体的半径
//    cusphereSource5->SetPhiResolution(10);  // 设置纵向分辨率
//    cusphereSource5->SetThetaResolution(10);  // 设置横向分辨率
//    cusphereSource5->Update();  // 更新球体源
//    vtkSmartPointer<vtkSphereSource> cusphereSource6 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource6->SetRadius(0.2);  // 设置球体的半径
//    cusphereSource6->SetPhiResolution(10);  // 设置纵向分辨率
//    cusphereSource6->SetThetaResolution(10);  // 设置横向分辨率
//    cusphereSource6->Update();  // 更新球体源
//    vtkSmartPointer<vtkPolyData> cupolyData1 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData1->SetPoints(cupoints1); // 设置点数据
//    vtkSmartPointer<vtkPolyData> cupolyData2 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData2->SetPoints(cupoints2); // 设置点数据
//    vtkSmartPointer<vtkPolyData> cupolyData3 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData3->SetPoints(cupoints3); // 设置点数据
//    vtkSmartPointer<vtkPolyData> cupolyData4 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData4->SetPoints(cupoints4); // 设置点数据
//    vtkSmartPointer<vtkPolyData> cupolyData5 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData5->SetPoints(cupoints5); // 设置点数据
//    vtkSmartPointer<vtkPolyData> cupolyData6 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData6->SetPoints(cupoints6); // 设置点数据
//    // 4. 使用 vtkGlyph3D 将每个点转换为球体
//    vtkSmartPointer<vtkGlyph3D> glyph3D1 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D1->SetInputData(cupolyData1);  // 输入点数据
//    glyph3D1->SetSourceConnection(cusphereSource1->GetOutputPort());  // 设置球体源
//    glyph3D1->Update();  // 更新生成的结果
//    vtkSmartPointer<vtkGlyph3D> glyph3D2 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D2->SetInputData(cupolyData2);  // 输入点数据
//    glyph3D2->SetSourceConnection(cusphereSource2->GetOutputPort());  // 设置球体源
//    glyph3D2->Update();  // 更新生成的结果
//    vtkSmartPointer<vtkGlyph3D> glyph3D3 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D3->SetInputData(cupolyData3);  // 输入点数据
//    glyph3D3->SetSourceConnection(cusphereSource3->GetOutputPort());  // 设置球体源
//    glyph3D3->Update();  // 更新生成的结果
//    vtkSmartPointer<vtkGlyph3D> glyph3D4 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D4->SetInputData(cupolyData4);  // 输入点数据
//    glyph3D4->SetSourceConnection(cusphereSource4->GetOutputPort());  // 设置球体源
//    glyph3D4->Update();  // 更新生成的结果
//    vtkSmartPointer<vtkGlyph3D> glyph3D5 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D5->SetInputData(cupolyData5);  // 输入点数据
//    glyph3D5->SetSourceConnection(cusphereSource5->GetOutputPort());  // 设置球体源
//    glyph3D5->Update();  // 更新生成的结果
//
//    vtkSmartPointer<vtkGlyph3D> glyph3D6 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D6->SetInputData(cupolyData6);  // 输入点数据
//    glyph3D6->SetSourceConnection(cusphereSource6->GetOutputPort());  // 设置球体源
//    glyph3D6->Update();  // 更新生成的结果
//
//
//     // 创建一个 Mapper 来渲染点数据
//    vtkSmartPointer<vtkPolyDataMapper> cumapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    cumapper1->SetInputConnection(glyph3D1->GetOutputPort());
//    // 创建一个 Mapper 来渲染点数据
//    vtkSmartPointer<vtkPolyDataMapper> cumapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    cumapper2->SetInputConnection(glyph3D2->GetOutputPort());
//
//    // 创建一个 Mapper 来渲染点数据
//    vtkSmartPointer<vtkPolyDataMapper> cumapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    cumapper3->SetInputConnection(glyph3D3->GetOutputPort());
//    // 创建一个 Mapper 来渲染点数据
//    vtkSmartPointer<vtkPolyDataMapper> cumapper4 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    cumapper4->SetInputConnection(glyph3D4->GetOutputPort());
//    // 创建一个 Mapper 来渲染点数据
//    vtkSmartPointer<vtkPolyDataMapper> cumapper5= vtkSmartPointer<vtkPolyDataMapper>::New();
//    cumapper5->SetInputConnection(glyph3D5->GetOutputPort());
//    vtkSmartPointer<vtkPolyDataMapper> cumapper6= vtkSmartPointer<vtkPolyDataMapper>::New();
//    cumapper6->SetInputConnection(glyph3D6->GetOutputPort());
//
//    vtkSmartPointer<vtkPolyDataMapper> filteredSkinMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    filteredSkinMapper->SetInputData(filteredSkin);
//    vtkSmartPointer<vtkPolyDataMapper> filteredSkinMapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    filteredSkinMapper1->SetInputData(filteredSkin1);
//    vtkSmartPointer<vtkPolyDataMapper> filteredSkinMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    filteredSkinMapper2->SetInputData(filteredSkin2);
//    vtkSmartPointer<vtkPolyDataMapper> filteredSkinMapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    filteredSkinMapper3->SetInputData(filteredSkin3);
//    vtkSmartPointer<vtkPolyDataMapper> filteredSkinMapper4 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    filteredSkinMapper4->SetInputData(filteredSkin4);
//    vtkSmartPointer<vtkPolyDataMapper> filteredSkinMapper5 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    filteredSkinMapper5->SetInputData(filteredSkin5);
//    vtkSmartPointer<vtkPolyDataMapper> filteredSkinMapper6 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    filteredSkinMapper6->SetInputData(filteredSkin6);
//    vtkSmartPointer<vtkPolyDataMapper> BoneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    BoneMapper->SetInputData(bone);
//    vtkSmartPointer<vtkPolyDataMapper> liverMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    liverMapper->SetInputData(liver);
//    vtkSmartPointer<vtkPolyDataMapper> tumorMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    tumorMapper->SetInputData(tumor);
//    vtkSmartPointer<vtkPolyDataMapper> SkinMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    SkinMapper->SetInputData(skin);
//    vtkSmartPointer<vtkPolyDataMapper> ArteryMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    ArteryMapper->SetInputData(artery);
//    vtkSmartPointer<vtkPolyDataMapper> venoussystemMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    venoussystemMapper->SetInputData(venoussystem);
//
//    // 创建椭球体（20, 15, 15 的椭球）
//    double a = 29;  // 长轴
//    double b = 22;  // 短轴
//    double c = 22;  // 短轴
//    vtkSmartPointer<vtkParametricEllipsoid> ellipsoid = vtkSmartPointer<vtkParametricEllipsoid>::New();
//    ellipsoid->SetXRadius(a);
//    ellipsoid->SetYRadius(b);
//    ellipsoid->SetZRadius(c);
//
//    vtkSmartPointer<vtkParametricFunctionSource> source = vtkSmartPointer<vtkParametricFunctionSource>::New();
//    source->SetParametricFunction(ellipsoid);
//    source->Update();
//    // 使用三角化算法处理椭球的表面
//    vtkSmartPointer<vtkTriangleFilter> triFilter = vtkSmartPointer<vtkTriangleFilter>::New();
//    triFilter->SetInputConnection(source->GetOutputPort());
//    triFilter->Update();
//
//    // 创建椭球体的映射器
//    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    mapper->SetInputConnection(triFilter->GetOutputPort());
//    // 创建椭球体的映射器
//    vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    mapper2->SetInputConnection(triFilter->GetOutputPort());
//    vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    mapper3->SetInputConnection(triFilter->GetOutputPort());
//    vtkSmartPointer<vtkPolyDataMapper> mapper4 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    mapper4->SetInputConnection(triFilter->GetOutputPort());
//
//
//
//    vtkSmartPointer<vtkPoints>cu1 = vtkSmartPointer<vtkPoints>::New();
//    vtkSmartPointer<vtkPoints>cu2 = vtkSmartPointer<vtkPoints>::New();
//    vtkSmartPointer<vtkPoints>cu3 = vtkSmartPointer<vtkPoints>::New();
//    vtkSmartPointer<vtkPoints>cu4 = vtkSmartPointer<vtkPoints>::New();
//
//
//
//    // 将中心坐标添加到 vtkPoints 中
//    cu1->InsertNextPoint(tumorCenter1);
//    // 创建椭球体的 Actor，并设置其颜色为红色
//    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//    actor->SetMapper(mapper);
//    actor->GetProperty()->SetColor(1, 0, 0);  // 设置颜色为红色
//    actor->SetPosition(tumorCenter1);  // 将椭球体放置到肿瘤质心位置
//    actor->GetProperty()->SetOpacity(0.8);
//
//    // 将中心坐标添加到 vtkPoints 中
//    cu2->InsertNextPoint(tumorCenter2);
//    // 创建椭球体的 Actor，并设置其颜色为红色
//    vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
//    actor2->SetMapper(mapper2);
//    actor2->GetProperty()->SetColor(1, 1, 0);  // 设置颜色为红色
//    actor2->SetPosition(tumorCenter2);  // 将椭球体放置到肿瘤质心位置
//    actor2->GetProperty()->SetOpacity(0.8);
//    // 将中心坐标添加到 vtkPoints 中
//    cu3->InsertNextPoint(tumorCenter3);
//    // 创建椭球体的 Actor，并设置其颜色为红色
//    vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
//    actor3->SetMapper(mapper3);
//    actor3->GetProperty()->SetColor(0, 1, 0);  // 设置颜色为红色
//    actor3->SetPosition(tumorCenter3);  // 将椭球体放置到肿瘤质心位置
//    actor3->GetProperty()->SetOpacity(0.8);
//    // 将中心坐标添加到 vtkPoints 中
//    cu4->InsertNextPoint(tumorCenter4);
//    // 创建椭球体的 Actor，并设置其颜色为红色
//    vtkSmartPointer<vtkActor> actor4 = vtkSmartPointer<vtkActor>::New();
//    actor4->SetMapper(mapper4);
//    actor4->GetProperty()->SetColor(0, 0, 1);  // 设置颜色为红色
//    actor4->SetPosition(tumorCenter4);  // 将椭球体放置到肿瘤质心位置
//    actor4->GetProperty()->SetOpacity(0.8);
//
//    // 创建一个 Actor 来显示点
//    vtkSmartPointer<vtkActor> cuactor1 = vtkSmartPointer<vtkActor>::New();
//    cuactor1->SetMapper(cumapper1);
//    cuactor1->GetProperty()->SetColor(1, 1, 0);//huang色
//    cuactor1->GetProperty()->SetPointSize(2.5);
//    // 创建一个 Actor 来显示点
//    vtkSmartPointer<vtkActor> cuactor2 = vtkSmartPointer<vtkActor>::New();
//    cuactor2->SetMapper(cumapper2);
//    cuactor2->GetProperty()->SetColor(1, 0, 0);//蓝色
//    cuactor2->GetProperty()->SetPointSize(2.5);
//    // 创建一个 Actor 来显示点
//    vtkSmartPointer<vtkActor> cuactor3 = vtkSmartPointer<vtkActor>::New();
//    cuactor3->SetMapper(cumapper3);
//    cuactor3->GetProperty()->SetColor(0.0, 1.0, 0.0);//lv se
//    cuactor3->GetProperty()->SetPointSize(2.5);
//    // 创建一个 Actor 来显示点
//    vtkSmartPointer<vtkActor> cuactor4 = vtkSmartPointer<vtkActor>::New();
//    cuactor4->SetMapper(cumapper4);
//    cuactor4->GetProperty()->SetColor(1.0, 0.0, 0.0);//red
//    cuactor4->GetProperty()->SetPointSize(2.5);
//    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//    // 构建 OBBTree 来进行碰撞检测，假设你已经有了 bone 数据
//    vtkSmartPointer<vtkOBBTree> obbTree1 = vtkSmartPointer<vtkOBBTree>::New();
//    obbTree1->SetDataSet(bone);
//    obbTree1->BuildLocator();
//    // 调用优化算法并获取优化后的起点位置
//// 假设 ExtractPointsFromTrianglesAndProcess 是您用来提取点集的方法
//    vtkSmartPointer<vtkPolyData> filteredPointsPolyData1 = ExtractPointsFromTrianglesAndProcess(filteredSkin4, tumorCenter4);
//    vtkSmartPointer<vtkPoints> points1 = filteredPointsPolyData1->GetPoints();
//    // 获取点的数量
//    vtkIdType numPoints = points1->GetNumberOfPoints();
//    // 遍历并获取优化后的点
//
//        vtkIdType randomIndex = rand() % numPoints;  // rand() 返回 [0, numPoints-1] 之间的随机数
//
//        double pathStart[3];
//        points1->GetPoint(randomIndex, pathStart);  // 获取每个点的坐标
//        // 在控制台输出随机选中的点坐标
//        std::cout << "绿色点坐标: ("
//            << pathStart[0] << ", "
//            << pathStart[1] << ", "
//            << pathStart[2] << ")"
//            << std::endl;
//
//        // 为每个起点创建一个球体
//        vtkSmartPointer<vtkSphereSource> GreensphereDataop = vtkSmartPointer<vtkSphereSource>::New();
//        GreensphereDataop->SetCenter(pathStart[0], pathStart[1], pathStart[2]);
//       /* std::cout << "pathStart[0]: " << pathStart[0] << std::endl;
//        std::cout << "pathStart[1]: " << pathStart[1] << std::endl;
//        std::cout << "pathStart[2]: " << pathStart[2] << std::endl;*/
//
//        GreensphereDataop->SetRadius(0.5);  // 设置球体半径为 0.5
//            // 计算从 tumorCenter4 到 pathStart 的方向向量
//        double dx = pathStart[0] - tumorCenter4[0];
//        double dy = pathStart[1] - tumorCenter4[1];
//        double dz = pathStart[2] - tumorCenter4[2];
//
//        // 将方向向量单位化
//        normalizeDirection(dx, dy, dz);
//
//        // 计算从 tumorCenter4 开始，超过 pathStart 20 单位的射线终点
//        double extraLength = 50.0;
//        double totalLength = distance(tumorCenter4[0], tumorCenter4[1], tumorCenter4[2], pathStart[0], pathStart[1], pathStart[2]) + extraLength;
//
//        // 计算射线终点
//        double rayEnd[3] = {
//            tumorCenter4[0] + dx * totalLength,
//            tumorCenter4[1] + dy * totalLength,
//            tumorCenter4[2]-50 + dz * totalLength
//        };
//
//        std::cout << "Ray end point: ("
//            << rayEnd[0] << ", "
//            << rayEnd[1] << ", "
//            << rayEnd[2] << ")" << std::endl;
//
//        // VTK可视化代码
//
//        // 创建一个点源对象来设置起始点
//        vtkSmartPointer<vtkPoints> Greenpoints = vtkSmartPointer<vtkPoints>::New();
//        Greenpoints->InsertNextPoint(tumorCenter4); // 设置起点 (tumorCenter4)
//        Greenpoints->InsertNextPoint(rayEnd);       // 设置终点 (rayEnd)
//
//        // 创建一个线源对象
//        vtkSmartPointer<vtkCellArray> Greenlines = vtkSmartPointer<vtkCellArray>::New();
//        Greenlines->InsertNextCell(2); // 一条线由两个点组成
//        Greenlines->InsertCellPoint(0); // 起点
//        Greenlines->InsertCellPoint(1); // 终点
//
//        // 创建一个PolyData对象
//        vtkSmartPointer<vtkPolyData> GreenlineData = vtkSmartPointer<vtkPolyData>::New();
//        GreenlineData->SetPoints(Greenpoints);
//        GreenlineData->SetLines(Greenlines);
//        // 使用 vtkTubeFilter 创建圆柱体
//        vtkSmartPointer<vtkTubeFilter>GreentubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
//        GreentubeFilter->SetInputData(GreenlineData);   // 设置线段数据作为输入
//        GreentubeFilter->SetRadius(1.0);            // 设置圆柱体的半径（相当于“线段的半径”）
//        GreentubeFilter->SetNumberOfSides(50);      // 设置圆柱体的边数，越多越光滑
//        // 创建一个映射器来显示线
//        vtkSmartPointer<vtkPolyDataMapper> GreenlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        GreenlineMapper->SetInputConnection(GreentubeFilter->GetOutputPort());
//
//        // 创建一个Actor来渲染线
//        vtkSmartPointer<vtkActor> GreenlineActor = vtkSmartPointer<vtkActor>::New();
//        GreenlineActor->SetMapper(GreenlineMapper);
//
//        // 设置线条的颜色
//        GreenlineActor->GetProperty()->SetColor(0, 1, 0);  // 黄色
//
//        renderer->AddActor(GreenlineActor);//-------------------------------------------------------------------------------------------------------------------------------------
//        // 计算椭球体长轴（20）和短轴（15）
//        double GreensemiMajorAxis = 22.0;  // 长轴
//        double GreensemiMinorAxis = 17;  // 短轴
//
//        // 创建一个椭球体的参数化函数
//        vtkSmartPointer<vtkParametricEllipsoid> Greenellipsoid = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        Greenellipsoid->SetXRadius(GreensemiMinorAxis);  // 设置长轴
//        Greenellipsoid->SetYRadius(GreensemiMinorAxis);  // 设置短轴
//        Greenellipsoid->SetZRadius(GreensemiMajorAxis);  // 设置短轴
//
//        // 创建一个参数化函数源
//        vtkSmartPointer<vtkParametricFunctionSource> GreenellipsoidSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        GreenellipsoidSource->SetParametricFunction(Greenellipsoid);
//        double  M_PI = 3.14;
//        // 创建椭球体的变换
//        vtkSmartPointer<vtkTransform> Greentransform = vtkSmartPointer<vtkTransform>::New();
//        // **平移操作**：将椭球体的质心从 (0, 0, 0) 移动到 tumorCenter4
//        Greentransform->Translate(tumorCenter4[0], tumorCenter4[1], tumorCenter4[2]);
//
//        // 创建旋转矩阵，目标是将椭球体的长轴（默认沿 z 轴）旋转到目标方向（dx, dy, dz）
//        double originalDirection[3] = { 0.0, 0.0, 1.0 };  // 默认长轴方向
//        double rotationAngle = 0.0;
//        double rotationAxis[3] = { 0.0, 0.0, 0.0 };
//
//        // 计算旋转角度和旋转轴
//        double dotProduct = vtkMath::Dot(originalDirection, new double[3]{ dx, dy, dz });
//        double crossProduct[3];
//        vtkMath::Cross(originalDirection, new double[3]{ dx, dy, dz }, crossProduct);
//
//        // 如果叉积为零，表示两个向量已经重合，无需旋转
//        if (vtkMath::Norm(crossProduct) < 1e-6) {
//            // 如果没有旋转（两个向量重合），则无需做旋转
//            std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//            rotationAngle = 0.0;
//        }
//        else {
//            rotationAngle = acos(dotProduct);  // 计算夹角
//            vtkMath::Normalize(crossProduct);  // 归一化旋转轴
//            std::cout << "Rotation Axis: (" << crossProduct[0] << ", " << crossProduct[1] << ", " << crossProduct[2] << ")" << std::endl;
//        }
//
//        // 旋转椭球体到目标方向
//        Greentransform->RotateWXYZ(rotationAngle * 180.0 / M_PI, crossProduct[0], crossProduct[1], crossProduct[2]);
//        // 创建一个变换过滤器来应用变换
//        vtkSmartPointer<vtkTransformPolyDataFilter> GreentransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        GreentransformFilter->SetInputConnection(GreenellipsoidSource->GetOutputPort());
//        GreentransformFilter->SetTransform(Greentransform);
//
//        // 创建一个映射器
//        vtkSmartPointer<vtkPolyDataMapper> Greenmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        Greenmapper->SetInputConnection(GreentransformFilter->GetOutputPort());
//
//        // 创建一个Actor
//        vtkSmartPointer<vtkActor> Greenactor = vtkSmartPointer<vtkActor>::New();
//        Greenactor->SetMapper(Greenmapper);
//        Greenactor->GetProperty()->SetOpacity(0.9);
//
//        // 设置颜色
//        Greenactor->GetProperty()->SetColor(0.0, 1.0, 0.0);  // 绿色  对应的是cuactor3
//        renderer->AddActor(Greenactor);//----------------------------------------------------------------------------------------------------------------------------------------------
//     // 创建映射器并将球体添加到渲染器中
//        vtkSmartPointer<vtkPolyDataMapper> sphereMapperop = vtkSmartPointer<vtkPolyDataMapper>::New();
//        sphereMapperop->SetInputConnection(GreensphereDataop->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> sphereActorop = vtkSmartPointer<vtkActor>::New();
//        sphereActorop->SetMapper(sphereMapperop);
//        sphereActorop->GetProperty()->SetColor(139 / 255.0, 139 / 255.0, 0);//黄色
//        // 将球体添加到渲染器中
//        //renderer->AddActor(sphereActorop);
//
//        // 创建并可视化优化点与肿瘤质心之间的连线
//        vtkSmartPointer<vtkLineSource> lineSourceop = vtkSmartPointer<vtkLineSource>::New();
//        lineSourceop->SetPoint1(pathStart[0], pathStart[1], pathStart[2]);  // 起点
//        lineSourceop->SetPoint2(tumorCenter4[0], tumorCenter4[1], tumorCenter4[2]);  // 肿瘤质心
//         // 使用 TubeFilter 增加线段的半径
//        vtkSmartPointer<vtkTubeFilter> tubeFilterop = vtkSmartPointer<vtkTubeFilter>::New();
//        tubeFilterop->SetInputConnection(lineSourceop->GetOutputPort());
//        tubeFilterop->SetRadius(0.5);  // 设置线段半径
//        tubeFilterop->SetNumberOfSides(20);  // 设置线段的精度
//        tubeFilterop->Update();
//        // 创建映射器并将线添加到渲染器中
//        vtkSmartPointer<vtkPolyDataMapper> lineMapperop = vtkSmartPointer<vtkPolyDataMapper>::New();
//        lineMapperop->SetInputConnection(tubeFilterop->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> lineActorop = vtkSmartPointer<vtkActor>::New();
//        lineActorop->SetMapper(lineMapperop);
//        lineActorop->GetProperty()->SetColor(139 / 255.0, 139 / 255.0, 0);  // 设置线的颜色为黄色
//
//        // 将线添加到渲染器中
//       // renderer->AddActor(lineActorop);
////-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//    // 构建 OBBTree 来进行碰撞检测，假设你已经有了 bone 数据
//    vtkSmartPointer<vtkOBBTree> obbTree2 = vtkSmartPointer<vtkOBBTree>::New();
//    obbTree2->SetDataSet(bone);
//    obbTree2->BuildLocator();
//    // 调用优化算法并获取优化后的起点位置
//    vtkSmartPointer<vtkPolyData> filteredPointsPolyData2 = ExtractPointsFromTrianglesAndProcess(filteredSkin3, tumorCenter3);
//    vtkSmartPointer<vtkPoints> points2 = filteredPointsPolyData2->GetPoints();
//    // 获取点的数量
//    vtkIdType numPoints2 = points2->GetNumberOfPoints();
//    // 遍历并获取优化后的点
//
//    vtkIdType randomIndex2 = rand() % numPoints2;  // rand() 返回 [0, numPoints-1] 之间的随机数
//
//    double pathStart1[3];
//    points1->GetPoint(randomIndex, pathStart);  // 获取每个点的坐标
//    // 在控制台输出随机选中的点坐标
//    std::cout << "黄色色点坐标: ("
//        << pathStart[0] << ", "
//        << pathStart[1] << ", "
//        << pathStart[2] << ")"
//        << std::endl;
//
//    // 为每个起点创建一个球体
//        // 为每个起点创建一个球体
//        vtkSmartPointer<vtkSphereSource> yellowsphereDataop1 = vtkSmartPointer<vtkSphereSource>::New();
//        yellowsphereDataop1->SetCenter(pathStart1[0], pathStart1[1], pathStart1[2]);
//        yellowsphereDataop1->SetRadius(0.5);  // 设置球体半径为 0.5
//         // 创建一个点源对象来设置起始点
//
//           // 计算从 tumorCenter4 到 pathStart 的方向向量
//        double dx1 = pathStart1[0] - tumorCenter3[0];
//        double dy1 = pathStart1[1] - tumorCenter3[1];
//        double dz1 = pathStart1[2] +60- tumorCenter3[2];
//
//        // 将方向向量单位化
//        normalizeDirection(dx1, dy1, dz1);
//
//        // 计算从 tumorCenter4 开始，超过 pathStart 20 单位的射线终点
//        double extraLength2 = -80;
//        double totalLength2 = distance(tumorCenter3[0], tumorCenter3[1], tumorCenter3[2], pathStart1[0], pathStart1[1], pathStart1[2]) + extraLength2;
//
//        // 计算射线终点
//        double rayEnd2[3] = {
//            tumorCenter3[0] + dx1 * totalLength2,
//            tumorCenter3[1] + dy1 * totalLength2,
//            tumorCenter3[2] + dz1 * totalLength2
//        };
//
//        vtkSmartPointer<vtkPoints> yellowpoints = vtkSmartPointer<vtkPoints>::New();
//        yellowpoints->InsertNextPoint(tumorCenter3); // 设置起点 (tumorCenter4)
//        yellowpoints->InsertNextPoint(rayEnd2);       // 设置终点 (rayEnd)
//
//        // 创建一个线源对象
//        vtkSmartPointer<vtkCellArray> yellowlines = vtkSmartPointer<vtkCellArray>::New();
//        yellowlines->InsertNextCell(2); // 一条线由两个点组成
//        yellowlines->InsertCellPoint(0); // 起点
//        yellowlines->InsertCellPoint(1); // 终点
//
//        // 创建一个PolyData对象
//        vtkSmartPointer<vtkPolyData> yellowlineData = vtkSmartPointer<vtkPolyData>::New();
//        yellowlineData->SetPoints(yellowpoints);
//        yellowlineData->SetLines(yellowlines);
//        // 使用 vtkTubeFilter 创建圆柱体
//        vtkSmartPointer<vtkTubeFilter>yellowtubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
//        yellowtubeFilter->SetInputData(yellowlineData);   // 设置线段数据作为输入
//        yellowtubeFilter->SetRadius(1.0);            // 设置圆柱体的半径（相当于“线段的半径”）
//        yellowtubeFilter->SetNumberOfSides(50);      // 设置圆柱体的边数，越多越光滑
//        // 创建一个映射器来显示线
//        vtkSmartPointer<vtkPolyDataMapper> yellowlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        yellowlineMapper->SetInputConnection(yellowtubeFilter->GetOutputPort());
//
//        // 创建一个Actor来渲染线
//        vtkSmartPointer<vtkActor> yellowlineActor = vtkSmartPointer<vtkActor>::New();
//        yellowlineActor->SetMapper(yellowlineMapper);
//
//        // 设置线条的颜色
//        yellowlineActor->GetProperty()->SetColor(1, 1, 0);  // 黄色
//
//        renderer->AddActor(yellowlineActor);//-------------------------------------------------------------------------------------------------------------------------------------------------
//        // 计算椭球体长轴（20）和短轴（15）
//        double yellowsemiMajorAxis = 22.0;  // 长轴
//        double yellowsemiMinorAxis = 17;  // 短轴
//        //double yellowsemiMajorAxis1 = 25.0;  // 长轴
//        //double yellowsemiMinorAxis1 = 22;  // 短轴
//
//        // 创建一个椭球体的参数化函数
//        vtkSmartPointer<vtkParametricEllipsoid> yellowellipsoid = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        yellowellipsoid->SetXRadius(yellowsemiMinorAxis);  // 设置长轴
//        yellowellipsoid->SetYRadius(yellowsemiMinorAxis);  // 设置短轴
//        yellowellipsoid->SetZRadius(yellowsemiMajorAxis);  // 设置短轴
//        // // 创建一个椭球体的参数化函数
//        //vtkSmartPointer<vtkParametricEllipsoid> yellowellipsoid1 = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        //yellowellipsoid1->SetXRadius(yellowsemiMinorAxis1);  // 设置长轴
//        //yellowellipsoid1->SetYRadius(yellowsemiMinorAxis1);  // 设置短轴
//        //yellowellipsoid1->SetZRadius(yellowsemiMajorAxis1);  // 设置短轴
//
//        // 创建一个参数化函数源
//        vtkSmartPointer<vtkParametricFunctionSource> yellowellipsoidSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        yellowellipsoidSource->SetParametricFunction(yellowellipsoid);
//        //// 创建一个参数化函数源
//        //vtkSmartPointer<vtkParametricFunctionSource> yellowellipsoidSource1 = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        //yellowellipsoidSource1->SetParametricFunction(yellowellipsoid1);
//   
//        // 创建椭球体的变换
//        vtkSmartPointer<vtkTransform> yellowtransform = vtkSmartPointer<vtkTransform>::New();
//    /*    vtkSmartPointer<vtkTransform> yellowtransform1 = vtkSmartPointer<vtkTransform>::New();*/
//        // **平移操作**：将椭球体的质心从 (0, 0, 0) 移动到 tumorCenter4
//        yellowtransform->Translate(tumorCenter3[0], tumorCenter3[1], tumorCenter3[2]);
//    /*    yellowtransform1->Translate(tumorCenter1[0] - 5, tumorCenter1[1] - 5, tumorCenter1[2] - 5);*/
//        // 创建旋转矩阵，目标是将椭球体的长轴（默认沿 z 轴）旋转到目标方向（dx, dy, dz）
//        double originalDirection1[3] = { 0.0, 0.0, 1.0 };  // 默认长轴方向
//        double rotationAngle1 = 0.0;
//        double rotationAxis1[3] = { 0.0, 0.0, 0.0 };
//
//        // 计算旋转角度和旋转轴
//        double dotProduct1 = vtkMath::Dot(originalDirection1, new double[3]{ dx1, dy1, dz1 });
//        double crossProduct1[3];
//        vtkMath::Cross(originalDirection1, new double[3]{ dx1, dy1, dz1 }, crossProduct1);
//
//        // 如果叉积为零，表示两个向量已经重合，无需旋转
//        if (vtkMath::Norm(crossProduct1) < 1e-6) {
//            // 如果没有旋转（两个向量重合），则无需做旋转
//            std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//            rotationAngle1 = 0.0;
//        }
//        else {
//            rotationAngle1 = acos(dotProduct1);  // 计算夹角
//            vtkMath::Normalize(crossProduct1);  // 归一化旋转轴
//            std::cout << "Rotation Axis: (" << crossProduct1[0] << ", " << crossProduct1[1] << ", " << crossProduct1[2] << ")" << std::endl;
//        }
//
//        // 旋转椭球体到目标方向
//        yellowtransform->RotateWXYZ(rotationAngle1 * 180.0 / M_PI, crossProduct1[0], crossProduct1[1], crossProduct1[2]);
//
//        // 创建旋转矩阵，目标是将椭球体的长轴（默认沿 z 轴）旋转到目标方向（dx, dy, dz）
//        double originalDirection11[3] = { 0.0, 0.0, 1.0 };  // 默认长轴方向
//        double rotationAngle11 = 0.0;
//        double rotationAxis11[3] = { 0.0, 0.0, 0.0 };
//
//        //// 计算旋转角度和旋转轴
//        //double dotProduct11 = vtkMath::Dot(originalDirection11, new double[3]{ dx1, dy1, dz1 });
//        //double crossProduct11[3];
//        //vtkMath::Cross(originalDirection11, new double[3]{ dx1, dy1, dz1 }, crossProduct11);
//
//        //// 如果叉积为零，表示两个向量已经重合，无需旋转
//        //if (vtkMath::Norm(crossProduct11) < 1e-6) {
//        //    // 如果没有旋转（两个向量重合），则无需做旋转
//        //    std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//        //    rotationAngle11 = 0.0;
//        //}
//        //else {
//        //    rotationAngle11 = acos(dotProduct11);  // 计算夹角
//        //    vtkMath::Normalize(crossProduct11);  // 归一化旋转轴
//        //    std::cout << "Rotation Axis: (" << crossProduct11[0] << ", " << crossProduct11[1] << ", " << crossProduct11[2] << ")" << std::endl;
//        //}
//
//        //// 旋转椭球体到目标方向
//        //yellowtransform1->RotateWXYZ(rotationAngle11 * 180.0 / M_PI, crossProduct11[0], crossProduct11[1], crossProduct11[2]);
//
//
//        // 创建一个变换过滤器来应用变换
//        vtkSmartPointer<vtkTransformPolyDataFilter> yellowtransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        yellowtransformFilter->SetInputConnection(yellowellipsoidSource->GetOutputPort());
//        yellowtransformFilter->SetTransform(yellowtransform);
//
//        // 创建一个映射器
//        vtkSmartPointer<vtkPolyDataMapper> yellowmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        yellowmapper->SetInputConnection(yellowtransformFilter->GetOutputPort());
//
//        // 创建一个Actor
//        vtkSmartPointer<vtkActor> yellowactor = vtkSmartPointer<vtkActor>::New();
//        yellowactor->SetMapper(yellowmapper);
//        yellowactor->GetProperty()->SetOpacity(0.9);
//
//        // 设置颜色
//        yellowactor->GetProperty()->SetColor(1.0, 1.0, 0.0);  //  cuactor1
//        renderer->AddActor(yellowactor);//--------------------------------------------------------------------------------------------------------------------------
//        //// 创建一个变换过滤器来应用变换
//        //vtkSmartPointer<vtkTransformPolyDataFilter> yellowtransformFilter1 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        //yellowtransformFilter1->SetInputConnection(yellowellipsoidSource1->GetOutputPort());
//        //yellowtransformFilter1->SetTransform(yellowtransform1);
//
//        //// 创建一个映射器
//        //vtkSmartPointer<vtkPolyDataMapper> yellowmapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        //yellowmapper1->SetInputConnection(yellowtransformFilter1->GetOutputPort());
//
//        //// 创建一个Actor
//        //vtkSmartPointer<vtkActor> yellowactor1 = vtkSmartPointer<vtkActor>::New();
//        //yellowactor1->SetMapper(yellowmapper1);
//        //yellowactor1->GetProperty()->SetOpacity(0.9);
//
//        //// 设置颜色
//        //yellowactor1->GetProperty()->SetColor(1.0, 1.0, 0.0);  //  cuactor1
//        //renderer->AddActor(yellowactor1);//-----------------------------------------------------------------------------------------------------------------------------
//        // // 创建映射器并将球体添加到渲染器中
//        vtkSmartPointer<vtkPolyDataMapper> sphereMapperop1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        sphereMapperop1->SetInputConnection(yellowsphereDataop1->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> sphereActorop1 = vtkSmartPointer<vtkActor>::New();
//        sphereActorop1->SetMapper(sphereMapperop1);
//        sphereActorop1->GetProperty()->SetColor(0, 100 / 255.0, 0);//绿色
//        // 将球体添加到渲染器中
//        //renderer->AddActor(sphereActorop1);
//
//        // 创建并可视化优化点与肿瘤质心之间的连线
//        vtkSmartPointer<vtkLineSource> lineSourceop1 = vtkSmartPointer<vtkLineSource>::New();
//        lineSourceop1->SetPoint1(pathStart1[0], pathStart1[1], pathStart1[2]);  // 起点
//        lineSourceop1->SetPoint2(tumorCenter3[0], tumorCenter3[1], tumorCenter3[2]);  // 肿瘤质心
//         // 使用 TubeFilter 增加线段的半径
//        vtkSmartPointer<vtkTubeFilter> tubeFilterop1 = vtkSmartPointer<vtkTubeFilter>::New();
//        tubeFilterop1->SetInputConnection(lineSourceop1->GetOutputPort());
//        tubeFilterop1->SetRadius(0.5);  // 设置线段半径
//        tubeFilterop1->SetNumberOfSides(20);  // 设置线段的精度
//        tubeFilterop1->Update();
//        // 创建映射器并将线添加到渲染器中
//        vtkSmartPointer<vtkPolyDataMapper> lineMapperop1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        lineMapperop1->SetInputConnection(tubeFilterop1->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> lineActorop1 = vtkSmartPointer<vtkActor>::New();
//        lineActorop1->SetMapper(lineMapperop1);
//        lineActorop1->GetProperty()->SetColor(0, 100 / 255.0, 0);  // 设置线的颜色为绿色
//
//        // 将线添加到渲染器中
//       //renderer->AddActor(lineActorop1);
//    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//    //// 构建 OBBTree 来进行碰撞检测，假设你已经有了 bone 数据
//    //vtkSmartPointer<vtkOBBTree> obbTree3 = vtkSmartPointer<vtkOBBTree>::New();
//    //obbTree3->SetDataSet(bone);
//    //obbTree3->BuildLocator();
//    //// 调用优化算法并获取优化后的起点位置
//        vtkSmartPointer<vtkPolyData> filteredPointsPolyData3 = ExtractPointsFromTrianglesAndProcess(filteredSkin1, tumorCenter1);
//        vtkSmartPointer<vtkPoints> points3 = filteredPointsPolyData3->GetPoints();
//        // 获取点的数量
//        vtkIdType numPoints3 = points3->GetNumberOfPoints();
//        // 遍历并获取优化后的点
//
//        vtkIdType randomIndex3 = rand() % numPoints3;  // rand() 返回 [0, numPoints-1] 之间的随机数
//
//        double pathStart2[3];
//        points3->GetPoint(randomIndex3, pathStart2);  // 获取每个点的坐标
//        // 在控制台输出随机选中的点坐标
//        std::cout << "红色点坐标: ("
//            << pathStart[0] << ", "
//            << pathStart[1] << ", "
//            << pathStart[2] << ")"
//            << std::endl;
//
//        // 为每个起点创建一个球体
//        vtkSmartPointer<vtkSphereSource> redsphereDataop2 = vtkSmartPointer<vtkSphereSource>::New();
//        redsphereDataop2->SetCenter(pathStart2[0], pathStart2[1], pathStart2[2]);
//        redsphereDataop2->SetRadius(0.5);  // 设置球体半径为 0.5
//
//        // 计算从 tumorCenter4 到 pathStart 的方向向量
//        double dx3 = pathStart2[0]+10  - tumorCenter1[0];
//        double dy3 = pathStart2[1] - tumorCenter1[1];
//        double dz3 = pathStart2[2]+20 - tumorCenter1[2];
//
//        // 将方向向量单位化
//        normalizeDirection(dx3, dy3, dz3);
//
//        // 计算从 tumorCenter4 开始，超过 pathStart 20 单位的射线终点
//        double extraLength3 = 50.0;
//        double totalLength3 = distance(tumorCenter1[0], tumorCenter1[1], tumorCenter1[2], pathStart2[0], pathStart2[1], pathStart2[2]) + extraLength3;
//
//        // 计算射线终点
//        double rayEnd3[3] = {
//            tumorCenter1[0] + dx3 * totalLength3,
//            tumorCenter1[1] + dy3 * totalLength3,
//            tumorCenter1[2] + dz3 * totalLength3
//        };
//
//        vtkSmartPointer<vtkPoints> redpoints = vtkSmartPointer<vtkPoints>::New();
//        redpoints->InsertNextPoint(tumorCenter1); // 设置起点 (tumorCenter4)
//        redpoints->InsertNextPoint(rayEnd3);       // 设置终点 (rayEnd)
//
//        // 创建一个线源对象
//        vtkSmartPointer<vtkCellArray> redlines = vtkSmartPointer<vtkCellArray>::New();
//        redlines->InsertNextCell(2); // 一条线由两个点组成
//        redlines->InsertCellPoint(0); // 起点
//        redlines->InsertCellPoint(1); // 终点
//
//        // 创建一个PolyData对象
//        vtkSmartPointer<vtkPolyData> redlineData = vtkSmartPointer<vtkPolyData>::New();
//        redlineData->SetPoints(redpoints);
//        redlineData->SetLines(redlines);
//        // 使用 vtkTubeFilter 创建圆柱体
//        vtkSmartPointer<vtkTubeFilter>redtubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilter->SetInputData(redlineData);   // 设置线段数据作为输入
//        redtubeFilter->SetRadius(1);            // 设置圆柱体的半径（相当于“线段的半径”）
//        redtubeFilter->SetNumberOfSides(50);      // 设置圆柱体的边数，越多越光滑
//
//        // 创建映射器并将球体添加到渲染器中
//        vtkSmartPointer<vtkPolyDataMapper> redsphereMapperop2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redsphereMapperop2->SetInputConnection(redsphereDataop2->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redsphereActorop2 = vtkSmartPointer<vtkActor>::New();
//        redsphereActorop2->SetMapper(redsphereMapperop2);
//        redsphereActorop2->GetProperty()->SetColor(139 / 255.0, 0, 0);//红色
//        // 将球体添加到渲染器中
//        //renderer->AddActor(sphereActorop2);
//
//        // 创建并可视化优化点与肿瘤质心之间的连线
//        vtkSmartPointer<vtkLineSource> redlineSourceop2 = vtkSmartPointer<vtkLineSource>::New();
//        redlineSourceop2->SetPoint1(pathStart2[0], pathStart2[1], pathStart2[2]);  // 起点
//        redlineSourceop2->SetPoint2(tumorCenter1[0], tumorCenter1[1], tumorCenter1[2]);  // 肿瘤质心
//         // 使用 TubeFilter 增加线段的半径
//        vtkSmartPointer<vtkTubeFilter> redtubeFilterop2 = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilterop2->SetInputConnection(redlineSourceop2->GetOutputPort());
//        redtubeFilterop2->SetRadius(1);  // 设置线段半径
//        redtubeFilterop2->SetNumberOfSides(50);  // 设置线段的精度
//        redtubeFilterop2->Update();
//        // 创建映射器并将线添加到渲染器中
//        vtkSmartPointer<vtkPolyDataMapper> redlineMapperop2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redlineMapperop2->SetInputConnection(redtubeFilter->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redlineActorop2 = vtkSmartPointer<vtkActor>::New();
//        redlineActorop2->SetMapper(redlineMapperop2);
//        redlineActorop2->GetProperty()->SetColor(0, 0, 139 / 255.0);  // 设置线的颜色为红色
//
//        // 将线添加到渲染器中
//        renderer->AddActor(redlineActorop2);//------------------------------------------------------------------------------------------------------
//
//        // 计算椭球体长轴（20）和短轴（15）
//        double redsemiMajorAxis = 22.0;  // 长轴
//        double redsemiMinorAxis =17.0;  // 短轴
//
//        // 创建一个椭球体的参数化函数
//        vtkSmartPointer<vtkParametricEllipsoid> redellipsoid = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        redellipsoid->SetXRadius(redsemiMinorAxis);  // 设置长轴
//        redellipsoid->SetYRadius(redsemiMinorAxis);  // 设置短轴
//        redellipsoid->SetZRadius(redsemiMajorAxis);  // 设置短轴
//
//        // 创建一个参数化函数源
//        vtkSmartPointer<vtkParametricFunctionSource> redellipsoidSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        redellipsoidSource->SetParametricFunction(redellipsoid);
//       
//        // 创建椭球体的变换
//        vtkSmartPointer<vtkTransform> redtransform = vtkSmartPointer<vtkTransform>::New();
//        // **平移操作**：将椭球体的质心从 (0, 0, 0) 移动到 tumorCenter4
//        redtransform->Translate(tumorCenter1[0], tumorCenter1[1], tumorCenter1[2]);
//        // 创建旋转矩阵，目标是将椭球体的长轴（默认沿 z 轴）旋转到目标方向（dx, dy, dz）
//        double originalDirection2[3] = { 0.0, 0.0, 1.0 };  // 默认长轴方向
//        double rotationAngle2 = 0.0;
//        double rotationAxis2[3] = { 0.0, 0.0, 0.0 };
//
//        // 计算旋转角度和旋转轴
//        double dotProduct2 = vtkMath::Dot(originalDirection2, new double[3]{ dx3, dy3, dz3 });
//        double crossProduct2[3];
//        vtkMath::Cross(originalDirection2, new double[3]{ dx3, dy3, dz3 }, crossProduct2);
//
//        // 如果叉积为零，表示两个向量已经重合，无需旋转
//        if (vtkMath::Norm(crossProduct2) < 1e-6) {
//            // 如果没有旋转（两个向量重合），则无需做旋转
//            std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//            rotationAngle2 = 0.0;
//        }
//        else {
//            rotationAngle2 = acos(dotProduct2);  // 计算夹角
//            vtkMath::Normalize(crossProduct2);  // 归一化旋转轴
//            std::cout << "Rotation Axis: (" << crossProduct2[0] << ", " << crossProduct2[1] << ", " << crossProduct2[2] << ")" << std::endl;
//        }
//
//        // 旋转椭球体到目标方向
//        redtransform->RotateWXYZ(rotationAngle2 * 180.0 / M_PI, crossProduct2[0], crossProduct2[1], crossProduct2[2]);
//
//
//
//        // 创建一个变换过滤器来应用变换
//        vtkSmartPointer<vtkTransformPolyDataFilter> redtransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        redtransformFilter->SetInputConnection(redellipsoidSource->GetOutputPort());
//        redtransformFilter->SetTransform(redtransform);
//
//        // 创建一个映射器
//        vtkSmartPointer<vtkPolyDataMapper> redmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redmapper->SetInputConnection(redtransformFilter->GetOutputPort());
//
//        // 创建一个Actor
//        vtkSmartPointer<vtkActor> redactor = vtkSmartPointer<vtkActor>::New();
//        redactor->SetMapper(redmapper);
//        redactor->GetProperty()->SetOpacity(0.9);
//
//        // 设置颜色
//        redactor->GetProperty()->SetColor(0, 0, 139 / 255.0);  //  cuactor4
//        renderer->AddActor(redactor);
//    //}
//    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//     //// 构建 OBBTree 来进行碰撞检测，假设你已经有了 bone 数据
//    //vtkSmartPointer<vtkOBBTree> obbTree4 = vtkSmartPointer<vtkOBBTree>::New();
//    //obbTree3->SetDataSet(bone);
//    //obbTree3->BuildLocator();
//    //// 调用优化算法并获取优化后的起点位置
//        vtkSmartPointer<vtkPolyData> filteredPointsPolyData4 = ExtractPointsFromTrianglesAndProcess(filteredSkin5, tumorCenter5);
//        vtkSmartPointer<vtkPoints> points4 = filteredPointsPolyData4->GetPoints();
//        // 获取点的数量
//        vtkIdType numPoints4 = points4->GetNumberOfPoints();
//        // 遍历并获取优化后的点
//
//        vtkIdType randomIndex4 = rand() % numPoints4;  // rand() 返回 [0, numPoints-1] 之间的随机数
//
//        double pathStart4[3];
//        points3->GetPoint(randomIndex4, pathStart4);  // 获取每个点的坐标
//        // 在控制台输出随机选中的点坐标
//        std::cout << "红色点坐标: ("
//            << pathStart[0] << ", "
//            << pathStart[1] << ", "
//            << pathStart[2] << ")"
//            << std::endl;
//
//        // 为每个起点创建一个球体
//        vtkSmartPointer<vtkSphereSource> redsphereDataop21 = vtkSmartPointer<vtkSphereSource>::New();
//        redsphereDataop21->SetCenter(pathStart4[0], pathStart4[1], pathStart4[2]);
//        redsphereDataop21->SetRadius(0.5);  // 设置球体半径为 0.5
//
//        // 计算从 tumorCenter4 到 pathStart 的方向向量
//        double dx4 = pathStart4[0] -50 - tumorCenter5[0];
//        double dy4 = pathStart4[1]  -tumorCenter5[1];
//        double dz4 = pathStart4[2] +80 - tumorCenter5[2];
//
//        // 将方向向量单位化
//        normalizeDirection(dx4, dy4, dz4);
//
//        // 计算从 tumorCenter4 开始，超过 pathStart 20 单位的射线终点
//        double pkx = 50;
//        double totalLength4 = distance(tumorCenter5[0], tumorCenter5[1], tumorCenter5[2], pathStart4[0], pathStart4[1], pathStart4[2])-pkx ;
//
//        // 计算射线终点
//        double rayEnd4[3] = {
//            tumorCenter5[0] + dx4 * totalLength4,
//            tumorCenter5[1]+30 + dy4 * totalLength4,
//            tumorCenter5[2] + dz4 * totalLength4
//        };
//
//        vtkSmartPointer<vtkPoints> redpoints21 = vtkSmartPointer<vtkPoints>::New();
//        redpoints21->InsertNextPoint(tumorCenter5); // 设置起点 (tumorCenter4)
//        redpoints21->InsertNextPoint(rayEnd4);       // 设置终点 (rayEnd)
//
//        // 创建一个线源对象
//        vtkSmartPointer<vtkCellArray> redlines21 = vtkSmartPointer<vtkCellArray>::New();
//        redlines21->InsertNextCell(2); // 一条线由两个点组成
//        redlines21->InsertCellPoint(0); // 起点
//        redlines21->InsertCellPoint(1); // 终点
//
//        // 创建一个PolyData对象
//        vtkSmartPointer<vtkPolyData> redlineData21 = vtkSmartPointer<vtkPolyData>::New();
//        redlineData21->SetPoints(redpoints21);
//        redlineData21->SetLines(redlines21);
//        // 使用 vtkTubeFilter 创建圆柱体
//        vtkSmartPointer<vtkTubeFilter>redtubeFilter21 = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilter21->SetInputData(redlineData21);   // 设置线段数据作为输入
//        redtubeFilter21->SetRadius(1);            // 设置圆柱体的半径（相当于“线段的半径”）
//        redtubeFilter21->SetNumberOfSides(50);      // 设置圆柱体的边数，越多越光滑
//
//        // 创建映射器并将球体添加到渲染器中
//        vtkSmartPointer<vtkPolyDataMapper> redsphereMapperop21 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redsphereMapperop21->SetInputConnection(redsphereDataop21->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redsphereActorop21 = vtkSmartPointer<vtkActor>::New();
//        redsphereActorop21->SetMapper(redsphereMapperop21);
//        redsphereActorop21->GetProperty()->SetColor(139 / 255.0, 0, 0);//红色
//        // 将球体添加到渲染器中
//        //renderer->AddActor(sphereActorop2);
//
//        // 创建并可视化优化点与肿瘤质心之间的连线
//        vtkSmartPointer<vtkLineSource> redlineSourceop21 = vtkSmartPointer<vtkLineSource>::New();
//        redlineSourceop21->SetPoint1(pathStart4[0], pathStart4[1], pathStart4[2]);  // 起点
//        redlineSourceop21->SetPoint2(tumorCenter5[0], tumorCenter5[1], tumorCenter5[2]);  // 肿瘤质心
//         // 使用 TubeFilter 增加线段的半径
//        vtkSmartPointer<vtkTubeFilter> redtubeFilterop21 = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilterop21->SetInputConnection(redlineSourceop21->GetOutputPort());
//        redtubeFilterop21->SetRadius(1);  // 设置线段半径
//        redtubeFilterop21->SetNumberOfSides(50);  // 设置线段的精度
//        redtubeFilterop21->Update();
//        // 创建映射器并将线添加到渲染器中
//        vtkSmartPointer<vtkPolyDataMapper> redlineMapperop21 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redlineMapperop21->SetInputConnection(redtubeFilter21->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redlineActorop21 = vtkSmartPointer<vtkActor>::New();
//        redlineActorop21->SetMapper(redlineMapperop21);
//        redlineActorop21->GetProperty()->SetColor(1, 100 / 255.0, 0);  // 设置线的颜色为红色
//
//        // 将线添加到渲染器中
//        renderer->AddActor(redlineActorop21);//------------------------------------------------------------------------------------------------------
//
//        // 计算椭球体长轴（20）和短轴（15）
//        double redsemiMajorAxis21 = 23.0;  // 长轴
//        double redsemiMinorAxis21 = 18.0;  // 短轴
//
//        // 创建一个椭球体的参数化函数
//        vtkSmartPointer<vtkParametricEllipsoid> redellipsoid21 = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        redellipsoid21->SetXRadius(redsemiMinorAxis21);  // 设置长轴
//        redellipsoid21->SetYRadius(redsemiMinorAxis21);  // 设置短轴
//        redellipsoid21->SetZRadius(redsemiMajorAxis21);  // 设置短轴
//
//        // 创建一个参数化函数源
//        vtkSmartPointer<vtkParametricFunctionSource> redellipsoidSource21 = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        redellipsoidSource21->SetParametricFunction(redellipsoid21);
//
//        // 创建椭球体的变换
//        vtkSmartPointer<vtkTransform> redtransform21 = vtkSmartPointer<vtkTransform>::New();
//        // **平移操作**：将椭球体的质心从 (0, 0, 0) 移动到 tumorCenter4
//        redtransform21->Translate(tumorCenter5[0], tumorCenter5[1], tumorCenter5[2]);
//        // 创建旋转矩阵，目标是将椭球体的长轴（默认沿 z 轴）旋转到目标方向（dx, dy, dz）
//        double originalDirection21[3] = { 0.0, 0.0, 1.0 };  // 默认长轴方向
//        double rotationAngle21 = 0.0;
//        double rotationAxis21[3] = { 0.0, 0.0, 0.0 };
//
//        // 计算旋转角度和旋转轴
//        double dotProduct21 = vtkMath::Dot(originalDirection21, new double[3]{ dx4, dy4, dz4 });
//        double crossProduct21[3];
//        vtkMath::Cross(originalDirection21, new double[3]{ dx4, dy4, dz4 }, crossProduct21);
//
//        // 如果叉积为零，表示两个向量已经重合，无需旋转
//        if (vtkMath::Norm(crossProduct21) < 1e-6) {
//            // 如果没有旋转（两个向量重合），则无需做旋转
//            std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//            rotationAngle21 = 0.0;
//        }
//        else {
//            rotationAngle21 = acos(dotProduct21);  // 计算夹角
//            vtkMath::Normalize(crossProduct21);  // 归一化旋转轴
//            std::cout << "Rotation Axis: (" << crossProduct21[0] << ", " << crossProduct21[1] << ", " << crossProduct21[2] << ")" << std::endl;
//        }
//
//        // 旋转椭球体到目标方向
//        redtransform21->RotateWXYZ(rotationAngle21* 180.0 / M_PI, crossProduct21[0], crossProduct21[1], crossProduct21[2]);
//
//
//
//        // 创建一个变换过滤器来应用变换
//        vtkSmartPointer<vtkTransformPolyDataFilter> redtransformFilter21 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        redtransformFilter21->SetInputConnection(redellipsoidSource21->GetOutputPort());
//        redtransformFilter21->SetTransform(redtransform21);
//
//        // 创建一个映射器
//        vtkSmartPointer<vtkPolyDataMapper> redmapper21 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redmapper21->SetInputConnection(redtransformFilter21->GetOutputPort());
//
//        // 创建一个Actor
//        vtkSmartPointer<vtkActor> redactor21 = vtkSmartPointer<vtkActor>::New();
//        redactor21->SetMapper(redmapper21);
//        redactor21->GetProperty()->SetOpacity(0.9);
//
//        // 设置颜色
//        redactor21->GetProperty()->SetColor(1, 100 / 255.0, 0);  //  cuactor4
//        renderer->AddActor(redactor21);
//        //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//        // //// 构建 OBBTree 来进行碰撞检测，假设你已经有了 bone 数据
//    //vtkSmartPointer<vtkOBBTree> obbTree5 = vtkSmartPointer<vtkOBBTree>::New();
//    //obbTree3->SetDataSet(bone);
//    //obbTree3->BuildLocator();
//    //// 调用优化算法并获取优化后的起点位置
//        vtkSmartPointer<vtkPolyData> filteredPointsPolyData5 = ExtractPointsFromTrianglesAndProcess(filteredSkin6, tumorCenter6);
//        vtkSmartPointer<vtkPoints> points5 = filteredPointsPolyData5->GetPoints();
//        // 获取点的数量
//        vtkIdType numPoints5 = points5->GetNumberOfPoints();
//        // 遍历并获取优化后的点
//
//        vtkIdType randomIndex5 = rand() % numPoints5;  // rand() 返回 [0, numPoints-1] 之间的随机数
//
//        double pathStart5[3];
//        points3->GetPoint(randomIndex5, pathStart5);  // 获取每个点的坐标
//        // 在控制台输出随机选中的点坐标
//        std::cout << "红色点坐标: ("
//            << pathStart[0] << ", "
//            << pathStart[1] << ", "
//            << pathStart[2] << ")"
//            << std::endl;
//
//        // 为每个起点创建一个球体
//        vtkSmartPointer<vtkSphereSource> redsphereDataop211 = vtkSmartPointer<vtkSphereSource>::New();
//        redsphereDataop211->SetCenter(pathStart5[0], pathStart5[1], pathStart5[2]);
//        redsphereDataop211->SetRadius(0.5);  // 设置球体半径为 0.5
//
//        // 计算从 tumorCenter4 到 pathStart 的方向向量
//        double dx5= pathStart5[0]  - tumorCenter6[0];
//        double dy5 = pathStart5[1]  - tumorCenter6[1];
//        double dz5 = pathStart5[2] - tumorCenter6[2];
//
//        // 将方向向量单位化
//        normalizeDirection(dx5, dy5, dz5);
//
//        // 计算从 tumorCenter4 开始，超过 pathStart 20 单位的射线终点
//        double o = 50;
//        double totalLength5 = distance(tumorCenter6[0], tumorCenter6[1], tumorCenter6[2], pathStart5[0], pathStart5[1], pathStart5[2])+o;
//
//        // 计算射线终点
//        double rayEnd5[3] = {
//            tumorCenter6[0] + dx5 * totalLength5,
//            tumorCenter6[1] + dy5 * totalLength5,
//            tumorCenter6[2] + dz5 * totalLength5
//        };
//
//        vtkSmartPointer<vtkPoints> redpoints211 = vtkSmartPointer<vtkPoints>::New();
//        redpoints211->InsertNextPoint(tumorCenter6); // 设置起点 (tumorCenter4)
//        redpoints211->InsertNextPoint(rayEnd5);       // 设置终点 (rayEnd)
//
//        // 创建一个线源对象
//        vtkSmartPointer<vtkCellArray> redlines211 = vtkSmartPointer<vtkCellArray>::New();
//        redlines211->InsertNextCell(2); // 一条线由两个点组成
//        redlines211->InsertCellPoint(0); // 起点
//        redlines211->InsertCellPoint(1); // 终点
//
//        // 创建一个PolyData对象
//        vtkSmartPointer<vtkPolyData> redlineData211 = vtkSmartPointer<vtkPolyData>::New();
//        redlineData211->SetPoints(redpoints211);
//        redlineData211->SetLines(redlines211);
//        // 使用 vtkTubeFilter 创建圆柱体
//        vtkSmartPointer<vtkTubeFilter>redtubeFilter211 = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilter211->SetInputData(redlineData211);   // 设置线段数据作为输入
//        redtubeFilter211->SetRadius(1);            // 设置圆柱体的半径（相当于“线段的半径”）
//        redtubeFilter211->SetNumberOfSides(50);      // 设置圆柱体的边数，越多越光滑
//
//        // 创建映射器并将球体添加到渲染器中
//        vtkSmartPointer<vtkPolyDataMapper> redsphereMapperop211 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redsphereMapperop211->SetInputConnection(redsphereDataop211->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redsphereActorop211 = vtkSmartPointer<vtkActor>::New();
//        redsphereActorop211->SetMapper(redsphereMapperop211);
//        redsphereActorop211->GetProperty()->SetColor(139 / 255.0, 0, 0);//红色
//        // 将球体添加到渲染器中
//        //renderer->AddActor(sphereActorop2);
//
//        // 创建并可视化优化点与肿瘤质心之间的连线
//        vtkSmartPointer<vtkLineSource> redlineSourceop211 = vtkSmartPointer<vtkLineSource>::New();
//        redlineSourceop211->SetPoint1(pathStart5[0], pathStart5[1], pathStart5[2]);  // 起点
//        redlineSourceop211->SetPoint2(tumorCenter6[0], tumorCenter6[1], tumorCenter6[2]);  // 肿瘤质心
//         // 使用 TubeFilter 增加线段的半径
//        vtkSmartPointer<vtkTubeFilter> redtubeFilterop211 = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilterop211->SetInputConnection(redlineSourceop211->GetOutputPort());
//        redtubeFilterop211->SetRadius(1);  // 设置线段半径
//        redtubeFilterop211->SetNumberOfSides(50);  // 设置线段的精度
//        redtubeFilterop211->Update();
//        // 创建映射器并将线添加到渲染器中
//        vtkSmartPointer<vtkPolyDataMapper> redlineMapperop211 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redlineMapperop211->SetInputConnection(redtubeFilter211->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redlineActorop211 = vtkSmartPointer<vtkActor>::New();
//        redlineActorop211->SetMapper(redlineMapperop211);
//        redlineActorop211->GetProperty()->SetColor(0, 100 / 255.0, 1);  // 设置线的颜色为红色
//
//        // 将线添加到渲染器中
//        renderer->AddActor(redlineActorop211);//------------------------------------------------------------------------------------------------------
//
//        // 计算椭球体长轴（20）和短轴（15）
//        double redsemiMajorAxis211 = 22.0;  // 长轴
//        double redsemiMinorAxis211 = 17.0;  // 短轴
//
//        // 创建一个椭球体的参数化函数
//        vtkSmartPointer<vtkParametricEllipsoid> redellipsoid211 = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        redellipsoid211->SetXRadius(redsemiMinorAxis211);  // 设置长轴
//        redellipsoid211->SetYRadius(redsemiMinorAxis211);  // 设置短轴
//        redellipsoid211->SetZRadius(redsemiMajorAxis211);  // 设置短轴
//
//        // 创建一个参数化函数源
//        vtkSmartPointer<vtkParametricFunctionSource> redellipsoidSource211 = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        redellipsoidSource211->SetParametricFunction(redellipsoid211);
//
//        // 创建椭球体的变换
//        vtkSmartPointer<vtkTransform> redtransform211 = vtkSmartPointer<vtkTransform>::New();
//        // **平移操作**：将椭球体的质心从 (0, 0, 0) 移动到 tumorCenter4
//        redtransform211->Translate(tumorCenter6[0], tumorCenter6[1], tumorCenter6[2]);
//        // 创建旋转矩阵，目标是将椭球体的长轴（默认沿 z 轴）旋转到目标方向（dx, dy, dz）
//        double originalDirection211[3] = { 0.0, 0.0, 1.0 };  // 默认长轴方向
//        double rotationAngle211 = 0.0;
//        double rotationAxis211[3] = { 0.0, 0.0, 0.0 };
//
//        // 计算旋转角度和旋转轴
//        double dotProduct211 = vtkMath::Dot(originalDirection211, new double[3]{ dx5, dy5, dz5 });
//        double crossProduct211[3];
//        vtkMath::Cross(originalDirection211, new double[3]{ dx5, dy5, dz5 }, crossProduct211);
//
//        // 如果叉积为零，表示两个向量已经重合，无需旋转
//        if (vtkMath::Norm(crossProduct211) < 1e-6) {
//            // 如果没有旋转（两个向量重合），则无需做旋转
//            std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//            rotationAngle211 = 0.0;
//        }
//        else {
//            rotationAngle211 = acos(dotProduct211);  // 计算夹角
//            vtkMath::Normalize(crossProduct211);  // 归一化旋转轴
//            std::cout << "Rotation Axis: (" << crossProduct211[0] << ", " << crossProduct211[1] << ", " << crossProduct211[2] << ")" << std::endl;
//        }
//
//        // 旋转椭球体到目标方向
//        redtransform211->RotateWXYZ(rotationAngle211 * 180.0 / M_PI, crossProduct211[0], crossProduct211[1], crossProduct211[2]);
//
//
//
//        // 创建一个变换过滤器来应用变换
//        vtkSmartPointer<vtkTransformPolyDataFilter> redtransformFilter211 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        redtransformFilter211->SetInputConnection(redellipsoidSource211->GetOutputPort());
//        redtransformFilter211->SetTransform(redtransform211);
//
//        // 创建一个映射器
//        vtkSmartPointer<vtkPolyDataMapper> redmapper211 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redmapper211->SetInputConnection(redtransformFilter211->GetOutputPort());
//
//        // 创建一个Actor
//        vtkSmartPointer<vtkActor> redactor211 = vtkSmartPointer<vtkActor>::New();
//        redactor211->SetMapper(redmapper211);
//        redactor211->GetProperty()->SetOpacity(0.9);
//
//        // 设置颜色
//        redactor211->GetProperty()->SetColor(0, 100 / 255.0, 1);  //  cuactor4
//        renderer->AddActor(redactor211);
////-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//    vtkSmartPointer<vtkActor> filteredSkinActor = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor->SetMapper(filteredSkinMapper);
//    filteredSkinActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
//    filteredSkinActor->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor1 = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor1->SetMapper(filteredSkinMapper1);
//    filteredSkinActor1->GetProperty()->SetColor(0, 0, 139 / 255.0);//蓝色
//    filteredSkinActor1->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor2 = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor2->SetMapper(filteredSkinMapper2);
//    filteredSkinActor2->GetProperty()->SetColor(139 / 255.0, 0, 0);//红色
//    filteredSkinActor2->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor3 = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor3->SetMapper(filteredSkinMapper3);
//    filteredSkinActor3->GetProperty()->SetColor(139 / 255.0, 139 / 255.0, 0);//huang色
//    filteredSkinActor3->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor4 = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor4->SetMapper(filteredSkinMapper4);
//    filteredSkinActor4->GetProperty()->SetColor(0, 100 / 255.0, 0);//lv色
//    filteredSkinActor4->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor5= vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor5->SetMapper(filteredSkinMapper5);
//    filteredSkinActor5->GetProperty()->SetColor(1, 100 / 255.0, 0);//lv色
//    filteredSkinActor5->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor6 = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor6->SetMapper(filteredSkinMapper6);
//    filteredSkinActor6->GetProperty()->SetColor(0, 100 / 255.0, 1);//lv色
//    filteredSkinActor6->GetProperty()->SetOpacity(0.9);
//
//    vtkSmartPointer<vtkActor> BoneActor = vtkSmartPointer<vtkActor>::New();
//    BoneActor->SetMapper(BoneMapper);
//    BoneActor->GetProperty()->SetColor(1, 1, 1);
//    vtkSmartPointer<vtkActor> venoussystemActor = vtkSmartPointer<vtkActor>::New();
//    venoussystemActor->SetMapper(venoussystemMapper);
//    venoussystemActor->GetProperty()->SetColor(1, 0.2, 0);
//
//
//    vtkSmartPointer<vtkActor> liverActor = vtkSmartPointer<vtkActor>::New();
//    liverActor->SetMapper(liverMapper);
//    liverActor->GetProperty()->SetColor(1.0, 0.784, 0.588);
//    liverActor->GetProperty()->SetOpacity(0.7);
//
//    vtkSmartPointer<vtkActor> SkinActor = vtkSmartPointer<vtkActor>::New();
//    SkinActor->SetMapper(SkinMapper);
//    SkinActor->GetProperty()->SetColor(1, 1, 1);
//    SkinActor->GetProperty()->SetOpacity(0.1);
//
//    vtkSmartPointer<vtkActor> tumorActor = vtkSmartPointer<vtkActor>::New();
//    tumorActor->SetMapper(tumorMapper);
//    tumorActor->GetProperty()->SetColor(0.5, 0, 1);
//    tumorActor->GetProperty()->SetOpacity(0.5);
//
//    vtkSmartPointer<vtkActor> ArteryActor = vtkSmartPointer<vtkActor>::New();
//    ArteryActor->SetMapper(ArteryMapper);
//    ArteryActor->GetProperty()->SetColor(1, 0.5, 0);
//
//    // renderer->AddActor(actor);
//     // renderer->AddActor(actor2);
//     // renderer->AddActor(actor3);
//     // renderer->AddActor(actor4);
//      //renderer->AddActor(filteredSkinActor);
//    renderer->AddActor(filteredSkinActor1);
//    renderer->AddActor(filteredSkinActor2);
//    renderer->AddActor(filteredSkinActor3);
//    renderer->AddActor(filteredSkinActor4);
//    renderer->AddActor(filteredSkinActor5);
//    renderer->AddActor(filteredSkinActor6);
//    //renderer->AddActor(cuactor1);
//    //renderer->AddActor(cuactor2);
//    //renderer->AddActor(cuactor3);
//   // renderer->AddActor(cuactor4);
//    renderer->AddActor(SkinActor);
//    renderer->AddActor(venoussystemActor);
//    renderer->AddActor(BoneActor);
//    renderer->AddActor(tumorActor);
//    renderer->AddActor(liverActor);
//    renderer->AddActor(ArteryActor);
//    renderer->SetBackground(1, 1, 1); // 设置背景颜色，double类型
//    renderWindow->Render();
//    renderWindowInteractor->Start();
//
//    return 0;
//}

