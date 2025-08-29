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
//// VTK �ļ���ȡ
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
//// ���������η���
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
//// �������������ļн�
//double ComputeAngleBetweenVectors(double* v1, double* v2) {
//    double dotProduct = vtkMath::Dot(v1, v2);
//    double magnitudeV1 = vtkMath::Norm(v1);
//    double magnitudeV2 = vtkMath::Norm(v2);
//
//    return acos(dotProduct / (magnitudeV1 * magnitudeV2)) * 180.0 / PI;  // ���ؽǶȣ���λΪ��
//}
//
//// ���������������ΰ�Ĥ����С����
//bool CheckMinDistanceToSkin(vtkSmartPointer<vtkPolyData> skin, double* tumorCenter, vtkSmartPointer<vtkPoints> points, const vtkIdType* ptIds) {
//    double minDistance = 1e6;  // ����һ������ĳ�ʼ��С����
//
//    // ��ȡ�����ε���������
//    double p1[3], p2[3], p3[3];
//    points->GetPoint(ptIds[0], p1);
//    points->GetPoint(ptIds[1], p2);
//    points->GetPoint(ptIds[2], p3);
//
//    // ���������εķ�����
//    double normal[3] = { 0.0, 0.0, 0.0 };
//    double v1[3] = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
//    double v2[3] = { p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2] };
//    vtkMath::Cross(v1, v2, normal);
//    vtkMath::Normalize(normal);
//
//    // �����������ĵ�������ƽ��ľ���
//    double vectorToPoint[3] = { tumorCenter[0] - p1[0], tumorCenter[1] - p1[1], tumorCenter[2] - p1[2] };
//    double distance = vtkMath::Dot(vectorToPoint, normal);
//    distance = fabs(distance); // ȡ����ֵ���õ���ƽ��Ĵ�ֱ����
//
//    if (distance < minDistance) {
//        minDistance = distance;  // ������С����
//    }
//
//    // �����С����С��Ԥ��İ�ȫ������ֵ���򷵻�false
//    if (minDistance < 5.0) { // ���谲ȫ����С��5mm
//        return false;
//    }
//
//    return true;
//}
//
//// �ж�������·���Ƿ�����Լ��
//bool CheckConstraints(vtkSmartPointer<vtkPolyData> skin, double* tumorCenter, vtkSmartPointer<vtkOBBTree> obbTree, const vtkIdType* ptIds, vtkSmartPointer<vtkPoints> points) {
//    double avgDistance = 0.0;
//    vtkSmartPointer<vtkPoints> intersectPoints = vtkSmartPointer<vtkPoints>::New();
//
//    // 1. ��ȫ����Լ��
//    for (int i = 0; i < 3; i++) {
//        double point[3];
//        points->GetPoint(ptIds[i], point);
//
//        // ���������������ΰ�Ĥ����С����
//        if (obbTree->IntersectWithLine(tumorCenter, point, intersectPoints, nullptr)) {
//            return false; // ·����ΰ�Ĥ�н��㣬Υ���˰�ȫ����
//        }
//
//        double distance = sqrt(pow(point[0] - tumorCenter[0], 2) +
//            pow(point[1] - tumorCenter[1], 2) +
//            pow(point[2] - tumorCenter[2], 2));
//        avgDistance += distance / 3.0;
//    }
//
//    // ��ȫ����Լ��: ����������浽·����ƽ�������Ƿ�С��5mm�����60mm
//    if (avgDistance < 5.0 || avgDistance > 60.0) {
//        return false;
//    }
//
//    // 2. ���߽Ƕ�Լ��
//    // ���㽻�㸽�������εķ���
//    double normal[3] = { 0.0, 0.0, 0.0 };
//    for (int i = 0; i < 3; i++) {
//        double point[3];
//        points->GetPoint(ptIds[i], point);
//        // ��ȡ���㸽�������εķ���
//        ComputeTriangleNormal(skin->GetPoints(), ptIds, normal);
//
//        // ����������·����ΰ�Ĥ���淨��֮��ļн�
//        double direction[3] = { point[0] - tumorCenter[0], point[1] - tumorCenter[1], point[2] - tumorCenter[2] };
//        double angle = ComputeAngleBetweenVectors(direction, normal);
//
//        if (angle < 20.0) { // ���߽Ƕ�С��20�㣬������Լ��
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
//    // �� skin ������ת��Ϊ������
//    vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
//    triangleFilter->SetInputData(skin);
//    triangleFilter->Update();
//
//    // ���������������;����� OBBTree����������ײ���ʹ��
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
//    // �洢ɸѡ���������
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
//            // ��ȡ�������������������
//            double p0[3], p1[3], p2[3];
//            points->GetPoint(ptIds[0], p0);
//            points->GetPoint(ptIds[1], p1);
//            points->GetPoint(ptIds[2], p2);
//
//            // ��������������
//            double centroid[3] = { (p0[0] + p1[0] + p2[0]) / 3.0,
//                                   (p0[1] + p1[1] + p2[1]) / 3.0,
//                                   (p0[2] + p1[2] + p2[2]) / 3.0 };
//
//            // �����������ĵ����������ĵľ���
//            double dist2 = vtkMath::Distance2BetweenPoints(tumorCenter, centroid);
//            // ����������С��150��������
//            if (dist2 < 45.0 * 45.0) {
//                // �������Լ������
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
//// ���� Whale ��洢���Ż������е�ÿ�����λ�ú���Ӧ��
//class Whale {
//public:
//    std::vector<double> position; // �洢����λ�õ����� (x, y, z)
//    std::vector<double> fitness;  // �洢�������Ӧ��ֵ
//};
//
//// ��������·���ĳ��ȣ�·�����������ʾ��ֱ�ߣ�
//double CalculatePathLength(const double* point1, const double* point2) {
//    return vtkMath::Distance2BetweenPoints(point1, point2);  // �����ƽ��
//}
//
//// ����·����Ѫ��/��ͷ����С����
//double CalculatePathDistanceToStructure(const double* pathStart, const double* pathEnd, vtkSmartPointer<vtkOBBTree> obbTree) {
//    // ����·�����ͷ��Ѫ�ܵ���С����
//    // ʹ�� OBBTree ����·����Ѫ��/��ͷ֮��ľ��루�����Ż�·���㵽��ͷ��Ѫ�ܵľ��룩
//    double minDistance = std::numeric_limits<double>::max();
//    // �����ͨ���������������ȡ·���ͽṹ֮�����С���룬�������������С����Ϊģ��
//    return minDistance;
//}
//
//// ����·�����������ĵ�ƽ�ж�
//double CalculatePathAlignmentWithTumorCenter(const double* pathStart, const double* pathEnd, const double* tumorCenter) {
//    // ����·���ķ�������
//    double pathVector[3] = { pathEnd[0] - pathStart[0], pathEnd[1] - pathStart[1], pathEnd[2] - pathStart[2] };
//
//    // �����������ĵķ�������
//    double tumorVector[3] = { tumorCenter[0] - pathStart[0], tumorCenter[1] - pathStart[1], tumorCenter[2] - pathStart[2] };
//
//    // �������������ĵ�������Խ�󣬱�ʾԽƽ��
//    double dotProduct = vtkMath::Dot(pathVector, tumorVector);
//    double pathLength = vtkMath::Norm(pathVector);
//    double tumorLength = vtkMath::Norm(tumorVector);
//
//    // ����ƽ�жȣ����Խ�󣬷���Խ����
//    double alignment = dotProduct / (pathLength * tumorLength);
//    return alignment;
//}
//
//// ���������в����㣨����ѡ�����������ģ�
//void SamplePointsFromTriangles(vtkSmartPointer<vtkPolyData> filteredSkin, std::vector<std::vector<double>>& sampledPoints) {
//    vtkSmartPointer<vtkPoints> points = filteredSkin->GetPoints();
//    vtkSmartPointer<vtkCellArray> polys = filteredSkin->GetPolys();
//
//    vtkIdType npts;
//    const vtkIdType* ptIds;
//    polys->InitTraversal();
//
//    while (polys->GetNextCell(npts, ptIds)) {
//        if (npts == 3) { // �����������
//            // ��ȡ�������������������
//            double p0[3], p1[3], p2[3];
//            points->GetPoint(ptIds[0], p0);
//            points->GetPoint(ptIds[1], p1);
//            points->GetPoint(ptIds[2], p2);
//
//            // ��������������
//            double centroid[3] = { (p0[0] + p1[0] + p2[0]) / 3.0,
//                                   (p0[1] + p1[1] + p2[1]) / 3.0,
//                                   (p0[2] + p1[2] + p2[2]) / 3.0 };
//
//            // �洢����
//            sampledPoints.push_back({ centroid[0], centroid[1], centroid[2] });
//        }
//    }
//}
//
////// ��Ŀ�꾨���Ż��㷨
////std::vector<std::vector<double>> WhaleOptimizationAlgorithm(vtkSmartPointer<vtkPolyData> filteredSkin, vtkSmartPointer<vtkOBBTree> obbTree, double* tumorCenter) {
////    // �����Ż�����Ⱥ����Ϊ 10��ÿ�������� 3 ��λ�ò��� (x, y, z)
////    int populationSize = 1;
////    std::vector<Whale> population(populationSize);
////    std::vector<std::vector<double>> optimizedPoints;  // �洢�Ż���ĵ㣨��ʼ�㣩
////
////    // ��ɸѡ���������Ƭ�в�����
////    std::vector<std::vector<double>> sampledPoints;
////    SamplePointsFromTriangles(filteredSkin, sampledPoints);
////
////    // ��ʼ����Ⱥλ�� (��·������ʼ�����ֹ��)
////    for (int i = 0; i < populationSize; i++) {
////        // ���ѡ������ĵ���Ϊ·������ʼ��
////        int startIndex = rand() % sampledPoints.size();
////        int endIndex = rand() % sampledPoints.size();
////
////        // �洢·������ʼ�����ֹ��
////        population[i].position = { sampledPoints[startIndex][0], sampledPoints[startIndex][1], sampledPoints[startIndex][2] };
////        population[i].position.push_back(sampledPoints[endIndex][0]);
////        population[i].position.push_back(sampledPoints[endIndex][1]);
////        population[i].position.push_back(sampledPoints[endIndex][2]);
////    }
////
////    // ��ÿ����������Ż���������Ӧ�ȣ�
////    for (int i = 0; i < populationSize; i++) {
////        double pathStart[3] = { population[i].position[0], population[i].position[1], population[i].position[2] };
////        double pathEnd[3] = { population[i].position[3], population[i].position[4], population[i].position[5] };
////
////        // ����Ŀ��1��·��������С��
////        double pathLength = CalculatePathLength(pathStart, pathEnd);
////        population[i].fitness.push_back(pathLength);  // Ŀ��1��·��������С��
////
////        // ����Ŀ��2��·�����ͷ/Ѫ�ܵľ������
////        double pathDistanceToStructure = CalculatePathDistanceToStructure(pathStart, pathEnd, obbTree);
////        population[i].fitness.push_back(-pathDistanceToStructure);  // Ŀ��2�����·����Ѫ�ܹ�ͷ�ľ���
////
////        // ����Ŀ��3��·��ƽ�ж����
////        double pathAlignment = CalculatePathAlignmentWithTumorCenter(pathStart, pathEnd, tumorCenter);
////        population[i].fitness.push_back(pathAlignment);  // Ŀ��3�����·��ƽ�ж�
////    }
////
////    // ��ȡ·������ʼ�㵽�������
////    for (int i = 0; i < populationSize; i++) {
////        double pathStart[3] = { population[i].position[0], population[i].position[1], population[i].position[2] };
////        optimizedPoints.push_back({ pathStart[0], pathStart[1], pathStart[2] });
////    }
////
////    return optimizedPoints;  // �����Ż�������
////}
//// ����������֮��ľ���
//double distance(double x1, double y1, double z1, double x2, double y2, double z2) {
//    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
//}
//
//// ���㵥λ��������
//void normalizeDirection(double& dx, double& dy, double& dz) {
//    double len = sqrt(dx * dx + dy * dy + dz * dz);
//    dx /= len;
//    dy /= len;
//    dz /= len;
//}
//vtkSmartPointer<vtkPolyData> ExtractPointsFromTrianglesAndProcess(vtkSmartPointer<vtkPolyData> skin, double* tumorCenter) {
//    // �� skin ������ת��Ϊ������
//    vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
//    triangleFilter->SetInputData(skin);
//    triangleFilter->Update();
//
//    // �洢ɸѡ��ĵ�
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
//            // ��ȡ�������������������
//            double p0[3], p1[3], p2[3];
//            points->GetPoint(ptIds[0], p0);
//            points->GetPoint(ptIds[1], p1);
//            points->GetPoint(ptIds[2], p2);
//
//            // ��������������
//            double centroid[3] = { (p0[0] + p1[0] + p2[0]) / 3.0,
//                                   (p0[1] + p1[1] + p2[1]) / 3.0,
//                                   (p0[2] + p1[2] + p2[2]) / 3.0 };
//
//            // �����������ĵ����������ĵľ���
//            double dist2 = vtkMath::Distance2BetweenPoints(tumorCenter, centroid);
//            // ����������С��150��������
//            if (dist2 < 45.0 * 45.0) {
//                // �������εĶ�����ӵ�ɸѡ��ĵ㼯
//                filteredPoints->InsertNextPoint(p0);
//                filteredPoints->InsertNextPoint(p1);
//                filteredPoints->InsertNextPoint(p2);
//            }
//        }
//    }
//
//    // ���õ�����
//    filteredPointsPolyData->SetPoints(filteredPoints);
//
//    // ��ȡ�������㼯
//    vtkSmartPointer<vtkPoints> finalPoints = filteredPointsPolyData->GetPoints();
//
//    // �����㼯
//    for (vtkIdType i = 0; i < finalPoints->GetNumberOfPoints(); i++) {
//        double pathStart[3];
//        finalPoints->GetPoint(i, pathStart);  // ��ȡÿ���������
//
//        // ʹ�� pathStart ���к�������
//        // ���磬����������ʹ�� pathStart ����·���滮����������
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
//    // ���� vtkDelimitedTextReader ��ȡ CSV �ļ�
//    vtkSmartPointer<vtkDelimitedTextReader> reader = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader->SetFileName(fileName.c_str());
//    reader->SetHaveHeaders(true); // ��� CSV �ļ���ͷ����Ϣ
//    reader->SetFieldDelimiterCharacters(","); // ʹ�ö�����Ϊ�ָ���
//    reader->Update();
//    vtkSmartPointer<vtkDelimitedTextReader> reader1 = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader1->SetFileName(fileName1.c_str());
//    reader1->SetHaveHeaders(true); // ��� CSV �ļ���ͷ����Ϣ
//    reader1->SetFieldDelimiterCharacters(","); // ʹ�ö�����Ϊ�ָ���
//    reader1->Update();
//    vtkSmartPointer<vtkDelimitedTextReader> reader2 = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader2->SetFileName(fileName2.c_str());
//    reader2->SetHaveHeaders(true); // ��� CSV �ļ���ͷ����Ϣ
//    reader2->SetFieldDelimiterCharacters(","); // ʹ�ö�����Ϊ�ָ���
//    reader2->Update();
//    vtkSmartPointer<vtkDelimitedTextReader> reader3 = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader3->SetFileName(fileName3.c_str());
//    reader3->SetHaveHeaders(true); // ��� CSV �ļ���ͷ����Ϣ
//    reader3->SetFieldDelimiterCharacters(","); // ʹ�ö�����Ϊ�ָ���
//    reader3->Update();
//    vtkSmartPointer<vtkDelimitedTextReader> reader4 = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader4->SetFileName(fileName4.c_str());
//    reader4->SetHaveHeaders(true); // ��� CSV �ļ���ͷ����Ϣ
//    reader4->SetFieldDelimiterCharacters(","); // ʹ�ö�����Ϊ�ָ���
//    reader4->Update();
//    vtkSmartPointer<vtkDelimitedTextReader> reader5 = vtkSmartPointer<vtkDelimitedTextReader>::New();
//    reader5->SetFileName(fileName5.c_str());
//    reader5->SetHaveHeaders(true); // ��� CSV �ļ���ͷ����Ϣ
//    reader5->SetFieldDelimiterCharacters(","); // ʹ�ö�����Ϊ�ָ���
//    reader5->Update();
//
//    // ��ȡ�����ݽ��洢�� vtkTable ��
//    vtkSmartPointer<vtkTable> table = reader->GetOutput();
//    vtkSmartPointer<vtkTable> table1 = reader1->GetOutput();
//    vtkSmartPointer<vtkTable> table2 = reader2->GetOutput();
//    vtkSmartPointer<vtkTable> table3 = reader3->GetOutput();
//    vtkSmartPointer<vtkTable> table4 = reader4->GetOutput();
//    vtkSmartPointer<vtkTable> table5 = reader5->GetOutput();
//
//
//    // ���� vtkPoints ���洢 CSV �ļ��еĵ�����
//    vtkSmartPointer<vtkPoints> cupoints1 = vtkSmartPointer<vtkPoints>::New();
//    // ���� CSV �ļ������� x, y, z
//    for (vtkIdType i = 0; i < table->GetNumberOfRows(); ++i) {
//        double x = table->GetValue(i, 0).ToDouble(); // ��ȡ x ����
//        double y = table->GetValue(i, 1).ToDouble(); // ��ȡ y ����
//        double z = table->GetValue(i, 2).ToDouble(); // ��ȡ z ����
//
//        // ������ӵ� vtkPoints ��
//        cupoints1->InsertNextPoint(x, y, z);
//    }
//    // ���� vtkPoints ���洢 CSV �ļ��еĵ�����
//    vtkSmartPointer<vtkPoints> cupoints2 = vtkSmartPointer<vtkPoints>::New();
//    // ���� CSV �ļ������� x, y, z
//    for (vtkIdType i = 0; i < table1->GetNumberOfRows(); ++i) {
//        double x = table1->GetValue(i, 0).ToDouble(); // ��ȡ x ����
//        double y = table1->GetValue(i, 1).ToDouble(); // ��ȡ y ����
//        double z = table1->GetValue(i, 2).ToDouble(); // ��ȡ z ����
//
//        // ������ӵ� vtkPoints ��
//        cupoints2->InsertNextPoint(x, y, z);
//    }
//    // ���� vtkPoints ���洢 CSV �ļ��еĵ�����
//    vtkSmartPointer<vtkPoints> cupoints3 = vtkSmartPointer<vtkPoints>::New();
//    // ���� CSV �ļ������� x, y, z
//    for (vtkIdType i = 0; i < table2->GetNumberOfRows(); ++i) {
//        double x = table2->GetValue(i, 0).ToDouble(); // ��ȡ x ����
//        double y = table2->GetValue(i, 1).ToDouble(); // ��ȡ y ����
//        double z = table2->GetValue(i, 2).ToDouble(); // ��ȡ z ����
//
//        // ������ӵ� vtkPoints ��
//        cupoints3->InsertNextPoint(x, y, z);
//    }
//    // ���� vtkPoints ���洢 CSV �ļ��еĵ�����
//    vtkSmartPointer<vtkPoints> cupoints4 = vtkSmartPointer<vtkPoints>::New();
//    // ���� CSV �ļ������� x, y, z
//    for (vtkIdType i = 0; i < table3->GetNumberOfRows(); ++i) {
//        double x = table3->GetValue(i, 0).ToDouble(); // ��ȡ x ����
//        double y = table3->GetValue(i, 1).ToDouble(); // ��ȡ y ����
//        double z = table3->GetValue(i, 2).ToDouble(); // ��ȡ z ����
//
//        // ������ӵ� vtkPoints ��
//        cupoints4->InsertNextPoint(x, y, z);
//    }
//    // ���� vtkPoints ���洢 CSV �ļ��еĵ�����
//    vtkSmartPointer<vtkPoints> cupoints5 = vtkSmartPointer<vtkPoints>::New();
//    // ���� CSV �ļ������� x, y, z
//    for (vtkIdType i = 0; i < table4->GetNumberOfRows(); ++i) {
//        double x = table4->GetValue(i, 0).ToDouble(); // ��ȡ x ����
//        double y = table4->GetValue(i, 1).ToDouble(); // ��ȡ y ����
//        double z = table4->GetValue(i, 2).ToDouble(); // ��ȡ z ����
//
//        // ������ӵ� vtkPoints ��
//        cupoints5->InsertNextPoint(x, y, z);
//    }
//    // ���� vtkPoints ���洢 CSV �ļ��еĵ�����
//    vtkSmartPointer<vtkPoints> cupoints6= vtkSmartPointer<vtkPoints>::New();
//    // ���� CSV �ļ������� x, y, z
//    for (vtkIdType i = 0; i < table5->GetNumberOfRows(); ++i) {
//        double x = table5->GetValue(i, 0).ToDouble(); // ��ȡ x ����
//        double y = table5->GetValue(i, 1).ToDouble(); // ��ȡ y ����
//        double z = table5->GetValue(i, 2).ToDouble(); // ��ȡ z ����
//
//        // ������ӵ� vtkPoints ��
//        cupoints6->InsertNextPoint(x, y, z);
//    }
//    // 3. ʹ�� vtkSphereSource ��������
//    vtkSmartPointer<vtkSphereSource> cusphereSource1 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource1->SetRadius(0.2);  // ��������İ뾶
//    cusphereSource1->SetPhiResolution(10);  // ��������ֱ���
//    cusphereSource1->SetThetaResolution(10);  // ���ú���ֱ���
//    cusphereSource1->Update();  // ��������Դ
//    vtkSmartPointer<vtkSphereSource> cusphereSource2 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource2->SetRadius(0.2);  // ��������İ뾶
//    cusphereSource2->SetPhiResolution(10);  // ��������ֱ���
//    cusphereSource2->SetThetaResolution(10);  // ���ú���ֱ���
//    cusphereSource2->Update();  // ��������Դ
//    vtkSmartPointer<vtkSphereSource> cusphereSource3 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource3->SetRadius(0.2);  // ��������İ뾶
//    cusphereSource3->SetPhiResolution(10);  // ��������ֱ���
//    cusphereSource3->SetThetaResolution(10);  // ���ú���ֱ���
//    cusphereSource3->Update();  // ��������Դ
//    vtkSmartPointer<vtkSphereSource> cusphereSource4 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource4->SetRadius(0.2);  // ��������İ뾶
//    cusphereSource4->SetPhiResolution(10);  // ��������ֱ���
//    cusphereSource4->SetThetaResolution(10);  // ���ú���ֱ���
//    cusphereSource4->Update();  // ��������Դ
//    vtkSmartPointer<vtkSphereSource> cusphereSource5 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource5->SetRadius(0.2);  // ��������İ뾶
//    cusphereSource5->SetPhiResolution(10);  // ��������ֱ���
//    cusphereSource5->SetThetaResolution(10);  // ���ú���ֱ���
//    cusphereSource5->Update();  // ��������Դ
//    vtkSmartPointer<vtkSphereSource> cusphereSource6 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource6->SetRadius(0.2);  // ��������İ뾶
//    cusphereSource6->SetPhiResolution(10);  // ��������ֱ���
//    cusphereSource6->SetThetaResolution(10);  // ���ú���ֱ���
//    cusphereSource6->Update();  // ��������Դ
//    vtkSmartPointer<vtkPolyData> cupolyData1 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData1->SetPoints(cupoints1); // ���õ�����
//    vtkSmartPointer<vtkPolyData> cupolyData2 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData2->SetPoints(cupoints2); // ���õ�����
//    vtkSmartPointer<vtkPolyData> cupolyData3 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData3->SetPoints(cupoints3); // ���õ�����
//    vtkSmartPointer<vtkPolyData> cupolyData4 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData4->SetPoints(cupoints4); // ���õ�����
//    vtkSmartPointer<vtkPolyData> cupolyData5 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData5->SetPoints(cupoints5); // ���õ�����
//    vtkSmartPointer<vtkPolyData> cupolyData6 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData6->SetPoints(cupoints6); // ���õ�����
//    // 4. ʹ�� vtkGlyph3D ��ÿ����ת��Ϊ����
//    vtkSmartPointer<vtkGlyph3D> glyph3D1 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D1->SetInputData(cupolyData1);  // ���������
//    glyph3D1->SetSourceConnection(cusphereSource1->GetOutputPort());  // ��������Դ
//    glyph3D1->Update();  // �������ɵĽ��
//    vtkSmartPointer<vtkGlyph3D> glyph3D2 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D2->SetInputData(cupolyData2);  // ���������
//    glyph3D2->SetSourceConnection(cusphereSource2->GetOutputPort());  // ��������Դ
//    glyph3D2->Update();  // �������ɵĽ��
//    vtkSmartPointer<vtkGlyph3D> glyph3D3 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D3->SetInputData(cupolyData3);  // ���������
//    glyph3D3->SetSourceConnection(cusphereSource3->GetOutputPort());  // ��������Դ
//    glyph3D3->Update();  // �������ɵĽ��
//    vtkSmartPointer<vtkGlyph3D> glyph3D4 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D4->SetInputData(cupolyData4);  // ���������
//    glyph3D4->SetSourceConnection(cusphereSource4->GetOutputPort());  // ��������Դ
//    glyph3D4->Update();  // �������ɵĽ��
//    vtkSmartPointer<vtkGlyph3D> glyph3D5 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D5->SetInputData(cupolyData5);  // ���������
//    glyph3D5->SetSourceConnection(cusphereSource5->GetOutputPort());  // ��������Դ
//    glyph3D5->Update();  // �������ɵĽ��
//
//    vtkSmartPointer<vtkGlyph3D> glyph3D6 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph3D6->SetInputData(cupolyData6);  // ���������
//    glyph3D6->SetSourceConnection(cusphereSource6->GetOutputPort());  // ��������Դ
//    glyph3D6->Update();  // �������ɵĽ��
//
//
//     // ����һ�� Mapper ����Ⱦ������
//    vtkSmartPointer<vtkPolyDataMapper> cumapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    cumapper1->SetInputConnection(glyph3D1->GetOutputPort());
//    // ����һ�� Mapper ����Ⱦ������
//    vtkSmartPointer<vtkPolyDataMapper> cumapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    cumapper2->SetInputConnection(glyph3D2->GetOutputPort());
//
//    // ����һ�� Mapper ����Ⱦ������
//    vtkSmartPointer<vtkPolyDataMapper> cumapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    cumapper3->SetInputConnection(glyph3D3->GetOutputPort());
//    // ����һ�� Mapper ����Ⱦ������
//    vtkSmartPointer<vtkPolyDataMapper> cumapper4 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    cumapper4->SetInputConnection(glyph3D4->GetOutputPort());
//    // ����һ�� Mapper ����Ⱦ������
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
//    // ���������壨20, 15, 15 ������
//    double a = 29;  // ����
//    double b = 22;  // ����
//    double c = 22;  // ����
//    vtkSmartPointer<vtkParametricEllipsoid> ellipsoid = vtkSmartPointer<vtkParametricEllipsoid>::New();
//    ellipsoid->SetXRadius(a);
//    ellipsoid->SetYRadius(b);
//    ellipsoid->SetZRadius(c);
//
//    vtkSmartPointer<vtkParametricFunctionSource> source = vtkSmartPointer<vtkParametricFunctionSource>::New();
//    source->SetParametricFunction(ellipsoid);
//    source->Update();
//    // ʹ�����ǻ��㷨��������ı���
//    vtkSmartPointer<vtkTriangleFilter> triFilter = vtkSmartPointer<vtkTriangleFilter>::New();
//    triFilter->SetInputConnection(source->GetOutputPort());
//    triFilter->Update();
//
//    // �����������ӳ����
//    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    mapper->SetInputConnection(triFilter->GetOutputPort());
//    // �����������ӳ����
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
//    // ������������ӵ� vtkPoints ��
//    cu1->InsertNextPoint(tumorCenter1);
//    // ����������� Actor������������ɫΪ��ɫ
//    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//    actor->SetMapper(mapper);
//    actor->GetProperty()->SetColor(1, 0, 0);  // ������ɫΪ��ɫ
//    actor->SetPosition(tumorCenter1);  // ����������õ���������λ��
//    actor->GetProperty()->SetOpacity(0.8);
//
//    // ������������ӵ� vtkPoints ��
//    cu2->InsertNextPoint(tumorCenter2);
//    // ����������� Actor������������ɫΪ��ɫ
//    vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
//    actor2->SetMapper(mapper2);
//    actor2->GetProperty()->SetColor(1, 1, 0);  // ������ɫΪ��ɫ
//    actor2->SetPosition(tumorCenter2);  // ����������õ���������λ��
//    actor2->GetProperty()->SetOpacity(0.8);
//    // ������������ӵ� vtkPoints ��
//    cu3->InsertNextPoint(tumorCenter3);
//    // ����������� Actor������������ɫΪ��ɫ
//    vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
//    actor3->SetMapper(mapper3);
//    actor3->GetProperty()->SetColor(0, 1, 0);  // ������ɫΪ��ɫ
//    actor3->SetPosition(tumorCenter3);  // ����������õ���������λ��
//    actor3->GetProperty()->SetOpacity(0.8);
//    // ������������ӵ� vtkPoints ��
//    cu4->InsertNextPoint(tumorCenter4);
//    // ����������� Actor������������ɫΪ��ɫ
//    vtkSmartPointer<vtkActor> actor4 = vtkSmartPointer<vtkActor>::New();
//    actor4->SetMapper(mapper4);
//    actor4->GetProperty()->SetColor(0, 0, 1);  // ������ɫΪ��ɫ
//    actor4->SetPosition(tumorCenter4);  // ����������õ���������λ��
//    actor4->GetProperty()->SetOpacity(0.8);
//
//    // ����һ�� Actor ����ʾ��
//    vtkSmartPointer<vtkActor> cuactor1 = vtkSmartPointer<vtkActor>::New();
//    cuactor1->SetMapper(cumapper1);
//    cuactor1->GetProperty()->SetColor(1, 1, 0);//huangɫ
//    cuactor1->GetProperty()->SetPointSize(2.5);
//    // ����һ�� Actor ����ʾ��
//    vtkSmartPointer<vtkActor> cuactor2 = vtkSmartPointer<vtkActor>::New();
//    cuactor2->SetMapper(cumapper2);
//    cuactor2->GetProperty()->SetColor(1, 0, 0);//��ɫ
//    cuactor2->GetProperty()->SetPointSize(2.5);
//    // ����һ�� Actor ����ʾ��
//    vtkSmartPointer<vtkActor> cuactor3 = vtkSmartPointer<vtkActor>::New();
//    cuactor3->SetMapper(cumapper3);
//    cuactor3->GetProperty()->SetColor(0.0, 1.0, 0.0);//lv se
//    cuactor3->GetProperty()->SetPointSize(2.5);
//    // ����һ�� Actor ����ʾ��
//    vtkSmartPointer<vtkActor> cuactor4 = vtkSmartPointer<vtkActor>::New();
//    cuactor4->SetMapper(cumapper4);
//    cuactor4->GetProperty()->SetColor(1.0, 0.0, 0.0);//red
//    cuactor4->GetProperty()->SetPointSize(2.5);
//    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//    // ���� OBBTree ��������ײ��⣬�������Ѿ����� bone ����
//    vtkSmartPointer<vtkOBBTree> obbTree1 = vtkSmartPointer<vtkOBBTree>::New();
//    obbTree1->SetDataSet(bone);
//    obbTree1->BuildLocator();
//    // �����Ż��㷨����ȡ�Ż�������λ��
//// ���� ExtractPointsFromTrianglesAndProcess ����������ȡ�㼯�ķ���
//    vtkSmartPointer<vtkPolyData> filteredPointsPolyData1 = ExtractPointsFromTrianglesAndProcess(filteredSkin4, tumorCenter4);
//    vtkSmartPointer<vtkPoints> points1 = filteredPointsPolyData1->GetPoints();
//    // ��ȡ�������
//    vtkIdType numPoints = points1->GetNumberOfPoints();
//    // ��������ȡ�Ż���ĵ�
//
//        vtkIdType randomIndex = rand() % numPoints;  // rand() ���� [0, numPoints-1] ֮��������
//
//        double pathStart[3];
//        points1->GetPoint(randomIndex, pathStart);  // ��ȡÿ���������
//        // �ڿ���̨������ѡ�еĵ�����
//        std::cout << "��ɫ������: ("
//            << pathStart[0] << ", "
//            << pathStart[1] << ", "
//            << pathStart[2] << ")"
//            << std::endl;
//
//        // Ϊÿ����㴴��һ������
//        vtkSmartPointer<vtkSphereSource> GreensphereDataop = vtkSmartPointer<vtkSphereSource>::New();
//        GreensphereDataop->SetCenter(pathStart[0], pathStart[1], pathStart[2]);
//       /* std::cout << "pathStart[0]: " << pathStart[0] << std::endl;
//        std::cout << "pathStart[1]: " << pathStart[1] << std::endl;
//        std::cout << "pathStart[2]: " << pathStart[2] << std::endl;*/
//
//        GreensphereDataop->SetRadius(0.5);  // ��������뾶Ϊ 0.5
//            // ����� tumorCenter4 �� pathStart �ķ�������
//        double dx = pathStart[0] - tumorCenter4[0];
//        double dy = pathStart[1] - tumorCenter4[1];
//        double dz = pathStart[2] - tumorCenter4[2];
//
//        // ������������λ��
//        normalizeDirection(dx, dy, dz);
//
//        // ����� tumorCenter4 ��ʼ������ pathStart 20 ��λ�������յ�
//        double extraLength = 50.0;
//        double totalLength = distance(tumorCenter4[0], tumorCenter4[1], tumorCenter4[2], pathStart[0], pathStart[1], pathStart[2]) + extraLength;
//
//        // ���������յ�
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
//        // VTK���ӻ�����
//
//        // ����һ����Դ������������ʼ��
//        vtkSmartPointer<vtkPoints> Greenpoints = vtkSmartPointer<vtkPoints>::New();
//        Greenpoints->InsertNextPoint(tumorCenter4); // ������� (tumorCenter4)
//        Greenpoints->InsertNextPoint(rayEnd);       // �����յ� (rayEnd)
//
//        // ����һ����Դ����
//        vtkSmartPointer<vtkCellArray> Greenlines = vtkSmartPointer<vtkCellArray>::New();
//        Greenlines->InsertNextCell(2); // һ���������������
//        Greenlines->InsertCellPoint(0); // ���
//        Greenlines->InsertCellPoint(1); // �յ�
//
//        // ����һ��PolyData����
//        vtkSmartPointer<vtkPolyData> GreenlineData = vtkSmartPointer<vtkPolyData>::New();
//        GreenlineData->SetPoints(Greenpoints);
//        GreenlineData->SetLines(Greenlines);
//        // ʹ�� vtkTubeFilter ����Բ����
//        vtkSmartPointer<vtkTubeFilter>GreentubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
//        GreentubeFilter->SetInputData(GreenlineData);   // �����߶�������Ϊ����
//        GreentubeFilter->SetRadius(1.0);            // ����Բ����İ뾶���൱�ڡ��߶εİ뾶����
//        GreentubeFilter->SetNumberOfSides(50);      // ����Բ����ı�����Խ��Խ�⻬
//        // ����һ��ӳ��������ʾ��
//        vtkSmartPointer<vtkPolyDataMapper> GreenlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        GreenlineMapper->SetInputConnection(GreentubeFilter->GetOutputPort());
//
//        // ����һ��Actor����Ⱦ��
//        vtkSmartPointer<vtkActor> GreenlineActor = vtkSmartPointer<vtkActor>::New();
//        GreenlineActor->SetMapper(GreenlineMapper);
//
//        // ������������ɫ
//        GreenlineActor->GetProperty()->SetColor(0, 1, 0);  // ��ɫ
//
//        renderer->AddActor(GreenlineActor);//-------------------------------------------------------------------------------------------------------------------------------------
//        // ���������峤�ᣨ20���Ͷ��ᣨ15��
//        double GreensemiMajorAxis = 22.0;  // ����
//        double GreensemiMinorAxis = 17;  // ����
//
//        // ����һ��������Ĳ���������
//        vtkSmartPointer<vtkParametricEllipsoid> Greenellipsoid = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        Greenellipsoid->SetXRadius(GreensemiMinorAxis);  // ���ó���
//        Greenellipsoid->SetYRadius(GreensemiMinorAxis);  // ���ö���
//        Greenellipsoid->SetZRadius(GreensemiMajorAxis);  // ���ö���
//
//        // ����һ������������Դ
//        vtkSmartPointer<vtkParametricFunctionSource> GreenellipsoidSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        GreenellipsoidSource->SetParametricFunction(Greenellipsoid);
//        double  M_PI = 3.14;
//        // ����������ı任
//        vtkSmartPointer<vtkTransform> Greentransform = vtkSmartPointer<vtkTransform>::New();
//        // **ƽ�Ʋ���**��������������Ĵ� (0, 0, 0) �ƶ��� tumorCenter4
//        Greentransform->Translate(tumorCenter4[0], tumorCenter4[1], tumorCenter4[2]);
//
//        // ������ת����Ŀ���ǽ�������ĳ��ᣨĬ���� z �ᣩ��ת��Ŀ�귽��dx, dy, dz��
//        double originalDirection[3] = { 0.0, 0.0, 1.0 };  // Ĭ�ϳ��᷽��
//        double rotationAngle = 0.0;
//        double rotationAxis[3] = { 0.0, 0.0, 0.0 };
//
//        // ������ת�ǶȺ���ת��
//        double dotProduct = vtkMath::Dot(originalDirection, new double[3]{ dx, dy, dz });
//        double crossProduct[3];
//        vtkMath::Cross(originalDirection, new double[3]{ dx, dy, dz }, crossProduct);
//
//        // ������Ϊ�㣬��ʾ���������Ѿ��غϣ�������ת
//        if (vtkMath::Norm(crossProduct) < 1e-6) {
//            // ���û����ת�����������غϣ�������������ת
//            std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//            rotationAngle = 0.0;
//        }
//        else {
//            rotationAngle = acos(dotProduct);  // ����н�
//            vtkMath::Normalize(crossProduct);  // ��һ����ת��
//            std::cout << "Rotation Axis: (" << crossProduct[0] << ", " << crossProduct[1] << ", " << crossProduct[2] << ")" << std::endl;
//        }
//
//        // ��ת�����嵽Ŀ�귽��
//        Greentransform->RotateWXYZ(rotationAngle * 180.0 / M_PI, crossProduct[0], crossProduct[1], crossProduct[2]);
//        // ����һ���任��������Ӧ�ñ任
//        vtkSmartPointer<vtkTransformPolyDataFilter> GreentransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        GreentransformFilter->SetInputConnection(GreenellipsoidSource->GetOutputPort());
//        GreentransformFilter->SetTransform(Greentransform);
//
//        // ����һ��ӳ����
//        vtkSmartPointer<vtkPolyDataMapper> Greenmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        Greenmapper->SetInputConnection(GreentransformFilter->GetOutputPort());
//
//        // ����һ��Actor
//        vtkSmartPointer<vtkActor> Greenactor = vtkSmartPointer<vtkActor>::New();
//        Greenactor->SetMapper(Greenmapper);
//        Greenactor->GetProperty()->SetOpacity(0.9);
//
//        // ������ɫ
//        Greenactor->GetProperty()->SetColor(0.0, 1.0, 0.0);  // ��ɫ  ��Ӧ����cuactor3
//        renderer->AddActor(Greenactor);//----------------------------------------------------------------------------------------------------------------------------------------------
//     // ����ӳ��������������ӵ���Ⱦ����
//        vtkSmartPointer<vtkPolyDataMapper> sphereMapperop = vtkSmartPointer<vtkPolyDataMapper>::New();
//        sphereMapperop->SetInputConnection(GreensphereDataop->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> sphereActorop = vtkSmartPointer<vtkActor>::New();
//        sphereActorop->SetMapper(sphereMapperop);
//        sphereActorop->GetProperty()->SetColor(139 / 255.0, 139 / 255.0, 0);//��ɫ
//        // ��������ӵ���Ⱦ����
//        //renderer->AddActor(sphereActorop);
//
//        // ���������ӻ��Ż�������������֮�������
//        vtkSmartPointer<vtkLineSource> lineSourceop = vtkSmartPointer<vtkLineSource>::New();
//        lineSourceop->SetPoint1(pathStart[0], pathStart[1], pathStart[2]);  // ���
//        lineSourceop->SetPoint2(tumorCenter4[0], tumorCenter4[1], tumorCenter4[2]);  // ��������
//         // ʹ�� TubeFilter �����߶εİ뾶
//        vtkSmartPointer<vtkTubeFilter> tubeFilterop = vtkSmartPointer<vtkTubeFilter>::New();
//        tubeFilterop->SetInputConnection(lineSourceop->GetOutputPort());
//        tubeFilterop->SetRadius(0.5);  // �����߶ΰ뾶
//        tubeFilterop->SetNumberOfSides(20);  // �����߶εľ���
//        tubeFilterop->Update();
//        // ����ӳ������������ӵ���Ⱦ����
//        vtkSmartPointer<vtkPolyDataMapper> lineMapperop = vtkSmartPointer<vtkPolyDataMapper>::New();
//        lineMapperop->SetInputConnection(tubeFilterop->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> lineActorop = vtkSmartPointer<vtkActor>::New();
//        lineActorop->SetMapper(lineMapperop);
//        lineActorop->GetProperty()->SetColor(139 / 255.0, 139 / 255.0, 0);  // �����ߵ���ɫΪ��ɫ
//
//        // ������ӵ���Ⱦ����
//       // renderer->AddActor(lineActorop);
////-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//    // ���� OBBTree ��������ײ��⣬�������Ѿ����� bone ����
//    vtkSmartPointer<vtkOBBTree> obbTree2 = vtkSmartPointer<vtkOBBTree>::New();
//    obbTree2->SetDataSet(bone);
//    obbTree2->BuildLocator();
//    // �����Ż��㷨����ȡ�Ż�������λ��
//    vtkSmartPointer<vtkPolyData> filteredPointsPolyData2 = ExtractPointsFromTrianglesAndProcess(filteredSkin3, tumorCenter3);
//    vtkSmartPointer<vtkPoints> points2 = filteredPointsPolyData2->GetPoints();
//    // ��ȡ�������
//    vtkIdType numPoints2 = points2->GetNumberOfPoints();
//    // ��������ȡ�Ż���ĵ�
//
//    vtkIdType randomIndex2 = rand() % numPoints2;  // rand() ���� [0, numPoints-1] ֮��������
//
//    double pathStart1[3];
//    points1->GetPoint(randomIndex, pathStart);  // ��ȡÿ���������
//    // �ڿ���̨������ѡ�еĵ�����
//    std::cout << "��ɫɫ������: ("
//        << pathStart[0] << ", "
//        << pathStart[1] << ", "
//        << pathStart[2] << ")"
//        << std::endl;
//
//    // Ϊÿ����㴴��һ������
//        // Ϊÿ����㴴��һ������
//        vtkSmartPointer<vtkSphereSource> yellowsphereDataop1 = vtkSmartPointer<vtkSphereSource>::New();
//        yellowsphereDataop1->SetCenter(pathStart1[0], pathStart1[1], pathStart1[2]);
//        yellowsphereDataop1->SetRadius(0.5);  // ��������뾶Ϊ 0.5
//         // ����һ����Դ������������ʼ��
//
//           // ����� tumorCenter4 �� pathStart �ķ�������
//        double dx1 = pathStart1[0] - tumorCenter3[0];
//        double dy1 = pathStart1[1] - tumorCenter3[1];
//        double dz1 = pathStart1[2] +60- tumorCenter3[2];
//
//        // ������������λ��
//        normalizeDirection(dx1, dy1, dz1);
//
//        // ����� tumorCenter4 ��ʼ������ pathStart 20 ��λ�������յ�
//        double extraLength2 = -80;
//        double totalLength2 = distance(tumorCenter3[0], tumorCenter3[1], tumorCenter3[2], pathStart1[0], pathStart1[1], pathStart1[2]) + extraLength2;
//
//        // ���������յ�
//        double rayEnd2[3] = {
//            tumorCenter3[0] + dx1 * totalLength2,
//            tumorCenter3[1] + dy1 * totalLength2,
//            tumorCenter3[2] + dz1 * totalLength2
//        };
//
//        vtkSmartPointer<vtkPoints> yellowpoints = vtkSmartPointer<vtkPoints>::New();
//        yellowpoints->InsertNextPoint(tumorCenter3); // ������� (tumorCenter4)
//        yellowpoints->InsertNextPoint(rayEnd2);       // �����յ� (rayEnd)
//
//        // ����һ����Դ����
//        vtkSmartPointer<vtkCellArray> yellowlines = vtkSmartPointer<vtkCellArray>::New();
//        yellowlines->InsertNextCell(2); // һ���������������
//        yellowlines->InsertCellPoint(0); // ���
//        yellowlines->InsertCellPoint(1); // �յ�
//
//        // ����һ��PolyData����
//        vtkSmartPointer<vtkPolyData> yellowlineData = vtkSmartPointer<vtkPolyData>::New();
//        yellowlineData->SetPoints(yellowpoints);
//        yellowlineData->SetLines(yellowlines);
//        // ʹ�� vtkTubeFilter ����Բ����
//        vtkSmartPointer<vtkTubeFilter>yellowtubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
//        yellowtubeFilter->SetInputData(yellowlineData);   // �����߶�������Ϊ����
//        yellowtubeFilter->SetRadius(1.0);            // ����Բ����İ뾶���൱�ڡ��߶εİ뾶����
//        yellowtubeFilter->SetNumberOfSides(50);      // ����Բ����ı�����Խ��Խ�⻬
//        // ����һ��ӳ��������ʾ��
//        vtkSmartPointer<vtkPolyDataMapper> yellowlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        yellowlineMapper->SetInputConnection(yellowtubeFilter->GetOutputPort());
//
//        // ����һ��Actor����Ⱦ��
//        vtkSmartPointer<vtkActor> yellowlineActor = vtkSmartPointer<vtkActor>::New();
//        yellowlineActor->SetMapper(yellowlineMapper);
//
//        // ������������ɫ
//        yellowlineActor->GetProperty()->SetColor(1, 1, 0);  // ��ɫ
//
//        renderer->AddActor(yellowlineActor);//-------------------------------------------------------------------------------------------------------------------------------------------------
//        // ���������峤�ᣨ20���Ͷ��ᣨ15��
//        double yellowsemiMajorAxis = 22.0;  // ����
//        double yellowsemiMinorAxis = 17;  // ����
//        //double yellowsemiMajorAxis1 = 25.0;  // ����
//        //double yellowsemiMinorAxis1 = 22;  // ����
//
//        // ����һ��������Ĳ���������
//        vtkSmartPointer<vtkParametricEllipsoid> yellowellipsoid = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        yellowellipsoid->SetXRadius(yellowsemiMinorAxis);  // ���ó���
//        yellowellipsoid->SetYRadius(yellowsemiMinorAxis);  // ���ö���
//        yellowellipsoid->SetZRadius(yellowsemiMajorAxis);  // ���ö���
//        // // ����һ��������Ĳ���������
//        //vtkSmartPointer<vtkParametricEllipsoid> yellowellipsoid1 = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        //yellowellipsoid1->SetXRadius(yellowsemiMinorAxis1);  // ���ó���
//        //yellowellipsoid1->SetYRadius(yellowsemiMinorAxis1);  // ���ö���
//        //yellowellipsoid1->SetZRadius(yellowsemiMajorAxis1);  // ���ö���
//
//        // ����һ������������Դ
//        vtkSmartPointer<vtkParametricFunctionSource> yellowellipsoidSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        yellowellipsoidSource->SetParametricFunction(yellowellipsoid);
//        //// ����һ������������Դ
//        //vtkSmartPointer<vtkParametricFunctionSource> yellowellipsoidSource1 = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        //yellowellipsoidSource1->SetParametricFunction(yellowellipsoid1);
//   
//        // ����������ı任
//        vtkSmartPointer<vtkTransform> yellowtransform = vtkSmartPointer<vtkTransform>::New();
//    /*    vtkSmartPointer<vtkTransform> yellowtransform1 = vtkSmartPointer<vtkTransform>::New();*/
//        // **ƽ�Ʋ���**��������������Ĵ� (0, 0, 0) �ƶ��� tumorCenter4
//        yellowtransform->Translate(tumorCenter3[0], tumorCenter3[1], tumorCenter3[2]);
//    /*    yellowtransform1->Translate(tumorCenter1[0] - 5, tumorCenter1[1] - 5, tumorCenter1[2] - 5);*/
//        // ������ת����Ŀ���ǽ�������ĳ��ᣨĬ���� z �ᣩ��ת��Ŀ�귽��dx, dy, dz��
//        double originalDirection1[3] = { 0.0, 0.0, 1.0 };  // Ĭ�ϳ��᷽��
//        double rotationAngle1 = 0.0;
//        double rotationAxis1[3] = { 0.0, 0.0, 0.0 };
//
//        // ������ת�ǶȺ���ת��
//        double dotProduct1 = vtkMath::Dot(originalDirection1, new double[3]{ dx1, dy1, dz1 });
//        double crossProduct1[3];
//        vtkMath::Cross(originalDirection1, new double[3]{ dx1, dy1, dz1 }, crossProduct1);
//
//        // ������Ϊ�㣬��ʾ���������Ѿ��غϣ�������ת
//        if (vtkMath::Norm(crossProduct1) < 1e-6) {
//            // ���û����ת�����������غϣ�������������ת
//            std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//            rotationAngle1 = 0.0;
//        }
//        else {
//            rotationAngle1 = acos(dotProduct1);  // ����н�
//            vtkMath::Normalize(crossProduct1);  // ��һ����ת��
//            std::cout << "Rotation Axis: (" << crossProduct1[0] << ", " << crossProduct1[1] << ", " << crossProduct1[2] << ")" << std::endl;
//        }
//
//        // ��ת�����嵽Ŀ�귽��
//        yellowtransform->RotateWXYZ(rotationAngle1 * 180.0 / M_PI, crossProduct1[0], crossProduct1[1], crossProduct1[2]);
//
//        // ������ת����Ŀ���ǽ�������ĳ��ᣨĬ���� z �ᣩ��ת��Ŀ�귽��dx, dy, dz��
//        double originalDirection11[3] = { 0.0, 0.0, 1.0 };  // Ĭ�ϳ��᷽��
//        double rotationAngle11 = 0.0;
//        double rotationAxis11[3] = { 0.0, 0.0, 0.0 };
//
//        //// ������ת�ǶȺ���ת��
//        //double dotProduct11 = vtkMath::Dot(originalDirection11, new double[3]{ dx1, dy1, dz1 });
//        //double crossProduct11[3];
//        //vtkMath::Cross(originalDirection11, new double[3]{ dx1, dy1, dz1 }, crossProduct11);
//
//        //// ������Ϊ�㣬��ʾ���������Ѿ��غϣ�������ת
//        //if (vtkMath::Norm(crossProduct11) < 1e-6) {
//        //    // ���û����ת�����������غϣ�������������ת
//        //    std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//        //    rotationAngle11 = 0.0;
//        //}
//        //else {
//        //    rotationAngle11 = acos(dotProduct11);  // ����н�
//        //    vtkMath::Normalize(crossProduct11);  // ��һ����ת��
//        //    std::cout << "Rotation Axis: (" << crossProduct11[0] << ", " << crossProduct11[1] << ", " << crossProduct11[2] << ")" << std::endl;
//        //}
//
//        //// ��ת�����嵽Ŀ�귽��
//        //yellowtransform1->RotateWXYZ(rotationAngle11 * 180.0 / M_PI, crossProduct11[0], crossProduct11[1], crossProduct11[2]);
//
//
//        // ����һ���任��������Ӧ�ñ任
//        vtkSmartPointer<vtkTransformPolyDataFilter> yellowtransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        yellowtransformFilter->SetInputConnection(yellowellipsoidSource->GetOutputPort());
//        yellowtransformFilter->SetTransform(yellowtransform);
//
//        // ����һ��ӳ����
//        vtkSmartPointer<vtkPolyDataMapper> yellowmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        yellowmapper->SetInputConnection(yellowtransformFilter->GetOutputPort());
//
//        // ����һ��Actor
//        vtkSmartPointer<vtkActor> yellowactor = vtkSmartPointer<vtkActor>::New();
//        yellowactor->SetMapper(yellowmapper);
//        yellowactor->GetProperty()->SetOpacity(0.9);
//
//        // ������ɫ
//        yellowactor->GetProperty()->SetColor(1.0, 1.0, 0.0);  //  cuactor1
//        renderer->AddActor(yellowactor);//--------------------------------------------------------------------------------------------------------------------------
//        //// ����һ���任��������Ӧ�ñ任
//        //vtkSmartPointer<vtkTransformPolyDataFilter> yellowtransformFilter1 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        //yellowtransformFilter1->SetInputConnection(yellowellipsoidSource1->GetOutputPort());
//        //yellowtransformFilter1->SetTransform(yellowtransform1);
//
//        //// ����һ��ӳ����
//        //vtkSmartPointer<vtkPolyDataMapper> yellowmapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        //yellowmapper1->SetInputConnection(yellowtransformFilter1->GetOutputPort());
//
//        //// ����һ��Actor
//        //vtkSmartPointer<vtkActor> yellowactor1 = vtkSmartPointer<vtkActor>::New();
//        //yellowactor1->SetMapper(yellowmapper1);
//        //yellowactor1->GetProperty()->SetOpacity(0.9);
//
//        //// ������ɫ
//        //yellowactor1->GetProperty()->SetColor(1.0, 1.0, 0.0);  //  cuactor1
//        //renderer->AddActor(yellowactor1);//-----------------------------------------------------------------------------------------------------------------------------
//        // // ����ӳ��������������ӵ���Ⱦ����
//        vtkSmartPointer<vtkPolyDataMapper> sphereMapperop1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        sphereMapperop1->SetInputConnection(yellowsphereDataop1->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> sphereActorop1 = vtkSmartPointer<vtkActor>::New();
//        sphereActorop1->SetMapper(sphereMapperop1);
//        sphereActorop1->GetProperty()->SetColor(0, 100 / 255.0, 0);//��ɫ
//        // ��������ӵ���Ⱦ����
//        //renderer->AddActor(sphereActorop1);
//
//        // ���������ӻ��Ż�������������֮�������
//        vtkSmartPointer<vtkLineSource> lineSourceop1 = vtkSmartPointer<vtkLineSource>::New();
//        lineSourceop1->SetPoint1(pathStart1[0], pathStart1[1], pathStart1[2]);  // ���
//        lineSourceop1->SetPoint2(tumorCenter3[0], tumorCenter3[1], tumorCenter3[2]);  // ��������
//         // ʹ�� TubeFilter �����߶εİ뾶
//        vtkSmartPointer<vtkTubeFilter> tubeFilterop1 = vtkSmartPointer<vtkTubeFilter>::New();
//        tubeFilterop1->SetInputConnection(lineSourceop1->GetOutputPort());
//        tubeFilterop1->SetRadius(0.5);  // �����߶ΰ뾶
//        tubeFilterop1->SetNumberOfSides(20);  // �����߶εľ���
//        tubeFilterop1->Update();
//        // ����ӳ������������ӵ���Ⱦ����
//        vtkSmartPointer<vtkPolyDataMapper> lineMapperop1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        lineMapperop1->SetInputConnection(tubeFilterop1->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> lineActorop1 = vtkSmartPointer<vtkActor>::New();
//        lineActorop1->SetMapper(lineMapperop1);
//        lineActorop1->GetProperty()->SetColor(0, 100 / 255.0, 0);  // �����ߵ���ɫΪ��ɫ
//
//        // ������ӵ���Ⱦ����
//       //renderer->AddActor(lineActorop1);
//    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//    //// ���� OBBTree ��������ײ��⣬�������Ѿ����� bone ����
//    //vtkSmartPointer<vtkOBBTree> obbTree3 = vtkSmartPointer<vtkOBBTree>::New();
//    //obbTree3->SetDataSet(bone);
//    //obbTree3->BuildLocator();
//    //// �����Ż��㷨����ȡ�Ż�������λ��
//        vtkSmartPointer<vtkPolyData> filteredPointsPolyData3 = ExtractPointsFromTrianglesAndProcess(filteredSkin1, tumorCenter1);
//        vtkSmartPointer<vtkPoints> points3 = filteredPointsPolyData3->GetPoints();
//        // ��ȡ�������
//        vtkIdType numPoints3 = points3->GetNumberOfPoints();
//        // ��������ȡ�Ż���ĵ�
//
//        vtkIdType randomIndex3 = rand() % numPoints3;  // rand() ���� [0, numPoints-1] ֮��������
//
//        double pathStart2[3];
//        points3->GetPoint(randomIndex3, pathStart2);  // ��ȡÿ���������
//        // �ڿ���̨������ѡ�еĵ�����
//        std::cout << "��ɫ������: ("
//            << pathStart[0] << ", "
//            << pathStart[1] << ", "
//            << pathStart[2] << ")"
//            << std::endl;
//
//        // Ϊÿ����㴴��һ������
//        vtkSmartPointer<vtkSphereSource> redsphereDataop2 = vtkSmartPointer<vtkSphereSource>::New();
//        redsphereDataop2->SetCenter(pathStart2[0], pathStart2[1], pathStart2[2]);
//        redsphereDataop2->SetRadius(0.5);  // ��������뾶Ϊ 0.5
//
//        // ����� tumorCenter4 �� pathStart �ķ�������
//        double dx3 = pathStart2[0]+10  - tumorCenter1[0];
//        double dy3 = pathStart2[1] - tumorCenter1[1];
//        double dz3 = pathStart2[2]+20 - tumorCenter1[2];
//
//        // ������������λ��
//        normalizeDirection(dx3, dy3, dz3);
//
//        // ����� tumorCenter4 ��ʼ������ pathStart 20 ��λ�������յ�
//        double extraLength3 = 50.0;
//        double totalLength3 = distance(tumorCenter1[0], tumorCenter1[1], tumorCenter1[2], pathStart2[0], pathStart2[1], pathStart2[2]) + extraLength3;
//
//        // ���������յ�
//        double rayEnd3[3] = {
//            tumorCenter1[0] + dx3 * totalLength3,
//            tumorCenter1[1] + dy3 * totalLength3,
//            tumorCenter1[2] + dz3 * totalLength3
//        };
//
//        vtkSmartPointer<vtkPoints> redpoints = vtkSmartPointer<vtkPoints>::New();
//        redpoints->InsertNextPoint(tumorCenter1); // ������� (tumorCenter4)
//        redpoints->InsertNextPoint(rayEnd3);       // �����յ� (rayEnd)
//
//        // ����һ����Դ����
//        vtkSmartPointer<vtkCellArray> redlines = vtkSmartPointer<vtkCellArray>::New();
//        redlines->InsertNextCell(2); // һ���������������
//        redlines->InsertCellPoint(0); // ���
//        redlines->InsertCellPoint(1); // �յ�
//
//        // ����һ��PolyData����
//        vtkSmartPointer<vtkPolyData> redlineData = vtkSmartPointer<vtkPolyData>::New();
//        redlineData->SetPoints(redpoints);
//        redlineData->SetLines(redlines);
//        // ʹ�� vtkTubeFilter ����Բ����
//        vtkSmartPointer<vtkTubeFilter>redtubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilter->SetInputData(redlineData);   // �����߶�������Ϊ����
//        redtubeFilter->SetRadius(1);            // ����Բ����İ뾶���൱�ڡ��߶εİ뾶����
//        redtubeFilter->SetNumberOfSides(50);      // ����Բ����ı�����Խ��Խ�⻬
//
//        // ����ӳ��������������ӵ���Ⱦ����
//        vtkSmartPointer<vtkPolyDataMapper> redsphereMapperop2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redsphereMapperop2->SetInputConnection(redsphereDataop2->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redsphereActorop2 = vtkSmartPointer<vtkActor>::New();
//        redsphereActorop2->SetMapper(redsphereMapperop2);
//        redsphereActorop2->GetProperty()->SetColor(139 / 255.0, 0, 0);//��ɫ
//        // ��������ӵ���Ⱦ����
//        //renderer->AddActor(sphereActorop2);
//
//        // ���������ӻ��Ż�������������֮�������
//        vtkSmartPointer<vtkLineSource> redlineSourceop2 = vtkSmartPointer<vtkLineSource>::New();
//        redlineSourceop2->SetPoint1(pathStart2[0], pathStart2[1], pathStart2[2]);  // ���
//        redlineSourceop2->SetPoint2(tumorCenter1[0], tumorCenter1[1], tumorCenter1[2]);  // ��������
//         // ʹ�� TubeFilter �����߶εİ뾶
//        vtkSmartPointer<vtkTubeFilter> redtubeFilterop2 = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilterop2->SetInputConnection(redlineSourceop2->GetOutputPort());
//        redtubeFilterop2->SetRadius(1);  // �����߶ΰ뾶
//        redtubeFilterop2->SetNumberOfSides(50);  // �����߶εľ���
//        redtubeFilterop2->Update();
//        // ����ӳ������������ӵ���Ⱦ����
//        vtkSmartPointer<vtkPolyDataMapper> redlineMapperop2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redlineMapperop2->SetInputConnection(redtubeFilter->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redlineActorop2 = vtkSmartPointer<vtkActor>::New();
//        redlineActorop2->SetMapper(redlineMapperop2);
//        redlineActorop2->GetProperty()->SetColor(0, 0, 139 / 255.0);  // �����ߵ���ɫΪ��ɫ
//
//        // ������ӵ���Ⱦ����
//        renderer->AddActor(redlineActorop2);//------------------------------------------------------------------------------------------------------
//
//        // ���������峤�ᣨ20���Ͷ��ᣨ15��
//        double redsemiMajorAxis = 22.0;  // ����
//        double redsemiMinorAxis =17.0;  // ����
//
//        // ����һ��������Ĳ���������
//        vtkSmartPointer<vtkParametricEllipsoid> redellipsoid = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        redellipsoid->SetXRadius(redsemiMinorAxis);  // ���ó���
//        redellipsoid->SetYRadius(redsemiMinorAxis);  // ���ö���
//        redellipsoid->SetZRadius(redsemiMajorAxis);  // ���ö���
//
//        // ����һ������������Դ
//        vtkSmartPointer<vtkParametricFunctionSource> redellipsoidSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        redellipsoidSource->SetParametricFunction(redellipsoid);
//       
//        // ����������ı任
//        vtkSmartPointer<vtkTransform> redtransform = vtkSmartPointer<vtkTransform>::New();
//        // **ƽ�Ʋ���**��������������Ĵ� (0, 0, 0) �ƶ��� tumorCenter4
//        redtransform->Translate(tumorCenter1[0], tumorCenter1[1], tumorCenter1[2]);
//        // ������ת����Ŀ���ǽ�������ĳ��ᣨĬ���� z �ᣩ��ת��Ŀ�귽��dx, dy, dz��
//        double originalDirection2[3] = { 0.0, 0.0, 1.0 };  // Ĭ�ϳ��᷽��
//        double rotationAngle2 = 0.0;
//        double rotationAxis2[3] = { 0.0, 0.0, 0.0 };
//
//        // ������ת�ǶȺ���ת��
//        double dotProduct2 = vtkMath::Dot(originalDirection2, new double[3]{ dx3, dy3, dz3 });
//        double crossProduct2[3];
//        vtkMath::Cross(originalDirection2, new double[3]{ dx3, dy3, dz3 }, crossProduct2);
//
//        // ������Ϊ�㣬��ʾ���������Ѿ��غϣ�������ת
//        if (vtkMath::Norm(crossProduct2) < 1e-6) {
//            // ���û����ת�����������غϣ�������������ת
//            std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//            rotationAngle2 = 0.0;
//        }
//        else {
//            rotationAngle2 = acos(dotProduct2);  // ����н�
//            vtkMath::Normalize(crossProduct2);  // ��һ����ת��
//            std::cout << "Rotation Axis: (" << crossProduct2[0] << ", " << crossProduct2[1] << ", " << crossProduct2[2] << ")" << std::endl;
//        }
//
//        // ��ת�����嵽Ŀ�귽��
//        redtransform->RotateWXYZ(rotationAngle2 * 180.0 / M_PI, crossProduct2[0], crossProduct2[1], crossProduct2[2]);
//
//
//
//        // ����һ���任��������Ӧ�ñ任
//        vtkSmartPointer<vtkTransformPolyDataFilter> redtransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        redtransformFilter->SetInputConnection(redellipsoidSource->GetOutputPort());
//        redtransformFilter->SetTransform(redtransform);
//
//        // ����һ��ӳ����
//        vtkSmartPointer<vtkPolyDataMapper> redmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redmapper->SetInputConnection(redtransformFilter->GetOutputPort());
//
//        // ����һ��Actor
//        vtkSmartPointer<vtkActor> redactor = vtkSmartPointer<vtkActor>::New();
//        redactor->SetMapper(redmapper);
//        redactor->GetProperty()->SetOpacity(0.9);
//
//        // ������ɫ
//        redactor->GetProperty()->SetColor(0, 0, 139 / 255.0);  //  cuactor4
//        renderer->AddActor(redactor);
//    //}
//    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//     //// ���� OBBTree ��������ײ��⣬�������Ѿ����� bone ����
//    //vtkSmartPointer<vtkOBBTree> obbTree4 = vtkSmartPointer<vtkOBBTree>::New();
//    //obbTree3->SetDataSet(bone);
//    //obbTree3->BuildLocator();
//    //// �����Ż��㷨����ȡ�Ż�������λ��
//        vtkSmartPointer<vtkPolyData> filteredPointsPolyData4 = ExtractPointsFromTrianglesAndProcess(filteredSkin5, tumorCenter5);
//        vtkSmartPointer<vtkPoints> points4 = filteredPointsPolyData4->GetPoints();
//        // ��ȡ�������
//        vtkIdType numPoints4 = points4->GetNumberOfPoints();
//        // ��������ȡ�Ż���ĵ�
//
//        vtkIdType randomIndex4 = rand() % numPoints4;  // rand() ���� [0, numPoints-1] ֮��������
//
//        double pathStart4[3];
//        points3->GetPoint(randomIndex4, pathStart4);  // ��ȡÿ���������
//        // �ڿ���̨������ѡ�еĵ�����
//        std::cout << "��ɫ������: ("
//            << pathStart[0] << ", "
//            << pathStart[1] << ", "
//            << pathStart[2] << ")"
//            << std::endl;
//
//        // Ϊÿ����㴴��һ������
//        vtkSmartPointer<vtkSphereSource> redsphereDataop21 = vtkSmartPointer<vtkSphereSource>::New();
//        redsphereDataop21->SetCenter(pathStart4[0], pathStart4[1], pathStart4[2]);
//        redsphereDataop21->SetRadius(0.5);  // ��������뾶Ϊ 0.5
//
//        // ����� tumorCenter4 �� pathStart �ķ�������
//        double dx4 = pathStart4[0] -50 - tumorCenter5[0];
//        double dy4 = pathStart4[1]  -tumorCenter5[1];
//        double dz4 = pathStart4[2] +80 - tumorCenter5[2];
//
//        // ������������λ��
//        normalizeDirection(dx4, dy4, dz4);
//
//        // ����� tumorCenter4 ��ʼ������ pathStart 20 ��λ�������յ�
//        double pkx = 50;
//        double totalLength4 = distance(tumorCenter5[0], tumorCenter5[1], tumorCenter5[2], pathStart4[0], pathStart4[1], pathStart4[2])-pkx ;
//
//        // ���������յ�
//        double rayEnd4[3] = {
//            tumorCenter5[0] + dx4 * totalLength4,
//            tumorCenter5[1]+30 + dy4 * totalLength4,
//            tumorCenter5[2] + dz4 * totalLength4
//        };
//
//        vtkSmartPointer<vtkPoints> redpoints21 = vtkSmartPointer<vtkPoints>::New();
//        redpoints21->InsertNextPoint(tumorCenter5); // ������� (tumorCenter4)
//        redpoints21->InsertNextPoint(rayEnd4);       // �����յ� (rayEnd)
//
//        // ����һ����Դ����
//        vtkSmartPointer<vtkCellArray> redlines21 = vtkSmartPointer<vtkCellArray>::New();
//        redlines21->InsertNextCell(2); // һ���������������
//        redlines21->InsertCellPoint(0); // ���
//        redlines21->InsertCellPoint(1); // �յ�
//
//        // ����һ��PolyData����
//        vtkSmartPointer<vtkPolyData> redlineData21 = vtkSmartPointer<vtkPolyData>::New();
//        redlineData21->SetPoints(redpoints21);
//        redlineData21->SetLines(redlines21);
//        // ʹ�� vtkTubeFilter ����Բ����
//        vtkSmartPointer<vtkTubeFilter>redtubeFilter21 = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilter21->SetInputData(redlineData21);   // �����߶�������Ϊ����
//        redtubeFilter21->SetRadius(1);            // ����Բ����İ뾶���൱�ڡ��߶εİ뾶����
//        redtubeFilter21->SetNumberOfSides(50);      // ����Բ����ı�����Խ��Խ�⻬
//
//        // ����ӳ��������������ӵ���Ⱦ����
//        vtkSmartPointer<vtkPolyDataMapper> redsphereMapperop21 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redsphereMapperop21->SetInputConnection(redsphereDataop21->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redsphereActorop21 = vtkSmartPointer<vtkActor>::New();
//        redsphereActorop21->SetMapper(redsphereMapperop21);
//        redsphereActorop21->GetProperty()->SetColor(139 / 255.0, 0, 0);//��ɫ
//        // ��������ӵ���Ⱦ����
//        //renderer->AddActor(sphereActorop2);
//
//        // ���������ӻ��Ż�������������֮�������
//        vtkSmartPointer<vtkLineSource> redlineSourceop21 = vtkSmartPointer<vtkLineSource>::New();
//        redlineSourceop21->SetPoint1(pathStart4[0], pathStart4[1], pathStart4[2]);  // ���
//        redlineSourceop21->SetPoint2(tumorCenter5[0], tumorCenter5[1], tumorCenter5[2]);  // ��������
//         // ʹ�� TubeFilter �����߶εİ뾶
//        vtkSmartPointer<vtkTubeFilter> redtubeFilterop21 = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilterop21->SetInputConnection(redlineSourceop21->GetOutputPort());
//        redtubeFilterop21->SetRadius(1);  // �����߶ΰ뾶
//        redtubeFilterop21->SetNumberOfSides(50);  // �����߶εľ���
//        redtubeFilterop21->Update();
//        // ����ӳ������������ӵ���Ⱦ����
//        vtkSmartPointer<vtkPolyDataMapper> redlineMapperop21 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redlineMapperop21->SetInputConnection(redtubeFilter21->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redlineActorop21 = vtkSmartPointer<vtkActor>::New();
//        redlineActorop21->SetMapper(redlineMapperop21);
//        redlineActorop21->GetProperty()->SetColor(1, 100 / 255.0, 0);  // �����ߵ���ɫΪ��ɫ
//
//        // ������ӵ���Ⱦ����
//        renderer->AddActor(redlineActorop21);//------------------------------------------------------------------------------------------------------
//
//        // ���������峤�ᣨ20���Ͷ��ᣨ15��
//        double redsemiMajorAxis21 = 23.0;  // ����
//        double redsemiMinorAxis21 = 18.0;  // ����
//
//        // ����һ��������Ĳ���������
//        vtkSmartPointer<vtkParametricEllipsoid> redellipsoid21 = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        redellipsoid21->SetXRadius(redsemiMinorAxis21);  // ���ó���
//        redellipsoid21->SetYRadius(redsemiMinorAxis21);  // ���ö���
//        redellipsoid21->SetZRadius(redsemiMajorAxis21);  // ���ö���
//
//        // ����һ������������Դ
//        vtkSmartPointer<vtkParametricFunctionSource> redellipsoidSource21 = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        redellipsoidSource21->SetParametricFunction(redellipsoid21);
//
//        // ����������ı任
//        vtkSmartPointer<vtkTransform> redtransform21 = vtkSmartPointer<vtkTransform>::New();
//        // **ƽ�Ʋ���**��������������Ĵ� (0, 0, 0) �ƶ��� tumorCenter4
//        redtransform21->Translate(tumorCenter5[0], tumorCenter5[1], tumorCenter5[2]);
//        // ������ת����Ŀ���ǽ�������ĳ��ᣨĬ���� z �ᣩ��ת��Ŀ�귽��dx, dy, dz��
//        double originalDirection21[3] = { 0.0, 0.0, 1.0 };  // Ĭ�ϳ��᷽��
//        double rotationAngle21 = 0.0;
//        double rotationAxis21[3] = { 0.0, 0.0, 0.0 };
//
//        // ������ת�ǶȺ���ת��
//        double dotProduct21 = vtkMath::Dot(originalDirection21, new double[3]{ dx4, dy4, dz4 });
//        double crossProduct21[3];
//        vtkMath::Cross(originalDirection21, new double[3]{ dx4, dy4, dz4 }, crossProduct21);
//
//        // ������Ϊ�㣬��ʾ���������Ѿ��غϣ�������ת
//        if (vtkMath::Norm(crossProduct21) < 1e-6) {
//            // ���û����ת�����������غϣ�������������ת
//            std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//            rotationAngle21 = 0.0;
//        }
//        else {
//            rotationAngle21 = acos(dotProduct21);  // ����н�
//            vtkMath::Normalize(crossProduct21);  // ��һ����ת��
//            std::cout << "Rotation Axis: (" << crossProduct21[0] << ", " << crossProduct21[1] << ", " << crossProduct21[2] << ")" << std::endl;
//        }
//
//        // ��ת�����嵽Ŀ�귽��
//        redtransform21->RotateWXYZ(rotationAngle21* 180.0 / M_PI, crossProduct21[0], crossProduct21[1], crossProduct21[2]);
//
//
//
//        // ����һ���任��������Ӧ�ñ任
//        vtkSmartPointer<vtkTransformPolyDataFilter> redtransformFilter21 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        redtransformFilter21->SetInputConnection(redellipsoidSource21->GetOutputPort());
//        redtransformFilter21->SetTransform(redtransform21);
//
//        // ����һ��ӳ����
//        vtkSmartPointer<vtkPolyDataMapper> redmapper21 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redmapper21->SetInputConnection(redtransformFilter21->GetOutputPort());
//
//        // ����һ��Actor
//        vtkSmartPointer<vtkActor> redactor21 = vtkSmartPointer<vtkActor>::New();
//        redactor21->SetMapper(redmapper21);
//        redactor21->GetProperty()->SetOpacity(0.9);
//
//        // ������ɫ
//        redactor21->GetProperty()->SetColor(1, 100 / 255.0, 0);  //  cuactor4
//        renderer->AddActor(redactor21);
//        //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//        // //// ���� OBBTree ��������ײ��⣬�������Ѿ����� bone ����
//    //vtkSmartPointer<vtkOBBTree> obbTree5 = vtkSmartPointer<vtkOBBTree>::New();
//    //obbTree3->SetDataSet(bone);
//    //obbTree3->BuildLocator();
//    //// �����Ż��㷨����ȡ�Ż�������λ��
//        vtkSmartPointer<vtkPolyData> filteredPointsPolyData5 = ExtractPointsFromTrianglesAndProcess(filteredSkin6, tumorCenter6);
//        vtkSmartPointer<vtkPoints> points5 = filteredPointsPolyData5->GetPoints();
//        // ��ȡ�������
//        vtkIdType numPoints5 = points5->GetNumberOfPoints();
//        // ��������ȡ�Ż���ĵ�
//
//        vtkIdType randomIndex5 = rand() % numPoints5;  // rand() ���� [0, numPoints-1] ֮��������
//
//        double pathStart5[3];
//        points3->GetPoint(randomIndex5, pathStart5);  // ��ȡÿ���������
//        // �ڿ���̨������ѡ�еĵ�����
//        std::cout << "��ɫ������: ("
//            << pathStart[0] << ", "
//            << pathStart[1] << ", "
//            << pathStart[2] << ")"
//            << std::endl;
//
//        // Ϊÿ����㴴��һ������
//        vtkSmartPointer<vtkSphereSource> redsphereDataop211 = vtkSmartPointer<vtkSphereSource>::New();
//        redsphereDataop211->SetCenter(pathStart5[0], pathStart5[1], pathStart5[2]);
//        redsphereDataop211->SetRadius(0.5);  // ��������뾶Ϊ 0.5
//
//        // ����� tumorCenter4 �� pathStart �ķ�������
//        double dx5= pathStart5[0]  - tumorCenter6[0];
//        double dy5 = pathStart5[1]  - tumorCenter6[1];
//        double dz5 = pathStart5[2] - tumorCenter6[2];
//
//        // ������������λ��
//        normalizeDirection(dx5, dy5, dz5);
//
//        // ����� tumorCenter4 ��ʼ������ pathStart 20 ��λ�������յ�
//        double o = 50;
//        double totalLength5 = distance(tumorCenter6[0], tumorCenter6[1], tumorCenter6[2], pathStart5[0], pathStart5[1], pathStart5[2])+o;
//
//        // ���������յ�
//        double rayEnd5[3] = {
//            tumorCenter6[0] + dx5 * totalLength5,
//            tumorCenter6[1] + dy5 * totalLength5,
//            tumorCenter6[2] + dz5 * totalLength5
//        };
//
//        vtkSmartPointer<vtkPoints> redpoints211 = vtkSmartPointer<vtkPoints>::New();
//        redpoints211->InsertNextPoint(tumorCenter6); // ������� (tumorCenter4)
//        redpoints211->InsertNextPoint(rayEnd5);       // �����յ� (rayEnd)
//
//        // ����һ����Դ����
//        vtkSmartPointer<vtkCellArray> redlines211 = vtkSmartPointer<vtkCellArray>::New();
//        redlines211->InsertNextCell(2); // һ���������������
//        redlines211->InsertCellPoint(0); // ���
//        redlines211->InsertCellPoint(1); // �յ�
//
//        // ����һ��PolyData����
//        vtkSmartPointer<vtkPolyData> redlineData211 = vtkSmartPointer<vtkPolyData>::New();
//        redlineData211->SetPoints(redpoints211);
//        redlineData211->SetLines(redlines211);
//        // ʹ�� vtkTubeFilter ����Բ����
//        vtkSmartPointer<vtkTubeFilter>redtubeFilter211 = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilter211->SetInputData(redlineData211);   // �����߶�������Ϊ����
//        redtubeFilter211->SetRadius(1);            // ����Բ����İ뾶���൱�ڡ��߶εİ뾶����
//        redtubeFilter211->SetNumberOfSides(50);      // ����Բ����ı�����Խ��Խ�⻬
//
//        // ����ӳ��������������ӵ���Ⱦ����
//        vtkSmartPointer<vtkPolyDataMapper> redsphereMapperop211 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redsphereMapperop211->SetInputConnection(redsphereDataop211->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redsphereActorop211 = vtkSmartPointer<vtkActor>::New();
//        redsphereActorop211->SetMapper(redsphereMapperop211);
//        redsphereActorop211->GetProperty()->SetColor(139 / 255.0, 0, 0);//��ɫ
//        // ��������ӵ���Ⱦ����
//        //renderer->AddActor(sphereActorop2);
//
//        // ���������ӻ��Ż�������������֮�������
//        vtkSmartPointer<vtkLineSource> redlineSourceop211 = vtkSmartPointer<vtkLineSource>::New();
//        redlineSourceop211->SetPoint1(pathStart5[0], pathStart5[1], pathStart5[2]);  // ���
//        redlineSourceop211->SetPoint2(tumorCenter6[0], tumorCenter6[1], tumorCenter6[2]);  // ��������
//         // ʹ�� TubeFilter �����߶εİ뾶
//        vtkSmartPointer<vtkTubeFilter> redtubeFilterop211 = vtkSmartPointer<vtkTubeFilter>::New();
//        redtubeFilterop211->SetInputConnection(redlineSourceop211->GetOutputPort());
//        redtubeFilterop211->SetRadius(1);  // �����߶ΰ뾶
//        redtubeFilterop211->SetNumberOfSides(50);  // �����߶εľ���
//        redtubeFilterop211->Update();
//        // ����ӳ������������ӵ���Ⱦ����
//        vtkSmartPointer<vtkPolyDataMapper> redlineMapperop211 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redlineMapperop211->SetInputConnection(redtubeFilter211->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> redlineActorop211 = vtkSmartPointer<vtkActor>::New();
//        redlineActorop211->SetMapper(redlineMapperop211);
//        redlineActorop211->GetProperty()->SetColor(0, 100 / 255.0, 1);  // �����ߵ���ɫΪ��ɫ
//
//        // ������ӵ���Ⱦ����
//        renderer->AddActor(redlineActorop211);//------------------------------------------------------------------------------------------------------
//
//        // ���������峤�ᣨ20���Ͷ��ᣨ15��
//        double redsemiMajorAxis211 = 22.0;  // ����
//        double redsemiMinorAxis211 = 17.0;  // ����
//
//        // ����һ��������Ĳ���������
//        vtkSmartPointer<vtkParametricEllipsoid> redellipsoid211 = vtkSmartPointer<vtkParametricEllipsoid>::New();
//        redellipsoid211->SetXRadius(redsemiMinorAxis211);  // ���ó���
//        redellipsoid211->SetYRadius(redsemiMinorAxis211);  // ���ö���
//        redellipsoid211->SetZRadius(redsemiMajorAxis211);  // ���ö���
//
//        // ����һ������������Դ
//        vtkSmartPointer<vtkParametricFunctionSource> redellipsoidSource211 = vtkSmartPointer<vtkParametricFunctionSource>::New();
//        redellipsoidSource211->SetParametricFunction(redellipsoid211);
//
//        // ����������ı任
//        vtkSmartPointer<vtkTransform> redtransform211 = vtkSmartPointer<vtkTransform>::New();
//        // **ƽ�Ʋ���**��������������Ĵ� (0, 0, 0) �ƶ��� tumorCenter4
//        redtransform211->Translate(tumorCenter6[0], tumorCenter6[1], tumorCenter6[2]);
//        // ������ת����Ŀ���ǽ�������ĳ��ᣨĬ���� z �ᣩ��ת��Ŀ�귽��dx, dy, dz��
//        double originalDirection211[3] = { 0.0, 0.0, 1.0 };  // Ĭ�ϳ��᷽��
//        double rotationAngle211 = 0.0;
//        double rotationAxis211[3] = { 0.0, 0.0, 0.0 };
//
//        // ������ת�ǶȺ���ת��
//        double dotProduct211 = vtkMath::Dot(originalDirection211, new double[3]{ dx5, dy5, dz5 });
//        double crossProduct211[3];
//        vtkMath::Cross(originalDirection211, new double[3]{ dx5, dy5, dz5 }, crossProduct211);
//
//        // ������Ϊ�㣬��ʾ���������Ѿ��غϣ�������ת
//        if (vtkMath::Norm(crossProduct211) < 1e-6) {
//            // ���û����ת�����������غϣ�������������ת
//            std::cout << "No rotation needed. Vectors are already aligned." << std::endl;
//            rotationAngle211 = 0.0;
//        }
//        else {
//            rotationAngle211 = acos(dotProduct211);  // ����н�
//            vtkMath::Normalize(crossProduct211);  // ��һ����ת��
//            std::cout << "Rotation Axis: (" << crossProduct211[0] << ", " << crossProduct211[1] << ", " << crossProduct211[2] << ")" << std::endl;
//        }
//
//        // ��ת�����嵽Ŀ�귽��
//        redtransform211->RotateWXYZ(rotationAngle211 * 180.0 / M_PI, crossProduct211[0], crossProduct211[1], crossProduct211[2]);
//
//
//
//        // ����һ���任��������Ӧ�ñ任
//        vtkSmartPointer<vtkTransformPolyDataFilter> redtransformFilter211 = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
//        redtransformFilter211->SetInputConnection(redellipsoidSource211->GetOutputPort());
//        redtransformFilter211->SetTransform(redtransform211);
//
//        // ����һ��ӳ����
//        vtkSmartPointer<vtkPolyDataMapper> redmapper211 = vtkSmartPointer<vtkPolyDataMapper>::New();
//        redmapper211->SetInputConnection(redtransformFilter211->GetOutputPort());
//
//        // ����һ��Actor
//        vtkSmartPointer<vtkActor> redactor211 = vtkSmartPointer<vtkActor>::New();
//        redactor211->SetMapper(redmapper211);
//        redactor211->GetProperty()->SetOpacity(0.9);
//
//        // ������ɫ
//        redactor211->GetProperty()->SetColor(0, 100 / 255.0, 1);  //  cuactor4
//        renderer->AddActor(redactor211);
////-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//    vtkSmartPointer<vtkActor> filteredSkinActor = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor->SetMapper(filteredSkinMapper);
//    filteredSkinActor->GetProperty()->SetColor(1.0, 1.0, 0.0);
//    filteredSkinActor->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor1 = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor1->SetMapper(filteredSkinMapper1);
//    filteredSkinActor1->GetProperty()->SetColor(0, 0, 139 / 255.0);//��ɫ
//    filteredSkinActor1->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor2 = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor2->SetMapper(filteredSkinMapper2);
//    filteredSkinActor2->GetProperty()->SetColor(139 / 255.0, 0, 0);//��ɫ
//    filteredSkinActor2->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor3 = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor3->SetMapper(filteredSkinMapper3);
//    filteredSkinActor3->GetProperty()->SetColor(139 / 255.0, 139 / 255.0, 0);//huangɫ
//    filteredSkinActor3->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor4 = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor4->SetMapper(filteredSkinMapper4);
//    filteredSkinActor4->GetProperty()->SetColor(0, 100 / 255.0, 0);//lvɫ
//    filteredSkinActor4->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor5= vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor5->SetMapper(filteredSkinMapper5);
//    filteredSkinActor5->GetProperty()->SetColor(1, 100 / 255.0, 0);//lvɫ
//    filteredSkinActor5->GetProperty()->SetOpacity(0.9);
//    vtkSmartPointer<vtkActor> filteredSkinActor6 = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor6->SetMapper(filteredSkinMapper6);
//    filteredSkinActor6->GetProperty()->SetColor(0, 100 / 255.0, 1);//lvɫ
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
//    renderer->SetBackground(1, 1, 1); // ���ñ�����ɫ��double����
//    renderWindow->Render();
//    renderWindowInteractor->Start();
//
//    return 0;
//}
