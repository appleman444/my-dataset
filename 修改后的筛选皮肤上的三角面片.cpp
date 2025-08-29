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
//
//using namespace std;
//const double PI = 3.141592653589793;
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
//// 筛选满足要求的三角面片
//vtkSmartPointer<vtkPolyData> FilterValidSkinTriangles(vtkSmartPointer<vtkPolyData> skin, vtkSmartPointer<vtkPolyData> bone, double* tumorCenter) {
//    vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
//    triangleFilter->SetInputData(skin);
//    triangleFilter->Update();
//
//    vtkSmartPointer<vtkOBBTree> obbTree = vtkSmartPointer<vtkOBBTree>::New();
//    obbTree->SetDataSet(bone);
//    obbTree->BuildLocator();
//
//    vtkSmartPointer<vtkPolyData> filteredTriangles = vtkSmartPointer<vtkPolyData>::New();
//    vtkSmartPointer<vtkCellArray> newCells = vtkSmartPointer<vtkCellArray>::New();
//    vtkSmartPointer<vtkPoints> points = triangleFilter->GetOutput()->GetPoints();
//    vtkSmartPointer<vtkCellArray> polys = triangleFilter->GetOutput()->GetPolys();
//
//    vtkIdType npts;
//    const vtkIdType* ptIds;
//    polys->InitTraversal();
//    while (polys->GetNextCell(npts, ptIds)) {
//        if (npts == 3) {
//            bool valid = true;
//            for (int i = 0; i < 3; i++) {
//                double point[3];
//                points->GetPoint(ptIds[i], point);
//                vtkNew<vtkPoints> intersectPoints;
//                if (obbTree->IntersectWithLine(tumorCenter, point, intersectPoints, nullptr)) {
//                    valid = false;
//                    break;
//                }
//            }
//            if (valid) {
//                newCells->InsertNextCell(3, ptIds);
//            }
//        }
//    }
//    filteredTriangles->SetPoints(points);
//    filteredTriangles->SetPolys(newCells);
//    return filteredTriangles;
//}
//
//int main() {
//    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//    renderWindow->AddRenderer(renderer);
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    renderWindowInteractor->SetRenderWindow(renderWindow);
//
//    vtkSmartPointer<vtkPolyData> bone = Reader_VTK("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\", "bone.vtk");
//    vtkSmartPointer<vtkPolyData> tumor = Reader_VTK("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\", "livertumor01.vtk");
//    vtkSmartPointer<vtkPolyData> skin = Reader_VTK("D:\\dataSet\\3Dircadb1\\3Dircadb1.1\\MESHES_VTK\\", "skin.vtk");
//
//    double* tumorCenter = tumor->GetCenter();
//    vtkSmartPointer<vtkPolyData> filteredSkin = FilterValidSkinTriangles(skin, bone, tumorCenter);
//
//    vtkSmartPointer<vtkPolyDataMapper> filteredSkinMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    filteredSkinMapper->SetInputData(filteredSkin);
//    vtkSmartPointer<vtkPolyDataMapper> BoneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    BoneMapper->SetInputData(bone);
//    vtkSmartPointer<vtkPolyDataMapper> tumorMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//   tumorMapper->SetInputData(tumor);
//    vtkSmartPointer<vtkPolyDataMapper> SkinMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//  SkinMapper->SetInputData(skin);
//    vtkSmartPointer<vtkActor> filteredSkinActor = vtkSmartPointer<vtkActor>::New();
//    filteredSkinActor->SetMapper(filteredSkinMapper);
//    filteredSkinActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
//    filteredSkinActor->GetProperty()->SetOpacity(0.4);
//    vtkSmartPointer<vtkActor> BoneActor = vtkSmartPointer<vtkActor>::New();
//    BoneActor->SetMapper(BoneMapper);
//    BoneActor->GetProperty()->SetColor(1,1,1);
//    vtkSmartPointer<vtkActor> SkinActor = vtkSmartPointer<vtkActor>::New();
//    SkinActor->SetMapper(SkinMapper);
//    SkinActor->GetProperty()->SetColor(1, 1, 1);
//    SkinActor->GetProperty()->SetOpacity(0.1);
//    vtkSmartPointer<vtkActor> tumorActor = vtkSmartPointer<vtkActor>::New();
//    tumorActor->SetMapper(tumorMapper);
//    tumorActor->GetProperty()->SetColor(1, 0, 0);
//   
//    renderer->AddActor(filteredSkinActor);
//    renderer->AddActor(SkinActor);
//    renderer->AddActor(BoneActor);
//    renderer->AddActor(tumorActor);
//    renderer->SetBackground(1, 1, 1);//设置背景颜色，double类型
//    renderWindow->Render();
//    renderWindowInteractor->Start();
//    return 0;
//}
