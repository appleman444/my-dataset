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
//
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
//   
//    double tumorCenter1[] = { 101.90, 99.48, 105.09 };
//    double tumorCenter2[] = { 105.40, 75.65, 104.77 };
//    double tumorCenter3[] = { 87.99, 86.50, 119.05 };
//    double tumorCenter4[] = { 75.90, 97.84, 98.02 };
//    // Create a points object to store the tumor centers
//    vtkSmartPointer<vtkPoints>  centerpoints1 = vtkSmartPointer<vtkPoints>::New();
//    centerpoints1->InsertNextPoint(tumorCenter1);
//    vtkSmartPointer<vtkPoints>  centerpoints2 = vtkSmartPointer<vtkPoints>::New();
//    centerpoints2->InsertNextPoint(tumorCenter2);
//    vtkSmartPointer<vtkPoints>  centerpoints3 = vtkSmartPointer<vtkPoints>::New();
//    centerpoints3->InsertNextPoint(tumorCenter3);
//    vtkSmartPointer<vtkPoints>  centerpoints4 = vtkSmartPointer<vtkPoints>::New();
//    centerpoints4->InsertNextPoint(tumorCenter4);
//
//    // Create a sphere source
//    vtkSmartPointer<vtkSphereSource> centersphereSource1 = vtkSmartPointer<vtkSphereSource>::New();
//    centersphereSource1->SetRadius(1.0);  // Set the radius of the sphere to 1
//    centersphereSource1->SetPhiResolution(20);  // Resolution of the sphere (smoothness)
//    centersphereSource1->SetThetaResolution(20);
//    vtkSmartPointer<vtkSphereSource> centersphereSource2 = vtkSmartPointer<vtkSphereSource>::New();
//    centersphereSource2->SetRadius(1.0);  // Set the radius of the sphere to 1
//    centersphereSource2->SetPhiResolution(20);  // Resolution of the sphere (smoothness)
//    centersphereSource2->SetThetaResolution(20);
//    vtkSmartPointer<vtkSphereSource> centersphereSource3 = vtkSmartPointer<vtkSphereSource>::New();
//    centersphereSource3->SetRadius(1.0);  // Set the radius of the sphere to 1
//    centersphereSource3->SetPhiResolution(20);  // Resolution of the sphere (smoothness)
//    centersphereSource3->SetThetaResolution(20);
//    vtkSmartPointer<vtkSphereSource> centersphereSource4 = vtkSmartPointer<vtkSphereSource>::New();
//    centersphereSource4->SetRadius(1.0);  // Set the radius of the sphere to 1
//    centersphereSource4->SetPhiResolution(20);  // Resolution of the sphere (smoothness)
//    centersphereSource4->SetThetaResolution(20);
//    vtkSmartPointer<vtkPolyData> centerpolyData1 = vtkSmartPointer<vtkPolyData>::New();
//    centerpolyData1->SetPoints(centerpoints1); // ���õ�����
//    vtkSmartPointer<vtkPolyData> centerpolyData2 = vtkSmartPointer<vtkPolyData>::New();
//    centerpolyData2->SetPoints(centerpoints2); // ���õ�����
//    vtkSmartPointer<vtkPolyData> centerpolyData3 = vtkSmartPointer<vtkPolyData>::New();
//    centerpolyData3->SetPoints(centerpoints3); // ���õ�����
//    vtkSmartPointer<vtkPolyData> centerpolyData4 = vtkSmartPointer<vtkPolyData>::New();
//
//    centerpolyData4->SetPoints(centerpoints4); // ���õ�����
//    // 4. ʹ�� vtkGlyph3D ��ÿ����ת��Ϊ����
//    vtkSmartPointer<vtkGlyph3D> cglyph3D1 = vtkSmartPointer<vtkGlyph3D>::New();
//    cglyph3D1->SetInputData(centerpolyData1);  // ���������
//    cglyph3D1->SetSourceConnection(centersphereSource1->GetOutputPort());  // ��������Դ
//    cglyph3D1->Update();  // �������ɵĽ��
//    vtkSmartPointer<vtkGlyph3D> cglyph3D2 = vtkSmartPointer<vtkGlyph3D>::New();
//    cglyph3D2->SetInputData(centerpolyData2);  // ���������
//    cglyph3D2->SetSourceConnection(centersphereSource2->GetOutputPort());  // ��������Դ
//    cglyph3D2->Update();  // �������ɵĽ��
//    vtkSmartPointer<vtkGlyph3D> cglyph3D3 = vtkSmartPointer<vtkGlyph3D>::New();
//    cglyph3D3->SetInputData(centerpolyData3);  // ���������
//    cglyph3D3->SetSourceConnection(centersphereSource3->GetOutputPort());  // ��������Դ
//    cglyph3D3->Update();  // �������ɵĽ��
//    vtkSmartPointer<vtkGlyph3D> cglyph3D4 = vtkSmartPointer<vtkGlyph3D>::New();
//    cglyph3D4->SetInputData(centerpolyData4);  // ���������
//    cglyph3D4->SetSourceConnection(centersphereSource4->GetOutputPort());  // ��������Դ
//    cglyph3D4->Update();  // �������ɵĽ��
//
//    // Create a mapper and actor for each tumor center
//    vtkSmartPointer<vtkPolyDataMapper>  centersphereMapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    centersphereMapper1->SetInputConnection(cglyph3D1->GetOutputPort());
//    vtkSmartPointer<vtkPolyDataMapper>  centersphereMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    centersphereMapper2->SetInputConnection(cglyph3D2->GetOutputPort());
//    vtkSmartPointer<vtkPolyDataMapper>  centersphereMapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    centersphereMapper3->SetInputConnection(cglyph3D3->GetOutputPort());
//    vtkSmartPointer<vtkPolyDataMapper>  centersphereMapper4 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    centersphereMapper4->SetInputConnection(cglyph3D4->GetOutputPort());
//
//    vtkSmartPointer<vtkActor>  centersphereActor1 = vtkSmartPointer<vtkActor>::New();
//    centersphereActor1->SetMapper(centersphereMapper1);
//    centersphereActor1->GetProperty()->SetColor(1, 0, 0);  // ������ɫΪ��ɫ
//    vtkSmartPointer<vtkActor>  centersphereActor2 = vtkSmartPointer<vtkActor>::New();
//    centersphereActor2->SetMapper(centersphereMapper2);
//    centersphereActor2->GetProperty()->SetColor(0, 0, 1);  // ������ɫΪ��ɫ
//    vtkSmartPointer<vtkActor>  centersphereActor3 = vtkSmartPointer<vtkActor>::New();
//    centersphereActor3->SetMapper(centersphereMapper3);
//    centersphereActor3->GetProperty()->SetColor(1, 1, 0);  // ������ɫΪ��ɫ
//    vtkSmartPointer<vtkActor>  centersphereActor4 = vtkSmartPointer<vtkActor>::New();
//    centersphereActor4->SetMapper(centersphereMapper4);
//    centersphereActor4->GetProperty()->SetColor(0, 1, 0);  // ������ɫΪ��ɫ
//
//   
//  
//
//        renderer->AddActor(centersphereActor1);
//        renderer->AddActor(centersphereActor2);
//        renderer->AddActor(centersphereActor3);
//        renderer->AddActor(centersphereActor4);
//  
//   
//    std::string fileName = "D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\CSV\\PKX\\tumor_cluster_1.csv";
//    std::string fileName1 = "D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\CSV\\PKX\\tumor_cluster_2.csv";
//    std::string fileName2 = "D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\CSV\\PKX\\tumor_cluster_3.csv";
//    std::string fileName3 = "D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\CSV\\PKX\\tumor_cluster_4.csv";
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
//
//    // ��ȡ�����ݽ��洢�� vtkTable ��
//    vtkSmartPointer<vtkTable> table = reader->GetOutput();
//    vtkSmartPointer<vtkTable> table1 = reader1->GetOutput();
//    vtkSmartPointer<vtkTable> table2 = reader2->GetOutput();
//    vtkSmartPointer<vtkTable> table3 = reader3->GetOutput();
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
//    // 3. ʹ�� vtkSphereSource ��������
//    vtkSmartPointer<vtkSphereSource> cusphereSource1 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource1->SetRadius(0.5);  // ��������İ뾶
//    cusphereSource1->SetPhiResolution(10);  // ��������ֱ���
//    cusphereSource1->SetThetaResolution(10);  // ���ú���ֱ���
//    cusphereSource1->Update();  // ��������Դ
//    vtkSmartPointer<vtkSphereSource> cusphereSource2 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource2->SetRadius(0.5);  // ��������İ뾶
//    cusphereSource2->SetPhiResolution(10);  // ��������ֱ���
//    cusphereSource2->SetThetaResolution(10);  // ���ú���ֱ���
//    cusphereSource2->Update();  // ��������Դ
//    vtkSmartPointer<vtkSphereSource> cusphereSource3 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource3->SetRadius(0.5);  // ��������İ뾶
//    cusphereSource3->SetPhiResolution(10);  // ��������ֱ���
//    cusphereSource3->SetThetaResolution(10);  // ���ú���ֱ���
//    cusphereSource3->Update();  // ��������Դ
//    vtkSmartPointer<vtkSphereSource> cusphereSource4 = vtkSmartPointer<vtkSphereSource>::New();
//    cusphereSource4->SetRadius(0.5);  // ��������İ뾶
//    cusphereSource4->SetPhiResolution(10);  // ��������ֱ���
//    cusphereSource4->SetThetaResolution(10);  // ���ú���ֱ���
//    cusphereSource4->Update();  // ��������Դ
//    vtkSmartPointer<vtkPolyData> cupolyData1 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData1->SetPoints(cupoints1); // ���õ�����
//    vtkSmartPointer<vtkPolyData> cupolyData2 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData2->SetPoints(cupoints2); // ���õ�����
//    vtkSmartPointer<vtkPolyData> cupolyData3 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData3->SetPoints(cupoints3); // ���õ�����
//    vtkSmartPointer<vtkPolyData> cupolyData4 = vtkSmartPointer<vtkPolyData>::New();
//    cupolyData4->SetPoints(cupoints4); // ���õ�����
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
//
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
//
//  
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
//    cuactor1->GetProperty()->SetPointSize(2);
//    // ����һ�� Actor ����ʾ��
//    vtkSmartPointer<vtkActor> cuactor2 = vtkSmartPointer<vtkActor>::New();
//    cuactor2->SetMapper(cumapper2);
//    cuactor2->GetProperty()->SetColor(1, 0, 0);//��ɫ
//    cuactor2->GetProperty()->SetPointSize(2);
//    // ����һ�� Actor ����ʾ��
//    vtkSmartPointer<vtkActor> cuactor3 = vtkSmartPointer<vtkActor>::New();
//    cuactor3->SetMapper(cumapper3);
//    cuactor3->GetProperty()->SetColor(0.0, 1.0, 0.0);//lv se
//    cuactor3->GetProperty()->SetPointSize(2);
//    // ����һ�� Actor ����ʾ��
//    vtkSmartPointer<vtkActor> cuactor4 = vtkSmartPointer<vtkActor>::New();
//    cuactor4->SetMapper(cumapper4);
//    cuactor4->GetProperty()->SetColor(0.0, 0.0, 1.0);//red
//    cuactor4->GetProperty()->SetPointSize(2);
//
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
//    tumorActor->GetProperty()->SetColor(1, 1, 1);
//    tumorActor->GetProperty()->SetOpacity(0.5);
//
//    vtkSmartPointer<vtkActor> ArteryActor = vtkSmartPointer<vtkActor>::New();
//    ArteryActor->SetMapper(ArteryMapper);
//    ArteryActor->GetProperty()->SetColor(1, 0.5, 0);
//
//     renderer->AddActor(actor);
//    renderer->AddActor(actor2);
//     renderer->AddActor(actor3);
//      renderer->AddActor(actor4);
//      //renderer->AddActor(filteredSkinActor);
//     // renderer->AddActor(filteredSkinActor1);
//      //renderer->AddActor(filteredSkinActor2);
//      //renderer->AddActor(filteredSkinActor3);
//     // renderer->AddActor(filteredSkinActor4);
//      renderer->AddActor(cuactor1);
//      renderer->AddActor(cuactor2);
//      renderer->AddActor(cuactor3);
//      renderer->AddActor(cuactor4);
//    renderer->AddActor(SkinActor);
//    renderer->AddActor(venoussystemActor);
//    renderer->AddActor(BoneActor);
//    //renderer->AddActor(tumorActor);
//    renderer->AddActor(liverActor);
//    renderer->AddActor(ArteryActor);
//    renderer->SetBackground(1, 1, 1); // ���ñ�����ɫ��double����
//    renderWindow->Render();
//    renderWindowInteractor->Start();
//
//    return 0;
//}
