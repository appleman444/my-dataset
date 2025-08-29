//#include "vtkAutoInit.h"
//VTK_MODULE_INIT(vtkRenderingOpenGL2);
//VTK_MODULE_INIT(vtkInteractionStyle);
//
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <vtkPoints.h>
//#include <vtkPolyDataReader.h>
//#include <vtkSelectEnclosedPoints.h>
//#include <vtkIntArray.h>
//#include <vtkPointData.h>
//#include <vtkSphereSource.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkNamedColors.h>
//#include <vtkProperty.h>
//// ����㼯�� CSV �ļ�
//void SavePointsToCSV(vtkSmartPointer<vtkPoints> points, const std::string& fileName) {
//    // ���ļ������ļ��������򴴽�
//    std::ofstream outFile(fileName);
//
//    if (!outFile.is_open()) {
//        std::cerr << "�޷����ļ�����д��: " << fileName << std::endl;
//        return;
//    }
//
//    // д�� CSV ��ͷ
//    outFile << "X,Y,Z" << std::endl;
//
//    // �����㣬д��ÿ���������
//    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
//        double* point = points->GetPoint(i);  // ��ȡ��ǰ�������
//        outFile << point[0] << "," << point[1] << "," << point[2] << std::endl;
//    }
//
//    // �ر��ļ�
//    outFile.close();
//    std::cout << "���ѱ��浽�ļ�: " << fileName << std::endl;
//}
//
//
//
//int main() {
//    // ��ȡ��������
//    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
//    reader->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.17//MESHES_VTK//livertumor2.vtk");  // �����������������ļ�
//    reader->Update();
//    vtkSmartPointer<vtkPolyData> poly_tumor = reader->GetOutput();
//
//    // ��ȡ�����İ�Χ��
//    double tumorBounds[6];
//    poly_tumor->GetBounds(tumorBounds);
//    double tumorMinX = tumorBounds[0], tumorMaxX = tumorBounds[1];
//    double tumorMinY = tumorBounds[2], tumorMaxY = tumorBounds[3];
//    double tumorMinZ = tumorBounds[4], tumorMaxZ = tumorBounds[5];
//
//    // ���ɾ��ȷֲ��� sample_pointsPolydata (�������߽緶Χ��)
//    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//    vtkSmartPointer<vtkPolyData> sample_pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
//    int numPoints = 20000;  // �������� 1000 ����
//
//    // �������ɷ�Χ�޶��������İ�Χ�з�Χ��
//    for (int i = 0; i < numPoints; i++) {
//        double x = tumorMinX + (tumorMaxX - tumorMinX) * (rand() / (RAND_MAX + 1.0));  // ������� x ����
//        double y = tumorMinY + (tumorMaxY - tumorMinY) * (rand() / (RAND_MAX + 1.0));  // ������� y ����
//        double z = tumorMinZ + (tumorMaxZ - tumorMinZ) * (rand() / (RAND_MAX + 1.0));  // ������� z ����
//        points->InsertNextPoint(x, y, z);
//    }
//
//    // ���㼯��ӵ� sample_pointsPolydata ��
//    sample_pointsPolydata->SetPoints(points);
//
//    // ʹ�� vtkSelectEnclosedPoints �жϵ��Ƿ��������ڲ�
//    vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
//    enclosedPoints->SetInputData(sample_pointsPolydata);  // ����Ҫ�жϵĵ�����
//    enclosedPoints->SetSurfaceData(poly_tumor);           // ���������ı�������
//    enclosedPoints->Update();
//
//    // ��ȡÿ�����Ƿ��������ڲ�
//    vtkSmartPointer<vtkDataArray> insideArray = enclosedPoints->GetOutput()->GetPointData()->GetArray("SelectedPoints");
//
//    // ��������ڲ��������
//    std::vector<vtkIdType> insidePoints;  // �洢�����ڲ��ĵ������
//    std::vector<vtkIdType> outPoints;  // �洢�����ڲ��ĵ������
//    // ���� vtkPoints ���������������ڲ��ĵ�
//    vtkSmartPointer<vtkPoints> insidePointsData = vtkSmartPointer<vtkPoints>::New();
//    for (vtkIdType i = 0; i < sample_pointsPolydata->GetNumberOfPoints(); i++) {
//        // ʹ�� GetComponent ����ȡ�Ƿ��������ڵ���Ϣ
//        bool insideTumor = insideArray->GetComponent(i, 0) == 1;  // �жϵ��Ƿ���������
//        double inpoint[3];
//        sample_pointsPolydata->GetPoint(i, inpoint);  // ��ȡ�������
//
//        if (insideTumor) {
//            insidePoints.push_back(i);  // ������������ڣ��������������б�
//            insidePointsData->InsertNextPoint(inpoint);  // ���õ���뵽 vtkPoints ������
//        }
//        else {
//            outPoints.push_back(i);  // ������������ڣ��������������б�
//        }
//    }
//
//    std::cout << "Points inside the tumor: " << insidePoints.size() << std::endl;
//    vtkSmartPointer<vtkPolyData> insidePolyData = vtkSmartPointer<vtkPolyData>::New();
//    insidePolyData->SetPoints(insidePointsData);
//
//    // �� main �������޸��ļ�·��
//    std::string fileName = "D:/Data Disk/LW/Data/3Dircadb1.1/CSV/tumorinpoints/livertumor2-17.csv";  // �����ļ�·��
//    SavePointsToCSV(insidePolyData->GetPoints(), fileName);  // ���������ڲ��ĵ㵽ָ��·��
//
//    // ���ӻ������ڲ����С��
//    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//    renderWindow->AddRenderer(renderer);
//
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    renderWindowInteractor->SetRenderWindow(renderWindow);
//
//    vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
//
//    // Ϊ�����ڲ���ÿ���㴴��С��
//    for (size_t i = 0; i < insidePoints.size(); i++) {
//        double point[3];
//        sample_pointsPolydata->GetPoint(insidePoints[i], point);
//
//        // ����С��
//        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
//        sphereSource->SetCenter(point);
//        sphereSource->SetRadius(0.2);  // ����С��뾶Ϊ0.5
//
//        vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
//        sphereActor->SetMapper(sphereMapper);
//
//        // ����С����ɫ
//        sphereActor->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());
//
//        // ��С����ӵ���Ⱦ��
//        renderer->AddActor(sphereActor);
//    }
//   // Ϊ�����ⲿ��ÿ���㴴��С��
//    //for (size_t i = 0; i < outPoints.size(); i++) {
//    //    double point1[3];
//    //    sample_pointsPolydata->GetPoint(outPoints[i], point1);
//
//    //    // ����С��
//    //    vtkSmartPointer<vtkSphereSource> sphereSource1 = vtkSmartPointer<vtkSphereSource>::New();
//    //    sphereSource1->SetCenter(point1);
//    //    sphereSource1->SetRadius(0.15);  // ����С��뾶Ϊ0.5
//
//    //    vtkSmartPointer<vtkPolyDataMapper> sphereMapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    //    sphereMapper1->SetInputConnection(sphereSource1->GetOutputPort());
//
//    //    vtkSmartPointer<vtkActor> sphereActor1= vtkSmartPointer<vtkActor>::New();
//    //    sphereActor1->SetMapper(sphereMapper1);
//
//    //    // ����С����ɫ
//    //    sphereActor1->GetProperty()->SetColor(colors->GetColor3d("Green").GetData());
//
//    //    // ��С����ӵ���Ⱦ��
//    //    renderer->AddActor(sphereActor1);
//    //}
//
//    // ���ӻ�������������ɫ��͸����0.4��
//    vtkSmartPointer<vtkPolyDataMapper> tumorMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    tumorMapper->SetInputData(poly_tumor);
//
//    vtkSmartPointer<vtkActor> tumorActor = vtkSmartPointer<vtkActor>::New();
//    tumorActor->SetMapper(tumorMapper);
//
//    // ������������ɫ��͸����
//    tumorActor->GetProperty()->SetColor(1, 1, 1);  // ��ɫ
//    tumorActor->GetProperty()->SetOpacity(0.4);  // ͸����0.4
//
//    // ��������ӵ���Ⱦ��
//   renderer->AddActor(tumorActor);
//
//    // ���ñ�����ɫ
//    renderer->SetBackground(colors->GetColor3d("White").GetData());
//
//    // ������Ⱦ
//    renderWindow->Render();
//    renderWindowInteractor->Start();
//
//    return 0;
//}
