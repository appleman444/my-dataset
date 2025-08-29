//#include <vtkSTLReader.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <vtkPoints.h>
//#include <vtkBoundingBox.h>
//#include <iostream>
//
//int main() {
//    // ��ȡSTL�ļ���ȷ��·����ȷ
//    vtkSmartPointer<vtkSTLReader> stlReader = vtkSmartPointer<vtkSTLReader>::New();
//    stlReader->SetFileName("D:/Data Disk/LW/Data/3Dircadb1.1/STL/livertumor01.stl"); // ʹ����б��
//    stlReader->Update();
//
//    // ��ȡģ�͵� polydata
//    vtkSmartPointer<vtkPolyData> tumorData = stlReader->GetOutput();
//
//    // ��ȡ������
//    vtkSmartPointer<vtkPoints> points = tumorData->GetPoints();
//
//    // ����һ����Χ�ж���
//    vtkBoundingBox boundingBox; // ʹ�ù��캯����������
//
//    // �����е���ӵ���Χ����
//    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
//        double point[3];
//        points->GetPoint(i, point);
//        boundingBox.AddPoint(point);
//    }
//
//    // ��ȡ��Χ�е�������С����
//    double minBounds[3], maxBounds[3];
//    boundingBox.GetMinPoint(minBounds);
//    boundingBox.GetMaxPoint(maxBounds);
//
//    // ���������ĳߴ�
//    double a = maxBounds[0] - minBounds[0]; // X����ĳߴ�
//    double b = maxBounds[1] - minBounds[1]; // Y����ĳߴ�
//    double c = maxBounds[2] - minBounds[2]; // Z����ĳߴ�
//
//    // ��������ĳߴ�
//    std::cout << "Tumor size (a �� b �� c): "
//        << a << " �� " << b << " �� " << c << std::endl;
//
//    return 0;
//}
