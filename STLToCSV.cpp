//#include <vtkSmartPointer.h>
//#include <vtkSTLReader.h>
//#include <vtkPoints.h>
//#include <vtkCellArray.h>
//#include <vtkPolyData.h>
//#include <fstream>
//#include <iostream>
//
//int main() {
//    // ��ȡ STL �ļ�
//    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
//    reader->SetFileName("D://Data Disk//LW//Data//3Dircadb1.1//STL//livertumor17.stl");
//    reader->Update();
//
//    // ��ȡ������
//    vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();
//    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
//
//    // ���� CSV �ļ�
//    std::ofstream csvFile("D://Data Disk//LW//Data//3Dircadb1.1//CSV//livertumor17.csv");
//    if (!csvFile.is_open()) {
//        std::cerr << "Failed to open CSV file." << std::endl;
//        return -1;
//    }
//
//    // д�� CSV �ļ�ͷ
//    csvFile << "X,Y,Z" << std::endl;
//
//    // ����ÿ���㲢д�� CSV �ļ�
//    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
//        double point[3];
//        points->GetPoint(i, point);
//        csvFile << point[0] << "," << point[1] << "," << point[2] << std::endl;
//    }
//
//    // �ر��ļ�
//    csvFile.close();
//    std::cout << "CSV file saved successfully!" << std::endl;
//
//    return 0;
//}
