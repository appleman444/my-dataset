//#include <vtkSmartPointer.h>
//#include <vtkSTLReader.h>
//#include <vtkPoints.h>
//#include <vtkCellArray.h>
//#include <vtkPolyData.h>
//#include <fstream>
//#include <iostream>
//
//int main() {
//    // 读取 STL 文件
//    vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
//    reader->SetFileName("D://Data Disk//LW//Data//3Dircadb1.1//STL//livertumor17.stl");
//    reader->Update();
//
//    // 获取点数据
//    vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();
//    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
//
//    // 创建 CSV 文件
//    std::ofstream csvFile("D://Data Disk//LW//Data//3Dircadb1.1//CSV//livertumor17.csv");
//    if (!csvFile.is_open()) {
//        std::cerr << "Failed to open CSV file." << std::endl;
//        return -1;
//    }
//
//    // 写入 CSV 文件头
//    csvFile << "X,Y,Z" << std::endl;
//
//    // 遍历每个点并写入 CSV 文件
//    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
//        double point[3];
//        points->GetPoint(i, point);
//        csvFile << point[0] << "," << point[1] << "," << point[2] << std::endl;
//    }
//
//    // 关闭文件
//    csvFile.close();
//    std::cout << "CSV file saved successfully!" << std::endl;
//
//    return 0;
//}
