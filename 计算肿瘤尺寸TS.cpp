//#include <vtkSTLReader.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <vtkPoints.h>
//#include <vtkBoundingBox.h>
//#include <iostream>
//
//int main() {
//    // 读取STL文件，确保路径正确
//    vtkSmartPointer<vtkSTLReader> stlReader = vtkSmartPointer<vtkSTLReader>::New();
//    stlReader->SetFileName("D:/Data Disk/LW/Data/3Dircadb1.1/STL/livertumor01.stl"); // 使用正斜杠
//    stlReader->Update();
//
//    // 获取模型的 polydata
//    vtkSmartPointer<vtkPolyData> tumorData = stlReader->GetOutput();
//
//    // 获取点数据
//    vtkSmartPointer<vtkPoints> points = tumorData->GetPoints();
//
//    // 创建一个包围盒对象
//    vtkBoundingBox boundingBox; // 使用构造函数创建对象
//
//    // 将所有点添加到包围盒中
//    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
//        double point[3];
//        points->GetPoint(i, point);
//        boundingBox.AddPoint(point);
//    }
//
//    // 获取包围盒的最大和最小坐标
//    double minBounds[3], maxBounds[3];
//    boundingBox.GetMinPoint(minBounds);
//    boundingBox.GetMaxPoint(maxBounds);
//
//    // 计算肿瘤的尺寸
//    double a = maxBounds[0] - minBounds[0]; // X方向的尺寸
//    double b = maxBounds[1] - minBounds[1]; // Y方向的尺寸
//    double c = maxBounds[2] - minBounds[2]; // Z方向的尺寸
//
//    // 输出肿瘤的尺寸
//    std::cout << "Tumor size (a × b × c): "
//        << a << " × " << b << " × " << c << std::endl;
//
//    return 0;
//}
