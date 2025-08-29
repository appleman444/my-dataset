//#include <vtkSTLReader.h>
//#include <vtkMassProperties.h>
//#include <vtkPolyData.h>
//#include <vtkSmartPointer.h>
//#include <vtkTriangleFilter.h>
//#include <iostream>
//
//int main() {
//    // 1. 读取STL文件
//    vtkSmartPointer<vtkSTLReader> stlReader = vtkSmartPointer<vtkSTLReader>::New();
//    stlReader->SetFileName("D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\STL\\livertumor01.stl");
//    stlReader->Update();
//
//    // 2. 如果STL文件没有三角形面片形式，可以转换为三角形面片
//    vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
//    triangleFilter->SetInputData(stlReader->GetOutput());
//    triangleFilter->Update();
//
//    // 3. 使用 vtkMassProperties 计算体积
//    vtkSmartPointer<vtkMassProperties> massProperties = vtkSmartPointer<vtkMassProperties>::New();
//    massProperties->SetInputData(triangleFilter->GetOutput());
//
//    // 4. 获取体积
//    double volume = massProperties->GetVolume();
//    std::cout << "Tumor volume: " << volume << std::endl;
//
//    return 0;
//}
