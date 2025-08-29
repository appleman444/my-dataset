//#include <vtkAutoInit.h>
//#include <vtkPolyDataReader.h>
//#include <vtkSTLWriter.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <iostream>
//
//VTK_MODULE_INIT(vtkRenderingOpenGL2);
//VTK_MODULE_INIT(vtkRenderingFreeType);
//VTK_MODULE_INIT(vtkInteractionStyle);
//
//int main() {
//    // 读取 VTK 文件
//    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
//    reader->SetFileName("D://Data Disk//LW//Data//3Dircadb1.1//MESHES_VTK//livertumor1.vtk");
//    reader->Update();
//
//    vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();
//
//    // 创建 STL 写入器
//    vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
//    writer->SetFileName("D://Data Disk//LW//Data//3Dircadb1.1//STL//livertumor17.stl");
//
//    // 设置输入数据并写入 STL 文件
//    writer->SetInputData(polyData);
//    writer->Write();
//
//    std::cout << "STL file saved successfully!" << std::endl;
//
//    return 0;
//}
//
//
