//#include <vtkSTLReader.h>
//#include <vtkMassProperties.h>
//#include <vtkPolyData.h>
//#include <vtkSmartPointer.h>
//#include <vtkTriangleFilter.h>
//#include <iostream>
//
//int main() {
//    // 1. ��ȡSTL�ļ�
//    vtkSmartPointer<vtkSTLReader> stlReader = vtkSmartPointer<vtkSTLReader>::New();
//    stlReader->SetFileName("D:\\Data Disk\\LW\\Data\\3Dircadb1.1\\STL\\livertumor01.stl");
//    stlReader->Update();
//
//    // 2. ���STL�ļ�û����������Ƭ��ʽ������ת��Ϊ��������Ƭ
//    vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
//    triangleFilter->SetInputData(stlReader->GetOutput());
//    triangleFilter->Update();
//
//    // 3. ʹ�� vtkMassProperties �������
//    vtkSmartPointer<vtkMassProperties> massProperties = vtkSmartPointer<vtkMassProperties>::New();
//    massProperties->SetInputData(triangleFilter->GetOutput());
//
//    // 4. ��ȡ���
//    double volume = massProperties->GetVolume();
//    std::cout << "Tumor volume: " << volume << std::endl;
//
//    return 0;
//}
