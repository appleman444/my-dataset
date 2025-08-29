//#include "vtkAutoInit.h"
//VTK_MODULE_INIT(vtkRenderingOpenGL2);
//VTK_MODULE_INIT(vtkInteractionStyle);
//
//#include <vtkPolyDataReader.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkNamedColors.h>
//#include <vtkCellArray.h>
//#include <iostream>
//#include <vtkPolyData.h>
//#include <vtkProperty.h>
//#include <vtkPoints.h>
//#include <vtkLine.h>
//int main()
//{
//    vtkSmartPointer<vtkRenderer> aRenderer = vtkSmartPointer<vtkRenderer>::New();
//    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
//    renWin->AddRenderer(aRenderer);
//
//    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    iren->SetRenderWindow(renWin);
//    //// 基于三维空间中两点画直线
//    double p0[3] = { 102.564,192.021,125.602 };
//    double p1[3] = { 100.0, 100.0, 100.0 };
//
//    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
//    pts->InsertNextPoint(p0);
//    pts->InsertNextPoint(p1);
//
//    vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
//    line0->GetPointIds()->SetId(0, 0);
//    line0->GetPointIds()->SetId(1, 1);
//
//    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
//    lines->InsertNextCell(line0);
//
//    vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
//    linesPolyData->SetPoints(pts);
//    linesPolyData->SetLines(lines);
//
//    // 创建路径的Mapper和Actor
//    vtkSmartPointer<vtkPolyDataMapper> lineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    lineMapper->SetInputData(linesPolyData);
//
//    vtkSmartPointer<vtkActor> lineActor = vtkSmartPointer<vtkActor>::New();
//    lineActor->SetMapper(lineMapper);
//    lineActor->GetProperty()->SetColor(1.0, 0.0, 0.0);  // 红色显示路径
//    lineActor->GetProperty()->SetLineWidth(3);
//
//    // 读取 liver 数据
//    vtkSmartPointer<vtkPolyDataReader> vtkReaderLiver = vtkSmartPointer<vtkPolyDataReader>::New();
//    vtkReaderLiver->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.1//MESHES_VTK//MESHES_VTK//livertumor02.vtk");
//    vtkReaderLiver->Update();
//    vtkSmartPointer<vtkPolyDataReader> vtkReaderLiver2 = vtkSmartPointer<vtkPolyDataReader>::New();
//    vtkReaderLiver2->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.1//MESHES_VTK//MESHES_VTK//livertumor03.vtk");
//    vtkReaderLiver2->Update();
//    vtkSmartPointer<vtkPolyDataReader> vtkReaderLiver3 = vtkSmartPointer<vtkPolyDataReader>::New();
//    vtkReaderLiver3->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.1//MESHES_VTK//MESHES_VTK//livertumor04.vtk");
//    vtkReaderLiver3->Update();
//    vtkSmartPointer<vtkPolyDataReader> vtkReaderLiver4 = vtkSmartPointer<vtkPolyDataReader>::New();
//    vtkReaderLiver4->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.1//MESHES_VTK//MESHES_VTK//livertumor05.vtk");
//    vtkReaderLiver4->Update();
//
//    // 读取 artery 数据
//    vtkSmartPointer<vtkPolyDataReader> vtkReaderArtery = vtkSmartPointer<vtkPolyDataReader>::New();
//    vtkReaderArtery->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.1//MESHES_VTK//MESHES_VTK//artery.vtk");
//    vtkReaderArtery->Update();
//
//    vtkSmartPointer<vtkPolyDataReader> vtkReaderSkin = vtkSmartPointer<vtkPolyDataReader>::New();
//    vtkReaderSkin->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.1//MESHES_VTK//MESHES_VTK//skin.vtk");
//    vtkReaderSkin->Update();
//    // 读取 bone 数据
//    vtkSmartPointer<vtkPolyDataReader> vtkReaderBone = vtkSmartPointer<vtkPolyDataReader>::New();
//    vtkReaderBone->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.1//MESHES_VTK//MESHES_VTK//bone.vtk");
//    vtkReaderBone->Update();
//
//    // 创建数据映射器
//    vtkSmartPointer<vtkPolyDataMapper> liverMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    liverMapper->SetInputConnection(vtkReaderLiver->GetOutputPort());
//   // liverMapper->ScalarVisibilityOff(); // 不显示颜色
//    vtkSmartPointer<vtkPolyDataMapper> liverMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    liverMapper2->SetInputConnection(vtkReaderLiver2->GetOutputPort());
//    //liverMapper2->ScalarVisibilityOff(); // 不显示颜色
//    vtkSmartPointer<vtkPolyDataMapper> liverMapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    liverMapper3->SetInputConnection(vtkReaderLiver3->GetOutputPort());
//    //liverMapper3->ScalarVisibilityOff(); // 不显示颜色
//    vtkSmartPointer<vtkPolyDataMapper> liverMapper4 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    liverMapper4->SetInputConnection(vtkReaderLiver4->GetOutputPort());
//    //liverMapper4->ScalarVisibilityOff(); // 不显示颜色
//
//    vtkSmartPointer<vtkPolyDataMapper> skinMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    skinMapper->SetInputConnection(vtkReaderSkin->GetOutputPort());
//    skinMapper->ScalarVisibilityOff(); // 不显示颜色 
//
//    vtkSmartPointer<vtkPolyDataMapper> arteryMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    arteryMapper->SetInputConnection(vtkReaderArtery->GetOutputPort());
//   // arteryMapper->ScalarVisibilityOff(); // 不显示颜色
//
//    vtkSmartPointer<vtkPolyDataMapper> boneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    boneMapper->SetInputConnection(vtkReaderBone->GetOutputPort());
//    boneMapper->ScalarVisibilityOff(); // 不显示颜色
//
//    // 创建演员
//    vtkSmartPointer<vtkActor> liverActor = vtkSmartPointer<vtkActor>::New();
//
//    liverActor->SetMapper(liverMapper);
//   liverActor->GetProperty()->SetOpacity(0.4);
//   liverActor->GetProperty()->SetColor(0.5, 0, 1);
//   vtkSmartPointer<vtkActor> liverActor2 = vtkSmartPointer<vtkActor>::New();
//
//   liverActor2->SetMapper(liverMapper2);
//   liverActor2->GetProperty()->SetOpacity(0.4);
//   liverActor2->GetProperty()->SetColor(0.5, 0, 1);
//   vtkSmartPointer<vtkActor> liverActor3 = vtkSmartPointer<vtkActor>::New();
//
//   liverActor3->SetMapper(liverMapper3);
//   liverActor3->GetProperty()->SetOpacity(0.4);
//   liverActor3->GetProperty()->SetColor(0.5, 0, 1);
//   vtkSmartPointer<vtkActor> liverActor4 = vtkSmartPointer<vtkActor>::New();
//
//   liverActor4->SetMapper(liverMapper4);
//   liverActor4->GetProperty()->SetOpacity(0.4);
//   liverActor4->GetProperty()->SetColor(0.5, 0, 1);
//
//   // 创建演员
//   vtkSmartPointer<vtkActor> skinActor = vtkSmartPointer<vtkActor>::New();
//   skinActor->SetMapper(skinMapper);
//   skinActor->GetProperty()->SetOpacity(0.3);
//
//    vtkSmartPointer<vtkActor> arteryActor = vtkSmartPointer<vtkActor>::New();
//    arteryActor->SetMapper(arteryMapper);
//    arteryActor->GetProperty()->SetColor(1, 0.5, 0);
//
//    vtkSmartPointer<vtkActor> boneActor = vtkSmartPointer<vtkActor>::New();
//    boneActor->SetMapper(boneMapper);
//
//    // 设置摄像机
//    vtkSmartPointer<vtkCamera> aCamera = vtkSmartPointer<vtkCamera>::New();
//    aCamera->SetViewUp(0, 0, 1);
//    aCamera->SetPosition(0, 0, 2);  // 改为适当的摄像机位置
//    aCamera->SetFocalPoint(0, 0, 0);
//    aCamera->ComputeViewPlaneNormal();
//    aCamera->Azimuth(45.0);
//    aCamera->Elevation(30.0);
//    aCamera->Dolly(1.5);
//
//    // 将所有演员添加到渲染器中
//    aRenderer->AddActor(lineActor);
//    aRenderer->AddActor(liverActor);
//    aRenderer->AddActor(liverActor2);
//    aRenderer->AddActor(liverActor3);
//    aRenderer->AddActor(liverActor4);
//    aRenderer->AddActor(arteryActor);
//    aRenderer->AddActor(boneActor);
//    aRenderer->AddActor(skinActor);
//
//    // 设置渲染器属性
//    aRenderer->SetActiveCamera(aCamera);
//    aRenderer->ResetCamera();
//    aRenderer->SetBackground(1.0, 1.0,1.0); // 背景色
//    aRenderer->ResetCameraClippingRange();
//
//    // 渲染并启动交互
//    renWin->Render();
//    iren->Initialize();
//    iren->Start();
//
//    return EXIT_SUCCESS;
//}
