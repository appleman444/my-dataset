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
//    //// ������ά�ռ������㻭ֱ��
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
//    // ����·����Mapper��Actor
//    vtkSmartPointer<vtkPolyDataMapper> lineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    lineMapper->SetInputData(linesPolyData);
//
//    vtkSmartPointer<vtkActor> lineActor = vtkSmartPointer<vtkActor>::New();
//    lineActor->SetMapper(lineMapper);
//    lineActor->GetProperty()->SetColor(1.0, 0.0, 0.0);  // ��ɫ��ʾ·��
//    lineActor->GetProperty()->SetLineWidth(3);
//
//    // ��ȡ liver ����
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
//    // ��ȡ artery ����
//    vtkSmartPointer<vtkPolyDataReader> vtkReaderArtery = vtkSmartPointer<vtkPolyDataReader>::New();
//    vtkReaderArtery->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.1//MESHES_VTK//MESHES_VTK//artery.vtk");
//    vtkReaderArtery->Update();
//
//    vtkSmartPointer<vtkPolyDataReader> vtkReaderSkin = vtkSmartPointer<vtkPolyDataReader>::New();
//    vtkReaderSkin->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.1//MESHES_VTK//MESHES_VTK//skin.vtk");
//    vtkReaderSkin->Update();
//    // ��ȡ bone ����
//    vtkSmartPointer<vtkPolyDataReader> vtkReaderBone = vtkSmartPointer<vtkPolyDataReader>::New();
//    vtkReaderBone->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.1//MESHES_VTK//MESHES_VTK//bone.vtk");
//    vtkReaderBone->Update();
//
//    // ��������ӳ����
//    vtkSmartPointer<vtkPolyDataMapper> liverMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    liverMapper->SetInputConnection(vtkReaderLiver->GetOutputPort());
//   // liverMapper->ScalarVisibilityOff(); // ����ʾ��ɫ
//    vtkSmartPointer<vtkPolyDataMapper> liverMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    liverMapper2->SetInputConnection(vtkReaderLiver2->GetOutputPort());
//    //liverMapper2->ScalarVisibilityOff(); // ����ʾ��ɫ
//    vtkSmartPointer<vtkPolyDataMapper> liverMapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    liverMapper3->SetInputConnection(vtkReaderLiver3->GetOutputPort());
//    //liverMapper3->ScalarVisibilityOff(); // ����ʾ��ɫ
//    vtkSmartPointer<vtkPolyDataMapper> liverMapper4 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    liverMapper4->SetInputConnection(vtkReaderLiver4->GetOutputPort());
//    //liverMapper4->ScalarVisibilityOff(); // ����ʾ��ɫ
//
//    vtkSmartPointer<vtkPolyDataMapper> skinMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    skinMapper->SetInputConnection(vtkReaderSkin->GetOutputPort());
//    skinMapper->ScalarVisibilityOff(); // ����ʾ��ɫ 
//
//    vtkSmartPointer<vtkPolyDataMapper> arteryMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    arteryMapper->SetInputConnection(vtkReaderArtery->GetOutputPort());
//   // arteryMapper->ScalarVisibilityOff(); // ����ʾ��ɫ
//
//    vtkSmartPointer<vtkPolyDataMapper> boneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    boneMapper->SetInputConnection(vtkReaderBone->GetOutputPort());
//    boneMapper->ScalarVisibilityOff(); // ����ʾ��ɫ
//
//    // ������Ա
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
//   // ������Ա
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
//    // ���������
//    vtkSmartPointer<vtkCamera> aCamera = vtkSmartPointer<vtkCamera>::New();
//    aCamera->SetViewUp(0, 0, 1);
//    aCamera->SetPosition(0, 0, 2);  // ��Ϊ�ʵ��������λ��
//    aCamera->SetFocalPoint(0, 0, 0);
//    aCamera->ComputeViewPlaneNormal();
//    aCamera->Azimuth(45.0);
//    aCamera->Elevation(30.0);
//    aCamera->Dolly(1.5);
//
//    // ��������Ա��ӵ���Ⱦ����
//    aRenderer->AddActor(lineActor);
//    aRenderer->AddActor(liverActor);
//    aRenderer->AddActor(liverActor2);
//    aRenderer->AddActor(liverActor3);
//    aRenderer->AddActor(liverActor4);
//    aRenderer->AddActor(arteryActor);
//    aRenderer->AddActor(boneActor);
//    aRenderer->AddActor(skinActor);
//
//    // ������Ⱦ������
//    aRenderer->SetActiveCamera(aCamera);
//    aRenderer->ResetCamera();
//    aRenderer->SetBackground(1.0, 1.0,1.0); // ����ɫ
//    aRenderer->ResetCameraClippingRange();
//
//    // ��Ⱦ����������
//    renWin->Render();
//    iren->Initialize();
//    iren->Start();
//
//    return EXIT_SUCCESS;
//}
