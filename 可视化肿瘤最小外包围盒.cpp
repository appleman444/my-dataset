//#include <vtkAutoInit.h>
//VTK_MODULE_INIT(vtkRenderingOpenGL2);
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingFreeType);
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataReader.h>
//#include <vtkBoundingBox.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkSphereSource.h>
//#include <vtkGlyph3D.h>
//#include <vtkPoints.h>
//#include <vtkPolyData.h>
//#include <vtkBox.h>
//#include <vtkProperty.h>
//#include <vtkCellArray.h>
//#include <vtkCamera.h>
//#include <vtkTextActor.h>
//#include <vtkTextProperty.h>
//
//int main(int argc, char* argv[])
//{
//    // 1. ��ȡ������VTK�ļ��������ļ���Ϊlivertumor01.vtk��
//    std::string inputFileName = "D:/Data Disk/LW/Data/3Dircadb1.1/MESHES_VTK/livertumor03.vtk";  // ����Ϊʵ���ļ�·��
//    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
//    reader->SetFileName(inputFileName.c_str());
//    reader->Update();
//
//    // 2. ����Ƿ�ɹ�����VTK�ļ�
//    if (reader->GetOutput()->GetNumberOfCells() == 0) {
//        std::cerr << "Error: No data loaded from the VTK file." << std::endl;
//        return -1;
//    }
//
//    // 3. ������С���Χ�У�ʹ��vtkBoundingBox��
//    vtkBoundingBox boundingBox(reader->GetOutput()->GetBounds());
//
//    // 4. ��ȡ��Сֵ�� P1 �� ���ֵ�� P2
//    double P1[3] = { boundingBox.GetMinPoint()[0], boundingBox.GetMinPoint()[1], boundingBox.GetMinPoint()[2] };
//    double P2[3] = { boundingBox.GetMaxPoint()[0], boundingBox.GetMaxPoint()[1], boundingBox.GetMaxPoint()[2] };
//
//    std::cout << "��Сֵ�� P1: (" << P1[0] << ", " << P1[1] << ", " << P1[2] << ")" << std::endl;
//    std::cout << "���ֵ�� P2: (" << P2[0] << ", " << P2[1] << ", " << P2[2] << ")" << std::endl;
//
//    // 5. ���ӻ���������
//    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    mapper->SetInputData(reader->GetOutput());
//    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//    actor->SetMapper(mapper);
//    actor->GetProperty()->SetColor(1.0, 0.784, 0.588);
//    actor->GetProperty()->SetOpacity(0.7);
//
//    // 6. ���ӻ���Сֵ������ֵ��
//    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
//    sphereSource->SetRadius(0.1);
//
//    // ���ӻ���Сֵ��
//    vtkSmartPointer<vtkGlyph3D> glyph1 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph1->SetSourceConnection(sphereSource->GetOutputPort());
//    vtkSmartPointer<vtkPoints> points1 = vtkSmartPointer<vtkPoints>::New();
//    points1->InsertNextPoint(P1);
//    vtkSmartPointer<vtkPolyData> polyData1 = vtkSmartPointer<vtkPolyData>::New();
//    polyData1->SetPoints(points1);
//    glyph1->SetInputData(polyData1);
//
//    // ����һ��Actor����ʾ��Сֵ��
//    vtkSmartPointer<vtkPolyDataMapper> glyphMapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    glyphMapper1->SetInputConnection(glyph1->GetOutputPort());
//    vtkSmartPointer<vtkActor> glyphActor1 = vtkSmartPointer<vtkActor>::New();
//    glyphActor1->SetMapper(glyphMapper1);
//
//    // ���ӻ����ֵ��
//    vtkSmartPointer<vtkGlyph3D> glyph2 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph2->SetSourceConnection(sphereSource->GetOutputPort());
//    vtkSmartPointer<vtkPoints> points2 = vtkSmartPointer<vtkPoints>::New();
//    points2->InsertNextPoint(P2);
//    vtkSmartPointer<vtkPolyData> polyData2 = vtkSmartPointer<vtkPolyData>::New();
//    polyData2->SetPoints(points2);
//    glyph2->SetInputData(polyData2);
//
//    // ����һ��Actor����ʾ���ֵ��
//    vtkSmartPointer<vtkPolyDataMapper> glyphMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    glyphMapper2->SetInputConnection(glyph2->GetOutputPort());
//    vtkSmartPointer<vtkActor> glyphActor2 = vtkSmartPointer<vtkActor>::New();
//    glyphActor2->SetMapper(glyphMapper2);
//
//    // 7. �������Χ�У�ʹ��vtkBox��
//    vtkSmartPointer<vtkPolyData> boxPolyData = vtkSmartPointer<vtkPolyData>::New();
//    vtkSmartPointer<vtkCellArray> boxCells = vtkSmartPointer<vtkCellArray>::New();
//
//    // ����8���ǵ�
//    double min[3], max[3];
//    boundingBox.GetMinPoint(min);
//    boundingBox.GetMaxPoint(max);
//
//    vtkIdType pts[8] = {
//        0, 1, 2, 3, 4, 5, 6, 7
//    };
//
//    // �������Χ�е�8����
//    vtkSmartPointer<vtkPoints> boxPoints = vtkSmartPointer<vtkPoints>::New();
//    boxPoints->InsertNextPoint(min[0], min[1], min[2]);
//    boxPoints->InsertNextPoint(max[0], min[1], min[2]);
//    boxPoints->InsertNextPoint(max[0], max[1], min[2]);
//    boxPoints->InsertNextPoint(min[0], max[1], min[2]);
//    boxPoints->InsertNextPoint(min[0], min[1], max[2]);
//    boxPoints->InsertNextPoint(max[0], min[1], max[2]);
//    boxPoints->InsertNextPoint(max[0], max[1], max[2]);
//    boxPoints->InsertNextPoint(min[0], max[1], max[2]);
//
//    boxPolyData->SetPoints(boxPoints);
//
//    // �������Χ�е�12����
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(0); boxCells->InsertCellPoint(1);
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(1); boxCells->InsertCellPoint(2);
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(2); boxCells->InsertCellPoint(3);
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(3); boxCells->InsertCellPoint(0);
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(4); boxCells->InsertCellPoint(5);
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(5); boxCells->InsertCellPoint(6);
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(6); boxCells->InsertCellPoint(7);
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(7); boxCells->InsertCellPoint(4);
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(0); boxCells->InsertCellPoint(4);
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(1); boxCells->InsertCellPoint(5);
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(2); boxCells->InsertCellPoint(6);
//    boxCells->InsertNextCell(2); boxCells->InsertCellPoint(3); boxCells->InsertCellPoint(7);
//
//    boxPolyData->SetLines(boxCells);
//
//    // �������Χ�е�Actor��������ɫ
//    vtkSmartPointer<vtkPolyDataMapper> boxMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    boxMapper->SetInputData(boxPolyData);
//    vtkSmartPointer<vtkActor> boxActor = vtkSmartPointer<vtkActor>::New();
//    boxActor->SetMapper(boxMapper);
//    boxActor->GetProperty()->SetColor(0.0, 0.0, 0.0); // ��ɫ���Χ��
//    boxActor->GetProperty()->SetLineWidth(3);  // �������Χ�е��߿�
//
//    // 8. �����ı���Ա���������ı�����
//    vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
//    vtkSmartPointer<vtkTextProperty> textProperty = vtkSmartPointer<vtkTextProperty>::New();
//    textProperty->SetFontSize(24);
//    textProperty->SetColor(1.0, 0.0, 0.0); // ����������ɫΪ��ɫ
//    textActor->SetTextProperty(textProperty);
//
//    // �����ı�����
//   // textActor->SetText("��С���Χ�� P1 �� P2");
//
//    // �����ı�λ��
//    textActor->SetPosition(10, 10);  // �����ڴ����е�λ��
//
//    // 9. ������Ⱦ�������ںͽ�����
//    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//    renderer->AddActor(actor);
//    renderer->AddActor(glyphActor1);
//    renderer->AddActor(glyphActor2);
//    renderer->AddActor(boxActor);
//    renderer->AddActor(textActor); // ����ı�����Ⱦ��
//    renderer->SetBackground(1,1,1);//���ñ�����ɫ��double����
//    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//    renderWindow->AddRenderer(renderer);
//
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    renderWindowInteractor->SetRenderWindow(renderWindow);
//
//    // ������Ⱦ
//    renderWindow->Render();
//    renderWindowInteractor->Start();
//
//    return 0;
//}
