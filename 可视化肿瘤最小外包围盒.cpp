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
//    // 1. 读取肿瘤的VTK文件（假设文件名为livertumor01.vtk）
//    std::string inputFileName = "D:/Data Disk/LW/Data/3Dircadb1.1/MESHES_VTK/livertumor03.vtk";  // 更新为实际文件路径
//    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
//    reader->SetFileName(inputFileName.c_str());
//    reader->Update();
//
//    // 2. 检查是否成功加载VTK文件
//    if (reader->GetOutput()->GetNumberOfCells() == 0) {
//        std::cerr << "Error: No data loaded from the VTK file." << std::endl;
//        return -1;
//    }
//
//    // 3. 计算最小外包围盒（使用vtkBoundingBox）
//    vtkBoundingBox boundingBox(reader->GetOutput()->GetBounds());
//
//    // 4. 获取最小值点 P1 和 最大值点 P2
//    double P1[3] = { boundingBox.GetMinPoint()[0], boundingBox.GetMinPoint()[1], boundingBox.GetMinPoint()[2] };
//    double P2[3] = { boundingBox.GetMaxPoint()[0], boundingBox.GetMaxPoint()[1], boundingBox.GetMaxPoint()[2] };
//
//    std::cout << "最小值点 P1: (" << P1[0] << ", " << P1[1] << ", " << P1[2] << ")" << std::endl;
//    std::cout << "最大值点 P2: (" << P2[0] << ", " << P2[1] << ", " << P2[2] << ")" << std::endl;
//
//    // 5. 可视化肿瘤网格
//    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    mapper->SetInputData(reader->GetOutput());
//    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//    actor->SetMapper(mapper);
//    actor->GetProperty()->SetColor(1.0, 0.784, 0.588);
//    actor->GetProperty()->SetOpacity(0.7);
//
//    // 6. 可视化最小值点和最大值点
//    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
//    sphereSource->SetRadius(0.1);
//
//    // 可视化最小值点
//    vtkSmartPointer<vtkGlyph3D> glyph1 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph1->SetSourceConnection(sphereSource->GetOutputPort());
//    vtkSmartPointer<vtkPoints> points1 = vtkSmartPointer<vtkPoints>::New();
//    points1->InsertNextPoint(P1);
//    vtkSmartPointer<vtkPolyData> polyData1 = vtkSmartPointer<vtkPolyData>::New();
//    polyData1->SetPoints(points1);
//    glyph1->SetInputData(polyData1);
//
//    // 创建一个Actor来显示最小值点
//    vtkSmartPointer<vtkPolyDataMapper> glyphMapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    glyphMapper1->SetInputConnection(glyph1->GetOutputPort());
//    vtkSmartPointer<vtkActor> glyphActor1 = vtkSmartPointer<vtkActor>::New();
//    glyphActor1->SetMapper(glyphMapper1);
//
//    // 可视化最大值点
//    vtkSmartPointer<vtkGlyph3D> glyph2 = vtkSmartPointer<vtkGlyph3D>::New();
//    glyph2->SetSourceConnection(sphereSource->GetOutputPort());
//    vtkSmartPointer<vtkPoints> points2 = vtkSmartPointer<vtkPoints>::New();
//    points2->InsertNextPoint(P2);
//    vtkSmartPointer<vtkPolyData> polyData2 = vtkSmartPointer<vtkPolyData>::New();
//    polyData2->SetPoints(points2);
//    glyph2->SetInputData(polyData2);
//
//    // 创建一个Actor来显示最大值点
//    vtkSmartPointer<vtkPolyDataMapper> glyphMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    glyphMapper2->SetInputConnection(glyph2->GetOutputPort());
//    vtkSmartPointer<vtkActor> glyphActor2 = vtkSmartPointer<vtkActor>::New();
//    glyphActor2->SetMapper(glyphMapper2);
//
//    // 7. 创建外包围盒（使用vtkBox）
//    vtkSmartPointer<vtkPolyData> boxPolyData = vtkSmartPointer<vtkPolyData>::New();
//    vtkSmartPointer<vtkCellArray> boxCells = vtkSmartPointer<vtkCellArray>::New();
//
//    // 定义8个角点
//    double min[3], max[3];
//    boundingBox.GetMinPoint(min);
//    boundingBox.GetMaxPoint(max);
//
//    vtkIdType pts[8] = {
//        0, 1, 2, 3, 4, 5, 6, 7
//    };
//
//    // 创建外包围盒的8个点
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
//    // 定义外包围盒的12条边
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
//    // 创建外包围盒的Actor并设置颜色
//    vtkSmartPointer<vtkPolyDataMapper> boxMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    boxMapper->SetInputData(boxPolyData);
//    vtkSmartPointer<vtkActor> boxActor = vtkSmartPointer<vtkActor>::New();
//    boxActor->SetMapper(boxMapper);
//    boxActor->GetProperty()->SetColor(0.0, 0.0, 0.0); // 红色外包围盒
//    boxActor->GetProperty()->SetLineWidth(3);  // 设置外包围盒的线宽
//
//    // 8. 创建文本演员，并设置文本属性
//    vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
//    vtkSmartPointer<vtkTextProperty> textProperty = vtkSmartPointer<vtkTextProperty>::New();
//    textProperty->SetFontSize(24);
//    textProperty->SetColor(1.0, 0.0, 0.0); // 设置字体颜色为白色
//    textActor->SetTextProperty(textProperty);
//
//    // 设置文本内容
//   // textActor->SetText("最小外包围盒 P1 和 P2");
//
//    // 设置文本位置
//    textActor->SetPosition(10, 10);  // 设置在窗口中的位置
//
//    // 9. 创建渲染器、窗口和交互器
//    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//    renderer->AddActor(actor);
//    renderer->AddActor(glyphActor1);
//    renderer->AddActor(glyphActor2);
//    renderer->AddActor(boxActor);
//    renderer->AddActor(textActor); // 添加文本到渲染器
//    renderer->SetBackground(1,1,1);//设置背景颜色，double类型
//    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//    renderWindow->AddRenderer(renderer);
//
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    renderWindowInteractor->SetRenderWindow(renderWindow);
//
//    // 启动渲染
//    renderWindow->Render();
//    renderWindowInteractor->Start();
//
//    return 0;
//}
