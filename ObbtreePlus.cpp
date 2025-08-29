//#include <vtkActor.h>
//#include <vtkCallbackCommand.h>
//#include <vtkCommand.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkMath.h>
//#include <vtkNamedColors.h>
//#include <vtkNew.h>
//#include <vtkOBBTree.h>
//#include <vtkPolyData.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkProperty.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkSliderRepresentation2D.h>
//#include <vtkSliderWidget.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataReader.h>
//
//namespace {
//    class vtkSliderCallback : public vtkCommand
//    {
//    public:
//        static vtkSliderCallback* New()
//        {
//            return new vtkSliderCallback;
//        }
//
//        vtkSliderCallback() : OBBTree(0), Level(0), PolyData(0), Renderer(0) {}
//
//        virtual void Execute(vtkObject* caller, unsigned long, void*)
//        {
//            vtkSliderWidget* sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
//            this->Level = vtkMath::Round(
//                static_cast<vtkSliderRepresentation*>(sliderWidget->GetRepresentation())
//                ->GetValue());
//            this->OBBTree->GenerateRepresentation(this->Level, this->PolyData);
//            this->Renderer->Render();
//        }
//
//        vtkOBBTree* OBBTree;
//        int Level;
//        vtkPolyData* PolyData;
//        vtkRenderer* Renderer;
//    };
//}
//
//int main(int argc, char* argv[])
//{
//    // 读取肝脏数据
//    vtkSmartPointer<vtkPolyDataReader> liverReader =
//        vtkSmartPointer<vtkPolyDataReader>::New();
//    liverReader->SetFileName("D://Data Disk//LW//Data//3Dircadb1.1//MyDo_VTK//bone.vtk");
//    liverReader->Update();
//    auto liverPolyData = liverReader->GetOutput();
//
//    // 读取肿瘤数据
//    vtkSmartPointer<vtkPolyDataReader> tumorReader =
//        vtkSmartPointer<vtkPolyDataReader>::New();
//    tumorReader->SetFileName("D://Data Disk//LW//Data//3Dircadb1.1//MyDo_VTK//livertumor01.vtk");
//    tumorReader->Update();
//    auto tumorPolyData = tumorReader->GetOutput();
//
//    vtkNew<vtkNamedColors> colors;
//
//    // 创建肝脏的Mapper和Actor
//    vtkNew<vtkPolyDataMapper> liverMapper;
//    liverMapper->SetInputData(liverPolyData);
//    liverMapper->ScalarVisibilityOff();
//
//    vtkNew<vtkActor> liverActor;
//    liverActor->SetMapper(liverMapper);
//    liverActor->GetProperty()->SetInterpolationToFlat();
//    liverActor->GetProperty()->SetColor(colors->GetColor4d("Yellow").GetData());
//    liverActor->GetProperty()->SetOpacity(.3);
//
//    // 创建肿瘤的Mapper和Actor
//    vtkNew<vtkPolyDataMapper> tumorMapper;
//    tumorMapper->SetInputData(tumorPolyData);
//    tumorMapper->ScalarVisibilityOff();
//
//    vtkNew<vtkActor> tumorActor;
//    tumorActor->SetMapper(tumorMapper);
//    tumorActor->GetProperty()->SetInterpolationToFlat();
//    tumorActor->GetProperty()->SetColor(colors->GetColor4d("Red").GetData());
//    tumorActor->GetProperty()->SetOpacity(.5);
//
//    int maxLevel = 5;
//
//    // 为肝脏创建OBB树
//    vtkNew<vtkOBBTree> liverOBBTree;
//    liverOBBTree->SetDataSet(liverPolyData);
//    liverOBBTree->SetMaxLevel(maxLevel);
//    liverOBBTree->BuildLocator();
//
//    // 为肿瘤创建OBB树
//    vtkNew<vtkOBBTree> tumorOBBTree;
//    tumorOBBTree->SetDataSet(tumorPolyData);
//    tumorOBBTree->SetMaxLevel(maxLevel);
//    tumorOBBTree->BuildLocator();
//
//    // 获取肝脏的包围盒信息
//    double liverCorner[3] = { 0.0, 0.0, 0.0 };
//    double liverMax[3] = { 0.0, 0.0, 0.0 };
//    double liverMid[3] = { 0.0, 0.0, 0.0 };
//    double liverMin[3] = { 0.0, 0.0, 0.0 };
//    double liverSize[3] = { 0.0, 0.0, 0.0 };
//
//    liverOBBTree->ComputeOBB(liverPolyData, liverCorner, liverMax, liverMid, liverMin, liverSize);
//
//    // 获取肿瘤的包围盒信息
//    double tumorCorner[3] = { 0.0, 0.0, 0.0 };
//    double tumorMax[3] = { 0.0, 0.0, 0.0 };
//    double tumorMid[3] = { 0.0, 0.0, 0.0 };
//    double tumorMin[3] = { 0.0, 0.0, 0.0 };
//    double tumorSize[3] = { 0.0, 0.0, 0.0 };
//
//    tumorOBBTree->ComputeOBB(tumorPolyData, tumorCorner, tumorMax, tumorMid, tumorMin, tumorSize);
//
//    std::cout << "Liver OBB Corner: " << liverCorner[0] << ", " << liverCorner[1] << ", " << liverCorner[2] << std::endl;
//    std::cout << "Tumor OBB Corner: " << tumorCorner[0] << ", " << tumorCorner[1] << ", " << tumorCorner[2] << std::endl;
//
//    // 初始化OBB树的Representation
//    vtkNew<vtkPolyData> liverOBBPolyData;
//    liverOBBTree->GenerateRepresentation(0, liverOBBPolyData);
//
//    vtkNew<vtkPolyDataMapper> liverOBBTreeMapper;
//    liverOBBTreeMapper->SetInputData(liverOBBPolyData);
//
//    vtkNew<vtkActor> liverOBBTreeActor;
//    liverOBBTreeActor->SetMapper(liverOBBTreeMapper);
//    liverOBBTreeActor->GetProperty()->SetInterpolationToFlat();
//    liverOBBTreeActor->GetProperty()->SetOpacity(.5);
//    liverOBBTreeActor->GetProperty()->EdgeVisibilityOn();
//    liverOBBTreeActor->GetProperty()->SetColor(colors->GetColor4d("SpringGreen").GetData());
//
//    vtkNew<vtkPolyData> tumorOBBPolyData;
//    tumorOBBTree->GenerateRepresentation(0, tumorOBBPolyData);
//
//    vtkNew<vtkPolyDataMapper> tumorOBBTreeMapper;
//    tumorOBBTreeMapper->SetInputData(tumorOBBPolyData);
//
//    vtkNew<vtkActor> tumorOBBTreeActor;
//    tumorOBBTreeActor->SetMapper(tumorOBBTreeMapper);
//    tumorOBBTreeActor->GetProperty()->SetInterpolationToFlat();
//    tumorOBBTreeActor->GetProperty()->SetOpacity(.5);
//    tumorOBBTreeActor->GetProperty()->EdgeVisibilityOn();
//    tumorOBBTreeActor->GetProperty()->SetColor(colors->GetColor4d("MediumPurple").GetData());
//
//    // 创建渲染器和窗口
//    vtkNew<vtkRenderer> renderer;
//    vtkNew<vtkRenderWindow> renderWindow;
//    renderWindow->AddRenderer(renderer);
//
//    vtkNew<vtkInteractorStyleTrackballCamera> style;
//    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
//    renderWindowInteractor->SetRenderWindow(renderWindow);
//    renderWindowInteractor->SetInteractorStyle(style);
//
//    // 添加Actors到渲染器
//    renderer->AddActor(liverActor);
//    renderer->AddActor(tumorActor);
//    renderer->AddActor(liverOBBTreeActor);
//    renderer->AddActor(tumorOBBTreeActor);
//    renderer->SetBackground(colors->GetColor3d("MidnightBlue").GetData());
//    renderer->UseHiddenLineRemovalOn();
//
//    // 渲染窗口
//    renderWindow->SetWindowName("VisualizeLiverAndTumor");
//    renderWindow->SetSize(600, 600);
//    renderWindow->Render();
//
//    vtkNew<vtkSliderRepresentation2D> sliderRep;
//    sliderRep->SetMinimumValue(0);
//    sliderRep->SetMaximumValue(maxLevel);
//    sliderRep->SetValue(maxLevel / 2);
//    sliderRep->SetTitleText("Level");
//    sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
//    sliderRep->GetPoint1Coordinate()->SetValue(.2, .2);
//    sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
//    sliderRep->GetPoint2Coordinate()->SetValue(.8, .2);
//    sliderRep->SetSliderLength(0.075);
//    sliderRep->SetSliderWidth(0.05);
//    sliderRep->SetEndCapLength(0.05);
//
//    vtkNew<vtkSliderWidget> sliderWidget;
//    sliderWidget->SetInteractor(renderWindowInteractor);
//    sliderWidget->SetRepresentation(sliderRep);
//    sliderWidget->SetAnimationModeToAnimate();
//    sliderWidget->EnabledOn();
//
//    vtkNew<vtkSliderCallback> callback;
//    callback->OBBTree = liverOBBTree;
//    callback->PolyData = liverOBBPolyData;
//    callback->Renderer = renderer;
//
//    sliderWidget->AddObserver(vtkCommand::InteractionEvent, callback);
//
//    renderWindowInteractor->Initialize();
//    renderWindow->Render();
//    renderWindowInteractor->Start();
//
//    return EXIT_SUCCESS;
//}
