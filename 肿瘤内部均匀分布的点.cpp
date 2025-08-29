//#include "vtkAutoInit.h"
//VTK_MODULE_INIT(vtkRenderingOpenGL2);
//VTK_MODULE_INIT(vtkInteractionStyle);
//
//#include <vtkSmartPointer.h>
//#include <vtkPolyData.h>
//#include <vtkPoints.h>
//#include <vtkPolyDataReader.h>
//#include <vtkSelectEnclosedPoints.h>
//#include <vtkIntArray.h>
//#include <vtkPointData.h>
//#include <vtkSphereSource.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkNamedColors.h>
//#include <vtkProperty.h>
//// 保存点集到 CSV 文件
//void SavePointsToCSV(vtkSmartPointer<vtkPoints> points, const std::string& fileName) {
//    // 打开文件，若文件不存在则创建
//    std::ofstream outFile(fileName);
//
//    if (!outFile.is_open()) {
//        std::cerr << "无法打开文件进行写入: " << fileName << std::endl;
//        return;
//    }
//
//    // 写入 CSV 表头
//    outFile << "X,Y,Z" << std::endl;
//
//    // 遍历点，写入每个点的坐标
//    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
//        double* point = points->GetPoint(i);  // 获取当前点的坐标
//        outFile << point[0] << "," << point[1] << "," << point[2] << std::endl;
//    }
//
//    // 关闭文件
//    outFile.close();
//    std::cout << "点已保存到文件: " << fileName << std::endl;
//}
//
//
//
//int main() {
//    // 读取肿瘤数据
//    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
//    reader->SetFileName("D://dataSet//3Dircadb1//3Dircadb1.17//MESHES_VTK//livertumor2.vtk");  // 假设肿瘤数据来自文件
//    reader->Update();
//    vtkSmartPointer<vtkPolyData> poly_tumor = reader->GetOutput();
//
//    // 获取肿瘤的包围盒
//    double tumorBounds[6];
//    poly_tumor->GetBounds(tumorBounds);
//    double tumorMinX = tumorBounds[0], tumorMaxX = tumorBounds[1];
//    double tumorMinY = tumorBounds[2], tumorMaxY = tumorBounds[3];
//    double tumorMinZ = tumorBounds[4], tumorMaxZ = tumorBounds[5];
//
//    // 生成均匀分布的 sample_pointsPolydata (在肿瘤边界范围内)
//    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//    vtkSmartPointer<vtkPolyData> sample_pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
//    int numPoints = 20000;  // 假设生成 1000 个点
//
//    // 将点生成范围限定在肿瘤的包围盒范围内
//    for (int i = 0; i < numPoints; i++) {
//        double x = tumorMinX + (tumorMaxX - tumorMinX) * (rand() / (RAND_MAX + 1.0));  // 随机生成 x 坐标
//        double y = tumorMinY + (tumorMaxY - tumorMinY) * (rand() / (RAND_MAX + 1.0));  // 随机生成 y 坐标
//        double z = tumorMinZ + (tumorMaxZ - tumorMinZ) * (rand() / (RAND_MAX + 1.0));  // 随机生成 z 坐标
//        points->InsertNextPoint(x, y, z);
//    }
//
//    // 将点集添加到 sample_pointsPolydata 中
//    sample_pointsPolydata->SetPoints(points);
//
//    // 使用 vtkSelectEnclosedPoints 判断点是否在肿瘤内部
//    vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
//    enclosedPoints->SetInputData(sample_pointsPolydata);  // 设置要判断的点数据
//    enclosedPoints->SetSurfaceData(poly_tumor);           // 设置肿瘤的表面数据
//    enclosedPoints->Update();
//
//    // 获取每个点是否在肿瘤内部
//    vtkSmartPointer<vtkDataArray> insideArray = enclosedPoints->GetOutput()->GetPointData()->GetArray("SelectedPoints");
//
//    // 输出肿瘤内部点的数量
//    std::vector<vtkIdType> insidePoints;  // 存储肿瘤内部的点的索引
//    std::vector<vtkIdType> outPoints;  // 存储肿瘤内部的点的索引
//    // 创建 vtkPoints 对象来保存肿瘤内部的点
//    vtkSmartPointer<vtkPoints> insidePointsData = vtkSmartPointer<vtkPoints>::New();
//    for (vtkIdType i = 0; i < sample_pointsPolydata->GetNumberOfPoints(); i++) {
//        // 使用 GetComponent 来获取是否在肿瘤内的信息
//        bool insideTumor = insideArray->GetComponent(i, 0) == 1;  // 判断点是否在肿瘤内
//        double inpoint[3];
//        sample_pointsPolydata->GetPoint(i, inpoint);  // 获取点的坐标
//
//        if (insideTumor) {
//            insidePoints.push_back(i);  // 如果点在肿瘤内，将其索引加入列表
//            insidePointsData->InsertNextPoint(inpoint);  // 将该点插入到 vtkPoints 对象中
//        }
//        else {
//            outPoints.push_back(i);  // 如果点在肿瘤内，将其索引加入列表
//        }
//    }
//
//    std::cout << "Points inside the tumor: " << insidePoints.size() << std::endl;
//    vtkSmartPointer<vtkPolyData> insidePolyData = vtkSmartPointer<vtkPolyData>::New();
//    insidePolyData->SetPoints(insidePointsData);
//
//    // 在 main 函数中修改文件路径
//    std::string fileName = "D:/Data Disk/LW/Data/3Dircadb1.1/CSV/tumorinpoints/livertumor2-17.csv";  // 保存文件路径
//    SavePointsToCSV(insidePolyData->GetPoints(), fileName);  // 保存肿瘤内部的点到指定路径
//
//    // 可视化肿瘤内部点的小球
//    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//    renderWindow->AddRenderer(renderer);
//
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    renderWindowInteractor->SetRenderWindow(renderWindow);
//
//    vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
//
//    // 为肿瘤内部的每个点创建小球
//    for (size_t i = 0; i < insidePoints.size(); i++) {
//        double point[3];
//        sample_pointsPolydata->GetPoint(insidePoints[i], point);
//
//        // 创建小球
//        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
//        sphereSource->SetCenter(point);
//        sphereSource->SetRadius(0.2);  // 设置小球半径为0.5
//
//        vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//        sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
//
//        vtkSmartPointer<vtkActor> sphereActor = vtkSmartPointer<vtkActor>::New();
//        sphereActor->SetMapper(sphereMapper);
//
//        // 设置小球颜色
//        sphereActor->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());
//
//        // 将小球添加到渲染器
//        renderer->AddActor(sphereActor);
//    }
//   // 为肿瘤外部的每个点创建小球
//    //for (size_t i = 0; i < outPoints.size(); i++) {
//    //    double point1[3];
//    //    sample_pointsPolydata->GetPoint(outPoints[i], point1);
//
//    //    // 创建小球
//    //    vtkSmartPointer<vtkSphereSource> sphereSource1 = vtkSmartPointer<vtkSphereSource>::New();
//    //    sphereSource1->SetCenter(point1);
//    //    sphereSource1->SetRadius(0.15);  // 设置小球半径为0.5
//
//    //    vtkSmartPointer<vtkPolyDataMapper> sphereMapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//    //    sphereMapper1->SetInputConnection(sphereSource1->GetOutputPort());
//
//    //    vtkSmartPointer<vtkActor> sphereActor1= vtkSmartPointer<vtkActor>::New();
//    //    sphereActor1->SetMapper(sphereMapper1);
//
//    //    // 设置小球颜色
//    //    sphereActor1->GetProperty()->SetColor(colors->GetColor3d("Green").GetData());
//
//    //    // 将小球添加到渲染器
//    //    renderer->AddActor(sphereActor1);
//    //}
//
//    // 可视化肿瘤（设置绿色、透明度0.4）
//    vtkSmartPointer<vtkPolyDataMapper> tumorMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//    tumorMapper->SetInputData(poly_tumor);
//
//    vtkSmartPointer<vtkActor> tumorActor = vtkSmartPointer<vtkActor>::New();
//    tumorActor->SetMapper(tumorMapper);
//
//    // 设置肿瘤的颜色和透明度
//    tumorActor->GetProperty()->SetColor(1, 1, 1);  // 绿色
//    tumorActor->GetProperty()->SetOpacity(0.4);  // 透明度0.4
//
//    // 将肿瘤添加到渲染器
//   renderer->AddActor(tumorActor);
//
//    // 设置背景颜色
//    renderer->SetBackground(colors->GetColor3d("White").GetData());
//
//    // 启动渲染
//    renderWindow->Render();
//    renderWindowInteractor->Start();
//
//    return 0;
//}
