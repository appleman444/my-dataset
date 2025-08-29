//#include<iostream>
//#include<cmath>
//using namespace std;
//
//class Point {
//public:
//    Point(double x = 0.0, double y = 0.0, double z = 0.0);
//    void disp();
//    friend class ComputeTools;
//private:
//    double x, y, z;
//};
//Point::Point(double x, double y, double z)
//{
//    this->x = x;
//    this->y = y;
//    this->z = z;
//}
//
//void Point::disp()
//{
//    cout << "point:(" << x << "," << y << "," << z << ")  ";
//}
//class Line
//{
//private:
//    double a, b, c, d;//ax+by+cz+d=0
//public:
//    Line(double a = 0.0, double b = 0.0, double c = 0.0, double d = 0.0);
//    void disp();
//    friend class ComputeTools;
//};
//Line::Line(double a, double b, double c, double d)
//{
//    this->a = a; this->b = b; this->c = c; this->d;
//}
//
//void Line::disp()
//{
//    cout << "line:" << a << "x+" << b << "y+" << c << "z" << "=" << d << "  ";
//}
////class Surface {
////private:
////    double a, b, c, d;
////public:
////    Surface(double a = 0.0, double b = 0.0,
////        double c = 0.0, double d = 0.0);
////    void disp();
////    friend class ComputeTools;
////};
////Surface::Surface(double a, double b, double c, double d)
////{
////    this->a = a; this->b = b; this->c = c; this->d = d;
////}
////void Surface::disp()
////{
////    cout << a << "x+" << b << "y+" << c << "z+" << d << "=0   ";
////}
//class ComputeTools {
//public:
//    //static double distance(Point p1, Point p2);
//    static double distance(Point p1, Line l1);
//    //static double distance(Point p1, Surface s1);
//    //static double distance(Line l1, Line l2);
//   // static double distance(Line l1, Surface s1);
//    //static double distance(Surface s1, Surface s2);
//};
////double ComputeTools::distance(Point p1, Point p2)
////{
////    return sqrt((p1.x - p2.x) * (p1.x - p2.x)
////        + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
////}
//double ComputeTools::distance(Point p1, Line l1)
//{
//    return fabs(l1.a * p1.x + l1.b * p1.y + l1.c * p1.z) / sqrt(l1.a * l1.a + l1.b * l1.b + l1.c * l1.c);
//}
////double ComputeTools::distance(Point p1, Surface s1)
////{
////    return fabs(s1.a * p1.x + s1.b * p1.y + s1.c * p1.z + s1.d) /
////        sqrt(s1.a * s1.a + s1.b * s1.b + s1.c * s1.c);
////}
////double ComputeTools::distance(Line l1, Line l2)
////{
////    return fabs(l1.c - l2.c) / sqrt(l1.a * l1.a + l1.b * l1.b  );
////}
//int main()
//{
//    Point p1(1.0, 1.0, 1.0), p2(0.0, 0.0, 0.0);
//    Line l1(1.0, 1.0, 1.0), 
//        l2(2.0, 2.0, 2.0);
//    //Surface s1(0.0, 0.0, 1.0);
//    //p1.disp();   p2.disp();
//   // cout << "  Distance: " << ComputeTools::distance(p1, p2) << endl;
//    p1.disp(); l1.disp();
//    cout << "  Distance: " << ComputeTools::distance(p1, l1) << endl;
//    //p1.disp(); s1.disp();
//   // cout << "  Distance: " << ComputeTools::distance(p1, s1) << endl;
//    //l1.disp(); l2.disp();
//   // cout << "  Distance: " << ComputeTools::distance(l1, l2) << endl;
//}