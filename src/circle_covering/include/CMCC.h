#pragma once
// STL
#include <cassert>
#include <list>
#include <array>

// My includes
#include "Point.h"

// MCC
inline double sqr(double r);
class MCC
{
public:
    Point center;
    double radius;
    MCC(std::vector<Point> subject);

private:
    std::vector<Point> mSubject;
    //...for the algorithms
    std::list<Point> mSupportPoints;
    std::_List_iterator<Point> support_end;
    int fsize; // number of forced points
    int ssize; // number of support points

    // ...for the ball updates
    Point mCenter;
    double mRadius;
    Point c[3];
    double sqr_r[3];

    // helper arrays
    Point q0;
    double z[3];
    double f[3];
    Point v[3];
    double a[3][2];

    void mtf_mb(std::_List_iterator<Point> n);
    void mtf_move_to_front(std::_List_iterator<Point> j);
    void computer_mcc();
    void pivot_move_to_front(Point p);
    double excess(Point p) const;
    void pop();
    bool push(Point p);
};

// CMCC
class CMCC
{
public:
    Point center;
    double radius;
    CMCC(std::vector<Point> subject, std::vector<Point> border);
    CMCC(std::vector<Point> subject, std::vector<Point> border, std::vector<std::vector<Point>> holes);

private:
    std::vector<Point> mSubject;
    std::vector<Point> mBorder;
    std::vector<std::vector<Point>> mHoles;

    // find convex hull
    std::vector<Point> convexhull(std::vector<Point> points);

    // find cmcc
    Point point_to_line_segment(Point a, Point b, Point point) const;
    Point bisect_line_segment(Point a, Point b, Point c, Point d) const;
    double find_y_on_line(Point a, Point b, double x) const;
    double find_x_on_line(Point a, Point b, double y) const;
    double dis_point_to_line(Point a, Point b, Point point) const;
    double dis_point_to_line_segment(const Point a, const Point b, const Point point);
    Point coord_point_to_line(Point a, Point b, Point point) const;
    Point coord_point_to_line_segment(Point a, Point b, Point point) const;
    std::pair<Point, double> Megiddo(Point a, Point b, std::vector<Point> subject);
    std::pair<Point, double> find_cmcc(std::vector<Point> subject, std::vector<Point> border, Point MCCcenter, double MCCradius);
};