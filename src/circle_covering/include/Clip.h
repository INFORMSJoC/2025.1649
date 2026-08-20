#pragma once
// STL
#include <math.h>
#include <iostream>
#include <chrono>

// My includes
#include "CMCC.h"
#include "clipper.h"

class Clip
{
public:
    Clip(std::vector<std::vector<Point>> subject, std::vector<Point> border);
    Clip(std::vector<std::vector<Point>> subject, std::vector<Point> border, std::vector<std::vector<Point>> holes);
    Clip(std::vector<std::vector<Point>> subject, std::vector<std::vector<std::vector<Point>>> borders, std::vector<std::vector<Point>> holes);
    std::vector<std::vector<Point>> clip();
    std::vector<std::vector<std::vector<Point>>> clip_border();
    std::vector<std::vector<Point>> clip_acc();

private:
    std::vector<std::vector<Point>> mSubject;
    std::vector<std::vector<std::vector<Point>>> mBorders;
    std::vector<Point> mBorder;
    std::vector<std::vector<Point>> mHoles;
};

// Simple Clip
enum Sign
{
    /** Value is negative */
    NEGATIVE = -1,
    /** Value is zero */
    ZERO = 0,
    /** Value is positive */
    POSITIVE = 1
};

bool are_points_collinear(const Point a, const Point b, const Point c);
std::vector<Point> convexhull(std::vector<Point> points);
static inline Sign point_is_in_half_plane(const Point &p, const Point &q1, const Point &q2);
static inline bool intersect_segments(const Point &p1, const Point &p2, const Point &q1, const Point &q2, Point &result);
void clip_polygon_by_half_plane(const std::vector<Point> &P, const Point &q1, const Point &q2, std::vector<Point> &result);
std::vector<Point> convex_clip_polygon(const std::vector<Point> &P, const std::vector<Point> &clip);
