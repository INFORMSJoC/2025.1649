#pragma once
// STL
#include <ostream>
#include <vector>
#include <algorithm>
#include <cmath>

// Declarations
class Point;
Point operator-(Point lhs, const Point &rhs);

// Implementations
class Point
{
public:
    double x;
    double y;

    Point(double x = 0.0, double y = 0.0);

    // Unary operators
    Point operator-() const;
    Point &operator+=(const Point &other);
    Point &operator-=(const Point &other);
    Point &operator*=(double t);
    bool operator==(const Point &other) const;
    bool operator!=(const Point &other) const;

    // Other operations
    Point getOrthogonal() const;
    double dot(const Point &other) const;
    double getNorm() const;
    double getDistance(const Point &other) const;
    double getDet(const Point &other) const;
    bool is_in_polygon(std::vector<Point> border);
};

// Binary operators
Point operator+(Point lhs, const Point &rhs);
Point operator-(Point lhs, const Point &rhs);
Point operator*(double t, Point vec);
Point operator*(Point vec, double t);
std::ostream &operator<<(std::ostream &os, const Point &vec);
