#include "Point.h"

Point::Point(double x, double y) : x(x), y(y)
{
}

// Unary operators
Point Point::operator-() const
{
    return Point(-x, -y);
}

Point &Point::operator+=(const Point &other)
{
    x += other.x;
    y += other.y;
    return *this;
}

Point &Point::operator-=(const Point &other)
{
    x -= other.x;
    y -= other.y;
    return *this;
}

Point &Point::operator*=(double t)
{
    x *= t;
    y *= t;
    return *this;
}

bool Point::operator==(const Point &other) const
{
    return x == other.x && y == other.y;
}

bool Point::operator!=(const Point &other) const
{
    return x != other.x || y != other.y;
}

// Other operations
Point Point::getOrthogonal() const
{
    return Point(-y, x);
}

double Point::dot(const Point &other) const
{
    return x * other.x + y * other.y;
}

double Point::getNorm() const
{
    return std::sqrt(x * x + y * y);
}

double Point::getDistance(const Point &other) const
{
    return (*this - other).getNorm();
}

double Point::getDet(const Point &other) const
{
    return x * other.y - y * other.x;
}

Point operator+(Point lhs, const Point &rhs)
{
    lhs += rhs;
    return lhs;
}

Point operator-(Point lhs, const Point &rhs)
{
    lhs -= rhs;
    return lhs;
}

Point operator*(double t, Point vec)
{
    vec *= t;
    return vec;
}

Point operator*(Point vec, double t)
{
    return t * vec;
}

std::ostream &operator<<(std::ostream &os, const Point &vec)
{
    os << "(" << vec.x << ", " << vec.y << ")";
    return os;
}

#define MIN(x, y) x < y ? x : y
#define MAX(x, y) x > y ? x : y

bool Point::is_in_polygon(std::vector<Point> border)
{
    int jds = 0;
    double delta = 1e-10;
    for (int i = 0; i < border.size(); i++)
    {
        double maxx = MAX(border[i].x, border[(i + 1) % border.size()].x);
        double minx = MIN(border[i].x, border[(i + 1) % border.size()].x);
        double maxy = MAX(border[i].y, border[(i + 1) % border.size()].y);
        double miny = MIN(border[i].y, border[(i + 1) % border.size()].y);
        if (fabs(maxy - miny) > delta)
        {
            if (x <= minx + delta)
            {
                if (y <= maxy + delta && y > miny - delta)
                    jds++;
            }
            else if (minx - delta < x && x < maxx + delta)
            {
                if (y <= maxy + delta && y > miny - delta)
                {
                    double kk = (border[(i + 1) % border.size()].y - border[i].y) / (border[(i + 1) % border.size()].x - border[i].x);
                    double xx = border[i].x + (y - border[i].y) / kk;
                    if (xx > x - delta)
                        jds++;
                }
            }
        }
    }
    return (jds % 2 == 1);
}
