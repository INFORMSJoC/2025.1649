#include "Clip.h"

Clip::Clip(std::vector<std::vector<Point>> subject, std::vector<Point> border) : mSubject(subject), mBorder(border)
{
}
Clip::Clip(std::vector<std::vector<Point>> subject, std::vector<Point> border, std::vector<std::vector<Point>> holes) : mSubject(subject), mBorder(border), mHoles(holes)
{
}
Clip::Clip(std::vector<std::vector<Point>> subject, std::vector<std::vector<std::vector<Point>>> borders, std::vector<std::vector<Point>> holes) : mSubject(subject), mBorders(borders), mHoles(holes)
{
}

std::vector<std::vector<Point>> Clip::clip()
{
    std::vector<std::vector<Point>> voronoi_new;
    Clipper2Lib::Paths clip_border, clip_holes;
    clip_border.push_back(Clipper2Lib::MakePath(mBorder));
    for (size_t i = 0; i < mHoles.size(); i++)
    {
        clip_holes.push_back(Clipper2Lib::MakePath(mHoles[i]));
    }
    for (size_t i = 0; i < mSubject.size(); i++)
    {
        Clipper2Lib::Paths subject_border, solution_border, solution_holes;
        subject_border.push_back(Clipper2Lib::MakePath(mSubject[i]));
        // clip border
        solution_border = Clipper2Lib::Intersect(subject_border, clip_border, Clipper2Lib::FillRule::NonZero);
        // clip holes
        solution_holes = Clipper2Lib::Difference(solution_border, clip_holes, Clipper2Lib::FillRule::NonZero);
        // save results
        std::vector<Point> result;
        auto start3 = std::chrono::steady_clock::now();
        for (size_t j = 0; j < solution_holes.size(); j++)
        {
            for (size_t k = 0; k < solution_holes[j].size(); k++)
            {
                result.push_back(solution_holes[j][k]);
            }
        }
        voronoi_new.push_back(result);
    }
    return voronoi_new;
}

std::vector<std::vector<std::vector<Point>>> Clip::clip_border()
{
    std::vector<std::vector<std::vector<Point>>> sub_boders;
    Clipper2Lib::Paths clip_border, clip_holes;
    clip_border.push_back(Clipper2Lib::MakePath(mBorder));
    for (size_t i = 0; i < mHoles.size(); i++)
    {
        clip_holes.push_back(Clipper2Lib::MakePath(mHoles[i]));
    }
    for (size_t i = 0; i < mSubject.size(); i++)
    {
        Clipper2Lib::Paths subject_border, solution_border, solution_holes;
        subject_border.push_back(Clipper2Lib::MakePath(mSubject[i]));
        solution_border = Clipper2Lib::Intersect(subject_border, clip_border, Clipper2Lib::FillRule::NonZero);
        solution_holes = Clipper2Lib::Difference(solution_border, clip_holes, Clipper2Lib::FillRule::NonZero);
        sub_boders.push_back(solution_holes);
    }
    return sub_boders;
}

std::vector<std::vector<Point>> Clip::clip_acc()
{
    std::vector<std::vector<Point>> voronoi_new;
    Clipper2Lib::Paths clip_border;
    for (size_t i = 0; i < mSubject.size(); i++)
    {
        // clip border
        clip_border.clear();
        for (size_t j = 0; j < mBorders[i].size(); j++)
            clip_border.push_back(Clipper2Lib::MakePath(mBorders[i][j]));
        Clipper2Lib::Paths subject_border, solution_border;
        subject_border.push_back(Clipper2Lib::MakePath(mSubject[i]));
        solution_border = Clipper2Lib::Intersect(subject_border, clip_border, Clipper2Lib::FillRule::NonZero);
        // save results
        std::vector<Point> result;
        auto start3 = std::chrono::steady_clock::now();
        for (size_t j = 0; j < solution_border.size(); j++)
        {
            for (size_t k = 0; k < solution_border[j].size(); k++)
                result.push_back(solution_border[j][k]);
        }
        voronoi_new.push_back(result);
    }
    return voronoi_new;
}

// Simple Clip
bool are_points_collinear(const Point a, const Point b, const Point c)
{
    return fabs((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)) < 1e-8;
}

std::vector<Point> convexhull(std::vector<Point> points)
{
    int n = points.size();
    if (n <= 3)
        return points;

    sort(points.begin(), points.end(), [](const Point &a, const Point &b)
         { return (a.x < b.x || (a.x == b.x && a.y < b.y)); });

    std::vector<Point> lower;
    for (auto it = points.begin(); it != points.end(); ++it)
    {
        // Pop off any points that make a convex angle with *it
        while (lower.size() >= 2 && (((*lower.rbegin()).x - (*(lower.rbegin() + 1)).x) * ((*it).y - (*(lower.rbegin() + 1)).y) - ((*lower.rbegin()).y - (*(lower.rbegin() + 1)).y) * ((*it).x - (*(lower.rbegin() + 1)).x)) <= 0)
            lower.pop_back();
        lower.push_back(*it);
    }

    std::vector<Point> upper;
    for (auto it = points.rbegin(); it != points.rend(); ++it)
    {
        // Pop off any points that make a convex angle with *it
        while (upper.size() >= 2 && (((*upper.rbegin()).x - (*(upper.rbegin() + 1)).x) * ((*it).y - (*(upper.rbegin() + 1)).y) - ((*upper.rbegin()).y - (*(upper.rbegin() + 1)).y) * ((*it).x - (*(upper.rbegin() + 1)).x)) <= 0)
            upper.pop_back();
        upper.push_back(*it);
    }

    std::vector<Point> hull;
    hull.insert(hull.end(), lower.begin(), lower.end());
    // Both hulls include both endpoints, so leave them out when we
    // append the upper hull.
    hull.insert(hull.end(), upper.begin() + 1, upper.end() - 1);

    return hull;
}

static inline Sign point_is_in_half_plane(const Point &p, const Point &q1, const Point &q2)
{
    Point a = q1 - p;
    Point b = q2 - p;
    double Delta = a.getDet(b);
    int int_tmp_result;
    double eps;
    double max1 = fabs(a.x);
    if ((max1 < fabs(a.y)))
        max1 = fabs(a.y);
    double max2 = fabs(b.x);
    if ((max2 < fabs(b.y)))
        max2 = fabs(b.y);
    double lower_bound_1;
    double upper_bound_1;
    lower_bound_1 = max1;
    upper_bound_1 = max1;
    if ((max2 < lower_bound_1))
        lower_bound_1 = max2;
    else
    {
        if ((max2 > upper_bound_1))
            upper_bound_1 = max2;
    }
    if ((lower_bound_1 < 5.00368081960964635413e-147))
        int_tmp_result = 0;
    else
    {
        if ((upper_bound_1 > 1.67597599124282407923e+153))
            int_tmp_result = 0;
        eps = (8.88720573725927976811e-16 * (max1 * max2));
        if ((Delta > eps))
            int_tmp_result = 1;
        else
        {
            if ((Delta < -eps))
                int_tmp_result = -1;
            else
                int_tmp_result = 0;
        }
    }
    return Sign(int_tmp_result);
}

static inline bool intersect_segments(const Point &p1, const Point &p2, const Point &q1, const Point &q2, Point &result)
{
    Point Vp = p2 - p1;
    Point Vq = q2 - q1;
    Point pq = q1 - p1;
    double a = Vp.x;
    double b = -Vq.x;
    double c = Vp.y;
    double d = -Vq.y;
    double delta = a * d - b * c;
    if (fabs(delta) <= 1e-10)
        return false;

    double tp = (d * pq.x - b * pq.y) / delta;
    result = Point((1.0 - tp) * p1.x + tp * p2.x, (1.0 - tp) * p1.y + tp * p2.y);

    return true;
}

void clip_polygon_by_half_plane(const std::vector<Point> &P, const Point &q1, const Point &q2, std::vector<Point> &result)
{
    result.clear();
    if (P.size() == 0)
        return;

    if (P.size() == 1)
    {
        if (point_is_in_half_plane(P[0], q1, q2))
            result.push_back(P[0]);
        return;
    }

    Point prev_p = P[P.size() - 1];
    Sign prev_status = point_is_in_half_plane(prev_p, q1, q2);

    for (size_t i = 0; i < P.size(); i++)
    {
        Point p = P[i];
        Sign status = point_is_in_half_plane(p, q1, q2);
        if (status != prev_status && status != ZERO && prev_status != ZERO)
        {
            Point intersect;
            if (intersect_segments(prev_p, p, q1, q2, intersect))
                result.push_back(intersect);
        }

        switch (status)
        {
        case NEGATIVE:
            break;
        case ZERO:
            break;
        case POSITIVE:
            result.push_back(p);
            break;
        }
        prev_p = p;
        prev_status = status;
    }
}

std::vector<Point> convex_clip_polygon(const std::vector<Point> &P, const std::vector<Point> &clip)
{
    std::vector<Point> result;
    std::vector<Point> tmp1 = P;
    std::vector<Point> tmp2;
    std::vector<Point> *src = &tmp1;
    std::vector<Point> *dst = &tmp2;
    for (unsigned int i = 0; i < clip.size(); i++)
    {
        unsigned int j = (i + 1) % clip.size();
        const Point &p1 = clip[i];
        const Point &p2 = clip[j];
        clip_polygon_by_half_plane(*src, p1, p2, *dst);
        std::swap(src, dst);
    }
    result = *src;
    return result;
}
