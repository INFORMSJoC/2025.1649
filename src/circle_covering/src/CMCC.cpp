#include "CMCC.h"

// MCC
inline double sqr(double r)
{
    return r * r;
}

// Class Definition
MCC::MCC(std::vector<Point> subject) : mSubject(subject)
{
    mRadius = 0;
    support_end = mSupportPoints.begin();
    fsize = 0;
    ssize = 0;

    // set initial center
    mCenter = Point(0.0, 0.0);

    // compute MCC
    computer_mcc();

    center = mCenter;
    radius = mRadius;
}

void MCC::mtf_mb(std::_List_iterator<Point> n)
{
    // Algorithm 1: mtf_mb (L_{n-1}, B), where L_{n-1} = [L.begin, n)
    // B: the set of forced points, defining the current ball
    // S: the superset of support points computed by the algorithm
    // --------------------------------------------------------------
    // from B. Gaertner, Fast and Robust Smallest Enclosing Balls, ESA 1999,
    // http://www.inf.ethz.ch/personal/gaertner/texts/own_work/esa99_final.pdf

    //   PRE: B = S
    assert(fsize == ssize);

    support_end = mSupportPoints.begin();
    if ((fsize) == 3)
        return;

    // incremental construction
    for (std::_List_iterator<Point> i = mSupportPoints.begin(); i != n;)
    {
        // INV: (support_end - L.begin() == |S|-|B|)
        assert(std::distance(mSupportPoints.begin(), support_end) == ssize - fsize);

        std::_List_iterator<Point> j = i++;
        if (excess(*j) > 0)
            if (push(*j))
            {              // B := B + p_i
                mtf_mb(j); // mtf_mb (L_{i-1}, B + p_i)
                pop();     // B := B - p_i
                mtf_move_to_front(j);
            }
    }
    // POST: the range [L.begin(), support_end) stores the set S\B
}

void MCC::mtf_move_to_front(std::_List_iterator<Point> j)
{
    if (support_end == j)
        support_end++;
    mSupportPoints.splice(mSupportPoints.begin(), mSupportPoints, j);
}

void MCC::computer_mcc()
{
    // Algorithm 2: pivot_mb (L_{n-1}), where L_{n-1} = [L.begin, n)
    // --------------------------------------------------------------
    // from B. Gaertner, Fast and Robust Smallest Enclosing Balls, ESA 1999,
    // http://www.inf.ethz.ch/personal/gaertner/texts/own_work/esa99_final.pdf
    Point current_p, current_c;
    double e, max_e, current_sqr_r;
    unsigned int loops = 0;
    double best_sqr_r = mRadius * mRadius;
    do
    {
        current_sqr_r = mRadius * mRadius;
        current_p = mSubject[0];
        max_e = 0.0;
        for (size_t i = 0; i < mSubject.size(); i++)
        {
            current_c = mCenter;
            e = sqr(mSubject[i].getDistance(current_c)) - current_sqr_r;
            if (e > max_e)
            {
                max_e = e;
                current_p = mSubject[i];
            }
        }

        if (current_sqr_r < -1e-7 || max_e > 1e-7)
        {
            // check if the pivot is already contained in the support set
            if (std::find(mSupportPoints.begin(), support_end, current_p) == support_end)
            {
                assert(fsize == 0);
                if (push(current_p))
                {
                    mtf_mb(support_end);
                    pop();
                    pivot_move_to_front(current_p);
                }
            }
        }
        if (best_sqr_r < mRadius * mRadius)
        {
            best_sqr_r = mRadius * mRadius;
            loops = 0;
        }
        else
            ++loops;
    } while (loops < 2);
}

void MCC::pivot_move_to_front(Point p)
{
    mSupportPoints.push_front(p);
    if (std::distance(mSupportPoints.begin(), support_end) == 4)
        support_end--;
}

inline double MCC::excess(Point p) const
{
    double e = sqr(p.getDistance(mCenter)) - mRadius * mRadius;
    return e;
}

void MCC::pop()
{
    --fsize;
}

bool MCC::push(Point p)
{
    int i, j;
    double eps = 1e-10;

    if (fsize == 0)
    {
        q0 = p;
        c[0] = q0;
        sqr_r[0] = 0;
    }
    else
    {
        // set v_fsize to Q_fsize
        v[fsize] = p - q0;

        // compute the a_{fsize,i}, i< fsize
        for (i = 1; i < fsize; ++i)
        {
            a[fsize][i] = v[i].x * v[fsize].x + v[i].y * v[fsize].y;
            a[fsize][i] *= (2 / z[i]);
        }

        // update v_fsize to Q_fsize-\bar{Q}_fsize
        for (i = 1; i < fsize; ++i)
        {
            v[fsize].x -= a[fsize][i] * v[i].x;
            v[fsize].y -= a[fsize][i] * v[i].y;
        }

        // compute z_fsize
        z[fsize] = v[fsize].dot(v[fsize]);
        z[fsize] *= 2;

        // reject push if z_fsize too small
        if (z[fsize] < sqr(eps * mRadius))
        {
            return false;
        }

        // update c, sqr_r
        double e = sqr(p.getDistance(c[fsize - 1])) - sqr_r[fsize - 1];
        f[fsize] = e / z[fsize];

        c[fsize].x = c[fsize - 1].x + f[fsize] * v[fsize].x;
        c[fsize].y = c[fsize - 1].y + f[fsize] * v[fsize].y;
        sqr_r[fsize] = sqr_r[fsize - 1] + e * f[fsize] / 2;
    }
    mCenter = c[fsize];
    mRadius = std::sqrt(sqr_r[fsize]);
    ssize = ++fsize;
    return true;
}

// CMCC
// Class Definition
CMCC::CMCC(std::vector<Point> subject, std::vector<Point> border) : mSubject(subject), mBorder(border)
{
    MCC mb(mSubject);
    mSubject = convexhull(subject);
    std::pair<Point, double> result = find_cmcc(mSubject, mBorder, mb.center, mb.radius);
    center = result.first;
    radius = result.second;
}

CMCC::CMCC(std::vector<Point> subject, std::vector<Point> border, std::vector<std::vector<Point>> holes) : mSubject(subject), mBorder(border), mHoles(holes)
{
    MCC mb(mSubject);
    int in_hole_number;
    bool in_hole = false;
    for (size_t i = 0; i < mHoles.size(); i++)
    {
        if (mb.center.is_in_polygon(mHoles[i]))
        {
            in_hole = true;
            in_hole_number = i;
            break;
        }
    }

    if (mb.center.is_in_polygon(mBorder) && !in_hole)
    {
        center = mb.center;
        radius = mb.radius;
    }
    else if (!mb.center.is_in_polygon(mBorder))
    {
        mSubject = convexhull(subject);
        std::pair<Point, double> result = find_cmcc(mSubject, mBorder, mb.center, mb.radius);
        center = result.first;
        radius = result.second;
    }
    else if (in_hole)
    {
        mSubject = convexhull(subject);
        std::pair<Point, double> result = find_cmcc(mSubject, mHoles[in_hole_number], mb.center, mb.radius);
        center = result.first;
        radius = result.second;
    }
}

// find convex hull
std::vector<Point> CMCC::convexhull(std::vector<Point> points)
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

// find cmcc
Point CMCC::point_to_line_segment(Point a, Point b, Point point) const
{
    Point result;
    // double distance;
    Point lineseg = b - a;
    Point a_point = point - a;
    double point_dot_lineseg = a_point.dot(lineseg);
    double lineseg_dot_lineseg = lineseg.dot(lineseg);
    double l = point_dot_lineseg / lineseg_dot_lineseg;
    if (l <= 0)
    {
        result.x = a.x;
        result.y = a.y;
    }
    else if (l >= 1)
    {
        result.x = b.x;
        result.y = b.y;
    }
    else
    {
        result.x = a.x + l * lineseg.x;
        result.y = a.y + l * lineseg.y;
    }
    // distance = result.getDistance(point);

    return result;
}

Point CMCC::bisect_line_segment(Point a, Point b, Point c, Point d) const
{
    Point result;
    const double toler = 1e-10;
    if (a.x - b.x > toler)
        std::swap(a, b);
    if (c.x - d.x > toler)
        std::swap(c, d);

    Point vertex_pair = b - a;
    Point line_segment = d - c;
    double distance;

    double f1_x = a.x - b.x;
    double f1_y = a.y - b.y;
    double const1 = (a.x - b.x) * (a.x + b.x) / 2 + (a.y - b.y) * (a.y + b.y) / 2;

    double f2_x = d.y - c.y;
    double f2_y = c.x - d.x;
    double const2 = c.x * d.y - d.x * c.y;

    double det = f1_x * f2_y - f2_x * f1_y;

    if (fabs(det) <= toler)
        return Point{-1, -1};

    double inv_det = 1.0 / det;
    result.x = (f2_y * const1 - f1_y * const2) * inv_det;
    result.y = (f1_x * const2 - f2_x * const1) * inv_det;
    return result;
}

double CMCC::find_y_on_line(Point a, Point b, double x) const
{
    if (fabs(x - a.x) > fabs(x - b.x))
        return ((x - a.x) * (b.y - a.y) / (b.x - a.x) + a.y);
    else
        return ((x - b.x) * (a.y - b.y) / (a.x - b.x) + b.y);
}

double CMCC::find_x_on_line(Point a, Point b, double y) const
{
    if (fabs(y - a.y) > fabs(y - b.y))
        return ((y - a.y) * (b.x - a.x) / (b.y - a.y) + a.x);
    else
        return ((y - b.y) * (a.x - b.x) / (a.y - b.y) + b.x);
}

double CMCC::dis_point_to_line(Point a, Point b, Point point) const
{
    double distance = fabs((b.x - a.x) * (a.y - point.y) - (a.x - point.x) * (b.y - a.y)) / ((b - a).getNorm());
    return distance;
}

double CMCC::dis_point_to_line_segment(const Point a, const Point b, const Point point)
{
    double result;
    Point lineseg = b - a;
    Point a_point = point - a;
    double point_dot_lineseg = a_point.dot(lineseg);
    double lineseg_dot_lineseg = lineseg.dot(lineseg);
    double l = point_dot_lineseg / lineseg_dot_lineseg;
    if (l <= 0)
        result = point.getDistance(a);
    else if (l >= 1)
        result = point.getDistance(b);
    else
        result = fabs((b.x - a.x) * (a.y - point.y) - (a.x - point.x) * (b.y - a.y)) / ((b - a).getNorm());
    return result;
}

Point CMCC::coord_point_to_line(Point a, Point b, Point point) const
{
    Point result;
    Point lineseg = b - a;
    Point a_point = point - a;
    double point_dot_lineseg = a_point.dot(lineseg);
    double lineseg_dot_lineseg = lineseg.dot(lineseg);
    double l = point_dot_lineseg / lineseg_dot_lineseg;
    result.x = a.x + l * lineseg.x;
    result.y = a.y + l * lineseg.y;
    return result;
}

Point CMCC::coord_point_to_line_segment(Point a, Point b, Point point) const
{
    Point result;
    Point lineseg = b - a;
    Point a_point = point - a;
    double point_dot_lineseg = a_point.dot(lineseg);
    double lineseg_dot_lineseg = lineseg.dot(lineseg);
    double l = point_dot_lineseg / lineseg_dot_lineseg;
    if (l <= 0)
    {
        result.x = a.x;
        result.y = a.y;
    }
    else if (l >= 1)
    {
        result.x = b.x;
        result.y = b.y;
    }
    else
    {
        result.x = a.x + l * lineseg.x;
        result.y = a.y + l * lineseg.y;
    }
    return result;
}

std::pair<Point, double> CMCC::Megiddo(Point a, Point b, std::vector<Point> subject)
{
    double toler = 1e-10;
    std::vector<Point> original_subject = subject; // Save a copy of the original before elimination
    int the_end_indx = -1;                         // initialize
    Point p0;
    Point p1;
    Point p;
    if (fabs(a.x - b.x) > fabs(a.y - b.y))
    {
        while (subject.size() > 0)
        {
            std::vector<int> indx_to_discard;
            int n_pairs = subject.size() / 2;
            int odd = subject.size() % 2;
            std::vector<Point> crits;
            std::vector<std::array<int, 2>> crits_pivots_indx; // the point whose projx is smaller goes first
            for (int i = 0; i < n_pairs; i++)
            {
                p0 = coord_point_to_line(a, b, subject[2 * i]);
                p1 = coord_point_to_line(a, b, subject[2 * i + 1]);
                if (fabs(p0.x - p1.x) < toler)
                {
                    double d0 = dis_point_to_line(a, b, subject[2 * i]);
                    double d1 = dis_point_to_line(a, b, subject[2 * i + 1]);
                    if (d0 < d1)
                        indx_to_discard.push_back(2 * i);
                    else
                        indx_to_discard.push_back(2 * i + 1);
                }
                else
                {
                    Point this_crit = bisect_line_segment(subject[2 * i], subject[2 * i + 1], a, b);
                    crits.push_back(this_crit);
                    if (p0.x < p1.x)
                        crits_pivots_indx.push_back({2 * i, 2 * i + 1});
                    else
                        crits_pivots_indx.push_back({2 * i + 1, 2 * i});
                }
            }
            if (odd == 1)
            {
                if (n_pairs > 0)
                {
                    p0 = coord_point_to_line(a, b, subject[2 * n_pairs]);
                    p1 = coord_point_to_line(a, b, subject[0]);
                    if (fabs(p0.x - p1.x) < toler)
                    {
                        double d0 = dis_point_to_line(a, b, subject[2 * n_pairs]);
                        double d1 = dis_point_to_line(a, b, subject[0]);
                        if (d0 < d1)
                            indx_to_discard.push_back(2 * n_pairs);
                        else
                            indx_to_discard.push_back(0);
                    }
                    else
                    {
                        Point this_crit = bisect_line_segment(subject[2 * n_pairs], subject[0], a, b);
                        crits.push_back(this_crit);
                        if (p0.x < p1.x)
                            crits_pivots_indx.push_back({2 * n_pairs, 0});
                        else
                            crits_pivots_indx.push_back({0, 2 * n_pairs});
                    }
                }
                else
                {
                    p = coord_point_to_line(a, b, subject[0]);
                    p = coord_point_to_line_segment(a, b, p);
                    double cover_radius = 0.0;
                    for (size_t i = 0; i < original_subject.size(); i++)
                    {
                        double this_dist = p.getDistance(original_subject[i]);
                        if (this_dist > cover_radius)
                            cover_radius = this_dist;
                    }
                    return std::make_pair(p, cover_radius);
                }
            }
            if (!crits.empty())
            {
                double xmedian;
                double ymedian;
                if (crits.size() % 2 == 0)
                {
                    int xmin_indx = 0;
                    double temp_x = crits[0].x;
                    for (size_t i = 0; i < crits.size(); i++)
                    {
                        if (crits[i].x < temp_x)
                        {
                            temp_x = crits[i].x;
                            xmin_indx = i;
                        }
                    }
                    std::vector<double> remaining_crits;
                    for (size_t i = 0; i < crits.size(); i++)
                    {
                        if (i != xmin_indx)
                            remaining_crits.push_back(crits[i].x);
                    }
                    std::sort(remaining_crits.begin(), remaining_crits.end());
                    int median_index = remaining_crits.size() / 2;
                    xmedian = remaining_crits[median_index];
                }
                else
                {
                    std::vector<double> crits_x;
                    for (size_t i = 0; i < crits.size(); i++)
                        crits_x.push_back(crits[i].x);
                    std::sort(crits_x.begin(), crits_x.end());
                    int median_index = crits.size() / 2;
                    xmedian = crits_x[median_index];
                }
                ymedian = find_y_on_line(a, b, xmedian);
                // compute g(xm) and find I
                double gxm = 0.0;
                std::vector<int> Iset;
                Point median_point = {xmedian, ymedian};
                for (size_t i = 0; i < subject.size(); i++)
                {
                    double this_one = median_point.getDistance(subject[i]);
                    if (this_one - gxm > -toler)
                    {
                        if (this_one - gxm > toler)
                            Iset.clear();
                        gxm = this_one;
                        Iset.push_back(i);
                    }
                }
                // decide
                std::vector<Point> proj_Iset;
                for (size_t i = 0; i < Iset.size(); i++)
                    proj_Iset.push_back(coord_point_to_line(a, b, subject[Iset[i]]));
                bool condition1 = true;
                bool condition2 = true;
                for (size_t i = 0; i < Iset.size(); i++)
                {
                    if (xmedian >= proj_Iset[i].x - toler)
                        condition1 = false;
                    if (xmedian <= proj_Iset[i].x + toler)
                        condition2 = false;
                }
                if (condition1)
                {
                    // xm < x*, then for each pair whose crit is < xm, discard the point whose projx is greater
                    if (xmedian > std::max(a.x, b.x))
                    {
                        the_end_indx = a.x > b.x ? 0 : 1;
                        break;
                    }

                    for (size_t i = 0; i < crits.size(); i++)
                    {
                        if (crits[i].x <= xmedian)
                            indx_to_discard.push_back(crits_pivots_indx[i][1]);
                    }
                }
                else if (condition2)
                {
                    // xm > x*, then for each pair whose crit is > xm, discard the point whose projx is smaller
                    if (xmedian < std::min(a.x, b.x))
                    {
                        the_end_indx = a.x < b.x ? 0 : 1;
                        break;
                    }
                    for (size_t i = 0; i < crits.size(); i++)
                    {
                        if (crits[i].x >= xmedian)
                            indx_to_discard.push_back(crits_pivots_indx[i][0]);
                    }
                }
                else
                {
                    p = coord_point_to_line_segment(a, b, Point{xmedian, ymedian});
                    // Find the cover radius: must check for all original pivots
                    double cover_radius = 0.0;
                    for (size_t i = 0; i < original_subject.size(); i++)
                    {
                        double this_dist = p.getDistance(original_subject[i]);
                        if (this_dist > cover_radius)
                            cover_radius = this_dist;
                    }
                    return std::make_pair(p, cover_radius);
                }
            }
            // Update pivots
            if (indx_to_discard.size() > 0)
            {
                std::sort(indx_to_discard.begin(), indx_to_discard.end(), std::greater<int>());
                for (int index : indx_to_discard)
                    subject.erase(subject.begin() + index);
            }
            else
            {
                if (subject.size() < 3)
                {
                    double d0 = dis_point_to_line(a, b, subject[0]);
                    double d1 = dis_point_to_line(a, b, subject[1]);
                    p = d0 < d1 ? p1 : p0;
                    p = coord_point_to_line_segment(a, b, p);
                    double cover_radius = 0.0;
                    for (size_t i = 0; i < original_subject.size(); i++)
                    {
                        double this_dist = p.getDistance(original_subject[i]);
                        if (this_dist > cover_radius)
                            cover_radius = this_dist;
                    }
                    return std::make_pair(p, cover_radius);
                }
            }
        }
    }
    else
    {
        while (subject.size() > 0)
        {
            std::vector<int> indx_to_discard;
            int n_pairs = subject.size() / 2;
            int odd = subject.size() % 2;
            std::vector<Point> crits;
            std::vector<std::array<int, 2>> crits_pivots_indx; // the point whose projx is smaller goes first
            for (int i = 0; i < n_pairs; i++)
            {
                p0 = coord_point_to_line(a, b, subject[2 * i]);
                p1 = coord_point_to_line(a, b, subject[2 * i + 1]);
                if (fabs(p0.y - p1.y) < toler)
                {
                    double d0 = dis_point_to_line(a, b, subject[2 * i]);
                    double d1 = dis_point_to_line(a, b, subject[2 * i + 1]);
                    if (d0 < d1)
                        indx_to_discard.push_back(2 * i);
                    else
                        indx_to_discard.push_back(2 * i + 1);
                }
                else
                {
                    Point this_crit = bisect_line_segment(subject[2 * i], subject[2 * i + 1], a, b);
                    crits.push_back(this_crit);
                    if (p0.y < p1.y)
                        crits_pivots_indx.push_back({2 * i, 2 * i + 1});
                    else
                        crits_pivots_indx.push_back({2 * i + 1, 2 * i});
                }
            }
            if (odd == 1)
            {
                if (n_pairs > 0)
                {
                    p0 = coord_point_to_line(a, b, subject[2 * n_pairs]);
                    p1 = coord_point_to_line(a, b, subject[0]);
                    if (fabs(p0.y - p1.y) < toler)
                    {
                        double d0 = dis_point_to_line(a, b, subject[2 * n_pairs]);
                        double d1 = dis_point_to_line(a, b, subject[0]);
                        if (d0 < d1)
                            indx_to_discard.push_back(2 * n_pairs);
                        else
                            indx_to_discard.push_back(0);
                    }
                    else
                    {
                        Point this_crit = bisect_line_segment(subject[2 * n_pairs], subject[0], a, b);
                        crits.push_back(this_crit);
                        if (p0.y < p1.y)
                            crits_pivots_indx.push_back({2 * n_pairs, 0});
                        else
                            crits_pivots_indx.push_back({0, 2 * n_pairs});
                    }
                }
                else
                {
                    p = coord_point_to_line(a, b, subject[0]);
                    p = coord_point_to_line_segment(a, b, p);
                    double cover_radius = 0.0;
                    for (size_t i = 0; i < original_subject.size(); i++)
                    {
                        double this_dist = p.getDistance(original_subject[i]);
                        if (this_dist > cover_radius)
                            cover_radius = this_dist;
                    }
                    return std::make_pair(p, cover_radius);
                }
            }

            if (!crits.empty())
            {
                double xmedian;
                double ymedian;
                if (crits.size() % 2 == 0)
                {
                    int ymin_indx = 0;
                    double temp_y = crits[0].y;
                    for (size_t i = 0; i < crits.size(); i++)
                    {
                        if (crits[i].y < temp_y)
                        {
                            temp_y = crits[i].y;
                            ymin_indx = i;
                        }
                    }
                    std::vector<double> remaining_crits;
                    for (size_t i = 0; i < crits.size(); i++)
                    {
                        if (i != ymin_indx)
                            remaining_crits.push_back(crits[i].y);
                    }
                    std::sort(remaining_crits.begin(), remaining_crits.end());
                    int median_index = remaining_crits.size() / 2;
                    ymedian = remaining_crits[median_index];
                }
                else
                {
                    std::vector<double> crits_y;
                    for (size_t i = 0; i < crits.size(); i++)
                        crits_y.push_back(crits[i].y);
                    std::sort(crits_y.begin(), crits_y.end());
                    int median_index = crits.size() / 2;
                    ymedian = crits_y[median_index];
                }
                xmedian = find_x_on_line(a, b, ymedian);
                // compute g(xm) and find I
                double gxm = 0;
                std::vector<int> Iset;
                Point median_point = {xmedian, ymedian};
                for (size_t i = 0; i < subject.size(); i++)
                {
                    double this_one = median_point.getDistance(subject[i]);
                    if (this_one - gxm > -toler)
                    {
                        if (this_one - gxm > toler)
                            Iset.clear();
                        gxm = this_one;
                        Iset.push_back(i);
                    }
                }
                // decide
                std::vector<Point> proj_Iset;
                for (size_t i = 0; i < Iset.size(); i++)
                    proj_Iset.push_back(coord_point_to_line(a, b, subject[Iset[i]]));
                bool condition1 = true;
                bool condition2 = true;
                for (size_t i = 0; i < Iset.size(); i++)
                {
                    if (ymedian >= proj_Iset[i].y - toler)
                        condition1 = false;
                    if (ymedian <= proj_Iset[i].y + toler)
                        condition2 = false;
                }
                if (condition1)
                {
                    // ym < y*, then for each pair whose crit is < ym, discard the point whose projy is greater
                    if (ymedian > std::max(a.y, b.y))
                    {
                        the_end_indx = a.y > b.y ? 0 : 1;
                        break;
                    }
                    for (size_t i = 0; i < crits.size(); i++)
                    {
                        if (crits[i].y <= ymedian)
                            indx_to_discard.push_back(crits_pivots_indx[i][1]);
                    }
                }
                else if (condition2)
                {
                    // xm > x*, then for each pair whose crit is > xm, discard the point whose projx is smaller
                    if (ymedian < std::min(a.y, b.y))
                    {
                        the_end_indx = a.y < b.y ? 0 : 1;
                        break;
                    }
                    for (size_t i = 0; i < crits.size(); i++)
                    {
                        if (crits[i].y >= ymedian)
                            indx_to_discard.push_back(crits_pivots_indx[i][0]);
                    }
                }
                else
                {
                    p = coord_point_to_line_segment(a, b, Point{xmedian, ymedian});
                    // Find the cover radius: must check for all original pivots
                    double cover_radius = 0.0;
                    for (size_t i = 0; i < original_subject.size(); i++)
                    {
                        double this_dist = p.getDistance(original_subject[i]);
                        if (this_dist > cover_radius)
                            cover_radius = this_dist;
                    }
                    return std::make_pair(p, cover_radius);
                }
            }
            // Update pivots
            if (indx_to_discard.size() > 0)
            {
                std::sort(indx_to_discard.begin(), indx_to_discard.end(), std::greater<int>());
                for (int index : indx_to_discard)
                    subject.erase(subject.begin() + index);
            }
            else
            {
                if (subject.size() < 3)
                {
                    double d0 = dis_point_to_line(a, b, subject[0]);
                    double d1 = dis_point_to_line(a, b, subject[1]);
                    p = d0 < d1 ? p1 : p0;
                    p = coord_point_to_line_segment(a, b, p);
                    double cover_radius = 0.0;
                    for (size_t i = 0; i < original_subject.size(); i++)
                    {
                        double this_dist = p.getDistance(original_subject[i]);
                        if (this_dist > cover_radius)
                            cover_radius = this_dist;
                    }
                    return std::make_pair(p, cover_radius);
                }
            }
        }
    }
    if (the_end_indx >= 0)
    {
        // return the early termination solution
        // Find the cover radius: must check for all original pivots
        p = the_end_indx == 0 ? a : b;
        double cover_radius = 0.0;
        for (size_t i = 0; i < original_subject.size(); i++)
        {
            double this_dist = p.getDistance(original_subject[i]);
            if (this_dist > cover_radius)
                cover_radius = this_dist;
        }
        return std::make_pair(p, cover_radius);
    }
    else
        return std::make_pair(Point(0, 0), -1);
}

std::pair<Point, double> CMCC::find_cmcc(std::vector<Point> subject, std::vector<Point> border, Point MCCcenter, double MCCradius)
{
    double toler = 1e-10; // toler must be a small number strictly greater than 0
    // only take the pivots that lie on the convex hull of all the pivots
    std::vector<std::pair<int, double>> d;
    for (size_t i = 0; i < border.size(); i++)
    {
        int j = (i + 1) % border.size();
        Point projected_point = coord_point_to_line_segment(border[i], border[j], MCCcenter);
        double projected_dist = projected_point.getDistance(MCCcenter);
        if (projected_dist - MCCradius >= toler)
            continue; // n_edges_skipped++;
        d.push_back({i, projected_dist});
    }
    std::sort(d.begin(), d.end(), [](const std::pair<int, double> &a, const std::pair<int, double> &b)
              { return a.second < b.second; });

    // Iterate the sorted list
    Point best_cover_center;
    double best_cover_radius = 100000;
    int di = 0;
    Point found_center;
    double found_radius;
    while (di < d.size())
    {
        if (d[di].second > sqrt(best_cover_radius * best_cover_radius - MCCradius * MCCradius))
        {
            di++;
            continue;
        }
        int this_segment_index = d[di].first; // the key in the (key, value) pair of the di-th entry of sorted_list_d
        int this_segment_index_ = (this_segment_index + 1) % border.size();
        std::vector<Point> the_two_vertices = {border[this_segment_index], border[this_segment_index_]};
        // First try method 1
        int found = 0;
        for (size_t ei = 0; ei < 2; ei++)
        {
            // #for the two vertices of the segment
            double max_dist_to_segment_end = 0;
            int pivot_indx = 0;
            for (size_t i = 0; i < subject.size(); i++)
            {
                double this_dist_to_segment_end = the_two_vertices[ei].getDistance(subject[i]);
                if (this_dist_to_segment_end - max_dist_to_segment_end >= toler)
                {
                    max_dist_to_segment_end = this_dist_to_segment_end;
                    pivot_indx = i;
                }
            }
            Point pivot_project_point = coord_point_to_line_segment(the_two_vertices[0], the_two_vertices[1], subject[pivot_indx]);
            if ((pivot_project_point - the_two_vertices[ei]).getNorm() <= toler)
            {
                found_center = the_two_vertices[ei];
                found_radius = max_dist_to_segment_end;
                found = 1;
                break;
            }
        }
        if (found)
        {
            // if found_radius < best_cover_radius:
            if (found_radius <= best_cover_radius + toler)
            {
                best_cover_radius = found_radius;
                best_cover_center = found_center;
            }
            di++;
            continue;
        }
        // Use Megiddo's algorithm
        std::pair<Point, double> result = Megiddo(the_two_vertices[0], the_two_vertices[1], subject);
        found_center = result.first;
        found_radius = result.second;
        // if found_radius < best_cover_radius:
        if (found_radius <= best_cover_radius + toler)
        {
            best_cover_radius = found_radius;
            best_cover_center = found_center;
        }
        di++;
    }
    return std::make_pair(best_cover_center, best_cover_radius);
}