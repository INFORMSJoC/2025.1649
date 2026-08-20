#pragma once

#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <climits>
#include <numeric>

#include "Point.h"

namespace Clipper2Lib
{

#if (defined(__cpp_exceptions) && __cpp_exceptions) || (defined(__EXCEPTIONS) && __EXCEPTIONS)

  class Clipper2Exception : public std::exception
  {
  public:
    explicit Clipper2Exception(const char *description) : m_descr(description) {}
    virtual const char *what() const throw() override { return m_descr.c_str(); }

  private:
    std::string m_descr;
  };

  static const char *precision_error =
      "Precision exceeds the permitted range";
  static const char *range_error =
      "Values exceed permitted range";
  static const char *scale_error =
      "Invalid scale (either 0 or too large)";
  static const char *non_pair_error =
      "There must be 2 values for each coordinate";
#endif

  // error codes (2^n)
  const int precision_error_i = 1; // non-fatal
  const int scale_error_i = 2;     // non-fatal
  const int non_pair_error_i = 4;  // non-fatal
  const int range_error_i = 64;

#ifndef PI
  static const double PI = 3.141592653589793238;
#endif
  static const int MAX_DECIMAL_PRECISION = 15; // see https://github.com/AngusJohnson/Clipper2/discussions/564
  static const int64_t MAX_COORD = INT64_MAX;
  static const int64_t MIN_COORD = -MAX_COORD;
  static const int64_t INVALID = INT64_MAX;
  const double max_coord = static_cast<double>(MAX_COORD);
  const double min_coord = static_cast<double>(MIN_COORD);

  static const double MAX_DBL = (std::numeric_limits<double>::max)();

  static void DoError(int error_code)
  {
#if (defined(__cpp_exceptions) && __cpp_exceptions) || (defined(__EXCEPTIONS) && __EXCEPTIONS)
    switch (error_code)
    {
    case precision_error_i:
      throw Clipper2Exception(precision_error);
    case scale_error_i:
      throw Clipper2Exception(scale_error);
    case non_pair_error_i:
      throw Clipper2Exception(non_pair_error);
    case range_error_i:
      throw Clipper2Exception(range_error);
    }
#else
    ++error_code; // only to stop compiler warning
#endif
  }

  // By far the most widely used filling rules for polygons are EvenOdd
  // and NonZero, sometimes called Alternate and Winding respectively.
  // https://en.wikipedia.org/wiki/Nonzero-rule
  enum class FillRule
  {
    EvenOdd,
    NonZero,
    Positive,
    Negative
  };

  // Point ------------------------------------------------------------------------
  // inclue Point.h

  // nb: using 'using' here (instead of typedef) as they can be used in templates

  using Path = std::vector<Point>;
  using Paths = std::vector<Path>;

  // Rect ------------------------------------------------------------------------
  struct Rect
  {
    double left;
    double top;
    double right;
    double bottom;

    Rect() : left(0), top(0), right(0), bottom(0)
    {
    }

    Rect(double l, double t, double r, double b) : left(l), top(t), right(r), bottom(b)
    {
    }

    Rect(bool is_valid)
    {
      if (is_valid)
      {
        left = right = top = bottom = 0;
      }
      else
      {
        left = top = (std::numeric_limits<double>::max)();
        right = bottom = -(std::numeric_limits<double>::max)();
      }
    }

    double Width() const { return right - left; }
    double Height() const { return bottom - top; }
    void Width(double width) { right = left + width; }
    void Height(double height) { bottom = top + height; }

    Point MidPoint() const
    {
      return Point((left + right) / 2, (top + bottom) / 2);
    }

    Path AsPath() const
    {
      Path result;
      result.reserve(4);
      result.push_back(Point(left, top));
      result.push_back(Point(right, top));
      result.push_back(Point(right, bottom));
      result.push_back(Point(left, bottom));
      return result;
    }

    bool Contains(const Point &pt) const
    {
      return pt.x > left && pt.x < right && pt.y > top && pt.y < bottom;
    }

    bool Contains(const Rect &rec) const
    {
      return rec.left >= left && rec.right <= right &&
             rec.top >= top && rec.bottom <= bottom;
    }

    void Scale(double scale)
    {
      left *= scale;
      top *= scale;
      right *= scale;
      bottom *= scale;
    }

    bool IsEmpty() const { return bottom <= top || right <= left; };

    bool Intersects(const Rect &rec) const
    {
      return ((std::max)(left, rec.left) <= (std::min)(right, rec.right)) &&
             ((std::max)(top, rec.top) <= (std::min)(bottom, rec.bottom));
    };

    bool operator==(const Rect &other) const
    {
      return left == other.left && right == other.right &&
             top == other.top && bottom == other.bottom;
    }
  };

  inline Rect ScaleRect(const Rect &rect, double scale)
  {
    Rect result;
    result.left = rect.left * scale;
    result.top = rect.top * scale;
    result.right = rect.right * scale;
    result.bottom = rect.bottom * scale;
    return result;
  }

  static const Rect MaxInvalidRect = Rect(MAX_DBL, MAX_DBL, -MAX_DBL, -MAX_DBL);

  static Rect GetBounds(const Path &path)
  {
    auto xmin = (std::numeric_limits<double>::max)();
    auto ymin = (std::numeric_limits<double>::max)();
    auto xmax = std::numeric_limits<double>::lowest();
    auto ymax = std::numeric_limits<double>::lowest();
    for (const auto &p : path)
    {
      if (p.x < xmin)
        xmin = p.x;
      if (p.x > xmax)
        xmax = p.x;
      if (p.y < ymin)
        ymin = p.y;
      if (p.y > ymax)
        ymax = p.y;
    }
    return Rect(xmin, ymin, xmax, ymax);
  }

  static Rect GetBounds(const Paths &paths)
  {
    auto xmin = (std::numeric_limits<double>::max)();
    auto ymin = (std::numeric_limits<double>::max)();
    auto xmax = std::numeric_limits<double>::lowest();
    auto ymax = std::numeric_limits<double>::lowest();
    for (const Path &path : paths)
      for (const Point &p : path)
      {
        if (p.x < xmin)
          xmin = p.x;
        if (p.x > xmax)
          xmax = p.x;
        if (p.y < ymin)
          ymin = p.y;
        if (p.y > ymax)
          ymax = p.y;
      }
    return Rect(xmin, ymin, xmax, ymax);
  }

  inline Path ScalePath(const Path &path, double scale_x, double scale_y, int &error_code)
  {
    Path result;
    if (scale_x == 0 || scale_y == 0)
    {
      error_code |= scale_error_i;
      DoError(scale_error_i);
      // if no exception, treat as non-fatal error
      if (scale_x == 0)
        scale_x = 1.0;
      if (scale_y == 0)
        scale_y = 1.0;
    }

    result.reserve(path.size());
    std::transform(path.begin(), path.end(), back_inserter(result),
                   [scale_x, scale_y](const auto &pt)
                   { return Point(round(pt.x * scale_x), round(pt.y * scale_y)); });
    return result;
  }

  inline Path ScalePath(const Path &path, double scale, int &error_code)
  {
    return ScalePath(path, scale, scale, error_code);
  }

  inline Paths ScalePaths(const Paths &paths, double scale_x, double scale_y, int &error_code)
  {
    Paths result;

    Rect r = GetBounds(paths);
    if ((r.left * scale_x) < min_coord || (r.right * scale_x) > max_coord || (r.top * scale_y) < min_coord || (r.bottom * scale_y) > max_coord)
    {
      error_code |= range_error_i;
      DoError(range_error_i);
      return result; // empty path
    }
    result.reserve(paths.size());
    std::transform(paths.begin(), paths.end(), back_inserter(result),
                   [=, &error_code](const auto &path)
                   { return ScalePath(path, scale_x, scale_y, error_code); });
    return result;
  }

  inline Paths ScalePaths(const Paths &paths, double scale, int &error_code)
  {
    return ScalePaths(paths, scale, scale, error_code);
  }

  inline double Sqr(double val)
  {
    return val * val;
  }

  inline bool NearEqual(const Point &p1, const Point &p2, double max_dist_sqrd)
  {
    return Sqr(p1.x - p2.x) + Sqr(p1.y - p2.y) < max_dist_sqrd;
  }

  inline Path StripNearEqual(const Path &path, double max_dist_sqrd, bool is_closed_path)
  {
    if (path.size() == 0)
      return Path();
    Path result;
    result.reserve(path.size());
    typename Path::const_iterator path_iter = path.cbegin();
    Point first_pt = *path_iter++, last_pt = first_pt;
    result.push_back(first_pt);
    for (; path_iter != path.cend(); ++path_iter)
    {
      if (!NearEqual(*path_iter, last_pt, max_dist_sqrd))
      {
        last_pt = *path_iter;
        result.push_back(last_pt);
      }
    }
    if (!is_closed_path)
      return result;
    while (result.size() > 1 &&
           NearEqual(result.back(), first_pt, max_dist_sqrd))
      result.pop_back();
    return result;
  }

  inline Paths StripNearEqual(const Paths &paths, double max_dist_sqrd, bool is_closed_path)
  {
    Paths result;
    result.reserve(paths.size());
    for (typename Paths::const_iterator paths_citer = paths.cbegin();
         paths_citer != paths.cend(); ++paths_citer)
    {
      result.push_back(StripNearEqual(*paths_citer, max_dist_sqrd, is_closed_path));
    }
    return result;
  }

  inline void StripDuplicates(Path &path, bool is_closed_path)
  {
    // https://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector#:~:text=Let%27s%20compare%20three%20approaches%3A
    path.erase(std::unique(path.begin(), path.end()), path.end());
    if (is_closed_path)
      while (path.size() > 1 && path.back() == path.front())
        path.pop_back();
  }

  inline void StripDuplicates(Paths &paths, bool is_closed_path)
  {
    for (typename Paths::iterator paths_citer = paths.begin();
         paths_citer != paths.end(); ++paths_citer)
    {
      StripDuplicates(*paths_citer, is_closed_path);
    }
  }

  // Miscellaneous ------------------------------------------------------------

  inline void CheckPrecision(int &precision, int &error_code)
  {
    if (precision >= -MAX_DECIMAL_PRECISION && precision <= MAX_DECIMAL_PRECISION)
      return;
    error_code |= precision_error_i; // non-fatal error
    DoError(precision_error_i);      // does nothing unless exceptions enabled
    precision = precision > 0 ? MAX_DECIMAL_PRECISION : -MAX_DECIMAL_PRECISION;
  }

  inline void CheckPrecision(int &precision)
  {
    int error_code = 0;
    CheckPrecision(precision, error_code);
  }

  inline double CrossProduct(const Point &pt1, const Point &pt2, const Point &pt3)
  {
    return ((pt2.x - pt1.x) * (pt3.y - pt2.y) - (pt2.y - pt1.y) * (pt3.x - pt2.x));
  }

  inline double CrossProduct(const Point &vec1, const Point &vec2)
  {
    return (vec1.y * vec2.x) - (vec2.y * vec1.x);
  }

  inline double DotProduct(const Point &pt1, const Point &pt2, const Point &pt3)
  {
    return ((pt2.x - pt1.x) * (pt3.x - pt2.x) + (pt2.y - pt1.y) * (pt3.y - pt2.y));
  }

  inline double DotProduct(const Point &vec1, const Point &vec2)
  {
    return (vec1.x * vec2.x) + (vec1.y * vec2.y);
  }

  inline double DistanceSqr(const Point pt1, const Point pt2)
  {
    return Sqr(pt1.x - pt2.x) + Sqr(pt1.y - pt2.y);
  }

  inline double DistanceFromLineSqrd(const Point &pt, const Point &ln1, const Point &ln2)
  {
    // perpendicular distance of point (x³,y³) = (Ax³ + By³ + C)/Sqrt(A² + B²)
    // see http://en.wikipedia.org/wiki/Perpendicular_distance
    double A = ln1.y - ln2.y;
    double B = ln2.x - ln1.x;
    double C = A * ln1.x + B * ln1.y;
    C = A * pt.x + B * pt.y - C;
    return (C * C) / (A * A + B * B);
  }

  inline double Area(const Path &path)
  {
    size_t cnt = path.size();
    if (cnt < 3)
      return 0.0;
    double a = 0.0;
    typename Path::const_iterator it1, it2 = path.cend() - 1, stop = it2;
    if (!(cnt & 1))
      ++stop;
    for (it1 = path.cbegin(); it1 != stop;)
    {
      a += (it2->y + it1->y) * (it2->x - it1->x);
      it2 = it1 + 1;
      a += (it1->y + it2->y) * (it1->x - it2->x);
      it1 += 2;
    }
    if (cnt & 1)
      a += (it2->y + it1->y) * (it2->x - it1->x);
    return a * 0.5;
  }

  inline double Area(const Paths &paths)
  {
    double a = 0.0;
    for (typename Paths::const_iterator paths_iter = paths.cbegin();
         paths_iter != paths.cend(); ++paths_iter)
    {
      a += Area(*paths_iter);
    }
    return a;
  }

  inline bool IsPositive(const Path &poly)
  {
    // A curve has positive orientation [and area] if a region 'R'
    // is on the left when traveling around the outside of 'R'.
    // https://mathworld.wolfram.com/CurveOrientation.html
    // nb: This statement is premised on using Cartesian coordinates
    return Area(poly) >= 0;
  }

  inline bool GetIntersectPoint(const Point &ln1a, const Point &ln1b, const Point &ln2a, const Point &ln2b, Point &ip)
  {
    // https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    double dx1 = ln1b.x - ln1a.x;
    double dy1 = ln1b.y - ln1a.y;
    double dx2 = ln2b.x - ln2a.x;
    double dy2 = ln2b.y - ln2a.y;

    double det = dy1 * dx2 - dy2 * dx1;
    if (det == 0.0)
      return false;
    double t = ((ln1a.x - ln2a.x) * dy2 - (ln1a.y - ln2a.y) * dx2) / det;
    if (t <= 0.0)
      ip = ln1a; // ?? check further (see also #568)
    else if (t >= 1.0)
      ip = ln1b; // ?? check further
    else
    {
      ip.x = ln1a.x + t * dx1;
      ip.y = ln1a.y + t * dy1;
    }
    return true;
  }

  inline bool SegmentsIntersect(const Point &seg1a, const Point &seg1b, const Point &seg2a, const Point &seg2b, bool inclusive = false)
  {
    if (inclusive)
    {
      double res1 = CrossProduct(seg1a, seg2a, seg2b);
      double res2 = CrossProduct(seg1b, seg2a, seg2b);
      if (res1 * res2 > 0)
        return false;
      double res3 = CrossProduct(seg2a, seg1a, seg1b);
      double res4 = CrossProduct(seg2b, seg1a, seg1b);
      if (res3 * res4 > 0)
        return false;
      return (res1 || res2 || res3 || res4); // ensures not collinear
    }
    else
    {
      return (CrossProduct(seg1a, seg2a, seg2b) *
                  CrossProduct(seg1b, seg2a, seg2b) <
              0) &&
             (CrossProduct(seg2a, seg1a, seg1b) *
                  CrossProduct(seg2b, seg1a, seg1b) <
              0);
    }
  }

  inline Point GetClosestPointOnSegment(const Point &offPt, const Point &seg1, const Point &seg2)
  {
    if (seg1.x == seg2.x && seg1.y == seg2.y)
      return seg1;
    double dx = seg2.x - seg1.x;
    double dy = seg2.y - seg1.y;
    double q = ((offPt.x - seg1.x) * dx + (offPt.y - seg1.y) * dy) / (Sqr(dx) + Sqr(dy));
    if (q < 0)
      q = 0;
    else if (q > 1)
      q = 1;
    return Point(seg1.x + nearbyint(q * dx), seg1.y + nearbyint(q * dy));
  }

  enum class PointInPolygonResult
  {
    IsOn,
    IsInside,
    IsOutside
  };

  inline PointInPolygonResult PointInPolygon(const Point &pt, const Path &polygon)
  {
    if (polygon.size() < 3)
      return PointInPolygonResult::IsOutside;

    int val = 0;
    typename Path::const_iterator cbegin = polygon.cbegin(), first = cbegin, curr, prev;
    typename Path::const_iterator cend = polygon.cend();

    while (first != cend && first->y == pt.y)
      ++first;
    if (first == cend) // not a proper polygon
      return PointInPolygonResult::IsOutside;

    bool is_above = first->y < pt.y, starting_above = is_above;
    curr = first + 1;
    while (true)
    {
      if (curr == cend)
      {
        if (cend == first || first == cbegin)
          break;
        cend = first;
        curr = cbegin;
      }

      if (is_above)
      {
        while (curr != cend && curr->y < pt.y)
          ++curr;
        if (curr == cend)
          continue;
      }
      else
      {
        while (curr != cend && curr->y > pt.y)
          ++curr;
        if (curr == cend)
          continue;
      }

      if (curr == cbegin)
        prev = polygon.cend() - 1; // nb: NOT cend (since might equal first)
      else
        prev = curr - 1;

      if (curr->y == pt.y)
      {
        if (curr->x == pt.x ||
            (curr->y == prev->y &&
             ((pt.x < prev->x) != (pt.x < curr->x))))
          return PointInPolygonResult::IsOn;
        ++curr;
        if (curr == first)
          break;
        continue;
      }

      if (pt.x < curr->x && pt.x < prev->x)
      {
        // we're only interested in edges crossing on the left
      }
      else if (pt.x > prev->x && pt.x > curr->x)
        val = 1 - val; // toggle val
      else
      {
        double d = CrossProduct(*prev, *curr, pt);
        if (d == 0)
          return PointInPolygonResult::IsOn;
        if ((d < 0) == is_above)
          val = 1 - val;
      }
      is_above = !is_above;
      ++curr;
    }

    if (is_above != starting_above)
    {
      cend = polygon.cend();
      if (curr == cend)
        curr = cbegin;
      if (curr == cbegin)
        prev = cend - 1;
      else
        prev = curr - 1;
      double d = CrossProduct(*prev, *curr, pt);
      if (d == 0)
        return PointInPolygonResult::IsOn;
      if ((d < 0) == is_above)
        val = 1 - val;
    }

    return (val == 0) ? PointInPolygonResult::IsOutside : PointInPolygonResult::IsInside;
  }

} // namespace
