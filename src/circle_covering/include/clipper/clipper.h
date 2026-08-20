#pragma once


#include <cstdlib>
#include <type_traits>
#include <vector>
 
#include "clipper.offset.h"

namespace Clipper2Lib
{

  inline Paths BooleanOp(ClipType cliptype, FillRule fillrule,
                         const Paths &subjects, const Paths &clips, int precision = 7)
  {
    int error_code = 0;
    CheckPrecision(precision, error_code);
    Paths result;
    if (error_code)
      return result;
    Clipper clipper(precision);
    clipper.AddSubject(subjects);
    clipper.AddClip(clips);
    clipper.Execute(cliptype, fillrule, result);
    return result;
  }

  inline void BooleanOp(ClipType cliptype, FillRule fillrule,
                        const Paths &subjects, const Paths &clips,
                        PolyTree &polytree, int precision = 7)
  {
    polytree.Clear();
    int error_code = 0;
    CheckPrecision(precision, error_code);
    if (error_code)
      return;
    Clipper clipper(precision);
    clipper.AddSubject(subjects);
    clipper.AddClip(clips);
    clipper.Execute(cliptype, fillrule, polytree);
  }

  inline Paths Intersect(const Paths &subjects, const Paths &clips, FillRule fillrule, int decimal_prec = 10)
  {
    return BooleanOp(ClipType::Intersection, fillrule, subjects, clips, decimal_prec);
  }

  inline Paths Union(const Paths &subjects, const Paths &clips, FillRule fillrule, int decimal_prec = 7)
  {
    return BooleanOp(ClipType::Union, fillrule, subjects, clips, decimal_prec);
  }

  inline Paths Union(const Paths &subjects, FillRule fillrule, int precision = 7)
  {
    Paths result;
    int error_code = 0;
    CheckPrecision(precision, error_code);
    if (error_code)
      return result;
    Clipper clipper(precision);
    clipper.AddSubject(subjects);
    clipper.Execute(ClipType::Union, fillrule, result);
    return result;
  }

  inline Paths Difference(const Paths &subjects, const Paths &clips, FillRule fillrule, int decimal_prec = 10)
  {
    return BooleanOp(ClipType::Difference, fillrule, subjects, clips, decimal_prec);
  }

  inline Paths Xor(const Paths &subjects, const Paths &clips, FillRule fillrule, int decimal_prec = 7)
  {
    return BooleanOp(ClipType::Xor, fillrule, subjects, clips, decimal_prec);
  }

  inline Paths InflatePaths(const Paths &paths, double delta,
                            JoinType jt, EndType et, double miter_limit = 2.0,
                            int precision = 7, double arc_tolerance = 0.0)
  {
    int error_code = 0;
    CheckPrecision(precision, error_code);
    if (!delta)
      return paths;
    if (error_code)
      return Paths();
    const double scale = std::pow(10, precision);
    ClipperOffset clip_offset(miter_limit, arc_tolerance);
    clip_offset.AddPaths(ScalePaths(paths, scale, error_code), jt, et);
    if (error_code)
      return Paths();
    Paths solution;
    clip_offset.Execute(delta * scale, solution);
    return ScalePaths(solution, 1 / scale, error_code);
  }

  inline Path TranslatePath(const Path &path, double dx, double dy)
  {
    Path result;
    result.reserve(path.size());
    std::transform(path.begin(), path.end(), back_inserter(result),
                   [dx, dy](const auto &pt)
                   { return Point(pt.x + dx, pt.y + dy); });
    return result;
  }
  inline Paths TranslatePaths(const Paths &paths, double dx, double dy)
  {
    Paths result;
    result.reserve(paths.size());
    std::transform(paths.begin(), paths.end(), back_inserter(result),
                   [dx, dy](const auto &path)
                   { return TranslatePath(path, dx, dy); });
    return result;
  }

  namespace details
  {

    inline void PolyPathToPaths(const PolyPath &polypath, Paths &paths)
    {
      paths.push_back(polypath.Polygon());
      for (const auto &child : polypath)
        PolyPathToPaths(*child, paths);
    }

    static void OutlinePolyPath(std::ostream &os,
                                size_t idx, bool isHole, size_t count, const std::string &preamble)
    {
      std::string plural = (count == 1) ? "." : "s.";
      if (isHole)
        os << preamble << "+- Hole (" << idx << ") contains " << count << " nested polygon" << plural << std::endl;
      else
        os << preamble << "+- Polygon (" << idx << ") contains " << count << " hole" << plural << std::endl;
    }

    static void OutlinePolyPath(std::ostream &os, const PolyPath &pp,
                                size_t idx, std::string preamble)
    {
      OutlinePolyPath(os, idx, pp.IsHole(), pp.Count(), preamble);
      for (size_t i = 0; i < pp.Count(); ++i)
        if (pp.Child(i)->Count())
          details::OutlinePolyPath(os, *pp.Child(i), i, preamble + "  ");
    }

  } // end details namespace

  inline std::ostream &operator<<(std::ostream &os, const PolyTree &pp)
  {
    std::string plural = (pp.Count() == 1) ? " polygon." : " polygons.";
    os << std::endl
       << "Polytree with " << pp.Count() << plural << std::endl;
    for (size_t i = 0; i < pp.Count(); ++i)
      if (pp.Child(i)->Count())
        details::OutlinePolyPath(os, *pp.Child(i), i, "  ");
    os << std::endl
       << std::endl;
    if (!pp.Level())
      os << std::endl;
    return os;
  }

  inline Paths PolyTreeToPathsD(const PolyTree &polytree)
  {
    Paths result;
    for (const auto &child : polytree)
      details::PolyPathToPaths(*child, result);
    return result;
  }

  inline Path MakePath(const std::vector<Point> &list)
  {
    Path result = list;
    return result;
  }

  inline Path TrimCollinear(const Path &path, int precision, bool is_open_path = false)
  {
    int error_code = 0;
    CheckPrecision(precision, error_code);
    if (error_code)
      return Path();
    const double scale = std::pow(10, precision);
    Path p = ScalePath(path, scale, error_code);
    if (error_code)
      return Path();
    p = TrimCollinear(p, is_open_path);
    return ScalePath(p, 1 / scale, error_code);
  }

  inline double Distance(const Point pt1, const Point pt2)
  {
    return std::sqrt(DistanceSqr(pt1, pt2));
  }

  inline double Length(const Path &path, bool is_closed_path = false)
  {
    double result = 0.0;
    if (path.size() < 2)
      return result;
    auto it = path.cbegin(), stop = path.end() - 1;
    for (; it != stop; ++it)
      result += Distance(*it, *(it + 1));
    if (is_closed_path)
      result += Distance(*stop, *path.cbegin());
    return result;
  }

  inline bool NearCollinear(const Point &pt1, const Point &pt2, const Point &pt3, double sin_sqrd_min_angle_rads)
  {
    double cp = std::abs(CrossProduct(pt1, pt2, pt3));
    return (cp * cp) / (DistanceSqr(pt1, pt2) * DistanceSqr(pt2, pt3)) < sin_sqrd_min_angle_rads;
  }

  inline Path Ellipse(const Point &center, double radiusX, double radiusY = 0, int steps = 0)
  {
    if (radiusX <= 0)
      return Path();
    if (radiusY <= 0)
      radiusY = radiusX;
    if (steps <= 2)
      steps = static_cast<int>(PI * sqrt((radiusX + radiusY) / 2));

    double si = std::sin(2 * PI / steps);
    double co = std::cos(2 * PI / steps);
    double dx = co, dy = si;
    Path result;
    result.reserve(steps);
    result.push_back(Point(center.x + radiusX, center.y));
    for (int i = 1; i < steps; ++i)
    {
      result.push_back(Point(center.x + radiusX * dx, center.y + radiusY * dy));
      double x = dx * co - dy * si;
      dy = dy * co + dx * si;
      dx = x;
    }
    return result;
  }

  inline double PerpendicDistFromLineSqrd(const Point &pt, const Point &line1, const Point &line2)
  {
    double a = static_cast<double>(pt.x - line1.x);
    double b = static_cast<double>(pt.y - line1.y);
    double c = static_cast<double>(line2.x - line1.x);
    double d = static_cast<double>(line2.y - line1.y);
    if (c == 0 && d == 0)
      return 0;
    return Sqr(a * d - c * b) / (c * c + d * d);
  }

  inline size_t GetNext(size_t current, size_t high, const std::vector<bool> &flags)
  {
    ++current;
    while (current <= high && flags[current])
      ++current;
    if (current <= high)
      return current;
    current = 0;
    while (flags[current])
      ++current;
    return current;
  }

  inline size_t GetPrior(size_t current, size_t high, const std::vector<bool> &flags)
  {
    if (current == 0)
      current = high;
    else
      --current;
    while (current > 0 && flags[current])
      --current;
    if (!flags[current])
      return current;
    current = high;
    while (flags[current])
      --current;
    return current;
  }

  inline Path SimplifyPath(const Path path, double epsilon, bool isClosedPath = true)
  {
    const size_t len = path.size(), high = len - 1;
    const double epsSqr = Sqr(epsilon);
    if (len < 4)
      return Path(path);

    std::vector<bool> flags(len);
    std::vector<double> distSqr(len);
    size_t prior = high, curr = 0, start, next, prior2, next2;
    if (isClosedPath)
    {
      distSqr[0] = PerpendicDistFromLineSqrd(path[0], path[high], path[1]);
      distSqr[high] = PerpendicDistFromLineSqrd(path[high], path[0], path[high - 1]);
    }
    else
    {
      distSqr[0] = MAX_DBL;
      distSqr[high] = MAX_DBL;
    }
    for (size_t i = 1; i < high; ++i)
      distSqr[i] = PerpendicDistFromLineSqrd(path[i], path[i - 1], path[i + 1]);

    for (;;)
    {
      if (distSqr[curr] > epsSqr)
      {
        start = curr;
        do
        {
          curr = GetNext(curr, high, flags);
        } while (curr != start && distSqr[curr] > epsSqr);
        if (curr == start)
          break;
      }

      prior = GetPrior(curr, high, flags);
      next = GetNext(curr, high, flags);
      if (next == prior)
        break;

      if (distSqr[next] < distSqr[curr])
      {
        flags[next] = true;
        next = GetNext(next, high, flags);
        next2 = GetNext(next, high, flags);
        distSqr[curr] = PerpendicDistFromLineSqrd(path[curr], path[prior], path[next]);
        if (next != high || isClosedPath)
          distSqr[next] = PerpendicDistFromLineSqrd(path[next], path[curr], path[next2]);
        curr = next;
      }
      else
      {
        flags[curr] = true;
        curr = next;
        next = GetNext(next, high, flags);
        prior2 = GetPrior(prior, high, flags);
        distSqr[curr] = PerpendicDistFromLineSqrd(path[curr], path[prior], path[next]);
        if (prior != 0 || isClosedPath)
          distSqr[prior] = PerpendicDistFromLineSqrd(path[prior], path[prior2], path[curr]);
      }
    }
    Path result;
    result.reserve(len);
    for (typename Path::size_type i = 0; i < len; ++i)
      if (!flags[i])
        result.push_back(path[i]);
    return result;
  }

  inline Paths SimplifyPaths(const Paths paths, double epsilon, bool isClosedPath = true)
  {
    Paths result;
    result.reserve(paths.size());
    for (const auto &path : paths)
      result.push_back(SimplifyPath(path, epsilon, isClosedPath));
    return result;
  }

  inline void RDP(const Path path, std::size_t begin, std::size_t end, double epsSqrd, std::vector<bool> &flags)
  {
    typename Path::size_type idx = 0;
    double max_d = 0;
    while (end > begin && path[begin] == path[end])
      flags[end--] = false;
    for (typename Path::size_type i = begin + 1; i < end; ++i)
    {
      // PerpendicDistFromLineSqrd - avoids expensive Sqrt()
      double d = PerpendicDistFromLineSqrd(path[i], path[begin], path[end]);
      if (d <= max_d)
        continue;
      max_d = d;
      idx = i;
    }
    if (max_d <= epsSqrd)
      return;
    flags[idx] = true;
    if (idx > begin + 1)
      RDP(path, begin, idx, epsSqrd, flags);
    if (idx < end - 1)
      RDP(path, idx, end, epsSqrd, flags);
  }

  inline Path RamerDouglasPeucker(const Path &path, double epsilon)
  {
    const typename Path::size_type len = path.size();
    if (len < 5)
      return Path(path);
    std::vector<bool> flags(len);
    flags[0] = true;
    flags[len - 1] = true;
    RDP(path, 0, len - 1, Sqr(epsilon), flags);
    Path result;
    result.reserve(len);
    for (typename Path::size_type i = 0; i < len; ++i)
      if (flags[i])
        result.push_back(path[i]);
    return result;
  }

  inline Paths RamerDouglasPeucker(const Paths &paths, double epsilon)
  {
    Paths result;
    result.reserve(paths.size());
    std::transform(paths.begin(), paths.end(), back_inserter(result),
                   [epsilon](const auto &path)
                   { return RamerDouglasPeucker(path, epsilon); });
    return result;
  }

} // end Clipper2Lib namespace
