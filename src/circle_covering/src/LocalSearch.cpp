#include "LocalSearch.h"

// GLOBAL
extern bool CONSTRAINT;
extern double R;
extern std::vector<Point> border;
extern std::vector<std::vector<Point>> holes;
extern std::vector<std::vector<int>> neighbors;
extern std::vector<std::vector<std::vector<Point>>> borders;

// LOCAL
double penalty_index;
std::vector<std::vector<int>> vertex_neighbors;
std::vector<Point> sub_points;
double d_star;

// gradient
Point coord_point_to_line_segment(Point a, Point b, Point point)
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

std::pair<Point, Point> find_line(Point p)
{
    Point a, b;
    for (size_t i = 0; i < border.size(); i++)
    {

        if (are_points_collinear(p, border[i], border[(i + 1) % border.size()]))
        {
            a = border[i];
            b = border[(i + 1) % border.size()];
        }
    }

    for (size_t i = 0; i < holes.size(); i++)
    {
        for (size_t j = 0; j < holes[i].size(); j++)
        {
            if (are_points_collinear(p, holes[i][j], holes[i][(j + 1) % holes[i].size()]))
            {
                a = holes[i][j];
                b = holes[i][(j + 1) % holes[i].size()];
            }
        }
    }

    return std::make_pair(a, b);
}

std::pair<double, double> duv_x(const std::vector<Point> &points, int n, int t, Point vertex)
{
    double x = points[n].x;
    double y = points[n].y;
    double k0, b0, k1, b1, dk1, db1;
    double du = 0.0;
    double dv = 0.0;
    // case 1
    if (vertex_neighbors[t].empty())
    {
        du = 0.0;
        dv = 0.0;
    }

    // case 2
    else if (vertex_neighbors[t].size() == 1)
    {
        double x1 = points[vertex_neighbors[t][0]].x;
        double y1 = points[vertex_neighbors[t][0]].y;
        // if Y = k1*X+b1 is vertical k1 does not exist
        if (fabs(x1 - x) < 1e-8)
        {
            du = 0.0;
            dv = 0.0;
        }
        else
        {
            // find the line
            Point a, b;
            std::pair<Point, Point> line = find_line(vertex);
            a = line.first;
            b = line.second;
            if (fabs(a.x - b.x) > 1e-8)
            {
                k0 = (a.y - b.y) / (a.x - b.x);
                b0 = a.y - k0 * a.x;
            }
            // if Y = k1*X+b1 is horizontal k1 does not exist
            if (fabs(y1 - y) < 1e-8)
            {
                du = 0.5;
                dv = k0 / 2;
            }
            else
            {
                k1 = (x1 - x) / (y - y1);
                b1 = (y + y1) / 2 - k1 * ((x + x1) / 2);
                dk1 = 1 / (y1 - y);
                db1 = x / (y - y1);
                // if k0 does not exist
                if (fabs(a.x - b.x) < 1e-8)
                {
                    du = 0.0;
                    dv = a.x * dk1 + db1;
                }
                else
                {
                    du = (1 / ((k1 - k0) * (k1 - k0))) * (-1 * db1 * (k1 - k0) - (b0 - b1) * dk1);
                    dv = (1 / ((k1 - k0) * (k1 - k0))) * ((b0 * dk1 - k0 * db1) * (k1 - k0) - (k1 * b0 - k0 * b1) * dk1);
                }
            }
        }
    }

    // case 3
    else if (vertex_neighbors[t].size() == 2)
    {

        double x1 = points[vertex_neighbors[t][0]].x;
        double y1 = points[vertex_neighbors[t][0]].y;
        double x2 = points[vertex_neighbors[t][1]].x;
        double y2 = points[vertex_neighbors[t][1]].y;
        double a = x - x1;
        double b = y - y1;
        double c = x - x2;
        double d = y - y2;
        double e = ((x * x - x1 * x1) - (y1 * y1 - y * y)) / 2;
        double f = ((x * x - x2 * x2) - (y2 * y2 - y * y)) / 2;
        double u = (d * e - b * f) / (a * d - b * c);
        double v = (a * f - c * e) / (a * d - b * c);

        du = (1 / ((a * d - b * c) * (a * d - b * c))) * (d - b) * (x * (a * d - b * c) - (d * e - b * f));
        dv = (1 / ((a * d - b * c) * (a * d - b * c))) * ((f - e + a * x - c * x) * (a * d - b * c) - (a * f - c * e) * (d - b));
    }

    return std::make_pair(du, dv);
}

std::pair<double, double> duv_y(const std::vector<Point> &points, int n, int t, Point vertex)
{
    double x = points[n].x;
    double y = points[n].y;
    double k0, b0, k1, b1, dk1, db1;
    double du = 0.0;
    double dv = 0.0;
    // case 1
    if (vertex_neighbors[t].empty())
    {
        du = 0.0;
        dv = 0.0;
    }

    // case 2
    else if (vertex_neighbors[t].size() == 1)
    {
        double x1 = points[vertex_neighbors[t][0]].x;
        double y1 = points[vertex_neighbors[t][0]].y;
        // if Y = k1*X+b1 is horizontal k1 does not exist
        if (fabs(y1 - y) < 1e-8)
        {
            du = 0.0;
            dv = 0.0;
        }
        else
        {
            // find the line
            Point a, b;
            std::pair<Point, Point> line = find_line(vertex);
            a = line.first;
            b = line.second;

            if (fabs(a.x - b.x) > 1e-8)
            {
                k0 = (a.y - b.y) / (a.x - b.x);
                b0 = a.y - k0 * a.x;
            }
            // if Y = k1*X+b1 is vertical k1 = 0
            if (fabs(x1 - x) < 1e-8)
            {
                du = fabs(a.x - b.x) > 1e-8 ? 1 / (2 * k0) : 0.0;
                dv = 0.5;
            }
            else
            {
                k1 = (x1 - x) / (y - y1);
                b1 = (y + y1) / 2 - k1 * ((x + x1) / 2);
                dk1 = (x - x1) / ((y - y1) * (y - y1));
                db1 = 0.5 - ((x * x - x1 * x1) / (2 * (y - y1) * (y - y1)));
                // if k0 does not exist
                if (fabs(a.x - b.x) < 1e-8)
                {
                    du = 0.0;
                    dv = a.x * dk1 + db1;
                }
                else
                {
                    du = (1 / ((k1 - k0) * (k1 - k0))) * (-1 * db1 * (k1 - k0) - (b0 - b1) * dk1);
                    dv = (1 / ((k1 - k0) * (k1 - k0))) * ((b0 * dk1 - k0 * db1) * (k1 - k0) - (k1 * b0 - k0 * b1) * dk1);
                }
            }
        }
    }

    // case 3
    else if (vertex_neighbors[t].size() == 2)
    {

        double x1 = points[vertex_neighbors[t][0]].x;
        double y1 = points[vertex_neighbors[t][0]].y;
        double x2 = points[vertex_neighbors[t][1]].x;
        double y2 = points[vertex_neighbors[t][1]].y;
        double a = x - x1;
        double b = y - y1;
        double c = x - x2;
        double d = y - y2;
        double e = ((x * x - x1 * x1) - (y1 * y1 - y * y)) / 2;
        double f = ((x * x - x2 * x2) - (y2 * y2 - y * y)) / 2;
        double u = (d * e - b * f) / (a * d - b * c);
        double v = (a * f - c * e) / (a * d - b * c);

        du = (1 / ((a * d - b * c) * (a * d - b * c))) * ((e - f + d * y - b * y) * (a * d - b * c) - (d * e - b * f) * (a - c));
        dv = (1 / ((a * d - b * c) * (a * d - b * c))) * (a - c) * (y * (a * d - b * c) - (a * f - c * e));
    }

    return std::make_pair(du, dv);
}

double computer_x_gradient(int n, const std::vector<Point> &points, double r, const std::vector<Point> &voronoi)
{
    double g = 0;
    double dist, e;
    double x_i, y_i, u_j, v_j;
    double du, dv;
    for (size_t i = 0; i < voronoi.size(); i++)
    {
        x_i = points[n].x;
        y_i = points[n].y;
        u_j = voronoi[i].x;
        v_j = voronoi[i].y;
        dist = points[n].getDistance(voronoi[i]);
        e = dist - r;
        if (e > 0)
        {
            std::pair<double, double> duvx = duv_x(points, n, i, voronoi[i]);
            du = duvx.first;
            dv = duvx.second;
            g += 2 * e * ((x_i - u_j) * (1 - du) - (y_i - v_j) * dv) / dist;
            for (size_t j = 0; j < vertex_neighbors[i].size(); j++)
            {
                int t = vertex_neighbors[i][j];
                x_i = points[t].x;
                y_i = points[t].y;
                g -= 2 * e * ((x_i - u_j) * du + (y_i - v_j) * dv) / dist;
            }
        }
    }

    return g;
}

double computer_y_gradient(int n, const std::vector<Point> &points, double r, const std::vector<Point> &voronoi)
{
    double g = 0;
    double dist, e;
    double x_i, y_i, u_j, v_j;
    double du, dv;
    for (size_t i = 0; i < voronoi.size(); i++)
    {
        x_i = points[n].x;
        y_i = points[n].y;
        u_j = voronoi[i].x;
        v_j = voronoi[i].y;
        dist = points[n].getDistance(voronoi[i]);
        e = dist - r;
        if (e > 0)
        {
            std::pair<double, double> duvy = duv_y(points, n, i, voronoi[i]);
            du = duvy.first;
            dv = duvy.second;
            g += 2 * e * (-1 * (x_i - u_j) * du + (y_i - v_j) * (1 - dv)) / dist;
            for (size_t j = 0; j < vertex_neighbors[i].size(); j++)
            {
                int t = vertex_neighbors[i][j];
                x_i = points[t].x;
                y_i = points[t].y;
                g -= 2 * e * ((x_i - u_j) * du + (y_i - v_j) * dv) / dist;
            }
        }
    }

    return g;
}

// descent
double myvalue(double *x, INT n)
{
    double f = 0.0;
    double dist, e;
    int N = n / 2; // the number of points

    std::vector<Point> points;
    for (size_t i = 0; i < N; i++)
        points.push_back(Point(x[2 * i], x[2 * i + 1]));

    VoronoiDiagram result = generate_diagram(points);
    std::vector<std::vector<Point>> voronoi = result.result();
    Clip res(voronoi, border, holes);
    voronoi = res.clip();

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < voronoi[i].size(); j++)
        {
            dist = points[i].getDistance(voronoi[i][j]);
            e = dist - R;
            if (e > 0)
                f += e * e;
        }
    }

    d_star = CONSTRAINT ? 1e-10 : R;
    double d;
    for (size_t i = 0; i < N; i++)
    {
        if (!points[i].is_in_polygon(border))
        {
            d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < border.size(); j++)
            {
                Point p_new = coord_point_to_line_segment(border[j], border[(j + 1) % border.size()], points[i]);
                double d_new = p_new.getDistance(points[i]);
                if (d_new < d)
                    d = d_new;
            }
            if (d > d_star)
                f += d * d;
        }
        else
        {
            for (size_t j = 0; j < holes.size(); j++)
            {
                if (points[i].is_in_polygon(holes[j]))
                {
                    d = std::numeric_limits<double>::infinity();
                    for (size_t k = 0; k < holes[j].size(); k++)
                    {
                        Point p_new = coord_point_to_line_segment(holes[j][k], holes[j][(k + 1) % holes[j].size()], points[i]);
                        double d_new = p_new.getDistance(points[i]);
                        if (d_new < d)
                            d = d_new;
                    }
                    if (d > d_star)
                        f += d * d;
                }
            }
        }
    }
    return f;
}

void mygrad(double *g, double *x, INT n)
{
    double dist, e;
    double x_i, y_i, u_j, v_j;
    int N = n / 2;

    std::vector<Point> points;
    for (size_t i = 0; i < N; i++)
        points.push_back(Point(x[2 * i], x[2 * i + 1]));

    VoronoiDiagram result = generate_diagram(points);
    std::vector<std::vector<Point>> voronoi = result.result();
    Clip res(voronoi, border, holes);
    voronoi = res.clip();

    for (size_t i = 0; i < n; i++)
        g[i] = 0.0;

    for (size_t i = 0; i < N; i++)
    {
        x_i = points[i].x;
        y_i = points[i].y;
        for (size_t j = 0; j < voronoi[i].size(); j++)
        {
            u_j = voronoi[i][j].x;
            v_j = voronoi[i][j].y;
            dist = points[i].getDistance(voronoi[i][j]);
            e = dist - R;
            if (e > 0)
            {
                g[2 * i] += 2.0 * e * (x_i - u_j) / dist;
                g[2 * i + 1] += 2.0 * e * (y_i - v_j) / dist;
            }
        }
    }

    d_star = CONSTRAINT ? 1e-10 : R;
    double d, d_x, d_y, x1, y1, x2, y2;
    for (size_t i = 0; i < N; i++)
    {
        x_i = points[i].x;
        y_i = points[i].y;
        d = 0.0;
        if (!points[i].is_in_polygon(border))
        {
            d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < border.size(); j++)
            {
                Point p_new = coord_point_to_line_segment(border[j], border[(j + 1) % border.size()], points[i]);
                double d_new = p_new.getDistance(points[i]);
                if (d_new < d)
                {
                    d = d_new;
                    d_x = p_new.x;
                    d_y = p_new.y;
                    x1 = border[j].x;
                    y1 = border[j].y;
                    x2 = border[(j + 1) % border.size()].x;
                    y2 = border[(j + 1) % border.size()].y;
                }
            }
        }
        else
        {
            for (size_t j = 0; j < holes.size(); j++)
            {
                if (points[i].is_in_polygon(holes[j]))
                {
                    d = std::numeric_limits<double>::infinity();
                    for (size_t k = 0; k < holes[j].size(); k++)
                    {
                        Point p_new = coord_point_to_line_segment(holes[j][k], holes[j][(k + 1) % holes[j].size()], points[i]);
                        double d_new = p_new.getDistance(points[i]);
                        if (d_new < d)
                        {
                            d = d_new;
                            d_x = p_new.x;
                            d_y = p_new.y;
                            x1 = holes[j][k].x;
                            y1 = holes[j][k].y;
                            x2 = holes[j][(k + 1) % holes[j].size()].x;
                            y2 = holes[j][(k + 1) % holes[j].size()].y;
                        }
                    }
                }
            }
        }

        if (d > d_star)
        {
            if (fabs(d_x - x1) < 1e-10 || fabs(d_x - x2) < 1e-10)
            {
                g[2 * i] += 2 * (x_i - d_x);
                g[2 * i + 1] += 2 * (y_i - d_y);
            }
            else
            {
                double A = y2 - y1;
                double B = x1 - x2;
                double C = x2 * y1 - x1 * y2;
                g[2 * i] += 2 * A * (A * x_i + B * y_i + C) / (A * A + B * B);
                g[2 * i + 1] += 2 * B * (A * x_i + B * y_i + C) / (A * A + B * B);
            }
        }
    }
}

void mygrad_h(double *g, double *x, INT n)
{
    double x_i, y_i;
    int N = n / 2;
    std::vector<Point> points;
    for (size_t i = 0; i < N; i++)
        points.push_back(Point(x[2 * i], x[2 * i + 1]));

    VoronoiDiagram result = generate_diagram(points);
    std::vector<std::vector<Point>> voronoi = result.result_with_neighbors();
    Clip res(voronoi, border, holes);
    voronoi = res.clip();

    for (size_t i = 0; i < n; i++)
        g[i] = 0.0;

    for (size_t i = 0; i < N; i++)
    {
        if (!voronoi[i].empty())
        {
            sub_points.clear();
            sub_points.push_back(points[i]);
            for (size_t j = 0; j < neighbors[i].size(); j++)
            {
                if (neighbors[i][j] >= 0)
                    sub_points.push_back(points[neighbors[i][j]]);
            }
            vertex_neighbors.clear();
            for (size_t j = 0; j < voronoi[i].size(); j++)
            {
                std::vector<int> neighbor;
                double d1 = points[i].getDistance(voronoi[i][j]);
                for (size_t k = 0; k < neighbors[i].size(); k++)
                {
                    if (neighbors[i][k] >= 0)
                    {
                        double d2 = points[neighbors[i][k]].getDistance(voronoi[i][j]);
                        if (fabs(d1 - d2) < 1e-9)
                            neighbor.push_back(neighbors[i][k]);
                    }
                }
                if (neighbor.size() > 2)
                {
                    std::vector<int> new_neighbor;
                    for (size_t k = 0; k < neighbor.size(); k++)
                    {
                        Point p = 0.5 * (points[i] + points[neighbor[k]]);
                        if (p.getDistance(voronoi[i][j]) > 1e-10)
                            new_neighbor.push_back(neighbor[k]);
                    }
                    neighbor = new_neighbor;
                }
                vertex_neighbors.push_back(neighbor);
            }
            g[2 * i] = computer_x_gradient(i, points, R, voronoi[i]);
            g[2 * i + 1] = computer_y_gradient(i, points, R, voronoi[i]);
        }
    }
    d_star = CONSTRAINT ? 1e-10 : R;
    double d, d_x, d_y, x1, y1, x2, y2;
    for (size_t i = 0; i < N; i++)
    {
        x_i = points[i].x;
        y_i = points[i].y;
        d = 0.0;
        if (!points[i].is_in_polygon(border))
        {
            d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < border.size(); j++)
            {
                Point p_new = coord_point_to_line_segment(border[j], border[(j + 1) % border.size()], points[i]);
                double d_new = p_new.getDistance(points[i]);
                if (d_new < d)
                {
                    d = d_new;
                    d_x = p_new.x;
                    d_y = p_new.y;
                    x1 = border[j].x;
                    y1 = border[j].y;
                    x2 = border[(j + 1) % border.size()].x;
                    y2 = border[(j + 1) % border.size()].y;
                }
            }
        }
        else
        {
            for (size_t j = 0; j < holes.size(); j++)
            {
                if (points[i].is_in_polygon(holes[j]))
                {
                    d = std::numeric_limits<double>::infinity();
                    for (size_t k = 0; k < holes[j].size(); k++)
                    {
                        Point p_new = coord_point_to_line_segment(holes[j][k], holes[j][(k + 1) % holes[j].size()], points[i]);
                        double d_new = p_new.getDistance(points[i]);
                        if (d_new < d)
                        {
                            d = d_new;
                            d_x = p_new.x;
                            d_y = p_new.y;
                            x1 = holes[j][k].x;
                            y1 = holes[j][k].y;
                            x2 = holes[j][(k + 1) % holes[j].size()].x;
                            y2 = holes[j][(k + 1) % holes[j].size()].y;
                        }
                    }
                }
            }
        }

        if (d > d_star)
        {
            if (fabs(d_x - x1) < 1e-10 || fabs(d_x - x2) < 1e-10)
            {
                g[2 * i] += 2 * (x_i - d_x);
                g[2 * i + 1] += 2 * (y_i - d_y);
            }
            else
            {
                double A = y2 - y1;
                double B = x1 - x2;
                double C = x2 * y1 - x1 * y2;
                g[2 * i] += 2 * A * (A * x_i + B * y_i + C) / (A * A + B * B);
                g[2 * i + 1] += 2 * B * (A * x_i + B * y_i + C) / (A * A + B * B);
            }
        }
    }
}

double myvalueR(double *x, INT n)
{
    double f = 0.0;
    double dist, e;
    int N = n / 2; // the number of points
    std::vector<Point> points;
    for (size_t i = 0; i < N; i++)
        points.push_back(Point(x[2 * i], x[2 * i + 1]));

    VoronoiDiagram result = generate_diagram(points);
    std::vector<std::vector<Point>> voronoi = result.result();
    Clip res(voronoi, border, holes);
    voronoi = res.clip();

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < voronoi[i].size(); j++)
        {
            dist = points[i].getDistance(voronoi[i][j]);
            e = dist - x[n - 1];
            if (e > 0)
                f += e * e;
        }
    }

    d_star = CONSTRAINT ? 1e-10 : x[n - 1];
    double d;
    for (size_t i = 0; i < N; i++)
    {
        if (!points[i].is_in_polygon(border))
        {
            d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < border.size(); j++)
            {
                Point p_new = coord_point_to_line_segment(border[j], border[(j + 1) % border.size()], points[i]);
                double d_new = p_new.getDistance(points[i]);
                if (d_new < d)
                    d = d_new;
            }
            if (d > d_star)
                f += d * d;
        }
        else
        {
            for (size_t j = 0; j < holes.size(); j++)
            {
                if (points[i].is_in_polygon(holes[j]))
                {
                    d = std::numeric_limits<double>::infinity();
                    for (size_t k = 0; k < holes[j].size(); k++)
                    {
                        Point p_new = coord_point_to_line_segment(holes[j][k], holes[j][(k + 1) % holes[j].size()], points[i]);
                        double d_new = p_new.getDistance(points[i]);
                        if (d_new < d)
                            d = d_new;
                    }
                    if (d > d_star)
                        f += d * d;
                }
            }
        }
    }
    f += x[n - 1] * x[n - 1] / penalty_index;
    return f;
}

void mygradR(double *g, double *x, INT n)
{
    double x_i, y_i;
    int N = n / 2;
    std::vector<Point> points;
    for (size_t i = 0; i < N; i++)
        points.push_back(Point(x[2 * i], x[2 * i + 1]));
    VoronoiDiagram result = generate_diagram(points);
    std::vector<std::vector<Point>> voronoi = result.result_with_neighbors();
    Clip res(voronoi, border, holes);
    voronoi = res.clip();
    for (size_t i = 0; i < n; i++)
        g[i] = 0.0;
    for (size_t i = 0; i < N; i++)
    {
        if (!voronoi[i].empty())
        {
            sub_points.clear();
            sub_points.push_back(points[i]);
            for (size_t j = 0; j < neighbors[i].size(); j++)
            {
                if (neighbors[i][j] >= 0)
                    sub_points.push_back(points[neighbors[i][j]]);
            }
            vertex_neighbors.clear();
            for (size_t j = 0; j < voronoi[i].size(); j++)
            {
                std::vector<int> neighbor;
                double d1 = points[i].getDistance(voronoi[i][j]);
                for (size_t k = 0; k < neighbors[i].size(); k++)
                {
                    if (neighbors[i][k] >= 0)
                    {
                        double d2 = points[neighbors[i][k]].getDistance(voronoi[i][j]);
                        if (fabs(d1 - d2) < 1e-9)
                            neighbor.push_back(neighbors[i][k]);
                    }
                }
                if (neighbor.size() > 2)
                {
                    std::vector<int> new_neighbor;
                    for (size_t k = 0; k < neighbor.size(); k++)
                    {
                        Point p = 0.5 * (points[i] + points[neighbor[k]]);
                        if (p.getDistance(voronoi[i][j]) > 1e-10)
                            new_neighbor.push_back(neighbor[k]);
                    }
                    neighbor = new_neighbor;
                }
                vertex_neighbors.push_back(neighbor);
            }
            g[2 * i] = computer_x_gradient(i, points, x[n - 1], voronoi[i]);
            g[2 * i + 1] = computer_y_gradient(i, points, x[n - 1], voronoi[i]);
        }
    }

    d_star = CONSTRAINT ? 1e-10 : x[n - 1];
    double d, d_x, d_y, x1, y1, x2, y2;
    for (size_t i = 0; i < N; i++)
    {
        x_i = points[i].x;
        y_i = points[i].y;
        d = 0.0;
        if (!points[i].is_in_polygon(border))
        {
            d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < border.size(); j++)
            {
                Point p_new = coord_point_to_line_segment(border[j], border[(j + 1) % border.size()], points[i]);
                double d_new = p_new.getDistance(points[i]);
                if (d_new < d)
                {
                    d = d_new;
                    d_x = p_new.x;
                    d_y = p_new.y;
                    x1 = border[j].x;
                    y1 = border[j].y;
                    x2 = border[(j + 1) % border.size()].x;
                    y2 = border[(j + 1) % border.size()].y;
                }
            }
        }
        else
        {
            for (size_t j = 0; j < holes.size(); j++)
            {
                if (points[i].is_in_polygon(holes[j]))
                {
                    d = std::numeric_limits<double>::infinity();
                    for (size_t k = 0; k < holes[j].size(); k++)
                    {
                        Point p_new = coord_point_to_line_segment(holes[j][k], holes[j][(k + 1) % holes[j].size()], points[i]);
                        double d_new = p_new.getDistance(points[i]);
                        if (d_new < d)
                        {
                            d = d_new;
                            d_x = p_new.x;
                            d_y = p_new.y;
                            x1 = holes[j][k].x;
                            y1 = holes[j][k].y;
                            x2 = holes[j][(k + 1) % holes[j].size()].x;
                            y2 = holes[j][(k + 1) % holes[j].size()].y;
                        }
                    }
                }
            }
        }
        if (d > d_star)
        {
            if (fabs(d_x - x1) < 1e-10 || fabs(d_x - x2) < 1e-10)
            {
                g[2 * i] += 2 * (x_i - d_x);
                g[2 * i + 1] += 2 * (y_i - d_y);
            }
            else
            {
                double A = y2 - y1;
                double B = x1 - x2;
                double C = x2 * y1 - x1 * y2;
                g[2 * i] += 2 * A * (A * x_i + B * y_i + C) / (A * A + B * B);
                g[2 * i + 1] += 2 * B * (A * x_i + B * y_i + C) / (A * A + B * B);
            }
        }
    }
    double dist, e;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < voronoi[i].size(); j++)
        {
            dist = points[i].getDistance(voronoi[i][j]);
            e = dist - x[n - 1];
            if (e > 0.0)
                g[n - 1] -= 2.0 * e;
        }
    }
    g[n - 1] += 2 * x[n - 1] / penalty_index;
}

// descent acc
double myvalue_acc(double *x, INT n)
{
    double f = 0.0;
    double dist, e;
    int N = n / 2; // the number of points

    std::vector<Point> points;
    for (size_t i = 0; i < N; i++)
        points.push_back(Point(x[2 * i], x[2 * i + 1]));
    VoronoiDiagram result = generate_diagram(points);
    std::vector<std::vector<Point>> voronoi = result.result();
    Clip res(voronoi, borders, holes);
    voronoi = res.clip_acc();

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < voronoi[i].size(); j++)
        {
            dist = points[i].getDistance(voronoi[i][j]);
            e = dist - R;
            if (e > 0)
                f += e * e;
        }
    }

    d_star = CONSTRAINT ? 1e-10 : R;
    double d;
    for (size_t i = 0; i < N; i++)
    {
        if (!points[i].is_in_polygon(border))
        {
            d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < border.size(); j++)
            {
                Point p_new = coord_point_to_line_segment(border[j], border[(j + 1) % border.size()], points[i]);
                double d_new = p_new.getDistance(points[i]);
                if (d_new < d)
                    d = d_new;
            }
            if (d > d_star)
                f += d * d;
        }
        else
        {
            for (size_t j = 0; j < holes.size(); j++)
            {
                if (points[i].is_in_polygon(holes[j]))
                {
                    d = std::numeric_limits<double>::infinity();
                    for (size_t k = 0; k < holes[j].size(); k++)
                    {
                        Point p_new = coord_point_to_line_segment(holes[j][k], holes[j][(k + 1) % holes[j].size()], points[i]);
                        double d_new = p_new.getDistance(points[i]);
                        if (d_new < d)
                            d = d_new;
                    }
                    if (d > d_star)
                        f += d * d;
                }
            }
        }
    }
    return f;
}

void mygrad_acc(double *g, double *x, INT n)
{
    double dist, e;
    double x_i, y_i, u_j, v_j;
    int N = n / 2;

    std::vector<Point> points;
    for (size_t i = 0; i < N; i++)
        points.push_back(Point(x[2 * i], x[2 * i + 1]));

    VoronoiDiagram result = generate_diagram(points);
    std::vector<std::vector<Point>> voronoi = result.result();
    Clip res(voronoi, borders, holes);
    voronoi = res.clip_acc();

    for (size_t i = 0; i < n; i++)
        g[i] = 0.0;

    for (size_t i = 0; i < N; i++)
    {
        x_i = points[i].x;
        y_i = points[i].y;

        for (size_t j = 0; j < voronoi[i].size(); j++)
        {
            u_j = voronoi[i][j].x;
            v_j = voronoi[i][j].y;
            dist = points[i].getDistance(voronoi[i][j]);
            e = dist - R;
            if (e > 0)
            {
                g[2 * i] += 2.0 * e * (x_i - u_j) / dist;
                g[2 * i + 1] += 2.0 * e * (y_i - v_j) / dist;
            }
        }
    }

    d_star = CONSTRAINT ? 1e-10 : R;
    double d, d_x, d_y, x1, y1, x2, y2;
    for (size_t i = 0; i < N; i++)
    {
        x_i = points[i].x;
        y_i = points[i].y;
        d = 0.0;
        if (!points[i].is_in_polygon(border))
        {
            d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < border.size(); j++)
            {
                Point p_new = coord_point_to_line_segment(border[j], border[(j + 1) % border.size()], points[i]);
                double d_new = p_new.getDistance(points[i]);
                if (d_new < d)
                {
                    d = d_new;
                    d_x = p_new.x;
                    d_y = p_new.y;
                    x1 = border[j].x;
                    y1 = border[j].y;
                    x2 = border[(j + 1) % border.size()].x;
                    y2 = border[(j + 1) % border.size()].y;
                }
            }
        }
        else
        {
            for (size_t j = 0; j < holes.size(); j++)
            {
                if (points[i].is_in_polygon(holes[j]))
                {
                    d = std::numeric_limits<double>::infinity();
                    for (size_t k = 0; k < holes[j].size(); k++)
                    {
                        Point p_new = coord_point_to_line_segment(holes[j][k], holes[j][(k + 1) % holes[j].size()], points[i]);
                        double d_new = p_new.getDistance(points[i]);
                        if (d_new < d)
                        {
                            d = d_new;
                            d_x = p_new.x;
                            d_y = p_new.y;
                            x1 = holes[j][k].x;
                            y1 = holes[j][k].y;
                            x2 = holes[j][(k + 1) % holes[j].size()].x;
                            y2 = holes[j][(k + 1) % holes[j].size()].y;
                        }
                    }
                }
            }
        }

        if (d > d_star)
        {
            if (fabs(d_x - x1) < 1e-10 || fabs(d_x - x2) < 1e-10)
            {
                g[2 * i] += 2 * (x_i - d_x);
                g[2 * i + 1] += 2 * (y_i - d_y);
            }
            else
            {
                double A = y2 - y1;
                double B = x1 - x2;
                double C = x2 * y1 - x1 * y2;
                g[2 * i] += 2 * A * (A * x_i + B * y_i + C) / (A * A + B * B);
                g[2 * i + 1] += 2 * B * (A * x_i + B * y_i + C) / (A * A + B * B);
            }
        }
    }
}

void mygrad_h_acc(double *g, double *x, INT n)
{
    double x_i, y_i;
    int N = n / 2;

    std::vector<Point> points;
    for (size_t i = 0; i < N; i++)
        points.push_back(Point(x[2 * i], x[2 * i + 1]));
    VoronoiDiagram result = generate_diagram(points);
    std::vector<std::vector<Point>> voronoi = result.result_with_neighbors();
    Clip res(voronoi, borders, holes);
    voronoi = res.clip_acc();

    for (size_t i = 0; i < n; i++)
        g[i] = 0.0;

    for (size_t i = 0; i < N; i++)
    {
        if (!voronoi[i].empty())
        {
            sub_points.clear();
            sub_points.push_back(points[i]);
            for (size_t j = 0; j < neighbors[i].size(); j++)
            {
                if (neighbors[i][j] >= 0)
                    sub_points.push_back(points[neighbors[i][j]]);
            }
            vertex_neighbors.clear();
            for (size_t j = 0; j < voronoi[i].size(); j++)
            {
                std::vector<int> neighbor;
                double d1 = points[i].getDistance(voronoi[i][j]);
                for (size_t k = 0; k < neighbors[i].size(); k++)
                {
                    if (neighbors[i][k] >= 0)
                    {
                        double d2 = points[neighbors[i][k]].getDistance(voronoi[i][j]);
                        if (fabs(d1 - d2) < 1e-9)
                            neighbor.push_back(neighbors[i][k]);
                    }
                }
                if (neighbor.size() > 2)
                {
                    std::vector<int> new_neighbor;
                    for (size_t k = 0; k < neighbor.size(); k++)
                    {
                        Point p = 0.5 * (points[i] + points[neighbor[k]]);
                        if (p.getDistance(voronoi[i][j]) > 1e-10)
                            new_neighbor.push_back(neighbor[k]);
                    }
                    neighbor = new_neighbor;
                }
                vertex_neighbors.push_back(neighbor);
            }
            g[2 * i] = computer_x_gradient(i, points, R, voronoi[i]);
            g[2 * i + 1] = computer_y_gradient(i, points, R, voronoi[i]);
        }
    }

    d_star = CONSTRAINT ? 1e-10 : R;
    double d, d_x, d_y, x1, y1, x2, y2;
    for (size_t i = 0; i < N; i++)
    {
        x_i = points[i].x;
        y_i = points[i].y;
        d = 0.0;
        if (!points[i].is_in_polygon(border))
        {
            d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < border.size(); j++)
            {
                Point p_new = coord_point_to_line_segment(border[j], border[(j + 1) % border.size()], points[i]);
                double d_new = p_new.getDistance(points[i]);
                if (d_new < d)
                {
                    d = d_new;
                    d_x = p_new.x;
                    d_y = p_new.y;
                    x1 = border[j].x;
                    y1 = border[j].y;
                    x2 = border[(j + 1) % border.size()].x;
                    y2 = border[(j + 1) % border.size()].y;
                }
            }
        }
        else
        {
            for (size_t j = 0; j < holes.size(); j++)
            {
                if (points[i].is_in_polygon(holes[j]))
                {
                    d = std::numeric_limits<double>::infinity();
                    for (size_t k = 0; k < holes[j].size(); k++)
                    {
                        Point p_new = coord_point_to_line_segment(holes[j][k], holes[j][(k + 1) % holes[j].size()], points[i]);
                        double d_new = p_new.getDistance(points[i]);
                        if (d_new < d)
                        {
                            d = d_new;
                            d_x = p_new.x;
                            d_y = p_new.y;
                            x1 = holes[j][k].x;
                            y1 = holes[j][k].y;
                            x2 = holes[j][(k + 1) % holes[j].size()].x;
                            y2 = holes[j][(k + 1) % holes[j].size()].y;
                        }
                    }
                }
            }
        }

        if (d > d_star)
        {
            if (fabs(d_x - x1) < 1e-10 || fabs(d_x - x2) < 1e-10)
            {
                g[2 * i] += 2 * (x_i - d_x);
                g[2 * i + 1] += 2 * (y_i - d_y);
            }
            else
            {
                double A = y2 - y1;
                double B = x1 - x2;
                double C = x2 * y1 - x1 * y2;
                g[2 * i] += 2 * A * (A * x_i + B * y_i + C) / (A * A + B * B);
                g[2 * i + 1] += 2 * B * (A * x_i + B * y_i + C) / (A * A + B * B);
            }
        }
    }
}

double myvalueR_acc(double *x, INT n)
{
    double f = 0.0;
    double dist, e;
    int N = n / 2; // the number of points
    std::vector<Point> points;
    for (size_t i = 0; i < N; i++)
        points.push_back(Point(x[2 * i], x[2 * i + 1]));
    VoronoiDiagram result = generate_diagram(points);
    std::vector<std::vector<Point>> voronoi = result.result();
    Clip res(voronoi, borders, holes);
    voronoi = res.clip_acc();
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < voronoi[i].size(); j++)
        {
            dist = points[i].getDistance(voronoi[i][j]);
            e = dist - x[n - 1];
            if (e > 0)
                f += e * e;
        }
    }
    d_star = CONSTRAINT ? 1e-10 : x[n - 1];
    double d;
    for (size_t i = 0; i < N; i++)
    {
        if (!points[i].is_in_polygon(border))
        {
            d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < border.size(); j++)
            {
                Point p_new = coord_point_to_line_segment(border[j], border[(j + 1) % border.size()], points[i]);
                double d_new = p_new.getDistance(points[i]);
                if (d_new < d)
                    d = d_new;
            }
            if (d > d_star)
                f += d * d;
        }
        else
        {
            for (size_t j = 0; j < holes.size(); j++)
            {
                if (points[i].is_in_polygon(holes[j]))
                {
                    d = std::numeric_limits<double>::infinity();
                    for (size_t k = 0; k < holes[j].size(); k++)
                    {
                        Point p_new = coord_point_to_line_segment(holes[j][k], holes[j][(k + 1) % holes[j].size()], points[i]);
                        double d_new = p_new.getDistance(points[i]);
                        if (d_new < d)
                            d = d_new;
                    }

                    if (d > d_star)
                        f += d * d;
                }
            }
        }
    }
    f += x[n - 1] * x[n - 1] / penalty_index;
    return f;
}

void mygradR_acc(double *g, double *x, INT n)
{
    double x_i, y_i;
    int N = n / 2;
    std::vector<Point> points;
    for (size_t i = 0; i < N; i++)
        points.push_back(Point(x[2 * i], x[2 * i + 1]));
    VoronoiDiagram result = generate_diagram(points);
    std::vector<std::vector<Point>> voronoi = result.result_with_neighbors();
    Clip res(voronoi, borders, holes);
    voronoi = res.clip_acc();
    for (size_t i = 0; i < n; i++)
        g[i] = 0.0;
    for (size_t i = 0; i < N; i++)
    {
        if (!voronoi[i].empty())
        {
            sub_points.clear();
            sub_points.push_back(points[i]);
            for (size_t j = 0; j < neighbors[i].size(); j++)
            {
                if (neighbors[i][j] >= 0)
                    sub_points.push_back(points[neighbors[i][j]]);
            }
            vertex_neighbors.clear();
            for (size_t j = 0; j < voronoi[i].size(); j++)
            {
                std::vector<int> neighbor;
                double d1 = points[i].getDistance(voronoi[i][j]);
                for (size_t k = 0; k < neighbors[i].size(); k++)
                {
                    if (neighbors[i][k] >= 0)
                    {
                        double d2 = points[neighbors[i][k]].getDistance(voronoi[i][j]);
                        if (fabs(d1 - d2) < 1e-9)
                            neighbor.push_back(neighbors[i][k]);
                    }
                }
                if (neighbor.size() > 2)
                {
                    std::vector<int> new_neighbor;
                    for (size_t k = 0; k < neighbor.size(); k++)
                    {
                        Point p = 0.5 * (points[i] + points[neighbor[k]]);
                        if (p.getDistance(voronoi[i][j]) > 1e-10)
                            new_neighbor.push_back(neighbor[k]);
                    }
                    neighbor = new_neighbor;
                }
                vertex_neighbors.push_back(neighbor);
            }
            g[2 * i] = computer_x_gradient(i, points, x[n - 1], voronoi[i]);
            g[2 * i + 1] = computer_y_gradient(i, points, x[n - 1], voronoi[i]);
        }
    }
    d_star = CONSTRAINT ? 1e-10 : x[n - 1];
    double d, d_x, d_y, x1, y1, x2, y2;
    for (size_t i = 0; i < N; i++)
    {
        x_i = points[i].x;
        y_i = points[i].y;
        d = 0.0;
        if (!points[i].is_in_polygon(border))
        {
            d = std::numeric_limits<double>::infinity();
            for (size_t j = 0; j < border.size(); j++)
            {
                Point p_new = coord_point_to_line_segment(border[j], border[(j + 1) % border.size()], points[i]);
                double d_new = p_new.getDistance(points[i]);
                if (d_new < d)
                {
                    d = d_new;
                    d_x = p_new.x;
                    d_y = p_new.y;
                    x1 = border[j].x;
                    y1 = border[j].y;
                    x2 = border[(j + 1) % border.size()].x;
                    y2 = border[(j + 1) % border.size()].y;
                }
            }
        }
        else
        {
            for (size_t j = 0; j < holes.size(); j++)
            {
                if (points[i].is_in_polygon(holes[j]))
                {
                    d = std::numeric_limits<double>::infinity();
                    for (size_t k = 0; k < holes[j].size(); k++)
                    {
                        Point p_new = coord_point_to_line_segment(holes[j][k], holes[j][(k + 1) % holes[j].size()], points[i]);
                        double d_new = p_new.getDistance(points[i]);
                        if (d_new < d)
                        {
                            d = d_new;
                            d_x = p_new.x;
                            d_y = p_new.y;
                            x1 = holes[j][k].x;
                            y1 = holes[j][k].y;
                            x2 = holes[j][(k + 1) % holes[j].size()].x;
                            y2 = holes[j][(k + 1) % holes[j].size()].y;
                        }
                    }
                }
            }
        }
        if (d > d_star)
        {
            if (fabs(d_x - x1) < 1e-10 || fabs(d_x - x2) < 1e-10)
            {
                g[2 * i] += 2 * (x_i - d_x);
                g[2 * i + 1] += 2 * (y_i - d_y);
            }
            else
            {
                double A = y2 - y1;
                double B = x1 - x2;
                double C = x2 * y1 - x1 * y2;
                g[2 * i] += 2 * A * (A * x_i + B * y_i + C) / (A * A + B * B);
                g[2 * i + 1] += 2 * B * (A * x_i + B * y_i + C) / (A * A + B * B);
            }
        }
    }
    double dist, e;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < voronoi[i].size(); j++)
        {
            dist = points[i].getDistance(voronoi[i][j]);
            e = dist - x[n - 1];
            if (e > 0.0)
                g[n - 1] -= 2.0 * e;
        }
    }
    g[n - 1] += 2 * x[n - 1] / penalty_index;
}

void update_borders(const std::vector<Point> &points, double alpha)
{
    std::vector<std::vector<Point>> sub_borders;
    for (size_t i = 0; i < points.size(); i++)
    {
        double x = points[i].x;
        double y = points[i].y;
        sub_borders.push_back({
            Point(x - alpha * R, y - alpha * R),
            Point(x + alpha * R, y - alpha * R),
            Point(x + alpha * R, y + alpha * R),
            Point(x - alpha * R, y + alpha * R),
        });
    }
    Clip res(sub_borders, border, holes);
    borders = res.clip_border();
}

// TPMBH
std::pair<std::vector<Point>, int> TPMBH(const std::vector<Point> &points, int maxRound)
{
    std::pair<std::vector<Point>, double> res;
    R *= 1 - 1e-3;
    // std::cout << "try to find better radius " << R << std::endl;
    res = MBH(points, maxRound);
    if (res.second < 1e-15)
        return std::make_pair(res.first, 1);
    R /= 1 - 1e-3;
    // std::cout << "try to find raw radius " << R << std::endl;
    res = MBH(res.first, maxRound);
    if (res.second < 1e-15)
        return std::make_pair(res.first, 1);

    return std::make_pair(res.first, 0);
}

std::pair<std::vector<Point>, double> MBH(std::vector<Point> points, int maxRound)
{
    uint64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0.5, 1.0);
    double alpha;
    double fx, best_fx;
    std::pair<std::vector<Point>, double> ERx;
    std::vector<Point> best_points;

    ERx = local_optimization(points);
    best_points = ERx.first;
    best_fx = ERx.second;
    // std::cout << "fx: " << best_fx << std::endl;
    if (best_fx < 1e-15)
        return std::make_pair(best_points, best_fx);

    int round = 0;
    while (round < maxRound)
    {
        // std::cout << "MBH search Round " << round << " : ";
        std::uniform_real_distribution<double> dis(0.0, 1.0);
        double choice = dis(generator);
        if (choice < 0.3)
            points = move_point_to_vertex(best_points, 0.05); // choose move point to vertex
        else
        {
            alpha = distribution(generator);
            points = perturb(best_points, R, alpha); // choose random perturbation
        }
        ERx = local_optimization(points);
        fx = ERx.second;
        // std::cout << "fx: " << fx << std::endl;
        if (fx < 1e-15)
            return std::make_pair(ERx.first, fx);
        if (fx < best_fx)
        {
            round = 0;
            best_points = ERx.first;
            best_fx = fx;
        }
        else
            round++;
    }

    return std::make_pair(best_points, best_fx);
}

std::pair<std::vector<Point>, double> local_optimization(std::vector<Point> points)
{
    INT n = points.size();
    double fx;
    Solution S;
    S.x = new double[2 * n];
    for (size_t i = 0; i < n; i++)
    {
        S.x[2 * i] = points[i].x;
        S.x[2 * i + 1] = points[i].y;
    }

    // stage 1
    CG_DESCENT::cg_descent(S.x, 2 * n, NULL, NULL, 1e-4, myvalue, mygrad, NULL, NULL);
    points.clear();
    for (size_t i = 0; i < n; i++)
        points.push_back(Point(S.x[2 * i], S.x[2 * i + 1]));

    // stage 2
    update_borders(points, 2.55);
    CG_DESCENT::cg_descent(S.x, 2 * n, NULL, NULL, 1e-10, myvalue_acc, mygrad_acc, NULL, NULL);
    points.clear();
    for (size_t i = 0; i < n; i++)
        points.push_back(Point(S.x[2 * i], S.x[2 * i + 1]));
    fx = myvalue(S.x, 2 * n);
    delete[] S.x;
    return std::make_pair(points, fx);
}

std::vector<Point> perturb(const std::vector<Point> &points, double radius, double perturb_index)
{
    std::vector<Point> new_points;
    uint64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    for (int i = 0; i < points.size(); i++)
    {
        Point point = points[i];
        std::uniform_real_distribution<double> rand_u(0.0, 1.0);
        std::uniform_real_distribution<double> rand_r(0.0, 1.0);
        double u = rand_u(generator);
        double r = rand_r(generator);
        double theta = 2.0 * M_PI * u;
        double rad = perturb_index * r * radius;
        point.x += rad * cos(theta);
        point.y += rad * sin(theta);
        new_points.push_back(point);
    }
    return new_points;
}

std::vector<Point> move_point_to_vertex(std::vector<Point> points, double alpha)
{
    int n = points.size();
    int t = n * alpha;
    double dist, e;
    std::vector<int> point_index;
    std::vector<double> point_loss;
    std::vector<Point> vertices;
    std::vector<double> vertex_loss;
    std::vector<int> good_points;
    std::vector<Point> bad_vertices;
    std::vector<std::vector<Point>> voronoi = generate_Voronoi(points);
    for (size_t i = 0; i < n; i++)
    {
        point_index.push_back(i);
        point_loss.push_back(0.0);
        for (size_t j = 0; j < voronoi[i].size(); j++)
        {
            dist = points[i].getDistance(voronoi[i][j]);
            e = dist - R;
            if (e > 0)
            {
                bool find = false;
                for (size_t k = 0; k < vertices.size(); k++)
                {
                    if (voronoi[i][j].getDistance(vertices[k]) < 1e-10)
                    {
                        find = true;
                        break;
                    }
                }
                if (!find)
                {
                    vertices.push_back(voronoi[i][j]);
                    vertex_loss.push_back(e);
                }
            }
            point_loss[i] += e;
        }
    }
    while (good_points.size() < 3 + t)
    {
        int best_point_index = std::min_element(point_loss.begin(), point_loss.end()) - point_loss.begin();
        good_points.push_back(point_index[best_point_index]);
        point_loss.erase(point_loss.begin() + best_point_index);
        point_index.erase(point_index.begin() + best_point_index);
    }
    while (bad_vertices.size() < 3 + t && !vertices.empty())
    {
        int worst_vertex_index = std::max_element(vertex_loss.begin(), vertex_loss.end()) - vertex_loss.begin();
        bad_vertices.push_back(vertices[worst_vertex_index]);
        vertices.erase(vertices.begin() + worst_vertex_index);
        vertex_loss.erase(vertex_loss.begin() + worst_vertex_index);
    }
    uint64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);
    std::uniform_int_distribution<> rand_point(0, good_points.size() - 1);
    std::uniform_int_distribution<> rand_vertex(0, bad_vertices.size() - 1);
    int a = rand_point(gen);
    int b = rand_vertex(gen);
    points[good_points[a]] = bad_vertices[b];
    return points;
}

// OERA
std::vector<Point> OERA(std::vector<Point> points)
{
    int n = points.size();
    std::pair<std::vector<Point>, double> Exr;
    std::pair<std::vector<Point>, int> TPMBH_result;
    double fx;
    int find = 0;
    do
    {
        /* SUMT */
        Exr = SUMT(points, R, n);
        std::pair<std::vector<Point>, double> final_result = local_search(Exr.first);
        if (final_result.second < R * (1 - 1e-6))
        {
            find = 1;
            points = final_result.first;
            R = final_result.second;
            // std::cout << "SUMT new R = " << R << std::endl;
        }
        else
        {
            if (final_result.second < R)
            {
                points = final_result.first;
                R = final_result.second;
            }
            // std::cout << "New R not found " << std::endl;
            break;
        }
        /* TPMHB_h */
        TPMBH_result = TPMBH_h(points, 3);
        if (TPMBH_result.second)
            points = TPMBH_result.first;
        else
            break;

    } while (find);

    return points;
}

// NOTE:
// For simple regions (e.g., C4), using a smaller start toler (1e-5) and alpha (5) is recommended to improve numerical accuracy.
// For more complex regions (e.g., C202H3), slightly larger values of start toler (1e-3) and alpha (10) are preferable for robustness and faster convergence.
std::pair<std::vector<Point>, double> SUMT(std::vector<Point> points, double r, int n)
{
    penalty_index = 1e2;
    double alpha = sqrt(10);
    double r_past = 0.0;
    double toler = 1e-5;
    double fx;

    Solution S;
    S.x = new double[2 * n + 1];
    for (size_t i = 0; i < n; i++)
    {
        S.x[2 * i] = points[i].x;
        S.x[2 * i + 1] = points[i].y;
    }
    S.x[2 * n] = r;

    CG_DESCENT::cg_descent(S.x, 2 * n + 1, NULL, NULL, toler, myvalueR, mygradR, NULL, NULL);
    r = S.x[2 * n];
    fx = penalty_index * myvalueR(S.x, 2 * n + 1) - S.x[2 * n] * S.x[2 * n];
    // std::cout << "pho: " << penalty_index << " fx: " << fx << " r: " << r << std::endl;
    penalty_index *= alpha;
    r_past = r;
    for (size_t i = 1; i < 20; i++)
    {
        toler = std::max(toler / 10, 1e-12);
        points.clear();
        for (size_t i = 0; i < n; i++)
            points.push_back(Point(S.x[2 * i], S.x[2 * i + 1]));
        update_borders(points, 1.55);
        CG_DESCENT::cg_descent(S.x, 2 * n + 1, NULL, NULL, toler, myvalueR_acc, mygradR_acc, NULL, NULL);
        r = S.x[2 * n];
        fx = penalty_index * myvalueR(S.x, 2 * n + 1) - S.x[2 * n] * S.x[2 * n];
        // std::cout << "pho: " << penalty_index << " fx: " << fx << " r: " << r << std::endl;
        penalty_index *= alpha;
        if ((fabs(r_past - r) < 1e-7 && fx < 1e-7) || fx < 1e-10)
            break;
        r_past = r;
    }
    points.clear();
    for (size_t i = 0; i < n; i++)
        points.push_back(Point(S.x[2 * i], S.x[2 * i + 1]));
    delete[] S.x;
    return std::make_pair(points, r);
}

std::pair<std::vector<Point>, int> TPMBH_h(const std::vector<Point> &points, int maxRound)
{
    std::pair<std::vector<Point>, double> res;
    R *= 1 - 1e-2;
    // std::cout << "try to find better radius " << R << std::endl;
    res = MBH_h(points, maxRound);
    if (res.second < 1e-15)
        return std::make_pair(res.first, 1);
    R /= 1 - 1e-2;
    // std::cout << "try to find raw radius " << R << std::endl;
    res = MBH_h(res.first, maxRound);
    if (res.second < 1e-15)
        return std::make_pair(res.first, 1);

    return std::make_pair(res.first, 0);
}

std::pair<std::vector<Point>, double> MBH_h(std::vector<Point> points, int maxRound)
{
    uint64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution(0.05, 0.35);
    double alpha;
    double fx, best_fx;
    std::pair<std::vector<Point>, double> ERx;
    std::vector<Point> best_points;
    ERx = local_optimization_h(points);
    best_points = ERx.first;
    best_fx = ERx.second;
    // std::cout << "fx: " << best_fx << std::endl;
    if (best_fx < 1e-15)
        return std::make_pair(best_points, best_fx);
    int round = 0;
    while (round < maxRound)
    {
        // std::cout << "MBH_h search Round " << round << " : ";
        alpha = distribution(generator);
        points = perturb(best_points, R, alpha);
        ERx = local_optimization_h(points);
        fx = ERx.second;
        // std::cout << "fx: " << fx << std::endl;
        if (fx < 1e-15)
            return std::make_pair(ERx.first, fx);
        if (fx < best_fx)
        {
            round = 0;
            best_points = ERx.first;
            best_fx = fx;
        }
        else
            round++;
    }

    return std::make_pair(best_points, best_fx);
}

std::pair<std::vector<Point>, double> local_optimization_h(std::vector<Point> points)
{
    INT n = points.size();
    double fx;
    std::vector<Point> best_points;
    double best_fx;
    Solution S;
    S.x = new double[2 * n];
    for (size_t i = 0; i < n; i++)
    {
        S.x[2 * i] = points[i].x;
        S.x[2 * i + 1] = points[i].y;
    }
    update_borders(points, 2.05);
    CG_DESCENT::cg_descent(S.x, 2 * n, NULL, NULL, 1e-10, myvalue_acc, mygrad_acc, NULL, NULL);
    fx = myvalue(S.x, 2 * n);
    best_points.clear();
    for (size_t i = 0; i < n; i++)
        best_points.push_back(Point(S.x[2 * i], S.x[2 * i + 1]));
    best_fx = fx;
    if (fx < 1e-15)
    {
        // std::cout << "mygrad found" << std::endl;
        delete[] S.x;
        return std::make_pair(best_points, fx);
    }

    points.clear();
    for (size_t i = 0; i < n; i++)
        points.push_back(Point(S.x[2 * i], S.x[2 * i + 1]));
    update_borders(points, 2.05);
    CG_DESCENT::cg_descent(S.x, 2 * n, NULL, NULL, 1e-10, myvalue_acc, mygrad_h_acc, NULL, NULL);
    fx = myvalue(S.x, 2 * n);
    if (fx < 1e-15)
    {
        // std::cout << "mygrad + mygrad_h found" << std::endl;
        best_points.clear();
        for (size_t i = 0; i < n; i++)
            best_points.push_back(Point(S.x[2 * i], S.x[2 * i + 1]));
        delete[] S.x;
        return std::make_pair(best_points, fx);
    }
    if (fx < best_fx)
    {
        best_fx = fx;
        best_points.clear();
        for (size_t i = 0; i < n; i++)
            best_points.push_back(Point(S.x[2 * i], S.x[2 * i + 1]));
    }
    delete[] S.x;
    return std::make_pair(best_points, best_fx);
}

// lattice-based generate points
std::pair<std::vector<Point>, double> initial_solution(int n, int maxRound)
{
    std::vector<Point> best_points;
    double R_min;
    std::vector<Point> points = generate_points(n);
    std::pair<std::vector<Point>, double> final_result = local_search(points);
    best_points = final_result.first;
    R = final_result.second;
    R_min = R;
    std::cout << "initial R = " << R << std::endl;
    int round = 0;
    while (round < maxRound)
    {
        points = generate_lattice_points(n, R);
        Solution S;
        S.x = new double[2 * n + 1];
        for (size_t i = 0; i < n; i++)
        {
            S.x[2 * i] = points[i].x;
            S.x[2 * i + 1] = points[i].y;
        }
        S.x[2 * n] = R;
        penalty_index = 10;
        CG_DESCENT::cg_descent(S.x, 2 * n + 1, NULL, NULL, 1e-2, myvalueR, mygradR, NULL, NULL);
        double fx = penalty_index * myvalueR(S.x, 2 * n + 1) - S.x[2 * n] * S.x[2 * n];
        // std::cout << "pho: " << penalty_index << " fx: " << fx << " r: " << S.x[2 * n] << std::endl;
        if (S.x[2 * n] < R_min)
        {
            R_min = S.x[2 * n];
            best_points = points;
        }
        delete[] S.x;
        round++;
    }
    std::pair<std::vector<Point>, double> res = SUMT(best_points, R, n);
    return res;
}

std::vector<Point> generate_lattice_points(int n, double r)
{
    std::vector<Point> points;
    int d = 0;
    while (points.size() != n)
    {
        if (d > 100)
        {
            r = 0.99 * r;
            d = 0;
        }
        points = find_lattice_points(n, r);
        d++;
    }
    return points;
}

std::vector<Point> find_lattice_points(int n, double r)
{
    std::vector<Point> points;
    Point p0 = generate_points(1)[0];
    points.push_back(p0);
    int number = 1;
    int I = 1;
    Point vr(3, sqrt(3));
    Point wr(3, -sqrt(3));
    uint64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> rand_u(0.0, 1.0);
    double u = rand_u(generator);
    double theta = 2.0 * M_PI * u;
    int k, l;
    Point p_new;
    double d, x, y;
    std::vector<Point> sub_points;
    std::vector<double> sub_ds;
    double t = holes.empty() ? 0.7 : 0.5;
    while (number < n)
    {
        int find_new_point = 0;
        for (size_t i = 0; i < 2 * I + 1; i++)
        {
            k = i - I;
            if (fabs(k) < I)
            {
                l = -I;
                p_new = k * r / 2 * vr + l * r / 2 * wr;
                x = p_new.x;
                y = p_new.y;
                p_new = p0 + Point(x * cos(theta) - y * sin(theta), y * cos(theta) + x * sin(theta));
                d = find_minimun_distance(p_new);
                if (fabs(d) < 1e-10)
                {
                    points.push_back(p_new);
                    number++;
                    find_new_point++;
                    if (number == n)
                        break;
                }
                else if (d < t * r)
                {
                    sub_points.push_back(p_new);
                    sub_ds.push_back(d);
                    find_new_point++;
                }
                l = I;
                p_new = k * r / 2 * vr + l * r / 2 * wr;
                x = p_new.x;
                y = p_new.y;
                p_new = p0 + Point(x * cos(theta) - y * sin(theta), y * cos(theta) + x * sin(theta));
                d = find_minimun_distance(p_new);
                if (fabs(d) < 1e-10)
                {
                    points.push_back(p_new);
                    number++;
                    find_new_point++;
                    if (number == n)
                        break;
                }
                else if (d < t * r)
                {
                    sub_points.push_back(p_new);
                    sub_ds.push_back(d);
                    find_new_point++;
                }
            }
            else
            {
                for (size_t j = 0; j < 2 * I + 1; j++)
                {
                    l = j - I;
                    p_new = k * r / 2 * vr + l * r / 2 * wr;
                    x = p_new.x;
                    y = p_new.y;
                    p_new = p0 + Point(x * cos(theta) - y * sin(theta), y * cos(theta) + x * sin(theta));
                    d = find_minimun_distance(p_new);
                    if (fabs(d) < 1e-10)
                    {
                        points.push_back(p_new);
                        number++;
                        find_new_point++;
                        if (number == n)
                            break;
                    }
                    else if (d < t * r)
                    {
                        sub_points.push_back(p_new);
                        sub_ds.push_back(d);
                        find_new_point++;
                    }
                }
                if (number == n)
                    break;
            }
        }
        if (!find_new_point)
            break;
        I++;
    }
    while (number < n && !sub_points.empty())
    {
        int d_min = std::min_element(sub_ds.begin(), sub_ds.end()) - sub_ds.begin();
        points.push_back(sub_points[d_min]);
        number++;
        sub_points.erase(sub_points.begin() + d_min);
        sub_ds.erase(sub_ds.begin() + d_min);
    }
    return points;
}

double find_minimun_distance(Point p)
{
    double d = 0.0;
    if (!p.is_in_polygon(border))
    {
        d = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < border.size(); i++)
        {
            Point p_new = coord_point_to_line_segment(border[i], border[(i + 1) % border.size()], p);
            double d_new = p_new.getDistance(p);
            if (d_new < d)
                d = d_new;
        }
    }
    else
    {
        for (size_t i = 0; i < holes.size(); i++)
        {
            if (p.is_in_polygon(holes[i]))
            {
                d = std::numeric_limits<double>::infinity();
                for (size_t j = 0; j < holes[i].size(); j++)
                {
                    Point p_new = coord_point_to_line_segment(holes[i][j], holes[i][(j + 1) % holes[i].size()], p);
                    double d_new = p_new.getDistance(p);
                    if (d_new < d)
                        d = d_new;
                }
                break;
            }
        }
    }
    return d;
}

//
extern std::string Region;
std::vector<Point> ALGORITHM(int n, int run)
{
    std::string filename = Region + "_" + std::to_string(n) + "_" + std::to_string(run) + ".txt";
    std::string path = "../results/" + Region + "/" + std::to_string(n) + "/" + filename; // please change the path according to your own directory structure
    std::vector<Point> best_points;
    double R_old;

    /* First Stage of the Search */
    std::pair<std::vector<Point>, double> IS = initial_solution(n, 10);
    std::vector<Point> points = IS.first;
    R = IS.second;
    std::cout << "Raw R: " << R << std::endl;
    best_points = points;
    // print(best_points);
    write(best_points, R, path);

    /* Second Stage of the Search */
    std::pair<std::vector<Point>, int> ERx;
    while (!time_out())
    {
        std::cout << "find R: " << R << " --------------------------------------------------" << std::endl;
        std::vector<Point> random_points = local_search(generate_points(n)).first;
        ERx = TPMBH(random_points, 10);
        if (ERx.second)
        {
            points = ERx.first;
            // std::cout << "ER(x) = 0 found, now the radius is " << R << std::endl;
            best_points = points;
            log(n, run, R, print_time());
            write(best_points, R, path);
            R_old = R;
            points = OERA(best_points);
            if (R < R_old * (1 - 1e-4))
            {
                // std::cout << "Adjust radius done, new R = " << R << std::endl;
                std::pair<std::vector<Point>, double> final_result = local_search(points);
                best_points = final_result.first;
                R = final_result.second;
                // std::cout << "after fix R: " << R << std::endl;
                log(n, run, R, print_time());
                write(best_points, R, path);
            }
        }
        // else
        //     std::cout << "ER(x) = 0 not found " << std::endl;
    }

    // std::cout << "final R: " << R << std::endl;
    // if (result_is_legal(best_points, R))
    //     std::cout << "Legal results. " << std::endl;
    // else
    //     std::cout << "Illegal results. " << std::endl;
    // rename
    rename(path, R);
    return best_points;
}

void check_result(std::vector<Point> points)
{
    double f = 0.0;
    double dist, e;
    std::cout << std::setprecision(10);
    // check the coverage violations
    std::vector<std::vector<Point>> voronoi = generate_Voronoi(points);

    std::cout << "Check the coverage violations:" << std::endl;
    for (size_t i = 0; i < points.size(); i++)
    {
        for (size_t j = 0; j < voronoi[i].size(); j++)
        {
            dist = points[i].getDistance(voronoi[i][j]);
            e = dist - R;
            if (e > 0)
                f += e * e;
        }
    }
    std::cout << "Coverage violation (sum of squared errors): " << f << std::endl;

    f = 0.0;
    if (CONSTRAINT)
    {
        // check the constraint violations
        std::cout << "Check the constraint violations:" << std::endl;
        double d;
        for (size_t i = 0; i < points.size(); i++)
        {
            if (!points[i].is_in_polygon(border))
            {
                d = std::numeric_limits<double>::infinity();
                for (size_t j = 0; j < border.size(); j++)
                {
                    Point p_new = coord_point_to_line_segment(border[j], border[(j + 1) % border.size()], points[i]);
                    double d_new = p_new.getDistance(points[i]);
                    if (d_new < d)
                        d = d_new;
                }
                f += d * d;
            }
            else
            {
                for (size_t j = 0; j < holes.size(); j++)
                {
                    if (points[i].is_in_polygon(holes[j]))
                    {
                        d = std::numeric_limits<double>::infinity();
                        for (size_t k = 0; k < holes[j].size(); k++)
                        {
                            Point p_new = coord_point_to_line_segment(holes[j][k], holes[j][(k + 1) % holes[j].size()], points[i]);
                            double d_new = p_new.getDistance(points[i]);
                            if (d_new < d)
                                d = d_new;
                        }
                        f += d * d;
                    }
                }
            }
        }
        std::cout << "Constraint violation (sum of squared errors): " << f << std::endl;
    }
}