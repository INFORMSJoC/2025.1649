#pragma once
// My includes
#include "Base.h"
#include "cg_blas.h"
#include "cg_user.h"
#include "cg_descent.h"

typedef struct Configuration
{
    double f;
    double *x;
} Solution; //  Configuration X

// gradient
Point coord_point_to_line_segment(Point a, Point b, Point point);
std::pair<Point, Point> find_line(Point p);
std::pair<double, double> duv_x(const std::vector<Point>& points, int n, int t, Point vertex);
std::pair<double, double> duv_y(const std::vector<Point>& points, int n, int t, Point vertex);
double computer_x_gradient(int n, const std::vector<Point>& points, double r, const std::vector<Point>& voronoi);
double computer_y_gradient(int n, const std::vector<Point>& points, double r, const std::vector<Point>& voronoi);

// descent
double myvalue(double *x, INT n);
void mygrad(double *g, double *x, INT n);
void mygrad_h(double *g, double *x, INT n);
double myvalueR(double *x, INT n);
void mygradR(double *g, double *x, INT n);

// descent acc
double myvalue_acc(double *x, INT n);
void mygrad_acc(double *g, double *x, INT n);
void mygrad_h_acc(double *g, double *x, INT n);
double myvalueR_acc(double *x, INT n);
void mygradR_acc(double *g, double *x, INT n);
void update_borders(const std::vector<Point> &points, double alpha);

// fixed radius
std::pair<std::vector<Point>, int> TPMBH(const std::vector<Point>& points, int maxRound);
std::pair<std::vector<Point>, double> MBH(std::vector<Point> points, int maxRound);
std::pair<std::vector<Point>, double> local_optimization(std::vector<Point> points);
std::vector<Point> perturb(const std::vector<Point>& points, double radius, double perturb_index);
std::vector<Point> move_point_to_vertex(std::vector<Point> points, double alpha);

// adjust radius
std::vector<Point> OERA(std::vector<Point> points);
std::pair<std::vector<Point>, double> SUMT(std::vector<Point> points, double r, int n);
std::pair<std::vector<Point>, int> TPMBH_h(const std::vector<Point>& points, int maxRound);
std::pair<std::vector<Point>, double> MBH_h(std::vector<Point> points, int maxRound);
std::pair<std::vector<Point>, double> local_optimization_h(std::vector<Point> points);

// lattice-based generate points
std::pair<std::vector<Point>, double> initial_solution(int n, int maxRound);
std::vector<Point> generate_lattice_points(int n, double r);
std::vector<Point> find_lattice_points(int n, double r);
double find_minimun_distance(Point p);

//
std::vector<Point> ALGORITHM(int n, int run);
void check_result(std::vector<Point> points);