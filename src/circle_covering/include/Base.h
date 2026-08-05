#pragma once
// STL
#include <sstream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <iomanip>
// My includes
#include "Voronoi.h"

// Initialization
void input_files(std::vector<Point> &border, std::vector<std::vector<Point>> &holes, std::string region);
std::vector<Point> generate_points(int n);
std::pair<std::vector<Point>, double> read_points(std::string path);

// Print
void print(std::vector<Point> points);
double print_time();

// Result
void write(std::vector<Point> points, double r, std::string path);
void rename(std::string path, double r);
void log(int n, int run, double r, double time);
void output_details(std::vector<Point> points, double r);
void output_summary(int n, double r_best, double r_avg, double r_worst, double sigma, double time_avg);
double deviation(std::vector<double> a, int n);

// Draw
void save_region(std::string path);
void save_points(std::vector<Point> points, std::string name, std::string path);
void save_result(std::vector<Point> points, double radius, std::string path);
void save_Voronoi(std::vector<Point> points, std::vector<std::vector<Point>> voronoi, std::string path);

// Check
bool result_is_legal(std::vector<Point> points, double r);
bool time_out();

// Voronoi
VoronoiDiagram generate_diagram(std::vector<Point> points);
std::vector<std::vector<Point>> generate_Voronoi(std::vector<Point> points);

// Local search
bool is_move(std::vector<Point> old_points, std::vector<Point> new_points, double r);
std::pair<std::vector<Point>, double> local_search(std::vector<Point> points);