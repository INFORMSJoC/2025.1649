#include "Base.h"

// GLOBAL
extern bool CONSTRAINT;
extern std::string Region;
extern std::vector<Point> border;
extern std::vector<std::vector<Point>> holes;
extern std::chrono::_V2::steady_clock::time_point TIME_START;
extern double TIME_MAX;

// LOCAL
double x_min = std::numeric_limits<double>::infinity();
double x_max = -std::numeric_limits<double>::infinity();
double y_min = std::numeric_limits<double>::infinity();
double y_max = -std::numeric_limits<double>::infinity();

// Initialization
void input_files(std::vector<Point> &border, std::vector<std::vector<Point>> &holes, std::string region)
{
    std::ifstream infile;
    infile.open("../data/" + region + ".txt", std::ios::in);
    if (!infile.is_open())
    {
        std::cout << "ERROR" << std::endl;
    }
    char dataType;
    int dataSize;

    while (infile >> dataType >> dataSize)
    {
        if (dataType == 'C')
        {
            for (size_t i = 0; i < dataSize; i++)
            {
                double tmp_x;
                double tmp_y;
                infile >> tmp_x >> tmp_y;
                border.push_back(Point(tmp_x, tmp_y));
            }
        }
        else if (dataType == 'H')
        {
            std::vector<Point> tmp_hole;
            for (size_t i = 0; i < dataSize; i++)
            {
                double tmp_x;
                double tmp_y;
                infile >> tmp_x >> tmp_y;
                tmp_hole.push_back(Point(tmp_x, tmp_y));
            }
            holes.push_back(tmp_hole);
        }
    }

    // reset the range
    for (size_t i = 0; i < border.size(); i++)
    {
        if (border[i].x > x_max)
            x_max = border[i].x;
        if (border[i].x < x_min)
            x_min = border[i].x;
        if (border[i].y > y_max)
            y_max = border[i].y;
        if (border[i].y < y_min)
            y_min = border[i].y;
    }
}

std::vector<Point> generate_points(int n)
{
    uint64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution_x(x_min, x_max);
    std::uniform_real_distribution<double> distribution_y(y_min, y_max);
    std::vector<Point> points;
    for (size_t i = 0; i < n; i++)
    {
        Point point = {distribution_x(generator), distribution_y(generator)};
        bool in_hole = false;
        for (size_t j = 0; j < holes.size(); j++)
        {
            if (point.is_in_polygon(holes[j]))
            {
                in_hole = true;
                break;
            }
        }

        if (point.is_in_polygon(border) && !in_hole)
            points.push_back(point);
        else
            i--;
    }
    return points;
}

std::pair<std::vector<Point>, double> read_points(std::string path)
{
    std::vector<Point> points;
    std::ifstream infile(path);
    if (!infile.is_open())
    {
        std::cout << "ERROR, no file found: " << path << std::endl;
        return std::make_pair(points, -1.0);
    }
    int n = 0;
    double r = 0.0;
    infile >> n >> r;
    double x, y;
    while (infile >> x >> y)
    {
        points.emplace_back(x, y);
    }
    return std::make_pair(points, r);
}

// Print
void print(std::vector<Point> points)
{
    for (size_t i = 0; i < points.size(); i++)
        std::cout << "Point" << points[i] << "," << std::endl;
}

double print_time()
{
    double a = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - TIME_START).count() / 1000.0;
    if (a < 1.0)
        std::cout << a * 1000.0 << " ms" << std::endl;
    else if (a < 60.0)
        std::cout << a << " s" << std::endl;
    else if (a < (60.0 * 60.0))
        std::cout << a / 60.0 << " min" << std::endl;
    else
        std::cout << a / (60.0 * 60.0) << " h" << std::endl;
    return a;
}

// Result
void write(std::vector<Point> points, double r, std::string path)
{
    std::ofstream ofs;
    ofs.open(path, std::ios::out);
    ofs << std::setprecision(15);
    ofs << "Time: " << print_time() << std::endl;
    ofs << "R: " << r << std::endl;
    for (size_t i = 0; i < points.size(); i++)
        ofs << points[i].x << " " << points[i].y << std::endl;
    ofs.close();
}

void rename(std::string path, double r)
{
    std::string file_old = path;
    std::string file_new = path.substr(0, path.size() - 4) + "_" + std::to_string(r) + ".txt";
    std::rename(file_old.c_str(), file_new.c_str());
}

void log(int n, int run, double r, double time)
{
    std::ofstream ofs;
    ofs.open("../results/" + Region + "/" + std::to_string(n) + "/log.txt", std::ios::app);
    ofs << std::setprecision(15);
    ofs << run << " " << r << " " << time << std::endl;
    ofs.close();
}

void output_details(std::vector<Point> points, double r)
{
    std::ofstream ofs;
    ofs.open("../results/" + Region + "/" + std::to_string(points.size()) + ".txt", std::ios::out);
    ofs << std::setprecision(15);
    ofs << points.size() << "   " << r << std::endl;
    for (size_t i = 0; i < points.size(); i++)
        ofs << points[i].x << " " << points[i].y << std::endl;
    ofs.close();
}

void output_summary(int n, double r_best, double r_avg, double r_worst, double sigma, double time_avg)
{
    std::ofstream ofs;
    ofs.open("../results/" + Region + "/ResultSummary.txt", std::ios::app);
    ofs << std::setprecision(10);
    ofs << n << " " << r_best << " " << r_avg << " " << r_worst << " " << sigma << " " << time_avg << std::endl;
    ofs.close();
}

double deviation(std::vector<double> a, int n)
{
    double sum = 0, tmp = 0, x_avg;
    for (size_t i = 0; i < n; i++)
        sum += a[i];
    x_avg = sum / n;
    for (size_t i = 0; i < n; i++)
        tmp += (a[i] - x_avg) * (a[i] - x_avg);
    return sqrt(tmp / n);
} // standard deviation

// Check
bool result_is_legal(std::vector<Point> points, double r)
{
    std::vector<std::vector<Point>> voronoi = generate_Voronoi(points);
    for (size_t i = 0; i < voronoi.size(); i++)
    {
        for (size_t j = 0; j < voronoi[i].size(); j++)
        {
            if (points[i].getDistance(voronoi[i][j]) > r + 1e-8)
                return false;
        }
    }
    return true;
}

bool time_out()
{
    double current_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - TIME_START).count();
    if (current_time > TIME_MAX - 60000)
        return true;
    else
        return false;
}

// Voronoi
VoronoiDiagram generate_diagram(std::vector<Point> points)
{
    // Construct diagram
    FortuneAlgorithm algorithm(points);
    algorithm.construct();
    // Bound the diagram
    algorithm.bound(Box{x_min, y_min, x_max, y_max});
    // Intersect the diagram with a box
    VoronoiDiagram diagram = algorithm.getDiagram();
    return diagram;
}

std::vector<std::vector<Point>> generate_Voronoi(std::vector<Point> points)
{
    VoronoiDiagram result = generate_diagram(points);
    std::vector<std::vector<Point>> voronoi = result.result();
    Clip res(voronoi, border, holes);
    voronoi = res.clip();
    return voronoi;
}

// local search
bool is_move(std::vector<Point> old_points, std::vector<Point> new_points, double r)
{
    double L = 0.0;
    for (size_t i = 0; i < old_points.size(); i++)
    {
        L += old_points[i].getDistance(new_points[i]);
    }
    if (L < r * 1e-7)
        return false;
    else
        return true;
}

std::pair<std::vector<Point>, double> local_search(std::vector<Point> points)
{
    int n = points.size();
    std::vector<Point> old_points;
    std::vector<Point> new_points;
    std::vector<Point> min_points;
    std::vector<double> old_radii;
    std::vector<double> new_radii;
    double old_radius;
    double new_radius;
    double min_radius;

    // initial configuration
    new_points = points;
    min_points = points;
    min_radius = std::numeric_limits<double>::infinity();
    int round = 1;
    bool still_move = true;

    // start
    while (still_move && round < 1000)
    {
        old_points = new_points;
        new_points.clear();
        new_radii.clear();
        // voronoi
        std::vector<std::vector<Point>> voronoi = generate_Voronoi(old_points);

        // If any Voronoi cell is empty, regenerate a new input set
        bool voronoi_empty = false;
        for (size_t i = 0; i < voronoi.size(); ++i)
        {
            if (voronoi[i].empty())
            {
                voronoi_empty = true;
                break;
            }
        }
        if (voronoi_empty)
        {
            std::cout << "Empty Voronoi cell found - regenerating initial points." << std::endl;
            new_points = generate_points(n);
            min_points = new_points;
            min_radius = std::numeric_limits<double>::infinity();
            still_move = true;
            round = 1;
            continue;
        }

        if (CONSTRAINT)
        {
            for (size_t i = 0; i < voronoi.size(); i++)
            {
                CMCC cmb(voronoi[i], border, holes);
                new_points.push_back(cmb.center);
                new_radii.push_back(cmb.radius);
            }
        }
        else
        {
            for (size_t i = 0; i < voronoi.size(); i++)
            {
                MCC mb(voronoi[i]);
                new_points.push_back(mb.center);
                new_radii.push_back(mb.radius);
            }
        }

        new_radius = *std::max_element(new_radii.begin(), new_radii.end());
        if (new_radius < min_radius)
        {
            min_points = new_points;
            min_radius = new_radius;
        }
        still_move = is_move(old_points, new_points, min_radius);
        round++;
    }
    // std::cout << "local search done. (" << round << ")" << std::endl;
    // check_result(min_points);
    return std::make_pair(min_points, min_radius);
}
