#include "LocalSearch.h"

// GLOBAL
double R;
std::vector<Point> border;
std::vector<std::vector<Point>> holes;
std::vector<std::vector<std::vector<Point>>> borders;
std::string Region = "C4";
bool CONSTRAINT = 1; // 0: unconstrained CRICMR , 1: constrained CRICMR

// time
std::chrono::_V2::steady_clock::time_point TIME_START;
double TIME;
double TIME_MAX = 60.0 * 60.0 * 1000.0; // min/s/ms

// save_Voronoi(points, voronoi, "../figure");
// checkResult(points, "../figure");
// save_result(points, R, "../figure");

int main()
{
    std::cout << std::setprecision(20);

    std::string inputRegion;
    int n = 0;
    int constraintInt = -1;

    /* ------ Input Region ------ */
    std::cout << "Please input Region (e.g. Square, C14H2): ";
    if (!(std::cin >> inputRegion))
    {
        std::cerr << "Error: failed to read Region." << std::endl;
        return 1;
    }

    /* ------ Input n ------ */
    std::cout << "Please input n (number of circles): ";
    if (!(std::cin >> n) || n <= 0)
    {
        std::cerr << "Error: n must be positive." << std::endl;
        return 1;
    }

    /* ------ Input CONSTRAINT ------ */
    std::cout << "Please input CONSTRAINT (0 = unconstrained, 1 = constrained): ";
    if (!(std::cin >> constraintInt) || (constraintInt != 0 && constraintInt != 1))
    {
        std::cerr << "Error: CONSTRAINT must be 0 or 1." << std::endl;
        return 1;
    }

    Region = inputRegion;
    CONSTRAINT = (constraintInt != 0);

    input_files(border, holes, Region);
    //
    std::vector<Point> points_best;
    std::vector<Point> points;
    double R_best = 1.0e8;
    double R_avg = 0.0;
    double R_worst = 0.0;
    std::vector<double> R_list;
    double sigma = 0.0;
    double Time_avg = 0.0;

    TIME_START = std::chrono::steady_clock::now();
    int runs = 10;
    for (size_t i = 0; i < runs; i++)
    {
        int run = i;
        std::cout << "RUN " << run << std::endl;
        R = 0.0;
        TIME_START = std::chrono::steady_clock::now();
        points = ALGORITHM(n, run);
        TIME = print_time();
        log(n, run, R, TIME);
        R_list.push_back(R);
        R_avg += R;
        if (R < R_best - 1.0e-8)
        {
            R_best = R;
            points_best = points;
        }
        if (R > R_worst)
            R_worst = R;

        Time_avg += TIME;
        //
    }

    R_avg /= runs;
    Time_avg /= runs;
    sigma = deviation(R_list, runs);
    output_details(points_best, R_best);
    output_summary(n, R_best, R_avg, R_worst, sigma, Time_avg);
}

// int main(int argc, char *argv[])
// {
//     Region = argv[1];
//     int n = atoi(argv[2]);
//     CONSTRAINT = atoi(argv[3]);
//     std::cout << std::setprecision(20);
//     input_files(border, holes, Region);
//     //
//     std::vector<Point> points_best;
//     std::vector<Point> points;
//     double R_best = 1.0e8;
//     double R_avg = 0.0;
//     double R_worst = 0.0;
//     std::vector<double> R_list;
//     double sigma = 0.0;
//     double Time_avg = 0.0;

//     TIME_START = std::chrono::steady_clock::now();
//     int runs = 10;
//     for (size_t i = 0; i < runs; i++)
//     {
//         int run = i;
//         std::cout << "RUN " << run << std::endl;
//         R = 0.0;
//         TIME_START = std::chrono::steady_clock::now();
//         points = ALGORITHM(n, run);
//         TIME = print_time();
//         log(n, run, R, TIME);
//         R_list.push_back(R);
//         R_avg += R;
//         if (R < R_best - 1.0e-8)
//         {
//             R_best = R;
//             points_best = points;
//         }
//         if (R > R_worst)
//             R_worst = R;

//         Time_avg += TIME;
//         //
//     }

//     R_avg /= runs;
//     Time_avg /= runs;
//     sigma = deviation(R_list, runs);
//     output_details(points_best, R_best);
//     output_summary(n, R_best, R_avg, R_worst, sigma, Time_avg);
// }

// result checking
// int main()
// {
//     std::cout << std::setprecision(20);
//     input_files(border, holes, Region);
//     int n = 100;
//     //
//     std::pair<std::vector<Point>, double> result = read_points("../results/100.txt");
//     print(result.first);
//     R = result.second;
//     std::cout << "R = " << R << std::endl;
//     std::cout << "is legal: " << result_is_legal(result.first, R) << std::endl;
//     check_result(result.first);
// }