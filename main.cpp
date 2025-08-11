#include <iostream>
#include <cmath>    // For std::log
#include <iomanip>  // for formatting string stream
#include <vector>
#include <stdexcept>
#include <limits>   // For std::numeric_limits
#include <optional> // For std::optional


double getKConst(double CL, double V) {
    return CL/V;
}

double get_half_life(double ke){
    return std::log(2.0) / ke; // t1/2 = ln(2) / ke
}

double conc_single_iv_bolus(double dose, double V, double ke, double t) {
    return (dose / V) * std::exp(-ke * t); // C(t) = (dose/V) * e^(-ke *t)
}

double ask_dose() {
    // std::cin >> std::ws; // This causes issue of asking input before showing prompt.
    
    double dose;
    std::cout << "\nEnter dose in mg: " << std::flush;
    std::cin >> dose; 

    if (std::cin.fail() || dose <= 0) {
        std::cout << "Invalid dose. Please enter a positive number." << std::endl;
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        return ask_dose();
    }
    std::cout << std::endl;

    return dose;
}

std::vector<double> build_time_vectors(
    double t_end,
    std::optional<double> dt
) {
    double step = dt.value_or(t_end / 200.0);

    if (t_end < 0.0)  throw std::invalid_argument("t_end must be >= 0.");
    if (step <= 0.0)    throw std::invalid_argument("dt must be > 0.");
    
    const std::size_t n = static_cast<std::size_t>(std::floor(t_end / step) + 1.0); // how many multiples of dt between 0 and t_end 

    std::vector<double> times; // Creates an empty dynamic array.
    times.reserve(n); // preallocate size.

    for (std::size_t i = 0; i < n; ++i) {
        times.push_back(static_cast<double>(i) * step);
    }
    return times;
}

int main() {
    const double CL = 3.0;
    const double V = 30.0;

    std::cout << "PKPD mini project\n\n";
    std::cout << "Using CL=" << CL << " L/h, V=" << V << " L\n";
    
    double ke = getKConst(CL, V); // h^-1
    double t12 = get_half_life(ke); // h
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Ke=" << ke << " 1/h, Half-Life=" << t12 << " h\n" << std::endl; 
    
    double dose = ask_dose();

    double t_end = 5 * t12; 
    double dt = t_end/10; // This could later be user-defined

    std::cout << "5 half-lifes of the drug, displayed at intervals of " << dt << ":"<< std::endl;

    std::vector<double> times = build_time_vectors(t_end, dt);

    for (std::size_t i = 0; i< times.size(); ++i) {
        double c = conc_single_iv_bolus(dose, V, ke, times[i]);
        std::cout << "t=" << times[i] << " h, C=" << c << " mg/L" << std::endl;
    }

    return 0;
}


// double c = conc_single_iv_bolus(dose, V, ke, t);
// std::cout << "For " << dose << " mg of dose, the concentration will be " << c << " mg after " << t << " hour." << std::endl;
