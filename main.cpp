#include <iostream>
#include <cmath>    // For std::log
#include <iomanip>  // for formatting string stream
#include <vector>
#include <stdexcept>
#include <limits>   // For std::numeric_limits
#include <optional> // For std::optional

struct pk_calcs {
    // elimination constant. first order kinetics (?)
    static double ke(double CL, double V) {return CL/V;}

    static double half_life(double ke) {
        return std::log(2.0) / ke; // t1/2 = ln(2) / ke
    }

    /// weight in kg, age in years, creatinine plasma conc. (creat_cp) in mg/dl. 
    /// Result in ml/min. 
    static double creatinine_cl(double age, double weight, double creat_cp, bool female = false) {
        double creat_cp_mod = 72.0; // modifier constant for creatinine concentration in plasma
        if (female) creat_cp_mod = 85.0;
        
        return ((140 - age) * weight)/ creat_cp_mod * creat_cp;
    }
};

struct Dose { 
    double time_h; 
    double amount_mg; 
};

// 1-compartment, first-order elimination
class one_compartment_pk {
private:
    double CL_; // L/h
    double V_;  // L

public:
    one_compartment_pk(double CL, double V)
        : CL_(CL), V_(V)
    {
        if (CL_ <= 0.0) throw std::invalid_argument("CL must be > 0.");
        if (V_  <= 0.0) throw std::invalid_argument("V must be > 0.");
    }
   
    double conc_iv_bolus_at(double dose, double time, double ke) const {
        return (dose / V_) * std::exp(- ke * (time)); // C(t) = (dose/V) * e^(-ke *t)
    }
    
    double conc_iv_bolus_at_with_prev(double dose, double last_dose_time, double ke, double tau) {
        return conc_iv_bolus_at(dose, last_dose_time, ke)/(1 - std::exp(- ke * tau)); // C(t) = dose/V * e^(-ke * t) / (1-e^(-ke*tau))
    }

    std::vector<double> build_c_at_times_single_bolus(double dose, std::vector<double> times) {
        std::vector<double> c_at_times;
        const double ke = pk_calcs::ke(CL_, V_);
        // const double tau = 12;
        
        for (std::size_t i = 0; i< times.size(); ++i) {
            double c = conc_iv_bolus_at(dose, times[i], ke);
            // double c = conc_iv_bolus_at_with_prev(dose, times[i], ke, tau);
            c_at_times.push_back(c);
        }

        return c_at_times;
    }

    std::vector<double> build_c_at_times_multiple_bolus(double dose, std::vector<double> times) {
        std::vector<double> c_at_times;
        const double ke = pk_calcs::ke(CL_, V_);
        // const double tau = 12;
        
        // double prev_c = dose / V_; // start with initial c at time 0. 
        
        for (std::size_t i = 0; i< times.size(); ++i) {
            // double t_since_last = times[i] -  
            double c = conc_iv_bolus_at(dose, times[i], ke);

            c_at_times.push_back(c);
            // prev_c = c;
        }

        return c_at_times;
    }

    double avg_c_ss(double dose, int tau) {
        // (dose / tau) / CL
        return (dose / tau) / CL_;
    }

    double min_c_ss() {
        // Cmaxss ​e−ke​τ
    }
};

// ------------------------------ Utilities ------------------------------
std::vector<double> build_time_steps(double t_end, std::optional<double> dt) {
    double step = dt.value_or(t_end / 200.0);

    if (t_end < 0.0)    throw std::invalid_argument("t_end must be >= 0.");
    if (step <= 0.0)    throw std::invalid_argument("dt must be > 0.");
    
    const std::size_t n = static_cast<std::size_t>(std::floor(t_end / step) + 1.0); // how many multiples of dt between 0 and t_end 

    std::vector<double> times; // Creates an empty dynamic array.
    times.reserve(n); // preallocate size.

    for (std::size_t i = 0; i < n; ++i) {
        times.push_back(static_cast<double>(i) * step);
    }
    return times;
}

double ask_dose() {
    // std::cin >> std::ws; // This causes issue of asking input before showing prompt.
    
    double dose;
    std::cout << "Enter dose in mg: " << std::flush;
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


int main() {
    const double CL = 3.0;
    const double V = 30.0;

    const int n_doses = 2;
    const int tau = 12; // (h) time between doses

    std::cout << "PKPD mini project\n\n";
    std::cout << "Using CL=" << CL << " L/h, V=" << V << " L\n";
    

    double ke = pk_calcs::ke(CL, V); // h^-1
    double t12 = pk_calcs::half_life(ke); // h
    double t_end = 5 * t12; 
    double dt = t_end/10; // 10 can later be changed to std::size_t plot_points optional parameter.

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Ke=" << ke << " 1/h, Half-Life=" << t12 << " h\n" << std::endl; 
    
    double dose = ask_dose();
    
    one_compartment_pk model(CL, V); 

    const double avgcss = model.avg_c_ss(dose, tau);

    std::cout 
        << "The concentration at steady state is " << avgcss << " with infusions at "
        // << ", reached after " << n_doses << " doses at " // this could be calculated later...
        << tau << " h intervals."  << std::endl;

    auto times = build_time_steps(t_end, dt);
    auto concs = model.build_c_at_times_single_bolus(dose, times);

    std::cout << "5 half-lifes of the drug, displayed at intervals of " << dt << ":"<< std::endl;
    std::cout << "\nTime (h), Conc (mg/L)\n";

    for (std::size_t i = 0; i < times.size(); ++i) {
        std::cout << std::setw(7) << std::setprecision(2) << std::fixed 
        << times[i] << "   " << std::setprecision(4) << concs[i] << "\n";
    }
    
    return 0;
}


// double c = conc_single_iv_bolus(dose, V, ke, t);
// std::cout << "For " << dose << " mg of dose, the concentration will be " << c << " mg after " << t << " hour." << std::endl;
