#include "pk_models/one_compartment.hpp"

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

// Used to define the timing and quantity of each dosage for injection, oral, or other treatment.
struct RegimenDose {double time_h, dose;};

// Used for the final data.
struct DisplayData {double time, conc;};

class IrregularIVBolusRegimen{
private:
    double CL_; // L/h
    double V_;  // L

public:
    IrregularIVBolusRegimen(double CL, double V) : CL_(CL), V_(V) {
        if (CL_ <= 0.0) throw std::invalid_argument("CL must be > 0.");
        if (V_  <= 0.0) throw std::invalid_argument("V must be > 0.");
    }

    static std::vector<RegimenDose> generate_regular_regimen(double dose, double tau, std::size_t n_doses) {
        std::vector<RegimenDose> regimen;
        regimen.reserve(n_doses);
        
        for (std::size_t i = 0; i<n_doses; ++i) {
            regimen.emplace_back(tau*i, dose);
        }

        return regimen;
    }

    double decay_at_t(double residual_c, double ke, double t, double t_last) { 
        return residual_c * std::exp(-ke *(t-t_last));  // C(t)=Cr*e^(-ke*(t-tlast))
    }
    double bolus_at_t(double residual_c, double dose) {
        return residual_c + dose / V_;  // C(t)=Cr + Dk / V
    }

    // Assumes that both vector params are already time sorted. 
    std::vector<double> merge_times(const std::vector<RegimenDose>& regimen, const std::vector<double>& time_steps, double eps = 1e-12) {
        std::vector<double> merged_times;
        merged_times.reserve(regimen.size() + time_steps.size());
        
        // --- generated with gpt-5 --- //
        // lambda expression with captured by reference (modify). takes parameter t of timestamp
        auto push_if_new = [&merged_times, &eps](double t) { 
            // if vector is empty or distance from the last append is greater than eps, push back. Otherwise skip (treat as duplicate)
            if (merged_times.empty() || std::abs(merged_times.back() - t) > eps)
                merged_times.push_back(t);
        };

        std::size_t si = 0, di = 0;  // indexes for steps and dose
        while (si < time_steps.size() || di < regimen.size()) {
            
            // Take from steps if: doses are exhausted, OR (both present and step <= dose)
            bool takeStep =
                (di == regimen.size()) ||
                (si < time_steps.size() && time_steps[si] <= regimen[di].time_h);

            // push step time, advance step index
            if (takeStep)   push_if_new(time_steps[si++]);
            // push dose time, advance dose index
            else            push_if_new(regimen[di++].time_h);
        }
        
        return merged_times;
    }

    // TODO: this has too many sidecases which are not assumed to occur:
    // - doses need to be sufficiently spaced between each in the regimen.
    // - doses cant be stacked at same time without merging beforehand.
    std::vector<DisplayData> generate_regimen_data(std::vector<RegimenDose> regimen, std::vector<double> time_steps) {
        const double ke = pk_calcs::ke(CL_, V_);
        std::vector<double> times = merge_times(regimen, time_steps); 
        
        std::vector<DisplayData> data;
        data.reserve(times.size());
        
        double C = 0;
        std::size_t r = 0;
        for (std::size_t i = 0; i < times.size(); ++i) {
            double t = times[i];
            double last_t = (i == 0) ? t : times[i-1];
            
            C = decay_at_t(C, ke, t, last_t);
            last_t = t;
        
            if (r < regimen.size() && regimen[r].time_h == t) {
                C = bolus_at_t(C, regimen[r].dose);
                ++r;
            }
            
            data.emplace_back(t, C);
        }

        return data;
    }
};

// ------------------------------ Utilities ------------------------------
class TimeStepBuilder {
private:
    double t_end;

public:
    double get_end_time_required(double t12, std::size_t n_doses = 1) {
        // The idea is to display 5*t12 of decay, and + 3*t12 for each dose > 1.
        double t_end = (3.0 * (n_doses - 1) + 5.0) * t12;
        return t_end;
    }

    std::vector<double> build_time_steps(double t_end, std::size_t steps_n = 30, std::optional<double> dt = std::nullopt) {
        double step = dt.value_or(t_end / steps_n);

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
};

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

void display_single_iv_bolus(double CL, double V, double dose, std::vector<double> times) {
    SingleIVBolus_1C model(CL, V); 
    auto concs = model.build_c_at_times(dose, times);

    std::cout << "\nTime (h), Conc (mg/L)\n";

    for (std::size_t i = 0; i < times.size(); ++i) {
        std::cout << std::setw(7) << std::setprecision(2) << std::fixed 
        << times[i] << "   " << std::setprecision(4) << concs[i] << "\n";
    }
}

// ------------------------------ main ------------------------------

int main() {
    const double CL = 3.0;
    const double V = 30.0;

    std::cout << "PKPD mini project\n\n";
    std::cout << "Using CL=" << CL << " L/h, V=" << V << " L\n";

    double ke = pk_calcs::ke(CL, V); // h^-1
    double t12 = pk_calcs::half_life(ke); // h
    std::cout << std::setprecision(4) << "Half-life=" << t12 << std::endl;
    
    const std::size_t n_doses = 3;
    const double tau = 12.0;

    TimeStepBuilder tbuild;
    double t_end = tbuild.get_end_time_required(t12, n_doses);
    std::vector<double> time_steps = tbuild.build_time_steps(t_end);
 
    double dose = ask_dose();
    IrregularIVBolusRegimen IIVBolus(CL, V);
    std::vector<RegimenDose> regimen = IIVBolus.generate_regular_regimen(dose, tau, n_doses);
    std::vector<DisplayData> data = IIVBolus.generate_regimen_data(regimen, time_steps);

    // /*
    std::cout << "time (h);C (mg/L)\n";
    for (std::size_t i = 0; i < data.size(); ++i) {
        // std::cout << "T=" << data[i].time 
        //     << std::setprecision(3) << "\tc=" << data[i].conc 
        //     << "\n";

        std::cout << data[i].time << std::setprecision(3) << ";" << data[i].conc << "\n";
    }
    // */

    return 0;
}

/* Example 1 compartment single IV bolus
    double ke = pk_calcs::ke(CL, V); // h^-1
    double t12 = pk_calcs::half_life(ke); // h
    double t_end = 5 * t12; 
    double dt = t_end/10; // hardcoded 10 data points to display

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Ke=" << ke << " 1/h, Half-Life=" << t12 << " h\n" << std::endl; 

    auto times = build_time_steps(t_end, dt);
    double dose = ask_dose();
    std::cout << "5 half-lifes of the drug, displayed at intervals of " << dt << ":"<< std::endl;
    display_single_iv_bolus(CL, V, dose, times);

*/