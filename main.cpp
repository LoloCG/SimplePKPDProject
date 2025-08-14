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



// ------------------------------ Utilities ------------------------------
struct DoseEvent {double time, dose;}; // Input
struct TimelinePoint {double time; std::optional<double> dose;}; // output from RegimenBuilder, input for the PK classes.
namespace RegimenBuilder {

    // Assumes that both vector params are already time sorted.
    std::vector<TimelinePoint> generate_regimen_timeline(
        const std::vector<DoseEvent>& regimen, 
        const std::vector<double>& time_steps, 
        double eps = 1e-12
    ) {
        
        std::vector<TimelinePoint> merged;
        merged.reserve(regimen.size() + time_steps.size());


        // Append a (time, dose?) pair; if same time as last (within eps), merge doses instead.
        auto push_or_merge = [&](double t, std::optional<double> dose) {
            if (merged.empty() || std::abs(merged.back().time - t) > eps
                ) {
                    merged.push_back(TimelinePoint{t, dose});
            } else {
                if (dose) {
                    if (merged.back().dose) merged.back().dose = *merged.back().dose + *dose;
                    else merged.back().dose = dose; // set dose on an existing time step
                }
            }
        };

        std::size_t t = 0, d = 0;
        while (t < time_steps.size() || d < regimen.size()) {
            
            // take from time_steps if no more dosages or when steps<=dose 
            bool takeTime = (
                d == regimen.size() || 
                (t < time_steps.size() && 
                time_steps[t] <= regimen[d].time)
            );
        
            if (takeTime) {
                push_or_merge(time_steps[t], std::nullopt);
                ++t;
            } else {
                push_or_merge(regimen[d].time, std::optional<double>(regimen[d].dose));
                ++d;
            }
        }

        return merged;
    }

    // Generates the dosage regimen based on a static schedule of time between events (tau) and the number of doses.
    std::vector<DoseEvent> regular_dose_regimen(double dose, double tau, std::size_t n_doses) {
        std::vector<DoseEvent> regimen;
        regimen.reserve(n_doses);
        
        for (std::size_t i = 0; i<n_doses; ++i) {
            regimen.emplace_back(tau*i, dose);
        }

        return regimen;
    }

    // Calculates the end time of the time to display based on the half life and number of dosages.
    // Each dose is 2*t12, while the last 5*t12 to display the decay. 
    double end_time_by_hl_of_doses(double t12, std::size_t n_doses, double decay_mod = 5, double dose_mod = 2) {
        return (dose_mod * (n_doses - 1) + decay_mod) * t12;
    }

    std::vector<double> time_steps_by_delta(double t_end, double dt) {
        // TODO: ask why const. 
        const std::size_t n = static_cast<std::size_t>(std::floor(t_end / dt) + 1.0);

        std::vector<double> times;
        times.reserve(n);

        for (std::size_t i = 0; i < n; ++i) {
            times.push_back(static_cast<double>(i) * dt);
        }
        return times;
    }

    std::vector<double> time_steps_by_n(double t_end, std::size_t steps_n = 30) {
        double dt = t_end / steps_n;
        return time_steps_by_delta(t_end, dt);
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

struct DisplayPoint {double time, conc;};
class IrregularIVBolusRegimen{
private:
    double CL_; // L/h
    double V_;  // L

public:
    IrregularIVBolusRegimen(double CL, double V) : CL_(CL), V_(V) {
        if (CL_ <= 0.0) throw std::invalid_argument("CL must be > 0.");
        if (V_  <= 0.0) throw std::invalid_argument("V must be > 0.");
    }

    double decay_at_t(double residual_c, double ke, double t, double t_last) { 
        return residual_c * std::exp(-ke *(t-t_last));  // C(t)=Cr*e^(-ke*(t-tlast))
    }
    double bolus_at_t(double residual_c, double dose) {
        return residual_c + dose / V_;  // C(t)=Cr + Dk / V
    }


    // TODO: this has too many sidecases which are not assumed to occur:
    // - doses need to be sufficiently spaced between each in the regimen.
    // - doses cant be stacked at same time without merging beforehand.
    std::vector<DisplayPoint> multiple_iv_bolus_data(const std::vector<TimelinePoint>& time_data) {
        const double ke = pk_calcs::ke(CL_, V_);
        
        std::vector<DisplayPoint> data;
        data.reserve(time_data.size());
        
        double C = 0;
        double last_t = time_data.front().time;
        for (const auto& pt : time_data) {
            const double t = pt.time;
            
            C = decay_at_t(C, ke, t, last_t);
            last_t = t;
            
            if (pt.dose) C = bolus_at_t(C, pt.dose.value()); 

            data.emplace_back(t, C);
        }

        return data;
    }
};

// ------------------------------ main ------------------------------

int main() {
    const double CL = 3.0;
    const double V = 30.0;
    const std::size_t n_doses = 3;
    const double tau = 12.0;
    const double dose = 150; 

    std::cout << "PKPD mini project\n\n";
    std::cout << "Using CL=" << CL << " L/h, V=" << V << " L\n";

    double ke = pk_calcs::ke(CL, V); // h^-1
    double t12 = pk_calcs::half_life(ke); // h
    std::cout << n_doses << " doses of" << dose << " mg every " << tau << " h" << std::setprecision(4) << ", half-life=" << t12 <<std::endl;
    // double dose = ask_dose();


    double end_t = RegimenBuilder::end_time_by_hl_of_doses(t12, n_doses);
    std::vector<double> time_steps = RegimenBuilder::time_steps_by_delta(end_t, 1);
    std::vector<DoseEvent> dosage_regimen = RegimenBuilder::regular_dose_regimen(dose, tau, n_doses);
    std::vector<TimelinePoint> timeLine = RegimenBuilder::generate_regimen_timeline(dosage_regimen, time_steps);

    IrregularIVBolusRegimen iv_bolus(CL, V);
    std::vector<DisplayPoint> data = iv_bolus.multiple_iv_bolus_data(timeLine);

    // /*
    std::cout << "Time (h);C (mg/L)\n";
    for (std::size_t i = 0; i < data.size(); ++i) {
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