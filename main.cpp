// #include "pk_models/one_compartment.hpp"
#include "pk_models/utils.hpp"

#include <iostream>
#include <cmath>    // For std::log
#include <iomanip>  // for formatting string stream
#include <vector>
// #include <stdexcept>
// #include <limits>   // For std::numeric_limits
// #include <optional> 
// #include <fstream>
// #include <sstream>
// #include <iomanip>
// #include <string>
// #include <locale>
// #include <filesystem>

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

// ------------------------------ PK classes ------------------------------
// struct DisplayPoint {double time, Ag, Ac;};
class CCompartment {
private:
    double CL_; // L/h
    double V_;  // L
    double ke_;
    double ka_;
    double F_; // bioavailability

    double e_exp(double delta_t) {
        return std::exp(-ke_ * delta_t);
    }

    double a_exp(double delta_t) {
        return std::exp(-ka_ * delta_t);
    }

public:
    CCompartment(double ke, double ka, double F) : ke_(ke), ka_(ka), F_(F) {
        if (ka <= ke_) throw std::invalid_argument("Flip-flop kinetics not added. Must have ka > ke");
        
        Ac = 0;
        Ag = 0;
    }

    double Ac; // state of current ammount in compartment.
    double Ag; // state of current ammount in extravascular compartment.

    void add_depot_dose(double dose) {
        Ag += dose * F_; 
    }

    double absorption_change(double delta_t) {
        return (ka_/(ka_-ke_)) * Ag *(e_exp(delta_t)-a_exp(delta_t));
    }

    // Returns the remaining amount of the central after applying first order decay.
    double first_order_decay_central(double eexp) {
        return Ac * eexp;
    }

    double first_order_decay_depot(double delta_t, double aexp){
        return Ag * aexp;
    }

    std::vector<DisplayPoint> propagate_distribution(const std::vector<TimelinePoint>& time_data) {
        std::vector<DisplayPoint> data;
        data.reserve(time_data.size());
        
        double last_t = time_data.front().time;
        
        {
            const auto& pt0 = time_data.front();
            if (pt0.dose) add_depot_dose(pt0.dose.value());
            data.emplace_back(pt0.time, Ag, Ac);
        }
        
        for (size_t i = 1; i < time_data.size(); ++i) {
            const auto& pt = time_data[i];
            const double t = pt.time;
            const double delta_t = t - last_t;

            if (delta_t == 0) {
                if (pt.dose) add_depot_dose(pt.dose.value());
                data.emplace_back(t, Ag, Ac);
                continue;
            }
            
            const double aexp = a_exp(delta_t);
            const double eexp = e_exp(delta_t);

            const double Ag_end = first_order_decay_depot(delta_t, aexp);
            
            double Acr_decay = 0;
            if (Ac > 0) Acr_decay = first_order_decay_central(eexp);
            const double Acr_abs = absorption_change(delta_t);

            Ag = Ag_end;
            Ac = Acr_decay + Acr_abs;
            last_t = t;

            if (pt.dose) add_depot_dose(pt.dose.value());

            data.emplace_back(t, Ag, Ac);
        }

        return data;
    }
};

// ------------------------------ abstractions ------------------------------

// ------------------------------ main ------------------------------

int main() {
    const double t12 = 26;
    const std::size_t n_doses = 5;
    const double tau = 24;
    const double dose = 50;
    const double steps = 1; // h
    const int F = 1;

    std::cout << "PKPD mini project\n\n";

    const double ke = std::log(2)/t12;
    const double ka = 1;

    std::cout << std::setprecision(3) << "Using half-life=" << t12 << " h, ka="
        << ka << " 1/h, ke=" << ke << " 1/h\n" << std::setprecision(1)
        << "Regimen of " << n_doses << " doses of " << std::setprecision(3) 
        << dose << " mg every " << std::setprecision(1) << tau << " h." << std::endl;

    double end_t = RegimenBuilder::end_time_by_hl_of_doses(t12, n_doses);
    std::vector<double> time_steps = RegimenBuilder::time_steps_by_delta(end_t, steps);
    std::vector<DoseEvent> dosage_regimen = RegimenBuilder::regular_dose_regimen(dose, tau, n_doses);
    std::vector<TimelinePoint> timeLine = RegimenBuilder::generate_regimen_timeline(dosage_regimen, time_steps);

    CCompartment compartment(ke, ka, F);
    const std::vector<DisplayPoint> data = compartment.propagate_distribution(timeLine);
    
    std::filesystem::path home =
    #ifdef _WIN32
        std::getenv("USERPROFILE");
    #else
        std::getenv("HOME");
    #endif
    auto out = home / "Desktop" / "pk_output.csv";   // adjust if your Desktop is different
    exportutil::save_for_excel(out, data);
    
    /*
    for (std::size_t i = 0; i < data.size(); ++i) {
        std::cout << data[i].time << ";\t" << std::setprecision(3) 
        << data[i].Ag << ";\t" << (data[i].Ac) << "\n";
    }
    */

    return 0;
}

/*
- Oral, immediate-release (small molecules): ka = 0.5–2 
- Oral, modified/extended-release: ka = 0.05–0.3
    - zero order or transit models may fit better 
- Sublingual/buccal/inhaled (fast): ka = 2–10
- IM/SC (solutions of small molecules): ka = 0.2–1.5
    - depends on site and formulation
- SC biologics / long-acting depots / transdermal: slow: ka = 0.005–0.1 h−1
*/