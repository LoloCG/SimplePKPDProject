#include "utils.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <locale>
#include <fstream>
#include <sstream>

/* ------------------------------ RegimenBuilder ------------------------------ */

// Assumes that both vector params are already time sorted.
std::vector<TimelinePoint> RegimenBuilder::generate_regimen_timeline(
    const std::vector<DoseEvent>& regimen, 
    const std::vector<double>& time_steps, 
    double eps
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
std::vector<DoseEvent> RegimenBuilder::regular_dose_regimen(double dose, double tau, std::size_t n_doses) {
    std::vector<DoseEvent> regimen;
    regimen.reserve(n_doses);
    
    for (std::size_t i = 0; i<n_doses; ++i) {
        regimen.emplace_back(tau*i, dose);
    }

    return regimen;
}

// Calculates the end time of the time to display based on the half life and number of dosages.
// Each dose is 2*t12, while the last 5*t12 to display the decay. 
double RegimenBuilder::end_time_by_hl_of_doses(double t12, std::size_t n_doses, double decay_mod, double dose_mod) {
    return (dose_mod * (n_doses - 1) + decay_mod) * t12;
}

std::vector<double> RegimenBuilder::time_steps_by_delta(double t_end, double dt) {
    // TODO: ask why const. 
    const std::size_t n = static_cast<std::size_t>(std::floor(t_end / dt) + 1.0);

    std::vector<double> times;
    times.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        times.push_back(static_cast<double>(i) * dt);
    }
    return times;
}

std::vector<double> RegimenBuilder::time_steps_by_n(double t_end, std::size_t steps_n) {
    double dt = t_end / steps_n;
    return time_steps_by_delta(t_end, dt);
}


/* ------------------------------ exportutil ------------------------------ */

// Generated with GPT-5
inline std::string exportutil::num_to_string(double x, int precision, bool decimal_comma) {
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(precision) << x;
    std::string s = oss.str();
    if (decimal_comma) {
        for (char& c : s) if (c == '.') c = ',';
    }
    return s;
}

inline bool exportutil::save_for_excel(
    const std::filesystem::path& out_path,
    const std::vector<DisplayPoint>& rows,
    char delimiter,
    int precision,
    bool decimal_comma,
    bool write_sep_hint,
    bool include_header) {
    // Create parent dir if provided
    std::error_code ec;
    auto parent = out_path.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent, ec); // ignore error; open will fail if truly unwritable

    std::ofstream ofs(out_path, std::ios::out | std::ios::trunc);
    if (!ofs) {
        std::cerr << "Failed to open " << out_path << ": " << std::strerror(errno) << "\n";
        return false;
    }

    if (write_sep_hint) ofs << "sep=" << delimiter << "\r\n";
    if (include_header) ofs << "time" << delimiter << "Ag" << delimiter << "Ac" << "\r";

    for (const auto& r : rows) {
        ofs  << num_to_string(r.time, precision, decimal_comma) << delimiter
                << num_to_string(r.Ag,   precision, decimal_comma) << delimiter
                << num_to_string(r.Ac,   precision, decimal_comma) << "\r";
    }
    ofs.flush();
    if (!ofs) {
        std::cerr << "Write failed for " << out_path << "\n";
        return false;
    }
    std::cout << "Wrote " << out_path << " (" << rows.size() << " rows)\n";
    return true;
}

/* ------------------------------ Helpers ------------------------------ */

double ask_dose() {
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
