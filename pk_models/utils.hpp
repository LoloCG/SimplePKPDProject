#pragma once

#include <vector>
#include <optional>
#include <string>
#include <filesystem>

struct DisplayPoint {double time, Ag, Ac;};

// output from RegimenBuilder, input for the PK classes.
struct TimelinePoint {double time; std::optional<double> dose;};
// Input for RegimenBuilder
struct DoseEvent {double time, dose;};

namespace RegimenBuilder{
    std::vector<TimelinePoint> generate_regimen_timeline(
        const std::vector<DoseEvent>& regimen, 
        const std::vector<double>& time_steps, 
        double eps = 1e-12
    );

    std::vector<DoseEvent> regular_dose_regimen(double dose, double tau, std::size_t n_doses);
    
    double end_time_by_hl_of_doses(double t12, std::size_t n_doses, 
        double decay_mod = 5, double dose_mod = 2);

    std::vector<double> time_steps_by_delta(double t_end, double dt);

    std::vector<double> time_steps_by_n(double t_end, std::size_t steps_n = 30);
};

// Generated with GPT-5
namespace exportutil {
    inline std::string num_to_string(double x, int precision, bool decimal_comma);

    inline bool save_for_excel(
        const std::filesystem::path& out_path,
        const std::vector<DisplayPoint>& rows,
        char delimiter = ';',
        int precision = 6,
        bool decimal_comma = true,
        bool write_sep_hint = true,
        bool include_header = true
    );
} // namespace exportutil

double ask_dose();