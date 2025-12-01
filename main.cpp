#include <iostream> // For std::cout and std::cerr
#include <cmath>    // For std::log, std::exp, std::abs, ...
#include <iomanip>  // for number formatting
#include <vector>
#include <string>
#include <optional>
#include <filesystem>
#include <fstream> // std::ofstream
#include <sstream> // std::ostringstream (exportutil)
#include <cstring>
#include <cerrno> // std::strerror

// ---------------------------------------- Helpers ----------------------------------------

double parse_double(const char* s) {
    try {
        return std::stod(s);
    } catch (...) {
        throw std::runtime_error(std::string("Invalid double: ") + s);
    }
}

int parse_int(const char* s) {
    try {
        return std::stoi(s);
    } catch (...) {
        throw std::runtime_error(std::string("Invalid integer: ") + s);
    }
}

std::filesystem::path get_desktop_path(const std::string filename="pk_output.csv"){
    std::filesystem::path home =
    #ifdef _WIN32
        std::getenv("USERPROFILE");
    #else
        std::getenv("HOME");
    #endif

    return home / "Desktop" / filename;
}

void print_usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " --dose VALUE [--t12 VALUE | --kel VALUE]\n"
        << "                [--ka VALUE] [--ndoses N] [--tau HOURS]\n"
        << "                [--f VALUE] [--step-size HOURS]\n"
        << "                [--out PATH]\n\n"

        << "Required:\n"
        << "  --dose VALUE        Dose amount (double > 0).\n"
        << "  --t12 VALUE         Half-life in hours (double > 0).\n"
        << "     OR\n"
        << "  --kel VALUE         Elimination rate constant in 1/h (double > 0).\n\n"

        << "Optional:\n"
        << "  --ka VALUE          Absorption rate constant ka (default 0.1).\n"
        << "  --ndoses N          Number of doses (default 1).\n"
        << "  --tau HOURS         Interval between doses; needed if ndoses > 1.\n"
        << "  --f VALUE           Bioavailability F (default 1).\n"
        << "  --step-size HOURS   Time step for simulation grid (default 1).\n"
        << "  --out PATH          Output CSV path. If omitted, writes to Desktop.\n\n"

        << "General:\n"
        << "  --help, -h          Show this message.\n";
}

// ---------------------------------------- structs ----------------------------------------

struct DisplayPoint {double time, Ag, Ac;};
struct DoseEvent {double time, dose;};
struct TimelinePoint {double time; std::optional<double> dose;};

struct SimParams {
    // required params
    double dose;    
    double t12;
    bool has_t12 = false;
    double kel;
    bool has_kel = false;

    // Optional with defaults
    double ka = 0.1;
    double f = 1;
    std::size_t n_doses = 1;
    double step_size = 1; // time
    double tau = 24.0;
    
    bool evroute = false;

    std::string out_path = "pk_output.csv";
    // std::string out_path;
};
SimParams parse_args(int argc, char** argv) {
    SimParams p;

    if (argc == 1) {
        throw std::runtime_error("No arguments provided");
    }

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        auto need_value = [&](const std::string& name) {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value for " + name);
            }
        };

        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else if (arg == "--dose") {
            need_value(arg);
            p.dose = parse_double(argv[++i]);
        } else if (arg == "--t12") {
            need_value(arg);
            p.t12 = parse_double(argv[++i]);
            p.has_t12 = true;
        } else if (arg == "--kel") {
            need_value(arg);
            p.kel = parse_double(argv[++i]);
            p.has_kel = true;
        } else if (arg == "--ka") {
            need_value(arg);
            p.ka = parse_double(argv[++i]);
        } else if (arg == "--ndoses") {
            need_value(arg);
            p.n_doses = parse_int(argv[++i]);
        } else if (arg == "--tau") {
            need_value(arg);
            p.tau = parse_double(argv[++i]);
        } else if (arg == "--f") {
            need_value(arg);
            p.f = parse_double(argv[++i]);
        } else if (arg == "--out" || arg == "-o") {
            need_value(arg);
            p.out_path = argv[++i];
        } else {
            throw std::runtime_error("Unknown argument: " + arg + "\n");
            print_usage(argv[0]);
        }
    }
    
    if (p.dose <= 0.0 || !p.dose) {
        throw std::runtime_error("dose must be > 0");
    }
    
    // ---------- t12 and kel checks ----------
    if (!p.has_t12 && !p.has_kel) {
        throw std::runtime_error("You must pass either --t12 or --kel");
    }
    
    const double ln2 = std::log(2.0);

    if (p.has_t12 && !p.has_kel) {
        if (p.t12 <= 0.0) {
            throw std::runtime_error("t12 must be > 0");
        }
        p.kel = ln2 / p.t12;
        p.has_kel = true;
    } else if (p.has_kel && !p.has_t12) {
        if (p.kel <= 0.0) {
            throw std::runtime_error("kel must be > 0");
        }
        p.t12 = ln2 / p.kel;
        p.has_t12 = true;
    } else if (p.has_t12 && p.has_kel) {
        double kel_from_t12 = ln2 / p.t12;
        double rel_diff = std::fabs(kel_from_t12 - p.kel) / kel_from_t12;
        if (rel_diff > 0.1) {
            throw std::runtime_error("Inconsistent t12 and kel");
        }
    }

    return p;
}

// ---------------------------------------- PK classes ----------------------------------------

class CCompartment {
    // Outputs post-dose and post-propagate drug amounts
private:
    // double CL_; // L/h
    // double V_;  // L
    double kel_;
    double ka_;
    double F_; // bioavailability
    bool evroute_;

    double e_exp(double delta_t) {
        return std::exp(-kel_ * delta_t);
    }

    double a_exp(double delta_t) {
        return std::exp(-ka_ * delta_t);
    }

public:
    CCompartment(const SimParams& p) : 
        kel_(p.kel), ka_(p.ka), F_(p.f),
        evroute_(p.evroute)
    {        
        Ac = 0;
        Ag = 0;
    }

    double Ac; // state of current ammount in compartment.
    double Ag; // state of current ammount in extravascular compartment.

    void add_depot_dose(double dose) {
        Ag += dose * F_; 
    }

    double first_order_abs_change(double delta_t) {
        return Ag * (ka_/(ka_-kel_)) * (e_exp(delta_t)-a_exp(delta_t));
    }

    double collapsed_abs_change(double delta_t) {
        return Ag * kel_ * e_exp(delta_t) * delta_t;
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

        // Dictates the threshold for (ka-ke)*dt to avoid equation collapse. 
        const double exp_limit = 1e-6; 
        
        for (size_t i = 1; i < time_data.size(); ++i) {
            const auto& pt = time_data[i];
            const double t = pt.time;
            const double delta_t = t - last_t;

            const double aexp = a_exp(delta_t);
            const double eexp = e_exp(delta_t);

            const double Ag_end = first_order_decay_depot(delta_t, aexp);

            const double Acr_decay = first_order_decay_central(eexp);            
            
            const double x = std::abs((ka_-kel_)*delta_t);
            const double Acr_abs = (x < exp_limit)
                ? collapsed_abs_change(delta_t)
                : first_order_abs_change(delta_t);
            
            Ag = Ag_end;
            Ac = Acr_decay + Acr_abs;
            last_t = t;
            
            if (pt.dose) add_depot_dose(pt.dose.value());

            data.emplace_back(t, Ag, Ac);
        }

        return data;
    }
};

namespace exportutil {
    std::string num_to_string(double x, int precision, bool decimal_comma) {
        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss << std::setprecision(precision) << x;
        std::string s = oss.str();
        if (decimal_comma) {
            for (char& c : s) if (c == '.') c = ',';
        }
        return s;
    }

    bool save_for_excel(
        const std::filesystem::path& out_path,
        const std::vector<DisplayPoint>& rows,
        char delimiter = ';',
        int precision = 6,
        bool decimal_comma = true,
        bool write_sep_hint = true,
        bool include_header = true
    ) {
        std::error_code ec;
        auto parent = out_path.parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent, ec); // ignore error; open will fail if truly unwritable

        std::ofstream ofs(out_path, std::ios::out | std::ios::trunc);
        if (!ofs) {
            std::cerr << "Failed to open " << out_path << ": " << std::strerror(errno) << "\n";
            return false;
        }

        if (write_sep_hint) ofs << "sep=" << delimiter << "\r";
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
};

// ---------------------------------------- Timepoint builder ----------------------------------------

namespace RegimenBuilder{
    /// Generates the dosage regimen based on a static schedule of time between events (tau) and the number of doses.
    std::vector<DoseEvent> regular_dose_regimen(double dose, double tau, std::size_t n_doses){
        std::vector<DoseEvent> regimen;
        regimen.reserve(n_doses);
        for (std::size_t i = 0; i<n_doses; ++i) {
            regimen.emplace_back(tau*i, dose);
        }
        return regimen;
    }


    /// Calculates the maximum time of the regimen timeline by using the 7x half-lives of last dose.
    /// Asumes regimen of same doses at regular intervals.
    double end_time_regular_dose_by_hl(
        double t12, 
        double tau,
        std::size_t n_doses, 
        double decay_mod = 7
    ) {
        double end_t;
        if (n_doses > 1) {
            end_t = ((n_doses-1) * tau) + decay_mod * t12;
        } else {
            end_t = decay_mod * t12;
        }
        std::cout << "calculated end_t=" << end_t << std::endl;
        return end_t;
    }


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


    std::vector<double> time_steps_by_delta(double t_end, double dt) {
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

// ---------------------------------------- main ----------------------------------------

std::vector<TimelinePoint> build_timepoints(const SimParams& p) {
    double end_t = RegimenBuilder::end_time_regular_dose_by_hl(p.t12, p.tau, p.n_doses);
    std::vector<double> time_steps = RegimenBuilder::time_steps_by_delta(end_t, p.step_size);
    std::vector<DoseEvent> dosage_regimen = RegimenBuilder::regular_dose_regimen(p.dose, p.tau, p.n_doses);
    std::vector<TimelinePoint> timeLine = RegimenBuilder::generate_regimen_timeline(dosage_regimen, time_steps);

    return timeLine;
}

int main(int argc, char *argv[]){
    try {
        SimParams p = parse_args(argc, argv);
        
        std::vector<TimelinePoint> timeLine = build_timepoints(p);

        CCompartment compartment(p);
        const std::vector<DisplayPoint> data = compartment.propagate_distribution(timeLine);
        
        std::filesystem::path out;
        if (p.out_path.empty()) {
            out = get_desktop_path("pk_output.csv");
        } else {
            out = p.out_path;
        }
        exportutil::save_for_excel(out, data);

        return 0;

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
}