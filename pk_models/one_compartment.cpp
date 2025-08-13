#include <stdexcept>
#include <cmath>
#include <vector>
#include "one_compartment.hpp"

SingleIVBolus_1C::SingleIVBolus_1C(double CL, double V) : CL_(CL), V_(V) {
  if (CL_ <= 0.0) throw std::invalid_argument("CL must be > 0.");
  if (V_  <= 0.0) throw std::invalid_argument("V must be > 0.");
}

double SingleIVBolus_1C::conc_iv_bolus_at(double dose, double time, double ke) const {
  return (dose / V_) * std::exp(-ke * time);
}

std::vector<double> SingleIVBolus_1C::build_c_at_times(double dose, std::vector<double> times) const {
        std::vector<double> c;
        c.reserve(times.size());
        
        const double ke = CL_/V_;
        
        for (std::size_t i = 0; i< times.size(); ++i) {
            double ct = conc_iv_bolus_at(dose, times[i], ke);
            c.push_back(ct);
        }

        return c;
}