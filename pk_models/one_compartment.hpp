#pragma once
#include <vector>

class SingleIVBolus_1C{
private:
    double CL_, V_;

public:
    SingleIVBolus_1C(double CL, double V);
    double conc_iv_bolus_at(double dose, double time, double ke) const;
    std::vector<double> build_c_at_times(double dose, std::vector<double> times) const ;
};

