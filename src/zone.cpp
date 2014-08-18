#include "zone.h"
#include <limits>
#include <cassert>
#include <string>
#include <fstream>
#define NAN std::numeric_limits<real>::quiet_NaN()

zone::zone(){
        rho = NAN;
        T_gas = NAN;
        Ye = NAN;
        H = NAN;
        e_rad = NAN;
        e_abs = NAN;
        l_abs = NAN;
        e_emit = NAN;
        l_emit = NAN;
        t_eabs = NAN;
        t_eemit = NAN;
        t_labs = NAN;
        t_lemit = NAN;
}

zone::zone(const int dimensionality){
        rho = NAN;
        T_gas = NAN;
        Ye = NAN;
        H = NAN;
        e_rad = NAN;
        e_abs = NAN;
        l_abs = NAN;
        e_emit = NAN;
        l_emit = NAN;
        t_eabs = NAN;
        t_eemit = NAN;
        t_labs = NAN;
        t_lemit = NAN;
}
