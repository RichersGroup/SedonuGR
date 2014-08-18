#include "zone.h"
#include <limits>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <string>
#include <fstream>
#define NaN std::numeric_limits<real>::quiet_NaN()

zone::zone(const int dimensionality = 0){
	v.resize(dimensionality);
	rho = NaN;
	T_gas = NaN;
	Ye = NaN;
	H = NaN;
	e_rad = NaN;
	e_abs = NaN;
	l_abs = NaN;
	e_emit = NaN;
	l_emit = NaN;
	t_eabs = NaN;
    t_eemit = NaN;
    t_labs = NaN;
    t_lemit = NaN;
}

void zone::write_header(const int dimensionality, ofstream& outf){
    outf << setprecision(4);
    outf << scientific;
    outf << "# ";
    for(int i=0; i<dimensionality; i++) outf << "r[" << i << "] ";
    outf << "e_rad rho T_gas Ye t_therm t_lep" << endl;
}

void zone::write_line(const vector<double>& r, ofstream& outf) const{
    assert(r.size()>=0);
    for(int i=0; i<r.size(); i++) outf << r[i] << " ";
    outf << e_rad << " ";
    outf << rho   << " ";
    outf << T_gas << " ";
    outf << Ye    << " ";
    //outf << t_eemit << " ";
    //outf << t_eabs  << " ";
    //outf << t_lemit << " ";
    //outf << t_labs  << " ";
    outf << 1.0 / fabs(1.0/t_eabs - 1.0/t_eemit) << "\t";
    outf << 1.0 / fabs(1.0/t_labs - 1.0/t_lemit) << "\t";
    outf << endl;
}
