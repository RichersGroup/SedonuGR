#include "zone.h"
#include <iomanip>
#include <cmath>
#include <fstream>
#include "global_options.h"

zone::zone(const int dimensionality){
	v.resize(dimensionality);
	rho = NaN;
	T = NaN;
	Ye = NaN;
	H_com = NaN;
	e_rad = NaN;
	e_abs = NaN;
	nue_abs = NaN;
	anue_abs = NaN;
	e_emit = NaN;
	l_emit = NaN;
	t_eabs = NaN;
	t_eemit = NaN;
	t_labs = NaN;
	t_lemit = NaN;
	Q_annihil = NaN;
}
