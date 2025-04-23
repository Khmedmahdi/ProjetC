#ifndef indicators_h
#define indicators_h

#include "Common.h"
//#include "IBMOLS.h"

double calcIndicatorValue(ind *p_ind_a, ind *p_ind_b, int indicator, float rho, int dim, double max_bound);

#endif
