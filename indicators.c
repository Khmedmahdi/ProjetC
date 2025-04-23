#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "indicators.h"

float rho;
int dim;
double max_bound;
/*range *bounds;*/
float max_value=1000000;

/*double calcHypervolume(ind *p_ind_a, ind *p_ind_b, int dim)
/* calculates the hypervolume of that portion of the objective space that
   is dominated by individual a but not by individual b */
/*{
    double a, b, r, max;
    double volume = 0;

    r = rho * (bounds[dim - 1].max - bounds[dim - 1].min);
    max = bounds[dim - 1].min + r;

    assert(p_ind_a != NULL);
    a = p_ind_a->f[dim - 1];
    if (p_ind_b == NULL)
	b = max;
    else
	b = p_ind_b->f[dim - 1];

    assert(dim > 0);
    if (dim == 1)
    {
	if (a < b)
	    volume = (b - a) / r;
	else
	    volume = 0;
    }
    else
    {
	if (a < b)
	{
	    volume = calcHypervolume(p_ind_a, NULL, dim - 1) *
		(b - a) / r;
	    volume += calcHypervolume(p_ind_a, p_ind_b, dim - 1) *
		(max - b) / r;
	}
	else
	{
	    volume = calcHypervolume(p_ind_a, p_ind_b, dim - 1) *
		(max - b) / r;
	}
    }

    return (volume);
}
*/
/*double calcHypervolumeIndicator(ind *p_ind_a, ind *p_ind_b, int dim){
  if (dominates(p_ind_a, p_ind_b)) return -calcHypervolume(p_ind_a, p_ind_b, dim);
  else return calcHypervolume(p_ind_b, p_ind_a, dim);
}
*/
double calcAddEpsIndicator(ind *p_ind_a, ind *p_ind_b)
/* calculates the maximum epsilon value by which individual a must be
   decreased in all objectives such that individual b is weakly dominated */
{
    int i;
/*printf("max bound dans calcul indicator %lf \n ",max_bound);*/
    double eps = 0;
    eps = ((p_ind_a->v[0]/max_bound)-(p_ind_b->v[0]/max_bound));
    /*printf(" eps %lf",eps);*/
    for (i = 1; i < dim; i++)
    {
	double temp_eps;

	temp_eps = ((p_ind_a->v[i]/max_bound)-(p_ind_b->v[i]/max_bound)) ;
	/*printf(" temp_eps %lf",temp_eps);*/
	if (temp_eps > eps)
	    eps = temp_eps;
    }
/*printf(" return %lf",eps);*/
    return (eps);
}

double calcBentleyIndicator(ind *a, ind *b){
 /* 1 if a dominates b 0 otherwise
     to be used with indicator_merge=0! */
  int i;
  float res=0.0;
  for (i=0; i<dim; i++)
    if (a->f[i]<b->f[i]) res-=1.0;
    else if (a->f[i]==b->f[i]) res-=0.5;
  return res;
}

double calcFonsecaIndicator(ind *a, ind *b){
  /* 1 if a dominates b 0 otherwise
     to be used with indicator_merge=0! */
  if (dominates(a,b)) return -1; else return 0;
}

double calcDebIndicator(ind *a, ind *b){
  /* fit(a)+1 if a dominates b fit(b) otherwise
     to be used with indicator_merge=2! */
  if (dominates(a,b)) return a->fitness-1; else return max_value;
}

double calcZitzlerIndicator(ind *a, ind *b){
  /* 1 if b dominates a 0 otherwise
     to be used with indicator_merge=0! */
  if (dominates(b,a)) return 0; else return -1;
}

double calclex1Indicator(ind *a, ind *b){
  if ((a->f[0]<b->f[0]) || ((a->f[0]==b->f[0])&&(a->f[1]<b->f[1]))) return -1;
  else return 0;
}

double calclex2Indicator(ind *a, ind *b){
  if ((a->f[1]<b->f[1]) || ((a->f[1]==b->f[1])&&(a->f[0]<b->f[0]))) return -1;
  else return 0;
}

double calcIndicatorValue(ind *p_ind_a, ind *p_ind_b, int indicator, float r, int d, double b){
    rho=r;
    dim=d;
    max_bound=b;
    if (indicator == 0) return calcAddEpsIndicator(p_ind_a,p_ind_b);
/*    else if (indicator == 1) return calcHypervolumeIndicator(p_ind_a,p_ind_b,dim);*/
    else if (indicator == 2) return calcBentleyIndicator(p_ind_a,p_ind_b);
    else if (indicator == 3) return calcFonsecaIndicator(p_ind_a,p_ind_b);
    else if (indicator == 4) return calcDebIndicator(p_ind_a,p_ind_b);
    else if (indicator == 5) return calcZitzlerIndicator(p_ind_a,p_ind_b);
    else if (indicator == 6) return calclex1Indicator(p_ind_a,p_ind_b);
    else if (indicator == 7) return calclex2Indicator(p_ind_a,p_ind_b);
    else return 0;
}



