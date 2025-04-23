#ifndef HBACO_h
#define HBACO_h

#include "Common.h"



/*typedef struct {
int *items_nonpris;
  int nombr_nonpris;
  int nombr;float fitness;
	int explored;
	double* f;
	int* d;
}ind;
*/
typedef struct pop_st  /* a population */
{
    int size;
    int maxsize;
    ind **ind_array;
} pop;

typedef struct genetic_op /* a genetic operator (2 point mutations)*/
{
  int p1;
  int p2;
} mut;

ind* ind_copy(ind *i);
void* chk_malloc(size_t size);
pop* create_pop(int maxsize, int dim);
ind* create_ind(int dim);
void complete_free_pop(pop *pp);
int dominates(ind *p_ind_a, ind *p_ind_b);
void init_pop(pop *SP,int size);
void random_init_pop(pop *SP, int size);
void random_init_ind(ind *x);
void evaluate(ind *x);
void init_fitness(ind *x);
void save_seed(char *s, int x);
void save_params(char *s, char *data);
void save_pop(pop *P, char* s, int t);
int max(int a, int b);
double Tchebycheff(int nb,ind *x);
void Update_Reference_Point (pop *solutions);
void compute_ind_fitness(ind *x, pop *SP);
double calcIndicatorValue(ind *p_ind_a, ind *p_ind_b);
double calcAddEpsIndicator(ind *p_ind_a, ind *p_ind_b);
double calcIndicator(ind *p_ind_a, ind *p_ind_b, int indicator);
double calcR2Indicator(ind *p_ind_a, ind *p_ind_b);
double R2Indicator_pop(pop *SP);

#endif
