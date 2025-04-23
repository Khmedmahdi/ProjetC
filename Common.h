#ifndef Common_h
#define Common_h

typedef struct {

 /* int *items_nonpris;
  int *items_pris;*/
  int nombr_nonpris;
  int nombr;
  int rank;
  int *Items;
  float fitnessbest;
    float fitness;
	int explored;
	double* f;
	double* capa;
    double* v;
	int* d;


}ind;

typedef struct{
  float min;
  float max;
}range;

#endif
