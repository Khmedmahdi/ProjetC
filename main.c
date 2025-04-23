#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <conio.h>

#include "Common.h"
#include "IBMOLS.h"

#define dimension 2
#define NBITEMS 250
/*#define NBANTS 100*/
/*#define NBF 12*/
#define FREQUANCY 200
#define L 5

#define LARGE 10e50
#define FALSE 0
#define TRUE 1

#ifndef M_PI
 #define M_PI 3.14159265358979323846
#endif

double capacities[dimension];
int weights[dimension][NBITEMS];
int profits[dimension][NBITEMS];

int nf,ni,cardP;
int nombr;

int paretoIni=28000;


int inter;
int iseed;

//////////************


pop *archive=NULL;
pop *P=NULL;
pop *solutions=NULL;
/*pop *LSarchive=NULL;*/

double perturbation_rate=0.05;

int  *triI;
double U2[NBITEMS];
FILE *Wfile;
double ReferencePoint[dimension];
double vector_weight[dimension];
double max_bound;
double OBJ_Weights[dimension][10000];
int nombreLIGNE=0;
int nextLn=0;
int inv=0;

float smallValue=0.0000001;
float max_value;
double kappa=0.05;


int alpha=10;  /* number of individuals in initial population */


range *bounds;

void loadMOKP(char *s){

  FILE* source;

int i,f;

  char cl[20];

  // Opening
   source=fopen(s,"r");

  fscanf(source, " %d %d  \n", &nf,&ni);
  printf( " %d %d  \n ", nf,ni);

  for (f=0;f<nf;f++)
	{	fscanf(source, "%lf  \n ",&capacities[f]);
		/*printf( " %lf  \n \n", capacities[f]);*/

		for (i=0;i<ni;i++)
		{
			fscanf(source, " %s \n", cl); /*printf( " %s  \n \n", cl);*/
			fscanf(source, " %d  \n", &weights[f][i]); /*printf( " %d \n \n", weights[f][i]);*/
			fscanf(source, "  %d  \n ", &profits[f][i]); /*printf( " %d \n \n", profits[f][i]);*/

		}
	}
  fclose(source);
}


  /********************************/
 /* Begin memory (des)allocation */
/********************************/

void* chk_malloc(size_t size)
/* Wrapper function for malloc(). Checks for failed allocations. */
{
    void *return_value = malloc(size);
    if(return_value == NULL)
	printf("Selector: Out of memory.");
    return (return_value);
}




pop* create_pop(int maxsize, int nf)
/* Allocates memory for a population. */
{
    int i;
    pop *pp;

    assert(nf >= 0);
    assert(maxsize >= 0);

    pp = (pop*) chk_malloc(sizeof(pop));
    pp->size = 0;
    pp->maxsize = maxsize;
    pp->ind_array = (ind**) chk_malloc(maxsize * sizeof(ind*));

    for (i = 0; i < maxsize; i++)
	pp->ind_array[i] = NULL;

    return (pp);
}

ind* create_ind(int nf)
/* Allocates memory for one individual. */
{
    ind *p_ind;

    assert(nf >= 0);
    assert(ni >=0 );

    p_ind = (ind*) chk_malloc(sizeof(ind));
    /*p_ind->items_pris=(int*) chk_malloc(ni * sizeof(int));
    p_ind->items_nonpris=(int*) chk_malloc(ni * sizeof(int));*/
    p_ind->nombr_nonpris=0;
    p_ind->nombr=0;
    p_ind->rank=0;
    p_ind->fitness = -1.0;
    p_ind->fitnessbest = -1.0;
    p_ind->explored = 0;
    p_ind->f = (double*) chk_malloc(nf * sizeof(double));
    p_ind->capa = (double*) chk_malloc(nf * sizeof(double));
    p_ind->v = (double*) chk_malloc(nf * sizeof(double));
    p_ind->d = (int*) chk_malloc(ni * sizeof(int));
    p_ind->Items = (int*) chk_malloc(ni * sizeof(int));

    return (p_ind);
}


ind* ind_copy(ind *i){
/* Allocates memory for one individual. */
    ind *p_ind=NULL;
    int k;

    p_ind=create_ind(nf);

    p_ind->nombr_nonpris=i->nombr_nonpris;
    p_ind->nombr=i->nombr;

   /* for (k=0;k<i->nombr;k++)
    {p_ind->items_pris[k]=i->items_pris[k];}


    for (k=0;k<i->nombr_nonpris;k++)
    {p_ind->items_nonpris[k]=i->items_nonpris[k];}*/
    p_ind->rank = i->rank;
    p_ind->fitnessbest = i->fitnessbest;
    p_ind->fitness = i->fitness;
    p_ind->explored = i->explored;

    for (k=0;k<nf;k++)
        {p_ind->f[k]=i->f[k];
         p_ind->v[k]=i->v[k];
         p_ind->capa[k]=i->capa[k];}

    for (k=0;k<ni;k++)
        {p_ind->d[k]=i->d[k];
         p_ind->Items[k]=i->Items[k];}

    return (p_ind);
}




void free_ind(ind *p_ind)
/* Frees memory for given individual. */
{
  assert(p_ind != NULL);
  /*free(p_ind->items_pris);
  free(p_ind->items_nonpris);*/
  free(p_ind->d);
  free(p_ind->f);
  free(p_ind->capa);
  free(p_ind->v);
  free(p_ind->Items);
  free(p_ind);
}

void free_pop(pop *pp)
/* Frees memory for given population. */
{
   if (pp != NULL)
   {
      free(pp->ind_array);
      free(pp);
   }
}

void complete_free_pop(pop *pp)
/* Frees memory for given population and for all individuals in the
   population. */
{
   int i = 0;
   if (pp != NULL)
   {
      if(pp->ind_array != NULL)
      {
         for (i = 0; i < pp->size; i++)
         {
            if (pp->ind_array[i] != NULL)
            {
               free_ind(pp->ind_array[i]);
               pp->ind_array[i] = NULL;
            }
         }
/*printf( " for" );*/
         free(pp->ind_array);
/*printf( " free ind" );*/
      }

      free(pp);
   }
   /*printf( " free " );*/
}



int dominates(ind *p_ind_a, ind *p_ind_b)
/* Determines if one individual dominates another.
   Minimizing fitness values. */
{
    int i;
    int a_is_worse = 0;
    int equal = 1;

     for (i = 0; i < nf && !a_is_worse; i++)
     {
	 a_is_worse = p_ind_a->f[i] > p_ind_b->f[i];
          equal = (p_ind_a->f[i] == p_ind_b->f[i]) && equal;
     }

     return (!equal && !a_is_worse);
}

int non_dominated(ind *p_ind_a, ind *p_ind_b)
/* Determines if one individual is non dominated according to another one.
   Minimizing fitness values.
   -1 dominated, 0 equal, 1 non-dominated */
{
    int i;
    int a_is_good = -1;
    int equal = 1;

     for (i = 0; i < nf; i++)
     {
       if (p_ind_a->f[i] > p_ind_b->f[i]) a_is_good=1;
         equal = (p_ind_a->f[i] == p_ind_b->f[i]) && equal;
     }
     if (equal) return 0;
     return a_is_good;
}


int otherResult(char* file_name,  double mpareto2[][nf]) {
  FILE* fd;

  int f;
  int cardP2 =0;

   if ( (fd=fopen(file_name, "r"))==NULL)
    {	printf("ERREUR: Verifiez le nom de fichier %s", file_name); return 0;	}

  do{

  for (f=0;f<nf;f++)
	{

	  fscanf(fd, "%lf  ",&mpareto2[cardP2][f]);

	}
  cardP2++;
  } while(!feof(fd));

  return cardP2;

}
// fonction de calcul de performance
double Cmesure(double mpareto1[][nf], double mpareto2[][nf], int cardP1, int cardP2)
{
	int i=0,efficace=1,existe=0,f, cdom,cegal,c=0,j;

		printf("carp1 %d cardp2 %d \n ", cardP1, cardP2);
				while (i<cardP2)
				{j=0;
					while (j<cardP1)
					{
					f=0;
					cegal=0;
					cdom=0; //pour detecter si sol cour est dominée
					while (f<nf)
					{
						if (mpareto2[i][f]<mpareto1[j][f])
						{
						cdom++;
						}
						if (mpareto2[i][f]==mpareto1[j][f])
						{
						cegal++;
						}


						f++;
					}
					j++;
					}
					if ((cdom>0)&&((cdom+cegal)==nf))//la sol cou est dominé
						c++;
					i++;
				}

				return ((double) c/ (double) cardP2);

}






/******************************/
 /* WEIGHTS functions */
/******************************/


void initfile_weights_log()
{
    int F=FREQUANCY/4,i,j,k;
    double lamda1,lamda2,lamda3,lamda4,tmp1,tmp2,tmp3;
    int T1,T2,M1,Zd;
    double C=exp(1);

  Wfile = fopen( "Weights.txt", "a+" );
  fflush(stdout);

if (nf==2)
{
for(i=0;i<=F;i++)
{

lamda1=log((4.0*i/FREQUANCY*C)+(cos(2.0*M_PI*i/FREQUANCY)));
lamda2=1.0-lamda1;

fprintf(Wfile,"%lf ",lamda1);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda2);
fflush(stdout);
fprintf(Wfile,"\n");
}
/*w1=1.0; w2=0.0;
fprintf(Wfile,"%lf ",w1);
fflush(stdout);
fprintf(Wfile,"%lf ",w2);
fflush(stdout);*/
}

if (nf==3)
{

for(i=0;i<F;i++)
{
lamda1=log((4.0*i/FREQUANCY*C)+(cos(2.0*M_PI*i/FREQUANCY)));

for(j=0;j<F;j++)
{

lamda2=(1.0-lamda1)*log((4.0*j/FREQUANCY*C)+(cos(2.0*M_PI*j/FREQUANCY)));

lamda3=1.0-lamda1-lamda2;

fprintf(Wfile,"%lf ",lamda1);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda2);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda3);
fflush(stdout);
fprintf(Wfile,"\n");
}
}
lamda1=1.0; lamda2=0.0; lamda3=0.0;
fprintf(Wfile,"%lf ",lamda1);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda2);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda3);
fflush(stdout);
}
if (nf==4)
{
for(i=0;i<F;i++)
{
lamda1=log((4.0*i/FREQUANCY*C)+(cos(2.0*M_PI*i/FREQUANCY)));

for(j=0;j<F;j++)
{
lamda2=(1.0-lamda1)*log((4.0*j/FREQUANCY*C)+(cos(2.0*M_PI*j/FREQUANCY)));

for(k=0;k<F;k++)
{
lamda3=(1.0-lamda1-lamda2)*log((4.0*k/FREQUANCY*C)+(cos(2.0*M_PI*k/FREQUANCY)));
lamda4=1.0-lamda1-lamda2-lamda3;

fprintf(Wfile,"%lf ",lamda1);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda2);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda3);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda4);
fflush(stdout);
fprintf(Wfile,"\n");
}
}
}
lamda1=1.0; lamda2=0.0; lamda3=0.0; lamda4=0.0;
fprintf(Wfile,"%lf ",lamda1);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda2);
fflush(stdout);
fprintf(Wfile,"%f ",lamda3);
fflush(stdout);
fprintf(Wfile,"%lf ",lamda4);
fflush(stdout);
}
fclose(Wfile);
}



void read_weights_file (char *s)
 {
     int i,j,caractereLu;
      char cl[20];

      Wfile = fopen( s, "r" );
/*
    do
    {
        caractereLu = fgetc(Wfile);
        printf("%c", caractereLu);
        if (caractereLu == '\n')
            nombreLIGNE++;
    } while(caractereLu != EOF);

       rewind(Wfile);
*/
if (nf==2) {nombreLIGNE=FREQUANCY/4;}
if (nf==3) {nombreLIGNE=(FREQUANCY/4)*(FREQUANCY/4);}
if (nf==4) {nombreLIGNE=FREQUANCY/4*(FREQUANCY/4)*(FREQUANCY/4);}

printf("nombre ligne %d \n", nombreLIGNE);

		for (i=0;i<nombreLIGNE+1;i++)
		{
		    for (j=0;j<nf;j++)
	     {
			fscanf(Wfile, " %lf ", &OBJ_Weights[j][i]); /*printf( " obj_W %lf ", OBJ_Weights[j][i]);*/
		 }
          /*printf("%d",i); printf("\n");*/
	    }

  fclose(Wfile);

 }



/*void dynamic_weight_allpop()
{
    int i;

if (inv==0){
   for(i=0;i<dimension;i++)
{
    vector_weight[i]=OBJ_Weights[i][nextLn];

}
 nextLn++;

}
if (inv==1) {

         for(i=0;i<dimension;i++)
{
    vector_weight[i]=OBJ_Weights[i][nextLn];

}

nextLn--;

}

if (nextLn==nombreLIGNE)
{
    inv=1;
}

}*/


void dynamic_weight_allpop()
{
    int i;

   for(i=0;i<dimension;i++)
{
    vector_weight[i]=OBJ_Weights[i][nextLn];

}


if(nextLn==nombreLIGNE)
        {nextLn=0;}
     else
         nextLn++;


}



void choose_weight ()
{

int i;
/*random_normalisated_weights();*/
/*random_weights();*/
/*dynamic_weight (SP,Sarchive);*/
dynamic_weight_allpop();
/*for(i=0;i<nf;i++){
    printf("poids %lf",vector_weight[i]);}*/

}



 /******************************/
 /* MOACO functions */
/******************************/
void seed(unsigned int seed)
 {iseed= seed; }



double monPow(double x, int y){
  double p;
  if (y==0) return 1;
  else if (y==1) return x;
  else if (y==2) return x*x;
  else if (y==3) return x*x*x;
  else if (y==4) {x *= x; return x*x;}
  else if (y==5) {p = x*x; return p*p*x;}
  else if (y==6){x = x*x*x; return x*x;}
  else {
    if (y%2==0){p = monPow(x,y/2); return p*p;}
    else{p = monPow(x,(y-1)/2); return p*p*x;};
  };
}



/* Generate a random double. */
double drand(double range)
{
     double j;
     j=(range * (double) rand() / (RAND_MAX + 1.0));
     return (j);
}


int irand(int range)
/* Generate a random integer. */
{
    int j;
    j=(int) ((double)range * (double) rand() / (RAND_MAX+1.0));
    return (j);

}

int rand_a_b(int a, int b){
    return rand()%(b-a) +a;
}

void init_pop(pop *SP,int size){
  int i;
  SP->size=size;
  random_init_pop(SP,size);
  /* init explored to 0! */
  for (i=0; i<SP->size; i++)

    SP->ind_array[i]->explored=0;
  /* We could add different scenarios for initialization */
  /* But, we have to add a constant which corresponds to the init algorithm */
}


void random_init_pop(pop *SP, int size){
/* Random initialisation of a population of permutations */
  int i;
SP->size=size;
  for (i=0; i<SP->size; i++){
    SP->ind_array[i]=create_ind(nf); /*printf("create ind");*/
    random_init_ind(SP->ind_array[i]); /*printf("random init");*/
    evaluate(SP->ind_array[i]); /*printf("init pop \n");*/
  }

}

void random_init_ind(ind *x){
  int j,r,tmp;

/* Random initialisation of a permutation */
  for (j=0; j<ni; j++)
      x->d[j]=j;
    for (j=0; j<ni; j++){
      r=irand(ni);
      tmp=x->d[r];
      x->d[r]=x->d[j];
      x->d[j]=tmp;
    }
}

/*void evaluate(ind *x)
{
    int j,i,l,k,faisable;

x->nombr=0; x->nombr_nonpris=0;
for(j=0;j<nf;j++)
{
x->capa[j]=0; x->f[j]=0;
}

for (j=0; j<ni; j++)
{
l=0; faisable=1;
do{
if (x->capa[l]+ weights[l][x->d[j]] > capacities[l])
faisable=0;
l++;
}while((l<nf)&&(faisable==1));

if (faisable==1)
{
for (k=0;k<nf;k++) {
x->capa[k] = x->capa[k]+weights[k][x->d[j]];
x->f[k] = x->f[k]+profits[k][x->d[j]];
}

x->items_pris[x->nombr]=x->d[j];

x->nombr++;
}
else if (faisable==0)
{
x->items_nonpris[x->nombr_nonpris]=x->d[j];

x->nombr_nonpris++;
}


}

}

*/
void evaluate(ind *x)
{
    int j,i,l,k,faisable;

x->nombr=0; x->nombr_nonpris=0;
for(j=0;j<nf;j++)
{
x->capa[j]=0; x->f[j]=0;
}

for (j=0; j<ni; j++)
{
l=0; faisable=1;
do{
if (x->capa[l]+ weights[l][x->d[j]] > capacities[l])
faisable=0;
l++;
}while((l<nf)&&(faisable==1));

if (faisable==1)
{
for (k=0;k<nf;k++) {
x->capa[k] = x->capa[k]+weights[k][x->d[j]];
x->f[k] = x->f[k]+profits[k][x->d[j]];
}

x->Items[x->d[j]]=1;
x->nombr++;
}
else if (faisable==0)
{
x->Items[x->d[j]]=0;

x->nombr_nonpris++;
}


}

}



void perturbation(ind *x){
  int i,j,t,k,l,faisable, pos,objet,objet1, pos1;
  int bruit_rate;

bruit_rate=perturbation_rate*x->nombr;
round(bruit_rate);
/*printf("bruit rate %d",bruit_rate );*/

for (j=0;j<bruit_rate;j++){

do{
     objet=irand(NBITEMS);
}while(x->Items[objet]==0);

x->Items[objet]=0;

x->nombr--;
x->nombr_nonpris++;

 for(k=0;k<nf;k++)
   {
  x->f[k]=x->f[k]-profits[k][objet];
  x->capa[k]=x->capa[k]-weights[k][objet];
   }

}



for(j=0;j<x->nombr_nonpris;j++)
{

do{
     objet1=irand(NBITEMS);
}while(x->Items[objet1]==1);

   l=0; faisable=1;
do{
if ((x->capa[l]+weights[l][objet1])>capacities[l])
faisable=0;
l++;
}while((l<nf)&&(faisable==1));
/*printf("faisable %d",faisable );*/

if  (faisable==1) {

x->Items[objet1]=1;

x->nombr++;
x->nombr_nonpris--;

for (t=0;t<nf;t++)
{
x->f[t]=x->f[t]+profits[t][objet1];
x->capa[t]=x->capa[t]+weights[t][objet1];
}

}

}

   /* for (t=0;t<nf;t++){
printf("prof %lf",Sarchive->ind_array[i]->f[t]);
     }
printf("nombr %d",Sarchive->ind_array[i]->nombr);
printf("nombr pris %d",Sarchive->ind_array[i]->nombr_nonpris);*/
     /*for(j=0;j<Sarchive->ind_array[i]->nombr;j++)
{printf("pris %d",Sarchive->ind_array[i]->items_pris[j]);}
for(j=0;j<Sarchive->ind_array[i]->nombr_nonpris;j++)
{printf("non pris %d",Sarchive->ind_array[i]->items_nonpris[j]);}*/


}



void P_init_pop(pop *SP,pop *Sarchive,int alpha){
  int i,x,tmp;
  int t=max(alpha,Sarchive->size);
  int shuffle[t];

  SP->size=alpha;
  for (i=0;i<t;i++) shuffle[i]=i;
  for (i=0;i<t;i++){
    x=irand(alpha);
    tmp=shuffle[i];
    shuffle[i]=shuffle[x];
    shuffle[x]=tmp;
  }

  if (Sarchive->size>alpha){ /* Pareto set sufficient to create all the solutions? */
        printf("size P %d",Sarchive->size );
     for (i=0;i<alpha;i++){
       //P->ind_array[i]=create_ind(dim);
       SP->ind_array[i]=ind_copy(Sarchive->ind_array[shuffle[i]]);
      /* perturbation(SP->ind_array[i]);*/
       /*evaluate_perturbation(SP->ind_array[i]);*/

     }
  }
  else{
         /* insert some random individuals */
    for (i=0;i<alpha;i++){
      //P->ind_array[i]=create_ind(dim);
      if (shuffle[i]<Sarchive->size){

	SP->ind_array[i]=ind_copy(Sarchive->ind_array[shuffle[i]]);
	/*perturbation(SP->ind_array[i]);*/
      }
      else {

	SP->ind_array[i]=create_ind(nf);
	random_init_ind(SP->ind_array[i]);
      evaluate(SP->ind_array[i]);
      }

    }
  }

}





int extractPtoArchive(pop *P, pop *archive){
  int i,j,dom;
  int t=archive->size+P->size;
  pop *archiveAndP;
  int convergence_rate=0;


  archiveAndP=create_pop(t,nf);
  for (i=0;i<archive->size;i++){
        /*printf("A %d",Sarchive->size);*/
    archiveAndP->ind_array[i]=archive->ind_array[i]; /* do not copy */
    //archive->ind_array[i]=NULL;  /*!*/
  }
  /*printf("here");*/
  //archive->size=0;
  for (i=0;i<P->size;i++){
    /*printf("P %d",SP->size);*/
    archiveAndP->ind_array[i+archive->size]=ind_copy(P->ind_array[i]);
  }
  archiveAndP->size=t;
  archive->size=0;

  /* initialize the future size for archive
		     anyway the solutions are lost (-> in archiveAndP) */
		      /*printf("encore here");*/
  for (i=0;i<t;i++){
    for (j=0;j<t;j++){
      if (i!=j) {
	dom=non_dominated(archiveAndP->ind_array[i],archiveAndP->ind_array[j]);
	if (dom==-1 || (dom==0 && i>j)) j=t+1;
      }
    }
    if(j==t) {
      archive->ind_array[archive->size++]=ind_copy(archiveAndP->ind_array[i]); //save the current
      if (i>=t-P->size) convergence_rate++;
    }
  }
/*printf("%d",Sarchive->size);*/
   /*printf(" pareto maj %d",verif_MAJ);*/
  complete_free_pop(archiveAndP);
/*printf("encore here");*/

  return convergence_rate;
}





/*void determiner (pop *solutions, int size){
int its[NBITEMS]; int prof[dimension];
int nit, existe,s,i,j,k;
nit=NBITEMS;
solutions->size=size;
for(i=0;i<ni;i++){its[i]=i;}

for(k=0;k<solutions->size;k++){
s=0;
solutions->ind_array[k]->nombr_nonpris =nit-solutions->ind_array[k]->nombr;
    for(i=0;i<ni;i++){

    for(j=0;j<solutions->ind_array[k]->nombr;j++) {

         if(its[i]!=solutions->ind_array[k]->d[j]){existe=0;}

         else if (its[i]==solutions->ind_array[k]->d[j]){
            existe=-1; break;}

    }
        if (existe==0)
      {solutions->ind_array[k]->items_nonpris[s]=its[i]; s++;}
    }

}
}

*/




/******************************/

 /* indicator_local search functions*/

/******************************/


void calcul_weight (pop *SP, int size)
{
    int i,j;


        for(i=0;i<SP->size;i++)
        {
            for(j=0;j<nf;j++)
        {
          SP->ind_array[i]->v[j]=SP->ind_array[i]->f[j]*vector_weight[j];
           printf("poids %lf",vector_weight[j]);

        }
        printf("\n");

        }
}


void calcMaxbound (pop* SP,int size)
{
    SP->size=size;
    int i, j;

	max_bound = SP->ind_array[0]->v[0];
/*printf("avant %lf ",max_bound);*/
for (i = 0; i < SP->size; i++){
	for (j = 0; j < nf; j++){
	   /* printf("en cours %lf \n",SP->ind_array[i]->v[j]);*/
	    if (max_bound < SP->ind_array[i]->v[j])
		{max_bound  = SP->ind_array[i]->v[j];}
	}
    }
/*printf("apres %lf \n",max_bound);*/
}


int delete_fitness(ind *x,double I){

 x->fitness+=exp(-I/kappa);

  return 0;
}


double update_fitness_return(double f,double I){

  return f-exp(-I/kappa);

}

void compute_all_fitness(pop *SP){
  int i;
  for (i=0;i<SP->size;i++)
  compute_ind_fitness(SP->ind_array[i],SP);
}

void init_fitness(ind *x){
  /*to complete...*/
  x->fitness=0.0;
}

void update_fitness(ind *x,double I){

  x->fitness-=exp(-I/kappa); /*printf(" \n indicator %lf", I); printf(" \n fitness %lf", x->fitness);*/

}

/* compute the fitness of x according to the population P */
void compute_ind_fitness(ind *x, pop *SP){
  int j;
  init_fitness(x);
  for (j=0;j<SP->size;j++){
    if (SP->ind_array[j]!=x){
      update_fitness(x,calcAddEpsIndicator(SP->ind_array[j],x));

    }
/*printf(" \n fitness %lf", SP->ind_array[j]);*/
  }
}


double calcAddEpsIndicator(ind *p_ind_a, ind *p_ind_b)
/* calculates the maximum epsilon value by which individual a must be
   decreased in all objectives such that individual b is weakly dominated */
{
    int i;
/*printf("max bound dans calcul indicator %lf \n ",max_bound);*/
    double eps = 0;
    eps = ((p_ind_a->v[0]/max_bound)-(p_ind_b->v[0]/max_bound));
    /*printf(" eps %lf",eps);*/
    for (i = 1; i < nf; i++)
    {
	double temp_eps;

	temp_eps = ((p_ind_a->v[i]/max_bound)-(p_ind_b->v[i]/max_bound)) ;
	/*printf(" temp_eps %lf",temp_eps);*/
	if (temp_eps > eps)
	    eps = temp_eps;
    }
/*printf("  eps return %lf",eps);*/
    return (eps);
}


/* compute the fitness of a new individual x     */
/* compute the updated fitness of the population */
/* if x is the worst individual, then return -1 */
/* else delte the worst individual               */
/* and update the population fitness             */
/* return the indice of the replaced solution    */
int compute_fitness_and_select(pop *SP, ind* x, int size){
  int i,j;
  int worst=-1;
  double worst_fit,fit_tmp;
SP->size=size;
  /* fitness for x */
  x->fitness=0;
  compute_ind_fitness(x,SP);
  worst_fit=x->fitness;
  /*printf(" \n fitness %f",x->fitness);*/

  /* update the other fitness according to x */
  for (i=0;i<SP->size;i++){
    fit_tmp=update_fitness_return(SP->ind_array[i]->fitness,calcAddEpsIndicator(x,SP->ind_array[i]));
    if (fit_tmp>worst_fit) {
      worst=i;
      worst_fit=fit_tmp;
    }
   //printf(" \n fitness %f", SP->ind_array[i]->fitness);
  //  for(j=0;j<nf;j++)
      //  {
       // printf(" %f",SP->ind_array[i]->f[j]);
         /*printf("vector,%lf",vector_weight[j]);*/
   // }
  }
//  printf(" \n ");
  /* save x fitness */
  fit_tmp=x->fitness;

  /* did we keep x into the population? x */
  if (worst==-1) return -1; /* No */
  else{  /* Yes: remove the first ones, and update fitness, then replace worst */
    /* fitness: P[i]=P[i]+x-worst / x=x-worst */
    for (i=0;i<SP->size;i++){
      if (delete_fitness(SP->ind_array[i],calcAddEpsIndicator(SP->ind_array[worst],SP->ind_array[i])))
	compute_ind_fitness(SP->ind_array[i],SP); /* special case... we have to re-compute all */
      update_fitness(SP->ind_array[i],calcAddEpsIndicator(x,SP->ind_array[i]));
    }
    if (delete_fitness(x,calcAddEpsIndicator(SP->ind_array[worst],x)))
      compute_ind_fitness(x,SP); /* special case... we have to re-compute all */
    free_ind(SP->ind_array[worst]);
    SP->ind_array[worst]=ind_copy(x);
    if (fit_tmp-worst_fit>smallValue) /* to avoid approximation errors */
      return worst;
    else return -1; /* perhaps the difference is due to approximation errors */
  }
}



void Indicator_local_search1(pop *SP,pop *Sarchive,int size){

ind *y;
ind *x;
int i,j,r,t,k,l,v,sol,mino,mp,maxp,consistant,pos, stop,convergence, ii;
int tmp_pris, tmp_nonpris, remplace[L], taille,feasible;
SP->size=size;

    extractPtoArchive(SP,Sarchive);
//printf("do");//
do{
  convergence=0;

  for (i=0;i<SP->size;i++) {
      if (!SP->ind_array[i]->explored){

      x=ind_copy(SP->ind_array[i]);


for(j=0;j<x->nombr;j++)
{


for(l=0;l<L;l++)
{remplace[l]=0;}

do{

     mino=irand(NBITEMS);

}while(x->Items[mino]==0);

	x->Items[mino]=0;
	x->nombr--;
	x->nombr_nonpris++;

	for (r=0;r<nf;r++)// maj capacit
    {x->capa[r]=x->capa[r]-weights[r][mino];
     x->f[r]=x->f[r]-profits[r][mino];}


// add some items
int IM=0;
stop=1;
taille=0;
do{

do{

     maxp=irand(NBITEMS);

}while(x->Items[maxp]==1);

if (maxp!=mino){
consistant=1;r=0;
do
{ if (x->capa[r]+weights[r][maxp]>capacities[r])
consistant=0;

r++;
}while ((r<nf)&&(consistant==1));


if (consistant==1)
{
feasible=1;r=0;
do
{ if (maxp==remplace[r])
feasible=0;
r++;
}while ((r<taille)&&(feasible));


 if (feasible==1) {

  remplace[taille]=maxp;
  taille++;

    x->Items[maxp]=1;
	x->nombr_nonpris--;
    x->nombr++;

	for (r=0;r<nf;r++)// maj capacit
    {x->capa[r]=x->capa[r]+weights[r][maxp];
     x->f[r]=x->f[r]+profits[r][maxp];}




 }

}

}

IM++;
}while(IM < L );
int tv;
/*free(maxI);*/
for(tv=0;tv<nf;tv++)
{x->v[tv]=x->f[tv]*vector_weight[tv];}

    calcMaxbound(SP,SP->size);
     sol=compute_fitness_and_select(SP,x,SP->size);
         /*printf("indice %d\n", sol);*/
       if (sol!=-1){ // keep it if good //
	    j=x->nombr+1; // a enlever //
         if (sol>i){
       // do not explore the neighbors of the new solution immediately //
		y=SP->ind_array[i+1];
		SP->ind_array[i+1]=SP->ind_array[sol];
		SP->ind_array[sol]=y;
		i++;

	      }
       }
       else if (sol==-1)
       {



    x->Items[mino]=1;
    x->nombr_nonpris--;
    x->nombr++;


  for (r=0;r<nf;r++)// maj capacit
    {x->capa[r]=x->capa[r]+weights[r][mino];
     x->f[r]=x->f[r]+profits[r][mino];}



if (taille>=1) {
for(r=0;r<taille;r++)
 {

    x->Items[remplace[r]]=0;
   x->nombr--;
    x->nombr_nonpris++;

 for (t=0;t<nf;t++)// maj capacit
    {x->capa[t]=x->capa[t]-weights[t][remplace[r]];
     x->f[t]=x->f[t]-profits[t][remplace[r]];
     x->v[t]=x->f[t]*vector_weight[t];}



}
}

       }



}

tmp_pris=x->nombr;
tmp_nonpris=x->nombr_nonpris;
free_ind(x);
}

if (j==tmp_pris) SP->ind_array[i]->explored=1;
  }
j=convergence;
convergence=extractPtoArchive(SP,Sarchive);

}while(convergence);

}



int main(int argc, char* argv[])
{
	double  duration, durmoy=0;
	clock_t start, finish;
    int it, NBL; int k, i,j;

double init_time;
double run_time;


   for(k=1;k<=10;k++)
 {

    duration=0;
    double mpareto1[28000][dimension];

init_time=0.0;
run_time=2400.0;
NBL=100;

nombreLIGNE=0; nextLn=0; inv=0;
inter=k;

	FILE *fpareto, *ftime;
	fpareto = fopen( "2502_Resulats.txt", "a+" );
	fprintf(fpareto,"\n");
	fflush(stdout);


seed(inter);
bounds = (range*) chk_malloc(nf * sizeof(range));

P=create_pop(paretoIni,nf);



  /*  int cardP2=0,j,cardP3=0;
	double pareto2[8000][dimension];//, pareto3[100];*/

	loadMOKP("250.2.txt");



/*initfile_weights_log(); printf("FICHIER W");*/
read_weights_file("Weights_2obj_FQ200.txt");


it=0;

do{

solutions=create_pop(alpha,nf);
archive=create_pop(paretoIni,nf);

//for nieghobr structure//

 choose_weight();

//end//



P_init_pop(solutions,P,alpha);


/*init_pop(solutions,alpha);*/

extractPtoArchive(solutions,P);

//indicator_local_search function //

/* epsilon indicator*/

calcul_weight(solutions,alpha);
calcMaxbound(solutions,alpha);
compute_all_fitness(solutions);

//end indicator_local_search functions //


start=clock();

Indicator_local_search1(solutions,archive,alpha);
extractPtoArchive(archive,P);

finish = clock(); /*printf( " time" );*/
duration = (double)(finish - start) / CLOCKS_PER_SEC;
it++;
printf("\n %d ligne",it);
printf("DURATION %f",duration);

init_time+=duration;

complete_free_pop(solutions);
complete_free_pop(archive);


	/*determiner(P);  printf( " determiner" );*/

printf("time %f",init_time);
}while(it <NBL);


	for(i=0;i<P->size;i++)
	{for(j=0;j<nf;j++)
		{fprintf(fpareto,"%f ",P->ind_array[i]->f[j]);
         /*fprintf(fpareto,"%f ",P->ind_array[i]->capa[j]);*/
		fflush(stdout);
		//mpareto1[i][j]=P->ind_array[j]->f[j];
		}
		fprintf(fpareto,"\n");

	}


/*fprintf(fpareto,"ACOeps \n");
cardP2= otherResult("reslt ACOEps\\ACOEps750.4\\750.4.30.txt",pareto2);
	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));*/

	/*fprintf(fpareto,"ACOHD \n");
cardP2= otherResult("reslt ACOHD\\ACOHD250.2\\250.2.30.txt",pareto2);

	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));

fprintf(fpareto,"IBEAeps \n");
cardP2= otherResult("C:\\Users\\imen\\Desktop\\projet\\ibea_knapsack\\knapsack\\bin\\Debug\\IBEAeps\\IBEA250.2\\knapsack_output30.txt",pareto2);
	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));

	fprintf(fpareto,"IBEAhd \n");
cardP2= otherResult("C:\\Users\\imen\\Desktop\\projet\\ibea_knapsack\\knapsack\\bin\\Debug\\IBEAHD\\IBEA250.2\\knapsack_output30.txt",pareto2);
	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));

	fprintf(fpareto,"spea2 \n");
cardP2= otherResult("C:\\Users\\imen\\Desktop\\projet\\SPEA2_knapsack\\knapsack\\bin\\Debug\\SPEA2250.2\\knapsack_output30.txt",pareto2);

	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesurerandom_normalisated_weights()(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));

	fprintf(fpareto,"HypE \n");
cardP2= otherResult("C:\\Users\\imen\\Desktop\\projet\\HypE_knapsack\\knapsack\\bin\\Debug\\HypE250.2\\knapsack_output30.txt",pareto2);
	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));

	/*fprintf(fpareto,"mACO \n");
cardP2= otherResult("result500.3\\result500.3.1.txt",pareto2);
	fprintf(fpareto,"C Mesure 1 %lf \n", Cmesure(mpareto1, pareto2, cardP, cardP2));
	fprintf(fpareto,"C Mesure 2 %lf \n", Cmesure(pareto2, mpareto1, cardP2, cardP));*/


	complete_free_pop(P);
  /* complete_free_pop(pp_offspring);
complete_free_pop(pp_Parents);*/


/*fprintf(fpareto,"Moyenne temps CPU %f \n", durmoy);*/

	fclose(fpareto);
/*printf( " \n before free" );*/

/*free(bounds); printf( " \nfree bounds" );*/



}
return(0);

}


int max(int a, int b){
  if (a>b) return a; else return b;
}
