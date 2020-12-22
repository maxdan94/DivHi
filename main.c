#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <strings.h>//to use "bzero"
#include <time.h>//to estimate the runing time

#include "bisection.h"
#include "struct.h"

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

//compute the maximum of three unsigned long
inline unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the edgelist from file
adjlist* readedgelist(char* input){
	unsigned long e1=NLINKS,i;
	adjlist *g=malloc(sizeof(adjlist));
	FILE *file;

	g->n=0;
	g->e=0;
	file=fopen(input,"r");
	g->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%lu %lu", &(g->edges[g->e].s), &(g->edges[g->e].t))==2) {
		g->n=max3(g->n,g->edges[g->e].s,g->edges[g->e].t);
		if (++(g->e)==e1) {
			e1+=NLINKS;
			g->edges=realloc(g->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->edges=realloc(g->edges,g->e*sizeof(edge));

	g->n++;
	g->nodes=malloc(g->n*sizeof(unsigned long));
	for (i=0;i<g->n;i++){
		g->nodes[i]=i;
	}
	return g;
}

//building the adjacency matrix
void mkadjlist(adjlist* g){
	unsigned long i,u,v,n;
	unsigned long *new,*d;

	n=0;
	g->map=malloc(g->n*sizeof(unsigned long));
	for (i=0;i<g->n;i++) {
		g->map[i]=g->nodes[i];
		n=(n>g->nodes[i])?n:g->nodes[i];
	}
	n++;

	new=malloc(n*sizeof(unsigned long));
	for (i=0;i<g->n;i++) {
		new[g->nodes[i]]=i;
		g->nodes[i]=i;
	}

	for (i=0;i<g->e;i++) {
		g->edges[i].s=new[g->edges[i].s];
		g->edges[i].t=new[g->edges[i].t];
	}

	bzero(new,g->n*sizeof(unsigned long));
	d=new;

	for (i=0;i<g->e;i++) {
		d[g->edges[i].s]++;
		d[g->edges[i].t]++;
	}

	g->cd=malloc((g->n+1)*sizeof(unsigned long));
	g->cd[0]=0;
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+d[i-1];
		d[i-1]=0;
	}

	g->adj=malloc(2*g->e*sizeof(unsigned long));

	for (i=0;i<g->e;i++) {
		u=g->edges[i].s;
		v=g->edges[i].t;
		g->adj[ g->cd[u] + d[u]++ ]=v;
		g->adj[ g->cd[v] + d[v]++ ]=u;
	}

	free(new);
}

//freeing memory
void free_adjlist(adjlist *g){
	free(g->nodes);
	free(g->edges);
	free(g->cd);
	free(g->adj);
	free(g->map);
	free(g);
}

table* alloctab(){
	table* tab=malloc(sizeof(table));
	tab->n=0;
	tab->nmax=NLINKS;
	tab->t=malloc(sizeof(unsigned));
	return tab;
}

//Make the 2 subgraphs of graph g using the labels "lab"
clusters mkkids(adjlist* g, bool* lab){
	unsigned long i,e1=g->e/2+1,e2=g->e/2+1;

	clusters clust;
	clust.n=2;
	clust.sg=malloc(clust.n*sizeof(adjlist*));
	clust.sg[0]=malloc(sizeof(adjlist));
	clust.sg[1]=malloc(sizeof(adjlist));
	clust.sg[0]->n=0;
	clust.sg[1]->n=0;
	clust.sg[0]->nodes=malloc(g->n*sizeof(unsigned long));
	clust.sg[1]->nodes=malloc(g->n*sizeof(unsigned long));
	for (i=0;i<g->n;i++){
		if (lab[i]==0){
			clust.sg[0]->nodes[(clust.sg[0]->n)++]=i;
		}
		else{
			clust.sg[1]->nodes[(clust.sg[1]->n)++]=i;
		}
	}
	clust.sg[0]->nodes=realloc(clust.sg[0]->nodes,clust.sg[0]->n*sizeof(unsigned long));
	clust.sg[1]->nodes=realloc(clust.sg[1]->nodes,clust.sg[1]->n*sizeof(unsigned long));

	clust.sg[0]->edges=malloc(e1*sizeof(edge));
	clust.sg[1]->edges=malloc(e2*sizeof(edge));
	
	clust.sg[0]->e=0;
	clust.sg[1]->e=0;
	for (i=0;i<g->e;i++) {
		if (lab[g->edges[i].s]==lab[g->edges[i].t]){
			if (lab[g->edges[i].s]==0){
				clust.sg[0]->edges[clust.sg[0]->e]=g->edges[i];
				if (++(clust.sg[0]->e)==e1) {
					e1*=2;
					clust.sg[0]->edges=realloc(clust.sg[0]->edges,e1*sizeof(edge));
				}
			}
			else{
				clust.sg[1]->edges[clust.sg[1]->e]=g->edges[i];
				if (++(clust.sg[1]->e)==e2) {
					e2*=2;
					clust.sg[1]->edges=realloc(clust.sg[1]->edges,e2*sizeof(edge));
				}
			}
		}
	}
	clust.sg[0]->edges=realloc(clust.sg[0]->edges,clust.sg[0]->e*sizeof(edge));
	clust.sg[1]->edges=realloc(clust.sg[1]->edges,clust.sg[1]->e*sizeof(edge));
	mkadjlist(clust.sg[0]);
	for (i=0;i<clust.sg[0]->n;i++){
		clust.sg[0]->map[i]=g->map[clust.sg[0]->map[i]];
	}
	mkadjlist(clust.sg[1]);
	for (i=0;i<clust.sg[1]->n;i++){
		clust.sg[1]->map[i]=g->map[clust.sg[1]->map[i]];
	}
	return clust;
}

bool* severaltimes(bisection bisec,unsigned k,adjlist* g){
	static double totcut=0;
	double val1,val2;
	bool *lab1=bisec(g,&val1),*lab2;
	unsigned i;
	for (i=0;i<k;i++){
		lab2=bisec(g,&val2);
		if (val2<val1){
			val1=val2;
			free(lab1);
			lab1=lab2;
		}
		else{
			free(lab2);
		}
	}
	totcut+=val1;
	printf("sum of cuts = %lf\n",totcut);
	return lab1;
}

//recursive function
void recurs(bisection bisec, adjlist* g, unsigned h, FILE* file){
	clusters clust;
	bool* lab;
	unsigned long i;
	if (g->e>0){
		lab=severaltimes(bisec,5,g);////////////////
		clust=mkkids(g,lab);
		fprintf(file,"%u 2\n",h);
		free(lab);
		free_adjlist(g);
		recurs(bisec, clust.sg[0], h+1, file);
		recurs(bisec, clust.sg[1], h+1, file);
	}
	else{
		fprintf(file,"%u 1 %lu",h,g->n);
		for (i=0;i<g->n;i++){
			fprintf(file," %lu",g->map[i]);
		}
		fprintf(file,"\n");
		free_adjlist(g);
	}
	
}

//main function
int main(int argc,char** argv){
	adjlist* g;
	bisection bisec;
	table* hier;
	FILE* file;
	time_t 	t1=time(NULL),t2;
	srand(time(NULL));

	if (argc==3)
		bisec=choosebisection("1");
	else if (argc==4)
		bisec=choosebisection(argv[3]);
	else{
		exit(1);
	}

	printf("Reading edgelist from file %s\n",argv[1]);
	g=readedgelist(argv[1]);

	printf("Number of nodes: %lu\n",g->n);
	printf("Number of edges: %lu\n",g->e);

	printf("Building adjacency array\n");
	mkadjlist(g);

	printf("Will print resulting hierarchy in file %s\n",argv[2]);
	file=fopen(argv[2],"w");

	printf("Starting recursive bisections\n");
	recurs(bisec, g, 0, file);
	fclose(file);

	t2=time(NULL);

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

	return 0;
}

