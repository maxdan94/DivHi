#include "bisection.h"


//generating n random labels (boolean values)
bool *init(adjlist *g){
	unsigned long i,n=g->n,s=0;
	bool *lab=malloc(n*sizeof(bool));

	//random side for each node:
	for (i=0;i<n;i++){
		lab[i]=rand() > (RAND_MAX / 2);
		if (lab[i]==0)
			s++;
	}

	//avoids all nodes on one side:
	if (s==n)
		lab[0]=1;
	else if (s==0)
		lab[0]=0;

	return lab;
}

//A greedy heuristic for sparsest cut
bool *greedySparsestcut(adjlist *g) {
	unsigned long n=g->n,u,i,j,s=0;
	double cut=0;
	unsigned iter=0;
	long int k,d;
	bool *lab=init(g);
	bool b;

	for (i=0;i<n;i++){
		if (lab[i]==0)
			s++;
	}

	//computing cut
	for (i=0;i<g->e;i++){
		if (lab[g->edges[i].s]!=lab[g->edges[i].t])
			cut++;
	}

	do {//greedy heuristics
		iter++;
		b=0;
		for (u=0;u<n;u++) {
			k=0;
			d=g->cd[u+1]-g->cd[u];
			for (j=g->cd[u];j<g->cd[u+1];j++){
				if (lab[g->adj[j]]==0){
					k++;
				}
			}
			if (lab[u]==0 && s>1 && cut/(s*(n-s)) > (cut+2*k-d)/((s-1)*(n-s+1))){
					s--;
					cut+=2*k-d;
					b=1;
					lab[u]=1;
			}
			else if (lab[u]==1 && n-s>1 && cut/(s*(n-s)) > (cut+d-2*k)/((s+1)*(n-s-1))){
					s++;
					cut+=d-2*k;
					b=1;
					lab[u]=0;
			}
		}
	} while (b);

	//printf("Number of iterations = %u\n",iter);

	return lab;
}

//A greedy heuristic for densest cut
bool *greedyDensestcut(adjlist *g) {
	unsigned long n=g->n,u,i,j,s=0;
	double cut=0;
	unsigned iter=0;
	long int k,d;
	bool *lab=init(g);
	bool b;

	for (i=0;i<n;i++){
		if (lab[i]==0)
			s++;
	}

	//computing cut
	for (i=0;i<g->e;i++){
		if (lab[g->edges[i].s]!=lab[g->edges[i].t])
			cut++;
	}

	do {//greedy heuristics
		iter++;
		b=0;
		for (u=0;u<n;u++) {
			k=0;
			d=g->cd[u+1]-g->cd[u];
			for (j=g->cd[u];j<g->cd[u+1];j++){
				if (lab[g->adj[j]]==0){
					k++;
				}
			}
			if (lab[u]==0 && s>1 && cut/(s*(n-s)) < (cut+2*k-d)/((s-1)*(n-s+1))){
					s--;
					cut+=2*k-d;
					b=1;
					lab[u]=1;
			}
			else if (lab[u]==1 && n-s>1 && cut/(s*(n-s)) < (cut+d-2*k)/((s+1)*(n-s-1))){
					s++;
					cut+=d-2*k;
					b=1;
					lab[u]=0;
			}
		}
	} while (b);

	//printf("Number of iterations = %u\n",iter);

	return lab;
}

//A greedy heuristic for maxcut
bool *greedyMaxcut(adjlist *g) {
	unsigned long n=g->n,u,j;
	unsigned iter=0;
	int k;
	bool *lab=init(g);
	bool b;
	
	do {
		iter++;
		b=0;
		for (u=0;u<n;u++) {
			k=0;
			for (j=g->cd[u];j<g->cd[u+1];j++){
				if (lab[g->adj[j]]){
					k++;
				}
				else{
					k--;
				}
			}
			if (k>0){
				if (lab[u]==1){
					b=1;
					lab[u]=0;
				}
			}
			else if (k<0){
				if (lab[u]==0){
					b=1;
					lab[u]=1;
				}
			}
			//else{//leads to a better cut
			//		lab[u]=rand() > (RAND_MAX / 2);
			//}
		}
	} while (b);

	//printf("Number of iterations = %u\n",iter);

	return lab;
}

bisection choosebisection(char *c){
	printf("Chosen bisection algorithm: ");
	if (strcmp(c,"0")==0){
		printf("Random Bisection\n");
		return init;
	}
	if (strcmp(c,"1")==0){
		printf("Greedy Sparsest Cut\n");
		return greedySparsestcut;
	}
	if (strcmp(c,"2")==0){
		printf("Greedy Densest Cut\n");
		return greedyDensestcut;
	}
	if (strcmp(c,"3")==0){
		printf("Greedy Maxcut\n");
		return greedyMaxcut;
	}
	printf("unknown\n");
	exit(1);
}


