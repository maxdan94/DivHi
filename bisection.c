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
bool *greedySparsestcut(adjlist *g, double* val) {
	unsigned long n=g->n,u,v,i,s=0;
	double cut=0;
	unsigned iter=0;
	bool *lab=init(g);
	bool b;
	int *d0=calloc(g->n,sizeof(int)),*d1=calloc(g->n,sizeof(int));

	for (u=0;u<n;u++){
		if (lab[u]==0){
			s++;
		}
	}

	for (i=0;i<g->e;i++){
		u=g->edges[i].s;
		v=g->edges[i].t;
		if (lab[u]!=lab[v]){
			cut++;
		}
		if (lab[v]==0){
			d0[u]++;
		}
		else{
			d1[u]++;
		}
		if (lab[u]==0){
			d0[v]++;
		}
		else{
			d1[v]++;
		}
	}

	do {//greedy heuristics
		iter++;
		b=0;
		for (u=0;u<n;u++) {
			if (lab[u]==0 && s>1 && cut/(s*(n-s)) > (cut+d0[u]-d1[u])/((s-1)*(n-s+1))){
				b=1;
				s--;
				lab[u]=1;
				cut+=d0[u]-d1[u];
				for (i=g->cd[u];i<g->cd[u+1];i++){
					v=g->adj[i];
					d0[v]--;
					d1[v]++;
				}
			}
			else if (lab[u]==1 && n-s>1 && cut/(s*(n-s)) > (cut+d1[u]-d0[u])/((s+1)*(n-s-1))){
				b=1;
				s++;
				lab[u]=0;
				cut+=d1[u]-d0[u];
				for (i=g->cd[u];i<g->cd[u+1];i++){
					v=g->adj[i];
					d0[v]++;
					d1[v]--;
				}
			}
		}
	} while (b);

	//printf("Number of iterations = %u\n",iter);
	free(d0);
	free(d1);
	*val=cut/(s*(n-s));
	return lab;
}

//A greedy heuristic for densest cut
bool *greedyDensestcut(adjlist *g,double* val) {
	unsigned long n=g->n,u,v,i,s=0;
	double cut=0;
	unsigned iter=0;
	bool *lab=init(g);
	bool b;
	int *d0=calloc(g->n,sizeof(int)),*d1=calloc(g->n,sizeof(int));

	for (u=0;u<n;u++){
		if (lab[u]==0){
			s++;
		}
	}

	for (i=0;i<g->e;i++){
		u=g->edges[i].s;
		v=g->edges[i].t;
		if (lab[u]!=lab[v]){
			cut++;
		}
		if (lab[v]==0){
			d0[u]++;
		}
		else{
			d1[u]++;
		}
		if (lab[u]==0){
			d0[v]++;
		}
		else{
			d1[v]++;
		}
	}

	do {//greedy heuristics
		iter++;
		b=0;
		for (u=0;u<n;u++) {
			if (lab[u]==0 && s>1 && cut/(s*(n-s)) < (cut+d0[u]-d1[u])/((s-1)*(n-s+1))){
				b=1;
				s--;
				lab[u]=1;
				cut+=d0[u]-d1[u];
				for (i=g->cd[u];i<g->cd[u+1];i++){
					v=g->adj[i];
					d0[v]--;
					d1[v]++;
				}
			}
			else if (lab[u]==1 && n-s>1 && cut/(s*(n-s)) < (cut+d1[u]-d0[u])/((s+1)*(n-s-1))){
				b=1;
				s++;
				lab[u]=0;
				cut+=d1[u]-d0[u];
				for (i=g->cd[u];i<g->cd[u+1];i++){
					v=g->adj[i];
					d0[v]++;
					d1[v]--;
				}
			}
		}
	} while (b);

	//printf("Number of iterations = %u\n",iter);

	free(d0);
	free(d1);
	*val=-cut/(s*(n-s));
	return lab;
}



//A greedy heuristic for maxcut
bool *greedyMaxcut(adjlist *g,double *val) {
	unsigned long n=g->n,u,v,i;
	double cut=0;
	unsigned iter=0;
	bool *lab=init(g);
	bool b;
	int *d0=calloc(g->n,sizeof(int)),*d1=calloc(g->n,sizeof(int));

	for (i=0;i<g->e;i++){
		u=g->edges[i].s;
		v=g->edges[i].t;
		if (lab[u]!=lab[v]){
			cut++;
		}
		if (lab[v]==0){
			d0[u]++;
		}
		else{
			d1[u]++;
		}
		if (lab[u]==0){
			d0[v]++;
		}
		else{
			d1[v]++;
		}
	}

	do {//greedy heuristics
		iter++;
		b=0;
		for (u=0;u<n;u++) {
			if (lab[u]==0 && d1[u] < d0[u] ){
				b=1;
				lab[u]=1;
				cut+=d0[u]-d1[u];
				for (i=g->cd[u];i<g->cd[u+1];i++){
					v=g->adj[i];
					d0[v]--;
					d1[v]++;
				}
			}
			else if (lab[u]==1 && d0[u] < d1[u] ){
				b=1;
				lab[u]=0;
				cut+=d1[u]-d0[u];
				for (i=g->cd[u];i<g->cd[u+1];i++){
					v=g->adj[i];
					d0[v]++;
					d1[v]--;
				}
			}
		}
	} while (b);

	//printf("Number of iterations = %u\n",iter);

	free(d0);
	free(d1);
	*val=cut;
	return lab;
}



//A greedy heuristic for min cut of the complement graph
bool *greedyMincut(adjlist *g,double *val) {
	unsigned long n=g->n,u,v,i,s=0;
	double cut=0;
	unsigned iter=0;
	bool *lab=init(g);
	bool b;
	int *d0=calloc(g->n,sizeof(int)),*d1=calloc(g->n,sizeof(int));

	for (u=0;u<n;u++){
		if (lab[u]==0){
			s++;
		}
	}

	for (i=0;i<g->e;i++){
		u=g->edges[i].s;
		v=g->edges[i].t;
		if (lab[u]!=lab[v]){
			cut++;
		}
		if (lab[v]==0){
			d0[u]++;
		}
		else{
			d1[u]++;
		}
		if (lab[u]==0){
			d0[v]++;
		}
		else{
			d1[v]++;
		}
	}

	do {//greedy heuristics
		iter++;
		b=0;
		for (u=0;u<n;u++) {
			if (lab[u]==0 && s>1 && cut > (cut+d0[u]-d1[u])){
				b=1;
				s--;
				lab[u]=1;
				cut+=d0[u]-d1[u];
				for (i=g->cd[u];i<g->cd[u+1];i++){
					v=g->adj[i];
					d0[v]--;
					d1[v]++;
				}
			}
			else if (lab[u]==1 && n-s>1 && cut > (cut+d1[u]-d0[u])){
				b=1;
				s++;
				lab[u]=0;
				cut+=d1[u]-d0[u];
				for (i=g->cd[u];i<g->cd[u+1];i++){
					v=g->adj[i];
					d0[v]++;
					d1[v]--;
				}
			}
		}
	} while (b);

	//printf("Number of iterations = %u\n",iter);
	//printf("cut, s, n-s = %lf %lu %lu\n",cut,s,n-s);
	free(d0);
	free(d1);
	*val=cut;
	return lab;
}



//A greedy heuristic for min cut of the complement graph
bool *greedyCompmincut(adjlist *g,double* val) {
	unsigned long n=g->n,u,v,i,s=0;
	double cut=0;
	unsigned iter=0;
	bool *lab=init(g);
	bool b;
	int *d0=calloc(g->n,sizeof(int)),*d1=calloc(g->n,sizeof(int));

	for (u=0;u<n;u++){
		if (lab[u]==0){
			s++;
		}
	}

	for (i=0;i<g->e;i++){
		u=g->edges[i].s;
		v=g->edges[i].t;
		if (lab[u]!=lab[v]){
			cut++;
		}
		if (lab[v]==0){
			d0[u]++;
		}
		else{
			d1[u]++;
		}
		if (lab[u]==0){
			d0[v]++;
		}
		else{
			d1[v]++;
		}
	}

	do {//greedy heuristics
		iter++;
		b=0;
		for (u=0;u<n;u++) {
			if (lab[u]==0 && s>1 && (s*(n-s)-cut) > ((s-1)*(n-s+1)-(cut+d0[u]-d1[u]))){
				b=1;
				s--;
				lab[u]=1;
				cut+=d0[u]-d1[u];
				for (i=g->cd[u];i<g->cd[u+1];i++){
					v=g->adj[i];
					d0[v]--;
					d1[v]++;
				}
			}
			else if (lab[u]==1 && n-s>1 && cut-(s*(n-s)) > ((s+1)*(n-s-1)-(cut+d1[u]-d0[u]))){
				b=1;
				s++;
				lab[u]=0;
				cut+=d1[u]-d0[u];
				for (i=g->cd[u];i<g->cd[u+1];i++){
					v=g->adj[i];
					d0[v]++;
					d1[v]--;
				}
			}
		}
	} while (b);

	//printf("Number of iterations = %u\n",iter);
	//printf("cut, s, n-s, s*(n-s)-cut = %lf %lu %lu %lf\n",cut,s,n-s,s*(n-s)-cut);
	free(d0);
	free(d1);
	*val=cut-(s*(n-s));
	return lab;
}



bool *greedyMinmin(adjlist *g,double* val) {
	static double totcut1=0,totcut2=0;
	unsigned long n=g->n,u,v,i,s1=0,s2=0;
	double cut1=0,cut2=0;
	double val1,val2;
	bool *lab1=greedyMincut(g,&val1);
	bool *lab2=greedyCompmincut(g,&val2);
	static FILE *f=NULL;
	if (f==NULL)
		f=fopen("graph.txt","w");

	for (u=0;u<n;u++){
		if (lab1[u]==0){
			s1++;
		}
		if (lab2[u]==0){
			s2++;
		}
	}
	for (i=0;i<g->e;i++){
		u=g->edges[i].s;
		v=g->edges[i].t;
		if (lab1[u]!=lab1[v]){
			cut1++;
		}
		if (lab2[u]!=lab2[v]){
			cut2++;
		}
	}
	if (cut1<s2*(n-s2)-cut2){
		//printf("s, n-s, cut = %lu %lu %lf\n",s1,n-s1,cut1);
		//totcut1+=cut1;
		//printf("number of deletions = %lf\n",totcut1);
		//printf("number of additions = %lf\n",totcut2);
		//printf("number of editions = %lf\n",totcut1+totcut2);
/*////////////
	for (i=0;i<g->e;i++){
		u=g->edges[i].s;
		v=g->edges[i].t;
		if (lab1[u]!=lab1[v]){
			fprintf(f,"%lu %lu\n",g->map[u],g->map[v]);
		}
	}
*////////////
		free(lab2);
		*val=cut1;
		return lab1;
	}

	//printf("s, n-s, s*(n-s)-cut  = %lu %lu %lf\n",s2,n-s2,s2*(n-s2)-cut2 );
	totcut2+=s2*(n-s2)-cut2;
	//printf("number of deletions = %lf\n",totcut1);
	//printf("number of additions = %lf\n",totcut2);
	//printf("number of editions = %lf\n",totcut1+totcut2);
	free(lab1);
	*val=s2*(n-s2)-cut2;
	return lab2;
}


bool *greedyMinmin2(adjlist *g,double* val) {
	static double totcut1=0;
	unsigned long n=g->n,u,v,i,s=0;
	double cut1=0,cut2=0;
	double val1,val2;
	bool *lab1=greedyMincut(g,&val1);
	bool *lab2=greedyCompmincut(g,&val2);
	static FILE *f=NULL;
	if (f==NULL)
		f=fopen("graph.txt","w");

	for (u=0;u<n;u++){
		if (lab2[u]==0){
			s++;
		}
	}
	for (i=0;i<g->e;i++){
		u=g->edges[i].s;
		v=g->edges[i].t;
		if (lab1[u]!=lab1[v]){
			cut1++;
		}
		if (lab2[u]==lab2[v]){
			cut2++;
		}
	}
	if (cut2<1){
		*val=cut2;
		return lab2;
	}


	printf("s, n-s, cut = %lu %lu %lf\n",s,n-s,cut1);
	totcut1+=cut1;
	printf("number of deletions = %lf\n",totcut1);
/////////////
	for (i=0;i<g->e;i++){
		u=g->edges[i].s;
		v=g->edges[i].t;
		if (lab1[u]!=lab1[v]){
			fprintf(f,"%lu %lu\n",g->map[u],g->map[v]);
		}
	}
////////////
	*val=cut1;
	return lab1;

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
	if (strcmp(c,"4")==0){
		printf("Greedy Mincut\n");
		return greedyMincut;
	}
	if (strcmp(c,"5")==0){
		printf("Greedy CompMincut\n");
		return greedyCompmincut;
	}
	if (strcmp(c,"6")==0){
		printf("Greedy minimum between Mincut and CompMincut\n");
		return greedyMinmin;
	}
	if (strcmp(c,"7")==0){
		printf("Greedy minimum between Mincut and CompMincut (CompMincut=0 only)\n");
		return greedyMinmin2;
	}
	printf("unknown\n");
	exit(1);
}


