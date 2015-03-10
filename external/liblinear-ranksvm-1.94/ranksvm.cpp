#include "ranksvm.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
template <class T> static inline void swap(T& x, T& y) { T t=x; x=y; y=t; }
struct tree_node
{
	int size;
	double xv;
};

class selectiontree
{
public:
	selectiontree(int k); //k leaves
	~selectiontree();
	void insert_node(int key, double value);
	void larger(int key, int *l_plus_ret, double *gamma_plus_ret);
	void smaller(int key, int *l_minus_ret, double *gamma_minus_ret);
	double xv_larger(int key);
	double xv_smaller(int key);

private:
	int num_nonleaves;
	int num_leaves;
	tree_node *node;
};

selectiontree::selectiontree(int k)
{
	int i = 1;
	int j;
	while(i < k)
		i *= 2;
	this->num_leaves = k;
	this->num_nonleaves = i-1;
	node = new tree_node[i+k];
	for (j=0;j<i+k;j++)
	{
		node[j].size = 0;
		node[j].xv = 0;
	}
}

selectiontree::~selectiontree()
{
	delete[] node;
}

void selectiontree::insert_node(int key, double value)
{
	key += this->num_nonleaves;
	for (;key>0;key/=2)
	{
		node[key].xv += value;
		node[key].size++;
	}
}

void selectiontree::larger(int key, int *l_plus_ret, double *gamma_plus_ret)
{
	if (key >= this->num_leaves)
	{
		*l_plus_ret = 0;
		*gamma_plus_ret = 0;
		return;
	}
	int l_plus = 0;
	double gamma_plus = 0;
	key += num_nonleaves;
	for (;key>1;key/=2)
		if (key % 2 == 0)
		{
			l_plus += node[key+1].size;
			gamma_plus += node[key+1].xv;
		}
	*l_plus_ret = l_plus;
	*gamma_plus_ret = gamma_plus;
}

void selectiontree::smaller(int key, int *l_minus_ret, double *gamma_minus_ret)
{
	if (key <= 1)
	{
		*l_minus_ret = 0;
		*gamma_minus_ret = 0;
		return;
	}
	int l_minus = 0;
	double gamma_minus = 0;
	key += num_nonleaves;
	for (;key>1;key/=2)
		if (key % 2 == 1)
		{
			l_minus += node[key-1].size;
			gamma_minus += node[key-1].xv;
		}
	*l_minus_ret = l_minus;
	*gamma_minus_ret = gamma_minus;
}

double selectiontree::xv_smaller(int key)
{
	if (key <= 1)
		return 0;
	double gamma_minus = 0;
	key += num_nonleaves;
	for (;key>1;key/=2)
		if (key % 2 == 1)
			gamma_minus += node[key-1].xv;
	return gamma_minus;
}

double selectiontree::xv_larger(int key)
{
	if (key >= this->num_leaves)
		return 0;
	double gamma_plus = 0;
	key += num_nonleaves;
	for (;key>1;key/=2)
		if (key % 2 == 0)
			gamma_plus += node[key+1].xv;
	return gamma_plus;
}

static int compare_values(const void *a, const void *b)
{
	struct id_and_value *ia = (struct id_and_value *)a;
	struct id_and_value *ib = (struct id_and_value *)b;
	if(ia->value > ib->value)
		return 1;
	if(ia->value < ib->value)
		return -1;
	return 0;
}

// start: begin of each query, count: #data of queries, perm: indices to the original data
// perm, length l, must be allocated before calling this subroutine
static void group_queries(const int *query_id, int l, int *nr_query_ret, int **start_ret, int **count_ret, int *perm)
{
	int max_nr_query = 16;
	int nr_query = 0;
	int *query = Malloc(int,max_nr_query);
	int *count = Malloc(int,max_nr_query);
	int *data_query = Malloc(int,l);
	int i;

	for(i=0;i<l;i++)
	{
		int this_query = (int)query_id[i];
		int j;
		for(j=0;j<nr_query;j++)
		{
			if(this_query == query[j])
			{
				++count[j];
				break;
			}
		}
		data_query[i] = j;
		if(j == nr_query)
		{
			if(nr_query == max_nr_query)
			{
				max_nr_query *= 2;
				query = (int *)realloc(query,max_nr_query * sizeof(int));
				count = (int *)realloc(count,max_nr_query * sizeof(int));
			}
			query[nr_query] = this_query;
			count[nr_query] = 1;
			++nr_query;
		}
	}

	int *start = Malloc(int,nr_query);
	start[0] = 0;
	for(i=1;i<nr_query;i++)
		start[i] = start[i-1] + count[i-1];
	for(i=0;i<l;i++)
	{
		perm[start[data_query[i]]] = i;
		++start[data_query[i]];
	}
	start[0] = 0;
	for(i=1;i<nr_query;i++)
		start[i] = start[i-1] + count[i-1];

	*nr_query_ret = nr_query;
	*start_ret = start;
	*count_ret = count;
	free(query);
	free(data_query);
}

l2r_l2_ranksvm_fun::l2r_l2_ranksvm_fun(const problem *prob, double C)
{
	int q,i,k;
	int l=prob->l;
	double *y = prob->y;

	this->prob = prob;

	z = new double[l];
	l_plus = new int[l];
	l_minus = new int[l];
	gamma_plus = new double[l];
	gamma_minus = new double[l];
	ATAXw = new double[l];
	ATe = new double[l];
	int_y = new int[l];
	this->C = C;
// Checking number of relevance levels in each query
// Meanwhile transform the labels into 1,...,k
	perm = Malloc(int,l);
	group_queries(prob->query, l, &nr_query, &start, &count, perm);
	nr_class = new int[nr_query];
	order_perm = new id_and_value *[nr_query];
	for (q=0;q<nr_query;q++)
		order_perm[q] = new id_and_value[count[q]];
	for (q=0;q<nr_query;q++)
	{
		int *perm_q = &perm[start[q]];
		id_and_value *order_perm_q = order_perm[q];
		k = 1;
		for (i=0;i<count[q];i++)
		{
			order_perm_q[i].id = perm_q[i];
			order_perm_q[i].value = y[perm_q[i]];
		}
		qsort(order_perm_q, count[q], sizeof(id_and_value), compare_values);

		int_y[order_perm_q[0].id] = 1;
		for(i=1;i<count[q];i++)
		{
			if (order_perm_q[i-1].value<order_perm_q[i].value)
				k++;
			int_y[order_perm_q[i].id] = k;
		}
		nr_class[q] = k;
	}
}

l2r_l2_ranksvm_fun::~l2r_l2_ranksvm_fun()
{
	delete[] z;
	delete[] l_plus;
	delete[] l_minus;
	delete[] gamma_plus;
	delete[] gamma_minus;
	delete[] ATAXw;
	delete[] ATe;
	delete[] int_y;
	delete[] nr_class;
	delete[] start;
	delete[] count;
	delete[] perm;
	for (int q=0;q<nr_query;q++)
		delete[] order_perm[q];
	delete[] order_perm;
}

double l2r_l2_ranksvm_fun::fun(double *w)
{
	int q,i,j;
	double f = 0;
	int l = prob->l;
	int w_size = get_nr_variable();
	selectiontree *T;
	Xv(w,z);
	for (q=0;q<nr_query;q++)
	{
		int *perm_q = &perm[start[q]];
		id_and_value *order_perm_q = order_perm[q];
		for (i=0;i<count[q];i++)
		{
			order_perm_q[i].id = perm_q[i];
			order_perm_q[i].value = z[perm_q[i]];
		}
		qsort(order_perm_q, count[q], sizeof(id_and_value), compare_values);

		T = new selectiontree(nr_class[q]);
		j = 0;
		for (i=0;i<count[q];i++)
		{
			while (j<count[q] && (1 - order_perm_q[j].value + order_perm_q[i].value>0))
			{
				T->insert_node(int_y[order_perm_q[j].id], order_perm_q[j].value);
				j++;
			}
			T->larger(int_y[order_perm_q[i].id], &l_plus[order_perm_q[i].id], &gamma_plus[order_perm_q[i].id]);
		}
		delete T;

		j = count[q] - 1;
		T = new selectiontree(nr_class[q]);
		for (i=count[q]-1;i>=0;i--)
		{
			while (j>=0 && (1 - order_perm_q[i].value + order_perm_q[j].value>0))
			{
				T->insert_node(int_y[order_perm_q[j].id], order_perm_q[j].value);
				j--;
			}
			T->smaller(int_y[order_perm_q[i].id], &l_minus[order_perm_q[i].id], &gamma_minus[order_perm_q[i].id]);
		}
		delete T;
	}
	for(i=0;i<w_size;i++)
		f += w[i] * w[i];
	f /= 2.0;

	for (i=0;i<l;i++)
	{
		ATAXw[i] = (l_plus[i]+l_minus[i])*z[i] - gamma_plus[i] - gamma_minus[i];
		ATe[i] = l_minus[i] - l_plus[i];
	}
	for(i=0;i<l;i++)
		f += C * (z[i] * (ATAXw[i] - 2 * ATe[i]) + l_minus[i]);
	return(f);
}

void l2r_l2_ranksvm_fun::grad(double *w, double *g)
{
	int i;
	int l = prob->l;
	double *tmp_vector;
	tmp_vector = new double[l];
	int w_size = get_nr_variable();
	for (i=0;i<l;i++)
		tmp_vector[i] = ATAXw[i] - ATe[i];
	XTv(tmp_vector, g);
	for(i=0;i<w_size;i++)
		g[i] = w[i] + 2 * C * g[i];
	delete[] tmp_vector;
}

int l2r_l2_ranksvm_fun::get_nr_variable(void)
{
	return prob->n;
}

void l2r_l2_ranksvm_fun::Hv(double *s, double *Hs)
{
	int q,i,j;
	int w_size = get_nr_variable();
	int l = prob->l;
	double *wa = new double[l];
	selectiontree *T;
	double* gamma_plus_minus;
	gamma_plus_minus = new double[l];
	Xv(s, wa);
	for (q=0;q<nr_query;q++)
	{
		id_and_value *order_perm_q = order_perm[q];
		T = new selectiontree(nr_class[q]);
		j = 0;
		for (i=0;i<count[q];i++)
		{
			while (j<count[q] && (1 - order_perm_q[j].value + order_perm_q[i].value>0))
			{
				T->insert_node(int_y[order_perm_q[j].id],wa[order_perm_q[j].id]);
				j++;
			}
			gamma_plus_minus[order_perm_q[i].id] = T->xv_larger(int_y[order_perm_q[i].id]);
		}
		delete T;

		j = count[q] - 1;
		T = new selectiontree(nr_class[q]);
		for (i=count[q]-1;i>=0;i--)
		{
			while (j>=0 && (1 - order_perm_q[i].value + order_perm_q[j].value>0))
			{
				T->insert_node(int_y[order_perm_q[j].id], wa[order_perm_q[j].id]);
				j--;
			}
			gamma_plus_minus[order_perm_q[i].id] += T->xv_smaller(int_y[order_perm_q[i].id]);
		}
		delete T;
	}
	for (i=0;i<l;i++)
		wa[i] = wa[i] * (l_plus[i] + l_minus[i]) - gamma_plus_minus[i];
	delete[] gamma_plus_minus;
	XTv(wa, Hs);
	delete[] wa;
	for(i=0;i<w_size;i++)
		Hs[i] = s[i] + 2 * C * Hs[i];
}

void l2r_l2_ranksvm_fun::Xv(double *v, double *Xv)
{
	int i;
	int l=prob->l;
	feature_node **x=prob->x;

	for(i=0;i<l;i++)
	{
		feature_node *s=x[i];
		Xv[i]=0;
		while(s->index!=-1)
		{
			Xv[i]+=v[s->index-1]*s->value;
			s++;
		}
	}
}

void l2r_l2_ranksvm_fun::XTv(double *v, double *XTv)
{
	int i;
	int l=prob->l;
	int w_size=get_nr_variable();
	feature_node **x=prob->x;

	for(i=0;i<w_size;i++)
		XTv[i]=0;
	for(i=0;i<l;i++)
	{
		feature_node *s=x[i];
		while(s->index!=-1)
		{
			XTv[s->index-1]+=v[i]*s->value;
			s++;
		}
	}
}

void eval_list(double *label, double *target, int *query, int l, double *result_ret)
{
	int q,i,j,k;
	int nr_query;
	int *start = NULL;
	int *count = NULL;
	int *perm = Malloc(int, l);
	id_and_value *order_perm;
	int true_query;
	long long totalnc = 0, totalnd = 0;
	long long nc = 0;
	long long nd = 0;
	double tmp;
	double accuracy = 0;
	int *l_plus;
	int *int_y;
	int same_y = 0;
	double *ideal_dcg;
	double *dcg;
	double meanndcg = 0;
	double ndcg;
	selectiontree *T;
	group_queries(query, l, &nr_query, &start, &count, perm);
	true_query = nr_query;
	for (q=0;q<nr_query;q++)
	{
		//We use selection trees to compute pairwise accuracy
		nc = 0;
		nd = 0;
		l_plus = new int[count[q]];
		int_y = new int[count[q]];
		order_perm = new id_and_value[count[q]];
		int *perm_q = &perm[start[q]];
		for (i=0;i<count[q];i++)
		{
			order_perm[i].id = i;
			order_perm[i].value = label[perm_q[i]];
		}
		qsort(order_perm, count[q], sizeof(id_and_value), compare_values);
		int_y[order_perm[0].id] = 1;
		same_y = 0;
		k = 1;
		for(i=1;i<count[q];i++)
		{
			if (order_perm[i-1].value < order_perm[i].value)
			{
				same_y = 0;
				k++;
			}
			else
				same_y++;
			int_y[order_perm[i].id] = k;
			nc += (i - same_y);
		}
		for (i=0;i<count[q];i++)
		{
			order_perm[i].id = i;
			order_perm[i].value = target[perm_q[i]];
		}
		qsort(order_perm, count[q], sizeof(id_and_value), compare_values);
		//total pairs
		T = new selectiontree(k);
		j = 0;
		for (i=0;i<count[q];i++)
		{
			while (j<count[q] && ( order_perm[j].value < order_perm[i].value))
			{
				T->insert_node(int_y[order_perm[j].id], tmp);
				j++;
			}
			T->larger(int_y[order_perm[i].id], &l_plus[order_perm[i].id], &tmp);
		}
		delete T;

		for (i=0;i<count[q];i++)
			nd += l_plus[i];
		nc -= nd;
		if (nc != 0 || nd != 0)
			accuracy += double(nc)/double(nc+nd);
		else
			true_query--;
		totalnc += nc;
		totalnd += nd;
		delete[] l_plus;
		delete[] int_y;
		delete[] order_perm;
	}
	result_ret[0] = (double)totalnc/(double)(totalnc+totalnd);
	for (q=0;q<nr_query;q++)
	{
		//mean ndcg by the formulation of LETOR
		ideal_dcg = new double[count[q]];
		dcg = new double[count[q]];
		ndcg = 0;
		order_perm = new id_and_value[count[q]];
		int *perm_q = &perm[start[q]];
		for (i=0;i<count[q];i++)
		{
			order_perm[i].id = perm_q[i];
			order_perm[i].value = label[perm_q[i]];
		}
		qsort(order_perm, count[q], sizeof(id_and_value), compare_values);
		ideal_dcg[0] = pow(2.0,order_perm[count[q]-1].value) - 1;
		for (i=1;i<count[q];i++)
			ideal_dcg[i] = ideal_dcg[i-1] + (pow(2.0,order_perm[count[q]-1 - i].value) - 1) * log(2.0) / log(i+1.0);
		for (i=0;i<count[q];i++)
		{
			order_perm[i].id = perm_q[i];
			order_perm[i].value = target[perm_q[i]];
		}
		qsort(order_perm, count[q], sizeof(id_and_value), compare_values);
		dcg[0] = pow(2.0, label[order_perm[count[q] - 1].id]) - 1;
		for (i=1;i<count[q];i++)
			dcg[i] = dcg[i-1] + (pow(2.0, label[order_perm[count[q] - 1 - i].id]) - 1) * log(2.0) / log(i + 1.0);
		if (ideal_dcg[0]>0)
			for (i=0;i<count[q];i++)
				ndcg += dcg[i]/ideal_dcg[i];
		else
			ndcg = 0;
		meanndcg += ndcg/count[q];
		delete[] order_perm;
		delete[] ideal_dcg;
		delete[] dcg;
	}
	meanndcg /= nr_query;
	result_ret[1] = meanndcg;
	free(start);
	free(count);
	free(perm);
}

void rank_cross_validation(const problem *prob, const parameter *param, int nr_fold, double *result)
{
	int i,q;
	int *fold_start;
	int l = prob->l;
	int *query_set;
	double *target = Malloc(double,l);
	int nr_query;
	int *start = NULL;
	int *count = NULL;
	int *perm = Malloc(int,l);
	int *query_perm;
	group_queries(prob->query, l, &nr_query, &start, &count, perm);
	if (nr_query == 1)
	{
		if (nr_fold > prob->l / 2)
		{
			nr_fold = l / 2; // each fold should include at least 2 instances to form pairs
			fprintf(stderr,"WARNING: # folds > # data / 2. Will use # folds = # data / 2 instead (Every fold should contain 2 data to form a pair)\n");
		}
		nr_query = nr_fold;
// Treat each fold as a query in performance evaluation
// to avoid ranking inconsistency between models.
		start = (int *)realloc(start,nr_query * sizeof(int));
		count = (int *)realloc(count,nr_query * sizeof(int));
		query_set = Malloc(int,l);
		for(q=0;q<nr_query;q++)
		{
			count[q] = 0;
			start[q] = 0;
		}
		for(i=0;i<l;i++)
		{
			int j = rand() % nr_query;
			query_set[i] = j;
			count[j]++;
		}
		start[0] = 0;
		for(q=1;q<nr_query;q++)
			start[q] = start[q-1] + count[q-1];
		for(i=0;i<l;i++)
		{
			perm[start[query_set[i]]] = i;
			++start[query_set[i]];
		}
		start[0] = 0;
		for(q=1;q<nr_query;q++)
			start[q] = start[q-1] + count[q-1];
	}
	else
	{
		query_set = prob->query;
		if (nr_query < nr_fold)
		{
			nr_fold = nr_query;
			fprintf(stderr,"WARNING: # folds > # query. Will use # folds = # query instead.\n");
		}
	}
	fold_start = Malloc(int,nr_fold+1);
	query_perm = Malloc(int,nr_query);
	for(i=0;i<=nr_fold;i++)
		fold_start[i] = i * nr_query / nr_fold;
	for (q=0;q<nr_query;q++)
		query_perm[q] = q;
	for (q=0;q<nr_query;q++)
	{
		i = q + rand() % (nr_query-q);
		swap(query_perm[q], query_perm[i]);
	}

	for(i=0;i<nr_fold;i++)
	{
		int begin = fold_start[i];
		int end = fold_start[i+1];
		int j,m,counter;
		struct problem subprob;
		counter = 0;
		subprob.n = prob->n;

		for (q=begin;q<end;q++)
			counter += count[query_perm[q]];
		subprob.l = l - counter;
		subprob.x = Malloc(struct feature_node*,subprob.l);
		subprob.y = Malloc(double,subprob.l);
		subprob.query = Malloc(int,subprob.l);

		m = 0;
		for(q=0;q<begin;q++)
		{
			int *perm_q = &perm[start[query_perm[q]]];
			for (j=0;j<count[query_perm[q]];j++)
			{
				subprob.x[m] = prob->x[perm_q[j]];
				subprob.y[m] = prob->y[perm_q[j]];
				subprob.query[m] = prob->query[perm_q[j]];
				++m;
			}
		}
		for(q=end;q<nr_query;q++)
		{
			int *perm_q = &perm[start[query_perm[q]]];
			for (j=0;j<count[query_perm[q]];j++)
			{
				subprob.x[m] = prob->x[perm_q[j]];
				subprob.y[m] = prob->y[perm_q[j]];
				subprob.query[m] = prob->query[perm_q[j]];
				++m;
			}
		}
		struct model *submodel = train(&subprob,param);

		for(q=begin;q<end;q++)
		{
			int *perm_q = &perm[start[query_perm[q]]];
			for (j=0;j<count[query_perm[q]];j++)
				target[perm_q[j]] = predict(submodel,prob->x[perm_q[j]]);
		}

		free_and_destroy_model(&submodel);
		free(subprob.x);
		free(subprob.y);
		free(subprob.query);
	}
	eval_list(prob->y,target,query_set,l,result);
	free(fold_start);
	free(count);
	free(start);
	free(perm);
	free(target);
	free(query_perm);
	if (nr_query == 1)
		free(query_set);
}

