#ifndef _RANKSVM_H
#define _RANKSVM_H
#include "linear.h"
#include "tron.h"
#ifdef __cplusplus
extern "C" {
#endif
void eval_list(double *label, double *target, int *query, int l, double *result_ret);
void rank_cross_validation(const problem *prob, const parameter *param, int nr_fold, double *result);

struct id_and_value
{
	int id;
	double value;
};

class l2r_l2_ranksvm_fun: public function
{
public:
	l2r_l2_ranksvm_fun(const problem *prob, double C);
	~l2r_l2_ranksvm_fun();

	double fun(double *w);
	void grad(double *w, double *g);
	void Hv(double *s, double *Hs);

	int get_nr_variable(void);

private:
	void Xv(double *v, double *Xv);
	void XTv(double *v, double *XTv);

	double C;
	double *z;
	int *l_plus;
	int *l_minus;
	double *gamma_plus;
	double *gamma_minus;
	double *ATAXw;
	double *ATe;
	int nr_query;
	int *perm;
	int *start;
	int *count;
	id_and_value **order_perm;
	int *nr_class;
	int *int_y;
	const problem *prob;
};

#ifdef __cplusplus
}
#endif

#endif /* _RANKSVM_H */

