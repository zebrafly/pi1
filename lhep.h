//#ifndef _BN_ETX_H_
//#define _BN_ETX_H_
#include "bn_ext.h"
//#endif

//the public parameter for the lhe linearly homomorphic encryption scheme
typedef struct
{
	bn_t p; //relic
	fmpz_t pf; //flint
	bn_t q;
	fmpz_t qf;
	bn_t n;
	int ni;
	fmpz_t nf;
	bn_t r; //正太分布
	fmpz_t rf;// 转换类型
	fq_ctx_t ctx;// context manager
	fq_poly_t modf;//Rp Rq 除的fx
	fq_poly_t modp;
	g1_t g; // 指数
	g2_t h; // 指数
	gt_t gt;// e(x^g,y^h)=gt
} lhe_par;

//generate the parameters for the encryption scheme
int lhep_new(lhe_par *par);

//free the structure
int lhep_free(lhe_par *par);
