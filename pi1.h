#include "lhep.h"


typedef struct
{
  int m;
  int d;
  fq_mat_t F;
  g1_t *G;
  g1_t *T;
} vc_k; // VKf

typedef struct
{
  bn_t *x;
} vc_p; // vkx

int pi1_init();

int pi1_close();

int pi1_keygen(vc_k *vck, fq_mat_t F, lhe_par *par);

int fq_mat_randz(fq_mat_t F, int m, int d, lhe_par *par);

int pi1_eval(fq_t c, bn_t *C, fq_t *LL, int d, lhe_par *par);

int pi1_comp(fq_t *nv, vc_k *vck, vc_p *vcp, lhe_par *par);

int pi1_vrfy(vc_k *vck, vc_p *vcp, fq_t *nv, lhe_par *par);
