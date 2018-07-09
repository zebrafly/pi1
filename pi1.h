#include "lhep.h"
/*! \file pi1.h
    \brief File describing the pi1 scheme for verifying the computation

	The input x and function F of pi1 is not encrypted,so the client just pass (x,F)
  to the server and the server return an out y to the client.Before passing the
  (x,F) ,the client will generate a verification key to verify if y=Fx after
  the computation.

*/

/*! \struct vc_k
 * F is the function passing in.m and d is the row and column of the function
 * passing in.
 * G and T is the verification key .
 */

typedef struct
{
  int m;
  int d;
  fq_mat_t F;
  g1_t *G;
  g1_t *T;
} vc_k; // VKf

/*! \struct vc_p
 * x is the encrypted message.In pi1, the message is not encrypted,so x=message.
 */

typedef struct
{
  bn_t *x;
} vc_p; // vkx
/**
 * Initializing the pi1 scheme
 */
int pi1_init();
/**
 * stop the program
 */
int pi1_close();
/**
 * generate the verification key
 * @param[in] F 	passing the F to the vck->F
 * @param[in] par some useful parameters
 * @param[out]vck output the function vck->F and the verification key: vck->G ,vck ->T
 */
int pi1_keygen(vc_k *vck, const fq_mat_t F,lhe_par *par);
/**
 * generate a random matrix F stored fq_t
 * @param[in] m   the row of the random matrix
 * @param[in] d   the column of the random matrix
 * @param[out]F   the random matrix
 */
int fq_mat_randz(fq_mat_t F, const int m, const int d, lhe_par *par);
/**
 * compute the one row function matrix multiply the message matrix n given by
 * pi1_comp(). c=LL*C
 * @param[in] C   the message function with size d*1
 * @param[in] LL  one row of the function matrix with size 1*d
 * @param[in] d   the row number of C (the column number of LL)
 * @param[in] par some useful parameters
 * @param[out]c   the outcome of LL*C with size 1*1
 */
int pi1_eval(fq_t c, bn_t *C, const fq_t *LL, const int d, lhe_par *par);
/**
 * compute the function matrix multiply the message marix. nv=vck->F*vcp->x
 * @param[in] vck   the vck->F is needed here to compute
 * @param[in] vcp   the vcp->x is needed here to compute
 * @param[in] par   some useful parameters
 * @param[out]nv    the outcome of vck->F*vcp->x with size m*1
 */
int pi1_comp(fq_t *nv, const vc_k *vck, const vc_p *vcp, lhe_par *par);
/**
 * verify the outcome. output 0 if the outcome is true, otherwise 2
 * @param[in] vck   the vck->G and the vck->T is needed here to verify
 * @param[in] vcp   the input vcp->x is needed here to verify
 * @param[in] nv    the output nv is needed here to verify
 * @param[in] par   some useful parameters
 */
int pi1_vrfy(const vc_k *vck, const vc_p *vcp,fq_t *nv, lhe_par *par);
