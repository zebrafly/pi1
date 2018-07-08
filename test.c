#include "pi1.h"
#include <stdio.h>


int main()
{

//initialize the project /////////////////////////////////////////////////
        if (pi1_init())
        {
                printf("Testing FAILED\n");
                printf("Problem initializing the library\n");
               return 1;
        }
        lhe_par *par=malloc(sizeof(*par));
        lhep_new(par);

printf("\n--------------------------------------begin test-------------------------\n\n\n");

        //test the rand matrix algorithm
        int mat_m=100;
        int mat_d=100;
        fq_mat_t F;
        fq_mat_init(F,mat_m,mat_d,par->ctx);
        fq_mat_randz(F,mat_m,mat_d,par);

        //print
        printf("\n\n the matrix F is:\n\n");
        fq_mat_print_pretty(F,par->ctx);

  //
        //test the pi1_keygen algorithm
        vc_k *vck=malloc(sizeof(*vck));
        vck->m=mat_m;
        vck->d=mat_d;

        pi1_keygen(vck,F,par);
        //print
        printf("\n\n the integer vck->m=%d\n\n",vck->m);
        printf("\n\n the integer vck->d=%d\n\n",vck->d);
        printf("\n\n the matrix vck->F is: \n\n");
        fq_mat_print_pretty(vck->F,par->ctx);
        printf("\n\n the vck->G vector is:\n\n");
        for (int i=0;i<vck->m;i++)
        {
         g1_print(vck->G[i]);
         printf("\n\n");
        }
        printf("\n\n the vck->T vector is:\n\n");
        for (int j=0;j<vck->d;j++)
        {
         g1_print(vck->T[j]);
         printf("\n\n");
	      }

        //test the vc_pgen algorithm
        vc_p *vcp=malloc(sizeof(*vcp));
        //generate a random message vector of d elements
        bn_t x[vck->d];
        for (int j=0;j<vck->d;j++)
        {
          bn_rand_mod(x[j],par->p);
        }
        //compute the problem instance
        vcp->x = x;
        printf("\n\n the input x is:\n\n");
        for (int j=0;j<vck->d;j++)
        {
           printf("\n\n this is x[%d]:\n\n",j);
           bn_print(x[j]);
           printf("\n\n");
        }
        //test the pi1_comp algorithm
        fq_t *nv=malloc(sizeof(fq_t)*(vck->m));
        for (int j=0;j<vck->m;j++)
        {
          fq_init(nv[j],par->ctx);
        }
        pi1_comp(nv,vck,vcp,par);

        //print

        printf("\n\n the computation result is:\n\n");
        for (int i=0;i<vck->m;i++)
        {
           printf("\n\n this is nv[%d]:\n\n",i);
           fq_print(nv[i],par->ctx);
           printf("\n\n");
        }

       //test the pi1_vrfy algorithm
       int flag=pi1_vrfy(vck,vcp,nv,par);
       printf("\n\n the verification result is flag=%d\n\n",flag);

printf("\n\n\n--------------------------------------end  test------------------------\n\n\n");


        return 0;
}
