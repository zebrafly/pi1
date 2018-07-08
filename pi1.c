#include "pi1.h"

int pi1_init()
{
	// Set up the context
	if (core_init() != STS_OK)
	{
		core_clean();
		return 1;
	}
	if (pc_param_set_any() != STS_OK)
	{
		THROW(ERR_NO_CURVE);
		core_clean();
		return 1;
	}
	return 0;
}

int pi1_close()
{
	core_clean();

	return 0;
}


//the pi1_keygen algorithm: generating the evaluation key and public verification key of the scheme
int pi1_keygen(vc_k *vck, fq_mat_t F, lhe_par *par)
{
   fq_mat_init(vck->F,vck->m,vck->d,par->ctx);
   fq_mat_set(vck->F,F,par->ctx);

   vck->G=malloc(sizeof(g1_t)*vck->m);
   vck->T=malloc(sizeof(g1_t)*vck->d);

   fq_mat_t R;
   fq_mat_init(R,1,vck->m,par->ctx);
   bn_t zb;
   bn_new(zb);
   fmpz_t zf;
   fmpz_init(zf);
   fq_t zq;
   fq_init(zq,par->ctx);


   for (int i=0;i<vck->m;i++)
   {
      bn_rand_mod(zb,par->q);
      bn2fmpz(zf,zb);
      fq_set_fmpz(zq,zf,par->ctx);
      fq_mat_entry_set(R,0,i,zq,par->ctx);
      g1_new(vck->G[i]);
      g1_mul(vck->G[i],par->g,zb);
   }


   printf("\n\n this is the random vector R:\n\n");
   fq_mat_print_pretty(R,par->ctx);


   fq_mat_t RF;
   fq_mat_init(RF,1,vck->d,par->ctx);
   fq_mat_mul(RF,R,F,par->ctx);

   printf("\n\n this is the random vector RF:\n\n");
   fq_mat_print_pretty(RF,par->ctx);


   for (int j=0;j<vck->d;j++)
   {
      bn_read_str(zb,fq_get_str_pretty(fq_mat_entry(RF,0,j),par->ctx),strlen(fq_get_str_pretty(fq_mat_entry(RF,0,j),par->ctx)),10);
      g1_new(vck->T[j]);
      g1_mul(vck->T[j],par->g,zb);
   }

   //release the memory
   fq_mat_clear(R,par->ctx);
   fq_mat_clear(RF,par->ctx);
   bn_free(zb);
   fmpz_clear(zf);
   fq_clear(zq,par->ctx);
   return 0;
}


//generate a random matrix over fq
int fq_mat_randz(fq_mat_t F, int m, int d, lhe_par *par)
{
   bn_t zb;
   bn_new(zb);
   fmpz_t zf;
   fmpz_init(zf);
   fq_t zq;
   fq_init(zq,par->ctx);

   for (int i=0;i<m;i++)
   {
     for (int j=0;j<d;j++)
     {
         bn_rand_mod(zb,par->p);
         bn2fmpz(zf,zb);
         fq_set_fmpz(zq,zf,par->ctx);
         fq_mat_entry_set(F,i,j,zq,par->ctx);
     }
   }
   bn_free(zb);
   fmpz_clear(zf);
   fq_clear(zq,par->ctx);
   return 0;
}

int pi1_eval(fq_t c, bn_t *C, fq_t *LL, int d, lhe_par *par)
{

    fq_mat_t X;
    fq_mat_init(X,d,1,par->ctx);
    fmpz_t zf;
    fmpz_init(zf);
    fq_t zq;
    fq_init(zq,par->ctx);
    for (int i=0;i<d;i++)
    {
       bn2fmpz(zf,C[i]);
       fq_set_fmpz(zq,zf,par->ctx);
       fq_mat_entry_set(X,i,0,zq,par->ctx);
    }

    fq_mat_t F;
    fq_mat_init(F,1,d,par->ctx);

    for(int y=0;y<d;y++)
    {
      fq_mat_entry_set(F,0,y,LL[y],par->ctx);
    }

    fq_mat_t FX;
    fq_mat_init(FX,1,1,par->ctx);
    fq_mat_mul(FX,F,X,par->ctx);
    fq_set(c,fq_mat_entry(FX,0,0),par->ctx);
    // fq2bn(c,fq_mat_entry(FX,0,0),par->ctx);
    fq_mat_clear(X,par->ctx);
    fq_mat_clear(F,par->ctx);
    fq_mat_clear(FX,par->ctx);
    fmpz_clear(zf);
    fq_clear(zq,par->ctx);



   //   for (int i=0; i<d; i++)
   //   {
   //       fq_t temp;
   //       fq_init(temp,par->ctx);
   //       fq_t tempnew;
   //       fq_init(tempnew,par->ctx);
   //       fmpz_t zf;
   //       fmpz_init(zf);
   //       fq_t zq;
   //       fq_init(zq,par->ctx);
   //       bn2fmpz(zf,C[i]);
   //       fq_set_fmpz(zq,zf,par->ctx);
   //       // c=zq*L[i];
   //       fq_mul(temp,LL[i],zq,par->ctx);
   //       if(i==0)
   //       {
   //         printf("%s\n","aaa" );
   //         fq_set(c,temp,par->ctx);
   //         printf("%s\n","aaa" );
   //       }
   //       else
   //       {
   //         fq_set(tempnew,c,par->ctx);
   //         fq_add(c,tempnew,temp,par->ctx);
   //       }
   // }
    return 0;
 }


   //release the memory


//the pi1_comp algorithm
int pi1_comp(fq_t *nv, vc_k *vck, vc_p *vcp, lhe_par *par)
{
  fq_t *L=malloc(sizeof(fq_t)*vck->d);
  for (int j=0;j<vck->d;j++)
  {
    fq_init(L[j],par->ctx);
  }

  for (int i=0;i<vck->m;i++)
  {
     for (int j=0;j<vck->d;j++)
     {
       fq_set(L[j],fq_mat_entry(vck->F,i,j),par->ctx);
     }

     pi1_eval(nv[i],vcp->x,L,vck->d,par);

  }

  free(L);
  return 0;
}

//the pi1_vrfy algorithm
int pi1_vrfy(vc_k *vck, vc_p *vcp, fq_t *nv, lhe_par *par)
{
  int flag=-1;
  g1_t left;
  g1_new(left);
  g1_t right;
  g1_new(right);
  for(int i=0;i<vck->m;i++)
  {
    g1_t temp;
    g1_new(temp);
    bn_t aaa;
    bn_new(aaa);
    fq2bn(aaa,nv[i],par->ctx);
    g1_mul(temp,vck->G[i],aaa);
    if(i==0)
    {
      g1_copy(left,temp);
    }
    else
    {
      g1_add(left,left,temp);
    }
  }
  for(int j=0;j<vck->d;j++)
  {
    g1_t temp;
    g1_new(temp);
    g1_mul(temp,vck->T[j],vcp->x[j]);
    if(j==0)
    {
      g1_copy(right,temp);
    }
    else
    {
      g1_add(right,right,temp);
    }
  }
  printf("\n\n the verification tag right is: \n\n");
  g1_print(right);
  printf("\n\n the verification tag left is: \n\n");
  g1_print(left);
  flag=g1_cmp(left,right);
  if (flag == 2)
  {
     return 2;
  }

  return 0;
}
