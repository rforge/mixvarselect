/* for Abbas: some of unecessary libraries are deleted */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/*************List of global variables********/


void Ridge_MOE_EM(double *resp,double *myX,int *constant, 
double *myinitial_alpha0, double *myinitial_alpha, double *myinitial_beta1, 
double *myinitial_beta2, double *myinitial_sigma, double *myridgepen_prop, double *myeps, double *myalpha,double *mybeta)
{
double eps=myeps[0];
int nsize=constant[0], NCOV=constant[1], ONCOV1=constant[2], ONCOV2=constant[3],NMAX=constant[4],
NCOMP=constant[5],Totaliter=constant[6], EM_maxiter=constant[7],i1,i2;

//EM_maxiter  is defined by Vahid, previousely it was set to 15*//
	
int Total1 = (NCOMP*NCOV+NCOMP),MAXTotal1=Total1+1;

double b10[NCOV],b20[NCOV];
for (i1=0;i1<NCOV;i1++)
{
	b10[i1]=myinitial_beta1[i1];b20[i1]=myinitial_beta2[i1];
}



//*double b10[NCOV]={0.9,0.5,0.5,2.1,0.7,0.6,0.5,0.5,0.5,0.5},b20[NCOV]={-0.9,0.7,0.5,3.2,0.5,0.6,0.5,0.5,0.5,0.5},
//*b30[NCOV]={-0.9,0.7,0.5,3.2,0.5,0.6,0.5,0.5,0.5,0.5};
double pi0=0.4, pi[NCOMP],sigma1;
double Betahat[NCOV][NCOMP],initbeta[NCOV][NCOMP],optlam[NCOMP];
double Hatpi[NCOMP],Hatsigma,initpi[NCOMP],Oraclbeta[NCOV][NCOMP],
Oraclpi[NCOMP],initsigma,Orcsigma,STD[Total1],alpha0[NCOMP-1],alpha[NCOV][NCOMP-1],
  Ebteda_alpha0[NCOMP-1],Ebteda_alpha[NCOV][NCOMP-1],
initalpha0[NCOMP-1],initalpha[NCOV][NCOMP],Hatalpha[NCOV][NCOMP],Hatalpha0[NCOMP-1],optlam_mix[NCOMP-1],
  Glob_Mat_pi[nsize][NCOMP],C=(log(12))/nsize,NEWmultX[NCOMP+1][nsize][NCOV];
int selection[NCOV][NCOMP],Oraclevector[NCOMP][NCOV] ,
  Oraclealpha[NCOMP-1][NCOV+1] ,
selection2[NCOV][NCOMP-1],
  ONCOV[NCOMP] ,ONCOV_alpha[NCOMP-1] ;

  
sigma1=*myinitial_sigma; /*this initializes the common sigma to be used later in its estimation*/
pi[0]=pi0;pi[1]=1-pi0;


for(i1=0;i1<NCOV;i1++)
	{
		for(i2=0;i2<NCOMP;i2++)
		{
			Betahat[i1][i2]=0;
			Oraclevector[i2][i1]=0;
		}
	}

for(i2=0;i2<NCOMP;i2++)
		{
			Hatpi[i2]=0;
			initpi[i2]=0;
			ONCOV[i2]=2;
		}

for(i2=0;i2<(NCOMP-1);i2++)
		{
			alpha0[i2]=0.5;
			Ebteda_alpha0[i2]=myinitial_alpha0[i2];
			ONCOV_alpha[i2]=3;
		}

for(i1=0;i1<NCOV;i1++)
	{
		for(i2=0;i2<(NCOMP-1);i2++)
		{
			alpha[i1][i2]=0;
			Ebteda_alpha[i1][i2]=myinitial_alpha[i1]; /* here the initial value is just one vector*/
		}
	}

for(i1=0;i1<(NCOV+1);i1++)
	{
		for(i2=0;i2<(NCOMP-1);i2++)
		{
			Oraclealpha[i2][i1]=0;
		}
	}
	
  /*changes hapenned for initialization of vectors alpha, Ebteda_alpha,  Oraclevector, Oraclealpha,
  ONCOV, ONCOV_alpha
  
/* for Abbas: some functions which were necessary in Ridge_MOE_EM are inserted here   */

void maxabs(int row, int n, double a[][MAXTotal1], int *checksin, int *maxrow)
{
	double max =0;
	int loc_maxrow, loc_checksin;
	int i;

	loc_maxrow=row;
	loc_checksin = 1;
	for (i= row; i<n ;i++) {
		if ( fabs(a[i][row]) > max) {
			max = fabs(a[i][row]);
			loc_maxrow =i;
		}

		if (max == 0) {
			/* printf(" matrix is singular \n"); */
			loc_checksin = 0;
			break;
			/* exit(0); */
		}
	}

	maxrow[0]   = loc_maxrow;
	checksin[0] = loc_checksin;
}

void gauss(int n, double a[][MAXTotal1], int *checksin)
{
	int i,j,k;
	int k_vec[1];
	double temp;

	k = 0;
	for (i=0;i<n;i++) {
	  maxabs(i, n, a, checksin, k_vec);
		k=k_vec[0];

		if (checksin[0] == 0) {
		   break;
		}
		if ( k != i)
			for ( j=i;j<n+1;j++) {
				temp = a[i][j];
				a[i][j] = a[k][j];
				a[k][j] = temp;
			} /*for j */

		for ( j=i+1;j<n;j++) {
			temp= a[j][i]/a[i][i];
			for ( k=i;k<n+1;k++)
				a[j][k]=a[j][k]-1.0*temp*a[i][k];
		} /* for j */
	} /* for i */
}/* gauss */

void backsub(int n, double a[][MAXTotal1], double *y)
{
	int i,k;
	y[n-1]= a[n-1][n]/a[n-1][n-1];
	for(i=n-2;i>= 0;i--) {
		y[i]= a[i][n];
		for(k = n-1; k > i; k--)
			y[i] = y[i] - a[i][k]*y[k];
		y[i]=y[i]/a[i][i];
	}
}

void sol(int Nelem, double IS[][MAXTotal1], double *solution, int *checksin )
{
	int i,j,k;
	gauss(Nelem, IS, checksin);
	backsub(Nelem, IS, solution);
}
  
  
/* for Abbas: END of some functions which were necessary in Ridge_MOE_EM are inserted here   */


/* for Abbas: inserts myX values to multX. multX is actually the design matrix of your code */

double multX[nsize][NCOV];
/* inserts myX values to multX*/
int icount,jcount;
for (jcount=0;jcount<NCOV;jcount++)
	{
		for (icount=0;icount<nsize;icount++)
		{
		multX[icount][jcount]=myX[nsize*jcount+icount];
		}
	}
/* for Abbas: END of inserts myX values to multX. multX is actually the design matrix of your code */

  int i=0,j=0,k1,niter1=0,l,check1[1],u,check2[1],niter2=0,Tavan;
  double sumi=0,sumi1=0,sumi2=0,mui=0,deni,beta0[NCOV][NCOMP],W[nsize][NCOMP],
    phi[nsize][NCOMP],newbeta[NCOV][NCOMP],XTWY[NCOV],XTWX[NCOV][NCOV],miui[nsize][NCOMP], 
	jamconvg2 = 0.0,
    ComMat[Total1][Total1+1],solution1[Total1],jamconvg,newpi[NCOMP],sumwi[NCOMP],
    newsigma,pii[NCOMP],sig1,Mat_pi[nsize][NCOMP],newalpha0[NCOMP-1],newalpha[NCOV][NCOMP-1],
    one_X[nsize][NCOV+1],oneXTWY[NCOV+1],oneXTWX[NCOV+1][NCOV+1],oneComMat[Total1][Total1+1],
    onesolution1[Total1],oldalpha0[NCOMP-1],oldalpha[NCOV][NCOMP-1],Max_Log_like_logistic,
    Log_like_logistic,old_oldalpha[NCOV][NCOMP-1],ridge[NCOV+1];

    
  
		for(i2=0;i2<NCOMP;i2++)
		{
			sumwi[i2]=0;
		}
    
for(i1=0;i1<nsize;i1++)
	{
		for(i2=0;i2<NCOMP;i2++)
		{
			Mat_pi[i1][i2]=0;
		}
	}

    
    
  char convg1 = 'n',convg2 = 'n';

  for(k1=0;k1<NCOMP;k1++)
    sumwi[k1] = 0.0;

  for(j=0;j<NCOV+1;j++)
   ridge[j] = *myridgepen_prop;
  
  for(i=0;i<NCOV;i++)
    beta0[i][0] = b10[i];
  for(i=0;i<NCOV;i++)
    beta0[i][1] = b20[i];
 
  for(k1=0;k1<(NCOMP-1);k1++){
    oldalpha0[k1] = Ebteda_alpha0[k1]; 
    for(j=0;j<NCOV;j++){
      old_oldalpha[j][k1]  = oldalpha[j][k1] = Ebteda_alpha[j][k1];
	 } 
  }


  for(i=0;i<nsize;i++){
    one_X[i][0] = 1.0;
    for(k1=0;k1<NCOMP;k1++)
      Mat_pi[i][k1] = 0.0;
  }	

  sumi = 0.0;
  
  for(i=0;i<nsize;i++){
   sumi1 = 0.0;
   for(k1=0;k1<(NCOMP-1);k1++){
    mui = 0.0;
    for(j=0;j<NCOV;j++)
      mui += oldalpha[j][k1]*multX[i][j];
    mui += oldalpha0[k1];
    miui[i][k1] = exp(mui); 
    sumi1 += exp(mui);    
   }   
   for(k1=0;k1<(NCOMP-1);k1++){
     Mat_pi[i][k1] = miui[i][k1]/(1.0+sumi1);
     Mat_pi[i][NCOMP-1] += Mat_pi[i][k1];
   }
  }

  for(i=0;i<nsize;i++){
    Mat_pi[i][NCOMP-1] = 1.0-Mat_pi[i][NCOMP-1];
    for(j=1;j<(NCOV+1);j++)
      one_X[i][j] = multX[i][j-1];
  }

 
  sig1 = sigma1; Tavan=0;

  while((convg1 != 'y') && (niter1 <Totaliter)){

    /****Beginning of each iteration****/

    /*******The E-step of the EM********/

    for(i=0;i<nsize;i++){
      sumi = 0.0;
      for(k1=0;k1<NCOMP;k1++){
	mui = 0.0;
	for(j=0;j<NCOV;j++)
	  mui += multX[i][j]*beta0[j][k1];
	deni = gsl_ran_gaussian_pdf(resp[i]-mui,sig1);
	phi[i][k1] = Mat_pi[i][k1]*deni;
	sumi += phi[i][k1];
      }
      for(k1=0;k1<NCOMP;k1++){
    	W[i][k1] = phi[i][k1]/sumi;
	sumwi[k1] += W[i][k1];
      }  
  }
  
    /**********End of the E-step*******/

    /*********M-step of the EM*********/

    for(k1=0;k1<NCOMP;k1++){

    /******Constructing the weigthed vector XTWY***/

   for(j=0;j<NCOV;j++){
	sumi1 = 0;
	for(i=0;i<nsize;i++)
	  sumi1 += multX[i][j]*W[i][k1]*resp[i];
	XTWY[j] = sumi1;
      }

    /*****Constructing the weighted matrix XTWX***/

      for(i=0;i<NCOV;i++){
	for(j=i;j<NCOV;j++){
	  sumi2=0;
	  for(l=0;l<nsize;l++) 
	    sumi2 += multX[l][i]*W[l][k1]*multX[l][j];
	  XTWX[j][i] = XTWX[i][j] = sumi2;
	}
      }

    /***In a system Ax=b, adding b to A as its last column**/

      for(i=0;i<NCOV;i++)
	for(j=0;j<(NCOV+1);j++)
	  if(j!=NCOV)
	    ComMat[i][j] = XTWX[i][j];
          else
	    ComMat[i][j] = XTWY[i];

    /**************************************************************/
    /*Solving the system Ax=y to get betahat in the k-th component*/ 
    /**************************************************************/

      sol(NCOV,ComMat,solution1,check1);
      for(j=0;j<NCOV;j++)
    	newbeta[j][k1] = solution1[j];
		
    } //* End of the estimation of the component-wise regression coefficients 
      

//* Newton-Raphson algorithm for updating the regression coefficients involved the mixing probabilities 


    niter2=0; convg2='n'; jamconvg2 = 0.0;Tavan=1; 

    while((convg2 != 'y') && (niter2 <EM_maxiter)){

     for(k1=0;k1<(NCOMP-1);k1++){//* For each component

    /******Constructing the weigthed vector XTWY***/

   for(j=0;j<(NCOV+1);j++){
	sumi1 = 0;
	for(i=0;i<nsize;i++)
	  sumi1 += (one_X[i][j]*(W[i][k1]-Mat_pi[i][k1]));
	 if(j==0) 
	   oneXTWY[j] = sumi1 - ridge[j]*oldalpha0[k1];
	  else
	   oneXTWY[j] = sumi1 - ridge[j]*oldalpha[j-1][k1];
      }

    /*****Constructing the weighted matrix XTWX***/

   for(i=0;i<(NCOV+1);i++){
     for(j=0;j<(NCOV+1);j++){
       sumi2=0;
       for(l=0;l<nsize;l++)
	    sumi2 += one_X[l][i]*Mat_pi[l][k1]*(1-Mat_pi[l][k1])*one_X[l][j];
	   if(i==j)
		 oneXTWX[i][j] = sumi2 + ridge[j];
	   else		
        oneXTWX[i][j] = sumi2;
     }
   }

    /***In a system Ax=b, adding b to A as its last column**/

  for(i=0;i<(NCOV+1);i++)
    for(j=0;j<(NCOV+2);j++)
      if(j!=(NCOV+1))
	oneComMat[i][j] = oneXTWX[i][j];
      else
	oneComMat[i][j] = oneXTWY[i];

    /**************************************************************/
    /*Solving the system Ax=y to get betahat in the k-th component*/ 
    /**************************************************************/

      sol(NCOV+1,oneComMat,onesolution1,check2);
	  
      for(j=0;j<NCOV+1;j++){
	if(j==0)
	  newalpha0[k1] = oldalpha0[k1] + pow(0.25,Tavan)*onesolution1[j];
	else
	  newalpha[j-1][k1] = oldalpha[j-1][k1] + pow(0.25,Tavan)*onesolution1[j];
      }
     } //* End of each component

  //*cout<<"\n";

 sumi = 0.0;
  
 for(i=0;i<nsize;i++){
   Mat_pi[i][NCOMP-1] = 0.0;
   sumi1 = 0.0;
   for(k1=0;k1<(NCOMP-1);k1++){
    mui = 0.0;
    for(j=0;j<NCOV;j++)
      mui += newalpha[j][k1]*multX[i][j];
    mui += newalpha0[k1];
    miui[i][k1] = exp(mui);
    sumi1 += exp(mui);    
   }  
   
   for(k1=0;k1<(NCOMP-1);k1++){
     Mat_pi[i][k1] = miui[i][k1]/(1.0+sumi1);
     Mat_pi[i][NCOMP-1] += Mat_pi[i][k1];
   }	
 }	
 
    for(i=0;i<nsize;i++)
     Mat_pi[i][NCOMP-1] = 1.0-Mat_pi[i][NCOMP-1];   

   jamconvg2 = 0.0;
  
   for(k1=0;k1<(NCOMP-1);k1++){
      for(j=0;j<NCOV;j++)
	   jamconvg2 += pow(newalpha[j][k1]-oldalpha[j][k1],2);

      jamconvg2 += pow(newalpha0[k1] - oldalpha0[k1],2);
    }
  
  if((jamconvg2 <= eps) || (isnan(jamconvg2)))
      convg2='y';
	  
	  
	  
  Log_like_logistic = 0.0; 

  for(k1=0;k1<NCOMP;k1++)  
   for(i=0;i<nsize;i++)
    Log_like_logistic += W[i][k1]*log(Mat_pi[i][k1]);

  if((niter2==0) && (!isnan(Log_like_logistic))){
   Max_Log_like_logistic = Log_like_logistic;
   niter2++;
   for(k1=0;k1<(NCOMP-1);k1++){
    initalpha0[k1] = oldalpha0[k1] = newalpha0[k1]; 
    for(j=0;j<NCOV;j++)
      initalpha[j][k1] = oldalpha[j][k1]  = newalpha[j][k1];
    }

   }
  else if((Log_like_logistic < Max_Log_like_logistic) || (isnan(Log_like_logistic))){ //*
      niter2=0;
      Tavan++;
      sumi = 0.0;
		
      for(k1=0;k1<(NCOMP-1);k1++){
	initalpha0[k1] = oldalpha0[k1] = newalpha0[k1];
	for(j=0;j<NCOV;j++)
	  oldalpha[j][k1]  = old_oldalpha[j][k1];
      }
		
      for(i=0;i<nsize;i++)
	for(k1=0;k1<NCOMP;k1++)
	  Mat_pi[i][k1] = 0.0;
  
  for(i=0;i<nsize;i++){
   sumi1 = 0.0;
   for(k1=0;k1<(NCOMP-1);k1++){
    mui = 0.0;
    for(j=0;j<NCOV;j++)
      mui += old_oldalpha[j][k1]*multX[i][j];
    mui += oldalpha0[k1]; 
    miui[i][k1] = exp(mui); 
    sumi1 += exp(mui);    
   }   
   for(k1=0;k1<(NCOMP-1);k1++){
     Mat_pi[i][k1] = miui[i][k1]/(1.0+sumi1);
     Mat_pi[i][NCOMP-1] += Mat_pi[i][k1];
   }
  }

  for(i=0;i<nsize;i++)
    Mat_pi[i][NCOMP-1] = 1.0-Mat_pi[i][NCOMP-1];
  }
  else{ 
   niter2++; 
   for(k1=0;k1<(NCOMP-1);k1++){
    initalpha0[k1] = oldalpha0[k1] = newalpha0[k1]; 
    for(j=0;j<NCOV;j++)
      initalpha[j][k1] = oldalpha[j][k1]  = newalpha[j][k1];
    }
   }
   
   
} //* End of the Newton-Raphson algorithm   


    sumi = 0;
    for(i=0;i<nsize;i++){
      for(k1=0;k1<NCOMP;k1++){
	   mui = 0;
	   for(j=0;j<NCOV;j++)
	    mui += multX[i][j]*newbeta[j][k1];
	  sumi += W[i][k1]*pow(resp[i]-mui,2);  
      }    
    }     

    newsigma = sqrt(sumi/nsize);

    /*****End of the M-step of the EM*****/
	
    jamconvg = 0.0;
    niter1++;
	
    for(k1=0;k1<NCOMP;k1++){
      for(j=0;j<NCOV;j++)
	   jamconvg += pow(newbeta[j][k1]-beta0[j][k1],2);

     //* jamconvg += (newpi[k1]-pii[k1])*(newpi[k1]-pii[k1]);
    }

    jamconvg += pow(newsigma-sig1,2);

    jamconvg += jamconvg2;

    if((jamconvg <= eps) || (isnan(jamconvg)))
      convg1='y';

    for(k1=0;k1<NCOMP;k1++){
      for(j=0;j<NCOV;j++)
	initbeta[j][k1] = beta0[j][k1] = newbeta[j][k1];
      sumwi[k1] = 0;
    }

    initsigma = sig1 = newsigma;
    
    /*******End of each iteration*******/
  }
  
  sumi1 = 0.0;

/* for Abbas: I am putting the results in a pointer to bring them into R */

for (jcount=0;jcount<NCOV;jcount++)
{
	for (icount=0;icount<NCOMP;icount++)
	{
		mybeta[jcount*NCOMP+icount]=initbeta[jcount][icount];
		}
}		

for (jcount=0;jcount<NCOV;jcount++)
{
		myalpha[jcount]=initalpha[jcount][0];
}		

/* for Abbas: if never ever you wanted to print something
icount=1;
for (jcount=0;jcount<NCOV;jcount++)
	{
		printf("X[i=%d, j=%d]=%f \n",icount+1,jcount+1,multX[icount][jcount]);
	}

*/


/* for Abbas: END of I am putting the results in a pointer to bring them into R */


}//*End of the function

/* ***************************************************************** */
