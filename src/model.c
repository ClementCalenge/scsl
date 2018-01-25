#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>


typedef struct {
    int I;
    int J;
    int K;
    double nu;
    SEXP nuijl;
    int nT;
    SEXP Y;
    SEXP Ykj;
    SEXP Yik;
    SEXP Yij;
    SEXP Sj;
    SEXP lVj;
    SEXP Vj;
} Pouroptim;




void compteur(int i)
{
    if (i < 10)
	Rprintf("\b");
    if ((i > 9)&&(i<100))
	Rprintf("\b\b");
    if ((i > 99)&&(i<1000))
	Rprintf("\b\b\b");
    if ((i > 999)&&(i<10000))
	Rprintf("\b\b\b\b");
    if ((i > 9999)&&(i<100000))
	Rprintf("\b\b\b\b\b");
    if ((i > 99999)&&(i<1000000))
	Rprintf("\b\b\b\b\b\b");
    if ((i > 999999)&&(i<10000000))
	Rprintf("\b\b\b\b\b\b\b");
    if ((i > 9999999)&&(i<100000000))
	Rprintf("\b\b\b\b\b\b\b\b");
    Rprintf("%i", i+1);
}



double calcVrais(SEXP EpsViv, SEXP PiViv, SEXP Alp, int I, int J, int K, SEXP lVj, SEXP Y,
		 SEXP Sj, int nT)
{
    SEXP loglambdajki;
    int i, j, k;
    double sum1, sum2, Vrais;
    
    PROTECT(loglambdajki = allocVector(REALSXP, J*K*I));
    for (j = 0; j < J; j++) {
	for (i = 0; i < I; i++) {
	    REAL(loglambdajki)[j + (0 * J) + (i * K * J)] = REAL(Alp)[i + I*j] + REAL(lVj)[j];
	    REAL(loglambdajki)[j + (1 * J) + (i * K * J)] = REAL(Alp)[i + I*j] + 
		REAL(EpsViv)[j] + REAL(PiViv)[i];
	}
    }
    
    sum2 = 0.0;
    sum1 = 0.0;
    for (j = 0; j < J; j++) {
	for (k = 0; k < K; k++) {
	    for (i = 0; i < I; i++) {
		sum2 += (REAL(loglambdajki)[j+k*J+i*K*J] + log(REAL(Sj)[j])) * 
		    REAL(Y)[j+k*J+i*K*J];
		sum1 += exp(REAL(loglambdajki)[j+k*J+i*K*J])*REAL(Sj)[j];
	    }
	}
    }
    sum1 = sum1 * ((double) nT);
    Vrais = sum1-sum2;
    UNPROTECT(1);
    
    return(Vrais);
}


SEXP calcVraisr(SEXP EpsViv, SEXP PiViv, SEXP Alp, SEXP Ir, SEXP Jr, 
	       SEXP Kr, SEXP lVj, SEXP Y, SEXP Sj, SEXP nTr)

{
    SEXP resu;
    PROTECT(resu=allocVector(REALSXP, 1));
    
    REAL(resu)[0] = calcVrais(EpsViv, PiViv, Alp, INTEGER(Ir)[0], INTEGER(Jr)[0], 
			      INTEGER(Kr)[0], lVj, Y, Sj, INTEGER(nTr)[0]);
    UNPROTECT(1);
    return(resu);
}














SEXP EstimMaxVrais(SEXP AlpInit, SEXP PiInit, SEXP EpsInit,
		   int nT, int I, int J, SEXP Sj, SEXP lVj, 
		   SEXP Vj, SEXP Y,
		   SEXP Ykj, SEXP Yik, SEXP Yij, 
		   int maxIter, double stopCrit, double VraisCompl)
{
    double vrais1, vrais;
    int a, i, j;
    SEXP d, sumcold, sumrowd, resu, resuvrais, Alp, Pi, Eps;
    PROTECT(d = allocVector(REALSXP, I*J));
    PROTECT(Alp = allocVector(REALSXP, I*J));
    PROTECT(Pi = allocVector(REALSXP, I));
    PROTECT(Eps = allocVector(REALSXP, J));
    PROTECT(sumcold = allocVector(REALSXP, J));
    PROTECT(sumrowd = allocVector(REALSXP, I));
    PROTECT(resuvrais = allocVector(REALSXP, 1));
    PROTECT(resu = allocVector(VECSXP, 4));
    
    /* Copie des valeurs initiales */
    for (i=0; i< I*J; i++) {
	REAL(Alp)[i] = REAL(AlpInit)[i];
    }
    for (i=0; i< I; i++) {
	REAL(Pi)[i] = REAL(PiInit)[i];
    }
    for (i=0; i< J; i++) {
	REAL(Eps)[i] = REAL(EpsInit)[i];
    }

    /* Vraisemblance initiale */
    vrais1 = calcVrais(EpsInit, PiInit, AlpInit, I, J, 2, lVj, Y,
		       Sj, nT);
    vrais = vrais1;
    
    for (a = 0; a<maxIter; a++) {
	vrais = vrais1;
	/* Actualisation des EpsViv */
	for (j=0; j<J; j++) {
	    REAL(d)[j*I] = exp(REAL(Alp)[j*I])*REAL(Sj)[j];
	}
	for (i=1; i<I; i++) {
	    for (j=0; j<J; j++) {
		REAL(d)[i+j*I] = exp(REAL(Alp)[i+j*I]+REAL(Pi)[i])*REAL(Sj)[j];
	    }
	}
	for (j=0; j<J; j++) {
	    REAL(sumcold)[j] = 0;
	    for (i = 0; i<I; i++) {
		REAL(sumcold)[j] += REAL(d)[i+j*I];
	    }
	}
	for (j=0; j<J; j++) {
	    REAL(Eps)[j] = log(REAL(Ykj)[1+2*j]/(((double) nT)*REAL(sumcold)[j]));
	    if (REAL(Eps)[j] < -30.0)
		REAL(Eps)[j] = -30.0; /* xx ajout xx */

	}
	/* Actualisation des PiViv */
	for (i=0; i<I; i++) {
	    for (j=0; j<J; j++) {
		REAL(d)[i+j*I] = exp(REAL(Alp)[i+j*I]+REAL(Eps)[j])*REAL(Sj)[j];
	    }
	}
	for (i=0; i<I; i++) {
	    REAL(sumrowd)[i] = 0;
	    for (j = 0; j<J; j++) {
		REAL(sumrowd)[i] += REAL(d)[i+j*I];
	    }
	}
	for (i=0; i<I; i++) {
	    REAL(Pi)[i] = log(REAL(Yik)[i+1*I]/(((double) nT)*REAL(sumrowd)[i]));
	}
	REAL(Pi)[0] = 0;
	/* Actualisation des Alp */
	for (i=0; i<I; i++) {
	    for (j=0; j<J; j++) {
		REAL(d)[i+j*I] = (REAL(Vj)[j] + 
		    exp(REAL(Eps)[j])*exp(REAL(Pi)[i]))*REAL(Sj)[j];
	    }
	}
	for (i=0; i<I; i++) {
	    for (j=0; j<J; j++) {
		REAL(Alp)[i+j*I] = log(REAL(Yij)[i+I*j]/(((double) nT)*REAL(d)[i+j*I]));
		if (REAL(Yij)[i+I*j] < 0.5)
		    REAL(Alp)[i+j*I] = -30;
	    }
	}
	vrais1 = calcVrais(Eps, Pi, Alp, I, J, 2, lVj, Y,
			   Sj, nT);
	if ((vrais-vrais1)/fabs(vrais) < stopCrit)
	    break;
    }
    REAL(resuvrais)[0] = vrais;
    SET_VECTOR_ELT(resu, 0, resuvrais);
    SET_VECTOR_ELT(resu, 1, Alp);
    SET_VECTOR_ELT(resu, 2, Eps);
    SET_VECTOR_ELT(resu, 3, Pi);
        
    UNPROTECT(8);

    return(resu);
}



SEXP ModeliseBrut(SEXP Esp, SEXP Dep, SEXP Sta, SEXP Y, SEXP Vjt, 
		  SEXP Ir, SEXP Jr, SEXP Kr, SEXP Sjt, SEXP nTr,
		  SEXP maxIter, SEXP stopCrit)
{
    SEXP Vj, lVj, Sj, loglambdajki, resu, AlpInit, EpsInit, PiInit;
    SEXP Ykj, Yik, Yij;
    int i,j,k, I, J, K, nT;
    double sum2, sum1, VraisCompl;
    
    /* Définition des constantes */
    I = INTEGER(Ir)[0];
    J = INTEGER(Jr)[0];
    K = INTEGER(Kr)[0];
    nT = INTEGER(nTr)[0];
    
    /* Mise en forme des Vj et Sj */
    PROTECT(Vj = allocVector(REALSXP, J));
    PROTECT(Sj = allocVector(REALSXP, J));
    PROTECT(lVj = allocVector(REALSXP, J));    
    for (j = 0; j < J; j++) {
	REAL(Vj)[j] = REAL(Vjt)[j];
	REAL(Sj)[j] = REAL(Sjt)[j];
	REAL(lVj)[j] = log(REAL(Vjt)[j]);
    }
    
    /* Calcul des statistiques suffisantes */
    PROTECT(Ykj = allocVector(REALSXP, K*J));
    PROTECT(Yik = allocVector(REALSXP, I*K));
    PROTECT(Yij = allocVector(REALSXP, I*J));    
    /* Initialisation à zéro */
    for (j=0; j<J; j++){
	for (k=0; k<K; k++){
	    REAL(Ykj)[k+j*K] = 0;
	}
    }
    for (i=0; i<I; i++){
	for (k=0; k<K; k++){
	    REAL(Yik)[i+k*I] = 0;
	}
    }
    for (j=0; j<J; j++){
	for (i=0; i<I; i++){
	    REAL(Yij)[i+j*I] = 0;
	}
    }
    /* Calcul des totaux */
    for (i=0; i<I; i++){
	for (j=0; j<J; j++){
	    for (k=0; k<K; k++){
		REAL(Ykj)[k+j*K] += REAL(Y)[j+k*J+i*K*J];
		REAL(Yik)[i+k*I] += REAL(Y)[j+k*J+i*K*J];
		REAL(Yij)[i+j*I] += REAL(Y)[j+k*J+i*K*J];
	    }
	}
    }
	

    /* calcul de la vraisemblance pour le modèle complet */
    PROTECT(loglambdajki = allocVector(REALSXP, J*I*K));
    for (i = 0; i < (I*J*K); i++) {
	if (REAL(Y)[i]> 0.5) {
	    REAL(loglambdajki)[i] = log(REAL(Y)[i]);
	} else {
	    REAL(loglambdajki)[i] = -30.0;
	}
    }
    
    sum2 = 0.0;
    for (i = 0; i < (I*J*K); i++) {
	sum2 += REAL(loglambdajki)[i] * REAL(Y)[i];
    }
    sum1 = 0.0;
    for (j = 0; j < J; j++) {
	for (k = 0; k < K; k++) {
	    for (i = 0; i < I; i++) {
		sum1 += REAL(Y)[j+k*J+i*K*J]*REAL(Sj)[j];
	    }
	}
    }
    sum1 = sum1 * ((double) nT);
    VraisCompl = sum1-sum2;

    /* Valeurs initiales des paramètres */
    PROTECT(AlpInit = allocVector(REALSXP, J*I));
    PROTECT(EpsInit = allocVector(REALSXP, J));
    PROTECT(PiInit = allocVector(REALSXP, I));
    for (i = 0; i < I; i++) {
	for (j = 0; j < J; j++) {
	    REAL(AlpInit)[i+j*I] = REAL(loglambdajki)[j+K*J*i] - REAL(lVj)[j];
	}
    }
    for (j = 0; j < J; j++) {
	REAL(EpsInit)[j] = REAL(loglambdajki)[j+J] - REAL(AlpInit)[j*I];
    }
    REAL(PiInit)[0] = 0.0;
    for (i = 1; i < I; i++) {
	REAL(PiInit)[i] = 0.0;
	for (j = 0; j < J; j++) {
	    REAL(PiInit)[i] += REAL(loglambdajki)[j+J+K*J*i] - 
		REAL(AlpInit)[i+j*I] - REAL(EpsInit)[j];
	}
	REAL(PiInit)[i] = REAL(PiInit)[i]/((double) J);
    }
    
    PROTECT(resu = EstimMaxVrais(AlpInit, PiInit, EpsInit,
				 nT, I, J, Sj, lVj, 
				 Vj, Y,
				 Ykj, Yik, Yij, 
				 INTEGER(maxIter)[0], REAL(stopCrit)[0],
				 VraisCompl));
    

    UNPROTECT(11);
    return(resu);

}



double calcSum3(SEXP nuijl, int I, int J, SEXP Alp)
{
    double Sijl;
    int i, j, l;
    
    Sijl = 0.0;
    
    for (i = 0; i < I; i++) {
	for (j = 0; j < J; j++) {
	    for (l = 0; l <= j; l++) {
		Sijl += REAL(nuijl)[i+j*I+l*I*J]*R_pow(REAL(Alp)[i+j*I] - REAL(Alp)[i+l*I], 2.0);
	    }
	}
    }
    return(Sijl);
}



/* Renvoie moins la log-vraisemblance pénalisée (plus faible = mieux) */
double calcVraisPen(double nu, SEXP nuijl, int nT, SEXP Y, SEXP Ykj, 
		    SEXP Yik, SEXP Yij, int I, int J, int K, SEXP EpsViv,
		    SEXP PiViv, SEXP Alp, SEXP Sj, SEXP lVj)
{
    SEXP loglambdajki;
    int i, j, k;
    double sum2, sum1, sum3, Vrais;
    
    PROTECT(loglambdajki = allocVector(REALSXP, J*K*I));
    for (j = 0; j < J; j++) {
	for (i = 0; i < I; i++) {
	    REAL(loglambdajki)[j + (0 * J) + (i * K * J)] = REAL(Alp)[i + I*j] + REAL(lVj)[j];
	    REAL(loglambdajki)[j + (1 * J) + (i * K * J)] = REAL(Alp)[i + I*j] + 
		REAL(EpsViv)[j] + REAL(PiViv)[i];
	}
    }
    
    sum2 = 0.0;
    sum1 = 0.0;
    for (j = 0; j < J; j++) {
	for (k = 0; k < K; k++) {
	    for (i = 0; i < I; i++) {
		sum1 += exp(REAL(loglambdajki)[j+k*J+i*K*J])*REAL(Sj)[j];
		sum2 += (REAL(loglambdajki)[j+k*J+i*K*J] + log(REAL(Sj)[j])) * REAL(Y)[j+k*J+i*K*J];
	    }
	}
    }
    sum1 = sum1 * ((double) nT);
    sum3 = calcSum3(nuijl, I, J, Alp);

    Vrais = sum1-sum2+nu*sum3;
    UNPROTECT(1);
    return(Vrais);
}



SEXP calcSumd3(SEXP nuijl, int I, int J, SEXP Alp)
{
    SEXP Sijl;
    int i, j, l;
    
    PROTECT(Sijl = allocVector(REALSXP, I*J));
    for (l = 0; l < I*J; l++) {
	REAL(Sijl)[l] = 0.0;
    }
    
    for (l = 0; l < J; l++) {
 	for (i = 0; i < I; i++) {
	    for (j = 0; j < J; j++) {
		REAL(Sijl)[i+j*I] += REAL(nuijl)[i+j*I+l*I*J]*(REAL(Alp)[i+j*I] - REAL(Alp)[i+l*I]);
	    }
	}
    }
    UNPROTECT(1);
    return(Sijl);
}


SEXP scaleok(SEXP mat, int I, int J)
{
    double sdd, sdd2;
    SEXP mats;
    int i, j;
    
    PROTECT(mats=allocVector(REALSXP, I*J));

    for (j = 0; j < J; j++) {
	sdd = 0.0;
	for (i = 0; i < I; i++) {
	    sdd += (R_pow(REAL(mat)[i+j*I], 2.0));
	}
	sdd2 = sqrt(sdd);
//	Rprintf("suu:%f\n",sdd2);
	for (i = 0; i < I; i++) {
	    REAL(mats)[i+j*I] = REAL(mat)[i+j*I]/sdd2;
	}
    }
    
    UNPROTECT(1);
    return(mats);
}

SEXP transposeok(SEXP mat, int I, int J)
{
    int i,j;
    SEXP matt;
    PROTECT(matt=allocVector(REALSXP, I*J));
    for (j = 0; j < J; j++) {
	for (i = 0; i < I; i++) {
	    REAL(matt)[j+J*i] = REAL(mat)[i+j*I];
	}
    }
    UNPROTECT(1);
    return(matt);    
}




SEXP EstimMaxVraisPenl2(SEXP AlpInit, SEXP PiInit, SEXP EpsInit,
			int nT, int I, int J, SEXP Sj, SEXP lVj, 
			SEXP Vj, SEXP Y,
			SEXP Ykj, SEXP Yik, SEXP Yij, 
			int maxIter, double stopCrit, double nu,
			SEXP nuijl, double h, int verb)
{
    double vrais1, vrais, vraisPenInit;
    int a, i, j;
    SEXP d3, d, sumcold, sumrowd, resu, resuvrais, resupen, Alp, Pi, Eps, dd, dt, dt2;
    SEXP resuvraispeninit, NIter, resuh;

    PROTECT(d = allocVector(REALSXP, I*J));
    PROTECT(Alp = allocVector(REALSXP, I*J));
    PROTECT(Pi = allocVector(REALSXP, I));
    PROTECT(Eps = allocVector(REALSXP, J));
    PROTECT(sumcold = allocVector(REALSXP, J));
    PROTECT(sumrowd = allocVector(REALSXP, I));
    PROTECT(resuvrais = allocVector(REALSXP, 1));
    PROTECT(resuvraispeninit = allocVector(REALSXP, 1));
    PROTECT(resupen = allocVector(REALSXP, 1));
    PROTECT(resuh = allocVector(REALSXP, 1));
    PROTECT(NIter = allocVector(INTSXP, 1));
    PROTECT(resu = allocVector(VECSXP, 8));
    vrais=0.0;

    /* Copie des valeurs initiales */
    for (i=0; i< I*J; i++) {
	REAL(Alp)[i] = REAL(AlpInit)[i];
    }
    for (i=0; i< I; i++) {
	REAL(Pi)[i] = REAL(PiInit)[i];
    }
    REAL(Pi)[0] = 0;
    for (i=0; i< J; i++) {
	REAL(Eps)[i] = REAL(EpsInit)[i];
    }
    
    /* Vraisemblance initiale */
    vrais1 = calcVraisPen(nu, nuijl, nT, Y, Ykj, 
			  Yik, Yij, I, J, 2, Eps,
			  Pi, Alp, Sj, lVj);
    vraisPenInit = vrais1;
    
    if (verb)
	Rprintf("Initial Value of the penalized likelihood: %f\n", vrais1);
    
    for (a = 0; a < maxIter; a++) {
	if (verb)
	    Rprintf("Iteration %i: ", a);

	vrais = vrais1;
	R_CheckUserInterrupt();
	/* Actualisation des EpsViv */
	for (i=0; i<I; i++) {
	    for (j=0; j<J; j++) {
		REAL(d)[i+j*I] = exp(REAL(Alp)[i+j*I]+REAL(Pi)[i])*REAL(Sj)[j];
	    }
	}
	for (j=0; j<J; j++) {
	    REAL(sumcold)[j] = 0;
	    for (i = 0; i<I; i++) {
		REAL(sumcold)[j] += REAL(d)[i+j*I];
	    }
	}
	for (j=0; j<J; j++) {
	    REAL(Eps)[j] = log(REAL(Ykj)[1+2*j]/(((double) nT)*REAL(sumcold)[j]));
	    if (REAL(Eps)[j] < -30.0)
		REAL(Eps)[j] = -30.0; /* xx ajout xx */
	}
	/* Actualisation des PiViv */
	for (i=0; i<I; i++) {
	    for (j=0; j<J; j++) {
		REAL(d)[i+j*I] = exp(REAL(Alp)[i+j*I]+REAL(Eps)[j])*REAL(Sj)[j];
	    }
	}
	for (i=0; i<I; i++) {
	    REAL(sumrowd)[i] = 0;
	    for (j = 0; j<J; j++) {
		REAL(sumrowd)[i] += REAL(d)[i+j*I];
	    }
	}
	for (i=0; i<I; i++) {
	    REAL(Pi)[i] = log(REAL(Yik)[i+1*I]/(((double) nT)*REAL(sumrowd)[i]));
	}
	REAL(Pi)[0] = 0;
	/* Actualisation des Alp */
	PROTECT(d3=calcSumd3(nuijl, I, J, Alp));
	for (i=0; i<I; i++) {
	    for (j=0; j<J; j++) {
		REAL(d)[i+j*I] = ((double) nT)*REAL(Sj)[j]*exp(REAL(Alp)[i+j*I])*
		    (REAL(Vj)[j]+exp(REAL(Eps)[j])*exp(REAL(Pi)[i])) - 
		    REAL(Yij)[i+j*I] + 2.0*nu*REAL(d3)[i+j*I];
	    }
	}
	PROTECT(dt = transposeok(d, I, J));
	PROTECT(dt2 = scaleok(dt, J, I));
	PROTECT(dd = transposeok(dt2, J, I));
	for (i=0; i<I; i++) {
	    for (j=0; j<J; j++) {
		REAL(Alp)[i+j*I] = REAL(Alp)[i+j*I] - h * REAL(dd)[i+j*I];
	    }	
	}
	vrais1 = calcVraisPen(nu, nuijl, nT, Y, Ykj, 
			      Yik, Yij, I, J, 2, Eps,
			      Pi, Alp, Sj, lVj);

	if (verb)
	    Rprintf("%f\n", vrais1);
	
	if (vrais-vrais1 <0)
	    h = h/2.0;
	
//	Rprintf("vrais: %f %f\n", vrais, vrais1);
//	Rprintf("%f\n", h);

	UNPROTECT(4);
	if (((vrais-vrais1) >0)&&((vrais-vrais1)/fabs(vrais) < stopCrit)) {
//	    if (h > 0.0000001)
		break;
	}
    }
    if (verb) {
	if (a == maxIter) {
	    Rprintf("Failed to converge. Returning last value\n");
	} else {
	    Rprintf("Converged.\n");
	}
    }

    
    REAL(resuvrais)[0] = vrais;
    SET_VECTOR_ELT(resu, 0, resuvrais);
    SET_VECTOR_ELT(resu, 1, Alp);
    SET_VECTOR_ELT(resu, 2, Eps);
    SET_VECTOR_ELT(resu, 3, Pi);
    REAL(resupen)[0] = calcSum3(nuijl, I, J, Alp);
    SET_VECTOR_ELT(resu, 4, resupen);
    REAL(resuvraispeninit)[0] = vraisPenInit;
    SET_VECTOR_ELT(resu, 5, resuvraispeninit);
    INTEGER(NIter)[0] = a;
    REAL(resuh)[0] = h;
    SET_VECTOR_ELT(resu, 6, NIter);
    SET_VECTOR_ELT(resu, 7, resuh);
    
    UNPROTECT(12);

    return(resu);
}




SEXP EstimMaxVraisPenl2alt(SEXP AlpInit, SEXP PiInit, SEXP EpsInit,
			   int nT, int I, int J, SEXP Sj, SEXP lVj, 
			   SEXP Vj, SEXP Y,
			   SEXP Ykj, SEXP Yik, SEXP Yij, 
			   int maxIter, double stopCrit, double nu,
			   SEXP nuijl, double h, int verb)
{
    double vrais1, vrais, vraisPenInit, su;
    int a, i, j;
    SEXP d3, d, de, dp, sumcold, sumrowd, resu, resuvrais, resupen, Alp, Pi, Eps, dd, dt, dt2;
    SEXP resuvraispeninit, NIter, resuh;
    
    PROTECT(d = allocVector(REALSXP, I*J));
    PROTECT(de = allocVector(REALSXP, J));
    PROTECT(dp = allocVector(REALSXP, I));
    PROTECT(Alp = allocVector(REALSXP, I*J));
    PROTECT(Pi = allocVector(REALSXP, I));
    PROTECT(Eps = allocVector(REALSXP, J));
    PROTECT(sumcold = allocVector(REALSXP, J));
    PROTECT(sumrowd = allocVector(REALSXP, I));
    PROTECT(resuvrais = allocVector(REALSXP, 1));
    PROTECT(resuvraispeninit = allocVector(REALSXP, 1));
    PROTECT(resupen = allocVector(REALSXP, 1));
    PROTECT(resuh = allocVector(REALSXP, 1));
    PROTECT(NIter = allocVector(INTSXP, 1));
    PROTECT(resu = allocVector(VECSXP, 8));
    vrais = 0.0;
    
    /* Copie des valeurs initiales */
    for (i=0; i< I*J; i++) {
	REAL(Alp)[i] = REAL(AlpInit)[i];
    }
    for (i=0; i< I; i++) {
	REAL(Pi)[i] = REAL(PiInit)[i];
    }
    REAL(Pi)[0] = 0.0;
    
    for (i=0; i< J; i++) {
	REAL(Eps)[i] = REAL(EpsInit)[i];
    }
    
    /* Vraisemblance initiale */
    vrais1 = calcVraisPen(nu, nuijl, nT, Y, Ykj, 
			  Yik, Yij, I, J, 2, Eps,
			  Pi, Alp, Sj, lVj);
    vraisPenInit = vrais1;

    if (verb)
	Rprintf("Initial Value of the likelihood: %f\n", vrais1);

    
    for (a = 0; a < maxIter; a++) {
	if (verb)
	    Rprintf("Updating... ");

	vrais = vrais1;
	R_CheckUserInterrupt();

	

	/* Actualisation des EpsViv */
	for (j=0; j<J; j++) {
	    REAL(de)[j] = 0.0;
	}	
	for (j=0; j<J; j++) {
	    for (i=0; i<I; i++) {
		REAL(de)[j] += ((double) nT)*REAL(Sj)[j]*exp(REAL(Alp)[i+j*I] + REAL(Eps)[j] + REAL(Pi)[i]);
	    }
	    REAL(de)[j] = REAL(de)[j] - REAL(Ykj)[1+j*2];
	}
	su = 0.0;
	for (j=0; j<J; j++) {
	    su += R_pow(REAL(de)[j],2.0);
	}
	/* Si on recherche la solution du maximum de vraisemblance (nu = 0), alors le gradient
	   est nul et su, la norme de ce vecteur est nulle aussi. Dans ce cas, on
	   divise 0.0 par 0, et ça fait NaN. Donc dans ce cas de figure, on fixe la norme à 1
	   et pis voila
	*/
	if (su < 1e-10) {
	    su = 1;
	} else {
	    su = sqrt(su/((double) (J)));
	}
	for (j=0; j<J; j++) {
	    REAL(Eps)[j] = REAL(Eps)[j] - h * REAL(de)[j]/su;
	}	

	/* Actualisation des PiViv */
	for (i=0; i<I; i++) {
	    REAL(dp)[i] = 0.0;
	}
	for (i=0; i<I; i++) {
	    for (j=0; j<J; j++) {
		REAL(dp)[i] += ((double) nT)*REAL(Sj)[j]*exp(REAL(Alp)[i+j*I] + REAL(Eps)[j] + 
						  REAL(Pi)[i]);
	    }
	    REAL(dp)[i] = REAL(dp)[i] - REAL(Yik)[i+1*I];
	}
	su = 0.0;
	for (i=1; i<I; i++) {
	    su += R_pow(REAL(dp)[i],2.0);
	}
	/* Si on recherche la solution du maximum de vraisemblance (nu = 0), alors le gradient
	   est nul et su, la norme de ce vecteur est nulle aussi. Dans ce cas, on
	   divise 0.0 par 0, et ça fait NaN. Donc dans ce cas de figure, on fixe la norme à 1
	   et pis voila
	*/
	if (su < 1e-10) {
	    su = 1;
	} else {
	    su = sqrt(su/((double) (I-1)));
	}
	for (i=1; i<I; i++) {
	    REAL(Pi)[i] = REAL(Pi)[i] - h * REAL(dp)[i]/su;
	}
	REAL(Pi)[0] = 0;
	/* Actualisation des Alp */
	PROTECT(d3=calcSumd3(nuijl, I, J, Alp));
	for (i=0; i<I; i++) {
	    for (j=0; j<J; j++) {
		REAL(d)[i+j*I] = ((double) nT)*REAL(Sj)[j]*exp(REAL(Alp)[i+j*I])*
		    (REAL(Vj)[j]+exp(REAL(Eps)[j])*exp(REAL(Pi)[i])) - 
		    REAL(Yij)[i+j*I] + 2.0*nu*REAL(d3)[i+j*I];
	    }
	}
	PROTECT(dt = transposeok(d, I, J));
	PROTECT(dt2 = scaleok(dt, J, I));
	PROTECT(dd = transposeok(dt2, J, I));
	for (i=0; i<I; i++) {
	    for (j=0; j<J; j++) {
		REAL(Alp)[i+j*I] = REAL(Alp)[i+j*I] - h * (REAL(dd)[i+j*I])*sqrt(((double) (J-1)));
	    }	
	}

	
	vrais1 = calcVraisPen(nu, nuijl, nT, Y, Ykj, 
			      Yik, Yij, I, J, 2, Eps,
			      Pi, Alp, Sj, lVj);


	if (verb)
	    Rprintf("%f\n", vrais1);

	UNPROTECT(4);
	if (((vrais-vrais1) >= 0)&&((vrais-vrais1)/fabs(vrais) < stopCrit)) {
	    break;
	} 
	if (h < 1e-12)
	    break;
	
	if (vrais-vrais1 <0) {
	    h = h/2.0;
	}
	if (h < 1e-12)
	    break;

    }
    
    if (verb) {
	if (a == maxIter) {
	    Rprintf("Failed to converge. Returning last value\n");
	} else {
	    Rprintf("Converged.\n");
	}
    }    
    REAL(resuvrais)[0] = vrais;
    SET_VECTOR_ELT(resu, 0, resuvrais);
    SET_VECTOR_ELT(resu, 1, Alp);
    SET_VECTOR_ELT(resu, 2, Eps);
    SET_VECTOR_ELT(resu, 3, Pi);
    REAL(resupen)[0] = calcSum3(nuijl, I, J, Alp);
    SET_VECTOR_ELT(resu, 4, resupen);
    REAL(resuvraispeninit)[0] = vraisPenInit;
    SET_VECTOR_ELT(resu, 5, resuvraispeninit);
    INTEGER(NIter)[0] = a;
    REAL(resuh)[0] = h;
    SET_VECTOR_ELT(resu, 6, NIter);
    SET_VECTOR_ELT(resu, 7, resuh);
    
    UNPROTECT(14);

    return(resu);
}





double vraispenpropt(int n, double *par, void *ex)
{
    SEXP Eps, Pi, Alp;
    double res;
    int i, j, l;
    Pouroptim *ex2 = (Pouroptim *) ex;
    
    PROTECT(Alp = allocVector(REALSXP, (ex2->I)*(ex2->J)));
    PROTECT(Eps = allocVector(REALSXP, ex2->J));
    PROTECT(Pi = allocVector(REALSXP, ex2->I));
    
    l = 0;
    for (i = 0; i < ex2->I; i++) {
	for (j = 0; j < ex2->J; j++) {
	    REAL(Alp)[i+j*ex2->I] = par[l];
	    l++;
	}
    }
    for (j = 0; j < ex2->J; j++) {
	REAL(Eps)[j] = par[l];
	l++;
    }
    REAL(Pi)[0] = 0;
    for (i = 1; i < ex2->I; i++) {
	REAL(Pi)[i] = par[l];
	l++;
    }
    
    res = calcVraisPen(ex2->nu, ex2->nuijl, ex2->nT, ex2->Y, ex2->Ykj, 
		       ex2->Yik, ex2->Yij, ex2->I, ex2->J, ex2->K, Eps,
		       Pi, Alp, ex2->Sj, ex2->lVj);
    UNPROTECT(3);
    return(res);
}



void calculgradient(int n, double *par, double *gr, void *ex)
{
    SEXP Eps, Pi, Alp, d3, de, dp;
    int i, j, l;
    Pouroptim *ex2 = (Pouroptim *) ex;
    
    PROTECT(Alp = allocVector(REALSXP, (ex2->I)*(ex2->J)));
    PROTECT(Eps = allocVector(REALSXP, ex2->J));
    PROTECT(de = allocVector(REALSXP, ex2->J));
    PROTECT(Pi = allocVector(REALSXP, ex2->I));
    PROTECT(dp = allocVector(REALSXP, ex2->I));

    l = 0;
    for (i = 0; i < ex2->I; i++) {
	for (j = 0; j < ex2->J; j++) {
	    REAL(Alp)[i+j*ex2->I] = par[l];
	    l++;
	}
    }
    for (j = 0; j < ex2->J; j++) {
	REAL(Eps)[j] = par[l];
	l++;
    }
    REAL(Pi)[0] = 0;
    for (i = 1; i < ex2->I; i++) {
	REAL(Pi)[i] = par[l];
	l++;
    }
    
    /* Calcul du gradient */
    l = 0;
    PROTECT(d3=calcSumd3(ex2->nuijl, ex2->I, ex2->J, Alp));
    for (i=0; i < ex2->I; i++) {
	for (j=0; j< ex2->J; j++) {
	    gr[l] = ((double) ex2->nT)* REAL(ex2->Sj)[j]*exp(REAL(Alp)[i+j*(ex2->I)])*
		(REAL(ex2->Vj)[j]+exp(REAL(Eps)[j])*exp(REAL(Pi)[i])) - 
		REAL(ex2->Yij)[i+j*(ex2->I)] + 2.0*ex2->nu*REAL(d3)[i+j*(ex2->I)];
	    l++;
	}
    }
    
    for (j=0; j<ex2->J; j++) {
	REAL(de)[j] = 0.0;
    }	
    for (j=0; j<ex2->J; j++) {
	for (i=0; i< ex2->I; i++) {
	    REAL(de)[j] += ((double) ex2->nT)*REAL(ex2->Sj)[j]*exp(REAL(Alp)[i+j*ex2->I] + REAL(Eps)[j] + REAL(Pi)[i]);
	}
	REAL(de)[j] = REAL(de)[j] - REAL(ex2->Ykj)[1+j*2];
    }
    for (j = 0; j < ex2->J; j++) {
	gr[l] = REAL(de)[j];
	l++;
    }
    
    for (i=0; i< ex2->I; i++) {
	REAL(dp)[i] = 0.0;
    }
    for (i=0; i< ex2->I; i++) {
	for (j=0; j< ex2->J; j++) {
	    REAL(dp)[i] += ((double) ex2->nT)*REAL(ex2->Sj)[j]*exp(REAL(Alp)[i+j*(ex2->I)] + REAL(Eps)[j] + 
							REAL(Pi)[i]);
	}
	REAL(dp)[i] = REAL(dp)[i] - REAL(ex2->Yik)[i+1*(ex2->I)];
    }
    
    for (i = 1; i < ex2->I; i++) {
	gr[l] = REAL(dp)[i];
	l++;
    }
    
    UNPROTECT(6);
}


SEXP EstimMaxVraisPenl2alt2(SEXP AlpInit, SEXP PiInit, SEXP EpsInit,
			    int nT, int I, int J, SEXP Sj, SEXP lVj, 
			    SEXP Vj, SEXP Y,
			    SEXP Ykj, SEXP Yik, SEXP Yij, 
			    int maxIter, double stopCrit, double nu,
			    SEXP nuijl, double h, int verb)
{
    
    SEXP Alp, Pi, Eps, theta, resuvrais, resuvraispeninit, resu, resuh, resupen, NIter;
    double *thetar, vrais1, vraisPenInit, vrais;
    int nREPORT;
    double abstol = log(0.0);
    double reltol = stopCrit;
    int fncount, grcount, convergence;
    int *mask, i, l, j;
    Pouroptim ex;

    
    PROTECT(Alp = allocVector(REALSXP, I*J));
    PROTECT(Pi = allocVector(REALSXP, I));
    PROTECT(Eps = allocVector(REALSXP, J));
    PROTECT(theta = allocVector(REALSXP, I+J+(I*J)-1));
    PROTECT(resuvrais = allocVector(REALSXP, 1));
    PROTECT(resuvraispeninit = allocVector(REALSXP, 1));
    PROTECT(resupen = allocVector(REALSXP, 1));
    PROTECT(resuh = allocVector(INTSXP, 1));
    PROTECT(NIter = allocVector(INTSXP, 1));
    PROTECT(resu = allocVector(VECSXP, 8));

    thetar = REAL(theta);
    nREPORT = 1;
    
    ex.I = I;
    ex.J = J;
    ex.K = 2;
    ex.nu = nu;
    ex.nuijl = nuijl;
    ex.nT = nT;
    ex.Y = Y;
    ex.Ykj = Ykj;
    ex.Yik = Yik;
    ex.Yij = Yij;
    ex.Sj = Sj;
    ex.lVj = lVj;
    ex.Vj = Vj;

    /* Copie des valeurs initiales */
    l = 0;
    for (i = 0; i < I; i++) {
	for (j = 0; j < J; j++) {
	    thetar[l] = REAL(AlpInit)[i+j*I];
	    l++;
	}
    }
    for (j = 0; j < J; j++) {
	thetar[l] = REAL(EpsInit)[j];
	l++;
    }
    for (i = 1; i < I; i++) {
	thetar[l] = REAL(PiInit)[i];
	l++;
    }
    fncount = 0;
    grcount = 0;
    convergence = 0;

    /* Vraisemblance initiale */
    vrais1 = calcVraisPen(nu, nuijl, nT, Y, Ykj, 
			  Yik, Yij, I, J, 2, EpsInit,
			  PiInit, AlpInit, Sj, lVj);
    vraisPenInit = vrais1;
    
    mask = (int *) R_alloc(length(theta), sizeof(int));
    for (i = 0; i < length(theta); i++) mask[i] = 1;

    vmmin(length(theta), thetar, &vrais1,
	  vraispenpropt, calculgradient, maxIter, verb,
	  mask, abstol, reltol, nREPORT,
	  &ex, &fncount, &grcount, &convergence);
    /* Les arguments sont les suivants:
       - nombre de paramètres
       - Valeurs initiales des paramètres
       - Valeur de la vraisemblance qui devra être renvoyée pour le résultat
       - fonction permettant le calcul de la vraisemblance
       - fonction permettant le calcul du gradient
       - nombre maximum d'itérations
       - verbose?
       - un vecteur indiquant quels sont les éléments de theta à modifier dans
       cette minimisation (mask=1) ou les éléments à ne pas prendre en compte (mask=0) 
       - tolérance absolue (indique que la foction ne doit pas être inférieure à cette valeur)
       - tolérance relative (utilisée pour identifier la convergence)
       - afficher un rapport tous les combien d'itérations?
       - les données utilisées par les fonctions de vraisemblance et de gradient
       - le nombre de fois que optim fait appel à la fonction de vraisemblance
       - le nombre de fois que optim fait appel à la fonction de gradient
       - la convergence
    */
       
    l = 0;
    for (i = 0; i < I; i++) {
	for (j = 0; j < J; j++) {
	    REAL(Alp)[i+j*I] = thetar[l];
	    l++;
	}
    }
    for (j = 0; j < J; j++) {
	REAL(Eps)[j] = thetar[l];
	l++;
    }
    REAL(Pi)[0] = 0;
    for (i = 1; i < I; i++) {
	REAL(Pi)[i] = thetar[l];
	l++;
    }

    vrais = calcVraisPen(nu, nuijl, nT, Y, Ykj, 
			 Yik, Yij, I, J, 2, Eps,
			 Pi, Alp, Sj, lVj);
    
    REAL(resuvrais)[0] = vrais;
    SET_VECTOR_ELT(resu, 0, resuvrais);
    SET_VECTOR_ELT(resu, 1, Alp);
    SET_VECTOR_ELT(resu, 2, Eps);
    SET_VECTOR_ELT(resu, 3, Pi);
    REAL(resupen)[0] = calcSum3(nuijl, I, J, Alp);
    SET_VECTOR_ELT(resu, 4, resupen);
    REAL(resuvraispeninit)[0] = vraisPenInit;
    SET_VECTOR_ELT(resu, 5, resuvraispeninit);
    INTEGER(NIter)[0] = fncount;
    INTEGER(resuh)[0] = (double) convergence;
    SET_VECTOR_ELT(resu, 6, NIter);
    SET_VECTOR_ELT(resu, 7, resuh);
    UNPROTECT(10);
    

    return(resu);
}







SEXP ModeliseBien(SEXP Esp, SEXP Dep, SEXP Sta, SEXP Y, SEXP Vjt, 
		  SEXP Ir, SEXP Jr, SEXP Kr, SEXP Sjt, SEXP nTr,
		  SEXP maxIter, SEXP stopCrit, SEXP nuijl, SEXP nur, 
		  SEXP hr, SEXP typealgo, SEXP verbo, SEXP StartVal)
{
    SEXP Vj, lVj, Sj, loglambdajki, resu, AlpInit, EpsInit, PiInit;
    SEXP Ykj, Yik, Yij, resufin, resufin2;
    int i,j,k, I, J, K, nT, verb, deprot;
    double sum2, sum1, VraisCompl, nu, h;
    
    /* Définition des constantes */
    I = INTEGER(Ir)[0];
    J = INTEGER(Jr)[0];
    K = INTEGER(Kr)[0];
    nT = INTEGER(nTr)[0];
    verb = INTEGER(verbo)[0];
    nu = REAL(nur)[0];
    h = REAL(hr)[0];
    
    /* Mise en forme des Vj et Sj */
    PROTECT(Vj = allocVector(REALSXP, J));
    PROTECT(Sj = allocVector(REALSXP, J));
    PROTECT(lVj = allocVector(REALSXP, J));    
    for (j = 0; j < J; j++) {
	REAL(Vj)[j] = REAL(Vjt)[j];
	REAL(Sj)[j] = REAL(Sjt)[j];
	REAL(lVj)[j] = log(REAL(Vjt)[j]);
    }
    
    /* Calcul des statistiques suffisantes */
    PROTECT(Ykj = allocVector(REALSXP, K*J));
    PROTECT(Yik = allocVector(REALSXP, I*K));
    PROTECT(Yij = allocVector(REALSXP, I*J));    
    /* Initialisation à zéro */
    for (j=0; j<J; j++){
	for (k=0; k<K; k++){
	    REAL(Ykj)[k+j*K] = 0;
	}
    }
    for (i=0; i<I; i++){
	for (k=0; k<K; k++){
	    REAL(Yik)[i+k*I] = 0;
	}
    }
    for (j=0; j<J; j++){
	for (i=0; i<I; i++){
	    REAL(Yij)[i+j*I] = 0;
	}
    }
    /* Calcul des totaux */
    for (i=0; i<I; i++){
	for (j=0; j<J; j++){
	    for (k=0; k<K; k++){
		REAL(Ykj)[k+j*K] += REAL(Y)[j+k*J+i*K*J];
		REAL(Yik)[i+k*I] += REAL(Y)[j+k*J+i*K*J];
		REAL(Yij)[i+j*I] += REAL(Y)[j+k*J+i*K*J];
	    }
	}
    }
    
    if ((length(StartVal) == 1)||(REAL(nur)[0]< 1e-7)) {
	if (INTEGER(verbo)[0])
	    Rprintf("Starting values not provided, computing Maximum likelihood estimates\n");
	
	/* calcul de la vraisemblance pour le modèle complet */
	PROTECT(loglambdajki = allocVector(REALSXP, J*I*K));
	for (i = 0; i < (I*J*K); i++) {
	    if (REAL(Y)[i]> 0.5) {
		REAL(loglambdajki)[i] = log(REAL(Y)[i]);
	    } else {
		REAL(loglambdajki)[i] = -30.0;
	    }
	}
	
	sum2 = 0.0;
	for (i = 0; i < (I*J*K); i++) {
	    sum2 += REAL(loglambdajki)[i] * REAL(Y)[i];
	}
	sum1 = 0.0;
	for (j = 0; j < J; j++) {
	    for (k = 0; k < K; k++) {
		for (i = 0; i < I; i++) {
		    sum1 += REAL(Y)[j+k*J+i*K*J]*REAL(Sj)[j];
		}
	    }
	}
	sum1 = sum1 * ((double) nT);
	VraisCompl = sum1-sum2;
	
	/* Valeurs initiales des paramètres */
	
	PROTECT(AlpInit = allocVector(REALSXP, J*I));
	PROTECT(EpsInit = allocVector(REALSXP, J));
	PROTECT(PiInit = allocVector(REALSXP, I));
	for (i = 0; i < I; i++) {
	    for (j = 0; j < J; j++) {
		REAL(AlpInit)[i+j*I] = REAL(loglambdajki)[j+K*J*i] - REAL(lVj)[j];
	    }
	}
	for (j = 0; j < J; j++) {
	    REAL(EpsInit)[j] = REAL(loglambdajki)[j+J] - REAL(AlpInit)[j*I];
	}
	REAL(PiInit)[0] = 0.0;
	for (i = 1; i < I; i++) {
	    REAL(PiInit)[i] = 0.0;
	    for (j = 0; j < J; j++) {
		REAL(PiInit)[i] += REAL(loglambdajki)[j+J+K*J*i] - 
		    REAL(AlpInit)[i+j*I] - REAL(EpsInit)[j];
	    }
	    REAL(PiInit)[i] = REAL(PiInit)[i]/((double) J);
	}
		
	PROTECT(resu = EstimMaxVrais(AlpInit, PiInit, EpsInit,
				     nT, I, J, Sj, lVj, 
				     Vj, Y,
				     Ykj, Yik, Yij, 
				     INTEGER(maxIter)[0], REAL(stopCrit)[0],
				     VraisCompl));
    } else {
	resu = StartVal;
    }
    
    switch (INTEGER(typealgo)[0]) {
    case 1 :
	PROTECT(resufin = EstimMaxVraisPenl2alt(VECTOR_ELT(resu, 1), 
						VECTOR_ELT(resu, 3), VECTOR_ELT(resu, 2),
						nT, I, J, Sj, lVj, Vj, Y, Ykj, Yik, Yij, 
						INTEGER(maxIter)[0], REAL(stopCrit)[0], 
						nu, nuijl,  h, verb));
	break;
    case 2 : 
	PROTECT(resufin = EstimMaxVraisPenl2(VECTOR_ELT(resu, 1), 
					     VECTOR_ELT(resu, 3), VECTOR_ELT(resu, 2),
					     nT, I, J, Sj, lVj, Vj, Y, Ykj, Yik, Yij, 
					     INTEGER(maxIter)[0], REAL(stopCrit)[0], 
					     nu, nuijl,  h, verb));
	break;
    case 3 :
	PROTECT(resufin = EstimMaxVraisPenl2alt2(VECTOR_ELT(resu, 1), 
						 VECTOR_ELT(resu, 3), VECTOR_ELT(resu, 2),
						 nT, I, J, Sj, lVj, Vj, Y, Ykj, Yik, Yij, 
						 INTEGER(maxIter)[0], REAL(stopCrit)[0], 
						 nu, nuijl,  h, verb));
	break;
    default :
	PROTECT(resufin2 = EstimMaxVraisPenl2(VECTOR_ELT(resu, 1), 
					      VECTOR_ELT(resu, 3), VECTOR_ELT(resu, 2),
					      nT, I, J, Sj, lVj, Vj, Y, Ykj, Yik, Yij, 
					      INTEGER(maxIter)[0], 1e-6, 
					      nu, nuijl,  h, verb));
	PROTECT(resufin = EstimMaxVraisPenl2alt2(VECTOR_ELT(resufin2, 1), 
						 VECTOR_ELT(resufin2, 3), VECTOR_ELT(resufin2, 2),
						 nT, I, J, Sj, lVj, Vj, Y, Ykj, Yik, Yij, 
						 INTEGER(maxIter)[0], REAL(stopCrit)[0], 
						 nu, nuijl,  h, verb));
	
    }

    deprot = 7;
    if ((length(StartVal) == 1)||(REAL(nur)[0]< 1e-7)) {
	deprot+=5;
    }
    if (INTEGER(typealgo)[0]==4) 
	deprot++;
    
    UNPROTECT(deprot);
    
    return(resufin);

}




SEXP PreditBien(SEXP res, SEXP Sj, SEXP Vj, SEXP Yv, int I, int J, int K, double ng)
{
    SEXP aij, ejk, pik, ResultatsQ;
    int i, j, k, effe;
    double eff, pi, Q1, Q2, Q3, Q4, pred, obs;
    
    PROTECT(ResultatsQ = allocVector(REALSXP, 4));

    /* On récupère les coefficients */
    aij = VECTOR_ELT(res, 1);
    ejk = VECTOR_ELT(res, 2);
    pik = VECTOR_ELT(res, 3);
    

    /* On effectue la prédiction */
    Q1 = 0.0;
    Q2 = 0.0;
    Q3 = 0.0;
    Q4 = 0.0;
//    Rprintf("toto <- c(, ");
    effe=0;
    for (i = 0; i < I; i++) {
	for (k = 0; k < K; k++) {
	    for (j = 0; j < J; j++) {
		if (k == 0) {
		    eff = log(REAL(Vj)[j]);
		    pi = 0.0;
		} else {
		    eff = REAL(ejk)[j];
		    pi = REAL(pik)[i];
		}
		pred = REAL(aij)[i + j*I]  + eff + pi;
		if (pred < -30)
		    pred = -30;
		pred = REAL(Sj)[j+k*J+i*J*K] * exp(pred)/ng;
		obs = REAL(Yv)[j+k*J+i*J*K];
		if (pred > exp(-29.9999999999)) {
		    Q1 += -dpois(obs, pred, 1);
		    Q3 += R_pow(pred - obs, 2.0)/pred;
		    effe++;
		}		
		Q2 += R_pow(pred - obs, 2.0);
		Q4 += R_pow(pred - obs, 2.0)/(pred+1.0);
	    }
	}
    }
//    Rprintf("\n\n\n\n");
    
    
    REAL(ResultatsQ)[0] = Q1/((double) effe);
    REAL(ResultatsQ)[1] = Q2/((double) (I*J*K));
    REAL(ResultatsQ)[2] = Q3/((double) effe);
    REAL(ResultatsQ)[3] = Q4/((double) (I*J*K));

    UNPROTECT(1);

    /* resultat */
    return(ResultatsQ);
}




SEXP PreditBienParEsp(SEXP res, SEXP Sj, SEXP Vj, SEXP Yv, int I, int J, int K, double ng)
{
    SEXP aij, ejk, pik, ResultatsQ, Q1, Q2, Q3, Q4;
    int i, j, k;
    double eff, pi, pred, obs;
    
    PROTECT(ResultatsQ = allocVector(VECSXP, 4));
    PROTECT(Q1 = allocVector(REALSXP, I));
    PROTECT(Q2 = allocVector(REALSXP, I));
    PROTECT(Q3 = allocVector(REALSXP, I));
    PROTECT(Q4 = allocVector(REALSXP, I));
    
    for (i = 0; i < I; i++) {
	REAL(Q1)[i] = 0.0;
	REAL(Q2)[i] = 0.0;
	REAL(Q3)[i] = 0.0;
	REAL(Q4)[i] = 0.0;
    }

    /* On récupère les coefficients */
    aij = VECTOR_ELT(res, 1);
    ejk = VECTOR_ELT(res, 2);
    pik = VECTOR_ELT(res, 3);
    

    /* On effectue la prédiction */
    for (i = 0; i < I; i++) {
	for (k = 0; k < K; k++) {
	    for (j = 0; j < J; j++) {
		if (k == 0) {
		    eff = log(REAL(Vj)[j]);
		    pi = 0.0;
		} else {
		    eff = REAL(ejk)[j];
		    pi = REAL(pik)[i];
		}
		pred = REAL(aij)[i + j*I]  + eff + pi;
		if (pred < -30)
		    pred = -30;
		pred = REAL(Sj)[j+k*J+i*J*K] * exp(pred)/ng;
		obs = REAL(Yv)[j+k*J+i*J*K];
		REAL(Q1)[i] += -dpois(obs, pred, 1);
		REAL(Q3)[i] += R_pow(pred - obs, 2.0)/pred;
		REAL(Q2)[i] += R_pow(pred - obs, 2.0);
		REAL(Q4)[i] += R_pow(pred - obs, 2.0)/(pred+1.0);
	    }
	}
    }    
    
    for (i = 0; i < I; i++) {
	REAL(Q1)[i] = REAL(Q1)[i];
	REAL(Q2)[i] = REAL(Q2)[i];
	REAL(Q3)[i] = REAL(Q3)[i];
	REAL(Q4)[i] = REAL(Q4)[i];
    }

    SET_VECTOR_ELT(ResultatsQ, 0, Q1);
    SET_VECTOR_ELT(ResultatsQ, 1, Q2);
    SET_VECTOR_ELT(ResultatsQ, 2, Q3);
    SET_VECTOR_ELT(ResultatsQ, 3, Q4);

    UNPROTECT(5);

    /* resultat */
    return(ResultatsQ);
}






SEXP ValidationCroisee(SEXP obsEsp, SEXP obsDep, SEXP obsSta, SEXP obsGroup, SEXP ngr,
		       SEXP Vj, SEXP Ir, SEXP Jr, SEXP Kr, SEXP nTr, SEXP depArea,
		       SEXP maxIter, SEXP stopCrit, SEXP nuijl, SEXP nur, 
		       SEXP hr, SEXP typealgo, SEXP verbose, SEXP StartVal, SEXP tracer)
{
    int ngroup, g, i, j, k, l, n, I, J, K, trace;
    SEXP Esp, Dep, Sta, Y, Sjt, Yv, res, res2, critere, Vjt, gg;
    
    trace = INTEGER(tracer)[0];
    n = length(obsEsp);
    I = INTEGER(Ir)[0];
    J = INTEGER(Jr)[0];
    K = INTEGER(Kr)[0];
    ngroup = INTEGER(ngr)[0];
    
    /* allocation de mémoire */
    PROTECT(Esp=allocVector(INTSXP, I*J*K));
    PROTECT(Dep=allocVector(INTSXP, I*J*K));
    PROTECT(Sta=allocVector(INTSXP, I*J*K));
    PROTECT(Sjt=allocVector(REALSXP, I*J*K));
    PROTECT(Y=allocVector(REALSXP, I*J*K));
    PROTECT(Yv=allocVector(REALSXP, I*J*K));
    PROTECT(res2=allocVector(VECSXP, ngroup));    
    PROTECT(Vjt=allocVector(REALSXP, I*J*K));    
    
    /* Construction des vecteurs constants */
    l=0;
    for (i = 1; i <= I; i++) {
	for (k = 1; k <= K; k++) {
	    for (j = 1; j <= J; j++) {
		INTEGER(Esp)[l] = i;
		INTEGER(Dep)[l] = j;
		INTEGER(Sta)[l] = k;
		REAL(Sjt)[l] = REAL(depArea)[j-1];
		REAL(Vjt)[l] = REAL(Vj)[j-1];
		l++;
	    }
	}
    }
    

    PROTECT(gg = StartVal);

    /* et pour chaque groupe */
    for (g = 1; g <= ngroup; g++) {
	
	
	if (trace)
	    Rprintf("#### Group: %i\n", g);

	/* Initialisation de Y et Yv à 0 */
	for (l = 0; l < (I*J*K); l++) {
	    REAL(Y)[l] = 0.0;
	    REAL(Yv)[l] = 0.0;
	}
	
	/* on compte le nombre de détections de chaque
	   espèce, département, statut */
	for (l = 0; l<n; l++) {
	    if (INTEGER(obsGroup)[l]!=g) {
		REAL(Y)[(INTEGER(obsDep)[l]-1)+(INTEGER(obsSta)[l]-1)*J+
			(INTEGER(obsEsp)[l]-1)*K*J] += 1.0;
	    } else {
		REAL(Yv)[(INTEGER(obsDep)[l]-1)+(INTEGER(obsSta)[l]-1)*J+
			 (INTEGER(obsEsp)[l]-1)*K*J] += 1.0;
	    }
	} 
	
	R_CheckUserInterrupt();
		
	/* Ajustement du modèle */
	PROTECT(res = ModeliseBien(Esp, Dep, Sta, Y, Vjt, 
				   Ir, Jr, Kr, Sjt, nTr,
				   maxIter, stopCrit, nuijl, nur, 
				   hr, typealgo, verbose, gg));
	

	/* Vérif que convergence */
	if (INTEGER(VECTOR_ELT(res, 6))[0]==(INTEGER(maxIter)[0]))
	    Rprintf("Non convergence pour le modèle %i\n", g);
	
	/* Prediction: */
	PROTECT(critere = PreditBienParEsp(res, Sjt, Vj, Yv, I, J, K, (double) ngroup-1));
	SET_VECTOR_ELT(res2, g-1, critere);
	gg = res;
	UNPROTECT(2);
    }
    


    UNPROTECT(9);
    return(res2);
}






SEXP BootstrapModele(SEXP obsEsp, SEXP obsPra, SEXP obsSta, SEXP Vj, SEXP Ir, SEXP Jr, SEXP Kr, 
		     SEXP nTr, SEXP nur, SEXP nbootstrap, SEXP nuijl, SEXP maxIter, SEXP deparea,
		     SEXP stopCrit, SEXP hr, SEXP typealgo, SEXP verbose, SEXP tracer,
		     SEXP StartVal)
{
    int i, j, k, l, n, I, J, K;
    SEXP Yb, sam, Esp, Pra, Sta, Sjt, Vjt, res;
    
    n = length(obsEsp);
    I = INTEGER(Ir)[0];
    J = INTEGER(Jr)[0];
    K = INTEGER(Kr)[0];
    
    /* allocation de mémoire */
    PROTECT(sam=allocVector(INTSXP, n));
    PROTECT(Yb=allocVector(REALSXP, I*J*K));
    PROTECT(Esp=allocVector(INTSXP, I*J*K));
    PROTECT(Pra=allocVector(INTSXP, I*J*K));
    PROTECT(Sta=allocVector(INTSXP, I*J*K));
    PROTECT(Sjt=allocVector(REALSXP, I*J*K));
    PROTECT(Vjt=allocVector(REALSXP, I*J*K));
    
    /* Construction des vecteurs constants */
    l=0;
    for (i = 1; i <= I; i++) {
	for (k = 1; k <= K; k++) {
	    for (j = 1; j <= J; j++) {
		INTEGER(Esp)[l] = i;
		INTEGER(Pra)[l] = j;
		INTEGER(Sta)[l] = k;
		REAL(Sjt)[l] = REAL(deparea)[j-1];
		REAL(Vjt)[l] = REAL(Vj)[j-1];
		l++;
	    }
	}
    }
    
    
    /* Bootstrap */
    R_CheckUserInterrupt();
    
    /* Échantillonnage */
    GetRNGstate();
    for (i = 0; i<n; i++) {
	INTEGER(sam)[i] = ((int) floor(((double) n) * unif_rand()));
    }
    PutRNGstate();
    
    /* Initialisation de Y à 0 */
    for (l = 0; l < (I*J*K); l++) {
	REAL(Yb)[l] = 0.0;
    }
    
    /* on compte le nombre de détections de chaque
       espèce, PRA, statut */
    for (l = 0; l<n; l++) {
	REAL(Yb)[(INTEGER(obsPra)[INTEGER(sam)[l]]-1)+
		 (INTEGER(obsSta)[INTEGER(sam)[l]]-1)*J+
		 (INTEGER(obsEsp)[INTEGER(sam)[l]]-1)*K*J] += 1.0;
    }
//    Rprintf("%i\n", J);
    /* Pour chaque modèle */
    PROTECT(res = ModeliseBien(Esp, Pra, Sta, Yb, Vjt, 
			       Ir, Jr, Kr, Sjt, nTr,
			       maxIter, stopCrit, nuijl, nur, 
			       hr, typealgo, verbose, StartVal));
    /* Vérif que convergence */
    if (INTEGER(VECTOR_ELT(res, 6))[0]==(INTEGER(maxIter)[0]))
	Rprintf("Non convergence\n");    

    UNPROTECT(8);
    return(res);
}

SEXP Calculeffort(SEXP liepv, SEXP linumpra, SEXP inclutVeh, SEXP praArea, SEXP Jr)
{
    int j, l, m, J, nv;
    SEXP Vj, praa, prab;
    
    J = INTEGER(Jr)[0];
    nv = length(liepv);
    
    PROTECT(Vj = allocVector(REALSXP, J));

    /* initialisation de Vj à 0 */
    for (j = 0; j < J; j++) {
	REAL(Vj)[j] = 0;
    }

    for (l = 0; l<nv; l++) {
	if (INTEGER(inclutVeh)[l]) {
	    praa = VECTOR_ELT(liepv, l);
	    prab = VECTOR_ELT(linumpra, l);
	    for (m = 0; m < length(praa); m++) {
		REAL(Vj)[INTEGER(prab)[m]] += REAL(praa)[m];
	    }
	}
    }

    for (j = 0; j < J; j++) {
	REAL(Vj)[j] = REAL(Vj)[j]/REAL(praArea)[j];
    }
    

    UNPROTECT(1);
    return(Vj);
}
