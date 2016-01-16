#include "..\headers\my.h"

int Nit=0, Nz = _Nz, Nq = _Nq, Nd = _Nd, Nj = _Nj, Na = _Na, Nt, Nap, Ns, is_linearization, is_b_min_inheritance, is_b_split_inheritance;
double z_start[_Nz], zmin[_Nz], zmax[_Nz], **b_start, bmin, bmax, d[_Nd], dmin[_Nd], dmax[_Nd], **E, **A, *P, Rabs[_Nq], alpha[_Na];
double start_q[_Nq], start_minq[_Nq], start_maxq[_Nq];
double gamma;

//FILE *log;

//double z_start[_Nz] = {0.1};
//double zmin[_Nz]	= {0};
//double zmax[_Nz]	= {1};
//double bmin	= -1e01;
//double bmax	= 1e01;

struct T *ROPUD_pT = NULL;		//	Для передачи указателя обрабатываемой 
struct TK *ROPUD_pTK = NULL;	//	области и номера ограничения для работы 
//int ROPUD_ConstraintToMaximize;	//	upperBound, lowerBound, maximize


int read_z(double** z)	//	Создает массив стартовых z 
{
	int iz;
	*z = (double*)malloc(Nz*sizeof(double));
	if(!(*z)) return 0;
	for(iz=0; iz<Nz; iz++)
		(*z)[iz] = z_start[iz];	//	Начальные значения управляющих параметров
	return 1;
}

int set_z(struct TK *pT, double* q)	//	Пересчитывает управления z, для апроксимаций b и точек q
{
	int iz,iq;
	for(iz=0; iz<Nz; iz++)
	{
		pT->z[iz] = pT->b[iz][Nq];	//	Свободный член
		for(iq=0; iq<Nq; iq++)
			pT->z[iz] += pT->b[iz][iq]*q[iq];
	}
}

int reset_b(double*** b)	//	Создает массив стартовых b на основе глобальных z_start
						//	Работает с абсолютными значениями Q
{
	int iz,iq;
	double S = 1;	//	1 + сумма стартовых средних значений тетта
	//for(iq=0; iq<Nq; iq++)
	S = 25;
	*b = malloc(Nz*sizeof(double));
	if(!(*b)) return 0;
	for(iz=0; iz<Nz; iz++)
	{
		(*b)[iz] = malloc((Nq+1)*sizeof(double));
		if(!(*b)[iz]) return 0;
		for(iq=0; iq<Nq+1; iq++) 
			(*b)[iz][iq] = 0;z_start[iz] /S;
	}
	(*b)[0][Nq] = z_start[0];
	(*b)[1][Nq] = z_start[1];
	return 1;
}

void replant(struct BinTree *BT, struct TK *pTK)
{
	if(!BT->region)
	{
		replant(BT->lBranch, pTK);
		replant(BT->rBranch, pTK);
	}
	else 
		BT->region->relK = pTK;
}

void branchOut(struct BinTree *BT, struct T *pT1, struct T *pT2)
{
	if(!BT->region)
	{
		branchOut(BT->lBranch, pT1, pT2);
		branchOut(BT->rBranch, pT1, pT2);
	}
	else if(pT1 == BT->region)
	{
		BT->lBranch = malloc(sizeof(struct BinTree));
		BT->rBranch = malloc(sizeof(struct BinTree));
		BT->lBranch->region = pT1;
		BT->rBranch->region = pT2;
		BT->lBranch->lBranch = BT->lBranch->rBranch 
			= BT->rBranch->lBranch = BT->rBranch->rBranch
								= BT->region = NULL;
	}
}

double Fx (double xxi)
{
	double result, xx;
//----------------------------------------------------
	result=1;
	if (xxi<=-4) result=0;
	if (xxi<4 && xxi>-4)
		{
			result=5.000001928379E-1;
			xx=xxi;
			result+=3.987112373793E-1*xx;
			xx*=xxi;
			result+=6.551375743198E-9*xx;
			xx*=xxi;
			result+=-6.586542282821E-2*xx;
			xx*=xxi;
			result+=-6.36802163649E-9*xx;
			xx*=xxi;
			result+=9.468042848461E-3*xx;
			xx*=xxi;
			result+=2.390385414721E-9*xx;
			xx*=xxi;
			result+=-9.955656545868E-4*xx;
			xx*=xxi;
			result+=-4.331928868847E-10*xx;
			xx*=xxi;
			result+=7.399573280813E-5*xx;
			xx*=xxi;
			result+=4.195944540559E-11*xx;
			xx*=xxi;
			result+=-3.733098012039E-6*xx;
			xx*=xxi;
			result+=-2.222323537432E-12*xx;
			xx*=xxi;
			result+=1.203330765428E-7*xx;
			xx*=xxi;
			result+=6.059591583327E-14*xx;
			xx*=xxi;
			result+=-2.219774606587E-9*xx;
			xx*=xxi;
			result+=0*xx;
			xx*=xxi;
			result+=1.77458743861E-11*xx;
		}
	return result;
}

void maximize(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,
     int n,int m,int me,double ko[],double komod[],double teta[])
{
	double *lb, *rb, Q[_Nq];
	double D, Z;
	int i;

	switch(ROPUD_pT->isSoftCons)
	{
	case 0:
		lb = ROPUD_pT->lbound;
		rb = ROPUD_pT->rbound;
		break;
	case 1:
		lb = ROPUD_pT->slbound;
		rb = ROPUD_pT->srbound;
		break;
	}

	//D = *d;

	//	Приведение к абсолютному значению
		//	Все поисковые - тетта. Соответственно, x - массив, 
		//	содержащий относительное значение критической точки
	for(i=0; i<Nq; i++)
		Q[i] = x[i]*(rb[i]-lb[i]) + lb[i];

	set_z(ROPUD_pT->relK, Q);
	calcModel(d, ROPUD_pT->relK->z, Q);
	//Z = *ROPUD_pT->relK->z;

	*pfc = - calcConstraint(ROPUD_pT->NCons);
	//switch(ROPUD_pT->NCons)
	//{
	//		//	Мягкие
	//case 0:
	//	*pfc = -(-D-Z*Q[0]-2*Q[1]+1);	//-D-Z*Q[0]-2*Q[1]+2	
	//	break;
	//case 1:
	//	*pfc = -(-2*D-2*Q[0]-Z*Q[1]+2);	//-2*D-2*Q[0]-4*Z*Q[1]+5	
	//	break;

	//		//	Жесткие
	//case 2:
	//	*pfc = -(zmin[0] - Z);
	//	break;
	//case 3:
	//	*pfc = -(Z - zmax[0]);
	//	break;
	//}
}

void upperBound(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,
     int n,int m,int me,double ko[],double komod[],double teta[])
{
	double D, Z;
	double *lb, *rb, *Qcmp, Q[_Nq], **b;
	double *pE, ETK[_Nq];
	int *piE;
	double *dF_dQ, *MQi;
	int iz, iq;
	struct T *pT;
	struct TK *pTK;
	int indexq, indext, indexz, indexa,/* indexd,*/ icrit, i_komod, total_krit;


	double d[_Nd];
	int *NQcr;
	double U;
	int i,j,k,s,ij;
	double koef;
	double y[_Nq],sig[_Nq],kk;
	int nqcr,nqapcr;

	int N;
	double F;

	double *pF=NULL;


	*pfc=0;

	pF=malloc(Nap*sizeof(double));
	if(pF==NULL) exit(2);

	kk = 3.9;	//	Правило трех сигм

	for(i=0; i<Nq; i++)
		sig[i] = (start_q[i]-start_minq[i])/kk;


	ij = Nd;

	pT = ROPUD_pT;
	while(pT=pT->pnext)
		if(pT->isSoftCons)
		{
			for(iq=0; iq<Nq; iq++,ij++)
			{
				pT->slbound[iq] = x[ij];
				pT->srbound[iq] = x[ij+Nq];
			}
			ij+=Nq;
		};

		//	Копирование значений b
	pTK = ROPUD_pTK;
	while(pTK=pTK->pnext)
		for(iz=0; iz<Nz; iz++)
			for(iq=0; iq<Nq+1; iq++, ij++)
				pTK->b[iz][iq] = x[ij];

//ограничения
//ограничения Эф Круглое
	for(i=0; i<Na; i++)
			f[i] = alpha[i];

	pT = ROPUD_pT;
	while(pT=pT->pnext)
	{
		if(pT->isSoftCons)
		{
			for (j=0, F=1; j<Nq; j++)
				F *= (Fx((pT->srbound[j]-start_q[j])/sig[j]) - Fx((pT->slbound[j]-start_q[j])/sig[j]));
			f[pT->NCons] -= F;
		}
	}
			//f[0] *= 2;


//вероятностные
	i = Na;
	pT = ROPUD_pT;
	for(indext=0; indext<Nt; indext++)
	{
		pT = pT->pnext;

			//	Ограничения на области по каждому параметру (правая граница больше левой)
		if(pT->isSoftCons)
			for(j=0; j<Nq; j++,i++)
				f[i] = pT->slbound[j] - pT->srbound[j];

		for(icrit=0; icrit<pT->NQcr; icrit++)
		{
			Qcmp = pT->Qcr + icrit*Nq;
			lb = pT->slbound;
			rb = pT->srbound;

			for(iq=0; iq<Nq; iq++)
				Q[iq] = Qcmp[iq]*(rb[iq]-lb[iq]) + lb[iq];

			set_z(pT->relK, Q);
			calcModel(x, pT->relK->z, Q);

			f[i++] = calcConstraint(pT->NCons);
		}
	}


	j = 0;
	//pTK = ROPUD_pTK->pnext;
	//while(pTK = pTK->pnext)
	{
		//pE = pTK->E;
		//piE = pTK->iE;

		for(iq=0; iq<Nq; iq++)
			Q[iq] = (start_minq[iq] + start_maxq[iq]) /2;
		
		//set_z(pTK, Q);
		//calcModel(x, pTK->z, Q);
		//b = pTK->b;

//**************************************************************
		pF[j] = calcCriteria()*1;

		if(is_linearization)
		{
			dF_dQ = calcDerivative(x, pTK->z, pTK->b, Q);
			for(iq=0; iq<Nq; iq++)
				pF[j] += dF_dQ[iq]*(*pE++-Q[iq]*pTK->ai[iq]);	
		}

		//for(iz=0; iz<Nz*2; iz++)
		//	f[i++] = calcConstraintOnZ(iz);

		j++;
	}

	//if(Nap==1) *pfc = *pF;
	//else 
	{
		for(j=0;j<1;j++)
			*pfc+=pF[j];
	}

		//f[i++] = -*pfc;

	free(pF);
}

void maxCriteria(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,
     int n,int m,int me,double ko[],double komod[],double teta[])
{
	double a,*ai;
	int i;

	double *lb, *rb, Q[_Nq], Ql[_Nq], *dF_dQ, D, Z, **b;

	lb = ROPUD_pTK->lbound;
	rb = ROPUD_pTK->rbound;

	//	Приведение к абсолютному значению
		//	Все поисковые - тетта. Соответственно, x - массив, 
		//	содержащий относительное значение критической точки
	for(i=0; i<Nq; i++)
		Q[i] = x[i]*(rb[i]-lb[i]) + lb[i];

//f
	set_z(ROPUD_pTK, Q);
	calcModel(d, ROPUD_pTK->z, Q);

	*pfc = calcCriteria();

	f[0] = 0;-*pfc;

//f с чертой

	for(i=0; i<Nq; i++)
		Ql[i] = 0.5*(rb[i]-lb[i]) + lb[i];	//	0.5 - середина

	set_z(ROPUD_pTK, Ql);
	calcModel(d, ROPUD_pTK->z, Ql);

	*pfc -= calcCriteria();

	if(is_linearization)
	{
		dF_dQ = calcDerivative(d, ROPUD_pTK->z, ROPUD_pTK->b, Ql);
		for(i=0; i<Nq; i++)
			*pfc -= dF_dQ[i]*(Q[i] - Ql[i]);	
	}

	//*pfc=-pow(*pfc, 2);
	*pfc=-pow(*pfc, 2);
}


void calcDD()
{
	int i;
	double csb[9], csv1[9], csv2[9], csrt1[9], csrt2[9], dcsdv1[9], dcsdrt1[9], dcsdv2[9], dcsdrt2[9], dd = 0.0001;
	double psd[2], cr;
		//	---------------------------------
	psd[0] = 2.969;3.009429;
	psd[1] = 2.501;2.681535;

	z_start[0] = 0.55;
	z_start[1] = 0.5;
	calcModel(psd, z_start, start_q);
	cr = calcCriteria();
	for(i=1; i<9; i++)
		csb[i] = calcConstraint(i);
		//	---------------------------------
	psd[0] -= dd*psd[0];
	calcModel(psd, z_start, start_q);
	for(i=1; i<9; i++)
	{
		csv1[i] = calcConstraint(i);
		dcsdv1[i] = (csb[i] - csv1[i]) / dd;
	}
	psd[0] += dd*psd[0];
		//	---------------------------------
	psd[1] -= dd*psd[1];
	calcModel(psd, z_start, start_q);
	for(i=1; i<9; i++)
	{
		csv2[i] = calcConstraint(i);
		dcsdv2[i] = (csb[i] - csv2[i]) / dd;
	}
	psd[1] += dd*psd[1];
		//	---------------------------------
	z_start[0] -= dd*z_start[0];
	calcModel(psd, z_start, start_q);
	for(i=1; i<9; i++)
	{
		csrt1[i] = calcConstraint(i);
		dcsdrt1[i] = (csb[i] - csrt1[i]) / dd;
	}
	z_start[0] += dd*z_start[0];
		//	---------------------------------
	z_start[1] -= dd*z_start[1];
	calcModel(psd, z_start, start_q);
	for(i=1; i<9; i++)
	{
		csrt2[i] = calcConstraint(i);
		dcsdrt2[i] = (csb[i] - csrt2[i]) / dd;
	}
	z_start[1] += dd*z_start[1];
}

double varSpread(double bound, double spread, double leftBound, double rightBound, int isRight)
{
	spread *= rightBound - leftBound;
	if(isRight)
		if(rightBound < bound+spread)
			return rightBound;
		else
			return bound+spread;
	else
		if(leftBound > bound-spread)
			return leftBound;
		else
			return bound-spread;
}