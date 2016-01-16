#include "..\headers\my.h"

int MaximizeConstraint(								// ���������� ifail ������
					   FILE *F,						// ���� ������
					   struct T *pT,				// �������, � ������� ��������������� �����������
					   int ConstraintToMaximize,	// ����� ����������� ��� ������������
					   double *ConstraintValue,		// �������� ������������������ ����������
					   double *tetamax,				// ���������� ����� ���������
					   double U,
					   double gamma,
					   int iPrintFlag				// ������ �� ����� ������������� ����������
					   )
{
	int iq;


	int i_komod;
	static double *x,fc,*f,*xmin,*xmax,*multX,*multF,*lam,*ko,*komod,*teta;
	static int *mzamp; static int n,m,me;
	int can,parm;
	int parmget,nAl,mAl;
	char sz_FileName[80], sz_functrName[80], sz_FileShortName[80], sz_tmp1[80], sz_tmp2[80];
	int i,k;
	double *myx,*myxmin,*myxmax,*mykomod;

		//	----------------
	ROPUD_pT = pT;	//	�������� ��������� �������������� ������� ��� ������ maximize() 
		//	----------------


	myx		= (double*)malloc(Nq*sizeof(double));
	myxmin	= (double*)malloc(Nq*sizeof(double));
	myxmax	= (double*)malloc(Nq*sizeof(double));

//	������� ��������� �����
	for(iq=0; iq<_Nq; iq++)
	{
		myx[iq] = 0.001;
		myxmin[iq] = 0.0;
		myxmax[iq] = 1.0;
	}

	can=2;
	parm=2; 
	parmget=1;

		//	������� ����� ����� ��� ������� ������ ������������ �����������
	strcpy(sz_tmp1,"f");
	itoa(ConstraintToMaximize,sz_tmp2,10);
	strcat(sz_tmp1,sz_tmp2);
		//	������� ��������� ����� ��� ������ ����� getname()
	strcpy(sz_FileShortName,sz_tmp1);
	strcat(sz_tmp1,".dat");
	strcpy(sz_FileName,sz_tmp1);
		//	������� ����� ������� ��� ������� functr()	
	strcpy(sz_tmp1,"max");
	itoa(ConstraintToMaximize,sz_tmp2,10);
	strcat(sz_tmp1,sz_tmp2);
	strcpy(sz_functrName,sz_tmp1);

	if(!CreateDatFile(sz_FileName,"max",_Nq,0,0,0,0,2,myx,myxmin,myxmax,NULL))
	{
		getname("",sz_FileShortName,"dat","rez",can,parm);
		getmem(&x,&xmin,&xmax,&f,&multX,&multF,&lam,&mzamp,&ko,&komod,&teta,can,&nAl,&mAl);
		callropud1(can,parm,parmget,&x,&fc,&f,&xmin,&xmax,&multX,&multF,&lam,&mzamp,&ko,&komod,&teta,&n,&m,&me);
// ��� ��� ������ max(fi)=-min(-fi)
		*ConstraintValue=-fc;
		//pT->cv[ConstraintToMaximize-1]=-fc;

		nfunM+=(int)ko[31];
		nfungM+=(int)ko[36];

		if(iPrintFlag)
		{
			fprintf(F,"\nConstraint #%d for region with lower bound (",ConstraintToMaximize);
			for(i=0;i<_Nq;i++)
			{
				fprintf(F," %.8f",pT->lbound[i]);
			}
			fprintf(F," ) and upper bound (");
			for(i=0;i<_Nq;i++)
			{
				fprintf(F," %.8f",pT->rbound[i]);
			}
			fprintf(F," ) is %.8f, teta[0]=%.8f",fc,x[0]);
			for(i=1;i<_Nq;i++)
			{
				fprintf(F," teta[%d]=%.8f",i,x[i]);
			}
			fprintf(F,"\nNumber of function calculations is %.0f, number of direction calculations is %.0f.",ko[31],ko[36]);
		}

//	���������� ����� ���������� ����� ������������ ����������� �������� �����
		for(iq=0; iq<_Nq; iq++)
			tetamax[iq] = x[iq];
		i=(int)ko[75];
		put_end(can);
	}
	else
	{
		printf("\n\n!!!ERROR: error while creating %s file\n",sz_FileName);
		exit(10);
	}
	free(myx);
	free(myxmin);
	free(myxmax);
	return i;
}

int MaximizeAllConstraints(
						   FILE *F,						// ���� ������
						   struct T *proot,				// ��������� �� ������ ��������
						   double **S,					// ������ ����������� ����� (����������� Nq*nkrit (��������=NULL)
						   int *nCrit,					// ��������� ����� ���� ����� �� ���� �������� (=0)
						   double U,					// U
						   int Nd,
						   double gamma,
						   int iPrintFlag				// ������ �� ����� ������������� ����������
						   )
{
	struct T *pT;
	int j,q;
	double max;
	double *tetamax;
	int totallength=0;
	int ifail;
	char fn[10];

// �� ���� ��������
	pT = proot;
	*nCrit = 0;
	while(pT=pT->pnext)
	{
		//if(!pT->isLocked)
		{
			tetamax = (double*)malloc(Nq*sizeof(double));
			if(tetamax==NULL) exit(4);

			j = pT->NCons;
			if(j==12)
				j=j;

			pT->isActive = 0;

			ifail = MaximizeConstraint(F,pT,j,&max,tetamax,U,gamma,iPrintFlag);
			if(ifail>1) return ifail;
			if(max>=5e-7)
			{
				//if(j==0)
					for(q=0; q<Nq; q++)
						if(tetamax[q]<0.05)
							tetamax[q] = 0;
						else if(tetamax[q]>0.95)
							tetamax[q] = 1;
				AddQcr(pT, tetamax);
				if((*S)==NULL)
					(*S)=(double*)malloc(Nq*sizeof(double));
				else
					(*S)=(double*)realloc((*S),(totallength+Nq)*sizeof(double));
				if((*S)==NULL) exit(4);
				for(q=0;q<Nq;q++)
					(*S)[q+totallength]=tetamax[q];
				totallength+=q;
				if(!(isLocked(pT)))	//	����� ����� ������������� �������� �� ������ ���������� ���� ���
				{
					(*nCrit)++;
					for(q=0; q<Nq; q++)
						printf("\t%lf", tetamax[q]);
				}

			}
			if(fabs(max)<=EPS_MAXFI)
			{
				pT->isActive=1;
				//pT->ActiveConstraints[j]=1;
			}
		

			sprintf(fn, "f%d.dat", j);
			remove(fn);
			sprintf(fn, "f%d.rez", j);
			remove(fn);

			free(tetamax);
		}
	}// �� ���� ��������

	return ifail;
}

int MaximizeAllCriteria(
						   FILE *F,						// ���� ������
						   struct TK *proot,				// ��������� �� ������ ��������
						   double U,					// U
						   double **d,
						   int Nd,
						   double gamma,
						   int iPrintFlag				// ������ �� ����� ������������� ����������
						   )
{
	struct TK *pT;
	int j,q;
	double max;
	double *tetamax;
	int totallength=0;
	int ifail;

	int CritNum;
	double CritVal;
	int first=1;

	char fn[10];
	
// �� ���� ��������
	j=0;
	pT=proot;
	while(pT=pT->pnext)
	{
		tetamax = (double*)malloc(Nq*sizeof(double));
		if(!tetamax) exit(4);

		j++;

		pT->isActive=0;

		ifail=MaximizeCriteria(F,pT,j,&max,tetamax,U,gamma,d,Nd,iPrintFlag);
		sprintf(fn, "k%d.dat", j);
		remove(fn);
		sprintf(fn, "k%d.rez", j);
		remove(fn);

		if(first)
		{
			first=0;
			CritVal=max;
			CritNum=1;
		}
		
		if(max>CritVal)
		{
			CritVal=max;
			CritNum=j;
		}

		free(tetamax);
	}// �� ���� ��������
		
// �������� �����������
	pT=proot;
//	pT=pT->pnext;
	for(j=0;j<CritNum;j++)
	{
		pT=pT->pnext;
	}	

	pT->isActive=1;

	return ifail;
}

int MaximizeCriteria(								// ���������� ifail ������
					   FILE *F,						// ���� ������
					   struct TK *pT,				// �������, � ������� ��������������� �����������
					   int ConstraintToMaximize,	// ����� ����������� ��� ������������
					   double *ConstraintValue,		// �������� ������������������ ����������
					   double *tetamax,				// ���������� ����� ���������
					   double U,
					   double gamma,
					   double **d,
					   int Nd,
					   int iPrintFlag				// ������ �� ����� ������������� ����������
					   )
{
	int i_komod;
	static double *x,fc,*f,*xmin,*xmax,*multX,*multF,*lam,*ko,*komod,*teta;
	static int *mzamp; static int n,m,me;
	int can,parm;
	int parmget,nAl,mAl;
	char sz_FileName[80], sz_functrName[80], sz_FileShortName[80], sz_tmp1[80], sz_tmp2[80];
	int i,k;
	double *myx,*myxmin,*myxmax,*mykomod;

		//	----------------
	ROPUD_pTK = pT;	//	�������� ��������� �������������� ������� �������� ��� ������ �23()
		//	----------------



	myx=(double*)malloc(Nq*sizeof(double));
	myxmin=(double*)malloc(Nq*sizeof(double));
	myxmax=(double*)malloc(Nq*sizeof(double));

//	������� ��������� �����
	for(k=0; k<Nq; k++)
	{
		myx[k] = 0.001;
		myxmin[k] = 0.0;
		myxmax[k] = 1.0;
	}

	can=1;
	parm=2; 
	parmget=1;
			

//	������� ����� ����� ��� ������� ������ ������������ �����������
	strcpy(sz_tmp1,"k");
	itoa(ConstraintToMaximize,sz_tmp2,10);
	strcat(sz_tmp1,sz_tmp2);
//	������� ��������� ����� ��� ������ ����� getname()
	strcpy(sz_FileShortName,sz_tmp1);
	strcat(sz_tmp1,".dat");
	strcpy(sz_FileName,sz_tmp1);
//	������� ����� ������� ��� ������� functr()	
/*	strcpy(sz_tmp1,"max");
	itoa(ConstraintToMaximize,sz_tmp2,10);
	strcat(sz_tmp1,sz_tmp2);
	strcpy(sz_functrName,sz_tmp1);
*/
	if(!CreateDatFile(sz_FileName,"KU",Nq,1,0,0,0,3,myx,myxmin,myxmax,NULL))
	{
		getname("",sz_FileShortName,"dat","rez",can,parm);
		getmem(&x,&xmin,&xmax,&f,&multX,&multF,&lam,&mzamp,&ko,&komod,&teta,can,&nAl,&mAl);
		callropud1(can,parm,parmget,&x,&fc,&f,&xmin,&xmax,&multX,&multF,&lam,&mzamp,&ko,&komod,&teta,&n,&m,&me);
// ��� ��� ������ max(fi)=-min(-fi)
		*ConstraintValue=-fc;
//		pT->cv[ConstraintToMaximize-1]=-fc;
		
		if(iPrintFlag)
		{
			fprintf(F,"\nConstraint #%d for region with lower bound (",ConstraintToMaximize);
			for(i=0;i<Nq;i++)
			{
				fprintf(F," %.8f",pT->lbound[i]);
			}
			fprintf(F," ) and upper bound (");
			for(i=0;i<Nq;i++)
			{
				fprintf(F," %.8f",pT->rbound[i]);
			}
			fprintf(F," ) is %.8f, teta[0]=%.8f",fc,x[0]);
			for(i=1;i<Nq;i++)
			{
				fprintf(F," teta[%d]=%.8f",i,x[i]);
			}
			fprintf(F,"\nNumber of function calculations is %.0f, number of direction calculations is %.0f.",ko[31],ko[36]);
		}

//	���������� ����� ���������� ����� ������������ ����������� �������� �����
		for(k=0;k<Nq;k++)
			tetamax[k]=x[k];
		i=(int)ko[75];
		put_end(can);
	}
	else
	{
		printf("\n\n!!!ERROR: error while creating %s file\n",sz_FileName);
		exit(10);
	}
	free(myx);
	free(myxmin);
	free(myxmax);
	//free(mykomod);
	return i;
}