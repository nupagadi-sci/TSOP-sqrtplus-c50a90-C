#include "..\headers\my.h"

static double V1, V2;	//	Конструктивные
static double RT1, RT2;	//	Управляющие
static double E1, E2, K10, K20;	//	Неопределенные
static double K1, K2, K3, K4, CA1, CA2, CB1, CB2;	//	Прочие
double C;	//	Параметр



void calcModel(double *d, double *z, double *Q)
{
	V1 = d[0];
	V2 = d[1];
	RT1 = z[0];
	RT2 = z[1];
	E1 = Q[0];
	E2 = Q[1];
	K10 = Q[2];
	K20 = Q[3];

	K1 = K10*exp(-E1/RT1);
	K2 = K10*exp(-E1/RT2);
	K3 = K20*exp(-E2/RT1);
	K4 = K20*exp(-E2/RT2);

	CA1 = 1/ (1+K1*V1);
	CA2 = CA1/ (1+K2*V2);
	CB1 = (1-CA1) / (1+K3*V1);
	CB2 = (CB1+CA1-CA2) / (1+K4*V2);
}

double* calcDerivative(double *d, double *z, double **b, double *Q)
{
	double dF_dQ[_Nq];

	return dF_dQ;
}

double calcConstraint(int NCons)
{
	double f;
	switch(NCons)
	{
			//	Мягкие
	case 0:
		f = C - CB2;
		break;

			//	Жесткие
	case 1:
		f = -CA1;
		break;
	case 2:
		f = CA1 - 1;
		break;
	case 3:
		f = -CA2;
		break;
	case 4:
		f = CA2 - 1;
		break;
	case 5:
		f = -CB1;
		break;
	case 6:
		f = CB1 - 1;
		break;
	case 7:
		f = -CB2;
		break;
	case 8:
		f = CB2 - 1;
		break;
	case 9:
		f = zmin[0] - RT1;
		break;
	case 10:
		f = RT1 - zmax[0];
		break;
	case 11:
		f = zmin[1] - RT2;
		break;
	case 12:
		f = RT2 - zmax[1];
		break;
	}
	return f;
}

double calcConstraintOnZ(int NCons)
{
	double f;
	switch(NCons)
	{
	case 0:
		f = zmin[0] - RT1;
		break;
	case 1:
		f = RT1 - zmax[0];
		break;
	case 2:
		f = zmin[1] - RT2;
		break;
	case 3:
		f = RT2 - zmax[1];
		break;
	}
	return f;
}

double calcCriteria()
{
	return sqrt(V1) + sqrt(V2) + 0*(RT1+RT2);
}
