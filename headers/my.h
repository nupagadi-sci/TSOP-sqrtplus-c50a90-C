#ifndef __MY_H__
#define __MY_H__
#include "..\headers\ropud.h"
#define _Nj 13
#define _Nt 30
#define _Nq 4
#define _Nd 2
#define _Nz 2
#define _Na 1
#define EPS_MAXFI 1e-07
#define EPS_UL	1e-02
#define EPS_T 1e-02
#define	K	0.0001
#define KA 1000
#define RMIN 0.01

#define EPS_QCR 1e-06

extern int nfunU,nfungU,nfunM,nfungM,nfunL,nfungL;
// ��������������� ��� ������ ���������� ����������
//#define DEBUG 0

extern int Nit, Nz, Nq, Nd, Nj, Na, Nt, Nap, Ns, is_linearization, is_b_min_inheritance, is_b_split_inheritance;
extern double z_start[_Nz], zmin[_Nz], zmax[_Nz], **b_start, bmin, bmax, d[_Nd], dmin[_Nd], dmax[_Nd], **E, **A, *P, Rabs[_Nq], alpha[_Na];
extern double start_q[_Nq], start_minq[_Nq], start_maxq[_Nq];
extern double gamma, mltf[5];
extern struct T *ROPUD_pT;
extern struct TK *ROPUD_pTK;

extern double C;
//extern FILE *log;
//extern int ROPUD_ConstraintToMaximize;

struct T
{
	int id;
	//	//	������ ����������� ���������� 
	//double *z;
	//	//	������� ������������� ������������ ����������� ����������
	//	//	b[Nz][Nq+1]; b[i] - ������������ i-�� ������������ ���������
	//double **b;	//	b[i][Nq] - ��������� ����
//�������
	double *lbound;
	double *rbound;
//��������� �������� ������
	double *slbound;
	double *srbound;
	int isSoftCons;	//	������ �� �����������
//����������� �� ������� (���������)
	//double alpha;
//����� �����������, �������� ����������� ��� �������
	int NCons;
//	����������� �����, ����. �������� ���������� � �������, ���� ����������
	double *Qcr;
	int NQcr;
	//double *cv;
	int isActive;
	//int *ActiveConstraints;
	int isLocked;

	//int *fork;	//	���������
		//	���������� ������ ��������� ���-�. ���-�� �� ��-�� ���-�� ">0" ��� ��������
	//int NumExtraSplits;	//	��������� �� �������� - ������������, � ������� ����������� �� 1
//	������ ������� � ����� ����. �����	
	double r;
	int qr;
//	����������� ��� ������������ � ������� �������� (����� �������� ��� �������� � ����� ��������)
//	� �������� ��� ����������� =0 � �� �� ��� �� ������
	struct TK *relK;	//	������� ��������, ������� ����������� ������ ������� �����������
	struct T *pnext;
};

struct TK
{
	int id;
		//	������ ����������� ���������� 
	double *z;
		//	������� ������������� ������������ ����������� ����������
		//	b[Nz][Nq+1]; b[i] - ������������ i-�� ������������ ���������
	double **b;	//	b[i][Nq] - ��������� ����
		//	�������
	double *lbound;
	double *rbound;

	int iE[_Nq];	//	������ ����������� � ���������� ������� E
	double E[_Nq];

	int isActive;
		//	������ ������� � ����� ����. �����	
	double r;
	int qr;
		//	����������� ��� ������������ � ������� �������� (����� �������� ��� �������� � ����� ��������)
		//	� �������� ��� ����������� =0 � �� �� ��� �� ������
	double a;
	double ai[_Nq];
	struct TK *pnext;
};

struct BinTree
{
	struct T *region;
	struct BinTree *lBranch;
	struct BinTree *rBranch;
};

struct TSplit
{
	double *lbound;
	double *rbound;
	int iteration;
	struct TSplit *pnext;
};

	//	������������� ���� ������� ������ ����� ��������������� ������� pTK �� �������� 
void replant(struct BinTree *BT, struct TK *pTK);
	//	������������� ���������� z, ��� ������������ b � ����� q
int set_z(struct TK *pT, double* q);
	//	�������� ������ �� �������� id �������� ��������, ����������� ��� �������� ����� ����� ���� �������,
int getNewID(int *lid);		//	���������� �������� ������ �������

struct TSplit* CreateTS(double *lbound, double *rbound, int Nq, int i);

struct T* CreateT(int Nq, int Nj, double *lbound, double *rbound, double *Rabs);
struct TK* CreateTK(double *lbound, double *rbound);
void DeleteT(struct T **pT);
void CalculateSize(struct T *pT,double *Rabs);
void CalculateSize2(struct T *pT);
struct T* Split(struct T *pT1, double *Rabs);
int CreateNewQcr(int NQcr, double *pQcr, struct T *pT);
void DeleteList(struct T **proot);
int AssignConstraintsValues(int t, int Nj, double *cv, struct T *pT);
int AddQcr(struct T *pT, double *pQcr);
int DeleteQcr(struct T *pT);

int start();

int CreateDatFile(
				  char cdf_sz_name[], 
				  char cdf_sz_number[], 
				  int cdf_n, 
				  int cdf_m, 
				  int cdf_me, 
				  int cdf_nmod, 
				  int cdf_nteta, 
				  int cdf_prTask, 
				  double cdf_x[], 
				  double cdf__xmin[], 
				  double cdf_xmax[], 
				  double cdf_komod[]
				  );

int MaximizeConstraint(								// ���������� ifail ������
					   FILE *F,						// ���� ������
					   struct T *pT,				// �������, � ������� ��������������� �����������
					   int ConstraintToMaximize,	// ����� ����������� ��� ������������
					   double *ConstraintValue,		// �������� ������������������ ����������
					   double *tetamax,				// ���������� ����� ���������
					   double U,
					   double gamma,
					   int iPrintFlag				// ������ �� ����� ������������� ����������
					   );

int MaximizeAllConstraints(
						   FILE *F,						// ���� ������
						   struct T *proot,				// ��������� �� ������ ��������
						   double **S,					// ������ ����������� ����� (����������� Nq*nkrit (��������=NULL)
						   int *nCrit,					// ��������� ����� ���� ����� �� ���� �������� (=0)
						   double U,					// U
						   int Nd,
						   double gamma,
						   int iPrintFlag				// ������ �� ����� ������������� ����������
						   );

int mvadzo(											// ���������� ifail ������
					   FILE *F,						// ���� ������
					   struct T *proot,				// ������� ��� ������
					   struct TK *prootK,			// ������� ��������
					   double **S,					// ������������ ����������� �����
					   int *nCrit,					// ������������ ����� ����������� �����
					   double *U,					// ������������ �������� ��������
					   int Nd,
					   double gamma,
					   int Nap,
					   double *pQap,
					   int Na,
					   double *alpha,
					   int iPrintFlag			// ������ �� ����� ������������� ����������
					   );

int SolveMinUdzo(									// ���������� ifail ������
					   FILE *F,						// ���� ������
					   struct T *proot,				// ������� ��� ������
					   struct TK *prootK,			// ������� ��������
					   double *U,					// 
					   double *pQap,
					   double *alpha,
					   int iPrintFlag			// ������ �� ����� ������������� ����������
					   );


void SplitT(int Nt, int Ntsplit,int *V, int *qsplit, double *array_teta,
	   double *tetamin, double *tetamax, double *z1,
	   double *Qcr, int *NQcr);



int MaximizeAllCriteria(
						   FILE *F,						// ���� ������
						   struct TK *proot,				// ��������� �� ������ ��������
						   double U,					// U
						   double **d,
						   int Nd,
						   double gamma,
						   int iPrintFlag				// ������ �� ����� ������������� ����������
						   );


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
					   );

void upperBound(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,int n,int m,int me,double ko[],double komod[],double teta[]);
void maxCriteria(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,int n,int m,int me,double ko[],double komod[],double teta[]);
void maximize(double x[],double *pfc,double f[],double g[],DOUBLE **dfdx,int n,int m,int me,double ko[],double komod[],double teta[]);

double varSpread(double bound, double spread, double leftBound, double rightBound, int isRight);

double calcCriteria();
double calcConstraint(int NCons);
double calcConstraintOnZ(int NCons);

void DisplaySearchVarsNums(FILE* fptr);
void DisplayConstraintsNums(FILE* fptr, struct T* prootU);

#endif
