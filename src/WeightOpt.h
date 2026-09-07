/*******************************************************************************************************************************
Copyright (c) Xiaoqiang Huang

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/

#ifndef WEIGHT_OPT_H
#define WEIGHT_OPT_H

#define  NUM_ENERGY_TERM 72
#define  NUM_REFERENCE   20
#define  NUM_ENERGY_MON  46 //20 ref + 8 intra + 3 statistics + 15 interS
#define  NUM_ENERGY_PPI  15 //5terms + 1 ssbond + 9 HB terms
#define  NUM_ENERGY_LIG  11 //5terms + 6 HB terms

#define  DELTA_X        1e-4
#define  EPSILON        1e-3
#define  WOLFE_C1       1e-4
#define  WOLFE_C2       0.1
#define  MAX_STEP       10000


int isNumber(double x);
int isFiniteNumber(double d);


typedef struct _EnergyTermRot
{
  int rotTypeNum[2]; // two rotamer types: NAT and ROT
  double** energy[2];
}EnergyTermRot;

int EnergyTermRotCreate(EnergyTermRot* pTerm);
int EnergyTermRotDestroy(EnergyTermRot* pTerm);


typedef struct _EnergyTermsRot
{
  int posNum;
  EnergyTermRot* terms;
}EnergyTermsRot;

int EnergyTermsRotCreate(EnergyTermsRot* pTerms);
int EnergyTermsRotDestroy(EnergyTermsRot* pTerms);
int EnergyTermsRotRead(EnergyTermsRot* terms, char* pdblist);

int GradientUnit(double* grad, int n);
double GradientNorm(double* grad, int n);
double GradientMax(double* grad, int n);
int Xupdate(double* x, double* grad, int n, double alpha);
int Xshow(double* x, int n);

int PNATROT_BacktrackLineSearch(EnergyTermsRot* pTerms, double lossk, double* xk, double* gradk, double scalar, int n);
int PNATROT_LossFunction(EnergyTermsRot* pTerms, double* x, int n, double* loss);
int PNATROT_CalcGradientTwoSide(EnergyTermsRot* pTerms, double loss, double* x, double* grad, int n);
int PNATROT_CalcGradientForward(EnergyTermsRot* pTerms, double loss, double* x, double* grad, int n);
int PNATROT_CalcGradientBackword(EnergyTermsRot* pTerms, double loss, double* x, double* grad, int n);
int PNATROT_WeightOptByGradientDescent(char* pdblistfile);


#define TYPE_TWENTYTWO 22


typedef struct _EnergyTermAA
{
  //21 types: ACD...Y and NAT rotamer
  //index 0-19 are A,C,D,...,Y and index 21 is NAT
  int rotTypeNum[TYPE_TWENTYTWO];
  double** energy[TYPE_TWENTYTWO];
}EnergyTermAA;

int EnergyTermAllCreate(EnergyTermAA* pTerm);
int EnergyTermAllDestroy(EnergyTermAA* pTerm);

typedef struct _EnergyTermsAA
{
  int posNum;
  EnergyTermAA* terms;
}EnergyTermsAA;

int EnergyTermsAACreate(EnergyTermsAA* pTerms);
int EnergyTermsAADestroy(EnergyTermsAA* pTerms);
int EnergyTermsAARead(EnergyTermsAA* pTerms, char* pdblist);
int PNATAA_BacktrackLineSearch(EnergyTermsAA* pTerms, double lossk, double* xk, double* gradk, double scalar, int n);
int PNATAA_LossFunction(EnergyTermsAA* pTerms, double* x, int n, double* loss);
int PNATAA_CalcGradientForward(EnergyTermsAA* pTerms, double loss, double* x, double* grad, int n);
int PNATAA_WeightOptByGradientDescent(char* pdblistfile);

#endif