/*******************************************************************************************************************************
Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

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

#ifndef ENERGY_OPTIMIZATION_H
#define ENERGY_OPTIMIZATION_H

#include "ProteinDesign.h"
#include "DEE.h"
#include "EnergyMatrix.h"
#include "SmallMol.h"
#include "Sequence.h"

// simulated annealing
#define SA_CYCLE         3
#define	SA_TMAX          10.0
#define	SA_TMIN	         0.01
#define	SA_DECREASE_FAC  0.8
#define METROPOLIS_STEP  50000

int SitePairConsDeploy(CataConsSitePairArray *pSitePairArray, Structure *pStructure);
int EnergyMatrixUpdateForCataCons(EnergyMatrix *pMatrix, RotamerList* pList, CataConsSitePair *pSitePair, Structure *pStructure);
int EnergyMatrixUpdateForCataConsArray(EnergyMatrix *pMatrix, RotamerList* pList, CataConsSitePairArray *pSitePairArray, Structure *pStructure);
int DesignSiteShowRotamerTypeAndCount(RotamerList* pList, Structure* pStructure, StringArray** ppRotamerType, IntArray** ppRotamerCount);

int SequenceRandomSiteIndex(int *mutSiteIndex,int designSiteCount);
int SequenceRandRotamerIndex(Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int siteIndex,int* rotIndex);
int SequenceChangeSingleSite(Sequence* pThis, int mutationSiteIndex, int mutationRotamerIndex);

int SequenceGenerateRandomSeed(Sequence* pThis, RotamerList* pList);
int SequenceGenerateInitialSequenceSeed(Structure* pStructure, RotamerList* pList,Sequence* pThis);
int SequenceGenerateNativeSequenceSeed(Structure* pStructure, RotamerList* pList,Sequence* pSequence);

int StructureCalcSequenceTemplateEnergy(Structure* pStructure, Sequence* pSequence,double energyTerms[MAX_ENERGY_TERM],double energyTermsBind[MAX_ENERGY_TERM]);
int StructureCalcSequenceEnergy(Structure* pStructure,Sequence* pSequence);
int StructureCalcEnergyChangeUponSingleMutation(Structure* pStructure,Sequence* pSequence,int mutSiteIndex,int mutRotIndex,double* dtot,double* dphy,double* dbin,double *devo);
int MetropolisCriterion(Sequence* pOld,Sequence* pBest,Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int *seqIndex,double temp,int stepCount, FILE* fp,FILE *fp2);

int StructureCalcSequenceEnergyWithCataCons(Structure* pStructure,Sequence* pSequence,CataConsSitePairArray* pConsArray);
int StructureCalcEnergyChangeUponSingleMutationWithCataCons(Structure* pStructure,Sequence* pSequence,int mutSiteIndex,int mutRotIndex,double* dtot,double* dphy,double* dbin,double *devo,int *dcons,CataConsSitePairArray* pConsArray);
int MetropolisCriterionWithCataCons(Sequence* pOld,Sequence* pBest,Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int *seqIndex,double temp,int stepCount, FILE* pFileRot,FILE *pFileSeq,CataConsSitePairArray* pConsArray);

int SimulatedAnnealingOptimization(Structure* pStructure,RotamerList* pList);

#endif
