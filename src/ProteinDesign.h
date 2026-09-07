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

#ifndef PROTEIN_DESIGN_H
#define PROTEIN_DESIGN_H

#include "RotamerBuilder.h"
#include "SmallMol.h"
#include "Sequence.h"

int SitePairConsDeploy(CataConsSitePairArray* pConsArray, Structure* pStruct);
int DesignSiteShowRotamerTypeAndNum(RotamerList* pList, Structure* pStruct, StringArray** ppRotTypes, IntArray** ppRotNums);

int SequenceRandomSiteIndex(int* mutSiteNdx, int desSiteNum);
int SequenceRandRotamerIndex(Structure* pStruct, RotamerList* pList, StringArray** ppRotTypes, IntArray** ppRotNums, int siteNdx, int* rotNdx);
int SequenceUpdateSingleSite(Sequence* pSeq, int mutSiteNdx, int mutRotNdx);

int SequenceGenRandomSeed(Sequence* pSeq, RotamerList* pList);
int SequenceGenInitialSeqSeed(Sequence* pSeq, Structure* pStruct, RotamerList* pList);
int SequenceGenNativeSeqSeed(Sequence* pSeq, Structure* pStruct, RotamerList* pList);

int SequenceCalcTemplateEnergy(Sequence* pSeq, Structure* pStruct, double *energyTerms, double* energyTermsBind);
int SequenceCalcEnergy(Sequence* pSeq, Structure* pStruct, Sequence** ppSavedSeqs, int numSavedSeqs);
int SequenceCalcEnergyWithCataCons(Sequence* pSeq, Structure* pStruct, CataConsSitePairArray* pConsArray, Sequence** ppSavedSeqs, int numSavedSeqs);
int SequenceCalcEnergyChangeUponSingleMutation(Sequence* pSeq, Structure* pStruct, int mutSiteNdx, int mutRotNdx, double* dtot, double* dphy, double* dbin, double* devo, int* nmutsAfter, Sequence** ppSavedSeqs, int numSavedSeqs, BOOL* isSeqSameAsSaved);
int SequenceCalcEnergyChangeUponSingleMutationWithCataCons(Sequence* pSeq, Structure* pStruct, CataConsSitePairArray* pConsArray, int mutSiteNdx, int mutRotNdx, double* dtot, double* dphy, double* dbin, double* devo, int* dcons, int* nmutsAfter, Sequence** ppSavedSeqs, int numSavedSeqs, BOOL* isSeqSameAsSaved);

int SequenceCalcNumOfMutations(Sequence* pSeq, Structure* pStruct);
int SequenceCalcNumOfMutationsUponSingleMutation(Sequence* pSeq, Structure* pStruct, int mutSiteNdx, int mutRotNdx);
double EnergyPenaltyForOutOfRangeMutations(int count);
BOOL SequenceIsSameAsSavedSequences(Sequence* pSeq, Structure* pStruct, Sequence** ppSavedSeqs, int numSavedSeqs);
BOOL SequenceIsSameAsSavedSequencesUponSingleMutation(Sequence* pSeq, Structure* pStruct, int mutSiteNdx, int mutRotNdx, Sequence** ppSavedSeqs, int numSavedSeqs);

int Metropolis(Sequence* pSeq, Sequence* pBest, Structure* pStruct, RotamerList* pList, StringArray** ppRotTypes, IntArray** ppRotNums, double T, int stepNum, Sequence** ppSavedSeqs, int numSavedSeqs);
int MetropolisWithCataCons(Sequence* pSeq, Sequence* pBest, Structure* pStruct, CataConsSitePairArray* pConsArray, RotamerList* pList, StringArray** ppRotTypes, IntArray** ppRotNums, double T, int stepNum, Sequence** ppSavedSeqs, int numSavedSeqs);
int SimulatedAnnealing(Structure* pStruct, RotamerList* pList);

#endif
