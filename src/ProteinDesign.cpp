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

#pragma warning(disable:4244)
#pragma warning(disable:4305)
#include "ProteinDesign.h"
#include "Sequence.h"
#include "Evolution.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern BOOL FLAG_EVOLUTION;
extern BOOL FLAG_DESIGN_MON;
extern BOOL FLAG_DESIGN_PPI;
extern BOOL FLAG_DESIGN_PLI;
extern BOOL FLAG_DESIGN_ENZ;
extern BOOL FLAG_MOL2;
extern BOOL FLAG_EVOPHIPSI;
extern BOOL FLAG_DESIGN_FROM_NATAA;

extern char TAG_BEST_DESIGNS[MAX_LEN_FILE_NAME+1];
extern char TAG_BEST_SEQS[MAX_LEN_FILE_NAME+1];
extern char TAG_BEST_MODEL[MAX_LEN_FILE_NAME+1];
extern char TAG_BEST_LIG_MOL2[MAX_LEN_FILE_NAME+1];
extern char TAG_DES_ROTS[MAX_LEN_FILE_NAME+1];
extern char TAG_DES_SEQS[MAX_LEN_FILE_NAME+1];
extern char FILE_LIG_CATACONS[MAX_LEN_FILE_NAME+1];

extern int PROT_LEN_NORM;
extern int NTRAJ;
extern int NTRAJ_START_NDX;

extern char PDBID[MAX_LEN_FILE_NAME+1];
extern char DES_CHAINS[MAX_LEN_ONE_LINE_CONTENT+1];
extern double WGT_PROFILE;
extern float** PROT_PROFILE;
extern double WGT_BIND;

extern char PROGRAM_PATH[MAX_LEN_ONE_LINE_CONTENT+1];
extern char TGT_PRF[MAX_LEN_FILE_NAME+1];
extern char TGT_SA[MAX_LEN_FILE_NAME+1];
extern char TGT_SS[MAX_LEN_FILE_NAME+1];
extern char TGT_SEQ[MAX_LEN_FILE_NAME+1];
extern char TGT_PHIPSI[MAX_LEN_FILE_NAME+1];

extern int MIN_MUT_NUM;
extern int MAX_MUT_NUM;
extern BOOL FLAG_RESTRICT_NUM_OF_MUTATIONS;
extern BOOL FLAG_PENALIZE_IDENTICAL_SEQUENCES;

extern double SA_T_INITIAL;
extern double SA_T_FINAL;
extern int    SA_CYCLE;
extern double SA_DECREASE_FACTOR;
extern int    NUM_METROPOLIS_STEPS;

extern double WGT_CATA_CONS;
extern double WGT_MUT_NUM;
extern double WGT_SAME_SEQ;



int SitePairConsDeploy(CataConsSitePairArray* pConsArray, Structure* pStruct)
{
  for (int i = 0; i < CataConsSitePairArrayGetNum(pConsArray); i++) {
    CataConsSitePair* pCons = CataConsSitePairArrayGet(pConsArray, i);
    DesignSite* pFirst = StructureFindDesignSiteByChainName(pStruct, pCons->chnName1, pCons->pos1);
    DesignSite* pSecond = StructureFindDesignSiteByChainName(pStruct, pCons->chnName2, pCons->pos2);
    CataConsSitePairDeploy(pCons, DesignSiteGetRotamers(pFirst), DesignSiteGetRotamers(pSecond));
  }
  return Success;
}


int DesignSiteShowRotamerTypeAndNum(RotamerList* pList, Structure* pStruct, StringArray** ppRotTypes, IntArray** ppRotNums)
{
  StringArray* pTypeArray = *ppRotTypes;
  IntArray* pNumArray = *ppRotNums;

  printf("design site-specific rotamer details are as follows:\n");
  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++) {
    int rotNum = 0;
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    IntArrayCreate(&pNumArray[i], 0);
    StringArrayCreate(&pTypeArray[i]);
    for (int j = 0; j < pList->rotNum[i]; j++) {
      if (pList->remainFlag[i][j] == FALSE) continue;
      BOOL     isNewType = TRUE;
      Rotamer* pRot = RotamersGet(pSet, j);
      for (int k = 0; k < StringArrayGetNum(&pTypeArray[i]); k++) {
        if (strcmp(RotamerGetType(pRot), StringArrayGet(&pTypeArray[i], k)) == 0) {
          isNewType = FALSE;
          break;
        }
      }
      if (isNewType == TRUE) {
        int len = IntArrayGetLength(&pNumArray[i]);
        StringArrayAppend(&pTypeArray[i], RotamerGetType(pRot));
        len++;
        IntArrayResize(&pNumArray[i], len);
        IntArraySet(&pNumArray[i], len - 1, 1);
      }
      else {
        int pos, rotNum;
        StringArrayFind(&pTypeArray[i], RotamerGetType(pRot), &pos);
        rotNum = IntArrayGet(&pNumArray[i], pos);
        rotNum++;
        IntArraySet(&pNumArray[i], pos, rotNum);
      }
    }
    for (int j = 0; j < IntArrayGetLength(&pNumArray[i]); j++) rotNum += IntArrayGet(&pNumArray[i], j);
    printf("site %3d: %s %4d %3s, %6d rotamers:  ", i, ResidueGetChainName(pSite->pRes), ResidueGetPosInChain(pSite->pRes), ResidueGetName(pSite->pRes), rotNum);
    for (int j = 0; j < StringArrayGetNum(&pTypeArray[i]); j++) printf("%4d %s ", IntArrayGet(&pNumArray[i], j), StringArrayGet(&pTypeArray[i], j));
    printf("\n");
  }

  return Success;
}


int SequenceRandomSiteIndex(int* mutSiteNdx, int desSiteNum)
{
  *mutSiteNdx = rand() % desSiteNum;
  return Success;
}


int SequenceRandRotamerIndex(Structure* pStruct, RotamerList* pList, StringArray** ppRotTypes, IntArray** ppRotNums, int siteNdx, int* rotNdx)
{
  StringArray* pTypeArray = *ppRotTypes;
  IntArray* pNumArray = *ppRotNums;
  int typeNdx = rand() % IntArrayGetLength(&pNumArray[siteNdx]);
  int ndxInType = rand() % IntArrayGet(&pNumArray[siteNdx], typeNdx);
  int flag = 0;
  for (int i = 0; i < pList->rotNum[siteNdx]; i++) {
    Rotamer* pRot = RotamersGet(DesignSiteGetRotamers(StructureGetDesignSite(pStruct, siteNdx)), i);
    if (pList->remainFlag[siteNdx][i] == FALSE) continue;
    if (strcmp(StringArrayGet(&pTypeArray[siteNdx], typeNdx), RotamerGetType(pRot)) != 0) continue;
    if (flag == ndxInType) {
      *rotNdx = i;
      break;
    }
    flag++;
  }
  return Success;
}


int SequenceUpdateSingleSite(Sequence* pSeq, int siteNdx, int rotNdx)
{
  IntArraySet(&pSeq->rotNdxs, siteNdx, rotNdx);
  return Success;
}


int SequenceGenRandomSeed(Sequence* pSeq, RotamerList* pList)
{
  SequenceDestroy(pSeq);
  SequenceCreate(pSeq);
  pSeq->desSiteNum = pList->desSiteNum;
  pSeq->etot = 0;
  pSeq->ephy = 0;
  pSeq->eevo = 0;
  IntArrayResize(&pSeq->rotNdxs, pSeq->desSiteNum);
  for (int i = 0; i < pSeq->desSiteNum; i++) {
    int j = rand() % pList->remainRotNum[i];
    int trueJ;
    RotamerOriginalIndexGet(pList, i, j, &trueJ);
    IntArraySet(&pSeq->rotNdxs, i, trueJ);
  }

  return Success;
}


int SequenceGenInitialSeqSeed(Sequence* pSeq, Structure* pStruct, RotamerList* pList)
{
  SequenceDestroy(pSeq);
  SequenceCreate(pSeq);
  pSeq->desSiteNum = pList->desSiteNum;
  pSeq->etot = 0.0;
  IntArrayResize(&pSeq->rotNdxs, pSeq->desSiteNum);
  for (int i = 0; i < pSeq->desSiteNum; i++) {
    double minSelfEnergy = 1e8;
    int minSelfEnergyNdx = -1;
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    for (int j = 0;j < pList->rotNum[i];j++) {
      if (pList->remainFlag[i][j] == FALSE) continue;
      Rotamer* pRot = RotamersGet(pSet, j);
      if (pRot->selfEnergy < minSelfEnergy) {
        minSelfEnergy = pRot->selfEnergy;
        minSelfEnergyNdx = j;
      }
    }
    IntArraySet(&pSeq->rotNdxs, i, minSelfEnergyNdx);
  }

  return Success;
}


int SequenceGenNativeSeqSeed(Sequence* pSeq, Structure* pStruct, RotamerList* pList)
{
  SequenceDestroy(pSeq);
  SequenceCreate(pSeq);
  pSeq->desSiteNum = pList->desSiteNum;
  pSeq->etot = 0.0;
  IntArrayResize(&pSeq->rotNdxs, pSeq->desSiteNum);
  for (int i = 0; i < pSeq->desSiteNum; i++) {
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    for (int j = 0;j < RotamersGetNum(pSet);j++) {
      Rotamer* pRot = RotamersGet(pSet, j);
      RotamerRestore(pRot, pSet);
      if (RotamerAndResidueSidechainRMSD(pRot, pSite->pRes) < 0.05) {
        IntArraySet(&pSeq->rotNdxs, i, j);
        break;
      }
      RotamerExtract(pRot);
    }
  }

  return Success;
}


int SequenceCalcTemplateEnergy(Sequence* pSeq, Structure* pStruct, double* energyTerms, double* energyTermsBind)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
    Chain* pChn = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChn); j++) {
      Residue* pRes = ChainGetResidue(pChn, j);
      if (pRes->desType != Type_DesType_Fixed) {
        //the residue is a design site, set the coordinates of non-backbone atoms of this designable residue to FALSE
        for (int k = 0; k < ResidueGetAtomNum(pRes); k++) {
          Atom* pAtom = ResidueGetAtom(pRes, k);
          if (pAtom->isBBAtom == FALSE) pAtom->isXyzValid = FALSE;
        }
      }
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
    Chain* pChn = StructureGetChain(pStruct, i);
    for (int ir = 0; ir < ChainGetResidueNum(pChn); ir++) {
      Residue* pResIR = ChainGetResidue(pChn, ir);
      if (pResIR->desType == Type_DesType_Fixed && ChainGetType(pChn) == Type_Chain_Protein) { // calculate protein residue internal energy only
        energyTerms[91] += pResIR->aapp;
        energyTerms[92] += pResIR->rama;
        energyTerms[93] += ResidueGetDunbrack(pResIR);
        AminoAcidReferenceEnergy(ResidueGetName(pResIR), energyTerms);
        EnergyIntraResidue(pResIR, energyTerms);
      }
      for (int is = ir + 1; is < ChainGetResidueNum(pChn); is++) {
        Residue* pResIS = ChainGetResidue(pChn, is);
        if (ResidueGetPosInChain(pResIR) + 1 == ResidueGetPosInChain(pResIS)) EnergyResidueAndNextResidue(pResIR, pResIS, energyTerms);
        else EnergyResidueAndOtherResidueSameChain(pResIR, pResIS, energyTerms);
      }
      for (int k = i + 1; k < StructureGetChainNum(pStruct); k++) {
        Chain* pChnK = StructureGetChain(pStruct, k);
        for (int ks = 0; ks < ChainGetResidueNum(pChnK); ks++) {
          Residue* pResKS = ChainGetResidue(pChnK, ks);
          if (ChainGetType(pChn) == Type_Chain_SmallMol && ChainGetType(pChnK) != Type_Chain_SmallMol) {
            EnergyResidueAndLigandResidue(pResKS, pResIR, energyTerms);
            if (FLAG_DESIGN_PLI == TRUE || FLAG_DESIGN_ENZ == TRUE) EnergyResidueAndLigandResidue(pResKS, pResIR, energyTermsBind);
          }
          else if (ChainGetType(pChnK) == Type_Chain_SmallMol && ChainGetType(pChn) != Type_Chain_SmallMol) {
            EnergyResidueAndLigandResidue(pResIR, pResKS, energyTerms);
            if (FLAG_DESIGN_PLI == TRUE || FLAG_DESIGN_ENZ == TRUE) EnergyResidueAndLigandResidue(pResIR, pResKS, energyTermsBind);
          }
          else if (ChainGetType(pChnK) != Type_Chain_SmallMol && ChainGetType(pChn) != Type_Chain_SmallMol) {
            EnergyResidueAndOtherResidueDiffChain(pResIR, pResKS, energyTerms);
            if ((strstr(DES_CHAINS, ChainGetName(pChn)) != NULL && strstr(DES_CHAINS, ChainGetName(pChnK)) == NULL)
              || (strstr(DES_CHAINS, ChainGetName(pChn)) == NULL && strstr(DES_CHAINS, ChainGetName(pChnK)) != NULL))
            {
              if (FLAG_DESIGN_PPI == TRUE) EnergyResidueAndOtherResidueDiffChain(pResIR, pResKS, energyTermsBind);
            }
          }
        }
      }
    }
  }

  // restore the coordinates of non-backbone atoms of designable residues to TRUE
  for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
    Chain* pChn = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChn); j++) {
      Residue* pRes = ChainGetResidue(pChn, j);
      if (pRes->desType != Type_DesType_Fixed) {
        for (int k = 0; k < ResidueGetAtomNum(pRes); k++) {
          Atom* pAtom = ResidueGetAtom(pRes, k);
          if (pAtom->isBBAtom == FALSE) pAtom->isXyzValid = TRUE;
        }
      }
    }
  }
  return Success;
}


int SequenceCalcEnergy(Sequence* pSeq, Structure* pStruct, Sequence** ppSavedSeqs, int numSavedSeqs)
{
  pSeq->etot = 0, pSeq->eevo = 0, pSeq->ephy = 0, pSeq->ebin = 0;
  // 1. compute physical energy
  double energyTermsExBind[MAX_ENERGY_TERM] = { 0 };
  double energyTermsBind[MAX_ENERGY_TERM] = { 0 };
  SequenceCalcTemplateEnergy(pSeq, pStruct, energyTermsExBind, energyTermsBind);
  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++) {
    DesignSite* pSiteI = StructureGetDesignSite(pStruct, i);
    Rotamers* pSetI = DesignSiteGetRotamers(pSiteI);
    Rotamer* pRotI = RotamersGet(pSetI, IntArrayGet(&pSeq->rotNdxs, i));
    Chain* pChnI = StructureGetChain(pStruct, pSiteI->chnNdx);
    RotamerRestore(pRotI, pSetI);
    for (int k = i; k < StructureGetDesignSiteNum(pStruct); k++) {
      if (k == i) {
        pSeq->ephy += pRotI->selfEnergy;
        pSeq->ebin += pRotI->selfEnergyBin;
      }
      else {
        DesignSite* pSiteK = StructureGetDesignSite(pStruct, k);
        Rotamers* pSetK = DesignSiteGetRotamers(pSiteK);
        Rotamer* pRotK = RotamersGet(pSetK, IntArrayGet(&pSeq->rotNdxs, k));
        Chain* pChnK = StructureGetChain(pStruct, pSiteK->chnNdx);
        RotamerRestore(pRotK, pSetK);
        if (pSiteI->chnNdx == pSiteK->chnNdx) {// from the same chain
          if (ChainGetType(pChnI) == Type_Chain_Protein) EnergyRotamerAndRotamerSameChain(pRotI, pRotK, energyTermsExBind);
        }
        else
        {// from different chains
          if (ChainGetType(pChnI) == Type_Chain_SmallMol || ChainGetType(pChnK) == Type_Chain_SmallMol)
          {
            double energyTermsPL[MAX_ENERGY_TERM] = { 0 };
            EnergyRotamerAndLigRotamer(pRotI, pRotK, energyTermsPL);
            if ((strstr(DES_CHAINS, ChainGetName(pChnI)) != NULL && strstr(DES_CHAINS, ChainGetName(pChnK)) == NULL)
              || (strstr(DES_CHAINS, ChainGetName(pChnI)) == NULL && strstr(DES_CHAINS, ChainGetName(pChnK)) != NULL))
            {
              if (FLAG_DESIGN_PLI == TRUE || FLAG_DESIGN_ENZ == TRUE) {
                for (int index = 0; index < MAX_ENERGY_TERM; index++) energyTermsBind[index] += energyTermsPL[index];
              }
            }
          }
          else
          {
            double energyTermsDC[MAX_ENERGY_TERM] = { 0 };
            EnergyRotamerAndRotamerDiffChain(pRotI, pRotK, energyTermsDC);
            if ((strstr(DES_CHAINS, ChainGetName(pChnI)) != NULL && strstr(DES_CHAINS, ChainGetName(pChnK)) == NULL)
              || (strstr(DES_CHAINS, ChainGetName(pChnI)) == NULL && strstr(DES_CHAINS, ChainGetName(pChnK)) != NULL))
            {
              if (FLAG_DESIGN_PPI == TRUE) {
                for (int index = 0; index < MAX_ENERGY_TERM; index++) energyTermsBind[index] += energyTermsDC[index];
              }
            }
          }
        }
        RotamerExtract(pRotK);
      }
    }
    RotamerExtract(pRotI);
  }
  EnergyTermWeighting(energyTermsExBind);
  pSeq->ephy += energyTermsExBind[0];
  EnergyTermWeighting(energyTermsBind);
  pSeq->ebin += energyTermsBind[0];

  // 2. compute physical energy
  if (FLAG_EVOLUTION == TRUE) {
    char seqtxt[MAX_LEN_SEQ+1] = "";
    StructureGetWholeSequence(pSeq, pStruct, seqtxt);
    if (FLAG_EVOPHIPSI == TRUE) pSeq->eevo += EvolutionScoreAllFromSeq(seqtxt);
    else pSeq->eevo += EvolutionScoreFromPSSMWithoutAlignment(seqtxt);
  }

  pSeq->etot = WGT_PROFILE * pSeq->eevo + pSeq->ephy + WGT_BIND * pSeq->ebin;

  pSeq->mutNum = SequenceCalcNumOfMutations(pSeq, pStruct);
  if (FLAG_RESTRICT_NUM_OF_MUTATIONS == TRUE) pSeq->etot += EnergyPenaltyForOutOfRangeMutations(pSeq->mutNum);
  pSeq->seqExist = SequenceIsSameAsSavedSequences(pSeq, pStruct, ppSavedSeqs, numSavedSeqs);
  if (FLAG_PENALIZE_IDENTICAL_SEQUENCES == TRUE) pSeq->etot += WGT_SAME_SEQ * pSeq->seqExist;

  return Success;
}


int SequenceCalcEnergyWithCataCons(Sequence* pSeq, Structure* pStruct, CataConsSitePairArray* pConsArray, Sequence** ppSavedSeqs, int numSavedSeqs)
{
  pSeq->etot = 0, pSeq->eevo = 0, pSeq->ephy = 0, pSeq->ebin = 0;
  // 1. compute physical energy
  double energyTermsExBind[MAX_ENERGY_TERM] = { 0 };
  double energyTermsBind[MAX_ENERGY_TERM] = { 0 };
  SequenceCalcTemplateEnergy(pSeq, pStruct, energyTermsExBind, energyTermsBind);
  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++) 
  {
    DesignSite* pSiteI = StructureGetDesignSite(pStruct, i);
    Rotamers* pSetI = DesignSiteGetRotamers(pSiteI);
    Rotamer* pRotI = RotamersGet(pSetI, IntArrayGet(&pSeq->rotNdxs, i));
    Chain* pChainI = StructureGetChain(pStruct, pSiteI->chnNdx);
    RotamerRestore(pRotI, pSetI);
    for (int k = i; k < StructureGetDesignSiteNum(pStruct); k++) 
    {
      if (k == i) 
      {
        pSeq->ephy += pRotI->selfEnergy;
        pSeq->ebin += pRotI->selfEnergyBin;
      }
      else 
      {
        DesignSite* pSiteK = StructureGetDesignSite(pStruct, k);
        Rotamers* pSetK = DesignSiteGetRotamers(pSiteK);
        Rotamer* pRotK = RotamersGet(pSetK, IntArrayGet(&pSeq->rotNdxs, k));
        Chain* pChainK = StructureGetChain(pStruct, pSiteK->chnNdx);
        RotamerRestore(pRotK, pSetK);
        if (pSiteI->chnNdx == pSiteK->chnNdx) 
        { // from the same chain
          if (ChainGetType(pChainI) == Type_Chain_Protein)
          {
            EnergyRotamerAndRotamerSameChain(pRotI, pRotK, energyTermsExBind);
          }
        }
        else 
        { // from different chains
          if (ChainGetType(pChainI) == Type_Chain_SmallMol || ChainGetType(pChainK) == Type_Chain_SmallMol) 
          {
            if (pSiteI->pRes->desType == Type_DesType_Catalytic || pSiteK->pRes->desType == Type_DesType_Catalytic) 
            {
              CataConsSitePair* pCons1 = CataConsSitePairArrayFind(pConsArray, ResidueGetChainName(pSiteI->pRes), ResidueGetPosInChain(pSiteI->pRes), ResidueGetName(pSiteI->pRes),
                ResidueGetChainName(pSiteK->pRes), ResidueGetPosInChain(pSiteK->pRes), ResidueGetName(pSiteK->pRes));
              CataConsSitePair* pCons2 = CataConsSitePairArrayFind(pConsArray, ResidueGetChainName(pSiteK->pRes), ResidueGetPosInChain(pSiteK->pRes), ResidueGetName(pSiteK->pRes),
                ResidueGetChainName(pSiteI->pRes), ResidueGetPosInChain(pSiteI->pRes), ResidueGetName(pSiteI->pRes));
              if ((pCons1 != NULL && pCons1->pairConsType == Type_SitePairCons_Covalent) || (pCons2 != NULL && pCons2->pairConsType == Type_SitePairCons_Covalent))
              { // do not calculate energy between sites forming a covalent bond
                continue;
              }
            }
            double energyTermsPL[MAX_ENERGY_TERM] = { 0 };
            EnergyRotamerAndLigRotamer(pRotI, pRotK, energyTermsPL);
            if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) == NULL)
              || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) != NULL))
            {
              if (FLAG_DESIGN_PLI == TRUE || FLAG_DESIGN_ENZ == TRUE) 
              {
                for (int index = 0; index < MAX_ENERGY_TERM; index++) 
                {
                  energyTermsBind[index] += energyTermsPL[index];
                }
              }
            }
          }
          else
          {
            double energyTermsDC[MAX_ENERGY_TERM] = { 0 };
            EnergyRotamerAndRotamerDiffChain(pRotI, pRotK, energyTermsDC);
            if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) == NULL)
              || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) != NULL))
            {
              if (FLAG_DESIGN_PPI == TRUE)
              {
                for (int index = 0; index < MAX_ENERGY_TERM; index++)
                {
                  energyTermsBind[index] += energyTermsDC[index];
                }
              }
            }
          }
        }
        RotamerExtract(pRotK);
      }
    }
    RotamerExtract(pRotI);
  }
  EnergyTermWeighting(energyTermsExBind);
  pSeq->ephy += energyTermsExBind[0];
  EnergyTermWeighting(energyTermsBind);
  pSeq->ebin += energyTermsBind[0];

  // 2. compute evolutionary energy
  if (FLAG_EVOLUTION == TRUE) 
  {
    char seq[MAX_LEN_SEQ+1] = "";
    StructureGetWholeSequence(pSeq, pStruct, seq);
    if (FLAG_EVOPHIPSI == TRUE)
    {
      pSeq->eevo += EvolutionScoreAllFromSeq(seq);
    }
    else
    {
      pSeq->eevo += EvolutionScoreFromPSSMWithoutAlignment(seq);
    }
  }

  // 3. compute the number of unmet catalytic constriants
  if (FLAG_DESIGN_ENZ == TRUE)
  {
    for (int i = 0; i < CataConsSitePairArrayGetNum(pConsArray); i++)
    {
      CataConsSitePair* pCons = CataConsSitePairArrayGet(pConsArray, i);
      int siteIndex1 = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStruct, pCons->chnName1, pCons->pos1);
      int siteIndex2 = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStruct, pCons->chnName2, pCons->pos2);
      if (siteIndex1 == -1 || siteIndex2 == -1)
      {
        char errMsg[MAX_LEN_ERR_MSG + 1];
        sprintf(errMsg, "in file %s line %d, %s%d or %s%d has not been added as a catalytic design site", __FILE__, __LINE__,
          pCons->chnName1, pCons->pos1,
          pCons->chnName2, pCons->pos2);
        TraceError(errMsg, ValueError);
        continue;
      }
      DesignSite* pSiteI = StructureGetDesignSite(pStruct, siteIndex1);
      Rotamers* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotI = RotamersGet(pSetI, IntArrayGet(&pSeq->rotNdxs, siteIndex1));
      RotamerRestore(pRotI, pSetI);
      if (siteIndex1 == siteIndex2)
      {
        if (!CataConsSitePairCheck(pCons, pRotI, pRotI))
        {
          pSeq->unmetConsNum++;
        }
      }
      else
      {
        DesignSite* pSiteK = StructureGetDesignSite(pStruct, siteIndex2);
        Rotamers* pSetK = DesignSiteGetRotamers(pSiteK);
        Rotamer* pRotK = RotamersGet(pSetK, IntArrayGet(&pSeq->rotNdxs, siteIndex2));
        RotamerRestore(pRotK, pSetK);
        if (!CataConsSitePairCheck(pCons, pRotI, pRotK))
        {
          pSeq->unmetConsNum++;
        }
        RotamerExtract(pRotK);
      }
      RotamerExtract(pRotI);
    }
  }

  pSeq->etot = WGT_PROFILE * pSeq->eevo + pSeq->ephy + WGT_BIND * pSeq->ebin + WGT_CATA_CONS * pSeq->unmetConsNum;

  pSeq->mutNum = SequenceCalcNumOfMutations(pSeq, pStruct);
  if (FLAG_RESTRICT_NUM_OF_MUTATIONS == TRUE) 
  {
    pSeq->etot += EnergyPenaltyForOutOfRangeMutations(pSeq->mutNum);
  }
  pSeq->seqExist = SequenceIsSameAsSavedSequences(pSeq, pStruct, ppSavedSeqs, numSavedSeqs);
  if (FLAG_PENALIZE_IDENTICAL_SEQUENCES == TRUE) 
  {
    pSeq->etot += WGT_SAME_SEQ * pSeq->seqExist;
  }

  return Success;
}


int SequenceCalcEnergyChangeUponSingleMutation(Sequence* pSeq, Structure* pStruct, int mutSiteNdx, int mutRotNdx, double* dtot, double* dphy, double* dbin, double* devo, int *nmutsAfter, Sequence** ppSavedSeqs, int numSavedSeqs, BOOL *isSeqSameAsSaved)
{
  double phyBefore = 0, phyAfter = 0, binBefore = 0, binAfter = 0;
  // 1. compute physical energy
  DesignSite* pCurSite = StructureGetDesignSite(pStruct, mutSiteNdx);
  Rotamers* pCurSet = DesignSiteGetRotamers(pCurSite);
  Chain* pCurChain = StructureGetChain(pStruct, pCurSite->chnNdx);
  Rotamer* pCurRot = RotamersGet(pCurSet, IntArrayGet(&pSeq->rotNdxs, mutSiteNdx));
  Rotamer* pNewRot = RotamersGet(pCurSet, mutRotNdx);
  RotamerRestore(pCurRot, pCurSet);
  RotamerRestore(pNewRot, pCurSet);
  phyBefore += pCurRot->selfEnergy;
  phyAfter += pNewRot->selfEnergy;
  if (FLAG_DESIGN_MON == FALSE)
  {
    binBefore += pCurRot->selfEnergyBin;
    binAfter += pNewRot->selfEnergyBin;
  }
  double energyTermsExBind1[MAX_ENERGY_TERM] = { 0 };
  double energyTermsExBind2[MAX_ENERGY_TERM] = { 0 };
  double energyTermsBind1[MAX_ENERGY_TERM] = { 0 };
  double energyTermsBind2[MAX_ENERGY_TERM] = { 0 };
  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++)
  {
    if (i == mutSiteNdx)
    {
      continue;
    }
    DesignSite* pSiteI = StructureGetDesignSite(pStruct, i);
    Rotamers* pSetI = DesignSiteGetRotamers(pSiteI);
    Rotamer* pRotI = RotamersGet(pSetI, IntArrayGet(&pSeq->rotNdxs, i));
    Chain* pChainI = StructureGetChain(pStruct, pSiteI->chnNdx);
    RotamerRestore(pRotI, pSetI);
    if (pCurSite->chnNdx == pSiteI->chnNdx)
    {
      if (ChainGetType(pChainI) == Type_Chain_Protein)
      {
        EnergyRotamerAndRotamerSameChain(pRotI, pCurRot, energyTermsExBind1);
        EnergyRotamerAndRotamerSameChain(pRotI, pNewRot, energyTermsExBind2);
      }
    }
    else
    {
      if (ChainGetType(pChainI) == Type_Chain_SmallMol)
      {
        double energyTermsPL1[MAX_ENERGY_TERM] = { 0 };
        double energyTermsPL2[MAX_ENERGY_TERM] = { 0 };
        EnergyRotamerAndLigRotamer(pCurRot, pRotI, energyTermsPL1);
        EnergyRotamerAndLigRotamer(pNewRot, pRotI, energyTermsPL2);
        if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) == NULL)
          || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) != NULL))
        {
          if (FLAG_DESIGN_PLI == TRUE || FLAG_DESIGN_ENZ == TRUE)
          {
            for (int index = 0; index < MAX_ENERGY_TERM; index++)
            {
              energyTermsBind1[index] += energyTermsPL1[index];
              energyTermsBind2[index] += energyTermsPL2[index];
            }
          }
        }
      }
      else if (ChainGetType(pCurChain) == Type_Chain_SmallMol)
      {
        double energyTermsPL1[MAX_ENERGY_TERM] = { 0 };
        double energyTermsPL2[MAX_ENERGY_TERM] = { 0 };
        EnergyRotamerAndLigRotamer(pRotI, pCurRot, energyTermsPL1);
        EnergyRotamerAndLigRotamer(pRotI, pNewRot, energyTermsPL2);
        if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) == NULL)
          || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) != NULL))
        {
          if (FLAG_DESIGN_PLI == TRUE || FLAG_DESIGN_ENZ == TRUE)
          {
            for (int index = 0; index < MAX_ENERGY_TERM; index++)
            {
              energyTermsBind1[index] += energyTermsPL1[index];
              energyTermsBind2[index] += energyTermsPL2[index];
            }
          }
        }
      }
      else
      {
        double energyTermsDC1[MAX_ENERGY_TERM] = { 0 };
        double energyTermsDC2[MAX_ENERGY_TERM] = { 0 };
        EnergyRotamerAndRotamerDiffChain(pRotI, pCurRot, energyTermsDC1);
        EnergyRotamerAndRotamerDiffChain(pRotI, pNewRot, energyTermsDC2);
        if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) == NULL)
          || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) != NULL))
        {
          if (FLAG_DESIGN_PPI == TRUE)
          {
            for (int index = 0; index < MAX_ENERGY_TERM; index++)
            {
              energyTermsBind1[index] += energyTermsDC1[index];
              energyTermsBind2[index] += energyTermsDC2[index];
            }
          }
        }
      }
    }
    RotamerExtract(pRotI);
  }
  EnergyTermWeighting(energyTermsExBind1);
  EnergyTermWeighting(energyTermsExBind2);
  phyBefore += energyTermsExBind1[0];
  phyAfter += energyTermsExBind2[0];
  *dphy = phyAfter - phyBefore;

  if (FLAG_DESIGN_MON == FALSE)
  {
    EnergyTermWeighting(energyTermsBind1);
    EnergyTermWeighting(energyTermsBind2);
    binBefore += energyTermsBind1[0];
    binAfter += energyTermsBind2[0];
    *dbin = binAfter - binBefore;
  }
  RotamerExtract(pCurRot);
  RotamerExtract(pNewRot);

  // 2. compute evolutionary energy
  if (FLAG_EVOLUTION == TRUE)
  {
    double evoBefore = 0;
    double evoAfter = 0;
    DesignSite* pCurSite = StructureGetDesignSite(pStruct, mutSiteNdx);
    Rotamers* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRot = RotamersGet(pCurSet, IntArrayGet(&pSeq->rotNdxs, mutSiteNdx));
    Rotamer* pNewRot = RotamersGet(pCurSet, mutRotNdx);
    if (RotamerAndRotamerInSameType(pCurRot, pNewRot) == FALSE)
    {
      char seq1[MAX_LEN_SEQ];
      StructureGetWholeSequence(pSeq, pStruct, seq1);
      char seq2[MAX_LEN_SEQ];
      Sequence newSeq;
      SequenceCreate(&newSeq);
      SequenceCopy(&newSeq, pSeq);
      IntArraySet(&newSeq.rotNdxs, mutSiteNdx, mutRotNdx);
      StructureGetWholeSequence(&newSeq, pStruct, seq2);
      if (FLAG_EVOPHIPSI)
      {
        evoBefore += EvolutionScoreAllFromSeq(seq1);
        evoAfter += EvolutionScoreAllFromSeq(seq2);
      }
      else
      {
        evoBefore += EvolutionScoreFromPSSMWithoutAlignment(seq1);
        evoAfter += EvolutionScoreFromPSSMWithoutAlignment(seq2);
      }
      SequenceDestroy(&newSeq);
    }
    *devo = evoAfter - evoBefore;
  }

  *dtot = WGT_PROFILE * (*devo) + *dphy + WGT_BIND * (*dbin);

  *nmutsAfter = SequenceCalcNumOfMutationsUponSingleMutation(pSeq, pStruct, mutSiteNdx, mutRotNdx);
  if (FLAG_RESTRICT_NUM_OF_MUTATIONS == TRUE) {
    *dtot += EnergyPenaltyForOutOfRangeMutations(*nmutsAfter) - EnergyPenaltyForOutOfRangeMutations(pSeq->mutNum);
  }

  *isSeqSameAsSaved = SequenceIsSameAsSavedSequencesUponSingleMutation(pSeq, pStruct, mutSiteNdx, mutRotNdx, ppSavedSeqs, numSavedSeqs);
  if (FLAG_PENALIZE_IDENTICAL_SEQUENCES == TRUE) {
    *dtot += WGT_SAME_SEQ * (*isSeqSameAsSaved) - WGT_SAME_SEQ * pSeq->seqExist;
  }

  return Success;
}


int SequenceCalcEnergyChangeUponSingleMutationWithCataCons(Sequence* pSeq, Structure* pStruct, CataConsSitePairArray* pConsArray, int mutSiteNdx, int mutRotNdx, double* dtot, double* dphy, double* dbin, double* devo, int* dcons, int* nmutsAfter, Sequence** ppSavedSeqs, int numSavedSeqs, BOOL* isSeqSameAsSaved)
{
  double phyBefore = 0, phyAfter = 0, binBefore = 0, binAfter = 0, evoBefore = 0, evoAfter = 0, nConsBefore = 0, nConsAfter = 0;
  // 1. compute physical energy
  DesignSite* pCurSite = StructureGetDesignSite(pStruct, mutSiteNdx);
  Chain* pCurChain = StructureGetChain(pStruct, pCurSite->chnNdx);
  Rotamers* pCurSet = DesignSiteGetRotamers(pCurSite);
  Rotamer* pCurRot = RotamersGet(pCurSet, IntArrayGet(&pSeq->rotNdxs, mutSiteNdx));
  Rotamer* pNewRot = RotamersGet(pCurSet, mutRotNdx);
  RotamerRestore(pCurRot, pCurSet);
  RotamerRestore(pNewRot, pCurSet);
  phyBefore += pCurRot->selfEnergy;
  phyAfter += pNewRot->selfEnergy;
  binBefore += pCurRot->selfEnergyBin;
  binAfter += pNewRot->selfEnergyBin;
  double energyTermsExBind1[MAX_ENERGY_TERM] = { 0 };
  double energyTermsExBind2[MAX_ENERGY_TERM] = { 0 };
  double energyTermsBind1[MAX_ENERGY_TERM] = { 0 };
  double energyTermsBind2[MAX_ENERGY_TERM] = { 0 };

  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++)
  {
    if (i == mutSiteNdx)
    {
      continue;
    }
    DesignSite* pSiteI = StructureGetDesignSite(pStruct, i);
    Rotamers* pSetI = DesignSiteGetRotamers(pSiteI);
    Rotamer* pRotamerI = RotamersGet(pSetI, IntArrayGet(&pSeq->rotNdxs, i));
    Chain* pChainI = StructureGetChain(pStruct, pSiteI->chnNdx);
    RotamerRestore(pRotamerI, pSetI);
    if (pCurSite->chnNdx == pSiteI->chnNdx)
    {
      if (ChainGetType(pChainI) == Type_Chain_Protein)
      {
        EnergyRotamerAndRotamerSameChain(pRotamerI, pCurRot, energyTermsExBind1);
        EnergyRotamerAndRotamerSameChain(pRotamerI, pNewRot, energyTermsExBind2);
      }
    }
    else
    {
      if (ChainGetType(pChainI) == Type_Chain_SmallMol || ChainGetType(pCurChain) == Type_Chain_SmallMol)
      {
        double energyTermsPL1[MAX_ENERGY_TERM] = { 0 };
        double energyTermsPL2[MAX_ENERGY_TERM] = { 0 };
        EnergyRotamerAndLigRotamer(pRotamerI, pCurRot, energyTermsPL1);
        EnergyRotamerAndLigRotamer(pRotamerI, pNewRot, energyTermsPL2);
        if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) == NULL)
          || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) != NULL))
        {
          if (FLAG_DESIGN_PLI == TRUE || FLAG_DESIGN_ENZ == TRUE)
          {
            for (int index = 0; index < MAX_ENERGY_TERM; index++)
            {
              energyTermsBind1[index] += energyTermsPL1[index];
              energyTermsBind2[index] += energyTermsPL2[index];
            }
          }
        }
      }
      else
      {
        double energyTermsDC1[MAX_ENERGY_TERM] = { 0 };
        double energyTermsDC2[MAX_ENERGY_TERM] = { 0 };
        EnergyRotamerAndRotamerDiffChain(pRotamerI, pCurRot, energyTermsDC1);
        EnergyRotamerAndRotamerDiffChain(pRotamerI, pNewRot, energyTermsDC2);
        if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) == NULL)
          || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) != NULL))
        {
          if (FLAG_DESIGN_PPI == TRUE)
          {
            for (int index = 0; index < MAX_ENERGY_TERM; index++)
            {
              energyTermsBind1[index] += energyTermsDC1[index];
              energyTermsBind2[index] += energyTermsDC2[index];
            }
          }
        }
      }
    }
    RotamerExtract(pRotamerI);
  }
  EnergyTermWeighting(energyTermsExBind1);
  EnergyTermWeighting(energyTermsExBind2);
  EnergyTermWeighting(energyTermsBind1);
  EnergyTermWeighting(energyTermsBind2);
  phyBefore += energyTermsExBind1[0];
  phyAfter += energyTermsExBind2[0];
  *dphy = phyAfter - phyBefore;
  binBefore += energyTermsBind1[0];
  binAfter += energyTermsBind2[0];
  *dbin = binAfter - binBefore;
  RotamerExtract(pCurRot);
  RotamerExtract(pNewRot);

  // 2. compute evlutionary energy
  if (FLAG_EVOLUTION == TRUE)
  {
    DesignSite* pCurSite = StructureGetDesignSite(pStruct, mutSiteNdx);
    Rotamers* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRot = RotamersGet(pCurSet, IntArrayGet(&pSeq->rotNdxs, mutSiteNdx));
    Rotamer* pNewRot = RotamersGet(pCurSet, mutRotNdx);
    if (RotamerAndRotamerInSameType(pCurRot, pNewRot) == FALSE)
    {
      char seq1[MAX_LEN_SEQ];
      StructureGetWholeSequence(pSeq, pStruct, seq1);
      char seq2[MAX_LEN_SEQ];
      Sequence newSeq;
      SequenceCreate(&newSeq);
      SequenceCopy(&newSeq, pSeq);
      IntArraySet(&newSeq.rotNdxs, mutSiteNdx, mutRotNdx);
      StructureGetWholeSequence(&newSeq, pStruct, seq2);
      if (FLAG_EVOPHIPSI == TRUE)
      {
        evoBefore += EvolutionScoreAllFromSeq(seq1);
        evoAfter += EvolutionScoreAllFromSeq(seq2);
      }
      else
      {
        evoBefore += EvolutionScoreFromPSSMWithoutAlignment(seq1);
        evoAfter += EvolutionScoreFromPSSMWithoutAlignment(seq2);
      }
      SequenceDestroy(&newSeq);
    }
    *devo = evoAfter - evoBefore;
  }

  // 3. compute the number of unmet catalytic constraints
  if (FLAG_DESIGN_ENZ == TRUE)
  {
    DesignSite* pCurSite = StructureGetDesignSite(pStruct, mutSiteNdx);
    if (ResidueGetDesignType(pCurSite->pRes) == Type_DesType_Catalytic || ResidueGetDesignType(pCurSite->pRes) == Type_DesType_SmallMol)
    {
      Rotamers* pCurSet = DesignSiteGetRotamers(pCurSite);
      Rotamer* pCurRot = RotamersGet(pCurSet, IntArrayGet(&pSeq->rotNdxs, mutSiteNdx));
      Rotamer* pNewRot = RotamersGet(pCurSet, mutRotNdx);
      if (pNewRot != pCurRot)
      {
        RotamerRestore(pCurRot, pCurSet);
        RotamerRestore(pNewRot, pCurSet);
        for (int i = 0; i < CataConsSitePairArrayGetNum(pConsArray); i++)
        {
          CataConsSitePair* pCons = CataConsSitePairArrayGet(pConsArray, i);
          int siteIndex1 = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStruct, pCons->chnName1, pCons->pos1);
          int siteIndex2 = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStruct, pCons->chnName2, pCons->pos2);
          if (siteIndex1 == -1 || siteIndex2 == -1)
          {
            continue;
          }
          if (siteIndex1 == mutSiteNdx && siteIndex2 != mutSiteNdx)
          {
            DesignSite* pSiteK = StructureGetDesignSite(pStruct, siteIndex2);
            Rotamers* pSetK = DesignSiteGetRotamers(pSiteK);
            Rotamer* pRotK = RotamersGet(pSetK, IntArrayGet(&pSeq->rotNdxs, siteIndex2));
            RotamerRestore(pRotK, pSetK);
            if (CataConsSitePairCheck(pCons, pCurRot, pRotK) == FALSE)
            {
              nConsBefore++;
            }
            if (CataConsSitePairCheck(pCons, pNewRot, pRotK) == FALSE)
            {
              nConsAfter++;
            }
            RotamerExtract(pRotK);
          }
          else if (siteIndex1 != mutSiteNdx && siteIndex2 == mutSiteNdx)
          {
            DesignSite* pSiteK = StructureGetDesignSite(pStruct, siteIndex1);
            Rotamers* pSetK = DesignSiteGetRotamers(pSiteK);
            Rotamer* pRotK = RotamersGet(pSetK, IntArrayGet(&pSeq->rotNdxs, siteIndex1));
            RotamerRestore(pRotK, pSetK);
            if (CataConsSitePairCheck(pCons, pRotK, pCurRot) == FALSE)
            {
              nConsBefore++;
            }
            if (CataConsSitePairCheck(pCons, pRotK, pNewRot) == FALSE)
            {
              nConsAfter++;
            }
            RotamerExtract(pRotK);
          }
          else if (siteIndex1 == mutSiteNdx && siteIndex2 == mutSiteNdx)
          {
            if (CataConsSitePairCheck(pCons, pCurRot, pCurRot) == FALSE)
            {
              nConsBefore++;
            }
            if (CataConsSitePairCheck(pCons, pNewRot, pNewRot) == FALSE)
            {
              nConsAfter++;
            }
          }
        }
        RotamerExtract(pNewRot);
        RotamerExtract(pCurRot);
      }
    }
    *dcons = nConsAfter - nConsBefore;
  }

  *dtot = WGT_PROFILE * (*devo) + *dphy + WGT_BIND * (*dbin) + WGT_CATA_CONS * (*dcons);

  *nmutsAfter = SequenceCalcNumOfMutationsUponSingleMutation(pSeq, pStruct, mutSiteNdx, mutRotNdx);
  if (FLAG_RESTRICT_NUM_OF_MUTATIONS == TRUE)
  {
    *dtot += EnergyPenaltyForOutOfRangeMutations(*nmutsAfter) - EnergyPenaltyForOutOfRangeMutations(pSeq->mutNum);
  }

  *isSeqSameAsSaved = SequenceIsSameAsSavedSequencesUponSingleMutation(pSeq, pStruct, mutSiteNdx, mutRotNdx, ppSavedSeqs, numSavedSeqs);
  if (FLAG_PENALIZE_IDENTICAL_SEQUENCES == TRUE)
  {
    *dtot += WGT_SAME_SEQ * (*isSeqSameAsSaved) - WGT_SAME_SEQ * pSeq->seqExist;
  }

  return Success;
}


int Metropolis(Sequence* pSeq, Sequence* pBest, Structure* pStruct, RotamerList* pList, StringArray** ppRotTypes, IntArray** ppRotNums, double T, int stepNum, Sequence** ppSavedSeqs, int numSavedSeqs)
{
  int nacc = 0;
  for (int i = 0; i < stepNum; i++) 
  {
    int   mutSiteNdx;
    int   mutRotNdx;
    SequenceRandomSiteIndex(&mutSiteNdx, StructureGetDesignSiteNum(pStruct));
    SequenceRandRotamerIndex(pStruct, pList, ppRotTypes, ppRotNums, mutSiteNdx, &mutRotNdx);
    double dtot = 0, dphy = 0, devo = 0, dbin = 0;
    int nmutsAfter = 0;
    BOOL isSeqSameAsSaved = FALSE;
    SequenceCalcEnergyChangeUponSingleMutation(pSeq, pStruct, mutSiteNdx, mutRotNdx, &dtot, &dphy, &dbin, &devo, &nmutsAfter, ppSavedSeqs, numSavedSeqs, &isSeqSameAsSaved);
    if (exp(-1.0 * dtot / T) > (double)(rand() + 1.0) / (RAND_MAX + 1.0)) 
    {
      SequenceUpdateSingleSite(pSeq, mutSiteNdx, mutRotNdx);
      pSeq->mutNum = nmutsAfter;
      pSeq->seqExist = isSeqSameAsSaved;
      pSeq->etot += dtot;
      pSeq->ephy += dphy;
      pSeq->ebin += dbin;
      pSeq->eevo += devo;
      if (pSeq->etot < pBest->etot) SequenceCopy(pBest, pSeq);
      nacc++;
    }
  }
  float racc = (float)nacc / stepNum;
  printf("acceptance rate = %9.6f at T = %10.6f\n", racc, T);

  return Success;
}


int MetropolisWithCataCons(Sequence* pSeq, Sequence* pBest, Structure* pStruct, CataConsSitePairArray* pConsArray, RotamerList* pList, StringArray** ppRotTypes, IntArray** ppRotNums, double T, int stepNum, Sequence** ppSavedSeqs, int numSavedSeqs)
{
  int nacc = 0;
  for (int i = 0; i < stepNum; i++) 
  {
    int mutSiteNdx;
    int mutRotNdx;
    SequenceRandomSiteIndex(&mutSiteNdx, StructureGetDesignSiteNum(pStruct));
    SequenceRandRotamerIndex(pStruct, pList, ppRotTypes, ppRotNums, mutSiteNdx, &mutRotNdx);
    double dtot = 0, dphy = 0, devo = 0, dbin = 0;
    int dcons = 0, nmutsAfter = 0;
    BOOL isSeqSameAsSaved = FALSE;
    SequenceCalcEnergyChangeUponSingleMutationWithCataCons(pSeq, pStruct, pConsArray, mutSiteNdx, mutRotNdx, &dtot, &dphy, &dbin, &devo, &dcons, &nmutsAfter, ppSavedSeqs, numSavedSeqs, &isSeqSameAsSaved);
    if (exp(-1.0 * dtot / T) > (double)(rand() + 1.0) / (RAND_MAX + 1.0)) 
    {
      SequenceUpdateSingleSite(pSeq, mutSiteNdx, mutRotNdx);
      pSeq->mutNum = nmutsAfter;
      pSeq->seqExist = isSeqSameAsSaved;
      pSeq->etot += dtot;
      pSeq->ephy += dphy;
      pSeq->ebin += dbin;
      pSeq->eevo += devo;
      pSeq->unmetConsNum += dcons;
      if (pSeq->etot < pBest->etot)
      {
        SequenceCopy(pBest, pSeq);
      }
      nacc++;
    }
  }
  float racc = (float)nacc / stepNum;
  printf("acceptance rate = %9.6f at T = %10.6f\n", racc, T);

  return Success;
}


int SimulatedAnnealing(Structure* pStruct, RotamerList* pList)
{
  int result = Success;
  char errMsg[MAX_LEN_ERR_MSG+1];
  srand((unsigned int)time(NULL));
  StringArray* pRotTypes = (StringArray*)malloc(sizeof(StringArray) * pList->desSiteNum);
  IntArray* pRotNums = (IntArray*)malloc(sizeof(IntArray) * pList->desSiteNum);
  DesignSiteShowRotamerTypeAndNum(pList, pStruct, &pRotTypes, &pRotNums);
  int remainRotNum = 0;
  for (int i = 0; i < pList->desSiteNum; i++) 
  {
    remainRotNum += pList->remainRotNum[i];
  }
  printf("number of total remaining rotamers: %d\n", remainRotNum);
  remainRotNum = remainRotNum < NUM_METROPOLIS_STEPS ? remainRotNum : NUM_METROPOLIS_STEPS;
  printf("set the number of monte carlo moves at each temperature to %d\n", remainRotNum);

  CataConsSitePairArray consArray;
  if (FLAG_DESIGN_ENZ == TRUE) 
  {
    result = CataConsSitePairArrayCreate(&consArray, FILE_LIG_CATACONS);
    if (FAILED(result)) 
    {
      sprintf(errMsg, "in file %s line %d, failed to create catalytic constraints from file %s",
        __FILE__, __LINE__, FILE_LIG_CATACONS);
      result = ValueError;
      return result;
    }
    for (int i = 0;i < CataConsSitePairArrayGetNum(&consArray);i++) 
    {
      result = StructureDeployCataConsSitePair(pStruct, CataConsSitePairArrayGet(&consArray, i));
      if (FAILED(result)) 
      {
        sprintf(errMsg, "in file %s line %d, failed to deploy catalytic constraints",
          __FILE__, __LINE__);
        result = ValueError;
        return result;
      }
    }
  }

  printf("search low-energy sequences by monte-carlo simulated annealing ...\n");
  int ntrajs = NTRAJ - NTRAJ_START_NDX + 1;
  Sequence* allBestSeqs = (Sequence*)malloc(sizeof(Sequence)*ntrajs);
  for (int i = 0; i < ntrajs; i++) 
  {
    SequenceCreate(&allBestSeqs[i]);
  }
  // write saved lowest-energy sequences from all trajectories to file
  char FILE_BEST_DESIGNS[MAX_LEN_FILE_NAME+1], FILE_BEST_SEQS[MAX_LEN_FILE_NAME+1];
  sprintf(FILE_BEST_DESIGNS, "%s", TAG_BEST_DESIGNS);
  sprintf(FILE_BEST_SEQS, "%s", TAG_BEST_SEQS);
  FILE* pFileDesigns = NULL, *pFileSeqs = NULL;
  if (NTRAJ_START_NDX > 1) 
  {
    pFileDesigns = fopen(FILE_BEST_DESIGNS, "a");
    pFileSeqs = fopen(FILE_BEST_SEQS, "a");
  }
  else 
  {
    double wgt1 = 0.0, wgt2 = 0.0;
    if (FLAG_RESTRICT_NUM_OF_MUTATIONS == TRUE)
    {
      wgt1 = WGT_MUT_NUM;
    }
    if (FLAG_PENALIZE_IDENTICAL_SEQUENCES == TRUE) 
    {
      wgt2 = WGT_SAME_SEQ;
    }
    pFileDesigns = fopen(FILE_BEST_DESIGNS, "w");
    fprintf(pFileDesigns,
      "Model_pdb\t"
      "Sequence_recovery_rate\t"
      "Weighted_total_energy\t"
      "Raw_evolution_energy|Weight=%.1f\t"
      "Raw_stability_energy_without_binding|Weight=1.0\t"
      "Raw_binding_energy|Weight=%.1f\t"
      "Num_unmet_catalytic_constraints|Penalty for enzyme design|Weight=%.1f\t"
      "Num_mutations|Penalty if num of mutations is not in specified %d-%d|Weight=%.1f\t"
      "Sequence_already_saved|Penalty for already saved seq|Weight=%.1f\t"
      "Mutant_string\n"
      , WGT_PROFILE, WGT_BIND, WGT_CATA_CONS, MIN_MUT_NUM, MAX_MUT_NUM, wgt1, wgt2
    );
    fflush(pFileDesigns);

    pFileSeqs = fopen(FILE_BEST_SEQS, "w");
    fprintf(pFileSeqs,
      "Model_pdb\t"
      "Full_sequence\t"
      "Design_site_sequence\n"
    );
    fflush(pFileSeqs);
  }
  
  // Generate multiple independent design trajectories
  int nSavedTraj = 0;
  int traj = NTRAJ_START_NDX;
  int lastSavedTraj = traj;
  for (; traj <= NTRAJ; traj++) 
  {
    if(traj - lastSavedTraj > 10)
    {
      break;
    }
    printf("=== design trajectory %d ===\n", traj);
    Sequence oldSeq, bestSeq;
    SequenceCreate(&oldSeq);
    SequenceCreate(&bestSeq);
    if (FLAG_DESIGN_FROM_NATAA == TRUE)
    {
      SequenceGenNativeSeqSeed(&oldSeq, pStruct, pList);
    }
    else
    { 
      SequenceGenRandomSeed(&oldSeq, pList);
    }
    if (FLAG_DESIGN_ENZ == TRUE)
    {
      SequenceCalcEnergyWithCataCons(&oldSeq, pStruct, &consArray, &allBestSeqs, nSavedTraj);
    }
    else
    {
      SequenceCalcEnergy(&oldSeq, pStruct, &allBestSeqs, nSavedTraj);
    }
    SequenceCopy(&bestSeq, &oldSeq);

    for (int cycle = 1; cycle <= SA_CYCLE; cycle++) 
    {
      printf("simulated annealing cycle %d\n", cycle);
      double T = SA_T_INITIAL / cycle;
      while (T > SA_T_FINAL) 
      {
        if (FLAG_DESIGN_ENZ == TRUE) 
        {
          MetropolisWithCataCons(&oldSeq, &bestSeq, pStruct, &consArray, pList, &pRotTypes, &pRotNums, T, remainRotNum, &allBestSeqs, nSavedTraj);
        }
        else 
        { 
          Metropolis(&oldSeq, &bestSeq, pStruct, pList, &pRotTypes, &pRotNums, T, remainRotNum, &allBestSeqs, nSavedTraj);
        }
        T *= SA_DECREASE_FACTOR;
      }
    }

    if(bestSeq.mutNum <= MAX_MUT_NUM && bestSeq.mutNum >= MIN_MUT_NUM && bestSeq.seqExist == FALSE)
    {
      // write the lowest-energy structure to file
      char FILE_BEST_MODEL_TRAJ[MAX_LEN_FILE_NAME+1];
      sprintf(FILE_BEST_MODEL_TRAJ, "%s%04d.pdb", TAG_BEST_MODEL, traj);
      // remember to include this before save the bestSeq
      strcpy(bestSeq.pdbToSaveSeq, FILE_BEST_MODEL_TRAJ);
      SequenceCopy(&allBestSeqs[nSavedTraj], &bestSeq);
      SequenceWriteStructureInPDB(&bestSeq, pStruct, FILE_BEST_MODEL_TRAJ);

      // append the lowest-energy sequence (bestSeq) to FILE stream
      SequenceWriteDesignSummaryToFile(&bestSeq, pStruct, traj, pFileDesigns);
      SequenceWriteDesignSequenceToFile(&bestSeq, pStruct, traj, pFileSeqs);
      fflush(pFileDesigns);
      fflush(pFileSeqs);

      // write the lowest-energy-sequence design-sites (ALL) and corresponding lig pose to files
      if (FLAG_MOL2 == TRUE) {
        char FILE_BEST_LIG_MOL2_TRAJ[MAX_LEN_FILE_NAME+1];
        sprintf(FILE_BEST_LIG_MOL2_TRAJ, "%s%04d.mol2", TAG_BEST_LIG_MOL2, traj);
        SequenceWriteLigandInMol2(&bestSeq, pStruct, FILE_BEST_LIG_MOL2_TRAJ);
      }
      SequenceDestroy(&oldSeq);
      SequenceDestroy(&bestSeq);
      nSavedTraj++;
      lastSavedTraj = traj;
    }
    else if(bestSeq.mutNum > MAX_MUT_NUM || bestSeq.mutNum < MIN_MUT_NUM)
    {
      SequenceDestroy(&oldSeq);
      SequenceDestroy(&bestSeq);
      continue;
    }
    else if(bestSeq.seqExist == TRUE)
    {
      SequenceDestroy(&oldSeq);
      SequenceDestroy(&bestSeq);
      break;
    }
  }
  fclose(pFileDesigns);
  fclose(pFileSeqs);

  // release memory
  for (int i = 0; i < ntrajs; i++)
  {
    SequenceDestroy(&allBestSeqs[i]);
  }
  free(allBestSeqs);
  allBestSeqs = NULL;
  
  if (FLAG_DESIGN_ENZ == TRUE)
  {
    CataConsSitePairArrayDestroy(&consArray);
  }
  for (int i = 0; i < pList->desSiteNum; i++) 
  {
    StringArrayDestroy(&pRotTypes[i]);
    IntArrayDestroy(&pRotNums[i]);
  }
  free(pRotTypes);
  free(pRotNums);

  return Success;
}

int SequenceCalcNumOfMutations(Sequence* pSeq, Structure* pStruct)
{
  int nmuts = 0;
  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++) 
  {
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    int rotNdx = IntArrayGet(&pSeq->rotNdxs, i);
    Rotamer* pRot = RotamersGet(pSet, rotNdx);
    if (RotamerAndResidueInSameType(pRot, pSite->pRes) == FALSE)
    {
      nmuts++;
    }
  }

  return nmuts;
}

int SequenceCalcNumOfMutationsUponSingleMutation(Sequence* pSeq, Structure* pStruct, int mutSiteNdx, int mutRotNdx)
{
  int nmuts = 0;
  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++) {
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    Rotamer* pRot = NULL;
    if (i != mutSiteNdx) {
      int rotNdx = IntArrayGet(&pSeq->rotNdxs, i);
      pRot = RotamersGet(pSet, rotNdx);
    }
    else pRot = RotamersGet(pSet, mutRotNdx);
    if (RotamerAndResidueInSameType(pRot, pSite->pRes) == FALSE) nmuts++;
  }

  return nmuts;
}


double EnergyPenaltyForOutOfRangeMutations(int count)
{
  if (count > MAX_MUT_NUM)
  {
    return WGT_MUT_NUM * (count - MAX_MUT_NUM);
  }
  else if (count < MIN_MUT_NUM)
  {
    return WGT_MUT_NUM * (MIN_MUT_NUM - count);
  }
  return 0.0;
}


BOOL SequenceIsSameAsSavedSequences(Sequence* pSeq, Structure* pStruct, Sequence** ppSavedSeqs, int numSavedSeqs)
{
  for (int j = 0; j < numSavedSeqs; j++) 
  {
    int ndiff = 0;
    Sequence* pSavedSeq = (*ppSavedSeqs+j);
    for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++) 
    {
      DesignSite* pSite = StructureGetDesignSite(pStruct, i);
      Rotamers* pSet = DesignSiteGetRotamers(pSite);
      int rotNdx = IntArrayGet(&pSeq->rotNdxs, i);
      Rotamer* pRot = RotamersGet(pSet, rotNdx);
      int rotNdx2 = IntArrayGet(&pSavedSeq->rotNdxs, i);
      Rotamer* pRot2 = RotamersGet(pSet, rotNdx2);
      if (RotamerAndRotamerInSameType(pRot, pRot2) == FALSE)
      {
        ndiff++;
        break;
      }
    }
    if (ndiff == 0)
    {
      return TRUE;
    }
  }
  return FALSE;
}


BOOL SequenceIsSameAsSavedSequencesUponSingleMutation(Sequence* pSeq, Structure* pStruct, int mutSiteNdx, int mutRotNdx, Sequence** ppSavedSeqs, int numSavedSeqs)
{
  for (int j = 0; j < numSavedSeqs; j++) 
  {
    int ndiff = 0;
    Sequence* pSavedSeq = (*ppSavedSeqs+j);
    for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++) 
    {
      DesignSite* pSite = StructureGetDesignSite(pStruct, i);
      Rotamers* pSet = DesignSiteGetRotamers(pSite);
      int rotNdx2 = IntArrayGet(&pSavedSeq->rotNdxs, i);
      Rotamer* pRot2 = RotamersGet(pSet, rotNdx2);
      Rotamer* pRot = NULL;      
      if (i != mutSiteNdx) 
      {
        int rotNdx = IntArrayGet(&pSeq->rotNdxs, i);
        pRot = RotamersGet(pSet, rotNdx);
      }
      else pRot = RotamersGet(pSet, mutRotNdx);
      if (RotamerAndRotamerInSameType(pRot, pRot2) == FALSE)
      {
        ndiff++;
        break;
      }
    }
    if (ndiff == 0)
    {
      return TRUE;
    }
  }
  return FALSE;
}