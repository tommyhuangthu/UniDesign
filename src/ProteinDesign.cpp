/*******************************************************************************************************************************
Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com)

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
#include "RotamerBuilder.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern BOOL FLAG_EVOLUTION;
extern BOOL FLAG_PHYSICS;
extern BOOL FLAG_MONOMER;
extern BOOL FLAG_PPI;
extern BOOL FLAG_PROT_LIG;
extern BOOL FLAG_ENZYME;
extern BOOL FLAG_MOL2;
extern BOOL FLAG_EVOPHIPSI;
extern BOOL FLAG_DESIGN_FROM_NATAA;

extern char BEST_SEQ[MAX_LEN_FILE_NAME + 1];
extern char BEST_STRUCT[MAX_LEN_FILE_NAME + 1];
extern char BEST_DESSITE[MAX_LEN_FILE_NAME + 1];
extern char BEST_MOL2[MAX_LEN_FILE_NAME + 1];
extern char ROT_INDEX_FILE[MAX_LEN_FILE_NAME + 1];
extern char SEQ_FILE[MAX_LEN_FILE_NAME + 1];
extern char FILE_CATACONS[MAX_LEN_FILE_NAME + 1];

extern int PROT_LEN_NORM;
extern int NTRAJ;
extern int NTRAJ_START_NDX;

extern char PDBID[MAX_LEN_FILE_NAME + 1];
extern char DES_CHAINS[MAX_LEN_ONE_LINE_CONTENT + 1];
extern double WGT_PROFILE;
extern float** PROT_PROFILE;

extern double WGT_BIND;


extern char PROGRAM_PATH[MAX_LEN_ONE_LINE_CONTENT + 1];
extern char TGT_PRF[MAX_LEN_FILE_NAME + 1];
extern char TGT_SA[MAX_LEN_FILE_NAME + 1];
extern char TGT_SS[MAX_LEN_FILE_NAME + 1];
extern char TGT_SEQ[MAX_LEN_FILE_NAME + 1];
extern char TGT_PHIPSI[MAX_LEN_FILE_NAME + 1];

#define WGT_CATA_CONS  1000.0


int SitePairConsDeploy(CataConsSitePairArray* pSitePairArray, Structure* pStructure)
{
  for (int i = 0; i < CataConsSitePairArrayGetCount(pSitePairArray); i++)
  {
    CataConsSitePair* pSitePair = CataConsSitePairArrayGet(pSitePairArray, i);
    DesignSite* pFirstDesignSite = StructureFindDesignSiteByChainName(pStructure, pSitePair->chnName1, pSitePair->pos1);
    DesignSite* pSecondDesignSite = StructureFindDesignSiteByChainName(pStructure, pSitePair->chnName2, pSitePair->pos2);
    CataConsSitePairDeploy(pSitePair, DesignSiteGetRotamers(pFirstDesignSite), DesignSiteGetRotamers(pSecondDesignSite));
  }
  return Success;
}


int EnergyMatrixUpdateForCataCons(EnergyMatrix* pMatrix, RotamerList* pList, CataConsSitePair* pSitePair, Structure* pStructure)
{
  DesignSite* pFirstDesignSite = StructureFindDesignSiteByChainName(pStructure, pSitePair->chnName1, pSitePair->pos1);
  DesignSite* pSecondDesignSite = StructureFindDesignSiteByChainName(pStructure, pSitePair->chnName2, pSitePair->pos2);
  int firstDesignSiteIndex = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStructure, pSitePair->chnName1, pSitePair->pos1);
  int secondDesignSiteIndex = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStructure, pSitePair->chnName2, pSitePair->pos2);
  EnergyMatrixBlock* pBlockII = EnergyMatrixGetBlock(pMatrix, firstDesignSiteIndex, firstDesignSiteIndex);
  EnergyMatrixBlock* pBlockKK = EnergyMatrixGetBlock(pMatrix, secondDesignSiteIndex, secondDesignSiteIndex);
  for (int j = 0; j < pBlockII->RotamerCountSiteI; j++)
  {
    int trueIndexIJ;
    RotamerOriginalIndexGet(pList, firstDesignSiteIndex, j, &trueIndexIJ);
    Rotamer* pRotamerJ = RotamerSetGet(DesignSiteGetRotamers(pFirstDesignSite), trueIndexIJ);
    for (int s = 0; s < pBlockKK->RotamerCountSiteK; s++)
    {
      int trueIndexKS;
      RotamerOriginalIndexGet(pList, secondDesignSiteIndex, s, &trueIndexKS);
      Rotamer* pRotamerS = RotamerSetGet(DesignSiteGetRotamers(pSecondDesignSite), trueIndexKS);
      if (CataConsSitePairCheck(pSitePair, pRotamerJ, pRotamerS) == FALSE)
      {
        if (firstDesignSiteIndex < secondDesignSiteIndex)
        {
          *EnergyMatrixGet(pMatrix, firstDesignSiteIndex, secondDesignSiteIndex, j, s) += 1e4;
        }
        else
        {
          *EnergyMatrixGet(pMatrix, secondDesignSiteIndex, firstDesignSiteIndex, s, j) += 1e4;
        }
      }
    }
  }
  return Success;
}


int EnergyMatrixUpdateForCataConsArray(EnergyMatrix* pMatrix, RotamerList* pList, CataConsSitePairArray* pSitePairArray, Structure* pStructure)
{
  for (int i = 0; i < CataConsSitePairArrayGetCount(pSitePairArray); i++)
  {
    CataConsSitePair* pSitePair = CataConsSitePairArrayGet(pSitePairArray, i);
    EnergyMatrixUpdateForCataCons(pMatrix, pList, pSitePair, pStructure);
  }
  return Success;
}


int DesignSiteShowRotamerTypeAndCount(RotamerList* pList, Structure* pStructure, StringArray** ppRotamerType, IntArray** ppRotamerCount)
{
  StringArray* pStringArray = *ppRotamerType;
  IntArray* pIntArray = *ppRotamerCount;

  printf("show rotamer distribution at each design site\n");
  for (int i = 0; i < StructureGetDesignSiteCount(pStructure); i++)
  {
    int remainRotamerCount = 0;
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure, i);
    RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    IntArrayCreate(&pIntArray[i], 0);
    StringArrayCreate(&pStringArray[i]);
    for (int j = 0; j < pList->rotamerCount[i]; j++)
    {
      if (pList->remainFlag[i][j] == FALSE)
      {
        continue;
      }
      BOOL     isRotamerTypeNew = TRUE;
      Rotamer* pRotamer = RotamerSetGet(pRotamerSet, j);
      for (int k = 0; k < StringArrayGetCount(&pStringArray[i]); k++)
      {
        if (strcmp(RotamerGetType(pRotamer), StringArrayGet(&pStringArray[i], k)) == 0)
        {
          isRotamerTypeNew = FALSE;
          break;
        }
      }
      if (isRotamerTypeNew == TRUE)
      {
        int length = IntArrayGetLength(&pIntArray[i]);
        StringArrayAppend(&pStringArray[i], RotamerGetType(pRotamer));
        length++;
        IntArrayResize(&pIntArray[i], length);
        IntArraySet(&pIntArray[i], length - 1, 1);
      }
      else
      {
        int pos, rotamerCount;
        StringArrayFind(&pStringArray[i], RotamerGetType(pRotamer), &pos);
        rotamerCount = IntArrayGet(&pIntArray[i], pos);
        rotamerCount++;
        IntArraySet(&pIntArray[i], pos, rotamerCount);
      }
    }
    for (int j = 0; j < IntArrayGetLength(&pIntArray[i]); j++)
    {
      remainRotamerCount += IntArrayGet(&pIntArray[i], j);
    }
    printf("site %3d : %3s %s %4d, %6d rotamers:  ", i,
      ResidueGetName(pDesignSite->pRes),
      ResidueGetChainName(pDesignSite->pRes),
      ResidueGetPosInChain(pDesignSite->pRes),
      remainRotamerCount);
    for (int j = 0; j < StringArrayGetCount(&pStringArray[i]); j++)
    {
      printf("%4d %s ", IntArrayGet(&pIntArray[i], j), StringArrayGet(&pStringArray[i], j));
    }
    printf("\n");
  }

  return Success;
}


int SequenceRandomSiteIndex(int* mutSiteIndex, int designSiteCount)
{
  *mutSiteIndex = rand() % designSiteCount;
  return Success;
}


int SequenceRandRotamerIndex(Structure* pStruct, RotamerList* pList, StringArray** ppRotTypes, IntArray** ppRotCounts, int siteNdx, int* rotNdx)
{
  StringArray* pStringArray = *ppRotTypes;
  IntArray* pIntArray = *ppRotCounts;
  int typeIndex = rand() % IntArrayGetLength(&pIntArray[siteNdx]);
  int indexInType = rand() % IntArrayGet(&pIntArray[siteNdx], typeIndex);
  int typeflag = 0;
  for (int i = 0; i < pList->rotamerCount[siteNdx]; i++)
  {
    Rotamer* pRotamer = RotamerSetGet(DesignSiteGetRotamers(StructureGetDesignSite(pStruct, siteNdx)), i);
    if (pList->remainFlag[siteNdx][i] == FALSE)
    {
      continue;
    }
    if (strcmp(StringArrayGet(&pStringArray[siteNdx], typeIndex), RotamerGetType(pRotamer)) != 0)
    {
      continue;
    }
    if (typeflag == indexInType)
    {
      *rotNdx = i;
      break;
    }
    typeflag++;
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
  pSeq->desSiteCount = pList->desSiteCount;
  pSeq->etot = 0;
  pSeq->ephy = 0;
  pSeq->eevo = 0;
  IntArrayResize(&pSeq->rotNdxs, pSeq->desSiteCount);
  for (int i = 0; i < pSeq->desSiteCount; i++)
  {
    int j = rand() % pList->remainRotamerCount[i];
    int trueJ;
    RotamerOriginalIndexGet(pList, i, j, &trueJ);
    IntArraySet(&pSeq->rotNdxs, i, trueJ);
  }

  return Success;
}


int SequenceGenInitialSeqSeed(Structure* pStruct, RotamerList* pList, Sequence* pSeq)
{
  SequenceDestroy(pSeq);
  SequenceCreate(pSeq);
  pSeq->desSiteCount = pList->desSiteCount;
  pSeq->etot = 0.0;
  IntArrayResize(&pSeq->rotNdxs, pSeq->desSiteCount);
  for (int i = 0; i < pSeq->desSiteCount; i++)
  {
    double minSelfEnergy = 1e8;
    int minSelfEnergyIndex = -1;
    DesignSite* pDesignSite = StructureGetDesignSite(pStruct, i);
    RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    for (int j = 0;j < pList->rotamerCount[i];j++)
    {
      if (pList->remainFlag[i][j] == FALSE)
      {
        continue;
      }
      Rotamer* pRotamer = RotamerSetGet(pRotamerSet, j);
      if (pRotamer->selfEnergy < minSelfEnergy)
      {
        minSelfEnergy = pRotamer->selfEnergy;
        minSelfEnergyIndex = j;
      }
    }
    IntArraySet(&pSeq->rotNdxs, i, minSelfEnergyIndex);
  }

  return Success;
}


int SequenceGenNativeSeqSeed(Structure* pStruct, RotamerList* pList, Sequence* pSeq)
{
  SequenceDestroy(pSeq);
  SequenceCreate(pSeq);
  pSeq->desSiteCount = pList->desSiteCount;
  pSeq->etot = 0.0;
  IntArrayResize(&pSeq->rotNdxs, pSeq->desSiteCount);
  for (int i = 0; i < pSeq->desSiteCount; i++)
  {
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    RotamerSet* pSet = DesignSiteGetRotamers(pSite);
    for (int j = 0;j < RotamerSetGetCount(pSet);j++)
    {
      Rotamer* pRotamer = RotamerSetGet(pSet, j);
      RotamerRestore(pRotamer, pSet);
      if (RotamerAndResidueSidechainRMSD(pRotamer, pSite->pRes) < 0.05)
      {
        IntArraySet(&pSeq->rotNdxs, i, j);
        break;
      }
      RotamerExtract(pRotamer);
    }
  }

  return Success;
}


int SequenceTemplateEnergy(Structure* pStruct, Sequence* pSeq, double energyTerms[MAX_ENERGY_TERM], double energyTermsBind[MAX_ENERGY_TERM])
{
  for (int i = 0; i < StructureGetChainCount(pStruct); i++)
  {
    Chain* pChainI = StructureGetChain(pStruct, i);
    for (int ir = 0; ir < ChainGetResidueCount(pChainI); ir++)
    {
      Residue* pResIR = ChainGetResidue(pChainI, ir);
      if (pResIR->desType != Type_DesType_Fixed)
      {
        //the residue is a design site, set the coordinates of non-backbone atoms of designable residues to be FALSE
        for (int atom1 = 0; atom1 < ResidueGetAtomCount(pResIR); atom1++)
        {
          Atom* pAtom1 = ResidueGetAtom(pResIR, atom1);
          if (pAtom1->isBBAtom == FALSE)
          {
            pAtom1->isXyzValid = FALSE;
          }
        }
      }
    }
  }

  for (int i = 0; i < StructureGetChainCount(pStruct); i++)
  {
    Chain* pChainI = StructureGetChain(pStruct, i);
    for (int ir = 0; ir < ChainGetResidueCount(pChainI); ir++)
    {
      Residue* pResIR = ChainGetResidue(pChainI, ir);
      if (pResIR->desType == Type_DesType_Fixed && ChainGetType(pChainI) == Type_Chain_Protein)
      { // calculate protein residue internal energy only
        energyTerms[91] += pResIR->aapp;
        energyTerms[92] += pResIR->rama;
        energyTerms[93] += ResidueGetDunbrack(pResIR);
        AminoAcidReferenceEnergy(ResidueGetName(pResIR), energyTerms);
        EnergyIntraResidue(pResIR, energyTerms);
      }
      for (int is = ir + 1; is < ChainGetResidueCount(pChainI); is++)
      {
        Residue* pResIS = ChainGetResidue(pChainI, is);
        if (ResidueGetPosInChain(pResIR) + 1 == ResidueGetPosInChain(pResIS))
        {
          EnergyResidueAndNextResidue(pResIR, pResIS, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueSameChain(pResIR, pResIS, energyTerms);
        }
      }
      for (int k = i + 1; k < StructureGetChainCount(pStruct); k++)
      {
        Chain* pChainK = StructureGetChain(pStruct, k);
        for (int ks = 0; ks < ChainGetResidueCount(pChainK); ks++)
        {
          Residue* pResKS = ChainGetResidue(pChainK, ks);
          if (ChainGetType(pChainI) == Type_Chain_SmallMol && ChainGetType(pChainK) != Type_Chain_SmallMol)
          {
            EnergyResidueAndLigandResidue(pResKS, pResIR, energyTerms);
            if (FLAG_PROT_LIG == TRUE || FLAG_ENZYME == TRUE)
            {
              EnergyResidueAndLigandResidue(pResKS, pResIR, energyTermsBind);
            }
          }
          else if (ChainGetType(pChainK) == Type_Chain_SmallMol && ChainGetType(pChainI) != Type_Chain_SmallMol)
          {
            EnergyResidueAndLigandResidue(pResIR, pResKS, energyTerms);
            if (FLAG_PROT_LIG == TRUE || FLAG_ENZYME == TRUE)
            {
              EnergyResidueAndLigandResidue(pResIR, pResKS, energyTermsBind);
            }
          }
          //else if (ChainGetType(pChainK) != Type_Chain_DNA || ChainGetType(pChainI) != Type_Chain_DNA)
          else if (ChainGetType(pChainK) != Type_Chain_SmallMol && ChainGetType(pChainI) != Type_Chain_SmallMol)
          {
            EnergyResidueAndOtherResidueDiffChain(pResIR, pResKS, energyTerms);
            if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) == NULL)
              || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) != NULL))
            {
              if (FLAG_PPI == TRUE)
              {
                EnergyResidueAndOtherResidueDiffChain(pResIR, pResKS, energyTermsBind);
              }
            }
          }
          else
          {
            // do nothing if both chains are a small-molecule chain;
          }
        }
      }
    }
  }

  // restore the coordinates of non-backbone atoms of designable residues to be TRUE
  for (int i = 0; i < StructureGetChainCount(pStruct); i++)
  {
    Chain* pChainI = StructureGetChain(pStruct, i);
    for (int ir = 0; ir < ChainGetResidueCount(pChainI); ir++)
    {
      Residue* pResIR = ChainGetResidue(pChainI, ir);
      if (pResIR->desType != Type_DesType_Fixed)
      {
        for (int atom1 = 0; atom1 < ResidueGetAtomCount(pResIR); atom1++)
        {
          Atom* pAtom1 = ResidueGetAtom(pResIR, atom1);
          if (pAtom1->isBBAtom == FALSE)
          {
            pAtom1->isXyzValid = TRUE;
          }
        }
      }
    }
  }
  return Success;
}


int SequenceEnergy(Structure* pStruct, Sequence* pSeq)
{
  pSeq->etot = 0;
  pSeq->eevo = 0;
  pSeq->ephy = 0;
  pSeq->ebin = 0;
  if (FLAG_EVOLUTION == TRUE)
  {
    char seq[MAX_SEQ_LEN + 1] = "";
    StructureGetWholeSequence(pStruct, pSeq, seq);
    if (FLAG_EVOPHIPSI == TRUE)
    {
      pSeq->eevo += EvolutionScoreAllFromSeq(seq);
    }
    else
    {
      pSeq->eevo += EvolutionEnergyFromPSSMWithoutAlignment(seq);
    }
  }

  if (FLAG_PHYSICS == TRUE)
  {
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    double energyTermsBind[MAX_ENERGY_TERM] = { 0 };
    SequenceTemplateEnergy(pStruct, pSeq, energyTerms, energyTermsBind);
    for (int i = 0; i < StructureGetDesignSiteCount(pStruct); i++)
    {
      DesignSite* pSiteI = StructureGetDesignSite(pStruct, i);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotI = RotamerSetGet(pSetI, IntArrayGet(&pSeq->rotNdxs, i));
      Chain* pChainI = StructureGetChain(pStruct, pSiteI->chnNdx);
      RotamerRestore(pRotI, pSetI);
      for (int k = i; k < StructureGetDesignSiteCount(pStruct); k++)
      {
        if (k == i)
        {
          pSeq->ephy += pRotI->selfEnergy;
          pSeq->ebin += pRotI->selfEnergyBin;
        }
        else
        {
          DesignSite* pSiteK = StructureGetDesignSite(pStruct, k);
          RotamerSet* pSetK = DesignSiteGetRotamers(pSiteK);
          Rotamer* pRotK = RotamerSetGet(pSetK, IntArrayGet(&pSeq->rotNdxs, k));
          Chain* pChainK = StructureGetChain(pStruct, pSiteK->chnNdx);
          RotamerRestore(pRotK, pSetK);
          if (pSiteI->chnNdx == pSiteK->chnNdx)
          {// from the same chain
            if (ChainGetType(pChainI) == Type_Chain_Protein)
            {
              EnergyRotamerAndRotamerSameChain(pRotI, pRotK, energyTerms);
            }
          }
          else
          {// from different chains
            if (ChainGetType(pChainI) == Type_Chain_SmallMol || ChainGetType(pChainK) == Type_Chain_SmallMol)
            {
              double energyTermsPL[MAX_ENERGY_TERM] = { 0 };
              EnergyRotamerAndLigandRotamer(pRotI, pRotK, energyTermsPL);
              for (int index = 0; index < MAX_ENERGY_TERM; index++)
              {
                energyTerms[index] += energyTermsPL[index];
              }
              if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) == NULL)
                || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) != NULL))
              {
                if (FLAG_PROT_LIG == TRUE || FLAG_ENZYME == TRUE)
                {
                  for (int index = 0; index < MAX_ENERGY_TERM; index++)
                  {
                    energyTermsBind[index] += energyTermsPL[index];
                  }
                }
              }
            }
            //else if ( (ChainGetType(pChainI) == Type_Chain_Protein || ChainGetType(pChainI) == Type_Chain_DNA || ChainGetType(pChainI) == Type_Chain_RNA || ChainGetType(pChainI) == Type_Chain_Water)
            //  && (ChainGetType(pChainK) == Type_Chain_Protein || ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA || ChainGetType(pChainK) == Type_Chain_Water) )
            else
            {
              double energyTermsDC[MAX_ENERGY_TERM] = { 0 };
              EnergyRotamerAndRotamerDiffChain(pRotI, pRotK, energyTermsDC);
              for (int index = 0; index < MAX_ENERGY_TERM; index++)
              {
                energyTerms[index] += energyTermsDC[index];
              }
              if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) == NULL)
                || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) != NULL))
              {
                if (FLAG_PPI == TRUE)
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

    EnergyTermWeighting(energyTerms);
    pSeq->ephy += energyTerms[0];
    EnergyTermWeighting(energyTermsBind);
    pSeq->ebin += energyTermsBind[0];
  }

  pSeq->etot = WGT_PROFILE * pSeq->eevo + pSeq->ephy + (WGT_BIND - 1.0) * pSeq->ebin;

  return Success;
}


int EnergyDifferenceUponSingleMutation(Structure* pStruct, Sequence* pSeq, int mutSiteNdx, int mutRotNdx, double* dtot, double* dphy, double* dbin, double* devo)
{
  double phyBefore = 0;
  double phyAfter = 0;
  double binBefore = 0;
  double binAfter = 0;

  if (FLAG_PHYSICS == TRUE)
  {
    DesignSite* pCurSite = StructureGetDesignSite(pStruct, mutSiteNdx);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Chain* pCurChain = StructureGetChain(pStruct, pCurSite->chnNdx);
    Rotamer* pCurRot = RotamerSetGet(pCurSet, IntArrayGet(&pSeq->rotNdxs, mutSiteNdx));
    Rotamer* pNewRot = RotamerSetGet(pCurSet, mutRotNdx);
    RotamerRestore(pCurRot, pCurSet);
    RotamerRestore(pNewRot, pCurSet);
    phyBefore += pCurRot->selfEnergy;
    phyAfter += pNewRot->selfEnergy;
    if (FLAG_MONOMER == FALSE)
    {
      binBefore += pCurRot->selfEnergyBin;
      binAfter += pNewRot->selfEnergyBin;
    }
    double energyTerms1[MAX_ENERGY_TERM] = { 0 };
    double energyTerms2[MAX_ENERGY_TERM] = { 0 };
    double energyTermsBind1[MAX_ENERGY_TERM] = { 0 };
    double energyTermsBind2[MAX_ENERGY_TERM] = { 0 };
    for (int i = 0; i < StructureGetDesignSiteCount(pStruct); i++)
    {
      if (i == mutSiteNdx)
      {
        continue;
      }
      DesignSite* pSiteI = StructureGetDesignSite(pStruct, i);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotI = RotamerSetGet(pSetI, IntArrayGet(&pSeq->rotNdxs, i));
      Chain* pChainI = StructureGetChain(pStruct, pSiteI->chnNdx);
      RotamerRestore(pRotI, pSetI);
      if (pCurSite->chnNdx == pSiteI->chnNdx)
      {
        if (ChainGetType(pChainI) == Type_Chain_Protein)
        {
          EnergyRotamerAndRotamerSameChain(pRotI, pCurRot, energyTerms1);
          EnergyRotamerAndRotamerSameChain(pRotI, pNewRot, energyTerms2);
        }
      }
      else
      {
        if (ChainGetType(pChainI) == Type_Chain_SmallMol)
        {
          double energyTermsPL1[MAX_ENERGY_TERM] = { 0 };
          double energyTermsPL2[MAX_ENERGY_TERM] = { 0 };
          EnergyRotamerAndLigandRotamer(pCurRot, pRotI, energyTermsPL1);
          EnergyRotamerAndLigandRotamer(pNewRot, pRotI, energyTermsPL2);
          for (int index = 0; index < MAX_ENERGY_TERM; index++)
          {
            energyTerms1[index] += energyTermsPL1[index];
            energyTerms2[index] += energyTermsPL2[index];
          }
          if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) == NULL)
            || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) != NULL))
          {
            if (FLAG_PROT_LIG == TRUE || FLAG_ENZYME == TRUE)
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
          EnergyRotamerAndLigandRotamer(pRotI, pCurRot, energyTermsPL1);
          EnergyRotamerAndLigandRotamer(pRotI, pNewRot, energyTermsPL2);
          for (int index = 0; index < MAX_ENERGY_TERM; index++)
          {
            energyTerms1[index] += energyTermsPL1[index];
            energyTerms2[index] += energyTermsPL2[index];
          }
          if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) == NULL)
            || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) != NULL))
          {
            if (FLAG_PROT_LIG == TRUE || FLAG_ENZYME == TRUE)
            {
              for (int index = 0; index < MAX_ENERGY_TERM; index++)
              {
                energyTermsBind1[index] += energyTermsPL1[index];
                energyTermsBind2[index] += energyTermsPL2[index];
              }
            }
          }
        }
        //else if ( (ChainGetType(pChainI) == Type_Chain_Protein || ChainGetType(pChainI) == Type_Chain_DNA || ChainGetType(pChainI) == Type_Chain_RNA || ChainGetType(pChainI) == Type_Chain_Water)
        //  && (ChainGetType(pCurChain) == Type_Chain_Protein || ChainGetType(pCurChain) == Type_Chain_DNA || ChainGetType(pCurChain) == Type_Chain_RNA || ChainGetType(pCurChain) == Type_Chain_Water) )
        else
        {
          double energyTermsDC1[MAX_ENERGY_TERM] = { 0 };
          double energyTermsDC2[MAX_ENERGY_TERM] = { 0 };
          EnergyRotamerAndRotamerDiffChain(pRotI, pCurRot, energyTermsDC1);
          EnergyRotamerAndRotamerDiffChain(pRotI, pNewRot, energyTermsDC2);
          for (int index = 0; index < MAX_ENERGY_TERM; index++)
          {
            energyTerms1[index] += energyTermsDC1[index];
            energyTerms2[index] += energyTermsDC2[index];
          }
          if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) == NULL)
            || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) != NULL))
          {
            if (FLAG_PPI == TRUE)
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
    EnergyTermWeighting(energyTerms1);
    EnergyTermWeighting(energyTerms2);
    phyBefore += energyTerms1[0];
    phyAfter += energyTerms2[0];
    *dphy = phyAfter - phyBefore;

    if (FLAG_MONOMER == FALSE)
    {
      EnergyTermWeighting(energyTermsBind1);
      EnergyTermWeighting(energyTermsBind2);
      binBefore += energyTermsBind1[0];
      binAfter += energyTermsBind2[0];
      *dbin = binAfter - binBefore;
    }

    RotamerExtract(pCurRot);
    RotamerExtract(pNewRot);
  }

  if (FLAG_EVOLUTION == TRUE)
  {
    double evoBefore = 0;
    double evoAfter = 0;
    DesignSite* pCurSite = StructureGetDesignSite(pStruct, mutSiteNdx);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRot = RotamerSetGet(pCurSet, IntArrayGet(&pSeq->rotNdxs, mutSiteNdx));
    Rotamer* pNewRot = RotamerSetGet(pCurSet, mutRotNdx);
    if (RotamerAndRotamerInSameType(pCurRot, pNewRot) == FALSE)
    {
      char seq1[MAX_SEQ_LEN];
      StructureGetWholeSequence(pStruct, pSeq, seq1);
      char seq2[MAX_SEQ_LEN];
      Sequence newSeq;
      SequenceCreate(&newSeq);
      SequenceCopy(&newSeq, pSeq);
      IntArraySet(&newSeq.rotNdxs, mutSiteNdx, mutRotNdx);
      StructureGetWholeSequence(pStruct, &newSeq, seq2);
      if (FLAG_EVOPHIPSI)
      {
        evoBefore += EvolutionScoreAllFromSeq(seq1);
        evoAfter += EvolutionScoreAllFromSeq(seq2);
      }
      else
      {
        evoBefore += EvolutionEnergyFromPSSMWithoutAlignment(seq1);
        evoAfter += EvolutionEnergyFromPSSMWithoutAlignment(seq2);
      }
      SequenceDestroy(&newSeq);
    }
    *devo = evoAfter - evoBefore;
  }

  *dtot = WGT_PROFILE * (*devo) + *dphy + (WGT_BIND - 1.0) * (*dbin);

  return Success;
}


int Metropolis(Sequence* pSeq, Sequence* pBest, Structure* pStruct, RotamerList* pList, StringArray** ppRotType, IntArray** ppRotCount, int* seqNdx, double temp, int stepCount, FILE* pFileRot, FILE* pFileSeq)
{
  int nacc = 0;
  for (int i = 0; i < stepCount; i++)
  {
    int   mutSiteNdx;
    int   mutRotNdx;
    SequenceRandomSiteIndex(&mutSiteNdx, StructureGetDesignSiteCount(pStruct));
    SequenceRandRotamerIndex(pStruct, pList, ppRotType, ppRotCount, mutSiteNdx, &mutRotNdx);
    double dtot = 0, dphy = 0, devo = 0, dbin = 0;
    EnergyDifferenceUponSingleMutation(pStruct, pSeq, mutSiteNdx, mutRotNdx, &dtot, &dphy, &dbin, &devo);
    if (exp(-1.0 * dtot / temp) > (double)(rand() + 1.0) / (RAND_MAX + 1.0))
    {
      SequenceUpdateSingleSite(pSeq, mutSiteNdx, mutRotNdx);
      pSeq->etot += dtot;
      pSeq->ephy += dphy;
      pSeq->ebin += dbin;
      pSeq->eevo += devo;
      if (pSeq->etot < pBest->etot)
      {
        SequenceCopy(pBest, pSeq);
      }
      nacc++;
    }
  }
  *seqNdx += stepCount;
  float racc = 0;
  if (stepCount != 0)
  {
    racc = (float)nacc / stepCount;
  }
  else
  {
    racc = 0;
  }
  printf("acceptance ratio = %12.6f at temp = %12.6f\n", racc, temp);

  return Success;
}


int SequenceEnergyWithCons(Structure* pStruct, Sequence* pSeq, CataConsSitePairArray* pConsArray)
{
  pSeq->etot = 0;
  pSeq->eevo = 0;
  pSeq->ephy = 0;
  pSeq->ebin = 0;
  if (FLAG_EVOLUTION == TRUE)
  {
    char seq[MAX_SEQ_LEN + 1] = "";
    StructureGetWholeSequence(pStruct, pSeq, seq);
    if (FLAG_EVOPHIPSI == TRUE)
    {
      pSeq->eevo += EvolutionScoreAllFromSeq(seq);
    }
    else
    {
      pSeq->eevo += EvolutionEnergyFromPSSMWithoutAlignment(seq);
    }
  }

  if (FLAG_PHYSICS == TRUE)
  {
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    double energyTermsBind[MAX_ENERGY_TERM] = { 0 };
    SequenceTemplateEnergy(pStruct, pSeq, energyTerms, energyTermsBind);
    for (int i = 0; i < StructureGetDesignSiteCount(pStruct); i++)
    {
      DesignSite* pSiteI = StructureGetDesignSite(pStruct, i);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotI = RotamerSetGet(pSetI, IntArrayGet(&pSeq->rotNdxs, i));
      Chain* pChainI = StructureGetChain(pStruct, pSiteI->chnNdx);
      RotamerRestore(pRotI, pSetI);
      for (int k = i; k < StructureGetDesignSiteCount(pStruct); k++)
      {
        if (k == i)
        {
          pSeq->ephy += pRotI->selfEnergy;
          pSeq->ebin += pRotI->selfEnergyBin;
        }
        else
        {
          DesignSite* pSiteK = StructureGetDesignSite(pStruct, k);
          RotamerSet* pSetK = DesignSiteGetRotamers(pSiteK);
          Rotamer* pRotK = RotamerSetGet(pSetK, IntArrayGet(&pSeq->rotNdxs, k));
          Chain* pChainK = StructureGetChain(pStruct, pSiteK->chnNdx);
          RotamerRestore(pRotK, pSetK);
          if (pSiteI->chnNdx == pSiteK->chnNdx)
          { // from the same chain
            if (ChainGetType(pChainI) == Type_Chain_Protein)
            {
              EnergyRotamerAndRotamerSameChain(pRotI, pRotK, energyTerms);
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
              EnergyRotamerAndLigandRotamer(pRotI, pRotK, energyTermsPL);
              for (int index = 0; index < MAX_ENERGY_TERM; index++)
              {
                energyTerms[index] += energyTermsPL[index];
              }
              if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) == NULL)
                || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) != NULL))
              {
                if (FLAG_PROT_LIG == TRUE || FLAG_ENZYME == TRUE)
                {
                  for (int index = 0;index < MAX_ENERGY_TERM;index++)
                  {
                    energyTermsBind[index] += energyTermsPL[index];
                  }
                }
              }
            }
            //else if ((ChainGetType(pChainI) == Type_Chain_Protein || ChainGetType(pChainI) == Type_Chain_DNA || ChainGetType(pChainI) == Type_Chain_RNA || ChainGetType(pChainI) == Type_Chain_Water)
            //  && (ChainGetType(pChainK) == Type_Chain_Protein || ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA || ChainGetType(pChainK) == Type_Chain_Water))
            else
            {
              double energyTermsDC[MAX_ENERGY_TERM] = { 0 };
              EnergyRotamerAndRotamerDiffChain(pRotI, pRotK, energyTermsDC);
              for (int index = 0; index < MAX_ENERGY_TERM; index++)
              {
                energyTerms[index] += energyTermsDC[index];
              }
              if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) == NULL)
                || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pChainK)) != NULL))
              {
                if (FLAG_PPI == TRUE)
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

    EnergyTermWeighting(energyTerms);
    pSeq->ephy += energyTerms[0];
    EnergyTermWeighting(energyTermsBind);
    pSeq->ebin += energyTermsBind[0];
  }

  if (FLAG_ENZYME == TRUE)
  {
    for (int i = 0; i < CataConsSitePairArrayGetCount(pConsArray); i++)
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
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotI = RotamerSetGet(pSetI, IntArrayGet(&pSeq->rotNdxs, siteIndex1));
      RotamerRestore(pRotI, pSetI);
      if (siteIndex1 == siteIndex2)
      {
        if (!CataConsSitePairCheck(pCons, pRotI, pRotI))
        {
          pSeq->numOfUnsatisfiedCons++;
        }
      }
      else
      {
        DesignSite* pSiteK = StructureGetDesignSite(pStruct, siteIndex2);
        RotamerSet* pSetK = DesignSiteGetRotamers(pSiteK);
        Rotamer* pRotK = RotamerSetGet(pSetK, IntArrayGet(&pSeq->rotNdxs, siteIndex2));
        RotamerRestore(pRotK, pSetK);
        if (!CataConsSitePairCheck(pCons, pRotI, pRotK))
        {
          pSeq->numOfUnsatisfiedCons++;
        }
        RotamerExtract(pRotK);
      }
      RotamerExtract(pRotI);
    }
  }

  pSeq->etot = WGT_PROFILE * pSeq->eevo + pSeq->ephy + (WGT_BIND - 1.0) * pSeq->ebin + WGT_CATA_CONS * pSeq->numOfUnsatisfiedCons;
  return Success;
}


int EnergyChangeUponSingleMutationWithCataCons(Structure* pStruct, Sequence* pSeq, int mutSiteNdx, int mutRotNdx, double* dtot, double* dphy, double* dbin, double* devo, int* dcons, CataConsSitePairArray* pConsArray)
{
  double phyBefore = 0;
  double phyAfter = 0;
  double binBefore = 0;
  double binAfter = 0;
  double evoBefore = 0;
  double evoAfter = 0;
  int nConsBefore = 0;
  int nConsAfter = 0;

  if (FLAG_EVOLUTION == TRUE)
  {
    DesignSite* pCurSite = StructureGetDesignSite(pStruct, mutSiteNdx);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRot = RotamerSetGet(pCurSet, IntArrayGet(&pSeq->rotNdxs, mutSiteNdx));
    Rotamer* pNewRot = RotamerSetGet(pCurSet, mutRotNdx);
    if (RotamerAndRotamerInSameType(pCurRot, pNewRot) == FALSE)
    {
      char seq1[MAX_SEQ_LEN];
      StructureGetWholeSequence(pStruct, pSeq, seq1);
      char seq2[MAX_SEQ_LEN];
      Sequence newSeq;
      SequenceCreate(&newSeq);
      SequenceCopy(&newSeq, pSeq);
      IntArraySet(&newSeq.rotNdxs, mutSiteNdx, mutRotNdx);
      StructureGetWholeSequence(pStruct, &newSeq, seq2);
      if (FLAG_EVOPHIPSI == TRUE)
      {
        evoBefore += EvolutionScoreAllFromSeq(seq1);
        evoAfter += EvolutionScoreAllFromSeq(seq2);
      }
      else
      {
        evoBefore += EvolutionEnergyFromPSSMWithoutAlignment(seq1);
        evoAfter += EvolutionEnergyFromPSSMWithoutAlignment(seq2);
      }
      SequenceDestroy(&newSeq);
    }
    *devo = evoAfter - evoBefore;
  }

  if (FLAG_PHYSICS == TRUE)
  {
    DesignSite* pCurSite = StructureGetDesignSite(pStruct, mutSiteNdx);
    Chain* pCurChain = StructureGetChain(pStruct, pCurSite->chnNdx);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRot = RotamerSetGet(pCurSet, IntArrayGet(&pSeq->rotNdxs, mutSiteNdx));
    Rotamer* pNewRot = RotamerSetGet(pCurSet, mutRotNdx);
    RotamerRestore(pCurRot, pCurSet);
    RotamerRestore(pNewRot, pCurSet);

    phyBefore += pCurRot->selfEnergy;
    phyAfter += pNewRot->selfEnergy;
    binBefore += pCurRot->selfEnergyBin;
    binAfter += pNewRot->selfEnergyBin;

    double energyTerms1[MAX_ENERGY_TERM] = { 0 };
    double energyTerms2[MAX_ENERGY_TERM] = { 0 };
    double energyTermsBind1[MAX_ENERGY_TERM] = { 0 };
    double energyTermsBind2[MAX_ENERGY_TERM] = { 0 };

    for (int i = 0; i < StructureGetDesignSiteCount(pStruct); i++)
    {
      if (i == mutSiteNdx)
      {
        continue;
      }
      DesignSite* pSiteI = StructureGetDesignSite(pStruct, i);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotamerI = RotamerSetGet(pSetI, IntArrayGet(&pSeq->rotNdxs, i));
      Chain* pChainI = StructureGetChain(pStruct, pSiteI->chnNdx);
      RotamerRestore(pRotamerI, pSetI);
      if (pCurSite->chnNdx == pSiteI->chnNdx)
      {
        if (ChainGetType(pChainI) == Type_Chain_Protein)
        {
          EnergyRotamerAndRotamerSameChain(pRotamerI, pCurRot, energyTerms1);
          EnergyRotamerAndRotamerSameChain(pRotamerI, pNewRot, energyTerms2);
        }
      }
      else
      {
        if (ChainGetType(pChainI) == Type_Chain_SmallMol || ChainGetType(pCurChain) == Type_Chain_SmallMol)
        {
          double energyTermsPL1[MAX_ENERGY_TERM] = { 0 };
          double energyTermsPL2[MAX_ENERGY_TERM] = { 0 };
          EnergyRotamerAndLigandRotamer(pRotamerI, pCurRot, energyTermsPL1);
          EnergyRotamerAndLigandRotamer(pRotamerI, pNewRot, energyTermsPL2);
          for (int index = 0; index < MAX_ENERGY_TERM; index++)
          {
            energyTerms1[index] += energyTermsPL1[index];
            energyTerms2[index] += energyTermsPL2[index];
          }
          if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) == NULL)
            || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) != NULL))
          {
            if (FLAG_PROT_LIG == TRUE || FLAG_ENZYME == TRUE)
            {
              for (int index = 0; index < MAX_ENERGY_TERM; index++)
              {
                energyTermsBind1[index] += energyTermsPL1[index];
                energyTermsBind2[index] += energyTermsPL2[index];
              }
            }
          }
        }
        //else if ( (ChainGetType(pChainI) == Type_Chain_Protein || ChainGetType(pChainI) == Type_Chain_DNA || ChainGetType(pChainI) == Type_Chain_RNA || ChainGetType(pChainI) == Type_Chain_Water)
        //  && (ChainGetType(pCurChain) == Type_Chain_Protein || ChainGetType(pCurChain) == Type_Chain_DNA || ChainGetType(pCurChain) == Type_Chain_RNA || ChainGetType(pCurChain) == Type_Chain_Water) )
        else
        {
          double energyTermsDC1[MAX_ENERGY_TERM] = { 0 };
          double energyTermsDC2[MAX_ENERGY_TERM] = { 0 };
          EnergyRotamerAndRotamerDiffChain(pRotamerI, pCurRot, energyTermsDC1);
          EnergyRotamerAndRotamerDiffChain(pRotamerI, pNewRot, energyTermsDC2);
          for (int index = 0; index < MAX_ENERGY_TERM; index++)
          {
            energyTerms1[index] += energyTermsDC1[index];
            energyTerms2[index] += energyTermsDC2[index];
          }
          if ((strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) == NULL)
            || (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL && strstr(DES_CHAINS, ChainGetName(pCurChain)) != NULL))
          {
            if (FLAG_PPI == TRUE)
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
    EnergyTermWeighting(energyTerms1);
    EnergyTermWeighting(energyTerms2);
    EnergyTermWeighting(energyTermsBind1);
    EnergyTermWeighting(energyTermsBind2);
    phyBefore += energyTerms1[0];
    phyAfter += energyTerms2[0];
    *dphy = phyAfter - phyBefore;
    binBefore += energyTermsBind1[0];
    binAfter += energyTermsBind2[0];
    *dbin = binAfter - binBefore;
    RotamerExtract(pCurRot);
    RotamerExtract(pNewRot);
  }

  if (FLAG_ENZYME == TRUE)
  {
    DesignSite* pCurSite = StructureGetDesignSite(pStruct, mutSiteNdx);
    if (ResidueGetDesignType(pCurSite->pRes) == Type_DesType_Catalytic || ResidueGetDesignType(pCurSite->pRes) == Type_DesType_SmallMol)
    {
      RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
      Rotamer* pCurRot = RotamerSetGet(pCurSet, IntArrayGet(&pSeq->rotNdxs, mutSiteNdx));
      Rotamer* pNewRot = RotamerSetGet(pCurSet, mutRotNdx);
      if (pNewRot != pCurRot)
      {
        RotamerRestore(pCurRot, pCurSet);
        RotamerRestore(pNewRot, pCurSet);
        for (int i = 0; i < CataConsSitePairArrayGetCount(pConsArray); i++)
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
            RotamerSet* pSetK = DesignSiteGetRotamers(pSiteK);
            Rotamer* pRotK = RotamerSetGet(pSetK, IntArrayGet(&pSeq->rotNdxs, siteIndex2));
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
            RotamerSet* pSetK = DesignSiteGetRotamers(pSiteK);
            Rotamer* pRotK = RotamerSetGet(pSetK, IntArrayGet(&pSeq->rotNdxs, siteIndex1));
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

  *dtot = WGT_PROFILE * (*devo) + *dphy + (WGT_BIND - 1.0) * (*dbin) + WGT_CATA_CONS * (*dcons);
  return Success;
}



int MetropolisWithCataCons(Sequence* pOld, Sequence* pBest, Structure* pStructure, RotamerList* pList, StringArray** ppRotType, IntArray** ppRotCount, int* seqNdx, double temp, int stepCount, FILE* pFileRot, FILE* pFileSeq, CataConsSitePairArray* pConsArray)
{
  int nacc = 0;
  for (int i = 0; i < stepCount; i++)
  {
    int mutSiteNdx;
    int mutRotNdx;
    SequenceRandomSiteIndex(&mutSiteNdx, StructureGetDesignSiteCount(pStructure));
    SequenceRandRotamerIndex(pStructure, pList, ppRotType, ppRotCount, mutSiteNdx, &mutRotNdx);
    double dtot = 0, dphy = 0, devo = 0, dbin = 0;
    int dcons = 0;
    EnergyChangeUponSingleMutationWithCataCons(pStructure, pOld, mutSiteNdx, mutRotNdx, &dtot, &dphy, &dbin, &devo, &dcons, pConsArray);
    if (exp(-1.0 * dtot / temp) > (double)(rand() + 1.0) / (RAND_MAX + 1.0))
    {
      SequenceUpdateSingleSite(pOld, mutSiteNdx, mutRotNdx);
      pOld->etot += dtot;
      pOld->ephy += dphy;
      pOld->ebin += dbin;
      pOld->eevo += devo;
      pOld->numOfUnsatisfiedCons += dcons;
      if (pOld->etot < pBest->etot)
      {
        SequenceCopy(pBest, pOld);
      }
      nacc++;
    }
  }
  *seqNdx += stepCount;
  float racc = (float)nacc / stepCount;
  printf("acceptance ratio = %12.6f at temp = %12.6f\n", racc, temp);

  return Success;
}


int SimulatedAnnealing(Structure* pStruct, RotamerList* pList)
{
  int result = Success;
  char errMsg[MAX_LEN_ERR_MSG + 1];
  srand((unsigned int)time(NULL));

  StringArray* pRotTypes = (StringArray*)malloc(sizeof(StringArray) * pList->desSiteCount);
  IntArray* pRotCounts = (IntArray*)malloc(sizeof(IntArray) * pList->desSiteCount);
  DesignSiteShowRotamerTypeAndCount(pList, pStruct, &pRotTypes, &pRotCounts);
  int remainCount = 0;
  for (int i = 0; i < pList->desSiteCount; i++)
  {
    remainCount += pList->remainRotamerCount[i];
  }
  printf("total remained rotamer count: %d\n", remainCount);
  remainCount = remainCount < METROPOLIS_STEP ? remainCount : METROPOLIS_STEP;
  printf("set the number of monte carlo moves at each temperature to %d\n", remainCount);

  CataConsSitePairArray consArray;
  if (FLAG_ENZYME == TRUE)
  {
    result = CataConsSitePairArrayCreate(&consArray, FILE_CATACONS);
    if (FAILED(result))
    {
      sprintf(errMsg, "in file %s line %d, failed to create catalytic constraints from file %s",
        __FILE__, __LINE__, FILE_CATACONS);
      result = ValueError;
      return result;
    }
    for (int i = 0;i < CataConsSitePairArrayGetCount(&consArray);i++)
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

  printf("searching sequences using monte-carlo simulated annealing optimization\n");
  char BEST_SEQ_IDX[MAX_LEN_FILE_NAME + 1];
  sprintf(BEST_SEQ_IDX, "%s.txt", BEST_SEQ);
  FILE* pFileBestSeq = NULL;
  if (NTRAJ_START_NDX > 1)
  {
    pFileBestSeq = fopen(BEST_SEQ_IDX, "a");
  }
  else
  {
    pFileBestSeq = fopen(BEST_SEQ_IDX, "w");
    fprintf(pFileBestSeq, "# Comment lines start with the '#' symbol. Each designed sequence is written in one line with a few columns.\n"
      "# 1st column: the lowest-energy sequence\n"
      "# 2nd column: the index of the independent trajectory\n"
      "# 3rd column: the ratio of designable positions recapitulated to be the native amino acid\n"
      "# 4th column: the weighted total energy\n"
      "# 5th column: the unweighted evolutionary-profile energy calculated based on PSSM\n"
      "# 6th column: the unweighted physical energy\n"
      "# 7th column: the unweighted binding energy. This energy value is a part of the physical energy in the 5th column, and it is zero for non-interaction design)\n"
      "# 8th column: the number of unsatisfied catalytic constriants for enzyme design. For non-enzyme design, this value equals to zero\n"
    );
  }

  fflush(pFileBestSeq);
  for (int i = NTRAJ_START_NDX;i <= NTRAJ;i++)
  {
    printf("search for independent design trajectory #%d\n", i);
    Sequence oldSeq, bestSeq;
    SequenceCreate(&oldSeq);
    SequenceCreate(&bestSeq);
    if (FLAG_DESIGN_FROM_NATAA == TRUE)
    {
      SequenceGenNativeSeqSeed(pStruct, pList, &oldSeq);
    }
    else
    {
      SequenceGenRandomSeed(&oldSeq, pList);
    }
    if (FLAG_ENZYME == TRUE)
    {
      SequenceEnergyWithCons(pStruct, &oldSeq, &consArray);
    }
    else
    {
      SequenceEnergy(pStruct, &oldSeq);
    }
    SequenceCopy(&bestSeq, &oldSeq);

    int seqIndex = 0;
    // record decoy sequences and selected rots
    FILE* pFileRotDecoys = NULL;
    FILE* pFileSeqDecoys = NULL;
    //char ROT_FILE_ITER[MAX_LEN_FILE_NAME+1];
    //char SEQ_FILE_ITER[MAX_LEN_FILE_NAME+1];
    //sprintf(ROT_FILE_ITER,"%s%04d.txt",ROT_INDEX_FILE,i);
    //FILE* pFileRotDecoys=fopen(ROT_FILE_ITER,"w");
    //sprintf(SEQ_FILE_ITER,"%s%04d.txt",SEQ_FILE,i);
    //FILE* pFileSeqDecoys=fopen(SEQ_FILE_ITER,"w");
    //SequenceCopy(&bestSeq,&oldSeq);
    //SequenceWriteDesignRotamer(&bestSeq,pStruct,seqNdx,pFileRotDecoys);
    //SequenceWriteDesignFasta(&bestSeq,pStruct,seqNdx,pFileSeqDecoys);
    seqIndex++;
    //clock_t start = clock();
    for (int cycle = 1; cycle <= SA_CYCLE; cycle++)
    {
      printf("simulated annealing cycle %3d\n", cycle);
      //double t=SA_TMAX/pow(cycle,2.0);
      double t = SA_TMAX / cycle;
      while (t > SA_TMIN)
      {
        if (FLAG_ENZYME == TRUE)
        {
          MetropolisWithCataCons(&oldSeq, &bestSeq, pStruct, pList, &pRotTypes, &pRotCounts, &seqIndex, t, remainCount, pFileRotDecoys, pFileSeqDecoys, &consArray);
        }
        else
        {
          Metropolis(&oldSeq, &bestSeq, pStruct, pList, &pRotTypes, &pRotCounts, &seqIndex, t, remainCount, pFileRotDecoys, pFileSeqDecoys);
        }
        t *= SA_DECREASE_FAC;
      }
    }
    if (pFileRotDecoys != NULL)
    {
      fclose(pFileRotDecoys);
    }
    if (pFileSeqDecoys != NULL)
    {
      fclose(pFileSeqDecoys);
    }

    //clock_t end = clock();
    //printf("elapsed time for simulated annealing: %f sec\n", (double)(end - start) / CLOCKS_PER_SEC);

    // write the lowest-energy sequence
    char BEST_STRUCT_ITER[MAX_LEN_FILE_NAME + 1];
    sprintf(BEST_STRUCT_ITER, "%s%04d.pdb", BEST_STRUCT, i);
    strcpy(bestSeq.fileToSaveThisSeq, BEST_STRUCT_ITER);
    SequenceWriteDesignFasta(&bestSeq, pStruct, i, pFileBestSeq);
    fflush(pFileBestSeq);

    // write the lowest-energy protein structure, protein design sites, and corresponding ligand conformer
    DesignShowMinEnergyDesignStructure(pStruct, &bestSeq, BEST_STRUCT_ITER);
    char BEST_DESSITE_ITER[MAX_LEN_FILE_NAME + 1];
    sprintf(BEST_DESSITE_ITER, "%s%04d.pdb", BEST_DESSITE, i);
    DesignShowMinEnergyDesignSites(pStruct, &bestSeq, BEST_DESSITE_ITER);
    if (FLAG_MOL2 == TRUE)
    {
      char BEST_MOL2_ITER[MAX_LEN_FILE_NAME + 1];
      sprintf(BEST_MOL2_ITER, "%s%04d.mol2", BEST_MOL2, i);
      DesignShowMinEnergyDesignLigand(pStruct, &bestSeq, BEST_MOL2_ITER);
    }

    SequenceDestroy(&oldSeq);
    SequenceDestroy(&bestSeq);
  }
  fclose(pFileBestSeq);

  // release memory
  if (FLAG_ENZYME == TRUE)
  {
    CataConsSitePairArrayDestroy(&consArray);
  }
  for (int i = 0; i < pList->desSiteCount; i++)
  {
    StringArrayDestroy(&pRotTypes[i]);
    IntArrayDestroy(&pRotCounts[i]);
  }
  free(pRotTypes);
  free(pRotCounts);

  return Success;
}

