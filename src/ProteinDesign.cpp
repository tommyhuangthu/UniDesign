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

extern char BEST_SEQ[MAX_LENGTH_FILE_NAME+1];
extern char BEST_STRUCT[MAX_LENGTH_FILE_NAME+1];
extern char BEST_DESSITE[MAX_LENGTH_FILE_NAME+1];
extern char BEST_MOL2[MAX_LENGTH_FILE_NAME+1];
extern char ROT_INDEX_FILE[MAX_LENGTH_FILE_NAME+1];
extern char SEQ_FILE[MAX_LENGTH_FILE_NAME+1];
extern char FILE_CATALYTIC_CONSTRAINT[MAX_LENGTH_FILE_NAME+1];

extern int PROTEIN_LENGTH_NORMALIZATION;
extern int NTRAJ;

extern char PDBID[MAX_LENGTH_FILE_NAME+1];
extern char DES_CHAINS[MAX_LENGTH_ONE_LINE_IN_FILE+1];
extern double WGT_PROFILE;
extern float** PROT_PROFILE;


extern char PROGRAM_PATH[MAX_LENGTH_ONE_LINE_IN_FILE+1];
extern char TGT_PRF[MAX_LENGTH_FILE_NAME+1];
extern char TGT_SA[MAX_LENGTH_FILE_NAME+1];
extern char TGT_SS[MAX_LENGTH_FILE_NAME+1];
extern char TGT_SEQ[MAX_LENGTH_FILE_NAME+1];
extern char TGT_PHIPSI[MAX_LENGTH_FILE_NAME+1];

#define WGT_CATA_CONS  1000.0


int SitePairConsDeploy(CataConsSitePairArray *pSitePairArray,Structure *pStructure){
  for(int i=0; i<CataConsSitePairArrayGetCount(pSitePairArray); i++){
    CataConsSitePair *pSitePair   = CataConsSitePairArrayGet(pSitePairArray, i);
    DesignSite *pFirstDesignSite  = StructureFindDesignSiteByChainName(pStructure, pSitePair->chName1, pSitePair->pos1);
    DesignSite *pSecondDesignSite = StructureFindDesignSiteByChainName(pStructure, pSitePair->chName2, pSitePair->pos2);
    CataConsSitePairDeploy(pSitePair, DesignSiteGetRotamers(pFirstDesignSite), DesignSiteGetRotamers(pSecondDesignSite));
  }
  return Success;
}


int EnergyMatrixUpdateForCataCons(EnergyMatrix *pMatrix,RotamerList* pList,CataConsSitePair *pSitePair,Structure *pStructure){
  DesignSite *pFirstDesignSite  = StructureFindDesignSiteByChainName(pStructure, pSitePair->chName1, pSitePair->pos1);
  DesignSite *pSecondDesignSite = StructureFindDesignSiteByChainName(pStructure, pSitePair->chName2, pSitePair->pos2);
  int firstDesignSiteIndex      = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStructure, pSitePair->chName1, pSitePair->pos1);
  int secondDesignSiteIndex     = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStructure, pSitePair->chName2, pSitePair->pos2);
  EnergyMatrixBlock* pBlockII = EnergyMatrixGetBlock(pMatrix, firstDesignSiteIndex, firstDesignSiteIndex);
  EnergyMatrixBlock* pBlockKK = EnergyMatrixGetBlock(pMatrix, secondDesignSiteIndex, secondDesignSiteIndex);
  for(int j=0; j<pBlockII->RotamerCountSiteI; j++){
    int trueIndexIJ;
    RotamerOriginalIndexGet(pList, firstDesignSiteIndex, j, &trueIndexIJ);
    Rotamer* pRotamerJ = RotamerSetGet(DesignSiteGetRotamers(pFirstDesignSite), trueIndexIJ);
    for(int s=0; s<pBlockKK->RotamerCountSiteK; s++){
      int trueIndexKS;
      RotamerOriginalIndexGet(pList, secondDesignSiteIndex, s, &trueIndexKS);
      Rotamer* pRotamerS = RotamerSetGet(DesignSiteGetRotamers(pSecondDesignSite), trueIndexKS);
      if(CataConsSitePairCheck(pSitePair, pRotamerJ, pRotamerS) == FALSE){
        if(firstDesignSiteIndex < secondDesignSiteIndex) *EnergyMatrixGet(pMatrix, firstDesignSiteIndex, secondDesignSiteIndex, j, s) += 1e4;
        else *EnergyMatrixGet(pMatrix, secondDesignSiteIndex, firstDesignSiteIndex, s, j) += 1e4;
      }
    }
  }
  return Success;
}


int EnergyMatrixUpdateForCataConsArray(EnergyMatrix *pMatrix,RotamerList* pList,CataConsSitePairArray *pSitePairArray,Structure *pStructure){
  for(int i=0; i<CataConsSitePairArrayGetCount(pSitePairArray); i++){
    CataConsSitePair *pSitePair = CataConsSitePairArrayGet(pSitePairArray, i);
    EnergyMatrixUpdateForCataCons(pMatrix, pList, pSitePair, pStructure);
  }
  return Success;
}


int DesignSiteShowRotamerTypeAndCount(RotamerList* pList,Structure* pStructure,StringArray** ppRotamerType,IntArray** ppRotamerCount){
  StringArray* pStringArray = *ppRotamerType;
  IntArray*	 pIntArray	  = *ppRotamerCount;

  printf("show rotamer distribution at each design site\n");
  for(int i=0; i<StructureGetDesignSiteCount(pStructure); i++){
    int remainRotamerCount=0;
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure, i);
    RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    IntArrayCreate(&pIntArray[i], 0);
    StringArrayCreate(&pStringArray[i]);
    for(int j=0; j<pList->rotamerCount[i]; j++){
      if(pList->remainFlag[i][j]==FALSE) continue;
      BOOL     isRotamerTypeNew=TRUE;
      Rotamer* pRotamer=RotamerSetGet(pRotamerSet, j);
      for(int k=0; k<StringArrayGetCount(&pStringArray[i]); k++){
        if(strcmp(RotamerGetType(pRotamer),StringArrayGet(&pStringArray[i],k)) == 0){
          isRotamerTypeNew=FALSE;
          break;
        }
      }
      if(isRotamerTypeNew){
        int length=IntArrayGetLength(&pIntArray[i]);
        StringArrayAppend(&pStringArray[i],RotamerGetType(pRotamer));
        length++;
        IntArrayResize(&pIntArray[i],length);
        IntArraySet(&pIntArray[i],length-1,1);
      }else{
        int pos, rotamerCount;
        StringArrayFind(&pStringArray[i],RotamerGetType(pRotamer),&pos);
        rotamerCount=IntArrayGet(&pIntArray[i],pos);
        rotamerCount++;
        IntArraySet(&pIntArray[i],pos,rotamerCount);
      }
    }
    for(int j=0; j<IntArrayGetLength(&pIntArray[i]); j++){
      remainRotamerCount += IntArrayGet(&pIntArray[i],j);
    }
    printf("site %3d : %3s %s %4d, %6d rotamers:  ",i,
      ResidueGetName(pDesignSite->pResidue),
      ResidueGetChainName(pDesignSite->pResidue),
      ResidueGetPosInChain(pDesignSite->pResidue),
      remainRotamerCount);
    for(int j=0; j<StringArrayGetCount(&pStringArray[i]); j++){
      printf("%4d %s ",IntArrayGet(&pIntArray[i],j),StringArrayGet(&pStringArray[i],j));
    }
    printf("\n");
  }

  return Success;
}


int SequenceRandomSiteIndex(int *mutSiteIndex, int designSiteCount){
  *mutSiteIndex=rand()%designSiteCount;
  return Success;
}


int SequenceRandRotamerIndex(Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int siteIndex,int* rotIndex){
  StringArray* pStringArray = *ppRotamerType;
  IntArray*    pIntArray    = *ppRotamerCount;
  int typeIndex=rand()%IntArrayGetLength(&pIntArray[siteIndex]);
  int indexInType=rand()%IntArrayGet(&pIntArray[siteIndex],typeIndex);
  int typeflag = 0;
  for(int i = 0; i < pList->rotamerCount[siteIndex]; i++){
    Rotamer* pRotamer=RotamerSetGet(DesignSiteGetRotamers(StructureGetDesignSite(pStructure,siteIndex)),i);
    if(pList->remainFlag[siteIndex][i] == FALSE) continue;
    if(strcmp(StringArrayGet(&pStringArray[siteIndex],typeIndex),RotamerGetType(pRotamer)) != 0) continue;
    if(typeflag == indexInType){
      *rotIndex=i;
      break;
    }
    typeflag++;
  }
  return Success;
}


int SequenceChangeSingleSite(Sequence* pThis,int mutationSiteIndex,int mutationRotamerIndex){
  IntArraySet(&pThis->rotamerIndices,mutationSiteIndex,mutationRotamerIndex);
  return Success;
}


int SequenceGenerateRandomSeed(Sequence* pThis,RotamerList* pList){
  SequenceDestroy(pThis);
  SequenceCreate(pThis);
  pThis->designSiteCount=pList->designSiteCount;
  pThis->etot=0;
  pThis->ephy=0;
  pThis->eevo=0;
  IntArrayResize(&pThis->rotamerIndices, pThis->designSiteCount);
  for(int i=0; i<pThis->designSiteCount; i++){
    int j = rand()%pList->remainRotamerCount[i];
    int trueJ;
    RotamerOriginalIndexGet(pList, i, j, &trueJ);
    IntArraySet(&pThis->rotamerIndices, i, trueJ);
  }

  return Success;
}


int SequenceGenerateInitialSequenceSeed(Structure* pStructure,RotamerList* pList,Sequence* pSequence){
  SequenceDestroy(pSequence);
  SequenceCreate(pSequence);
  pSequence->designSiteCount=pList->designSiteCount;
  pSequence->etot=0.0;
  IntArrayResize(&pSequence->rotamerIndices, pSequence->designSiteCount);
  for(int i=0; i<pSequence->designSiteCount; i++){
    double minSelfEnergy = 1e8;
    int minSelfEnergyIndex = -1;
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure,i);
    RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    for(int j=0;j<pList->rotamerCount[i];j++){
      if(pList->remainFlag[i][j]==FALSE) continue;
      Rotamer* pRotamer=RotamerSetGet(pRotamerSet,j);
      if(pRotamer->selfenergy<minSelfEnergy){
        minSelfEnergy = pRotamer->selfenergy;
        minSelfEnergyIndex = j;
      }
    }
    IntArraySet(&pSequence->rotamerIndices, i, minSelfEnergyIndex);
  }

  return Success;
}


int SequenceGenerateNativeSequenceSeed(Structure* pStructure,RotamerList* pList,Sequence* pSequence){
  SequenceDestroy(pSequence);
  SequenceCreate(pSequence);
  pSequence->designSiteCount=pList->designSiteCount;
  pSequence->etot=0.0;
  IntArrayResize(&pSequence->rotamerIndices, pSequence->designSiteCount);
  for(int i=0; i<pSequence->designSiteCount; i++){
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure,i);
    RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    for(int j=0;j<RotamerSetGetCount(pRotamerSet);j++){
      Rotamer* pRotamer=RotamerSetGet(pRotamerSet,j);
      RotamerRestore(pRotamer,pRotamerSet);
      if(RotamerAndResidueSidechainRMSD(pRotamer,pDesignSite->pResidue)<0.05){
        IntArraySet(&pSequence->rotamerIndices, i, j);
        RotamerExtract(pRotamer);
        break;
      }
      RotamerExtract(pRotamer);
    }
  }

  return Success;
}


int StructureCalcSequenceTemplateEnergy(Structure* pStructure,Sequence* pSequence,double energyTerms[MAX_ENERGY_TERM],double energyTermsBind[MAX_ENERGY_TERM]){
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      if(pResIR->designType!=Type_ResidueDesignType_Fixed){
        //the residue is a design site, set the coordinates of non-backbone atoms to be FALSE
        for(int atom1=0; atom1<ResidueGetAtomCount(pResIR); atom1++){
          Atom* pAtom1 = ResidueGetAtom(pResIR, atom1);
          if(!pAtom1->isBBAtom) pAtom1->isXyzValid=FALSE;
        }
      }
    }
  }

  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain *pChainI=StructureGetChain(pStructure,i);
    for(int ir=0;ir<ChainGetResidueCount(pChainI);ir++){
      Residue *pResIR=ChainGetResidue(pChainI,ir);
      if(pResIR->designType==Type_ResidueDesignType_Fixed){
        energyTerms[91]+=pResIR->aapp;
        energyTerms[92]+=pResIR->rama;
        energyTerms[93]+=ResidueGetDunbrack(pResIR);
        AminoAcidReferenceEnergy(pResIR->name,energyTerms);
        EnergyIntraResidue(pResIR,energyTerms);
      }
      for(int is=ir+1;is<ChainGetResidueCount(pChainI);is++){
        Residue *pResIS=ChainGetResidue(pChainI,is);
        if(ResidueGetPosInChain(pResIR)+1==ResidueGetPosInChain(pResIS)) EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
        else EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
      }
      for(int k=i+1; k<StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        for(int ks=0; ks<ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          if(ChainGetType(pChainI)==Type_Chain_SmallMol){
            EnergyResidueAndLigResidue(pResKS,pResIR,energyTerms);
            if(FLAG_PROT_LIG || FLAG_ENZYME){
              EnergyResidueAndLigResidue(pResKS,pResIR,energyTermsBind);
            }
          }
          else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
            EnergyResidueAndLigResidue(pResIR,pResKS,energyTerms);
            if(FLAG_PROT_LIG || FLAG_ENZYME){
              EnergyResidueAndLigResidue(pResIR,pResKS,energyTermsBind);
            }
          }
          else{
            EnergyResidueAndOtherResidueDiffChain(pResIR,pResKS,energyTerms);
            //note that this code can easily cause a bug!!! (2/2/2020)
            if((strstr(DES_CHAINS,ChainGetName(pChainI))!=NULL && strstr(DES_CHAINS,ChainGetName(pChainK))==NULL) ||
              (strstr(DES_CHAINS,ChainGetName(pChainI))==NULL && strstr(DES_CHAINS,ChainGetName(pChainK))!=NULL)){
                if(FLAG_PPI){
                  EnergyResidueAndOtherResidueDiffChain(pResIR,pResKS,energyTermsBind);
                }
            }
          }
        }
      }
    }
  }

  //restore the coordinates
  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    Chain *pChainI=StructureGetChain(pStructure,i);
    for(int ir=0; ir<ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      if(pResIR->designType!=Type_ResidueDesignType_Fixed){
        for(int atom1=0; atom1<ResidueGetAtomCount(pResIR); atom1++){
          Atom* pAtom1=ResidueGetAtom(pResIR,atom1);
          if(pAtom1->isBBAtom==FALSE) pAtom1->isXyzValid=TRUE;
        }
      }
    }
  }
  return Success;
}


int StructureCalcSequenceEnergy(Structure* pStructure,Sequence* pSequence){
  pSequence->etot = 0;
  pSequence->eevo = 0;
  pSequence->ephy = 0;
  pSequence->ebin = 0;
  if(FLAG_EVOLUTION){
    char seq[MAX_SEQ_LEN+1]="";
    StructureGetWholeSequence(pStructure,pSequence,seq);
    if(FLAG_EVOPHIPSI) pSequence->eevo += EvolutionScoreAllFromSeq(seq);
    else pSequence->eevo += GetEvolutionScoreFromPrfWithoutAlignment(seq);
  }

  if(FLAG_PHYSICS){
    double energyTerms[MAX_ENERGY_TERM]={0};
    double energyTermsBind[MAX_ENERGY_TERM]={0};
    StructureCalcSequenceTemplateEnergy(pStructure,pSequence,energyTerms,energyTermsBind);
    for(int i=0;i<StructureGetDesignSiteCount(pStructure);i++){
      DesignSite* pSiteI = StructureGetDesignSite(pStructure,i);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotamerI = RotamerSetGet(pSetI,IntArrayGet(&pSequence->rotamerIndices,i));
      Chain* pChainI = StructureGetChain(pStructure,pSiteI->chainIndex);
      RotamerRestore(pRotamerI,pSetI);
      //RotamerShowAtomParameter(pRotamerI);
      //RotamerShowBondInformation(pRotamerI);
      for(int k=i;k<StructureGetDesignSiteCount(pStructure);k++){
        if(k==i){
          pSequence->ephy+=pRotamerI->selfenergy;
          pSequence->ebin+=pRotamerI->selfenergyBind;
        }
        else{
          DesignSite* pSiteK = StructureGetDesignSite(pStructure,k);
          RotamerSet* pSetK = DesignSiteGetRotamers(pSiteK);
          Rotamer* pRotamerK = RotamerSetGet(pSetK,IntArrayGet(&pSequence->rotamerIndices,k));
          Chain* pChainK = StructureGetChain(pStructure,pSiteK->chainIndex);
          RotamerRestore(pRotamerK,pSetK);
          if(pSiteI->chainIndex==pSiteK->chainIndex){
            if(ChainGetType(pChainI)==Type_Chain_Protein){
              EnergyRotamerAndRotamerSameChain(pRotamerI,pRotamerK,energyTerms);
            }
          }
          else{//different chains
            if(ChainGetType(pChainI)==Type_Chain_SmallMol || ChainGetType(pChainK)==Type_Chain_SmallMol){
                double energyTermsPL[MAX_ENERGY_TERM]={0};
                EnergyRotamerAndLigRotamer(pRotamerI,pRotamerK,energyTermsPL);
                for(int index=0;index<MAX_ENERGY_TERM;index++){
                  energyTerms[index]+=energyTermsPL[index];
                }
                if((strstr(DES_CHAINS,ChainGetName(pChainI))!=NULL && 
                  strstr(DES_CHAINS,ChainGetName(pChainK))==NULL) ||
                  (strstr(DES_CHAINS,ChainGetName(pChainI))==NULL && 
                  strstr(DES_CHAINS,ChainGetName(pChainK))!=NULL)){
                    if(FLAG_PROT_LIG || FLAG_ENZYME){
                      for(int index=0;index<MAX_ENERGY_TERM;index++){
                        energyTermsBind[index]+=energyTermsPL[index];
                      }
                    }
                }
            }
            else if((ChainGetType(pChainI)==Type_Chain_Protein || ChainGetType(pChainI)==Type_Chain_DNA || ChainGetType(pChainI)==Type_Chain_RNA) &&
              (ChainGetType(pChainK)==Type_Chain_Protein || ChainGetType(pChainK)==Type_Chain_DNA || ChainGetType(pChainK)==Type_Chain_RNA)){
                double energyTermsDC[MAX_ENERGY_TERM]={0};
                EnergyRotamerAndRotamerDiffChain(pRotamerI,pRotamerK,energyTermsDC);
                for(int index=0;index<MAX_ENERGY_TERM;index++){
                  energyTerms[index]+=energyTermsDC[index];
                  //energyTermsBind[index]+=energyTermsDC[index];
                }
                //note that this code can easily cause a bug!!! (2/2/2020)
                if((strstr(DES_CHAINS,ChainGetName(pChainI))!=NULL && 
                  strstr(DES_CHAINS,ChainGetName(pChainK))==NULL) ||
                  (strstr(DES_CHAINS,ChainGetName(pChainI))==NULL && 
                  strstr(DES_CHAINS,ChainGetName(pChainK))!=NULL)){
                    if(FLAG_PPI){
                      for(int index=0;index<MAX_ENERGY_TERM;index++){
                        energyTermsBind[index]+=energyTermsDC[index];
                      }
                    }
                }
            }
          }
          RotamerExtract(pRotamerK);
        }
      }
      RotamerExtract(pRotamerI);
    }

    EnergyTermWeighting(energyTerms);
    pSequence->ephy += energyTerms[0];
    EnergyTermWeighting(energyTermsBind);
    pSequence->ebin += energyTermsBind[0];
  }

  pSequence->etot = WGT_PROFILE*pSequence->eevo + pSequence->ephy;
  return Success;
}


int StructureCalcEnergyChangeUponSingleMutation(Structure* pStructure,Sequence* pSequence,int mutSiteIndex,int mutRotIndex,double* dtot,double* dphy,double* dbin,double *devo){
  double phyBefore = 0;
  double phyAfter = 0;
  double binBefore = 0;
  double binAfter = 0;
  double evoBefore=0;
  double evoAfter=0;

  if(FLAG_PHYSICS == TRUE){
    DesignSite* pCurSite=StructureGetDesignSite(pStructure,mutSiteIndex);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Chain* pCurChain = StructureGetChain(pStructure,pCurSite->chainIndex);
    Rotamer* pCurRotamer = RotamerSetGet(pCurSet,IntArrayGet(&pSequence->rotamerIndices,mutSiteIndex));
    Rotamer* pNewRotamer = RotamerSetGet(pCurSet,mutRotIndex);
    RotamerRestore(pCurRotamer,pCurSet);
    RotamerRestore(pNewRotamer,pCurSet);
    //get self energy (template)
    phyBefore += pCurRotamer->selfenergy;
    phyAfter  += pNewRotamer->selfenergy;
    binBefore += pCurRotamer->selfenergyBind;
    binAfter  += pNewRotamer->selfenergyBind;
    //get pairwise rotamer energy
    double energyTerms1[MAX_ENERGY_TERM]={0};
    double energyTerms2[MAX_ENERGY_TERM]={0};
    double energyTermsBind1[MAX_ENERGY_TERM]={0};
    double energyTermsBind2[MAX_ENERGY_TERM]={0};
    for(int i=0;i<pStructure->designSiteCount;i++){
      if(i==mutSiteIndex) continue;
      DesignSite* pSiteI = StructureGetDesignSite(pStructure,i);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotamerI = RotamerSetGet(pSetI,IntArrayGet(&pSequence->rotamerIndices,i));
      Chain* pChainI = StructureGetChain(pStructure,pSiteI->chainIndex);
      RotamerRestore(pRotamerI,pSetI);
      //get energy change
      if(pCurSite->chainIndex==pSiteI->chainIndex){
        if(ChainGetType(pChainI)==Type_Chain_Protein){
          EnergyRotamerAndRotamerSameChain(pRotamerI,pCurRotamer,energyTerms1);
          EnergyRotamerAndRotamerSameChain(pRotamerI,pNewRotamer,energyTerms2);
        }
      }
      else{
        if(ChainGetType(pChainI)==Type_Chain_SmallMol || ChainGetType(pCurChain)==Type_Chain_SmallMol){
            double energyTermsPL1[MAX_ENERGY_TERM]={0};
            double energyTermsPL2[MAX_ENERGY_TERM]={0};
            EnergyRotamerAndLigRotamer(pRotamerI,pCurRotamer,energyTermsPL1);
            EnergyRotamerAndLigRotamer(pRotamerI,pNewRotamer,energyTermsPL2);
            for(int index=0;index<MAX_ENERGY_TERM;index++){
              energyTerms1[index]+=energyTermsPL1[index];
              energyTerms2[index]+=energyTermsPL2[index];
            }
            if((strstr(DES_CHAINS,ChainGetName(pChainI))!=NULL && 
              strstr(DES_CHAINS,ChainGetName(pCurChain))==NULL) ||
              (strstr(DES_CHAINS,ChainGetName(pChainI))==NULL && 
              strstr(DES_CHAINS,ChainGetName(pCurChain))!=NULL)){
                if(FLAG_PROT_LIG || FLAG_ENZYME){
                  for(int index=0;index<MAX_ENERGY_TERM;index++){
                    energyTermsBind1[index]+=energyTermsPL1[index];
                    energyTermsBind2[index]+=energyTermsPL2[index];
                  }
                }
            }
        }
        else if((ChainGetType(pChainI)==Type_Chain_Protein || ChainGetType(pChainI)==Type_Chain_DNA || ChainGetType(pChainI)==Type_Chain_RNA) &&
          (ChainGetType(pCurChain)==Type_Chain_Protein || ChainGetType(pCurChain)==Type_Chain_DNA || ChainGetType(pCurChain)==Type_Chain_RNA)){
            double energyTermsDC1[MAX_ENERGY_TERM]={0};
            double energyTermsDC2[MAX_ENERGY_TERM]={0};
            EnergyRotamerAndRotamerDiffChain(pRotamerI,pCurRotamer,energyTermsDC1);
            EnergyRotamerAndRotamerDiffChain(pRotamerI,pNewRotamer,energyTermsDC2);
            for(int index=0;index<MAX_ENERGY_TERM;index++){
              energyTerms1[index]+=energyTermsDC1[index];
              energyTerms2[index]+=energyTermsDC2[index];
            }
            if((strstr(DES_CHAINS,ChainGetName(pChainI))!=NULL && 
              strstr(DES_CHAINS,ChainGetName(pCurChain))==NULL) ||
              (strstr(DES_CHAINS,ChainGetName(pChainI))==NULL && 
              strstr(DES_CHAINS,ChainGetName(pCurChain))!=NULL)){
                if(FLAG_PPI){
                  for(int index=0;index<MAX_ENERGY_TERM;index++){
                    energyTermsBind1[index]+=energyTermsDC1[index];
                    energyTermsBind2[index]+=energyTermsDC2[index];
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
    phyAfter  += energyTerms2[0];
    *dphy=phyAfter-phyBefore;
    binBefore += energyTermsBind1[0];
    binAfter  += energyTermsBind2[0];
    *dbin=binAfter-binBefore;
    RotamerExtract(pCurRotamer);
    RotamerExtract(pNewRotamer);
  }

  if(FLAG_EVOLUTION){
    DesignSite* pCurSite=StructureGetDesignSite(pStructure,mutSiteIndex);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRotamer = RotamerSetGet(pCurSet,IntArrayGet(&pSequence->rotamerIndices,mutSiteIndex));
    Rotamer* pNewRotamer = RotamerSetGet(pCurSet,mutRotIndex);
    if(RotamerAndRotamerInSameType(pCurRotamer,pNewRotamer)==FALSE){
      char seq1[MAX_SEQ_LEN];
      StructureGetWholeSequence(pStructure,pSequence,seq1);
      char seq2[MAX_SEQ_LEN];
      Sequence newSeq;
      SequenceCreate(&newSeq);
      SequenceCopy(&newSeq,pSequence);
      IntArraySet(&newSeq.rotamerIndices,mutSiteIndex,mutRotIndex);
      StructureGetWholeSequence(pStructure,&newSeq,seq2);
      if(FLAG_EVOPHIPSI){
        evoBefore += EvolutionScoreAllFromSeq(seq1);
        evoAfter += EvolutionScoreAllFromSeq(seq2);
      }
      else{
        evoBefore += GetEvolutionScoreFromPrfWithoutAlignment(seq1);
        evoAfter += GetEvolutionScoreFromPrfWithoutAlignment(seq2);
      }
      SequenceDestroy(&newSeq);
    }
    *devo=evoAfter-evoBefore;
  }

  *dtot = WGT_PROFILE*(*devo) + *dphy;
  return Success;
}


int MetropolisCriterion(Sequence* pOld,Sequence* pBest,Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int *seqIndex,double temp,int stepCount, FILE* fp,FILE *fp2){
  for(int i=0; i<stepCount; i++){
    int   mutationSiteIndex;
    int   mutationRotamerIndex;
    SequenceRandomSiteIndex(&mutationSiteIndex,StructureGetDesignSiteCount(pStructure));
    SequenceRandRotamerIndex(pStructure,pList,ppRotamerType,ppRotamerCount,mutationSiteIndex,&mutationRotamerIndex);
    double dtot=0,dphy=0,devo=0,dbin=0;
    StructureCalcEnergyChangeUponSingleMutation(pStructure,pOld,mutationSiteIndex,mutationRotamerIndex,&dtot,&dphy,&dbin,&devo);
    if(exp(-1.0*dtot/temp)>(double)(rand()+1.0)/(RAND_MAX+1.0)){
      SequenceChangeSingleSite(pOld,mutationSiteIndex,mutationRotamerIndex);
      pOld->etot += dtot;
      pOld->ephy += dphy;
      pOld->ebin += dbin;
      pOld->eevo += devo;
      //printf("step=%8d,temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => accepted\n",*seqIndex+i,temp,dtot,dphy,devo);
      //printf("temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => accepted\n",temp,pOld->tote,pOld->phye,pOld->evoe);
      //ouput all accepted sequences
      SequenceWriteDesignRotamer(pOld,pStructure,*seqIndex+i,fp);
      SequenceWriteDesignFasta(pOld,pStructure,*seqIndex+i,fp2);
      if(pOld->etot<pBest->etot){
        SequenceCopy(pBest,pOld);
      }
    }
    /*else{
      printf("step=%8d,temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => rejected\n",*seqIndex+i,temp,dtot,dphy,devo);
    }*/
  }
  *seqIndex+=stepCount;
  return Success;
}


int StructureCalcSequenceEnergyWithCataCons(Structure* pStructure,Sequence* pSequence,CataConsSitePairArray* pConsArray){
  pSequence->etot = 0;
  pSequence->eevo = 0;
  pSequence->ephy = 0;
  pSequence->ebin = 0;
  if(FLAG_EVOLUTION){
    char seq[MAX_SEQ_LEN+1]="";
    StructureGetWholeSequence(pStructure,pSequence,seq);
    if(FLAG_EVOPHIPSI) pSequence->eevo += EvolutionScoreAllFromSeq(seq);
    else pSequence->eevo += GetEvolutionScoreFromPrfWithoutAlignment(seq);
  }

  if(FLAG_PHYSICS){
    double energyTerms[MAX_ENERGY_TERM]={0};
    double energyTermsBind[MAX_ENERGY_TERM]={0};
    //calculate the template energy for the conformation-fixed residues
    //get template energy
    StructureCalcSequenceTemplateEnergy(pStructure,pSequence,energyTerms,energyTermsBind);
    //get rotamer pairwise energy
    for(int i=0;i<StructureGetDesignSiteCount(pStructure);i++){
      DesignSite* pSiteI = StructureGetDesignSite(pStructure,i);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotamerI = RotamerSetGet(pSetI,IntArrayGet(&pSequence->rotamerIndices,i));
      RotamerRestore(pRotamerI,pSetI);
      //RotamerShowAtomParameter(pRotamerI);
      //RotamerShowBondInformation(pRotamerI);
      for(int k=i;k<StructureGetDesignSiteCount(pStructure);k++){
        if(k==i){
          pSequence->ephy+=pRotamerI->selfenergy;
          pSequence->ebin+=pRotamerI->selfenergyBind;
        }
        else{
          DesignSite* pSiteK = StructureGetDesignSite(pStructure,k);
          RotamerSet* pSetK = DesignSiteGetRotamers(pSiteK);
          Rotamer* pRotamerK = RotamerSetGet(pSetK,IntArrayGet(&pSequence->rotamerIndices,k));
          RotamerRestore(pRotamerK,pSetK);
          if(pSiteI->chainIndex==pSiteK->chainIndex){
            if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_Protein){
              EnergyRotamerAndRotamerSameChain(pRotamerI,pRotamerK,energyTerms);
            }
          }
          else{//different chains
            if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_SmallMol ||
              ChainGetType(StructureGetChain(pStructure,pSiteK->chainIndex))==Type_Chain_SmallMol){
                if(pSiteI->pResidue->designType==Type_ResidueDesignType_Catalytic || pSiteK->pResidue->designType==Type_ResidueDesignType_Catalytic){
                  CataConsSitePair* pCons1 = CataConsSitePairArrayFind(pConsArray,ResidueGetChainName(pSiteI->pResidue),ResidueGetPosInChain(pSiteI->pResidue),ResidueGetName(pSiteI->pResidue),
                    ResidueGetChainName(pSiteK->pResidue),ResidueGetPosInChain(pSiteK->pResidue),ResidueGetName(pSiteK->pResidue));
                  CataConsSitePair* pCons2 = CataConsSitePairArrayFind(pConsArray,ResidueGetChainName(pSiteK->pResidue),ResidueGetPosInChain(pSiteK->pResidue),ResidueGetName(pSiteK->pResidue),
                    ResidueGetChainName(pSiteI->pResidue),ResidueGetPosInChain(pSiteI->pResidue),ResidueGetName(pSiteI->pResidue));
                  // do not calculate energy between sites forming a covalent bond
                  if((pCons1!=NULL && pCons1->pairConsType==Type_SitePairCons_Covalent) || (pCons2!=NULL && pCons2->pairConsType==Type_SitePairCons_Covalent)) continue;
                }
                double energyTermsPL[MAX_ENERGY_TERM]={0};
                EnergyRotamerAndLigRotamer(pRotamerI,pRotamerK,energyTermsPL);
                for(int index=0;index<MAX_ENERGY_TERM;index++){
                  energyTerms[index]+=energyTermsPL[index];
                }
                if((strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))!=NULL && 
                  strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteK->chainIndex)))==NULL) ||
                  (strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))==NULL && 
                  strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteK->chainIndex)))!=NULL)){
                    if(FLAG_PROT_LIG || FLAG_ENZYME){
                      for(int index=0;index<MAX_ENERGY_TERM;index++){
                        energyTermsBind[index]+=energyTermsPL[index];
                      }
                    }
                }
            }
            else if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_Protein &&
              ChainGetType(StructureGetChain(pStructure,pSiteK->chainIndex))==Type_Chain_Protein){
                double energyTermsDC[MAX_ENERGY_TERM]={0};
                EnergyRotamerAndRotamerDiffChain(pRotamerI,pRotamerK,energyTermsDC);
                for(int index=0;index<MAX_ENERGY_TERM;index++){
                  energyTerms[index]+=energyTermsDC[index];
                  //energyTermsBind[index]+=energyTermsDC[index];
                }
                //note that this code can easily cause a bug!!! (2/2/2020)
                if((strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))!=NULL && 
                  strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteK->chainIndex)))==NULL) ||
                  (strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))==NULL && 
                  strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteK->chainIndex)))!=NULL)){
                    if(FLAG_PPI){
                      for(int index=0;index<MAX_ENERGY_TERM;index++){
                        energyTermsBind[index]+=energyTermsDC[index];
                      }
                    }
                }
            }
          }
          RotamerExtract(pRotamerK);
        }
      }
      RotamerExtract(pRotamerI);
    }

    EnergyTermWeighting(energyTerms);
    pSequence->ephy += energyTerms[0];
    EnergyTermWeighting(energyTermsBind);
    pSequence->ebin += energyTermsBind[0];
  }

  if(FLAG_ENZYME){
    for(int i=0; i<CataConsSitePairArrayGetCount(pConsArray); i++){
      CataConsSitePair* pCons = CataConsSitePairArrayGet(pConsArray, i);
      int siteIndex1 = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStructure, pCons->chName1, pCons->pos1);
      int siteIndex2 = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStructure, pCons->chName2, pCons->pos2);
      if(siteIndex1 == -1 || siteIndex2 == -1){
        char errMsg[MAX_LENGTH_ERR_MSG+1];
        sprintf(errMsg, "in file %s line %d, %s%d or %s%d has not been added as a catalytic design site",__FILE__,__LINE__,
          pCons->chName1, pCons->pos1,
          pCons->chName2, pCons->pos2);
        TraceError(errMsg, ValueError);
        continue;
      }
      DesignSite* pSiteI = StructureGetDesignSite(pStructure,siteIndex1);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotamerI = RotamerSetGet(pSetI,IntArrayGet(&pSequence->rotamerIndices,siteIndex1));
      RotamerRestore(pRotamerI,pSetI);
      if(siteIndex1 == siteIndex2){
        if(!CataConsSitePairCheck(pCons, pRotamerI, pRotamerI)) pSequence->numOfUnsatisfiedCons++;
      }
      else{
        DesignSite* pSiteK = StructureGetDesignSite(pStructure,siteIndex2);
        RotamerSet* pSetK = DesignSiteGetRotamers(pSiteK);
        Rotamer* pRotamerK = RotamerSetGet(pSetK,IntArrayGet(&pSequence->rotamerIndices,siteIndex2));
        RotamerRestore(pRotamerK,pSetK);
        if(!CataConsSitePairCheck(pCons, pRotamerI, pRotamerK)) pSequence->numOfUnsatisfiedCons++;
        RotamerExtract(pRotamerK);
      }
      RotamerExtract(pRotamerI);
    }
  }

  //normalize by length
  //pSequence->eevo /= sqrt((double)TOT_SEQ_LEN);
  //pSequence->ephy /= sqrt((double)TOT_SEQ_LEN);
  //do not normalize binding energy
  //pSequence->ebin /= sqrt((double)TOT_SEQ_LEN);
  pSequence->etot = WGT_PROFILE*pSequence->eevo + 1.0*pSequence->ephy + WGT_CATA_CONS*pSequence->numOfUnsatisfiedCons;
  return Success;
}


int StructureCalcEnergyChangeUponSingleMutationWithCataCons(Structure* pStructure,Sequence* pSequence,int mutSiteIndex,int mutRotIndex,double* dtot,double* dphy,double* dbin,double *devo,int *dcons,CataConsSitePairArray* pConsArray){
  double phyBefore = 0;
  double phyAfter = 0;
  double binBefore = 0;
  double binAfter = 0;
  double evoBefore=0;
  double evoAfter=0;
  int nConsBefore=0;
  int nConsAfter=0;

  if(FLAG_EVOLUTION){
    DesignSite* pCurSite=StructureGetDesignSite(pStructure,mutSiteIndex);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRotamer = RotamerSetGet(pCurSet,IntArrayGet(&pSequence->rotamerIndices,mutSiteIndex));
    Rotamer* pNewRotamer = RotamerSetGet(pCurSet,mutRotIndex);
    if(RotamerAndRotamerInSameType(pCurRotamer,pNewRotamer)==FALSE){
      char seq1[MAX_SEQ_LEN];
      StructureGetWholeSequence(pStructure,pSequence,seq1);
      char seq2[MAX_SEQ_LEN];
      Sequence newSeq;
      SequenceCreate(&newSeq);
      SequenceCopy(&newSeq,pSequence);
      IntArraySet(&newSeq.rotamerIndices,mutSiteIndex,mutRotIndex);
      StructureGetWholeSequence(pStructure,&newSeq,seq2);
      if(FLAG_EVOPHIPSI){
        evoBefore += EvolutionScoreAllFromSeq(seq1);
        evoAfter += EvolutionScoreAllFromSeq(seq2);
      }
      else{
        evoBefore += GetEvolutionScoreFromPrfWithoutAlignment(seq1);
        evoAfter += GetEvolutionScoreFromPrfWithoutAlignment(seq2);
      }
      SequenceDestroy(&newSeq);
    }
    *devo=evoAfter-evoBefore;
  }

  if(FLAG_PHYSICS){
    DesignSite* pCurSite=StructureGetDesignSite(pStructure,mutSiteIndex);
    RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
    Rotamer* pCurRotamer = RotamerSetGet(pCurSet,IntArrayGet(&pSequence->rotamerIndices,mutSiteIndex));
    Rotamer* pNewRotamer = RotamerSetGet(pCurSet,mutRotIndex);
    RotamerRestore(pCurRotamer,pCurSet);
    RotamerRestore(pNewRotamer,pCurSet);
    //get self energy (template)
    phyBefore += pCurRotamer->selfenergy;
    phyAfter  += pNewRotamer->selfenergy;
    binBefore += pCurRotamer->selfenergyBind;
    binAfter  += pNewRotamer->selfenergyBind;
    //get pairwise rotamer energy
    double energyTerms1[MAX_ENERGY_TERM]={0};
    double energyTerms2[MAX_ENERGY_TERM]={0};
    double energyTermsBind1[MAX_ENERGY_TERM]={0};
    double energyTermsBind2[MAX_ENERGY_TERM]={0};
    for(int i=0;i<pStructure->designSiteCount;i++){
      if(i==mutSiteIndex) continue;
      DesignSite* pSiteI = StructureGetDesignSite(pStructure,i);
      RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
      Rotamer* pRotamerI = RotamerSetGet(pSetI,IntArrayGet(&pSequence->rotamerIndices,i));
      RotamerRestore(pRotamerI,pSetI);
      //get energy change
      if(pCurSite->chainIndex==pSiteI->chainIndex){
        if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_Protein){
          EnergyRotamerAndRotamerSameChain(pRotamerI,pCurRotamer,energyTerms1);
          EnergyRotamerAndRotamerSameChain(pRotamerI,pNewRotamer,energyTerms2);
        }
      }
      else{
        if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_SmallMol ||
          ChainGetType(StructureGetChain(pStructure,pCurSite->chainIndex))==Type_Chain_SmallMol){
            double energyTermsPL1[MAX_ENERGY_TERM]={0};
            double energyTermsPL2[MAX_ENERGY_TERM]={0};
            EnergyRotamerAndLigRotamer(pRotamerI,pCurRotamer,energyTermsPL1);
            EnergyRotamerAndLigRotamer(pRotamerI,pNewRotamer,energyTermsPL2);
            for(int index=0;index<MAX_ENERGY_TERM;index++){
              energyTerms1[index]+=energyTermsPL1[index];
              energyTerms2[index]+=energyTermsPL2[index];
            }
            if((strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))!=NULL && 
              strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pCurSite->chainIndex)))==NULL) ||
              (strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))==NULL && 
              strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pCurSite->chainIndex)))!=NULL)){
                if(FLAG_PROT_LIG || FLAG_ENZYME){
                  for(int index=0;index<MAX_ENERGY_TERM;index++){
                    energyTermsBind1[index]+=energyTermsPL1[index];
                    energyTermsBind2[index]+=energyTermsPL2[index];
                  }
                }
            }
        }
        else if(ChainGetType(StructureGetChain(pStructure,pSiteI->chainIndex))==Type_Chain_Protein &&
          ChainGetType(StructureGetChain(pStructure,pCurSite->chainIndex))==Type_Chain_Protein){
            double energyTermsDC1[MAX_ENERGY_TERM]={0};
            double energyTermsDC2[MAX_ENERGY_TERM]={0};
            EnergyRotamerAndRotamerDiffChain(pRotamerI,pCurRotamer,energyTermsDC1);
            EnergyRotamerAndRotamerDiffChain(pRotamerI,pNewRotamer,energyTermsDC2);
            for(int index=0;index<MAX_ENERGY_TERM;index++){
              energyTerms1[index]+=energyTermsDC1[index];
              energyTerms2[index]+=energyTermsDC2[index];
            }
            if((strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))!=NULL && 
              strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pCurSite->chainIndex)))==NULL) ||
              (strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pSiteI->chainIndex)))==NULL && 
              strstr(DES_CHAINS,ChainGetName(StructureGetChain(pStructure,pCurSite->chainIndex)))!=NULL)){
                if(FLAG_PPI){
                  for(int index=0;index<MAX_ENERGY_TERM;index++){
                    energyTermsBind1[index]+=energyTermsDC1[index];
                    energyTermsBind2[index]+=energyTermsDC2[index];
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
    phyAfter  += energyTerms2[0];
    *dphy=phyAfter-phyBefore;
    binBefore += energyTermsBind1[0];
    binAfter  += energyTermsBind2[0];
    *dbin=binAfter-binBefore;
    RotamerExtract(pCurRotamer);
    RotamerExtract(pNewRotamer);
  }

  if(FLAG_ENZYME){
    DesignSite* pCurSite=StructureGetDesignSite(pStructure,mutSiteIndex);
    if(ResidueGetDesignType(pCurSite->pResidue) == Type_ResidueDesignType_Catalytic || ResidueGetDesignType(pCurSite->pResidue) == Type_ResidueDesignType_SmallMol){
      RotamerSet* pCurSet = DesignSiteGetRotamers(pCurSite);
      Rotamer* pCurRotamer = RotamerSetGet(pCurSet,IntArrayGet(&pSequence->rotamerIndices,mutSiteIndex));
      Rotamer* pNewRotamer = RotamerSetGet(pCurSet,mutRotIndex);
      if(pNewRotamer!=pCurRotamer){
        RotamerRestore(pCurRotamer,pCurSet);
        RotamerRestore(pNewRotamer,pCurSet);
        for(int i=0; i<CataConsSitePairArrayGetCount(pConsArray); i++){
          CataConsSitePair* pCons = CataConsSitePairArrayGet(pConsArray, i);
          int siteIndex1 = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStructure, pCons->chName1, pCons->pos1);
          int siteIndex2 = StructureFindDesignSiteIndexByChainNameAndPosInChain(pStructure, pCons->chName2, pCons->pos2);
          if(siteIndex1 == -1 || siteIndex2 == -1) continue;
          if(siteIndex1 == mutSiteIndex && siteIndex2 != mutSiteIndex){
            DesignSite* pSiteK = StructureGetDesignSite(pStructure,siteIndex2);
            RotamerSet* pSetK = DesignSiteGetRotamers(pSiteK);
            Rotamer* pRotamerK = RotamerSetGet(pSetK,IntArrayGet(&pSequence->rotamerIndices,siteIndex2));
            RotamerRestore(pRotamerK,pSetK);
            if(!CataConsSitePairCheck(pCons, pCurRotamer, pRotamerK)) nConsBefore++;
            if(!CataConsSitePairCheck(pCons, pNewRotamer, pRotamerK)) nConsAfter++;
            RotamerExtract(pRotamerK);
          }
          else if(siteIndex1 != mutSiteIndex && siteIndex2 == mutSiteIndex){
            DesignSite* pSiteK = StructureGetDesignSite(pStructure,siteIndex1);
            RotamerSet* pSetK = DesignSiteGetRotamers(pSiteK);
            Rotamer* pRotamerK = RotamerSetGet(pSetK,IntArrayGet(&pSequence->rotamerIndices,siteIndex1));
            RotamerRestore(pRotamerK,pSetK);
            if(!CataConsSitePairCheck(pCons, pRotamerK, pCurRotamer)) nConsBefore++;
            if(!CataConsSitePairCheck(pCons, pRotamerK, pNewRotamer)) nConsAfter++;
            RotamerExtract(pRotamerK);
          }
          else if(siteIndex1 == mutSiteIndex && siteIndex2 == mutSiteIndex){
            if(CataConsSitePairCheck(pCons, pCurRotamer, pCurRotamer)==FALSE) nConsBefore++;
            if(CataConsSitePairCheck(pCons, pNewRotamer, pNewRotamer)==FALSE) nConsAfter++;
          }
        }
        RotamerExtract(pNewRotamer);
        RotamerExtract(pCurRotamer);
      }
    }
    *dcons=nConsAfter-nConsBefore;
  }

  *dtot = WGT_PROFILE*(*devo) + 1.0*(*dphy) + WGT_CATA_CONS*(*dcons);
  return Success;
}



int MetropolisCriterionWithCataCons(Sequence* pOld,Sequence* pBest,Structure* pStructure,RotamerList* pList,StringArray** ppRotamerType,IntArray** ppRotamerCount,int *seqIndex,double temp,int stepCount, FILE* pFileRot,FILE *pFileSeq,CataConsSitePairArray* pConsArray){
  for(int i=0; i<stepCount; i++){
    int   mutationSiteIndex;
    int   mutationRotamerIndex;
    SequenceRandomSiteIndex(&mutationSiteIndex,StructureGetDesignSiteCount(pStructure));
    SequenceRandRotamerIndex(pStructure,pList,ppRotamerType,ppRotamerCount,mutationSiteIndex,&mutationRotamerIndex);
    double dtot=0,dphy=0,devo=0,dbin=0;
    int dcons=0;
    StructureCalcEnergyChangeUponSingleMutationWithCataCons(pStructure,pOld,mutationSiteIndex,mutationRotamerIndex,&dtot,&dphy,&dbin,&devo,&dcons,pConsArray);
    if(exp(-1.0*dtot/temp)>(double)(rand()+1.0)/(RAND_MAX+1.0)){
      SequenceChangeSingleSite(pOld,mutationSiteIndex,mutationRotamerIndex);
      pOld->etot += dtot;
      pOld->ephy += dphy;
      pOld->ebin += dbin;
      pOld->eevo += devo;
      pOld->numOfUnsatisfiedCons += dcons;
      //printf("step=%8d,temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f, dcons=%d => accepted\n",*seqIndex+i,temp,dtot,dphy,devo,dcons);
      //printf("temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => accepted\n",temp,pOld->tote,pOld->phye,pOld->evoe);
      SequenceWriteDesignRotamer(pOld,pStructure,*seqIndex+i,pFileRot);
      SequenceWriteDesignFasta(pOld,pStructure,*seqIndex+i,pFileSeq);
      if(pOld->etot<pBest->etot){
        SequenceCopy(pBest,pOld);
      }
    }
    /*else{
      printf("step=%8d,temp=%12.6f,dtot=%12.6f, dphy=%12.6f, devo=%12.6f => rejected\n",*seqIndex+i,temp,dtot,dphy,devo);
    }*/
  }
  *seqIndex+=stepCount;
  return Success;
}


int SimulatedAnnealingOptimization(Structure* pStructure,RotamerList* pList){
  int result = Success;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  srand((unsigned int)time(NULL));

  // show rotamer types and counts at each design site
  StringArray *pRotTypes=(StringArray*)malloc(sizeof(StringArray)*pList->designSiteCount);
  IntArray *pRotCounts=(IntArray*)malloc(sizeof(IntArray)*pList->designSiteCount);
  DesignSiteShowRotamerTypeAndCount(pList,pStructure,&pRotTypes,&pRotCounts);
  int remainCount=0;
  for(int i=0; i<pList->designSiteCount; i++) remainCount += pList->remainRotamerCount[i];
  printf("total remained rotamer count: %d\n", remainCount);
  remainCount = remainCount<METROPOLIS_STEP ? remainCount:METROPOLIS_STEP;
  printf("set the number of monte carlo moves at each temperature to %d\n", remainCount);

  CataConsSitePairArray consArray;
  if(FLAG_ENZYME){
    result = CataConsSitePairArrayCreate(&consArray,FILE_CATALYTIC_CONSTRAINT);
    if(FAILED(result)){
      sprintf(errMsg,"in file %s line %d, failed to create catalytic constraints from file %s",
        __FILE__,__LINE__,FILE_CATALYTIC_CONSTRAINT);
      result = ValueError;
      return result;
    }
    for(int i=0;i<CataConsSitePairArrayGetCount(&consArray);i++){
      result = StructureDeployCataConsSitePair(pStructure,CataConsSitePairArrayGet(&consArray,i));
      if(FAILED(result)){
        sprintf(errMsg,"in file %s line %d, failed to deploy catalytic constraints",
          __FILE__,__LINE__);
        result = ValueError;
        return result;
      }
    }
  }

  printf("searching sequences using monte-carlo simulated annealing optimization\n");
  char BEST_SEQ_IDX[MAX_LENGTH_FILE_NAME+1];
  sprintf(BEST_SEQ_IDX,"%s.txt",BEST_SEQ);
  FILE* pFileBestSeq=fopen(BEST_SEQ_IDX,"w");
  fprintf(pFileBestSeq,"# Comment lines start with the '#' symbol. Each designed sequence is written in one line with a few columns.\n"
    "# 1st column: the lowest-energy sequence\n"
    "# 2nd column: the index of the independent trajectory\n"
    "# 3rd column: the ratio of designable positions recapitulated to be the native amino acid\n"
    "# 4th column: the weighted total energy\n"
    "# 5th column: the unweighted evolutionary-profile energy calculated based on PSSM\n"
    "# 6th column: the unweighted physical energy\n"
    "# 7th column: the unweighted binding energy. This energy value is a part of the physical energy in the 5th column, and it is zero for non-interaction design)\n"
    "# 8th column: the number of unsatisfied catalytic constriants for enzyme design. This valuezero for non-enzyme design\n"
    );
  fflush(pFileBestSeq);
  for(int i=1;i<=NTRAJ;i++){
    printf("search for independent design trajectory #%d\n",i);
    Sequence oldSequence,bestSequence;
    SequenceCreate(&oldSequence);
    SequenceCreate(&bestSequence);
    if(FLAG_DESIGN_FROM_NATAA) SequenceGenerateNativeSequenceSeed(pStructure,pList,&oldSequence);
    else SequenceGenerateRandomSeed(&oldSequence,pList);
    if(FLAG_ENZYME) StructureCalcSequenceEnergyWithCataCons(pStructure,&oldSequence,&consArray);
    else StructureCalcSequenceEnergy(pStructure,&oldSequence);
    int seqIndex=0;

    // record decoy sequences and selected rots
    char ROT_FILE_ITER[MAX_LENGTH_FILE_NAME+1];
    char SEQ_FILE_ITER[MAX_LENGTH_FILE_NAME+1];
    sprintf(ROT_FILE_ITER,"%s%04d.txt",ROT_INDEX_FILE,i);
    FILE* pFileRotDecoys=fopen(ROT_FILE_ITER,"w");
    sprintf(SEQ_FILE_ITER,"%s%04d.txt",SEQ_FILE,i);
    FILE* pFileSeqDecoys=fopen(SEQ_FILE_ITER,"w");
    SequenceCopy(&bestSequence,&oldSequence);
    SequenceWriteDesignRotamer(&bestSequence,pStructure,seqIndex,pFileRotDecoys);
    SequenceWriteDesignFasta(&bestSequence,pStructure,seqIndex,pFileSeqDecoys);
    seqIndex++;
    for(int cycle = 1; cycle <= SA_CYCLE; cycle++){
      printf("simulated annealing cycle %3d\n", cycle);
      //double t=SA_TMAX/pow(cycle,2.0);
      double t=SA_TMAX;
      while(t>SA_TMIN){
        if(FLAG_ENZYME) MetropolisCriterionWithCataCons(&oldSequence,&bestSequence,pStructure,pList,&pRotTypes,&pRotCounts,&seqIndex,t,remainCount,pFileRotDecoys,pFileSeqDecoys,&consArray);
        else MetropolisCriterion(&oldSequence,&bestSequence,pStructure,pList,&pRotTypes,&pRotCounts,&seqIndex,t,remainCount,pFileRotDecoys,pFileSeqDecoys);
        t *= SA_DECREASE_FAC;
      }
    }
    fclose(pFileRotDecoys);
    fclose(pFileSeqDecoys);

    // record the best design, writing the sequence, protein, and ligand structure if applicable
    // record the lowest-energy sequence
    SequenceWriteDesignFasta(&bestSequence,pStructure,i,pFileBestSeq);
    fflush(pFileBestSeq);
    // show lowest-energy protein and ligand structures
    char BEST_STRUCT_ITER[MAX_LENGTH_FILE_NAME+1];
    sprintf(BEST_STRUCT_ITER,"%s%04d.pdb",BEST_STRUCT,i);
    DesignShowMinEnergyDesignStructure(pStructure,&bestSequence,BEST_STRUCT_ITER);
    char BEST_DESSITE_ITER[MAX_LENGTH_FILE_NAME+1];
    sprintf(BEST_DESSITE_ITER,"%s%04d.pdb",BEST_DESSITE,i);
    DesignShowMinEnergyDesignSites(pStructure,&bestSequence,BEST_DESSITE_ITER);
    if(FLAG_MOL2){
      char BEST_MOL2_ITER[MAX_LENGTH_FILE_NAME+1];
      sprintf(BEST_MOL2_ITER,"%s%04d.mol2",BEST_MOL2,i);
      DesignShowMinEnergyDesignLigand(pStructure,&bestSequence,BEST_MOL2_ITER);
    }


    SequenceDestroy(&oldSequence);
    SequenceDestroy(&bestSequence);
  }
  fclose(pFileBestSeq);

  // release memory
  if(FLAG_ENZYME) CataConsSitePairArrayDestroy(&consArray);
  for(int i=0; i<pList->designSiteCount; i++){
    StringArrayDestroy(&pRotTypes[i]);
    IntArrayDestroy(&pRotCounts[i]);
  }
  free(pRotTypes);
  free(pRotCounts);

  return Success;
}

