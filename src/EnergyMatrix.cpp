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

#pragma warning(disable:4090)
#pragma warning(disable:4101)
#include "EnergyMatrix.h"
#include <string.h>
#include <time.h>

extern char DES_CHAINS[MAX_LENGTH_ONE_LINE_IN_FILE+1];
extern BOOL FLAG_MONOMER;
extern BOOL FLAG_PPI;
extern BOOL FLAG_PROT_LIG;
extern BOOL FLAG_ENZYME;

int EnergyMatrixBlockCreate(EnergyMatrixBlock* pThis){
  pThis->DesignSiteI = 0;
  pThis->DesignSiteK = 0;
  pThis->RotamerCountSiteI = 0;
  pThis->RotamerCountSiteK = 0;
  pThis->energyIK = NULL;
  return Success;
}


void EnergyMatrixBlockDestroy(EnergyMatrixBlock* pThis){
  if(pThis->energyIK!=NULL){
    free(pThis->energyIK);
    pThis->energyIK = NULL;
  }
}


int EnergyMatrixBlockCopy(EnergyMatrixBlock* pThis,EnergyMatrixBlock* pOther){
  EnergyMatrixBlockDestroy(pThis);
  EnergyMatrixBlockCreate(pThis);
  pThis->DesignSiteI = pOther->DesignSiteI;
  pThis->DesignSiteK = pOther->DesignSiteK;
  pThis->RotamerCountSiteI = pOther->RotamerCountSiteI;
  pThis->RotamerCountSiteK = pOther->RotamerCountSiteK;
  pThis->energyIK = (double*)malloc(sizeof(double)*pThis->RotamerCountSiteI*pThis->RotamerCountSiteK);
  memcpy(pThis->energyIK, pOther->energyIK, sizeof(double)*pThis->RotamerCountSiteI*pThis->RotamerCountSiteK);
  return Success;
}


double* EnergyMatrixBlockGet(EnergyMatrixBlock* pThis,int J,int S){
  return &pThis->energyIK[J*pThis->RotamerCountSiteK + S];
}


int EnergyMatrixBlockGenerate(EnergyMatrixBlock* pThis,Structure *pStructure,int designSiteI,int designSiteK,FILE* outputFile){
  DesignSite *pDesignSiteI = StructureGetDesignSite(pStructure, designSiteI);
  DesignSite *pDesignSiteK = StructureGetDesignSite(pStructure, designSiteK);
  if(designSiteI>designSiteK || !pDesignSiteI || !pDesignSiteK) return ValueError;

  pThis->DesignSiteI = designSiteI;
  pThis->DesignSiteK = designSiteK;
  pThis->RotamerCountSiteI = RotamerSetGetCount(DesignSiteGetRotamers(pDesignSiteI));
  pThis->RotamerCountSiteK = RotamerSetGetCount(DesignSiteGetRotamers(pDesignSiteK));
  pThis->energyIK = (double*)malloc(sizeof(double)*pThis->RotamerCountSiteI*pThis->RotamerCountSiteK);

  if(designSiteI==designSiteK){
    for(int j=0;j<pThis->RotamerCountSiteI;j++){
      Rotamer *pRotamerIJ  = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), j);
      RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));
      double* pEnergy = EnergyMatrixBlockGet(pThis,j,j);
      *pEnergy = 0.0;
      double energyTerms[MAX_ENERGY_TERM]={0};
      Chain *pChainI = StructureGetChain(pStructure, pDesignSiteI->chainIndex);
      if(ChainGetType(pChainI) == Type_Chain_Protein){
        AminoAcidReferenceEnergy(RotamerGetType(pRotamerIJ), energyTerms);
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(b==pDesignSiteI->resiIndex){
                EnergyIntraRotamer(pRotamerIJ, energyTerms);
              }
              else{
                if(ChainGetType(pChainA) != Type_Chain_SmallMol){
                  if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                    EnergyRotamerAndFixedResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                  else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                    pResidueAB->designType == Type_ResidueDesignType_Designable){
                      EnergyRotamerAndDesignResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                }
              }
            }
          }
          else{//different chains
            if(ChainGetType(pChainA) != Type_Chain_SmallMol){
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                  EnergyRotamerAndFixedResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                  pResidueAB->designType == Type_ResidueDesignType_Designable){
                    EnergyRotamerAndDesignResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
          }
        }
      }
      else if(ChainGetType(pChainI) == Type_Chain_SmallMol){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){
            //*pEnergy += pRotamerIJ->vdwInternal;
          }
          else{
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                EnergyRotamerAndFixedResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
              }
              else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                pResidueAB->designType == Type_ResidueDesignType_Designable){
                  EnergyRotamerAndDesignResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
              }
            }
          }
        }
      }

      EnergyTermWeighting(energyTerms);
      *pEnergy += energyTerms[0];
      fprintf(outputFile,"%f %d %d %d %d\n",*pEnergy,designSiteI,j,designSiteK,j);
      RotamerExtract(pRotamerIJ);
    }
    return Success;
  }

  for(int j=0; j<pThis->RotamerCountSiteI; j++){
    Rotamer *pRotamerIJ = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), j);
    RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));
    for(int s=0; s<pThis->RotamerCountSiteK; s++){
      Rotamer *pRotamerKS = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteK), s);
      RotamerRestore(pRotamerKS,DesignSiteGetRotamers(pDesignSiteK));
      double* pEnergy = EnergyMatrixBlockGet(pThis,j,s);
      *pEnergy = 0.0;

      double energyTerms[MAX_ENERGY_TERM]={0};
      if(pDesignSiteI->chainIndex == pDesignSiteK->chainIndex){
        EnergyRotamerAndRotamerSameChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      else{
        EnergyRotamerAndRotamerDiffChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      EnergyTermWeighting(energyTerms);
      for(int a=1; a<MAX_ENERGY_TERM;a++){
        *pEnergy += energyTerms[a];
      }
      fprintf(outputFile,"%f %d %d %d %d\n",*pEnergy,designSiteI,j,designSiteK,s);
      RotamerExtract(pRotamerKS);
    }
    RotamerExtract(pRotamerIJ);
  }
  return Success;
}


int EnergyMatrixCreate(EnergyMatrix* pThis){
  pThis->designSiteCount = 0;
  pThis->blocks = NULL;
  return Success;
}


void EnergyMatrixDestroy(EnergyMatrix* pThis){
  int i,j;
  for(i=0;i<pThis->designSiteCount;i++){
    for(j=i;j<pThis->designSiteCount;j++){
      EnergyMatrixBlockDestroy(&pThis->blocks[ i*pThis->designSiteCount + j ]);
    }
  }
  free(pThis->blocks);
  pThis->blocks = NULL;
}


int EnergyMatrixCopy(EnergyMatrix* pThis,EnergyMatrix* pOther){
  EnergyMatrixDestroy(pThis);
  EnergyMatrixCreate(pThis);
  pThis->designSiteCount = pOther->designSiteCount;
  pThis->blocks = (EnergyMatrixBlock*)malloc(sizeof(EnergyMatrixBlock)*pThis->designSiteCount*pThis->designSiteCount);
  for(int i=0;i<pThis->designSiteCount;i++){
    for(int j=i;j<pThis->designSiteCount;j++){
      EnergyMatrixBlockCreate(EnergyMatrixGetBlock(pThis,i,j));
      EnergyMatrixBlockCopy(EnergyMatrixGetBlock(pThis,i,j),EnergyMatrixGetBlock(pOther,i,j));
    }
  }
  return Success;
}


int EnergyMatrixRead(EnergyMatrix* pThis,char* filepath){
  int result = Success;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  FILE* pFile = fopen(filepath,"r");
  if(!pFile){
    result = IOError;
    sprintf(errMsg,"in file %s line %d, cannot read file %s",__FILE__,__LINE__,filepath);
    TraceError(errMsg,result);
    return result;
  }

  fseek(pFile,0,SEEK_END);
  int fileLength = ftell(pFile);
  fseek(pFile,0,SEEK_SET);

  printf("reading energy matrix from file %s\n",filepath);
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while( fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    double energy;
    int designSiteI,designSiteK,rotamerIJ,rotamerKS;
    EnergyMatrixBlock* energyIK;

    sscanf(line,"%lf %d %d %d %d",&energy,&designSiteI,&rotamerIJ,&designSiteK,&rotamerKS);
    if(designSiteI<0 || designSiteK<0 || rotamerIJ<0 || rotamerKS<0 || designSiteI>designSiteK){
      result = FormatError;
      sprintf(errMsg,"in file %s line %d, bad format at this line of file %s\n%s",__FILE__,__LINE__,filepath,line);
      TraceError(errMsg,result);
      return result;
    }

    if(designSiteI>=pThis->designSiteCount || designSiteK>=pThis->designSiteCount){
      int newDesignSiteCount = designSiteI>designSiteK? designSiteI+1 : designSiteK+1;
      EnergyMatrixBlock* newBlocks = (EnergyMatrixBlock*)malloc(
        sizeof(EnergyMatrixBlock) * newDesignSiteCount * newDesignSiteCount);
      for(int i=0;i<newDesignSiteCount;i++){
        for(int j=0;j<newDesignSiteCount;j++){
          EnergyMatrixBlockCreate(&newBlocks[i*newDesignSiteCount+j]);
        }
      }
      for(int i=0;i<pThis->designSiteCount;i++){
        for(int j=i;j<pThis->designSiteCount;j++){
          newBlocks[i*newDesignSiteCount+j] = *EnergyMatrixGetBlock(pThis,i,j);
        }
      }

      pThis->designSiteCount = newDesignSiteCount;
      free(pThis->blocks);
      pThis->blocks = newBlocks;
    }

    energyIK = EnergyMatrixGetBlock(pThis,designSiteI,designSiteK);
    if(energyIK->energyIK != NULL){
      result = FormatError;
      sprintf(errMsg,"in file %s line %d, energy block between site %d and %d has already been read, when reading file %s at this line:\n%s",__FILE__,__LINE__,designSiteI,designSiteK,filepath,line);
      TraceError(errMsg,result);
      return result;
    }
    if(rotamerIJ >= energyIK->RotamerCountSiteI){
      energyIK->RotamerCountSiteI = rotamerIJ+1;
    }
    if(rotamerKS >= energyIK->RotamerCountSiteK){
      energyIK->RotamerCountSiteK = rotamerKS+1;
    }

  }

  for(int i=0;i<pThis->designSiteCount;i++){
    for(int j=i;j<pThis->designSiteCount;j++){
      int k;
      EnergyMatrixBlock* pEnergyIK = EnergyMatrixGetBlock(pThis,i,j);
      if(pEnergyIK->RotamerCountSiteI && pEnergyIK->RotamerCountSiteK && !pEnergyIK->energyIK){
        pEnergyIK->energyIK = (double*)malloc(
          sizeof(double)*pEnergyIK->RotamerCountSiteI*pEnergyIK->RotamerCountSiteK);
      }
      for(k=0;k<pEnergyIK->RotamerCountSiteI*pEnergyIK->RotamerCountSiteK;k++){
        pEnergyIK->energyIK[k] = 0.0;
      }
      pEnergyIK->DesignSiteI = i;
      pEnergyIK->DesignSiteK = j;
    }
  }

  //read for a second time
  fseek(pFile,0,SEEK_SET);
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    double energy;
    int designSiteI,designSiteK,rotamerIJ,rotamerKS;
    sscanf(line,"%lf %d %d %d %d",&energy,&designSiteI,&rotamerIJ,&designSiteK,&rotamerKS);
    *EnergyMatrixGet(pThis,designSiteI,designSiteK,rotamerIJ,rotamerKS) = energy;
  }

  fclose(pFile);

  return Success;
}


int EnergyMatrixCheck(char* matrixFile){
  int MAX_DESIGN_COUNT = 1000;
  int siteCount;
  int rotamerCount[1000];
  int i,j,k,s;
  double value;
  int lineCounter;
  int fileLength;
  int result = Success;
  char errMsg[MAX_LENGTH_ERR_MSG+1];

  siteCount = 0;
  for(i=0; i<1000; i++) rotamerCount[i] = 0;

  FILE* pFile = fopen(matrixFile,"r");
  if(!pFile){
    sprintf(errMsg,"in file %s line %d, cannot read file %s",__FILE__,__LINE__,matrixFile);
    result = IOError;
    TraceError(errMsg,result);
    return result;
  }

  fseek(pFile,0,SEEK_END);
  fileLength = ftell(pFile);
  fseek(pFile,0,SEEK_SET);

  lineCounter = 0;
  printf("Counting design sites and rotamers...\n");
  while( fscanf(pFile,"%lf %d %d %d %d",&value,&i,&j,&k,&s)!=EOF ){
    lineCounter++;
    if(lineCounter%100000 == 0){
      printf("\r");ShowProgress(40,100.0*ftell(pFile)/fileLength);
    }
    if(i<0 || i>=MAX_DESIGN_COUNT || k<0 || k>=MAX_DESIGN_COUNT){
      printf("error in line %d, i and k must between 0 and %d :\n%f %d %d %d %d\n",lineCounter,MAX_DESIGN_COUNT,value,i,j,k,s);
      return FormatError;
    }
    if(i>=siteCount){
      siteCount = i+1;
    }
    if(k>=siteCount){
      siteCount = k+1;
    }
    if(j>=rotamerCount[i]){
      rotamerCount[i] = j+1;
    }
    if(s>=rotamerCount[k]){
      rotamerCount[k] = s+1;
    }
  }
  printf("\n totally %d sites, %d lines in file\n",siteCount,lineCounter);;
  for(i=0;i<siteCount;i++){
    printf("Design site #%d, %d rotamers\n",i,rotamerCount[i]);
  }

  fseek(pFile,0,SEEK_SET);
  lineCounter = 0;
  for(i=0;i<siteCount;i++){
    for(k=i;k<siteCount;k++){
      printf("Checking energy block i = %d, k = %d, size= %d*%d...",
        i,k,rotamerCount[i],rotamerCount[k]);
      for(j=0;j<rotamerCount[i];j++){
        for(s=0;s<rotamerCount[k];s++){
          int iread,jread,kread,sread,result;
          if(i==k && j!=s){
            continue;
          }

          lineCounter++;
          result = fscanf(pFile,"%lf %d %d %d %d",&value,&iread,&jread,&kread,&sread);
          if(result == EOF){
            printf("\nAt the end of file, energy item i=%d, j=%d, k=%d, s=%d is missing\n",
              i,j,k,s);
            return FormatError;
          }
          if(iread!=i || jread!=j || kread!=k || sread!=s){
            printf("\nAt line #%d, expect i=%d, j=%d, k=%d,s=%d, but this line is read:\n"
              "%f %d %d %d %d",lineCounter,i,j,k,s,value,iread,jread,kread,sread);
            return FormatError;
          }
          if(fabs(value)>1e8){
            printf("\nAt line #%d, energy value is too large:\n"
              "%f %d %d %d %d",lineCounter,value,i,j,k,s);
            return FormatError;
          }
        }
      }
      printf("OK\n");
    }
  }
  printf("all blocks of energy matrix file %s are okay\n",matrixFile);

  fclose(pFile);
  return Success;

}


int EnergyMatrixWrite(EnergyMatrix* pThis,char* filepath){
  int result = Success;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  FILE* pFile = fopen(filepath,"w");
  if(!pFile){
    result = IOError;
    sprintf(errMsg,"in file %s line %d, cannot write to file %s",__FILE__,__LINE__,filepath);
    TraceError(errMsg,result);
    return result;
  }
  printf("Write EnergyMatrix file %s\n",filepath);
  for(int i=0;i<pThis->designSiteCount;i++){
    for(int k=i;k<pThis->designSiteCount;k++){
      EnergyMatrixBlock* pEnergyBlock = EnergyMatrixGetBlock(pThis,i,k);
      for(int ij=0; ij<pEnergyBlock->RotamerCountSiteI; ij++){
        for(int ks=0; ks<pEnergyBlock->RotamerCountSiteK; ks++){
          if(i==k && ij!=ks) continue;
          fprintf(pFile,"%f  %d  %d  %d  %d\n",*EnergyMatrixBlockGet(pEnergyBlock,ij,ks),i,ij,k,ks);
        }
      }
    }
  }
  fclose(pFile);

  return Success;
}


EnergyMatrixBlock* EnergyMatrixGetBlock(EnergyMatrix* pThis,int designSiteI,int designSiteK){
  return &pThis->blocks[designSiteI*pThis->designSiteCount + designSiteK];
}


double* EnergyMatrixGet(EnergyMatrix* pThis,int designSiteI,int designSiteK,int rotamerIJ,int rotamerKS){
  if(designSiteI>designSiteK){
    int swap = designSiteI;
    designSiteI = designSiteK;
    designSiteK = swap;

    swap = rotamerIJ;
    rotamerIJ = rotamerKS;
    rotamerKS = swap;
  }
  return EnergyMatrixBlockGet(EnergyMatrixGetBlock(pThis,designSiteI,designSiteK),rotamerIJ,rotamerKS);
}


int  EnergyMatrixGetSiteCount(EnergyMatrix* pThis){
  return pThis->designSiteCount;
}


int  EnergyMatrixGetRotamerCount(EnergyMatrix* pThis,int designSiteI){
  return pThis->blocks[designSiteI*pThis->designSiteCount+designSiteI].RotamerCountSiteI;
}


int EnergyMatrixShow(EnergyMatrix* pThis){
  for(int i=0;i<pThis->designSiteCount;i++){
    for(int k=i;k<pThis->designSiteCount;k++){
      printf("energy block between site %d and site %d:\n",i,k);
      for(int j=0; j<EnergyMatrixGetRotamerCount(pThis,i);j++){
        for(int s=0; s<EnergyMatrixGetRotamerCount(pThis,k);s++){
          if(i==k && j!=s) continue;
          printf("%f ",*EnergyMatrixGet(pThis,i,k,j,s));
        }
        printf("\n");
      }
    }
  }
  return Success;
}


//'jobIndex' and 'totalJobCount' is designed for parallel computation,
//If not used in parallel computation environment, set both 'totalJobCount' and 'jobIndex' to be 1
int EnergyMatrixGenerate(Structure* pStructure,char* energyMatrixFilePath,int slotIndex,int slotCount){
  typedef struct{
    int computationAmount;
    int designSiteI;
    int designSiteK;
  } Job;

  typedef struct{
    Job* jobs;
    int jobCount;
    int totalComputationAmount;
  } Slot;

  int jobCount = pStructure->designSiteCount * (pStructure->designSiteCount+1) / 2;
  Job* jobs = (Job*)malloc(sizeof(Job)*jobCount);
  Slot* slots= (Slot*)malloc(sizeof(Slot)*slotCount);

  int i,j,counter;
  for(i=0;i<slotCount;i++){
    slots[i].jobCount = 0;
    slots[i].totalComputationAmount = 0;
    slots[i].jobs = NULL;
  }

  if(slotCount<=0 || slotIndex<=0 || slotIndex>slotCount){
    return IndexError;
  }
  slotIndex--;

  counter = 0;
  for(i=0;i<pStructure->designSiteCount;i++){
    for(j=i;j<pStructure->designSiteCount;j++){
      int rotamerCountOnSiteI = RotamerSetGetCount(&pStructure->designSites[i].rotamers);
      int rotamerCountOnSiteK = RotamerSetGetCount(&pStructure->designSites[j].rotamers);
      if(i==j){
        jobs[counter].computationAmount = rotamerCountOnSiteI;
      }
      else{
        jobs[counter].computationAmount = rotamerCountOnSiteI * rotamerCountOnSiteK;
      }
      jobs[counter].designSiteI = i;
      jobs[counter].designSiteK = j;
      counter++;
    }
  }

  //Sort 'jobs' by descending order ; Bubble sorting
  for(i=0;i<jobCount-1;i++){
    for(j=jobCount-1;j>i;j--){
      if(jobs[j].computationAmount > jobs[j-1].computationAmount){
        Job tempJob = jobs[j];
        jobs[j] = jobs[j-1];
        jobs[j-1] = tempJob;
      }
    }
  }

  //Assigning jobs to each slot
  for(i=0;i<jobCount;i++){
    Slot* pMinSlot = &slots[0];
    int minSlotAmount = slots[0].totalComputationAmount;
    for(j=0;j<slotCount;j++){
      if(slots[j].totalComputationAmount < minSlotAmount){
        pMinSlot = &slots[j];
        minSlotAmount = slots[j].totalComputationAmount;
      }
    }
    (pMinSlot->jobCount)++;
    pMinSlot->jobs = (Job*)realloc(pMinSlot->jobs,sizeof(Job)* pMinSlot->jobCount );
    pMinSlot->jobs[pMinSlot->jobCount-1] = jobs[i];
    pMinSlot->totalComputationAmount += jobs[i].computationAmount;
  }

  //Show the job assignment of each slot
  for(i=0;i<slotCount;i++){
    printf("Slot #%d, %8d ",i+1,slots[i].totalComputationAmount);
    for(j=0;j<slots[i].jobCount; j++){
      printf("(%2d,%2d) ",slots[i].jobs[j].designSiteI,slots[i].jobs[j].designSiteK);
    }
    printf("\n");
  }

  //Do the computational jobs assigned to current slot
  for(i=slots[slotIndex].jobCount-1;i>=0;i--){ 
    int result;
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    Job curJob = slots[slotIndex].jobs[i];
    EnergyMatrixBlock energyBlock;
    char energyBlockFileName[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(energyBlockFileName,"%s_%3.3d_%3.3d.txt",energyMatrixFilePath,curJob.designSiteI,curJob.designSiteK);
    printf("%s\n",energyBlockFileName);
    FILE* pFile = fopen(energyBlockFileName,"w");
    if(!pFile){
      sprintf(errMsg,"int file %s line %d, cannot write to file %s",__FILE__,__LINE__,energyBlockFileName);
      result = IOError;
      TraceError(errMsg,result);
      return result;        
    }

    EnergyMatrixBlockCreate(&energyBlock);
    result = EnergyMatrixBlockGenerate(&energyBlock,pStructure,
      curJob.designSiteI,curJob.designSiteK,pFile);
    if(FAILED(result)){
      return result;
    }

    EnergyMatrixBlockDestroy(&energyBlock);
    fclose(pFile);
  }

  for(i=0;i<slotCount;i++){
    free(slots[i].jobs);
  }
  free(jobs);
  free(slots);
  return Success;
}


//////////////////////////////////////////////////////////////////////////////////////////
//split the whole job of energy matrix computation into small jobs for parallel computing
//////////////////////////////////////////////////////////////////////////////////////////
int CompareJobAmount(const void *a, const void *b){
  double da = (*(Job*)a).computationAmount;
  double db = (*(Job*)b).computationAmount;
  return da < db ? -1 : 1;
}


int EnergyMatrixGenerateBasedOnPartition(Structure* pStructure, RotamerList* pList, char* energyMatrixFilePath,int slotIndex,int slotCount){
  int PARTITION_SIZE = PARTITION_SIZE_PARA;
  int jobCount;
  Job* jobs = NULL;
  Slot* slots= NULL;

  int i, j, k, s, counter;

  time_t totalTimeStart, totalTimeCurrent;
  int totalTimeElapsed;

  double cpuTimeStart, cpuTimeCurrent, cpuTimeElapsed;

  SitePartition* sitePartitions = (SitePartition*)malloc(sizeof(SitePartition)*StructureGetDesignSiteCount(pStructure));

  int * remainRotamerCount = (int*)malloc(sizeof(int)*pList->designSiteCount);
  double* energyInMemory = NULL;
  FILE* slotEnergyFile = NULL;
  char slotEnergyFileName[MAX_LENGTH_ONE_LINE_IN_FILE+1];

  for(i = 0; i < pList->designSiteCount; i++){
    remainRotamerCount[i] = 0;
    for(j = 0; j < pList->rotamerCount[i]; j++){
      if(pList->remainFlag[i][j] == TRUE) remainRotamerCount[i]++;
    }
  }
  for(i = 0; i < StructureGetDesignSiteCount(pStructure); i++){
    int rotamerIndex;
    DesignSite* pDesignSiteI = StructureGetDesignSite(pStructure, i);
    RotamerSet* pRotamerSetI = DesignSiteGetRotamers(pDesignSiteI);
    int partitionCountI = (int)(remainRotamerCount[i]/PARTITION_SIZE);
    if(remainRotamerCount[i]%PARTITION_SIZE != 0) partitionCountI++;
    printf("site: %d, partition count: %d\n", i, partitionCountI);
    sitePartitions[i].designSite = i;
    sitePartitions[i].partitionCount = partitionCountI;
    sitePartitions[i].rotamerCount = (int*)malloc(sizeof(int)*partitionCountI);
    sitePartitions[i].rotamerIndexOld = (IntArray*)malloc(sizeof(IntArray)*partitionCountI);
    sitePartitions[i].rotamerIndexNew = (IntArray*)malloc(sizeof(IntArray)*partitionCountI);

    for(j = 0; j < partitionCountI; j++){
      IntArrayCreate(&sitePartitions[i].rotamerIndexOld[j], 0);
      IntArrayCreate(&sitePartitions[i].rotamerIndexNew[j], 0);
      if(j == partitionCountI - 1){
        sitePartitions[i].rotamerCount[j] = remainRotamerCount[i] - j*PARTITION_SIZE;
        
      }
      else{
        sitePartitions[i].rotamerCount[j] = PARTITION_SIZE;
      }
      IntArrayResize(&sitePartitions[i].rotamerIndexOld[j], sitePartitions[i].rotamerCount[j]);
      IntArrayResize(&sitePartitions[i].rotamerIndexNew[j], sitePartitions[i].rotamerCount[j]);
    }

    rotamerIndex = 0;
    for(j = 0; j < RotamerSetGetCount(pRotamerSetI); j++){
      int partitionIndex;
      if(pList->remainFlag[i][j] == FALSE) continue;
      partitionIndex = rotamerIndex/PARTITION_SIZE;
      IntArraySet(&sitePartitions[i].rotamerIndexOld[partitionIndex], rotamerIndex - partitionIndex*PARTITION_SIZE, j);
      IntArraySet(&sitePartitions[i].rotamerIndexNew[partitionIndex], rotamerIndex - partitionIndex*PARTITION_SIZE, rotamerIndex);
      rotamerIndex++;
    }
  }

  jobCount = 0;
  for(i = 0; i < pStructure->designSiteCount; i++){
    int partitionCountOnSiteI = sitePartitions[i].partitionCount;
    for(k = i; k < pStructure->designSiteCount; k++){
      int partitionCountOnSiteK = sitePartitions[k].partitionCount;
      if(k == i){
        jobCount += partitionCountOnSiteI;
      }
      else{
        jobCount += partitionCountOnSiteI * partitionCountOnSiteK;
      }
    }
  }
  printf("jobCount: %d\n", jobCount);
  jobs = (Job*)malloc(sizeof(Job)*jobCount);
  slots= (Slot*)malloc(sizeof(Slot)*slotCount);


  for(i = 0; i < slotCount; i++){
    slots[i].jobCount = 0;
    slots[i].totalComputationAmount = 0.0;
    slots[i].jobs = NULL;
  }

  if(slotCount <= 0 || slotIndex <= 0 || slotIndex > slotCount){
    return ValueError;
  }
  slotIndex--;

  counter = 0;
  for(i = 0; i < pStructure->designSiteCount; i++){
    Chain* pChainI = StructureGetChain(pStructure, StructureGetDesignSite(pStructure, i)->chainIndex);
    int partitionCountOnSiteI = sitePartitions[i].partitionCount;
    for(k = i; k < pStructure->designSiteCount; k++){
      Chain* pChainK = StructureGetChain(pStructure, StructureGetDesignSite(pStructure, k)->chainIndex);
      int partitionCountOnSiteK = sitePartitions[k].partitionCount;
      for(j = 0; j < partitionCountOnSiteI; j++){
        if(k == i){
          jobs[counter].computationAmount = sitePartitions[i].rotamerCount[j];

          jobs[counter].designSiteI = i;
          jobs[counter].designSiteK = i;
          jobs[counter].rotamerPartitionIndexOnSiteI = j;
          jobs[counter].rotamerPartitionIndexOnSiteK = j;
          counter++;
        }
        else{
          for(s = 0; s < partitionCountOnSiteK; s++){
            jobs[counter].computationAmount = sitePartitions[i].rotamerCount[j] * sitePartitions[k].rotamerCount[s];
            
            jobs[counter].designSiteI = i;
            jobs[counter].designSiteK = k;
            jobs[counter].rotamerPartitionIndexOnSiteI = j;
            jobs[counter].rotamerPartitionIndexOnSiteK = s;
            counter++;
          }
        }
      }
    }
  }

  qsort(jobs, jobCount, sizeof(Job), CompareJobAmount);
  for(i = jobCount-1; i >= 0; i--){
    Slot* pMinSlot = &slots[0];
    double minSlotAmount = slots[0].totalComputationAmount;
    for(j = 0; j < slotCount; j++){
      if(slots[j].totalComputationAmount < minSlotAmount){
        pMinSlot = &slots[j];
        minSlotAmount = slots[j].totalComputationAmount;
      }
    }
    (pMinSlot->jobCount)++;
    pMinSlot->jobs = (Job*)realloc(pMinSlot->jobs,sizeof(Job)* pMinSlot->jobCount );
    pMinSlot->jobs[pMinSlot->jobCount-1] = jobs[i];
    pMinSlot->totalComputationAmount += jobs[i].computationAmount;
  }

  //Show the job assignment of each slot
  for(i=0;i<slotCount;i++){
    int index = 0;
    printf("Slot #%d, computational amount %.2f \n",i+1,slots[i].totalComputationAmount);
    for(j=0;j<slots[i].jobCount; j++){
      printf("(%2d,%2d,%4d,%4d) ",slots[i].jobs[j].designSiteI, slots[i].jobs[j].designSiteK,
        slots[i].jobs[j].rotamerPartitionIndexOnSiteI, slots[i].jobs[j].rotamerPartitionIndexOnSiteK);
      index++;
      if(index % 5 == 0){
        printf("\n");
      }
    }
    printf("\n");
  }

  //allocate memories to store the energy temporarily for each CPU, and finally output to files;
  counter = 0;
  for(i = slots[slotIndex].jobCount - 1; i >= 0; i--){
    int designSiteI;
    int designSiteK;
    int partitionIndexOnSiteI;
    int partitionIndexOnSiteK;
    Job curJob = slots[slotIndex].jobs[i];

    designSiteI = curJob.designSiteI;
    designSiteK = curJob.designSiteK;
    partitionIndexOnSiteI = curJob.rotamerPartitionIndexOnSiteI;
    partitionIndexOnSiteK = curJob.rotamerPartitionIndexOnSiteK;

    counter += sitePartitions[designSiteI].rotamerCount[partitionIndexOnSiteI]*sitePartitions[designSiteK].rotamerCount[partitionIndexOnSiteK];
  }

  energyInMemory = (double*)malloc(sizeof(double)*counter);
  for(i = 0; i < counter; i++){
    energyInMemory[i] = 0.0;
  }

  //Do the computational jobs assigned to current slot
  totalTimeStart = time(NULL);
  cpuTimeStart = clock();
  counter = 0;
  for(i = slots[slotIndex].jobCount - 1; i >= 0; i--){
    int result;
    int designSiteI;
    int designSiteK;
    int partitionIndexOnSiteI;
    int partitionIndexOnSiteK;
    time_t begin;
    time_t end;
    int elapsed;
    Job curJob = slots[slotIndex].jobs[i];
    char energyBlockFileName[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    begin = time(NULL);
    sprintf(energyBlockFileName,"%s_%3.3d_%3.3d_%4.4d_%4.4d.txt",energyMatrixFilePath, 
      curJob.designSiteI, curJob.designSiteK, curJob.rotamerPartitionIndexOnSiteI, curJob.rotamerPartitionIndexOnSiteK);
    printf("%s, ",energyBlockFileName);

    designSiteI = curJob.designSiteI;
    designSiteK = curJob.designSiteK;
    partitionIndexOnSiteI = curJob.rotamerPartitionIndexOnSiteI;
    partitionIndexOnSiteK = curJob.rotamerPartitionIndexOnSiteK;

    result = EnergyMatrixPartitionPairGenerateStoreInMemory(pStructure, designSiteI, designSiteK, 
      &sitePartitions[designSiteI].rotamerIndexOld[partitionIndexOnSiteI], &sitePartitions[designSiteK].rotamerIndexOld[partitionIndexOnSiteK], 
      &sitePartitions[designSiteI].rotamerIndexNew[partitionIndexOnSiteI], &sitePartitions[designSiteK].rotamerIndexNew[partitionIndexOnSiteK],
      energyInMemory, counter);
    counter += sitePartitions[designSiteI].rotamerCount[partitionIndexOnSiteI]*sitePartitions[designSiteK].rotamerCount[partitionIndexOnSiteK];
    if(FAILED(result)){
      return result;
    }

    end = time(NULL);
    elapsed = (int)(end - begin);
    printf("elapsed time: %d hr %d min %d sec\n", (int)(elapsed/3600), (int)(elapsed%3600/60), (int)(elapsed%3600%60));
  }

  // output energy values from memories to files;
  counter = 0;
  sprintf(slotEnergyFileName, "%s_CPU%03d.txt", energyMatrixFilePath, slotIndex+1);
  slotEnergyFile = fopen(slotEnergyFileName, "w");
  if(slotEnergyFile == NULL){
    printf("cannot open file %s for writing\n", slotEnergyFileName);
    return IOError;
  }
  for(i = slots[slotIndex].jobCount - 1; i >= 0; i--){
    int result;
    int designSiteI;
    int designSiteK;
    int partitionIndexOnSiteI;
    int partitionIndexOnSiteK;
    Job curJob = slots[slotIndex].jobs[i];

    designSiteI = curJob.designSiteI;
    designSiteK = curJob.designSiteK;
    partitionIndexOnSiteI = curJob.rotamerPartitionIndexOnSiteI;
    partitionIndexOnSiteK = curJob.rotamerPartitionIndexOnSiteK;

    result = EnergyMatrixPartitionPairWriteEnergyToFile(pStructure, designSiteI, designSiteK, 
      &sitePartitions[designSiteI].rotamerIndexOld[partitionIndexOnSiteI], &sitePartitions[designSiteK].rotamerIndexOld[partitionIndexOnSiteK], 
      &sitePartitions[designSiteI].rotamerIndexNew[partitionIndexOnSiteI], &sitePartitions[designSiteK].rotamerIndexNew[partitionIndexOnSiteK],
      energyInMemory, counter, slotEnergyFile);
    counter += sitePartitions[designSiteI].rotamerCount[partitionIndexOnSiteI]*sitePartitions[designSiteK].rotamerCount[partitionIndexOnSiteK];
    if(FAILED(result)){
      return result;
    }
  }
  fclose(slotEnergyFile);

  totalTimeCurrent = time(NULL);
  totalTimeElapsed = (int)(totalTimeCurrent - totalTimeStart);
  printf("elapsed time for calculating EnergyMatrixBlocks: %d hr %d min %d sec\n", totalTimeElapsed/3600, totalTimeElapsed%3600/60, totalTimeElapsed%3600%60);

  cpuTimeCurrent = clock();
  cpuTimeElapsed = cpuTimeCurrent - cpuTimeStart;
  printf("elapsed CPU time for calculating EnergyMatrixBlocks: %f secs\n", cpuTimeElapsed/CLOCKS_PER_SEC);


  for(i = 0; i < pStructure->designSiteCount; i++){
    free(sitePartitions[i].rotamerCount);
    for(j = 0; j < sitePartitions[i].partitionCount; j++){
      IntArrayDestroy(&sitePartitions[i].rotamerIndexOld[j]);
      IntArrayDestroy(&sitePartitions[i].rotamerIndexNew[j]);
    }
  }
  free(sitePartitions);

  for(i=0;i<slotCount;i++){
    free(slots[i].jobs);
  }
  free(jobs);
  free(slots);
  free(remainRotamerCount);
  free(energyInMemory);

  return Success;
}


int EnergyMatrixPartitionPairGenerate(Structure* pStructure, int designSiteI, int designSiteK, 
                     IntArray* pRotamerIndexOldOnSiteI, IntArray* pRotamerIndexOldOnSiteK, 
                     IntArray* pRotamerIndexNewOnSiteI, IntArray* pRotamerIndexNewOnSiteK, FILE* outputFile)
{
  DesignSite *pDesignSiteI = StructureGetDesignSite(pStructure, designSiteI);
  DesignSite *pDesignSiteK = StructureGetDesignSite(pStructure, designSiteK);
  Chain *pChainI = StructureGetChain(pStructure, pDesignSiteI->chainIndex);
  Chain* pChainK = StructureGetChain(pStructure, pDesignSiteK->chainIndex);
  double* partitionPairEnergy = NULL;

  if(designSiteI > designSiteK || pDesignSiteI == NULL || pDesignSiteK == NULL){
    return ValueError;
  }

  if(designSiteI == designSiteK){
    partitionPairEnergy = (double*)malloc(sizeof(double)*IntArrayGetLength(pRotamerIndexOldOnSiteI));
    for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
      Rotamer *pRotamerIJ  = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), IntArrayGet(pRotamerIndexOldOnSiteI, j));
      RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));

      double total = 0;
      double energyTerms[MAX_ENERGY_TERM]={0};
      if(ChainGetType(pChainI) == Type_Chain_Protein){
        AminoAcidReferenceEnergy(RotamerGetType(pRotamerIJ), energyTerms);
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(b==pDesignSiteI->resiIndex){
                EnergyIntraRotamer(pRotamerIJ, energyTerms);
              }
              else{
                if(ChainGetType(pChainA) != Type_Chain_SmallMol){
                  if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                    EnergyRotamerAndFixedResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                  else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                    pResidueAB->designType == Type_ResidueDesignType_Designable){
                      EnergyRotamerAndDesignResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                }
              }
            }
          }
          else{//different chains
            if(ChainGetType(pChainA) != Type_Chain_SmallMol){
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                  EnergyRotamerAndFixedResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                  pResidueAB->designType == Type_ResidueDesignType_Designable){
                    EnergyRotamerAndDesignResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
          }
        }
      }
      else if(ChainGetType(pChainI) == Type_Chain_SmallMol){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){
            total += pRotamerIJ->vdwInternal;
          }
          else{
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                EnergyRotamerAndFixedResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
              }
              else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                pResidueAB->designType == Type_ResidueDesignType_Designable){
                  EnergyRotamerAndDesignResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
              }
            }
          }
        }
      }

      EnergyTermWeighting(energyTerms);
      total += energyTerms[0];
      partitionPairEnergy[j] = total;
      RotamerExtract(pRotamerIJ);
    }

    for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
      fprintf(outputFile,"%f %d %d %d %d\n",partitionPairEnergy[j],designSiteI,IntArrayGet(pRotamerIndexNewOnSiteI, j),designSiteI,IntArrayGet(pRotamerIndexNewOnSiteI, j));
    }
    free(partitionPairEnergy);
    return Success;
  }

  partitionPairEnergy = (double*)malloc(sizeof(double)*IntArrayGetLength(pRotamerIndexOldOnSiteI)*IntArrayGetLength(pRotamerIndexOldOnSiteK));
  for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
    Rotamer *pRotamerIJ = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), IntArrayGet(pRotamerIndexOldOnSiteI, j));
    RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));
    for(int s = 0; s < IntArrayGetLength(pRotamerIndexOldOnSiteK); s++){
      Rotamer *pRotamerKS = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteK), IntArrayGet(pRotamerIndexOldOnSiteK, s));
      RotamerRestore(pRotamerKS,DesignSiteGetRotamers(pDesignSiteK));

      double energyTerms[MAX_ENERGY_TERM]={0};
      if(pDesignSiteI->chainIndex == pDesignSiteK->chainIndex){
        EnergyRotamerAndRotamerSameChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      else{
        EnergyRotamerAndRotamerDiffChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      EnergyTermWeighting(energyTerms);
      double total=0;
      for(int a=1; a<MAX_ENERGY_TERM;a++){
        total += energyTerms[a];
      }
      partitionPairEnergy[(j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+s)] = total;
      RotamerExtract(pRotamerKS);
    }
    RotamerExtract(pRotamerIJ);
  }

  for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
    for(int s = 0; s < IntArrayGetLength(pRotamerIndexOldOnSiteK); s++){
      fprintf(outputFile,"%f %d %d %d %d\n", 
        partitionPairEnergy[(j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+s)], 
        designSiteI, 
        IntArrayGet(pRotamerIndexNewOnSiteI, j), 
        designSiteK, 
        IntArrayGet(pRotamerIndexNewOnSiteK, j));
    }
  }
  free(partitionPairEnergy);

  return Success;
}


int EnergyMatrixPartitionPairGenerateStoreInMemory(Structure* pStructure, int designSiteI, int designSiteK, 
                    IntArray* pRotamerIndexOldOnSiteI, IntArray* pRotamerIndexOldOnSiteK, 
                    IntArray* pRotamerIndexNewOnSiteI, IntArray* pRotamerIndexNewOnSiteK, 
                    double* energyInMemory, int startIndexThisPartition)
{
  int j, s;
  double    solvNonPolar,vdwAttract,vdwRepul, hbond, elecDesolvation, elecScreenedCoulomb, entropy, total;
  DesignSite *pDesignSiteI = StructureGetDesignSite(pStructure, designSiteI);
  DesignSite *pDesignSiteK = StructureGetDesignSite(pStructure, designSiteK);
  Chain* pChainI = StructureFindChainByName(pStructure, ResidueGetChainName(pDesignSiteI->pResidue));
  Chain* pChainK = StructureFindChainByName(pStructure, ResidueGetChainName(pDesignSiteK->pResidue));
  if(designSiteI > designSiteK || pDesignSiteI == NULL || pDesignSiteK == NULL){
    return ValueError;
  }

  if(designSiteI == designSiteK){
    for(j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
      Rotamer *pRotamerIJ  = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), IntArrayGet(pRotamerIndexOldOnSiteI, j));
      RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));

      double total = 0;
      double energyTerms[MAX_ENERGY_TERM]={0};
      Chain *pChainI = StructureGetChain(pStructure, pDesignSiteI->chainIndex);
      if(ChainGetType(pChainI) == Type_Chain_Protein){
        AminoAcidReferenceEnergy(RotamerGetType(pRotamerIJ), energyTerms);
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(b==pDesignSiteI->resiIndex){
                EnergyIntraRotamer(pRotamerIJ, energyTerms);
              }
              else{
                if(ChainGetType(pChainA) != Type_Chain_SmallMol){
                  if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                    EnergyRotamerAndFixedResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                  else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                    pResidueAB->designType == Type_ResidueDesignType_Designable){
                      EnergyRotamerAndDesignResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                }
              }
            }
          }
          else{//different chains
            if(ChainGetType(pChainA) != Type_Chain_SmallMol){
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                  EnergyRotamerAndFixedResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                  pResidueAB->designType == Type_ResidueDesignType_Designable){
                    EnergyRotamerAndDesignResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
          }
        }
      }
      else if(ChainGetType(pChainI) == Type_Chain_SmallMol){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){
            total += pRotamerIJ->vdwInternal;
          }
          else{
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                EnergyRotamerAndFixedResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
              }
              else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                pResidueAB->designType == Type_ResidueDesignType_Designable){
                  EnergyRotamerAndDesignResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
              }
            }
          }
        }
      }

      EnergyTermWeighting(energyTerms);
      total += energyTerms[0];
      energyInMemory[startIndexThisPartition + (j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+j)] = total;
      RotamerExtract(pRotamerIJ);
    }

    return Success;
  }

  for(j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
    Rotamer *pRotamerIJ        = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteI), IntArrayGet(pRotamerIndexOldOnSiteI, j));
    RotamerRestore(pRotamerIJ,DesignSiteGetRotamers(pDesignSiteI));
    for(s = 0; s < IntArrayGetLength(pRotamerIndexOldOnSiteK); s++){
      Rotamer *pRotamerKS        = RotamerSetGet(DesignSiteGetRotamers(pDesignSiteK), IntArrayGet(pRotamerIndexOldOnSiteK, s));
      RotamerRestore(pRotamerKS,DesignSiteGetRotamers(pDesignSiteK));

      double energyTerms[MAX_ENERGY_TERM]={0};
      if(pDesignSiteI->chainIndex == pDesignSiteK->chainIndex){
        EnergyRotamerAndRotamerSameChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      else{
        EnergyRotamerAndRotamerDiffChain(pRotamerIJ, pRotamerKS, energyTerms);
      }
      EnergyTermWeighting(energyTerms);
      double total=0;
      for(int a=1; a<MAX_ENERGY_TERM;a++){
        total += energyTerms[a];
      }

      energyInMemory[startIndexThisPartition+(j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+s)] = total;
      RotamerExtract(pRotamerKS);
    }
    RotamerExtract(pRotamerIJ);
  }
  return Success;
}


int EnergyMatrixPartitionPairWriteEnergyToFile(Structure* pStructure, int designSiteI, int designSiteK, 
                    IntArray* pRotamerIndexOldOnSiteI, IntArray* pRotamerIndexOldOnSiteK, 
                    IntArray* pRotamerIndexNewOnSiteI, IntArray* pRotamerIndexNewOnSiteK, 
                    double* energyInMemory, int startIndexThisPartition,
                    FILE* outputFile)
{
  if(designSiteI == designSiteK){
    for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
      double total = energyInMemory[startIndexThisPartition + (j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+j)];
      fprintf(outputFile,"%f %d %d %d %d\n", total, designSiteI, IntArrayGet(pRotamerIndexNewOnSiteI, j), designSiteI, IntArrayGet(pRotamerIndexNewOnSiteI, j));
    }

    return Success;
  }

  for(int j = 0; j < IntArrayGetLength(pRotamerIndexOldOnSiteI); j++){
    for(int s = 0; s < IntArrayGetLength(pRotamerIndexOldOnSiteK); s++){
      double total = energyInMemory[startIndexThisPartition + (j*IntArrayGetLength(pRotamerIndexOldOnSiteK)+s)];
      fprintf(outputFile,"%f %d %d %d %d\n", total,designSiteI, IntArrayGet(pRotamerIndexNewOnSiteI, j), designSiteK, IntArrayGet(pRotamerIndexNewOnSiteK, s));
    }
  }

  return Success;
}


int SelfEnergyGenerate(Structure *pStructure,char* filepath){
  int result=Success;
  char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  FILE* outputFile = fopen(filepath,"w");
  if(outputFile==NULL){
    sprintf(errMsg,"in file %s line %d, cannot create file %s to write",__FILE__,__LINE__,filepath);
    result = IOError;
    TraceError(errMsg,result);
    return result;        
  }

  for(int i=0; i<StructureGetDesignSiteCount(pStructure); i++){
    DesignSite *pDesignSiteI = StructureGetDesignSite(pStructure,i);
    RotamerSet *pSetI = DesignSiteGetRotamers(pDesignSiteI);
    Chain *pChainI = StructureGetChain(pStructure, pDesignSiteI->chainIndex);
    for(int j=0; j<RotamerSetGetCount(pSetI); j++){
      Rotamer *pRotamerIJ  = RotamerSetGet(pSetI, j);
      RotamerRestore(pRotamerIJ,pSetI);
      double energyTerms[MAX_ENERGY_TERM]={0};
      if(ChainGetType(pChainI) == Type_Chain_Protein){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(b==pDesignSiteI->resiIndex){
                AminoAcidReferenceEnergy(RotamerGetType(pRotamerIJ), energyTerms);
                EnergyIntraRotamer(pRotamerIJ, energyTerms);
              }
              else{
                if(ChainGetType(pChainA) != Type_Chain_SmallMol){
                  if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                    EnergyRotamerAndFixedResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                  else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                    pResidueAB->designType == Type_ResidueDesignType_Designable ||
                    pResidueAB->designType == Type_ResidueDesignType_Catalytic){
                      EnergyRotamerAndDesignResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                }
              }
            }
          }
          else{//different chains
            if(ChainGetType(pChainA) != Type_Chain_SmallMol){
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                  EnergyRotamerAndFixedResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                  pResidueAB->designType == Type_ResidueDesignType_Designable ||
                  pResidueAB->designType == Type_ResidueDesignType_Catalytic){
                    EnergyRotamerAndDesignResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
            else{//fixed small molecule
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                  EnergyRotamerAndFixedLigResidue(pRotamerIJ,pResidueAB,energyTerms);
                }
                else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                  pResidueAB->designType == Type_ResidueDesignType_Designable ||
                  pResidueAB->designType == Type_ResidueDesignType_Catalytic){
                    EnergyRotamerAndDesignResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
                }
              }
            }
          }
        }
      }

      else if(ChainGetType(pChainI) == Type_Chain_SmallMol){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            pRotamerIJ->selfenergy += pRotamerIJ->vdwInternal;
          }
          else{//different chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                EnergyLigRotamerAndFixedResidue(pRotamerIJ,pResidueAB,energyTerms);
              }
              else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                pResidueAB->designType == Type_ResidueDesignType_Designable ||
                pResidueAB->designType == Type_ResidueDesignType_Catalytic){
                  EnergyLigRotamerAndDesignResidue(pRotamerIJ,pResidueAB,energyTerms);
              }
            }
          }
        }
      }
      EnergyTermWeighting(energyTerms);
      pRotamerIJ->selfenergy += energyTerms[0];
      fprintf(outputFile, "%d %d %f\n",i,j,pRotamerIJ->selfenergy);
      RotamerExtract(pRotamerIJ);
    }
  }

  if(FAILED(result)){
    return result;
  }
  fclose(outputFile);

  return Success;
}


int SelfEnergyGenerate2(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRamaTable,char* filepath){
  FILE* outputFile = fopen(filepath,"w");
  if(outputFile==NULL){
    char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(errMsg,"in file %s line %d, cannot write to file %s",__FILE__,__LINE__,filepath);
    TraceError(errMsg,IOError);
    return IOError;
  }

  for(int i=0; i<StructureGetDesignSiteCount(pStructure); i++){
    DesignSite *pDesignSiteI = StructureGetDesignSite(pStructure, i);
    RotamerSet *pSetI = DesignSiteGetRotamers(pDesignSiteI);
    Chain *pChainI = StructureGetChain(pStructure, pDesignSiteI->chainIndex);
    for(int j=0; j<RotamerSetGetCount(pSetI); j++){
      Rotamer *pRotamerIJ  = RotamerSetGet(pSetI, j);
      RotamerRestore(pRotamerIJ,pSetI);
      double energyTerms[MAX_ENERGY_TERM]={0};
      double energyTermsBind[MAX_ENERGY_TERM]={0};
      if(ChainGetType(pChainI)==Type_Chain_Protein || ChainGetType(pChainI)==Type_Chain_DNA || ChainGetType(pChainI)==Type_Chain_RNA){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(b==pDesignSiteI->resiIndex){
                AminoAcidReferenceEnergy(RotamerGetType(pRotamerIJ),energyTerms);
                EnergyIntraRotamer(pRotamerIJ,energyTerms);
                RotamerPropensityAndRamachandranEnergy(pRotamerIJ,pResidueAB,pAAppTable,pRamaTable,energyTerms);
                RotamerDunbrackEnergy(pRotamerIJ,energyTerms);
              }
              else{
                if(ChainGetType(pChainA)==Type_Chain_Protein || ChainGetType(pChainA)==Type_Chain_DNA || ChainGetType(pChainA)==Type_Chain_RNA){
                  if(pResidueAB->designType==Type_ResidueDesignType_Fixed){
                    EnergyRotamerAndFixedResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                  else if(pResidueAB->designType==Type_ResidueDesignType_Repacked ||
                    pResidueAB->designType==Type_ResidueDesignType_Designable ||
                    pResidueAB->designType==Type_ResidueDesignType_Catalytic){
                      EnergyRotamerAndDesignResidueSameChain(pRotamerIJ,pResidueAB,energyTerms);
                  }
                }
              }
            }
          }
          else{//different chains
            if(ChainGetType(pChainA)==Type_Chain_Protein || ChainGetType(pChainA)==Type_Chain_DNA || ChainGetType(pChainA)==Type_Chain_RNA){
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                  EnergyRotamerAndFixedResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
                  //note: this code can easily cause a bug!!! (2/2/2020)
                  if((strstr(DES_CHAINS,ChainGetName(pChainA))==NULL && strstr(DES_CHAINS,ChainGetName(pChainI))!=NULL) ||
                    (strstr(DES_CHAINS,ChainGetName(pChainA))!=NULL && strstr(DES_CHAINS,ChainGetName(pChainI))==NULL)){
                    if(FLAG_PPI){
                      EnergyRotamerAndFixedResidueDiffChain(pRotamerIJ,pResidueAB,energyTermsBind);
                    }
                  }
                }
                else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                  pResidueAB->designType == Type_ResidueDesignType_Designable ||
                  pResidueAB->designType == Type_ResidueDesignType_Catalytic){
                    EnergyRotamerAndDesignResidueDiffChain(pRotamerIJ,pResidueAB,energyTerms);
                    //note: this code can easily cause a bug!!! (2/2/2020)
                    if((strstr(DES_CHAINS,ChainGetName(pChainA))==NULL && strstr(DES_CHAINS,ChainGetName(pChainI))!=NULL) ||
                      (strstr(DES_CHAINS,ChainGetName(pChainA))!=NULL && strstr(DES_CHAINS,ChainGetName(pChainI))==NULL)){
                      if(FLAG_PPI){
                        EnergyRotamerAndDesignResidueDiffChain(pRotamerIJ,pResidueAB,energyTermsBind);
                      }
                    }
                }
              }
            }
            else if(ChainGetType(pChainA)==Type_Chain_SmallMol){//pChainA is small molecule
              for(int b=0; b<ChainGetResidueCount(pChainA); b++){
                Residue* pResidueAB = ChainGetResidue(pChainA, b);
                if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                  EnergyRotamerAndFixedLigResidue(pRotamerIJ,pResidueAB,energyTerms);
                  if(FLAG_PROT_LIG || FLAG_ENZYME){
                    EnergyRotamerAndFixedLigResidue(pRotamerIJ,pResidueAB,energyTermsBind);
                  }
                }
              }
            }
          }
        }
      }
      else if(ChainGetType(pChainI)==Type_Chain_SmallMol){
        for(int a=0; a<StructureGetChainCount(pStructure); a++){
          Chain* pChainA = StructureGetChain(pStructure, a);
          if(a == pDesignSiteI->chainIndex){//same chain
            pRotamerIJ->selfenergy += pRotamerIJ->vdwInternal;
            //note do not count vdwBackbone because it will be re-calculated
            //pRotamerIJ->selfenergyBind += pRotamerIJ->vdwBackbone;
          }
          else{//different chain
            for(int b=0; b<ChainGetResidueCount(pChainA); b++){
              Residue* pResidueAB = ChainGetResidue(pChainA, b);
              if(pResidueAB->designType == Type_ResidueDesignType_Fixed){
                EnergyLigRotamerAndFixedResidue(pRotamerIJ,pResidueAB,energyTerms);
                EnergyLigRotamerAndFixedResidue(pRotamerIJ,pResidueAB,energyTermsBind);
              }
              else if(pResidueAB->designType == Type_ResidueDesignType_Repacked ||
                pResidueAB->designType == Type_ResidueDesignType_Designable ||
                pResidueAB->designType == Type_ResidueDesignType_Catalytic){
                  EnergyLigRotamerAndDesignResidue(pRotamerIJ,pResidueAB,energyTerms);
                  EnergyLigRotamerAndDesignResidue(pRotamerIJ,pResidueAB,energyTermsBind);
              }
            }
          }
        }
      }
      EnergyTermWeighting(energyTerms);
      pRotamerIJ->selfenergy += energyTerms[0];
      EnergyTermWeighting(energyTermsBind);
      pRotamerIJ->selfenergyBind += energyTermsBind[0];
      fprintf(outputFile, "%d %d %f %f\n", i, j, pRotamerIJ->selfenergy,pRotamerIJ->selfenergyBind);
      RotamerExtract(pRotamerIJ);
    }
  }
  fclose(outputFile);

  return Success;
}


int SelfEnergyReadAndCheck(Structure* pStructure,RotamerList* pList,char* filepath){
  FILE* pFile = fopen(filepath, "r");
  if(pFile == NULL){
    char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(errMsg,"in file %s line %d, cannot write to file %s",__FILE__,__LINE__,filepath);
    TraceError(errMsg,IOError);
    return IOError;
  }

  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(fgets(buffer, MAX_LENGTH_ONE_LINE_IN_FILE, pFile)){
    int designSiteI=-1;
    int rotamerIJ=-1;
    double selfTot=1000.0;
    double selfBin=1000.0;
    sscanf(buffer, "%d %d %lf %lf\n",&designSiteI,&rotamerIJ,&selfTot,&selfBin);
    DesignSite* pSite=StructureGetDesignSite(pStructure,designSiteI);
    RotamerSet* pSet=DesignSiteGetRotamers(pSite);
    Rotamer* pRotamer=RotamerSetGet(pSet,rotamerIJ);
    pRotamer->selfenergy=selfTot;
    pRotamer->selfenergyBind=selfBin;
  }

  for(int i = 0; i < StructureGetDesignSiteCount(pStructure); i++){
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure, i);
    RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
    for(int j = 0; j < RotamerSetGetCount(pRotamerSet); j++){
      Rotamer* pRotamerJ = RotamerSetGet(pRotamerSet, j);
      if(pList->remainFlag[i][j] == FALSE) continue;
      for(int k = 0; k < RotamerSetGetCount(pRotamerSet); k++){
        Rotamer* pRotamerK = RotamerSetGet(pRotamerSet, k);
        if(pList->remainFlag[i][k] == FALSE) continue;
        if((pRotamerJ->selfenergy - pRotamerK->selfenergy > SELF_ENERGY_DIFFERENT_THRESHOLD && RotamerAndRotamerInSameType(pRotamerJ,pRotamerK)==FALSE) || 
          (pRotamerJ->selfenergy - pRotamerK->selfenergy > SELF_ENERGY_SAME_THRESHOLD && RotamerAndRotamerInSameType(pRotamerJ,pRotamerK))){
            pList->remainFlag[i][j] = FALSE;
            break;
        }
      }
    }
  }
  fclose(pFile);

  return Success;
}


int RotamerListCreateFromStructure(RotamerList* pThis,Structure* pStructure){
  pThis->designSiteCount = StructureGetDesignSiteCount(pStructure);
  pThis->rotamerCount = (int*)malloc(sizeof(int)*pThis->designSiteCount);
  pThis->remainFlag = (BOOL**)malloc(sizeof(BOOL*)*pThis->designSiteCount);
  for(int i=0;i<pThis->designSiteCount;i++){
    pThis->rotamerCount[i] = RotamerSetGetCount(DesignSiteGetRotamers(StructureGetDesignSite(pStructure,i)));
    pThis->remainFlag[i] = (BOOL*)malloc(sizeof(BOOL)*pThis->rotamerCount[i]);
    for(int j=0;j<pThis->rotamerCount[i];j++){
      pThis->remainFlag[i][j] = TRUE;
    }
  }
  pThis->remainRotamerCount = (int*)malloc(sizeof(int)*pThis->designSiteCount);
  for(int i=0;i<pThis->designSiteCount;i++){
    pThis->remainRotamerCount[i]=0;
    for(int j=0;j<pThis->rotamerCount[i];j++){
      if(pThis->remainFlag[i][j]){
        pThis->remainRotamerCount[i]++;
      }
    }
  }
  return Success;
}

int RotamerListCreateFromEnergyMatrix(RotamerList* pThis,EnergyMatrix* pEnergyMatrix){
  pThis->designSiteCount = pEnergyMatrix->designSiteCount;
  pThis->rotamerCount = (int*)malloc(sizeof(int)*EnergyMatrixGetSiteCount(pEnergyMatrix));
  pThis->remainFlag = (BOOL**)malloc(sizeof(BOOL*)*EnergyMatrixGetSiteCount(pEnergyMatrix));
  for(int i=0;i<pThis->designSiteCount;i++){
    pThis->rotamerCount[i] = EnergyMatrixGetRotamerCount(pEnergyMatrix,i);
    pThis->remainFlag[i] = (BOOL*)malloc(sizeof(BOOL)*pThis->rotamerCount[i]);
    for(int j=0;j<pThis->rotamerCount[i];j++){
      pThis->remainFlag[i][j] = TRUE;
    }
  }
  pThis->rotamerCount = (int*)malloc(sizeof(int)*pThis->designSiteCount);
  for(int i=0;i<pThis->designSiteCount;i++){
    pThis->remainRotamerCount[i]=0;
    for(int j=0;j<pThis->rotamerCount[i];j++){
      if(pThis->remainFlag[i][j]){
        pThis->remainRotamerCount[i]++;
      }
    }
  }
  return Success;
}

void RotamerListDestroy(RotamerList* pThis){
  for(int i=0;i<pThis->designSiteCount;i++) free(pThis->remainFlag[i]);
  free(pThis->remainFlag);
  free(pThis->rotamerCount);
  free(pThis->remainRotamerCount);
}


int RotamerListCopy(RotamerList* pThis, RotamerList* pOther){
  RotamerListDestroy(pThis);
  pThis->designSiteCount=pOther->designSiteCount;
  pThis->rotamerCount=(int*)malloc(pThis->designSiteCount*sizeof(int));
  memcpy(pThis->rotamerCount,pOther->rotamerCount,sizeof(int)*pThis->designSiteCount);
  pThis->remainFlag=(BOOL**)malloc(pThis->designSiteCount*sizeof(BOOL*));
  for(int i=0; i<pThis->designSiteCount; i++){
    pThis->remainFlag[i]=(BOOL*)malloc(pThis->rotamerCount[i]*sizeof(BOOL));
    memcpy(pThis->remainFlag[i],pOther->remainFlag[i],sizeof(BOOL)*pThis->rotamerCount[i]);
  }
  pThis->remainRotamerCount=(int*)malloc(pThis->designSiteCount*sizeof(int));
  memcpy(pThis->remainRotamerCount,pOther->remainRotamerCount,sizeof(int)*pThis->designSiteCount);
  return Success;
}


int RotamerListRead(RotamerList* pThis,char* filepath){
  int result;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  FILE* pFile = fopen(filepath,"r");
  if(pFile==NULL){
    sprintf(errMsg,"in file %s line %d, cannot read file %s",__FILE__,__LINE__,filepath);
    result = IOError;
    TraceError(errMsg,result);
    return result;
  }

  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    int i,j;
    char flag[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sscanf(line,"%d %d %s",&i,&j,flag);
    if(i>=pThis->designSiteCount){
      sprintf(errMsg,"in file %s line %d, structure does not have site #%d at this line:\n%s",__FILE__,__LINE__,i,line);
      result = ValueError;
      TraceError(errMsg,result);
      return result;
    }
    if(j>=pThis->rotamerCount[i]){
      sprintf(errMsg,"in file %s line %d, structure does not have rotamer #%d on site #%d at this line:\n%s",__FILE__,__LINE__,j,i,line);
      result = ValueError;
      TraceError(errMsg,result);
      return result;
    }
    if(strcmp(flag,"TRUE")==0) pThis->remainFlag[i][j] = TRUE;
    else if(strcmp(flag,"FALSE")==0) pThis->remainFlag[i][j] = FALSE;
    else{
      sprintf(errMsg,"in file %s line %d, incorrect file format at this line:\n%s",__FILE__,__LINE__,line);
      result = FormatError;
      TraceError(errMsg,result);
      return result;
    }
  }

  for(int i=0;i<pThis->designSiteCount;i++){
    pThis->remainRotamerCount[i]=0;
    for(int j=0;j<pThis->rotamerCount[i];j++){
      if(pThis->remainFlag[i][j]){
        pThis->remainRotamerCount[i]++;
      }
    }
  }
  return Success;
}


int RotamerListWrite(RotamerList* pThis,char* filepath){
  int result;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  FILE* pFile = fopen(filepath,"w");
  if(!pFile){
    result = IOError;
    sprintf(errMsg,"in file %s line %d, cannot read file %s",__FILE__,__LINE__,filepath);
    TraceError(errMsg,result);
    return result;
  }
  for(int i=0;i<pThis->designSiteCount;i++){
    for(int j=0;j<pThis->rotamerCount[i];j++){
      fprintf(pFile,"%d\t%d\t%s\n",i,j,  
        pThis->remainFlag[i][j]? "TRUE" : "FALSE"  );
    }
  }
  fclose(pFile);
  return Success;
}



int RotamerListShow(RotamerList* pThis){
  int totalRotamerCount=0;
  for(int i=0;i<pThis->designSiteCount;i++){
    int remainRotamerCount=0;
    printf("site %d, %d rotamers:\n",i,pThis->rotamerCount[i]);
    for(int j=0;j<pThis->rotamerCount[i];j++){
      if(pThis->remainFlag[i][j]){
        printf("rotamer %d\n",j);
        remainRotamerCount++;
      }
    }
    totalRotamerCount+=remainRotamerCount;
    printf("%d rotamers left\n",remainRotamerCount);
  }
  printf("total rotamer count:\t%d\n",totalRotamerCount);
  return Success;
}


int RotamerListAndEnergyMatrixDelete(RotamerList* pList,EnergyMatrix* pMatrix,IntArray* pDeletedRotamers){
  int result;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  BOOL** rotamerDeletedFlag = (BOOL**)malloc(sizeof(BOOL*)*pList->designSiteCount);

  for(int i=0;i<pList->designSiteCount;i++){
    rotamerDeletedFlag[i] = (BOOL*)malloc(sizeof(BOOL)*EnergyMatrixGetBlock(pMatrix,i,i)->RotamerCountSiteI);
    for(int j=0;j<EnergyMatrixGetBlock(pMatrix,i,i)->RotamerCountSiteI;j++){
      rotamerDeletedFlag[i][j] = FALSE;
    }
  }

  if(IntArrayGetLength(pDeletedRotamers)%2!=0){
    result = ValueError;
    sprintf(errMsg,"in file %s line %d, parameter #3 must be a int array with odd number length",__FILE__,__LINE__);
    TraceError(errMsg,result);
    return result;
  }
  for(int i=0;i<IntArrayGetLength(pDeletedRotamers);i+=2){
    int siteIndex = IntArrayGet(pDeletedRotamers,i);
    int rotamerIndex = IntArrayGet(pDeletedRotamers,i+1);
    if( siteIndex >= pList->designSiteCount ||
      rotamerIndex >= EnergyMatrixGetBlock(pMatrix,siteIndex,siteIndex)->RotamerCountSiteI ){
        result = ValueError;
        sprintf(errMsg,"in file %s line %d, rotamer #%d on site #%d does not exist",__FILE__,__LINE__,rotamerIndex,siteIndex);
        TraceError(errMsg,result);
        return result;
    }
    rotamerDeletedFlag[siteIndex][rotamerIndex] = TRUE;
  }

  for(int i=0;i<pList->designSiteCount;i++){
    int remainCount = 0;
    for(int j=0;j<pList->rotamerCount[i];j++){
      if(pList->remainFlag[i][j] == TRUE){
        if(rotamerDeletedFlag[i][remainCount]){
          pList->remainFlag[i][j] = FALSE;
        }
        remainCount++;      
      }
    }
    if(remainCount==0){
      result = ValueError;
      sprintf(errMsg,"in file %s line %d, all rotamers on design site %d have been deleted",__FILE__,__LINE__,i);
      TraceError(errMsg,result);
      return result;
    }
  }

  for(int i=0;i<pList->designSiteCount;i++){
    for(int j=i;j<pList->designSiteCount;j++){
      EnergyMatrixBlockUpdate(EnergyMatrixGetBlock(pMatrix,i,j),
        rotamerDeletedFlag[i],rotamerDeletedFlag[j]);
    }
  }

  for(int i=0;i<pList->designSiteCount;i++) free(rotamerDeletedFlag[i]);
  free(rotamerDeletedFlag);
  return Success;
}



int EnergyMatrixBlockUpdate(EnergyMatrixBlock* pEnergyBlock,BOOL* siteIRotamerDeletedFlag,BOOL* siteKRotamerDeletedFlag){
  int newRotamerCountSiteI = 0;
  int newRotamerCountSiteK = 0;
  double* newEnergyIK;
  double* pPositionInNewEnergyIK;

  for(int j=0;j<pEnergyBlock->RotamerCountSiteI;j++){
    if(siteIRotamerDeletedFlag[j]==FALSE) newRotamerCountSiteI++;
  }
  for(int s=0;s<pEnergyBlock->RotamerCountSiteK;s++){
    if(siteKRotamerDeletedFlag[s]==FALSE) newRotamerCountSiteK++;
  }
  newEnergyIK = (double*)malloc(sizeof(double)*newRotamerCountSiteI*newRotamerCountSiteK);
  pPositionInNewEnergyIK = &newEnergyIK[0];
  for(int j=0;j<pEnergyBlock->RotamerCountSiteI;j++){
    if(siteIRotamerDeletedFlag[j] == TRUE){
      continue;
    }
    for(int s=0;s<pEnergyBlock->RotamerCountSiteK;s++){
      if(siteKRotamerDeletedFlag[s] == TRUE) continue;
      *pPositionInNewEnergyIK = *EnergyMatrixBlockGet(pEnergyBlock,j,s);
      pPositionInNewEnergyIK++;
    }
  }

  pEnergyBlock->RotamerCountSiteI = newRotamerCountSiteI;
  pEnergyBlock->RotamerCountSiteK = newRotamerCountSiteK;
  free(pEnergyBlock->energyIK);
  pEnergyBlock->energyIK = newEnergyIK;

  return Success;
}


int EnergyMatrixUpdateByRotamerList(EnergyMatrix* pMatrix,RotamerList* pList){
  IntArray  dels;
  RotamerList rotamerList;
  IntArrayCreate(&dels, 0);
  RotamerListCreateFromEnergyMatrix(&rotamerList, pMatrix);
  for(int i=0; i<pList->designSiteCount; i++){
    for(int j=0; j<pList->rotamerCount[i]; j++){
      if(pList->remainFlag[i][j]==FALSE){
        IntArrayAppend(&dels,i);
        IntArrayAppend(&dels,j);
      }
    }
  }
  RotamerListAndEnergyMatrixDelete(&rotamerList, pMatrix, &dels);
  RotamerListDestroy(&rotamerList);
  IntArrayDestroy(&dels);
  return Success;
}


int StructureShowDesignSitesAfterRotamerDelete(Structure* pStructure,RotamerList* pList,FILE* pFile){
  int result = Success;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  double conformationSpace=0.0;
  int totalRotamerCount = 0;
  if(StructureGetDesignSiteCount(pStructure) != pList->designSiteCount){
    sprintf(errMsg,"in file %s line %d, number of design site in structure (%d) is not equal to that in rotamer list (%d)",
      __FILE__,__LINE__,StructureGetDesignSiteCount(pStructure),pList->designSiteCount);
    result = ValueError;
    TraceError(errMsg,result);
    return result;
  }
  for(int i = 0; i < StructureGetDesignSiteCount(pStructure); i++){
    DesignSite* pSiteI = StructureGetDesignSite(pStructure, i);
    RotamerSet* pSetI = DesignSiteGetRotamers(pSiteI);
    int remainRotamerCount = 0;

    if(RotamerSetGetCount(pSetI) != pList->rotamerCount[i]){
      sprintf(errMsg,"in file %s line %d, number of rotamers (%d) in site %d of structure is not equal to that in rotamer list (%d)",
        __FILE__,__LINE__,RotamerSetGetCount(pSetI),i,pList->rotamerCount[i]);
      result = ValueError;
      TraceError(errMsg,result);
      return result;
    }
    for(int j = 0; j < RotamerSetGetCount(pSetI); j++){
      Rotamer* pRotamer = RotamerSetGet(pSetI, j);
      if(pList->remainFlag[i][j] == TRUE){
        remainRotamerCount++;
      }
    }
    fprintf(stdout,"site %3d : %3s %s %4d, %6d rotamers : ",
      i,ResidueGetName(pSiteI->pResidue),ResidueGetChainName(pSiteI->pResidue),ResidueGetPosInChain(pSiteI->pResidue),remainRotamerCount);
    totalRotamerCount += remainRotamerCount;
    if(remainRotamerCount>0){
      conformationSpace += log((double)remainRotamerCount)/log(10.0);
    }
    switch(pSiteI->pResidue->designType){
      case Type_ResidueDesignType_Catalytic:
        fprintf(pFile,"catalytic\n"); break;
      case Type_ResidueDesignType_Fixed:
        fprintf(pFile,"fixed\n"); break;
      case Type_ResidueDesignType_Designable:
        fprintf(pFile,"mutated\n"); break;
      case Type_ResidueDesignType_SmallMol:
        fprintf(pFile,"smallmol\n"); break;
      case Type_ResidueDesignType_Repacked:
        fprintf(pFile,"rotameric\n"); break;
      default:
        break;
    }
  }
  fprintf(stdout, "total rotamer count: %d, conformation space: 1e^%.1f\n", totalRotamerCount, conformationSpace);

  return Success;
}


int RotamerDeleteBySelfEnergyCheck2(EnergyMatrix *pMatrix,Structure *pStructure,RotamerList* pList,IntArray *pDelete,EnergyMatrix *pFlag){
  IntArrayResize(pDelete, 0);
  for(int i=0; i<pMatrix->designSiteCount; i++){
    EnergyMatrixBlock *pBlock = EnergyMatrixGetBlock(pMatrix, i, i);
    for(int j=0; j<pBlock->RotamerCountSiteI; j++){
      int trueIndexJ;
      RotamerOriginalIndexGet(pList, i, j, &trueIndexJ);
      Rotamer *pRotamerIJ = RotamerSetGet(DesignSiteGetRotamers(StructureGetDesignSite(pStructure, i)), trueIndexJ);
      for(int u=0; u<pBlock->RotamerCountSiteI; u++){
        Rotamer *pRotamerIU;
        int trueIndexU;
        RotamerOriginalIndexGet(pList, i, u, &trueIndexU);
        pRotamerIU = RotamerSetGet(DesignSiteGetRotamers(StructureGetDesignSite(pStructure, i)), trueIndexU);
        if(*EnergyMatrixGet(pFlag, i, i, u, u) < 0.0) continue;
        //for different rotamer types
        if(*EnergyMatrixGet(pMatrix, i, i, j, j) > *EnergyMatrixGet(pMatrix, i, i, u, u) + SELF_ENERGY_DIFFERENT_THRESHOLD){
          IntArrayAppend(pDelete, i);
          IntArrayAppend(pDelete, j);
          *EnergyMatrixGet(pFlag, i, i, j, j) = -1.0;
          break;
        }
        //for same rotamer types
        else if(RotamerAndRotamerInSameType(pRotamerIJ,pRotamerIU) &&
          *EnergyMatrixGet(pMatrix, i, i, j, j) > *EnergyMatrixGet(pMatrix, i, i, u, u) + SELF_ENERGY_SAME_THRESHOLD){
            IntArrayAppend(pDelete, i);
            IntArrayAppend(pDelete, j);
            *EnergyMatrixGet(pFlag, i, i, j, j) = -1.0;
            break;
        }
      }
    }
  }

  return Success;
}


int DesignSiteGetRemainRotamerCount(int designSiteI, RotamerList* pList){
  int result=Success;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  if(designSiteI < 0 || designSiteI >= pList->designSiteCount){
    sprintf(errMsg,"in file %s line %d, invalid index\n",__FILE__,__LINE__);
    result = ValueError;
    TraceError(errMsg,ValueError);
    return 0;
  }
  int remainRotamerCount = 0;
  for(int j=0; j<pList->rotamerCount[designSiteI]; j++){
    if(pList->remainFlag[designSiteI][j]) remainRotamerCount++;
  }
  return remainRotamerCount;
}


int RotamerListUpdateByDeleteArray(RotamerList* pList, IntArray* pDeleteArray){
  for(int i=0; i<IntArrayGetLength(pDeleteArray); i+=2){
    int designSiteI = IntArrayGet(pDeleteArray, i);
    int rotamerIJ = IntArrayGet(pDeleteArray, i+1);
    if(designSiteI < 0 || designSiteI >= pList->designSiteCount|| rotamerIJ < 0 || rotamerIJ >= pList->rotamerCount[designSiteI]){
        char errMsg[MAX_LENGTH_ERR_MSG+1];
        sprintf(errMsg,"in file %s line %d, invalid parameter",__FILE__,__LINE__);
        TraceError(errMsg,ValueError);
        return ValueError;
    }
    pList->remainFlag[designSiteI][rotamerIJ] = FALSE;
  }

  return Success;
}


int RotamerOriginalIndexGet(RotamerList* pList, int designSiteI, int rotamerIJ, int *trueIndexIJ){
  int j, flag = 0;
  for(j=0; j<pList->rotamerCount[designSiteI]; j++){
    if(pList->remainFlag[designSiteI][j] == FALSE) continue;
    if(flag == rotamerIJ) break;
    flag++;
  }
  if(j == pList->rotamerCount[designSiteI]) return DataNotExistError;
  else *trueIndexIJ = j;
  return Success;
}


int RotamerReducedIndexGet(RotamerList* pList, int designSiteI, int trueIndexIJ, int *reducedIndex){
  int j, flag = 0;
  for(j=0; j<pList->rotamerCount[designSiteI]; j++){
    if(pList->remainFlag[designSiteI][j] == FALSE) continue;
    if(j == trueIndexIJ) break;
    flag++;
  }
  if(j == pList->rotamerCount[designSiteI]) return DataNotExistError;
  else *reducedIndex = flag;
  return Success;
}


int EnergyMatrixReadNew(EnergyMatrix* pThis, RotamerList* pList, char* energyMatrixFile){
  int *remainRotamerCount = (int*)malloc(sizeof(int)*pList->designSiteCount);
  for(int i=0; i<pList->designSiteCount; i++){
    remainRotamerCount[i] = 0;
    for(int j = 0; j < pList->rotamerCount[i]; j++){
      if(pList->remainFlag[i][j] == TRUE) remainRotamerCount[i]++;
    }
  }

  pThis->designSiteCount = pList->designSiteCount;
  pThis->blocks = (EnergyMatrixBlock*)malloc(sizeof(EnergyMatrixBlock)*pThis->designSiteCount*pThis->designSiteCount);
  for(int i=0; i<pThis->designSiteCount*pThis->designSiteCount; i++) EnergyMatrixBlockCreate(&pThis->blocks[i]);
  for(int i=0; i<pThis->designSiteCount; i++){
    for(int k=i; k<pThis->designSiteCount; k++){
      EnergyMatrixBlock* pBlockIK = EnergyMatrixGetBlock(pThis, i, k);
      pBlockIK->DesignSiteI = i;
      pBlockIK->DesignSiteK = k;
      pBlockIK->RotamerCountSiteI = remainRotamerCount[i];
      pBlockIK->RotamerCountSiteK = remainRotamerCount[k];
      pBlockIK->energyIK = (double*)malloc(sizeof(double)*pBlockIK->RotamerCountSiteI*pBlockIK->RotamerCountSiteK);
      for(int j = 0; j < pBlockIK->RotamerCountSiteI*pBlockIK->RotamerCountSiteK; j++) pBlockIK->energyIK[j] = 1000.0;
    }
  }

  FILE* pFile = fopen(energyMatrixFile, "r");
  if(!pFile){
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    sprintf(errMsg,"in file %s line %d, cannot open file %s",__FILE__,__LINE__,energyMatrixFile);
    TraceError(errMsg,IOError);
    return IOError;
  }

  long int lineCounter = 0;
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(fgets(buffer, MAX_LENGTH_ONE_LINE_IN_FILE, pFile)){
    double energy;
    int siteI, siteK, rotamerIJ, rotamerKS;
    sscanf(buffer, "%lf %d %d %d %d\n", &energy, &siteI, &rotamerIJ, &siteK, &rotamerKS);
    *EnergyMatrixGet(pThis, siteI, siteK, rotamerIJ, rotamerKS) = energy;
  }
  fclose(pFile);
  free(remainRotamerCount);

  return Success;
}


int DeleteArrayGenerateFromTwoRotamerList(IntArray* pDeleteRotamerArray, RotamerList* pOld, RotamerList* pNew){
  IntArrayResize(pDeleteRotamerArray,0);
  for(int i=0; i<pOld->designSiteCount; i++){
    int index = 0;
    for(int j=0; j<pOld->rotamerCount[i]; j++){
      if(pOld->remainFlag[i][j] == TRUE && pNew->remainFlag[i][j] == TRUE) index++;
      if(pOld->remainFlag[i][j] == TRUE && pNew->remainFlag[i][j] == FALSE){
        IntArrayAppend(pDeleteRotamerArray, i);
        IntArrayAppend(pDeleteRotamerArray, index);
        index++;
      }
    }
  }

  return Success;
}


int RotamerListCreateFromFile(RotamerList *pThis, char* filepath){
  int result;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  FILE* pFile = fopen(filepath,"r");
  if(!pFile){
    sprintf(errMsg,"in file %s line %d, cannot open file %s",__FILE__,__LINE__,filepath);
    result = IOError;
    TraceError(errMsg,result);
    return result;
  }

  int designSiteCount = 0;
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    int i,j;
    char flag[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    if(sscanf(line,"%d %d %s",&i,&j,flag) != 3){
      sprintf(errMsg,"in file %s line %d, wrong format in file %s line\n%s",__FILE__,__LINE__,filepath,line);
      result = IOError;
      TraceError(errMsg,result);
      return result;
    }
    if(i>=designSiteCount) designSiteCount++;
  }

  IntArray remainRotCountArray;
  IntArrayCreate(&remainRotCountArray, designSiteCount);
  for(int i=0; i<IntArrayGetLength(&remainRotCountArray); i++){
    IntArraySet(&remainRotCountArray,i,0);
  }

  fseek(pFile, 0, SEEK_SET);
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    int i,j;
    char flag[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sscanf(line,"%d %d %s",&i,&j,flag);
    if(j>=IntArrayGet(&remainRotCountArray, i)){
      IntArraySet(&remainRotCountArray, i, j+1);
    }
  }
  fseek(pFile, 0, SEEK_SET);

  pThis->designSiteCount = designSiteCount;
  pThis->rotamerCount = (int*)malloc(sizeof(int)*pThis->designSiteCount);
  pThis->remainFlag = (BOOL**)malloc(sizeof(BOOL*)*pThis->designSiteCount);
  for(int i=0; i<pThis->designSiteCount; i++){
    pThis->rotamerCount[i] = IntArrayGet(&remainRotCountArray, i);
    pThis->remainFlag[i] = (BOOL*)malloc(sizeof(BOOL)*pThis->rotamerCount[i]);
    for(int j=0;j<pThis->rotamerCount[i]; j++){
      pThis->remainFlag[i][j] = TRUE;
    }
  }
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    int i,j;
    char flag[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sscanf(line,"%d %d %s",&i,&j,flag);
    if(strcmp(flag, "TRUE") == 0) pThis->remainFlag[i][j] = TRUE;
    else if(strcmp(flag, "FALSE") == 0) pThis->remainFlag[i][j] = FALSE;
  }
  fclose(pFile);
  IntArrayDestroy(&remainRotCountArray);

  return Success;
}

