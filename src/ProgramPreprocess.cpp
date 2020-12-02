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

#include "ProgramPreprocess.h"
#include "Residue.h"
#include "SmallMolEEF1.h"
#include "ErrorTracker.h"
#include "EvoAminoName.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

extern char PROGRAM_PATH[MAX_LENGTH_ONE_LINE_IN_FILE+1];
extern char PROGRAM_NAME[MAX_LENGTH_FILE_NAME+1];
extern char PROGRAM_VERSION[MAX_LENGTH_FILE_NAME+1];
extern char PDBID[MAX_LENGTH_FILE_NAME+1];
extern char TGT_PRF[MAX_LENGTH_FILE_NAME+1];
extern char TGT_MSA[MAX_LENGTH_FILE_NAME+1];
extern char DES_CHAINS[MAX_LENGTH_ONE_LINE_IN_FILE+1];
extern float** PROT_PROFILE;

using namespace CharSeq;



int ExtractPathAndName(char* fullpath,char* path,char* name){
  BOOL slash=FALSE;
  for(int i=strlen(fullpath); i>=0; i--){
    if(fullpath[i]=='/' || fullpath[i]=='\\'){
      strcpy(name,fullpath+i+1);
      strncpy(path, fullpath, i);
      path[i]='\0';
      slash=TRUE;
      break;
    }
  }
  if(slash==FALSE){
    strcpy(name,fullpath);
    strcpy(path,".");
  }
  return Success;
}


int GetPDBID(char* pdbfile, char* pdbid){
  int result=Success;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
  int len=(int)strlen(pdbfile);
  if(len<=4 || tolower(pdbfile[len-1])!='b' || tolower(pdbfile[len-2])!='d' || tolower(pdbfile[len-3])!='p' || tolower(pdbfile[len-4])!='.'){
    strcpy(pdbid,"yourpdb");
    sprintf(errMsg,"in file %s line %d, failed to parse pdb ID from %s, designate the pdb ID as 'yourpdb'",__FILE__,__LINE__,pdbfile);
    result = Warning;
    TraceError(errMsg,Warning);
    return result;
  }
  else{
    strcpy(pdbid,pdbfile);
    pdbid[len-4]='\0';
  }

  return Success;
}


int GenerateProfile(char* pdb){
  char cmd[MAX_LENGTH_ONE_LINE_IN_FILE+1]="";
  sprintf(cmd,"perl %s/extbin/GenerateProfile.pl %s",PROGRAM_PATH,pdb);
  system(cmd);
  return Success;
}


int ReadProfile(Chain* pChain){
  //read profile
  FILE *pIn=fopen(TGT_PRF,"r");
  int len=ChainGetResidueCount(pChain);
  PROT_PROFILE = new float*[len];
  for(int i=0; i<len; i++){
    PROT_PROFILE[i] = new float[AMINOS];
    for(int j=0;j<AMINOS;j++){
      PROT_PROFILE[i][j]=0;
    }
  }
  AminoName map[AMINOS];
  char ami[2];
  for(int i=0; i<AMINOS; ++i) {
    fscanf(pIn, "%s", ami);
    map[i] = charToAmino(ami[0]);
  }
  fscanf(pIn, "%*s");
  for(int i=0; i<len; ++i) {
    fscanf(pIn, "%*s");
    for(int j=0; j<AMINOS; ++j) {
      fscanf(pIn, "%f", &PROT_PROFILE[i][map[j]]);
    }
    fscanf(pIn, "%*s");
  }
  fclose(pIn);
  return Success;
}


int GenerateSingleSequenceFeatures(char* pdb){
  char cmd[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  printf("secondary structure prediction of target will be done by DSSP\n");
  printf("running DSSP ... \n");
  sprintf(cmd,"%s/extbin/dsspcmbi %s > dssp.txt 2>/dev/null",PROGRAM_PATH,pdb);
  system(cmd);
  printf("done\n");
  printf("extract seondary structure information ... \n");
  sprintf(cmd,"%s/extbin/dssp-ver2hor dssp.txt > dssph.txt",PROGRAM_PATH);
  system(cmd);
  printf("done\n");
  printf("convert to I-TASSER like secondary structure format ... \n");
  sprintf(cmd,"%s/extbin/convSS.pl dssph.txt > ss.txt",PROGRAM_PATH);
  system(cmd);
  printf("done\n");
  printf("generated: secondary structure, surface's solvent-accessibility, and phi-psi files.\n");
  return Success;
}


int StructureDeployEvolutionInfo(Structure* pStructure){
  Chain* pChain=StructureFindChainByName(pStructure,DES_CHAINS);
  char TGT_CHAIN[MAX_LENGTH_FILE_NAME+1];
  sprintf(TGT_CHAIN,"%s_TGTCHN.pdb",PDBID);
  FILE* pFile=fopen(TGT_CHAIN,"w");
  if(pFile!=NULL){
    for(int i=0;i<ChainGetResidueCount(pChain);i++){
      Residue* pResi=ChainGetResidue(pChain,i);
      AtomArrayShowInPDBFormat(&pResi->atoms,"ATOM",ResidueGetName(pResi),ChainGetName(pChain),1,ResidueGetPosInChain(pResi),pFile);
    }
    fclose(pFile);
    printf("create an temp pdb %s for chain %s used for TM-align search\n",TGT_CHAIN,DES_CHAINS);
  }
  else{
    char errMsg[MAX_LENGTH_ERR_MSG+1];
    sprintf(errMsg,"in file %s line %d, cannot create a temp pdb %s for chain %s",__FILE__,__LINE__,TGT_CHAIN,DES_CHAINS);
    TraceError(errMsg,IOError);
    exit(IOError);
  }

  FileReader frPrf;
  if((!FAILED(FileReaderCreate(&frPrf,TGT_PRF)) && FileReaderGetLineCount(&frPrf)!=ChainGetResidueCount(pChain)) || FAILED(FileReaderCreate(&frPrf,TGT_PRF))){
    FileReader frMSA;
    if((!FAILED(FileReaderCreate(&frMSA,TGT_MSA)) && FileReaderGetLineCount(&frMSA)!=ChainGetResidueCount(pChain)) || FAILED(FileReaderCreate(&frMSA,TGT_MSA))){
      GenerateProfile(TGT_CHAIN);
    }
    else{
      char cmd[MAX_LENGTH_ONE_LINE_IN_FILE+1];
      printf("make profile from msa file %s\n", TGT_MSA);
      sprintf(cmd,"%s/extbin/profile/mkprf %s > %s 2>/dev/null",PROGRAM_PATH,TGT_MSA,TGT_PRF);
      system(cmd);
    }
    ReadProfile(pChain);
  }
  FileReaderDestroy(&frPrf);
  GenerateSingleSequenceFeatures(TGT_CHAIN);
  return Success;
}

