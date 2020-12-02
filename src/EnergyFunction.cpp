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

#include "EnergyFunction.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>


extern double CUT_EXCL_LOW_ROT_PROB;

double WEIGHTS[MAX_ENERGY_TERM];

#define DEAL_WITH_ENERGY_WEIGHTING
int EnergyTermInitialize(double *energyTerms){
  memset(energyTerms,0,sizeof(double)*MAX_ENERGY_TERM);
  return Success;
}

int EnergyTermWeighting(double *energyTerms){
  for(int i=1; i<MAX_ENERGY_TERM; i++){
    energyTerms[i] *= WEIGHTS[i];
    energyTerms[0]+=energyTerms[i];
  }
  return Success;
}


int EnergyTermShowMonomer(double *energyTerms){
  printf("reference_ALA         =            %12.6f\n", energyTerms[ 1]);
  printf("reference_CYS         =            %12.6f\n", energyTerms[ 2]);
  printf("reference_ASP         =            %12.6f\n", energyTerms[ 3]);
  printf("reference_GLU         =            %12.6f\n", energyTerms[ 4]);
  printf("reference_PHE         =            %12.6f\n", energyTerms[ 5]);
  printf("reference_GLY         =            %12.6f\n", energyTerms[ 6]);
  printf("reference_HIS         =            %12.6f\n", energyTerms[ 7]);
  printf("reference_ILE         =            %12.6f\n", energyTerms[ 8]);
  printf("reference_LYS         =            %12.6f\n", energyTerms[ 9]);
  printf("reference_LEU         =            %12.6f\n", energyTerms[10]);
  printf("reference_MET         =            %12.6f\n", energyTerms[11]);
  printf("reference_ASN         =            %12.6f\n", energyTerms[12]);
  printf("reference_PRO         =            %12.6f\n", energyTerms[13]);
  printf("reference_GLN         =            %12.6f\n", energyTerms[14]);
  printf("reference_ARG         =            %12.6f\n", energyTerms[15]);
  printf("reference_SER         =            %12.6f\n", energyTerms[16]);
  printf("reference_THR         =            %12.6f\n", energyTerms[17]);
  printf("reference_VAL         =            %12.6f\n", energyTerms[18]);
  printf("reference_TRP         =            %12.6f\n", energyTerms[19]);
  printf("reference_TYR         =            %12.6f\n", energyTerms[20]);
  printf("intraR_vdwatt         =            %12.6f\n", energyTerms[21]);
  printf("intraR_vdwrep         =            %12.6f\n", energyTerms[22]);
  printf("intraR_electr         =            %12.6f\n", energyTerms[23]);
  printf("intraR_deslvP         =            %12.6f\n", energyTerms[24]);
  printf("intraR_deslvH         =            %12.6f\n", energyTerms[25]);
  printf("intraR_hbscbb_dis     =            %12.6f\n", energyTerms[26]);
  printf("intraR_hbscbb_the     =            %12.6f\n", energyTerms[27]);
  printf("intraR_hbscbb_phi     =            %12.6f\n", energyTerms[28]);
  printf("interS_vdwatt         =            %12.6f\n", energyTerms[31]);
  printf("interS_vdwrep         =            %12.6f\n", energyTerms[32]);
  printf("interS_electr         =            %12.6f\n", energyTerms[33]);
  printf("interS_deslvP         =            %12.6f\n", energyTerms[34]);
  printf("interS_deslvH         =            %12.6f\n", energyTerms[35]);
  printf("interS_hbbbbb_dis     =            %12.6f\n", energyTerms[41]);
  printf("interS_hbbbbb_the     =            %12.6f\n", energyTerms[42]);
  printf("interS_hbbbbb_phi     =            %12.6f\n", energyTerms[43]);
  printf("interS_hbscbb_dis     =            %12.6f\n", energyTerms[44]);
  printf("interS_hbscbb_the     =            %12.6f\n", energyTerms[45]);
  printf("interS_hbscbb_phi     =            %12.6f\n", energyTerms[46]);
  printf("interS_hbscsc_dis     =            %12.6f\n", energyTerms[47]);
  printf("interS_hbscsc_the     =            %12.6f\n", energyTerms[48]);
  printf("interS_hbscsc_phi     =            %12.6f\n", energyTerms[49]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %12.6f\n\n", energyTerms[0]);

  return Success;
}


int EnergyTermShowComplex(double* energyTerms){
  printf("reference_ALA         =            %12.6f\n", energyTerms[ 1]);
  printf("reference_CYS         =            %12.6f\n", energyTerms[ 2]);
  printf("reference_ASP         =            %12.6f\n", energyTerms[ 3]);
  printf("reference_GLU         =            %12.6f\n", energyTerms[ 4]);
  printf("reference_PHE         =            %12.6f\n", energyTerms[ 5]);
  printf("reference_GLY         =            %12.6f\n", energyTerms[ 6]);
  printf("reference_HIS         =            %12.6f\n", energyTerms[ 7]);
  printf("reference_ILE         =            %12.6f\n", energyTerms[ 8]);
  printf("reference_LYS         =            %12.6f\n", energyTerms[ 9]);
  printf("reference_LEU         =            %12.6f\n", energyTerms[10]);
  printf("reference_MET         =            %12.6f\n", energyTerms[11]);
  printf("reference_ASN         =            %12.6f\n", energyTerms[12]);
  printf("reference_PRO         =            %12.6f\n", energyTerms[13]);
  printf("reference_GLN         =            %12.6f\n", energyTerms[14]);
  printf("reference_ARG         =            %12.6f\n", energyTerms[15]);
  printf("reference_SER         =            %12.6f\n", energyTerms[16]);
  printf("reference_THR         =            %12.6f\n", energyTerms[17]);
  printf("reference_VAL         =            %12.6f\n", energyTerms[18]);
  printf("reference_TRP         =            %12.6f\n", energyTerms[19]);
  printf("reference_TYR         =            %12.6f\n", energyTerms[20]);

  printf("intraR_vdwatt         =            %12.6f\n", energyTerms[21]);
  printf("intraR_vdwrep         =            %12.6f\n", energyTerms[22]);
  printf("intraR_electr         =            %12.6f\n", energyTerms[23]);
  printf("intraR_deslvP         =            %12.6f\n", energyTerms[24]);
  printf("intraR_deslvH         =            %12.6f\n", energyTerms[25]);
  printf("intraR_hbscbb_dis     =            %12.6f\n", energyTerms[26]);
  printf("intraR_hbscbb_the     =            %12.6f\n", energyTerms[27]);
  printf("intraR_hbscbb_phi     =            %12.6f\n", energyTerms[28]);
  printf("aapropensity          =            %12.6f\n", energyTerms[91]);
  printf("ramachandran          =            %12.6f\n", energyTerms[92]);
  printf("dunbrack              =            %12.6f\n", energyTerms[93]);

  printf("interS_vdwatt         =            %12.6f\n", energyTerms[31]);
  printf("interS_vdwrep         =            %12.6f\n", energyTerms[32]);
  printf("interS_electr         =            %12.6f\n", energyTerms[33]);
  printf("interS_deslvP         =            %12.6f\n", energyTerms[34]);
  printf("interS_deslvH         =            %12.6f\n", energyTerms[35]);
  printf("interS_ssbond         =            %12.6f\n", energyTerms[36]);
  printf("interS_hbbbbb_dis     =            %12.6f\n", energyTerms[41]);
  printf("interS_hbbbbb_the     =            %12.6f\n", energyTerms[42]);
  printf("interS_hbbbbb_phi     =            %12.6f\n", energyTerms[43]);
  printf("interS_hbscbb_dis     =            %12.6f\n", energyTerms[44]);
  printf("interS_hbscbb_the     =            %12.6f\n", energyTerms[45]);
  printf("interS_hbscbb_phi     =            %12.6f\n", energyTerms[46]);
  printf("interS_hbscsc_dis     =            %12.6f\n", energyTerms[47]);
  printf("interS_hbscsc_the     =            %12.6f\n", energyTerms[48]);
  printf("interS_hbscsc_phi     =            %12.6f\n", energyTerms[49]);

  printf("interD_vdwatt         =            %12.6f\n", energyTerms[51]);
  printf("interD_vdwrep         =            %12.6f\n", energyTerms[52]);
  printf("interD_electr         =            %12.6f\n", energyTerms[53]);
  printf("interD_deslvP         =            %12.6f\n", energyTerms[54]);
  printf("interD_deslvH         =            %12.6f\n", energyTerms[55]);
  printf("interD_ssbond         =            %12.6f\n", energyTerms[56]);
  printf("interD_hbbbbb_dis     =            %12.6f\n", energyTerms[61]);
  printf("interD_hbbbbb_the     =            %12.6f\n", energyTerms[62]);
  printf("interD_hbbbbb_phi     =            %12.6f\n", energyTerms[63]);
  printf("interD_hbscbb_dis     =            %12.6f\n", energyTerms[64]);
  printf("interD_hbscbb_the     =            %12.6f\n", energyTerms[65]);
  printf("interD_hbscbb_phi     =            %12.6f\n", energyTerms[66]);
  printf("interD_hbscsc_dis     =            %12.6f\n", energyTerms[67]);
  printf("interD_hbscsc_the     =            %12.6f\n", energyTerms[68]);
  printf("interD_hbscsc_phi     =            %12.6f\n", energyTerms[69]);

  printf("prolig_vdwatt         =            %12.6f\n", energyTerms[71]);
  printf("prolig_vdwrep         =            %12.6f\n", energyTerms[72]);
  printf("prolig_electr         =            %12.6f\n", energyTerms[73]);
  printf("prolig_deslvP         =            %12.6f\n", energyTerms[74]);
  printf("prolig_deslvH         =            %12.6f\n", energyTerms[75]);
  printf("prolig_hbscbb_dis     =            %12.6f\n", energyTerms[81]);
  printf("prolig_hbscbb_the     =            %12.6f\n", energyTerms[82]);
  printf("prolig_hbscbb_phi     =            %12.6f\n", energyTerms[83]);
  printf("prolig_hbscsc_dis     =            %12.6f\n", energyTerms[84]);
  printf("prolig_hbscsc_the     =            %12.6f\n", energyTerms[85]);
  printf("prolig_hbscsc_phi     =            %12.6f\n", energyTerms[86]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %12.6f\n\n", energyTerms[0]);
  return Success;
}


int EnergyWeightRead(char* weightfile){
  for(int i=0;i<MAX_ENERGY_TERM;i++){
    WEIGHTS[i]=1.0;
  }
  FILE* pf=fopen(weightfile,"r");
  if(pf!=NULL){
    char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pf)){
      char term[MAX_LENGTH_ONE_LINE_IN_FILE+1];
      double val=0.0;
      sscanf(line,"%s %lf",term,&val);
      if(!strcmp(term,"reference_ALA")) WEIGHTS[ 1]=val;
      else if(!strcmp(term,"reference_CYS")) WEIGHTS[ 2]=val;
      else if(!strcmp(term,"reference_ASP")) WEIGHTS[ 3]=val;
      else if(!strcmp(term,"reference_GLU")) WEIGHTS[ 4]=val;
      else if(!strcmp(term,"reference_PHE")) WEIGHTS[ 5]=val;
      else if(!strcmp(term,"reference_GLY")) WEIGHTS[ 6]=val;
      else if(!strcmp(term,"reference_HIS")) WEIGHTS[ 7]=val;
      else if(!strcmp(term,"reference_ILE")) WEIGHTS[ 8]=val;
      else if(!strcmp(term,"reference_LYS")) WEIGHTS[ 9]=val;
      else if(!strcmp(term,"reference_LEU")) WEIGHTS[10]=val;
      else if(!strcmp(term,"reference_MET")) WEIGHTS[11]=val;
      else if(!strcmp(term,"reference_ASN")) WEIGHTS[12]=val;
      else if(!strcmp(term,"reference_PRO")) WEIGHTS[13]=val;
      else if(!strcmp(term,"reference_GLN")) WEIGHTS[14]=val;
      else if(!strcmp(term,"reference_ARG")) WEIGHTS[15]=val;
      else if(!strcmp(term,"reference_SER")) WEIGHTS[16]=val;
      else if(!strcmp(term,"reference_THR")) WEIGHTS[17]=val;
      else if(!strcmp(term,"reference_VAL")) WEIGHTS[18]=val;
      else if(!strcmp(term,"reference_TRP")) WEIGHTS[19]=val;
      else if(!strcmp(term,"reference_TYR")) WEIGHTS[20]=val;

      else if(!strcmp(term,"intraR_vdwatt")) WEIGHTS[21]=val;
      else if(!strcmp(term,"intraR_vdwrep")) WEIGHTS[22]=val;
      else if(!strcmp(term,"intraR_electr")) WEIGHTS[23]=val;
      else if(!strcmp(term,"intraR_deslvP")) WEIGHTS[24]=val;
      else if(!strcmp(term,"intraR_deslvH")) WEIGHTS[25]=val;
      else if(!strcmp(term,"intraR_hbscbb_dis")) WEIGHTS[26]=val;
      else if(!strcmp(term,"intraR_hbscbb_the")) WEIGHTS[27]=val;
      else if(!strcmp(term,"intraR_hbscbb_phi")) WEIGHTS[28]=val;
      else if(!strcmp(term,"aapropensity")) WEIGHTS[91]=val;
      else if(!strcmp(term,"ramachandran")) WEIGHTS[92]=val;
      else if(!strcmp(term,"dunbrack"))     WEIGHTS[93]=val;

      else if(!strcmp(term,"interS_vdwatt")) WEIGHTS[31]=val;
      else if(!strcmp(term,"interS_vdwrep")) WEIGHTS[32]=val;
      else if(!strcmp(term,"interS_electr")) WEIGHTS[33]=val;
      else if(!strcmp(term,"interS_deslvP")) WEIGHTS[34]=val;
      else if(!strcmp(term,"interS_deslvH")) WEIGHTS[35]=val;
      else if(!strcmp(term,"interS_ssbond")) WEIGHTS[36]=val;
      else if(!strcmp(term,"interS_hbbbbb_dis")) WEIGHTS[41]=val;
      else if(!strcmp(term,"interS_hbbbbb_the")) WEIGHTS[42]=val;
      else if(!strcmp(term,"interS_hbbbbb_phi")) WEIGHTS[43]=val;
      else if(!strcmp(term,"interS_hbscbb_dis")) WEIGHTS[44]=val;
      else if(!strcmp(term,"interS_hbscbb_the")) WEIGHTS[45]=val;
      else if(!strcmp(term,"interS_hbscbb_phi")) WEIGHTS[46]=val;
      else if(!strcmp(term,"interS_hbscsc_dis")) WEIGHTS[47]=val;
      else if(!strcmp(term,"interS_hbscsc_the")) WEIGHTS[48]=val;
      else if(!strcmp(term,"interS_hbscsc_phi")) WEIGHTS[49]=val;

      else if(!strcmp(term,"interD_vdwatt")) WEIGHTS[51]=val;
      else if(!strcmp(term,"interD_vdwrep")) WEIGHTS[52]=val;
      else if(!strcmp(term,"interD_electr")) WEIGHTS[53]=val;
      else if(!strcmp(term,"interD_deslvP")) WEIGHTS[54]=val;
      else if(!strcmp(term,"interD_deslvH")) WEIGHTS[55]=val;
      else if(!strcmp(term,"interD_ssbond")) WEIGHTS[56]=val;
      else if(!strcmp(term,"interD_hbbbbb_dis")) WEIGHTS[61]=val;
      else if(!strcmp(term,"interD_hbbbbb_the")) WEIGHTS[62]=val;
      else if(!strcmp(term,"interD_hbbbbb_phi")) WEIGHTS[63]=val;
      else if(!strcmp(term,"interD_hbscbb_dis")) WEIGHTS[64]=val;
      else if(!strcmp(term,"interD_hbscbb_the")) WEIGHTS[65]=val;
      else if(!strcmp(term,"interD_hbscbb_phi")) WEIGHTS[66]=val;
      else if(!strcmp(term,"interD_hbscsc_dis")) WEIGHTS[67]=val;
      else if(!strcmp(term,"interD_hbscsc_the")) WEIGHTS[68]=val;
      else if(!strcmp(term,"interD_hbscsc_phi")) WEIGHTS[69]=val;

      else if(!strcmp(term,"ligand_vdwatt")) WEIGHTS[71]=val;
      else if(!strcmp(term,"ligand_vdwrep")) WEIGHTS[72]=val;
      else if(!strcmp(term,"ligand_electr")) WEIGHTS[73]=val;
      else if(!strcmp(term,"ligand_deslvP")) WEIGHTS[74]=val;
      else if(!strcmp(term,"ligand_deslvH")) WEIGHTS[75]=val;
      else if(!strcmp(term,"ligand_hbscbb_dis")) WEIGHTS[84]=val;
      else if(!strcmp(term,"ligand_hbscbb_the")) WEIGHTS[85]=val;
      else if(!strcmp(term,"ligand_hbscbb_phi")) WEIGHTS[86]=val;
      else if(!strcmp(term,"ligand_hbscsc_dis")) WEIGHTS[87]=val;
      else if(!strcmp(term,"ligand_hbscsc_the")) WEIGHTS[88]=val;
      else if(!strcmp(term,"ligand_hbscsc_phi")) WEIGHTS[89]=val;
    }
    fclose(pf);
  }
  else{
    char errMsg[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    sprintf(errMsg,"in file %s line %d, cannot open weight file %s",__FILE__,__LINE__,weightfile);
    TraceError(errMsg,IOError);
    return IOError;
  }

  return Success;
}

int EnergyWeightWrite(char* weightfile){
  FILE* fo=fopen(weightfile,"w");
  fprintf(fo,"reference_ALA      %8.3f\n",WEIGHTS[ 1]);
  fprintf(fo,"reference_CYS      %8.3f\n",WEIGHTS[ 2]);
  fprintf(fo,"reference_ASP      %8.3f\n",WEIGHTS[ 3]);
  fprintf(fo,"reference_GLU      %8.3f\n",WEIGHTS[ 4]);
  fprintf(fo,"reference_PHE      %8.3f\n",WEIGHTS[ 5]);
  fprintf(fo,"reference_GLY      %8.3f\n",WEIGHTS[ 6]);
  fprintf(fo,"reference_HIS      %8.3f\n",WEIGHTS[ 7]);
  fprintf(fo,"reference_ILE      %8.3f\n",WEIGHTS[ 8]);
  fprintf(fo,"reference_LYS      %8.3f\n",WEIGHTS[ 9]);
  fprintf(fo,"reference_LEU      %8.3f\n",WEIGHTS[10]);
  fprintf(fo,"reference_MET      %8.3f\n",WEIGHTS[11]);
  fprintf(fo,"reference_ASN      %8.3f\n",WEIGHTS[12]);
  fprintf(fo,"reference_PRO      %8.3f\n",WEIGHTS[13]);
  fprintf(fo,"reference_GLN      %8.3f\n",WEIGHTS[14]);
  fprintf(fo,"reference_ARG      %8.3f\n",WEIGHTS[15]);
  fprintf(fo,"reference_SER      %8.3f\n",WEIGHTS[16]);
  fprintf(fo,"reference_THR      %8.3f\n",WEIGHTS[17]);
  fprintf(fo,"reference_VAL      %8.3f\n",WEIGHTS[18]);
  fprintf(fo,"reference_TRP      %8.3f\n",WEIGHTS[19]);
  fprintf(fo,"reference_TYR      %8.3f\n",WEIGHTS[20]);

  fprintf(fo,"intraR_vdwatt      %8.3f\n",WEIGHTS[21]);
  fprintf(fo,"intraR_vdwrep      %8.3f\n",WEIGHTS[22]);
  fprintf(fo,"intraR_electr      %8.3f\n",WEIGHTS[23]);
  fprintf(fo,"intraR_deslvP      %8.3f\n",WEIGHTS[24]);
  fprintf(fo,"intraR_deslvH      %8.3f\n",WEIGHTS[25]);
  fprintf(fo,"intraR_hbscbb_dis  %8.3f\n",WEIGHTS[26]);
  fprintf(fo,"intraR_hbscbb_the  %8.3f\n",WEIGHTS[27]);
  fprintf(fo,"intraR_hbscbb_phi  %8.3f\n",WEIGHTS[28]);
  fprintf(fo,"aapropensity       %8.3f\n",WEIGHTS[91]);
  fprintf(fo,"ramachandran       %8.3f\n",WEIGHTS[92]);
  fprintf(fo,"dunbrack           %8.3f\n",WEIGHTS[93]);

  fprintf(fo,"interS_vdwatt      %8.3f\n",WEIGHTS[31]);
  fprintf(fo,"interS_vdwrep      %8.3f\n",WEIGHTS[32]);
  fprintf(fo,"interS_electr      %8.3f\n",WEIGHTS[33]);
  fprintf(fo,"interS_deslvP      %8.3f\n",WEIGHTS[34]);
  fprintf(fo,"interS_deslvH      %8.3f\n",WEIGHTS[35]);
  fprintf(fo,"interS_ssbond      %8.3f\n",WEIGHTS[36]);
  fprintf(fo,"interS_hbbbbb_dis  %8.3f\n",WEIGHTS[41]);
  fprintf(fo,"interS_hbbbbb_the  %8.3f\n",WEIGHTS[42]);
  fprintf(fo,"interS_hbbbbb_phi  %8.3f\n",WEIGHTS[43]);
  fprintf(fo,"interS_hbscbb_dis  %8.3f\n",WEIGHTS[44]);
  fprintf(fo,"interS_hbscbb_the  %8.3f\n",WEIGHTS[45]);
  fprintf(fo,"interS_hbscbb_phi  %8.3f\n",WEIGHTS[46]);
  fprintf(fo,"interS_hbscsc_dis  %8.3f\n",WEIGHTS[47]);
  fprintf(fo,"interS_hbscsc_the  %8.3f\n",WEIGHTS[48]);
  fprintf(fo,"interS_hbscsc_phi  %8.3f\n",WEIGHTS[49]);

  fprintf(fo,"interD_vdwatt      %8.3f\n",WEIGHTS[51]);
  fprintf(fo,"interD_vdwrep      %8.3f\n",WEIGHTS[52]);
  fprintf(fo,"interD_electr      %8.3f\n",WEIGHTS[53]);
  fprintf(fo,"interD_deslvP      %8.3f\n",WEIGHTS[54]);
  fprintf(fo,"interD_deslvH      %8.3f\n",WEIGHTS[55]);
  fprintf(fo,"interD_ssbond      %8.3f\n",WEIGHTS[56]);
  fprintf(fo,"interD_hbbbbb_dis  %8.3f\n",WEIGHTS[61]);
  fprintf(fo,"interD_hbbbbb_the  %8.3f\n",WEIGHTS[62]);
  fprintf(fo,"intraD_hbbbbb_phi  %8.3f\n",WEIGHTS[63]);
  fprintf(fo,"interD_hbscbb_dis  %8.3f\n",WEIGHTS[64]);
  fprintf(fo,"interD_hbscbb_the  %8.3f\n",WEIGHTS[65]);
  fprintf(fo,"interD_hbscbb_phi  %8.3f\n",WEIGHTS[66]);
  fprintf(fo,"interD_hbscsc_dis  %8.3f\n",WEIGHTS[67]);
  fprintf(fo,"interD_hbscsc_the  %8.3f\n",WEIGHTS[68]);
  fprintf(fo,"interD_hbscsc_phi  %8.3f\n",WEIGHTS[69]);

  fprintf(fo,"ligand_vdwatt      %8.3f\n",WEIGHTS[71]);
  fprintf(fo,"ligand_vdwrep      %8.3f\n",WEIGHTS[72]);
  fprintf(fo,"ligand_electr      %8.3f\n",WEIGHTS[73]);
  fprintf(fo,"ligand_deslvP      %8.3f\n",WEIGHTS[74]);
  fprintf(fo,"ligand_deslvH      %8.3f\n",WEIGHTS[75]);
  fprintf(fo,"ligand_hbscbb_dis  %8.3f\n",WEIGHTS[84]);
  fprintf(fo,"ligand_hbscbb_the  %8.3f\n",WEIGHTS[85]);
  fprintf(fo,"ligand_hbscbb_phi  %8.3f\n",WEIGHTS[86]);
  fprintf(fo,"ligand_hbscsc_dis  %8.3f\n",WEIGHTS[87]);
  fprintf(fo,"ligand_hbscsc_the  %8.3f\n",WEIGHTS[88]);
  fprintf(fo,"ligand_hbscsc_phi  %8.3f\n",WEIGHTS[89]);
  fclose(fo);
  return Success;
}


double CalcResidueBuriedRatio(Residue* pResi1){
  if(pResi1->nCbIn10A >= 15) return 1.0;
  else if(pResi1->nCbIn10A <= 8) return 0.0;
  else return (pResi1->nCbIn10A-8.0)/7.0;
}

double CalcAverageBuriedRatio(double ratio1, double ratio2){
  return (ratio1+ratio2)/2.0;
}

#define DEAL_WITH_BOND_CONNECTION
//these functions are used to check the 12, 13, 14 and 15 bond connectivity
BOOL ResidueIntraBond12Check(char *atom1, char *atom2, BondSet* pBondSet){
  for(int i = 0; i < pBondSet->count; i++){
    Bond *pBond = pBondSet->bonds+i;
    if( (strcmp(atom1, pBond->atomFromName) == 0 && strcmp(atom2, pBond->atomToName) == 0) ||
      (strcmp(atom2, pBond->atomFromName) == 0 && strcmp(atom1, pBond->atomToName) == 0) ){
      return TRUE;
    }
  }
  return FALSE;
}

BOOL ResidueIntraBond13Check(char *atom1, char *atom2, BondSet* pBondSet){
  for(int i = 0; i < pBondSet->count; i++){
    Bond *pBond = pBondSet->bonds+i;
    if( strcmp(atom1, pBond->atomFromName) == 0){
      if(ResidueIntraBond12Check(pBond->atomToName, atom2, pBondSet)){
        return TRUE;
      }
    }
    else if(strcmp(atom1, pBond->atomToName) == 0){
      if(ResidueIntraBond12Check(pBond->atomFromName, atom2, pBondSet)){
        return TRUE;
      }
    }
  }
  return FALSE;
}

BOOL ResidueIntraBond14Check(char *atom1, char *atom2, BondSet* pBondSet){
  for(int i = 0; i < pBondSet->count; i++){
    Bond *pBond = pBondSet->bonds+i;
    if( strcmp(atom1, pBond->atomFromName) == 0){
      if( ResidueIntraBond13Check(pBond->atomToName, atom2, pBondSet) ){
        return TRUE;
      }
    }
    else if( strcmp(atom1, pBond->atomToName) == 0 ){
      if( ResidueIntraBond13Check(pBond->atomFromName, atom2, pBondSet) ){
        return TRUE;
      }
    }
  }
  return FALSE;
}

int ResidueIntraBondConnectionCheck(char *atom1, char *atom2, BondSet* pBondSet){
  if(ResidueIntraBond12Check(atom1, atom2, pBondSet)){
    return 12;
  }
  else if(ResidueIntraBond13Check(atom1, atom2, pBondSet)){
    return 13;
  }
  else if(ResidueIntraBond14Check(atom1, atom2, pBondSet)){
    return 14;
  }
  else{
    return 15;
  }
}

int ResidueAndNextResidueInterBondConnectionCheck_charmm22(char *atomOnPreResi, char *atomOnNextResi, Residue *pPreResi, Residue *pNextResi){
  if( strcmp(atomOnPreResi, "C") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 12;
    else if( strcmp(atomOnNextResi, "CA") == 0 || strcmp(atomOnNextResi, "HN") == 0 || (strcmp(atomOnNextResi, "CD") == 0 && strcmp(pNextResi->name, "PRO") == 0) ){
      return 13;
    }
    else if( strcmp(atomOnNextResi, "CB") == 0|| strcmp(atomOnNextResi, "C") == 0 || 
      strcmp(atomOnNextResi, "HA") == 0 || strcmp(atomOnNextResi, "HA1") == 0 || strcmp(atomOnNextResi, "HA2") == 0 || 
      (strcmp(atomOnNextResi, "HD1") == 0 && strcmp(pNextResi->name, "PRO") == 0) || 
      (strcmp(atomOnNextResi, "HD2") == 0 && strcmp(pNextResi->name, "PRO") == 0)|| 
      (strcmp(atomOnNextResi, "CG") == 0 && strcmp(pNextResi->name, "PRO") == 0)){
      return 14;
    }
  }
  else if( strcmp(atomOnPreResi, "O") == 0 || strcmp(atomOnPreResi, "CA") == 0 ){
    if(strcmp(atomOnNextResi, "N") == 0) return 13;
    else if( strcmp(atomOnNextResi, "CA") == 0 || strcmp(atomOnNextResi, "HN") == 0|| 
      (strcmp(atomOnNextResi, "CD") == 0 && strcmp(pNextResi->name, "PRO") == 0) ){
      return 14;
    }
  }
  else if( strcmp(atomOnPreResi, "CB") == 0 || strcmp(atomOnPreResi, "HA") == 0 || 
    strcmp(atomOnPreResi, "HA1") == 0 || strcmp(atomOnPreResi, "HA2") == 0 || strcmp(atomOnPreResi, "N") == 0 ){
    if(strcmp(atomOnNextResi, "N") == 0) return 14;
  }
  return 15;
}

int ResidueAndNextResidueInterBondConnectionCheck_charmm19(char *atomOnPreResi, char *atomOnNextResi, char *nextResiName){
  if( strcmp(atomOnPreResi, "C") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 12;
    else if( strcmp(atomOnNextResi, "CA") == 0|| strcmp(atomOnNextResi, "H") == 0|| (strcmp(atomOnNextResi, "CD") == 0 && strcmp(nextResiName, "PRO") == 0) )
      return 13;
    else if( strcmp(atomOnNextResi, "CB") == 0|| strcmp(atomOnNextResi, "C") == 0 || (strcmp(atomOnNextResi, "CG") == 0 && strcmp(nextResiName, "PRO") == 0))
      return 14;
  }
  else if( strcmp(atomOnPreResi, "O") == 0 || strcmp(atomOnPreResi, "CA") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 13;
    else if( strcmp(atomOnNextResi, "CA") == 0|| strcmp(atomOnNextResi, "H") == 0||
      (strcmp(atomOnNextResi, "CD") == 0 && strcmp(nextResiName, "PRO") == 0) ){
      return 14;
    }
  }
  else if( strcmp(atomOnPreResi, "CB") == 0 || strcmp(atomOnPreResi, "N") == 0 ){
    if( strcmp(atomOnNextResi, "N") == 0 ) return 14;
  }
  return 15;
}




#define DEAL_WITH_BASIC_ENERGY_TERMS
//////////////////////////////////////////////////////////////
//Atomic pairwise energy for residue
//////////////////////////////////////////////////////////////
int VdwAttEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance, int bondType, double *vdwAtt){
  if(distance>=ENERGY_DISTANCE_CUTOFF) return Success;
  if(bondType==12||bondType==13) return Success;
  if(AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2)) return Success;
  double rmin = RADIUS_SCALE_FOR_VDW * (pAtom1->vdw_radius + pAtom2->vdw_radius);
  double ratio = distance/rmin;
  double energy=0.0;
  double scale=0.0;
  if(bondType==14) scale=ENERGY_SCALE_FACTOR_BOND_14;
  else if(bondType==15) scale=ENERGY_SCALE_FACTOR_BOND_15;

  if(ratio < 0.8909){
    energy=0.0;
  }
  else if(distance <= 5.0){
    double epsilon  = sqrt(pAtom1->vdw_epsilon*pAtom2->vdw_epsilon);
    double B6 = pow(1/ratio, 6.0);
    double A12 = B6*B6;
    energy = epsilon * (A12 - 2.0 * B6);
  }
  else if(distance > 5.0 && distance < ENERGY_DISTANCE_CUTOFF){
    double epsilon  = sqrt(pAtom1->vdw_epsilon*pAtom2->vdw_epsilon);
    double B6 = pow((double)rmin/5.0, 6.0);
    double A12 = B6 * B6;
    double M = epsilon * ( A12 - 2.0 * B6);
    double N = 2.4 * epsilon * (B6 - A12);
    double a = 2 * M + N;
    double b = -33 * M - 17 * N;
    double c = 180 * M + 96 * N;
    double d = -324 * M -180 * N;
    energy = a * distance * distance * distance + b * distance * distance + c * distance + d;
  }
  energy*=scale;
  *vdwAtt=energy;
  if(ENERGY_DEBUG_MODE_VDW_ATT){
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %f, ratio: %f, vdwAtt: %f\n", 
      AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1), 
      AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
      bondType,distance, ratio, energy);
  }
  return Success;
}



int VdwRepEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance,int bondType, double *vdwRep){
  //return Success;
  if(bondType==12||bondType==13) return Success;
  //if(AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2)) return Success;
  double rmin = RADIUS_SCALE_FOR_VDW * (pAtom1->vdw_radius + pAtom2->vdw_radius);
  double ratio = distance/rmin;
  double epsilon  = sqrt(pAtom1->vdw_epsilon*pAtom2->vdw_epsilon);
  double RATIO_CUTOFF = 0.70; // can be adjusted
  double energy=0.0;
  double scale=0.0;
  if(bondType==14) scale=ENERGY_SCALE_FACTOR_BOND_14;
  else if(bondType==15) scale=ENERGY_SCALE_FACTOR_BOND_15;

  if(ratio > 0.8909){ // 0.8909 ~ inf
    energy=0.0;
  }
  else if(ratio >= RATIO_CUTOFF){ // [0.70, 0.8909]
    double B6 = pow(1/ratio, 6.0);
    double A12 = B6*B6;
    energy = epsilon * (A12 - 2.0 * B6);
  }
  else{ // [0, 0.70]
    double B6_0 = pow(1/RATIO_CUTOFF, 6.0);
    double a = epsilon * (B6_0 * B6_0 - 2.0 * B6_0);
    double b = epsilon * 12.0 * (B6_0 / RATIO_CUTOFF - B6_0 * B6_0 / RATIO_CUTOFF);
    double y0 = a * epsilon;
    double k = b * epsilon;
    energy = k * (ratio - RATIO_CUTOFF) + y0;
  }
  energy*=scale;
  //set a cutoff for maximum clash
  //double MAX_CLASH=5.0*epsilon;
  //if(energy>MAX_CLASH) energy=MAX_CLASH;
  *vdwRep=energy;

  if(ENERGY_DEBUG_MODE_VDW_REP){
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %f, ratio: %f, vdwRep: %f\n", 
      AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1),
      AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
      bondType,distance, ratio, energy);
  }
  return Success;
}


int HBondEnergyAtomAndAtom(Atom *atomH, Atom *atomA, Atom*atomD, Atom *atomB,double distanceHA, int bondType, double *etotal, double *edist, double *etheta, double *ephi){
  //return Success;
  if(bondType==12||bondType==13) return Success;
  if(distanceHA > HBOND_DISTANCE_CUTOFF_MAX) return Success;
  XYZ xyzDH = XYZDifference(&atomD->xyz, &atomH->xyz);
  XYZ xyzHA = XYZDifference(&atomH->xyz, &atomA->xyz);
  XYZ xyzAAB = XYZDifference(&atomA->xyz, &atomB->xyz);
  double angleTheta = PI-XYZAngle(&xyzDH, &xyzHA);
  if(RadToDeg(angleTheta)<90) return Success;
  double anglePhi   = PI-XYZAngle(&xyzHA, &xyzAAB);
  if(RadToDeg(anglePhi)<80) return Success;

  double energyR=0.0;
  if(distanceHA<HBOND_OPTIMAL_DISTANCE){
    energyR=-1.0*HBOND_WELL_DEPTH*cos((distanceHA-HBOND_OPTIMAL_DISTANCE)*PI);
  }
  else{
    energyR=-0.5*cos(PI/(HBOND_DISTANCE_CUTOFF_MAX-HBOND_OPTIMAL_DISTANCE)*(distanceHA-HBOND_OPTIMAL_DISTANCE))-0.5;
  }
  if(energyR > 0.0) energyR = 0.0;

  double energyTheta = -1.0*cos(angleTheta)*cos(angleTheta)*cos(angleTheta)*cos(angleTheta);
  double energyPhi = 0.0;
  if(atomH->isBBAtom && atomA->isBBAtom){
    energyPhi=-1.0*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150));
  }
  else{
    if(atomA->hybridType == Type_AtomHybridType_SP3){
      energyPhi = -1.0*cos(anglePhi-DegToRad(135))*cos(anglePhi-DegToRad(135))*cos(anglePhi-DegToRad(135))*cos(anglePhi-DegToRad(135));
    }
    else if(atomA->hybridType == Type_AtomHybridType_SP2){
      energyPhi = -1.0*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150))*cos(anglePhi-DegToRad(150));
    }
  }
  
  // original function is error, we should add angle and restriction to calculate energy
  double energy = 0.0;
  if(RadToDeg(angleTheta) >= 90 && RadToDeg(anglePhi) >= 80 && distanceHA < HBOND_DISTANCE_CUTOFF_MAX){
    energy = energyR + energyTheta + energyPhi;
    if(energy>0.0) energy=0.0;
    *etotal = energy;
    *edist = energyR;
    *etheta = energyTheta;
    *ephi = energyPhi;
  }

  if(ENERGY_DEBUG_MODE_HBOND){
    printf("AtomH: %1s %4d %4s, AtomA: %1s %4d %4s, dist: %5.2f, theta: %5.1f, phi: %5.1f, edist: %5.2f, ethe: %5.2f, ephi: %5.2f\n", 
      AtomGetChainName(atomH), AtomGetPosInChain(atomH), AtomGetName(atomH),
      AtomGetChainName(atomA), AtomGetPosInChain(atomA), AtomGetName(atomA),
      distanceHA, RadToDeg(angleTheta),RadToDeg(anglePhi),energyR,energyTheta,energyPhi);
  }
  return Success;
}


int ElecEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance12,int bondType, double *elec){
  //return Success;
  if(bondType==12||bondType==13) return Success;
  if(distance12 > ELEC_DISTANCE_CUTOFF) return Success;
  if(fabs(pAtom1->charge)<1e-2 || fabs(pAtom2->charge)<1e-2) return Success;
  else if(distance12<0.8*(pAtom1->vdw_radius + pAtom2->vdw_radius)) distance12=0.8*(pAtom1->vdw_radius + pAtom2->vdw_radius);

  double energy = COULOMB_CONSTANT*pAtom1->charge*pAtom2->charge/distance12/distance12/40.0;
  double scale=0.0;
  if(bondType==14) scale=ENERGY_SCALE_FACTOR_BOND_14;
  else if(bondType==15) scale=ENERGY_SCALE_FACTOR_BOND_15;
  energy*=scale;
  *elec = energy;

  if(ENERGY_DEBUG_MODE_ELEC){
    if(fabs(energy) > 0.0 && distance12 < 4.0){
      printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %5.2f, elec: %5.2f\n", 
        AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1), 
        AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
        bondType, distance12, energy);
    }
  }
  return Success;
}

int LKDesolvationEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance,int bondType, double *energyP, double *energyH){
  //return Success;
  if(bondType==12||bondType==13) return Success;
  if(AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2)) return Success;
  if(distance>ENERGY_DISTANCE_CUTOFF) return Success;
  double volume1 = pAtom1->EEF1_volume;
  double volume2 = pAtom2->EEF1_volume;
  double dGFreeAtom1 = pAtom1->EEF1_freeDG;
  double dGFreeAtom2 = pAtom2->EEF1_freeDG;
  double coefficient = -0.089793561062582974; // 0.5/(pi^1.5)
  double r1=pAtom1->vdw_radius*RADIUS_SCALE_FOR_DESOLV;
  double r2=pAtom2->vdw_radius*RADIUS_SCALE_FOR_DESOLV;
  double r12 = r1+r2;

  distance = distance < r12 ? r12 : distance;
  double lamda1 = pAtom1->EEF1_lamda_ * distance * distance;
  double lamda2 = pAtom2->EEF1_lamda_ * distance * distance;
  double x1 = (distance - r1)/pAtom1->EEF1_lamda_;
  double x2 = (distance - r2)/pAtom2->EEF1_lamda_;

  double desolv12 = coefficient * volume2 * dGFreeAtom1 / lamda1;
  desolv12 *= exp( -1.0 * x1 * x1 );
  double desolv21 = coefficient * volume1 * dGFreeAtom2 / lamda2;
  desolv21 *= exp( -1.0 * x2 * x2 );
  if(pAtom1->polarity == Type_AtomPolarity_P || pAtom1->polarity == Type_AtomPolarity_C) *energyP += desolv12;
  else *energyH += desolv12;
  if(pAtom2->polarity == Type_AtomPolarity_P || pAtom2->polarity == Type_AtomPolarity_C) *energyP += desolv21;
  else *energyH += desolv21;

  if(ENERGY_DEBUG_MODE_DESOLV){
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, dist: %5.2f, desolv12: %7.3f, desolv21: %7.3f\n",
      AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1), 
      AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
      distance, desolv12, desolv21);
  }

  return Success;
}

///////////////////////////////////////////////////////////////////
//User can define new energy term here
///////////////////////////////////////////////////////////////////

int SSbondEnergyAtomAndAtom(Atom *pAtomS1,Atom *pAtomS2,Atom *pAtomCB1,Atom *pAtomCB2,Atom* pAtomCA1,Atom* pAtomCA2,double* sse){
  XYZ xyzC1S1=XYZDifference(&pAtomCB1->xyz,&pAtomS1->xyz);
  XYZ xyzS1S2=XYZDifference(&pAtomS1->xyz,&pAtomS2->xyz);
  XYZ xyzS2C2=XYZDifference(&pAtomS2->xyz,&pAtomCB2->xyz);
  double dSS=XYZDistance(&pAtomS1->xyz,&pAtomS2->xyz);
  double aC1S1S2=RadToDeg(PI-XYZAngle(&xyzS1S2,&xyzS2C2));
  double aC2S2S1=RadToDeg(PI-XYZAngle(&xyzC1S1,&xyzS1S2));
  double xC1S1S2C2=GetTorsionAngle(&pAtomCB1->xyz,&pAtomS1->xyz,&pAtomS2->xyz,&pAtomCB2->xyz);
  double xCA1CB1SG1SG2=GetTorsionAngle(&pAtomCA1->xyz,&pAtomCB1->xyz,&pAtomS1->xyz,&pAtomS2->xyz);
  double xCA2CB2SG2SG1=GetTorsionAngle(&pAtomCA2->xyz,&pAtomCB2->xyz,&pAtomS2->xyz,&pAtomS1->xyz);
  *sse=(0.8*(1-exp(-10.0*(dSS-SSBOND_DISTANCE)))*(1-exp(-10.0*(dSS-SSBOND_DISTANCE))))
    +0.005*(aC1S1S2-SSBOND_ANGLE)*(aC1S1S2-SSBOND_ANGLE)
    +0.005*(aC2S2S1-SSBOND_ANGLE)*(aC2S2S1-SSBOND_ANGLE)
    +cos(2.0*xC1S1S2C2)+1.0
    +1.25*sin(xCA1CB1SG1SG2+2.0*PI/3.0)-1.75
    +1.25*sin(xCA2CB2SG2SG1+2.0*PI/3.0)-1.75;
  if(*sse>0.0) *sse=0.0;
  return Success;
}



#define RESIDUE_PAIRWISE_ENERGY_COMPUTATION
int AminoAcidReferenceEnergy(char *AAname, double energyTerm[MAX_ENERGY_TERM]){
  if(strcmp(AAname, "ALA")  == 0)      energyTerm[ 1] += 1.0;
  else if(strcmp(AAname, "CYS")  == 0) energyTerm[ 2] += 1.0;
  else if(strcmp(AAname, "ASP")  == 0) energyTerm[ 3] += 1.0;
  else if(strcmp(AAname, "GLU")  == 0) energyTerm[ 4] += 1.0;
  else if(strcmp(AAname, "PHE")  == 0) energyTerm[ 5] += 1.0;
  else if(strcmp(AAname, "GLY")  == 0) energyTerm[ 6] += 1.0;
  else if(strcmp(AAname, "HIS")  == 0) energyTerm[ 7] += 1.0;
  else if(strcmp(AAname, "HSE")  == 0) energyTerm[ 7] += 1.0;
  else if(strcmp(AAname, "HSD")  == 0) energyTerm[ 7] += 1.0;
  else if(strcmp(AAname, "ILE")  == 0) energyTerm[ 8] += 1.0;
  else if(strcmp(AAname, "LYS")  == 0) energyTerm[ 9] += 1.0;
  else if(strcmp(AAname, "LEU")  == 0) energyTerm[10] += 1.0;
  else if(strcmp(AAname, "MET")  == 0) energyTerm[11] += 1.0;
  else if(strcmp(AAname, "ASN")  == 0) energyTerm[12] += 1.0;
  else if(strcmp(AAname, "PRO")  == 0) energyTerm[13] += 1.0;
  else if(strcmp(AAname, "GLN")  == 0) energyTerm[14] += 1.0;
  else if(strcmp(AAname, "ARG")  == 0) energyTerm[15] += 1.0;
  else if(strcmp(AAname, "SER")  == 0) energyTerm[16] += 1.0;
  else if(strcmp(AAname, "THR")  == 0) energyTerm[17] += 1.0;
  else if(strcmp(AAname, "VAL")  == 0) energyTerm[18] += 1.0;
  else if(strcmp(AAname, "TRP")  == 0) energyTerm[19] += 1.0;
  else if(strcmp(AAname, "TYR")  == 0) energyTerm[20] += 1.0;
  return Success;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// the following are energy between residue and residue
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int EnergyIntraResidue(Residue* pThis,double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); ++i){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    if(pAtom1->isXyzValid==FALSE)continue;
    for(int j=i+1; j<ResidueGetAtomCount(pThis);++j){
      Atom* pAtom2=ResidueGetAtom(pThis,j);
      if(pAtom2->isXyzValid==FALSE)continue;
      if(pAtom2->isBBAtom && pAtom1->isBBAtom) continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
        if(strcmp(ResidueGetName(pThis),"ILE")==0 || strcmp(ResidueGetName(pThis),"MET")==0 ||
          strcmp(ResidueGetName(pThis),"GLN")==0||strcmp(ResidueGetName(pThis),"GLU")==0||
          strcmp(ResidueGetName(pThis),"LYS")==0||strcmp(ResidueGetName(pThis),"ARG")==0){
            int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetBonds(pThis));
            if(bondType==12||bondType==13) continue;
            double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0;
            VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
            VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
            LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
            energyTerms[21]+=vdwAtt;
            energyTerms[22]+=vdwRep;
            energyTerms[24]+=desolvP;
            energyTerms[25]+=desolvH;
        }
      }
      else if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){
        if(strcmp(AtomGetName(pAtom1),"CB")==0 || strcmp(AtomGetName(pAtom2),"CB")==0) continue;
        int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetBonds(pThis));
        if(bondType==12||bondType==13) continue;
        double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0,ele=0;
        VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
        VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
        LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
        ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
        energyTerms[21]+=vdwAtt;
        energyTerms[22]+=vdwRep;
        energyTerms[23]+=ele;
        energyTerms[24]+=desolvP;
        energyTerms[25]+=desolvH;
        if(distance < HBOND_DISTANCE_CUTOFF_MAX){
          double hbtot=0,hbd=0,hbt=0,hbp=0;
          if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
            HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom2)),
              distance,bondType,&hbtot,&hbd,&hbt,&hbp);
          }
          else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
            HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
              distance,bondType,&hbtot,&hbd,&hbt,&hbp);
          }
          energyTerms[26]+=hbd;
          energyTerms[27]+=hbt;
          energyTerms[28]+=hbp;
        }
      }
    }

  }  
  return Success;
}

int EnergyResidueAndNextResidue(Residue* pThis,Residue* pOther,double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0;i<ResidueGetAtomCount(pThis);i++){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    if(pAtom1->isXyzValid==FALSE)continue;
    for(int j=0;j<ResidueGetAtomCount(pOther);j++){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      if(pAtom2->isXyzValid==FALSE)continue;
      if(pAtom1->isBBAtom && pAtom2->isBBAtom) continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
        int bondType=15;
        double att=0,rep=0,desP=0,desH=0,ele=0;
        VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
        VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
        LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
        ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
        energyTerms[31]+=att;
        energyTerms[32]+=rep;
        energyTerms[33]+=ele;
        energyTerms[34]+=desP;
        energyTerms[35]+=desH;
        if(distance<HBOND_DISTANCE_CUTOFF_MAX){
          double hbtot=0,hbd=0,hbt=0,hbp=0;
          if(pAtom1->isHBatomH && pAtom2->isHBatomA){
            HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
              distance,bondType,&hbtot,&hbd,&hbt,&hbp);
          }
          else if(pAtom2->isHBatomH && pAtom1->isHBatomA){
            HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
              distance,bondType,&hbtot,&hbd,&hbt,&hbp);
          }
          if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
            hbd *= HBOND_LOCAL_REDUCE;
            hbt *= HBOND_LOCAL_REDUCE;
            hbp *= HBOND_LOCAL_REDUCE;
          }
          if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
            energyTerms[41]+=hbd;
            energyTerms[42]+=hbt;
            energyTerms[43]+=hbp;
          }
          else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
            energyTerms[47]+=hbd;
            energyTerms[48]+=hbt;
            energyTerms[49]+=hbp;
          }
          else{
            energyTerms[44]+=hbd;
            energyTerms[45]+=hbt;
            energyTerms[46]+=hbp;
          }
        }
      }
      else if((pAtom1->isBBAtom && pAtom2->isBBAtom==FALSE) || (pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom)){
        int bondType=ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetName(pOther));
        if(bondType==12||bondType==13)continue;
        double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0,ele=0;
        VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
        VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
        LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
        ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
        energyTerms[31]+=vdwAtt;
        energyTerms[32]+=vdwRep;
        energyTerms[33]+=ele;
        energyTerms[34]+=desolvP;
        energyTerms[35]+=desolvH;
        if(distance < HBOND_DISTANCE_CUTOFF_MAX){
          double hbtot=0,hbd=0,hbt=0,hbp=0;
          if(pAtom1->isHBatomH && pAtom2->isHBatomA){
            HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
              distance,bondType,&hbtot,&hbd,&hbt,&hbp);
          }
          else if(pAtom2->isHBatomH && pAtom1->isHBatomA){
            HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
              distance,bondType,&hbtot,&hbd,&hbt,&hbp);
          }
          if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
            hbd *= HBOND_LOCAL_REDUCE;
            hbt *= HBOND_LOCAL_REDUCE;
            hbp *= HBOND_LOCAL_REDUCE;
          }
          if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
            energyTerms[41]+=hbd;
            energyTerms[42]+=hbt;
            energyTerms[43]+=hbp;
          }
          else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
            energyTerms[47]+=hbd;
            energyTerms[48]+=hbt;
            energyTerms[49]+=hbp;
          }
          else{
            energyTerms[44]+=hbd;
            energyTerms[45]+=hbt;
            energyTerms[46]+=hbp;
          }
        }
      }
    }
  }
  return Success;
}

int EnergyResidueAndOtherResidueSameChain(Residue* pThis,Residue* pOther,double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0;i<ResidueGetAtomCount(pThis);i++){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    if(pAtom1->isXyzValid==FALSE)continue;
    for(int j=0;j<ResidueGetAtomCount(pOther);j++){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      if(pAtom2->isXyzValid==FALSE)continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      int bondType=15;
      double vdwAtt=0,vdwRep=0,ele=0,desolvP=0,desolvH=0;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
      energyTerms[31]+=vdwAtt;
      energyTerms[32]+=vdwRep;
      energyTerms[33]+=ele;
      energyTerms[34]+=desolvP;
      energyTerms[35]+=desolvH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0,hbd=0,hbt=0,hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom2->isHBatomH && pAtom1->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
          hbd*=HBOND_LOCAL_REDUCE;
          hbt*=HBOND_LOCAL_REDUCE;
          hbp*=HBOND_LOCAL_REDUCE;
        }
        if(pAtom1->isBBAtom && pAtom2->isBBAtom){
          energyTerms[41]+=hbd;
          energyTerms[42]+=hbt;
          energyTerms[43]+=hbp;
        }
        else if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
          energyTerms[47]+=hbd;
          energyTerms[48]+=hbt;
          energyTerms[49]+=hbp;
        }
        else{
          energyTerms[44]+=hbd;
          energyTerms[45]+=hbt;
          energyTerms[46]+=hbp;
        }
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",ResidueGetName(pThis))==0 && strcmp("CYS",ResidueGetName(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<ResidueGetAtomCount(pThis);i++){
      Atom* pAtom=ResidueGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<ResidueGetAtomCount(pOther);i++){
      Atom* pAtom=ResidueGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL && pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[36]+=ssbond;
      }
    }
  }

  return Success;
}


int EnergyResidueAndOtherResidueDiffChain(Residue* pThis,Residue* pOther,double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<ResidueGetAtomCount(pThis); i++){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    if(!pAtom1->isXyzValid)continue;
    for(int j=0; j<ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      if(!pAtom2->isXyzValid)continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0,rep=0,ele=0,desP=0,desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[51]+=att;
      energyTerms[52]+=rep;
      energyTerms[53]+=ele;
      energyTerms[54]+=desP;
      energyTerms[55]+=desH;
      if(distance < HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0,hbd=0,hbt=0,hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        if(pAtom1->isBBAtom && pAtom2->isBBAtom){
          energyTerms[61]+=hbd;
          energyTerms[62]+=hbt;
          energyTerms[63]+=hbp;
        }
        else if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
          energyTerms[67]+=hbd;
          energyTerms[68]+=hbt;
          energyTerms[69]+=hbp;
        }
        else{
          energyTerms[64]+=hbd;
          energyTerms[65]+=hbt;
          energyTerms[66]+=hbp;
        }
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",ResidueGetName(pThis))==0 && strcmp("CYS",ResidueGetName(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<ResidueGetAtomCount(pThis);i++){
      Atom* pAtom=ResidueGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<ResidueGetAtomCount(pOther);i++){
      Atom* pAtom=ResidueGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL &&pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[56]+=ssbond;
      }
    }
  }

  return Success;
}


int EnergyResidueAndLigResidue(Residue* pProtein, Residue* pLigand,double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<ResidueGetAtomCount(pProtein); i++){
    Atom* pAtom1=ResidueGetAtom(pProtein,i);
    if(pAtom1->isXyzValid==FALSE)continue;
    for(int j=0; j<ResidueGetAtomCount(pLigand); j++){
      Atom* pAtom2=ResidueGetAtom(pLigand,j);
      if(pAtom2->isXyzValid==FALSE)continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      int bondType=15;
      double vdwAtt=0,vdwRep=0,ele=0,desolvP=0,desolvH=0;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
      energyTerms[71]+=vdwAtt;
      energyTerms[72]+=vdwRep;
      energyTerms[73]+=ele;
      energyTerms[74]+=desolvP;
      energyTerms[75]+=desolvH;
      if(distance < HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0,hbd=0,hbt=0,hbp=0;
        if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pProtein,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pLigand,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pLigand,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pProtein,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
          energyTerms[84]+=hbd;
          energyTerms[85]+=hbt;
          energyTerms[86]+=hbp;
        }
        else{
          energyTerms[81]+=hbd;
          energyTerms[82]+=hbt;
          energyTerms[83]+=hbp;
        }
      }
    }
  }
  return Success;
}




#define ROTAMER_PAIRWISE_ENERGY_COMPUTATION
////////////////////////////////////////////////////////////////////////////////////////////
//these functions are used to calculate energy for rotamers
////////////////////////////////////////////////////////////////////////////////////////////

int EnergyIntraRotamer(Rotamer* pThis,double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); ++i){
    Atom* pAtom1=RotamerGetAtom(pThis,i);
    if(pAtom1->isXyzValid==FALSE) continue;
    for(int j=i+1; j<RotamerGetAtomCount(pThis);++j){
      Atom* pAtom2=RotamerGetAtom(pThis,j);
      if(pAtom2->isXyzValid==FALSE) continue;
      if(pAtom2->isBBAtom && pAtom1->isBBAtom) continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
        if(strcmp(RotamerGetType(pThis),"ILE")==0 || strcmp(RotamerGetType(pThis),"MET")==0 ||
          strcmp(RotamerGetType(pThis),"GLN")==0||strcmp(RotamerGetType(pThis),"GLU")==0||
          strcmp(RotamerGetType(pThis),"LYS")==0||strcmp(RotamerGetType(pThis),"ARG")==0){
            int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),RotamerGetBonds(pThis));
            if(bondType==12||bondType==13) continue;
            double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0;
            VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
            VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
            LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
            energyTerms[21]+=vdwAtt;
            energyTerms[22]+=vdwRep;
            energyTerms[24]+=desolvP;
            energyTerms[25]+=desolvH;
        }
      }
      else if((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE)){
        if(strcmp(AtomGetName(pAtom1),"CB")==0 || strcmp(AtomGetName(pAtom2),"CB")==0) continue;
        int bondType=ResidueIntraBondConnectionCheck(AtomGetName(pAtom1),AtomGetName(pAtom2),RotamerGetBonds(pThis));
        if(bondType==12||bondType==13) continue;
        double vdwAtt=0,vdwRep=0,desolvP=0,desolvH=0,ele=0;
        VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwAtt);
        VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&vdwRep);
        LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desolvP,&desolvH);
        ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
        energyTerms[21]+=vdwAtt;
        energyTerms[22]+=vdwRep;
        energyTerms[23]+=ele;
        energyTerms[24]+=desolvP;
        energyTerms[25]+=desolvH;
        if(distance < HBOND_DISTANCE_CUTOFF_MAX){
          double hbtot=0,hbd=0,hbt=0,hbp=0;
          if(pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE){
            HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom2)),
              distance,bondType,&hbtot,&hbd,&hbt,&hbp);
          }
          else if(pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE){
            HBondEnergyAtomAndAtom(pAtom2,pAtom1,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
              distance,bondType,&hbtot,&hbd,&hbt,&hbp);
          }
          energyTerms[26]+=hbd;
          energyTerms[27]+=hbt;
          energyTerms[28]+=hbp;
        }
      }
    }
  }
  return Success;
}


int EnergyRotamerAndRotamerSameChain(Rotamer* pThis,Rotamer* pOther,double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0;i<RotamerGetAtomCount(pThis);i++){
    Atom* pAtom1=RotamerGetAtom(pThis,i);
    if(pAtom1->isXyzValid==FALSE) continue;
  	if(pAtom1->isBBAtom) continue;
    for(int j=0;j<RotamerGetAtomCount(pOther);j++){
      Atom* pAtom2=RotamerGetAtom(pOther,j);
      if(pAtom2->isXyzValid==FALSE) continue;
	    if(pAtom2->isBBAtom) continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      int bondType=15;
      double att=0,rep=0,ele=0,desP=0,desH=0;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[31]+=att;
      energyTerms[32]+=rep;
      energyTerms[33]+=ele;
      energyTerms[34]+=desP;
      energyTerms[35]+=desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0,hbd=0,hbt=0,hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
          hbd*=HBOND_LOCAL_REDUCE;
          hbt*=HBOND_LOCAL_REDUCE;
          hbp*=HBOND_LOCAL_REDUCE;
        }
        /*if(pAtom1->isBBAtom && pAtom2->isBBAtom){
          energyTerms[41]+=hbd;
          energyTerms[42]+=hbt;
          energyTerms[43]+=hbp;
        }
        else if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){*/
          energyTerms[47]+=hbd;
          energyTerms[48]+=hbt;
          energyTerms[49]+=hbp;
        /*}
        else{
          energyTerms[44]+=hbd;
          energyTerms[45]+=hbt;
          energyTerms[46]+=hbp;
        }*/
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",RotamerGetType(pThis))==0 && strcmp("CYS",RotamerGetType(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<RotamerGetAtomCount(pThis);i++){
      Atom* pAtom=RotamerGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<RotamerGetAtomCount(pOther);i++){
      Atom* pAtom=RotamerGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL && pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[36]+=ssbond;
      }
    }
  }

  return Success;
}

int EnergyRotamerAndRotamerDiffChain(Rotamer* pThis, Rotamer* pOther, double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isXyzValid==FALSE) continue;
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <RotamerGetAtomCount(pOther); j++){
      Atom* pAtom2 = RotamerGetAtom(pOther, j);
      if(pAtom2->isXyzValid==FALSE) continue;
      if(pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[51] += att;
      energyTerms[52] += rep;
      energyTerms[53] += ele;
      energyTerms[54] += desP;
      energyTerms[55] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0,hbd=0,hbt=0,hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        /*if(pAtom1->isBBAtom && pAtom2->isBBAtom){
          energyTerms[61]+=hbd;
          energyTerms[62]+=hbt;
          energyTerms[63]+=hbp;
        }
        else if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){*/
          energyTerms[67]+=hbd;
          energyTerms[68]+=hbt;
          energyTerms[69]+=hbp;
        /*}
        else{
          energyTerms[64]+=hbd;
          energyTerms[65]+=hbt;
          energyTerms[66]+=hbp;
        }*/
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",RotamerGetType(pThis))==0 && strcmp("CYS",RotamerGetType(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<RotamerGetAtomCount(pThis);i++){
      Atom* pAtom=RotamerGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<RotamerGetAtomCount(pOther);i++){
      Atom* pAtom=RotamerGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL && pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[56]+=ssbond;
      }
    }
  }
  return Success;
}



int EnergyRotamerAndLigRotamer(Rotamer* pThis,Rotamer* pOther,double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isXyzValid==FALSE) continue;
	if(pAtom1->isBBAtom) continue;
    for(int j=0; j <RotamerGetAtomCount(pOther); j++){
      Atom* pAtom2 = RotamerGetAtom(pOther, j);
      if(pAtom2->isXyzValid==FALSE) continue;
	  if(pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,RotamerGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        energyTerms[84]+=hbd;
        energyTerms[85]+=hbt;
        energyTerms[86]+=hbp;
      }
    }
  }

  return Success;
}



int EnergyRotamerAndFixedResidueSameChain(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM]){
  int neighborCheck = 0;
  if(RotamerGetPosInChain(pThis)+1==ResidueGetPosInChain(pOther)) neighborCheck=12;
  else if(RotamerGetPosInChain(pThis)-1==ResidueGetPosInChain(pOther)) neighborCheck=21;

  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isXyzValid==FALSE) continue;
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if(pAtom2->isXyzValid==FALSE) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      int bondType=15;
      if(neighborCheck==12) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetName(pOther));
      else if(neighborCheck==21) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom2),AtomGetName(pAtom1),RotamerGetType(pThis));
      if(bondType==12 || bondType==13) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[31] += att;
      energyTerms[32] += rep;
      energyTerms[33] += ele;
      energyTerms[34] += desP;
      energyTerms[35] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
          hbd *= HBOND_LOCAL_REDUCE;
          hbt *= HBOND_LOCAL_REDUCE;
          hbp *= HBOND_LOCAL_REDUCE;
        }
        /*if(pAtom1->isBBAtom && pAtom2->isBBAtom){
          energyTerms[41]+=hbd;
          energyTerms[42]+=hbt;
          energyTerms[43]+=hbp;
        }
        else */if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){
          energyTerms[47]+=hbd;
          energyTerms[48]+=hbt;
          energyTerms[49]+=hbp;
        }
        else{
          energyTerms[44]+=hbd;
          energyTerms[45]+=hbt;
          energyTerms[46]+=hbp;
        }
        
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",RotamerGetType(pThis))==0 && strcmp("CYS",ResidueGetName(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<RotamerGetAtomCount(pThis);i++){
      Atom* pAtom=RotamerGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<ResidueGetAtomCount(pOther);i++){
      Atom* pAtom=ResidueGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL && pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[36]+=ssbond;
      }
    }
  }


  return Success;
}

int EnergyRotamerAndFixedResidueDiffChain(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(!pAtom1->isXyzValid) continue;
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if(!pAtom2->isXyzValid) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0,rep=0,ele=0,desP=0,desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[51] += att;
      energyTerms[52] += rep;
      energyTerms[53] += ele;
      energyTerms[54] += desP;
      energyTerms[55] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        /*if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
          energyTerms[61]+=hbd;
          energyTerms[62]+=hbt;
          energyTerms[63]+=hbp;
        }
        else */if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
          energyTerms[67]+=hbd;
          energyTerms[68]+=hbt;
          energyTerms[69]+=hbp;
        }
        else{
          energyTerms[64]+=hbd;
          energyTerms[65]+=hbt;
          energyTerms[66]+=hbp;
        }
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",RotamerGetType(pThis))==0 && strcmp("CYS",ResidueGetName(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<RotamerGetAtomCount(pThis);i++){
      Atom* pAtom=RotamerGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<ResidueGetAtomCount(pOther);i++){
      Atom* pAtom=ResidueGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL && pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[56]+=ssbond;
      }
    }
  }


  return Success;
}


int EnergyRotamerAndFixedLigResidue(Rotamer* pThis, Residue* pLigand, double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isXyzValid==FALSE) continue;
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pLigand); j++){
      Atom* pAtom2 = ResidueGetAtom(pLigand, j);
      if(pAtom2->isXyzValid==FALSE) continue;
      //if(pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pLigand,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pLigand,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        //if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
          energyTerms[84]+=hbd;
          energyTerms[85]+=hbt;
          energyTerms[86]+=hbp;
        /*}
        else{
          energyTerms[81]+=hbd;
          energyTerms[82]+=hbt;
          energyTerms[83]+=hbp;
        }*/
      }
    }
  }

  return Success;
}



int EnergyLigRotamerAndFixedResidue(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isXyzValid==FALSE) continue;
    for(int j=0; j <ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if(pAtom2->isXyzValid==FALSE) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
          energyTerms[84]+=hbd;
          energyTerms[85]+=hbt;
          energyTerms[86]+=hbp;
        }
        else{
          energyTerms[81]+=hbd;
          energyTerms[82]+=hbt;
          energyTerms[83]+=hbp;
        }
      }
    }
  }
  return Success;
}



int EnergyRotamerAndDesignResidueSameChain(Rotamer* pThis,Residue* pOther,double energyTerms[MAX_ENERGY_TERM]){
  int neighborCheck = 0;
  if(RotamerGetPosInChain(pThis)+1==ResidueGetPosInChain(pOther)) neighborCheck=12;
  else if(RotamerGetPosInChain(pThis)-1==ResidueGetPosInChain(pOther)) neighborCheck=21;

  for(int i=0;i<RotamerGetAtomCount(pThis);i++){
    Atom* pAtom1=RotamerGetAtom(pThis,i);
    if(pAtom1->isXyzValid==FALSE) continue;
    if(pAtom1->isBBAtom) continue;
    for(int j=0;j<ResidueGetAtomCount(pOther);j++){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      if(pAtom2->isXyzValid==FALSE) continue;
      if(pAtom2->isBBAtom==FALSE) continue;
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      int bondType=15;
      if(neighborCheck==12) bondType=ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetName(pOther));
      else if(neighborCheck==21) bondType=ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom2),AtomGetName(pAtom1),RotamerGetType(pThis));
      if(bondType==12 || bondType==13) continue;
      double att=0,rep=0,ele=0,desP=0,desH=0;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[31] += att;
      energyTerms[32] += rep;
      energyTerms[33] += ele;
      energyTerms[34] += desP;
      energyTerms[35] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0,hbd=0,hbt=0,hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
          hbd *= HBOND_LOCAL_REDUCE;
          hbt *= HBOND_LOCAL_REDUCE;
          hbp *= HBOND_LOCAL_REDUCE;
        }
        energyTerms[44] += hbd;
        energyTerms[45] += hbt;
        energyTerms[46] += hbp;       
      }
    }
  }

  return Success;
}



int EnergyRotamerAndDesignResidueDiffChain(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis,i);
    if(!pAtom1->isXyzValid) continue;
    if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2 = ResidueGetAtom(pOther,j);
      if(!pAtom2->isXyzValid) continue;
      if(!pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0,rep=0,ele=0,desP=0,desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[51] += att;
      energyTerms[52] += rep;
      energyTerms[53] += ele;
      energyTerms[54] += desP;
      energyTerms[55] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0,hbd=0,hbt=0,hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        energyTerms[64] += hbd;
        energyTerms[65] += hbt;
        energyTerms[66] += hbp;
      }
    }
  }
  return Success;
}





int EnergyLigRotamerAndDesignResidue(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<RotamerGetAtomCount(pThis); i++){
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if(pAtom1->isXyzValid==FALSE) continue;
	if(pAtom1->isBBAtom) continue;
    for(int j=0; j <ResidueGetAtomCount(pOther); j++){
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if(pAtom2->isXyzValid==FALSE) continue;
      if(pAtom2->isBBAtom==FALSE) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      double att=0, rep=0, ele=0, desP=0, desH=0;
      int bondType=15;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0, hbd=0, hbt=0, hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),RotamerGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        energyTerms[84] += hbd;
        energyTerms[85] += hbt;
        energyTerms[86] += hbp;
      }
    }
  }

  return Success;
}


#define DEAL_WITH_STATISTICAL_ENERGY_TERMS
int RamaTableReadFromFile(RamaTable* pRama,char* ramafile){
  FILE* pFile=fopen(ramafile,"r");
  if(!pFile) return IOError;
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(fgets(buffer,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    if(buffer[0]==' '|| buffer[0]=='#') continue;
    int phipsi[2];
    double aae[20];
    sscanf(buffer,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
      &phipsi[0],&phipsi[1],&aae[0],&aae[1],&aae[2],&aae[3],&aae[4],&aae[5],&aae[6],&aae[7],&aae[8],&aae[9],
      &aae[10],&aae[11],&aae[12],&aae[13],&aae[14],&aae[15],&aae[16],&aae[17],&aae[18],&aae[19]);
    int phiindex=(phipsi[0]+180)/10;
    int psiindex=(phipsi[1]+180)/10;
    for(int j=0;j<20;j++){
      pRama->ramatable[phiindex][psiindex][j]=aae[j];
    }
  }
  fclose(pFile);
  return Success;
}

int AApropensityTableReadFromFile(AAppTable* pAAppTable,char* aappfile){
  FILE* pFile=fopen(aappfile,"r");
  if(!pFile) return IOError;
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  while(fgets(buffer,MAX_LENGTH_ONE_LINE_IN_FILE,pFile)){
    if(buffer[0]==' '|| buffer[0]=='#') continue;
    int phipsi[2];
    double aae[20];
    sscanf(buffer,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
      &phipsi[0],&phipsi[1],&aae[0],&aae[1],&aae[2],&aae[3],&aae[4],&aae[5],&aae[6],&aae[7],&aae[8],&aae[9],
      &aae[10],&aae[11],&aae[12],&aae[13],&aae[14],&aae[15],&aae[16],&aae[17],&aae[18],&aae[19]);
    int phiindex=(phipsi[0]+180)/10;
    int psiindex=(phipsi[1]+180)/10;
    for(int j=0;j<20;j++){
      pAAppTable->aapptable[phiindex][psiindex][j]=aae[j];
    }
  }
  fclose(pFile);
  return Success;
}

int AminoAcidPropensityAndRamachandranEnergy(Residue* pThis,AAppTable* pAAppTable,RamaTable* pRama){
  char aa1=ThreeLetterAAToOneLetterAA(ResidueGetName(pThis));
  int aaindex=AminoAcidGetIndex(aa1);
  int phiindex=((int)pThis->phipsi[0]+180)/10;
  int psiindex=((int)pThis->phipsi[1]+180)/10;
  pThis->aapp+=pAAppTable->aapptable[phiindex][psiindex][aaindex];
  pThis->rama+=pRama->ramatable[phiindex][psiindex][aaindex];
  return Success;
}


int RotamerPropensityAndRamachandranEnergy(Rotamer* pThis,Residue* pResidue,AAppTable* pAAppTable,RamaTable* pRama,double energyTerms[MAX_ENERGY_TERM]){
  char aa1=ThreeLetterAAToOneLetterAA(RotamerGetType(pThis));
  int aaindex=AminoAcidGetIndex(aa1);
  int phiindex=((int)pResidue->phipsi[0]+180)/10;
  int psiindex=((int)pResidue->phipsi[1]+180)/10;
  energyTerms[91]+=pAAppTable->aapptable[phiindex][psiindex][aaindex];
  energyTerms[92]+=pRama->ramatable[phiindex][psiindex][aaindex];
  return Success;
}


int AminoAcidDunbrackEnergy(Residue* pThis,BBdepRotamerLib* pBBdepRotLib){
  if(strcmp(ResidueGetName(pThis),"ALA")==0 || strcmp(ResidueGetName(pThis),"GLY")==0) return Success;
  if(pThis->isSCIntact){
    int binIdx=((int)(pThis->phipsi[0]+180)/10)*36+(int)(pThis->phipsi[1]+180)/10;
    RotLibPhiPsi* pRotLibPhiPsi=&pBBdepRotLib->rotlibphipsis[binIdx];
    int rotTypeIdx=-1;
    StringArrayFind(&pRotLibPhiPsi->rotTypes,ResidueGetName(pThis),&rotTypeIdx);
    DoubleArray* pTorsionsArray = pRotLibPhiPsi->torsions[rotTypeIdx];
    DoubleArray* pDeviationsArray = pRotLibPhiPsi->deviations[rotTypeIdx];
    int matchIdx=-1;
    double pMatch;
    double pMin=10;
    for(int i=0;i<IntArrayGet(&pRotLibPhiPsi->rotamerCounts,rotTypeIdx);i++){
      double p=DoubleArrayGet(&pRotLibPhiPsi->probability[rotTypeIdx],i);
      if(p<CUT_EXCL_LOW_ROT_PROB) break;
      if(p<pMin) pMin=p;
      DoubleArray* pTorsions=&pTorsionsArray[i];
      DoubleArray* pDeviations=&pDeviationsArray[i];

      BOOL match=TRUE;
      for(int j=0;j<DoubleArrayGetLength(&pThis->Xs);j++){
        //use a strict criteria
        //double min=DoubleArrayGet(pTorsions,j)-DegToRad(5.0);
        //double max=DoubleArrayGet(pTorsions,j)+DegToRad(5.0);
        double min=DoubleArrayGet(pTorsions,j)-DoubleArrayGet(pDeviations,j);
        double max=DoubleArrayGet(pTorsions,j)+DoubleArrayGet(pDeviations,j);
        double torsion=DoubleArrayGet(&pThis->Xs,j);
        double torsionm2pi=torsion-2*PI;
        double torsionp2pi=torsion+2*PI;
        double torsion2=torsion;
        if((strcmp(ResidueGetName(pThis),"PHE")==0 && j==1) || (strcmp(ResidueGetName(pThis),"TYR")==0 && j==1) ||
          (strcmp(ResidueGetName(pThis),"ASP")==0 && j==1)|| strcmp(ResidueGetName(pThis),"GLU")==0 && j==2){
            torsion2=torsion+PI;
            torsion2=torsion>0?torsion-PI:torsion2;
        }
        double torsion2m2pi=torsion2-2*PI;
        double torsion2p2pi=torsion2+2*PI;
        if(!((torsion<=max && torsion>=min) || (torsionm2pi<=max && torsionm2pi>=min) || (torsionp2pi<=max && torsionp2pi>=min) ||
          (torsion2<=max && torsion2>=min) || (torsion2m2pi<=max && torsion2m2pi>=min) || (torsion2p2pi<=max && torsion2p2pi>=min))){
            match=FALSE;
            break;
        }
      }
      if(match){
        matchIdx=i;
        pMatch=p;
        break;
      }
    }
    double pDelta=1e-7;
    if(matchIdx!=-1) ResidueSetDunbrack(pThis,-log(pMatch+pDelta));
    else ResidueSetDunbrack(pThis,-log(pMin+pDelta));
  }
  
  return Success;
}


int RotamerDunbrackEnergy(Rotamer* pThis,double energyTerms[MAX_ENERGY_TERM]){
  energyTerms[93]+=pThis->dunbrack;
  return Success;
}

int EnergyResidueAndResidueSameChain(Residue* pThis,Residue* pOther,double energyTerms[MAX_ENERGY_TERM]){
  int neighborCheck = 0;
  if(ResidueGetPosInChain(pThis)+1==ResidueGetPosInChain(pOther)) neighborCheck=12;
  else if(ResidueGetPosInChain(pThis)-1==ResidueGetPosInChain(pOther)) neighborCheck=21;

  for(int i=0;i<ResidueGetAtomCount(pThis);i++){
    Atom* pAtom1=ResidueGetAtom(pThis,i);
    if(pAtom1->isXyzValid==FALSE) continue;
    for(int j=0;j<ResidueGetAtomCount(pOther);j++){
      Atom* pAtom2=ResidueGetAtom(pOther,j);
      if(pAtom2->isXyzValid==FALSE) continue;
      if((neighborCheck==12 || neighborCheck==21) &&(pAtom1->isBBAtom && pAtom2->isBBAtom)){
        continue;
      }
      double distance=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
      if(distance>ENERGY_DISTANCE_CUTOFF) continue;
      int bondType=15;
      if(neighborCheck==12) bondType=ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1),AtomGetName(pAtom2),ResidueGetName(pOther));
      else if(neighborCheck==21) bondType=ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom2),AtomGetName(pAtom1),ResidueGetName(pThis));
      if(bondType==12 || bondType==13) continue;
      double att=0,rep=0,ele=0,desP=0,desH=0;
      VdwAttEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&att);
      VdwRepEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&rep);
      ElecEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1,pAtom2,distance,bondType,&desP,&desH);
      energyTerms[31]+=att;
      energyTerms[32]+=rep;
      energyTerms[33]+=ele;
      energyTerms[34]+=desP;
      energyTerms[35]+=desH;
      if(distance<HBOND_DISTANCE_CUTOFF_MAX){
        double hbtot=0,hbd=0,hbt=0,hbp=0;
        if(pAtom1->isHBatomH && pAtom2->isHBatomA){
          HBondEnergyAtomAndAtom(pAtom1,pAtom2,ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        else if(pAtom1->isHBatomA && pAtom2->isHBatomH){
          HBondEnergyAtomAndAtom(pAtom2,pAtom1,ResidueGetAtomByName(pOther,AtomGetHbDorB(pAtom2)),ResidueGetAtomByName(pThis,AtomGetHbDorB(pAtom1)),
            distance,bondType,&hbtot,&hbd,&hbt,&hbp);
        }
        if(pThis->posInChain-pOther->posInChain<=2 && pThis->posInChain-pOther->posInChain>=-2){
          hbd*=HBOND_LOCAL_REDUCE;
          hbt*=HBOND_LOCAL_REDUCE;
          hbp*=HBOND_LOCAL_REDUCE;
        }
        if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
          energyTerms[41]+=hbd;
          energyTerms[42]+=hbt;
          energyTerms[43]+=hbp;
        }
        else if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
          energyTerms[47]+=hbd;
          energyTerms[48]+=hbt;
          energyTerms[49]+=hbp;
        }
        else{
          energyTerms[44]+=hbd;
          energyTerms[45]+=hbt;
          energyTerms[46]+=hbp;
        }
      }
    }
  }

  //consider the SSBOND
  if(strcmp("CYS",ResidueGetName(pThis))==0 && strcmp("CYS",ResidueGetName(pOther))==0){
    Atom* pSG1=NULL,*pSG2=NULL,*pCB1=NULL,*pCB2=NULL,*pCA1=NULL,*pCA2=NULL;
    for(int i=0;i<ResidueGetAtomCount(pThis);i++){
      Atom* pAtom=ResidueGetAtom(pThis,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB1=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA1=pAtom;
    }
    for(int i=0;i<ResidueGetAtomCount(pOther);i++){
      Atom* pAtom=ResidueGetAtom(pOther,i);
      if(strcmp(AtomGetName(pAtom),"SG")==0) pSG2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CB")==0) pCB2=pAtom;
      else if(strcmp(AtomGetName(pAtom),"CA")==0) pCA2=pAtom;
    }
    if(pSG1!=NULL && pCB1!=NULL &&pSG2!=NULL && pCB2!=NULL && pCA1!=NULL && pCA2!=NULL){
      double dist=XYZDistance(&pSG1->xyz,&pSG2->xyz);
      if(dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN){
        double ssbond=0;
        SSbondEnergyAtomAndAtom(pSG1,pSG2,pCB1,pCB2,pCA1,pCA2,&ssbond);
        energyTerms[36]+=ssbond;
      }
    }
  }

  return Success;
}
