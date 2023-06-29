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

#include "EnergyFunction.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>


extern double CUT_EXCL_LOW_PROB_ROT;

double WEIGHTS[MAX_ENERGY_TERM];

#define DEAL_WITH_ENERGY_WEIGHTING
int EnergyTermInitialize(double* energyTerms)
{
  memset(energyTerms, 0, sizeof(double) * MAX_ENERGY_TERM);
  return Success;
}

int EnergyTermWeighting(double* energyTerms)
{
  for (int i = 1; i < MAX_ENERGY_TERM; i++)
  {
    energyTerms[i] *= WEIGHTS[i];
    energyTerms[0] += energyTerms[i];
  }
  return Success;
}


int EnergyTermShowMonomer(double* energyTerms)
{
  printf("reference_ALA         =            %12.6f\n", energyTerms[1]);
  printf("reference_CYS         =            %12.6f\n", energyTerms[2]);
  printf("reference_ASP         =            %12.6f\n", energyTerms[3]);
  printf("reference_GLU         =            %12.6f\n", energyTerms[4]);
  printf("reference_PHE         =            %12.6f\n", energyTerms[5]);
  printf("reference_GLY         =            %12.6f\n", energyTerms[6]);
  printf("reference_HIS         =            %12.6f\n", energyTerms[7]);
  printf("reference_ILE         =            %12.6f\n", energyTerms[8]);
  printf("reference_LYS         =            %12.6f\n", energyTerms[9]);
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


int EnergyTermShowComplex(double* energyTerms)
{
  printf("reference_ALA         =            %12.6f\n", energyTerms[1]);
  printf("reference_CYS         =            %12.6f\n", energyTerms[2]);
  printf("reference_ASP         =            %12.6f\n", energyTerms[3]);
  printf("reference_GLU         =            %12.6f\n", energyTerms[4]);
  printf("reference_PHE         =            %12.6f\n", energyTerms[5]);
  printf("reference_GLY         =            %12.6f\n", energyTerms[6]);
  printf("reference_HIS         =            %12.6f\n", energyTerms[7]);
  printf("reference_ILE         =            %12.6f\n", energyTerms[8]);
  printf("reference_LYS         =            %12.6f\n", energyTerms[9]);
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


int EnergyWeightRead(char* weightfile)
{
  for (int i = 0;i < MAX_ENERGY_TERM;i++)
  {
    WEIGHTS[i] = 1.0;
  }
  FILE* pFile = fopen(weightfile, "r");
  if (pFile != NULL)
  {
    char line[MAX_LEN_ONE_LINE_CONTENT + 1];
    while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pFile))
    {
      char term[MAX_LEN_ONE_LINE_CONTENT + 1];
      double val = 0.0;
      sscanf(line, "%s %lf", term, &val);
      if (!strcmp(term, "reference_ALA"))
      {
        WEIGHTS[1] = val;
      }
      else if (!strcmp(term, "reference_CYS"))
      {
        WEIGHTS[2] = val;
      }
      else if (!strcmp(term, "reference_ASP"))
      {
        WEIGHTS[3] = val;
      }
      else if (!strcmp(term, "reference_GLU"))
      {
        WEIGHTS[4] = val;
      }
      else if (!strcmp(term, "reference_PHE"))
      {
        WEIGHTS[5] = val;
      }
      else if (!strcmp(term, "reference_GLY"))
      {
        WEIGHTS[6] = val;
      }
      else if (!strcmp(term, "reference_HIS"))
      {
        WEIGHTS[7] = val;
      }
      else if (!strcmp(term, "reference_ILE"))
      {
        WEIGHTS[8] = val;
      }
      else if (!strcmp(term, "reference_LYS"))
      {
        WEIGHTS[9] = val;
      }
      else if (!strcmp(term, "reference_LEU"))
      {
        WEIGHTS[10] = val;
      }
      else if (!strcmp(term, "reference_MET"))
      {
        WEIGHTS[11] = val;
      }
      else if (!strcmp(term, "reference_ASN"))
      {
        WEIGHTS[12] = val;
      }
      else if (!strcmp(term, "reference_PRO"))
      {
        WEIGHTS[13] = val;
      }
      else if (!strcmp(term, "reference_GLN"))
      {
        WEIGHTS[14] = val;
      }
      else if (!strcmp(term, "reference_ARG"))
      {
        WEIGHTS[15] = val;
      }
      else if (!strcmp(term, "reference_SER"))
      {
        WEIGHTS[16] = val;
      }
      else if (!strcmp(term, "reference_THR"))
      {
        WEIGHTS[17] = val;
      }
      else if (!strcmp(term, "reference_VAL"))
      {
        WEIGHTS[18] = val;
      }
      else if (!strcmp(term, "reference_TRP"))
      {
        WEIGHTS[19] = val;
      }
      else if (!strcmp(term, "reference_TYR"))
      {
        WEIGHTS[20] = val;
      }

      else if (!strcmp(term, "intraR_vdwatt"))
      {
        WEIGHTS[21] = val;
      }
      else if (!strcmp(term, "intraR_vdwrep"))
      {
        WEIGHTS[22] = val;
      }
      else if (!strcmp(term, "intraR_electr"))
      {
        WEIGHTS[23] = val;
      }
      else if (!strcmp(term, "intraR_deslvP"))
      {
        WEIGHTS[24] = val;
      }
      else if (!strcmp(term, "intraR_deslvH"))
      {
        WEIGHTS[25] = val;
      }
      else if (!strcmp(term, "intraR_hbscbb_dis"))
      {
        WEIGHTS[26] = val;
      }
      else if (!strcmp(term, "intraR_hbscbb_the"))
      {
        WEIGHTS[27] = val;
      }
      else if (!strcmp(term, "intraR_hbscbb_phi"))
      {
        WEIGHTS[28] = val;
      }
      else if (!strcmp(term, "aapropensity"))
      {
        WEIGHTS[91] = val;
      }
      else if (!strcmp(term, "ramachandran"))
      {
        WEIGHTS[92] = val;
      }
      else if (!strcmp(term, "dunbrack"))
      {
        WEIGHTS[93] = val;
      }

      else if (!strcmp(term, "interS_vdwatt")) WEIGHTS[31] = val;
      else if (!strcmp(term, "interS_vdwrep")) WEIGHTS[32] = val;
      else if (!strcmp(term, "interS_electr")) WEIGHTS[33] = val;
      else if (!strcmp(term, "interS_deslvP")) WEIGHTS[34] = val;
      else if (!strcmp(term, "interS_deslvH")) WEIGHTS[35] = val;
      else if (!strcmp(term, "interS_ssbond")) WEIGHTS[36] = val;
      else if (!strcmp(term, "interS_hbbbbb_dis")) WEIGHTS[41] = val;
      else if (!strcmp(term, "interS_hbbbbb_the")) WEIGHTS[42] = val;
      else if (!strcmp(term, "interS_hbbbbb_phi")) WEIGHTS[43] = val;
      else if (!strcmp(term, "interS_hbscbb_dis")) WEIGHTS[44] = val;
      else if (!strcmp(term, "interS_hbscbb_the")) WEIGHTS[45] = val;
      else if (!strcmp(term, "interS_hbscbb_phi")) WEIGHTS[46] = val;
      else if (!strcmp(term, "interS_hbscsc_dis")) WEIGHTS[47] = val;
      else if (!strcmp(term, "interS_hbscsc_the")) WEIGHTS[48] = val;
      else if (!strcmp(term, "interS_hbscsc_phi")) WEIGHTS[49] = val;

      else if (!strcmp(term, "interD_vdwatt")) WEIGHTS[51] = val;
      else if (!strcmp(term, "interD_vdwrep")) WEIGHTS[52] = val;
      else if (!strcmp(term, "interD_electr")) WEIGHTS[53] = val;
      else if (!strcmp(term, "interD_deslvP")) WEIGHTS[54] = val;
      else if (!strcmp(term, "interD_deslvH")) WEIGHTS[55] = val;
      else if (!strcmp(term, "interD_ssbond")) WEIGHTS[56] = val;
      else if (!strcmp(term, "interD_hbbbbb_dis")) WEIGHTS[61] = val;
      else if (!strcmp(term, "interD_hbbbbb_the")) WEIGHTS[62] = val;
      else if (!strcmp(term, "interD_hbbbbb_phi")) WEIGHTS[63] = val;
      else if (!strcmp(term, "interD_hbscbb_dis")) WEIGHTS[64] = val;
      else if (!strcmp(term, "interD_hbscbb_the")) WEIGHTS[65] = val;
      else if (!strcmp(term, "interD_hbscbb_phi")) WEIGHTS[66] = val;
      else if (!strcmp(term, "interD_hbscsc_dis")) WEIGHTS[67] = val;
      else if (!strcmp(term, "interD_hbscsc_the")) WEIGHTS[68] = val;
      else if (!strcmp(term, "interD_hbscsc_phi")) WEIGHTS[69] = val;

      else if (!strcmp(term, "ligand_vdwatt")) WEIGHTS[71] = val;
      else if (!strcmp(term, "ligand_vdwrep")) WEIGHTS[72] = val;
      else if (!strcmp(term, "ligand_electr")) WEIGHTS[73] = val;
      else if (!strcmp(term, "ligand_deslvP")) WEIGHTS[74] = val;
      else if (!strcmp(term, "ligand_deslvH")) WEIGHTS[75] = val;
      else if (!strcmp(term, "ligand_hbscbb_dis")) WEIGHTS[84] = val;
      else if (!strcmp(term, "ligand_hbscbb_the")) WEIGHTS[85] = val;
      else if (!strcmp(term, "ligand_hbscbb_phi")) WEIGHTS[86] = val;
      else if (!strcmp(term, "ligand_hbscsc_dis")) WEIGHTS[87] = val;
      else if (!strcmp(term, "ligand_hbscsc_the")) WEIGHTS[88] = val;
      else if (!strcmp(term, "ligand_hbscsc_phi")) WEIGHTS[89] = val;
    }
    fclose(pFile);
  }
  else
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];
    sprintf(errMsg, "in file %s line %d, cannot open weight file %s", __FILE__, __LINE__, weightfile);
    TraceError(errMsg, IOError);
    return IOError;
  }

  return Success;
}

int EnergyWeightWrite(char* file)
{
  FILE* pFile = fopen(file, "w");
  fprintf(pFile, "reference_ALA      %8.3f\n", WEIGHTS[1]);
  fprintf(pFile, "reference_CYS      %8.3f\n", WEIGHTS[2]);
  fprintf(pFile, "reference_ASP      %8.3f\n", WEIGHTS[3]);
  fprintf(pFile, "reference_GLU      %8.3f\n", WEIGHTS[4]);
  fprintf(pFile, "reference_PHE      %8.3f\n", WEIGHTS[5]);
  fprintf(pFile, "reference_GLY      %8.3f\n", WEIGHTS[6]);
  fprintf(pFile, "reference_HIS      %8.3f\n", WEIGHTS[7]);
  fprintf(pFile, "reference_ILE      %8.3f\n", WEIGHTS[8]);
  fprintf(pFile, "reference_LYS      %8.3f\n", WEIGHTS[9]);
  fprintf(pFile, "reference_LEU      %8.3f\n", WEIGHTS[10]);
  fprintf(pFile, "reference_MET      %8.3f\n", WEIGHTS[11]);
  fprintf(pFile, "reference_ASN      %8.3f\n", WEIGHTS[12]);
  fprintf(pFile, "reference_PRO      %8.3f\n", WEIGHTS[13]);
  fprintf(pFile, "reference_GLN      %8.3f\n", WEIGHTS[14]);
  fprintf(pFile, "reference_ARG      %8.3f\n", WEIGHTS[15]);
  fprintf(pFile, "reference_SER      %8.3f\n", WEIGHTS[16]);
  fprintf(pFile, "reference_THR      %8.3f\n", WEIGHTS[17]);
  fprintf(pFile, "reference_VAL      %8.3f\n", WEIGHTS[18]);
  fprintf(pFile, "reference_TRP      %8.3f\n", WEIGHTS[19]);
  fprintf(pFile, "reference_TYR      %8.3f\n", WEIGHTS[20]);

  fprintf(pFile, "intraR_vdwatt      %8.3f\n", WEIGHTS[21]);
  fprintf(pFile, "intraR_vdwrep      %8.3f\n", WEIGHTS[22]);
  fprintf(pFile, "intraR_electr      %8.3f\n", WEIGHTS[23]);
  fprintf(pFile, "intraR_deslvP      %8.3f\n", WEIGHTS[24]);
  fprintf(pFile, "intraR_deslvH      %8.3f\n", WEIGHTS[25]);
  fprintf(pFile, "intraR_hbscbb_dis  %8.3f\n", WEIGHTS[26]);
  fprintf(pFile, "intraR_hbscbb_the  %8.3f\n", WEIGHTS[27]);
  fprintf(pFile, "intraR_hbscbb_phi  %8.3f\n", WEIGHTS[28]);
  fprintf(pFile, "aapropensity       %8.3f\n", WEIGHTS[91]);
  fprintf(pFile, "ramachandran       %8.3f\n", WEIGHTS[92]);
  fprintf(pFile, "dunbrack           %8.3f\n", WEIGHTS[93]);

  fprintf(pFile, "interS_vdwatt      %8.3f\n", WEIGHTS[31]);
  fprintf(pFile, "interS_vdwrep      %8.3f\n", WEIGHTS[32]);
  fprintf(pFile, "interS_electr      %8.3f\n", WEIGHTS[33]);
  fprintf(pFile, "interS_deslvP      %8.3f\n", WEIGHTS[34]);
  fprintf(pFile, "interS_deslvH      %8.3f\n", WEIGHTS[35]);
  fprintf(pFile, "interS_ssbond      %8.3f\n", WEIGHTS[36]);
  fprintf(pFile, "interS_hbbbbb_dis  %8.3f\n", WEIGHTS[41]);
  fprintf(pFile, "interS_hbbbbb_the  %8.3f\n", WEIGHTS[42]);
  fprintf(pFile, "interS_hbbbbb_phi  %8.3f\n", WEIGHTS[43]);
  fprintf(pFile, "interS_hbscbb_dis  %8.3f\n", WEIGHTS[44]);
  fprintf(pFile, "interS_hbscbb_the  %8.3f\n", WEIGHTS[45]);
  fprintf(pFile, "interS_hbscbb_phi  %8.3f\n", WEIGHTS[46]);
  fprintf(pFile, "interS_hbscsc_dis  %8.3f\n", WEIGHTS[47]);
  fprintf(pFile, "interS_hbscsc_the  %8.3f\n", WEIGHTS[48]);
  fprintf(pFile, "interS_hbscsc_phi  %8.3f\n", WEIGHTS[49]);

  fprintf(pFile, "interD_vdwatt      %8.3f\n", WEIGHTS[51]);
  fprintf(pFile, "interD_vdwrep      %8.3f\n", WEIGHTS[52]);
  fprintf(pFile, "interD_electr      %8.3f\n", WEIGHTS[53]);
  fprintf(pFile, "interD_deslvP      %8.3f\n", WEIGHTS[54]);
  fprintf(pFile, "interD_deslvH      %8.3f\n", WEIGHTS[55]);
  fprintf(pFile, "interD_ssbond      %8.3f\n", WEIGHTS[56]);
  fprintf(pFile, "interD_hbbbbb_dis  %8.3f\n", WEIGHTS[61]);
  fprintf(pFile, "interD_hbbbbb_the  %8.3f\n", WEIGHTS[62]);
  fprintf(pFile, "intraD_hbbbbb_phi  %8.3f\n", WEIGHTS[63]);
  fprintf(pFile, "interD_hbscbb_dis  %8.3f\n", WEIGHTS[64]);
  fprintf(pFile, "interD_hbscbb_the  %8.3f\n", WEIGHTS[65]);
  fprintf(pFile, "interD_hbscbb_phi  %8.3f\n", WEIGHTS[66]);
  fprintf(pFile, "interD_hbscsc_dis  %8.3f\n", WEIGHTS[67]);
  fprintf(pFile, "interD_hbscsc_the  %8.3f\n", WEIGHTS[68]);
  fprintf(pFile, "interD_hbscsc_phi  %8.3f\n", WEIGHTS[69]);

  fprintf(pFile, "ligand_vdwatt      %8.3f\n", WEIGHTS[71]);
  fprintf(pFile, "ligand_vdwrep      %8.3f\n", WEIGHTS[72]);
  fprintf(pFile, "ligand_electr      %8.3f\n", WEIGHTS[73]);
  fprintf(pFile, "ligand_deslvP      %8.3f\n", WEIGHTS[74]);
  fprintf(pFile, "ligand_deslvH      %8.3f\n", WEIGHTS[75]);
  fprintf(pFile, "ligand_hbscbb_dis  %8.3f\n", WEIGHTS[84]);
  fprintf(pFile, "ligand_hbscbb_the  %8.3f\n", WEIGHTS[85]);
  fprintf(pFile, "ligand_hbscbb_phi  %8.3f\n", WEIGHTS[86]);
  fprintf(pFile, "ligand_hbscsc_dis  %8.3f\n", WEIGHTS[87]);
  fprintf(pFile, "ligand_hbscsc_the  %8.3f\n", WEIGHTS[88]);
  fprintf(pFile, "ligand_hbscsc_phi  %8.3f\n", WEIGHTS[89]);
  fclose(pFile);
  return Success;
}


#define DEAL_WITH_BOND_CONNECTION
//these functions are used to check the 12, 13, 14 and 15 bond connectivity
BOOL ResidueIntraBond12Check(char* atom1, char* atom2, BondSet* pBondSet)
{
  for (int i = 0; i < pBondSet->count; i++)
  {
    Bond* pBond = pBondSet->bonds + i;
    if ((strcmp(atom1, pBond->atomFromName) == 0 && strcmp(atom2, pBond->atomToName) == 0) ||
      (strcmp(atom2, pBond->atomFromName) == 0 && strcmp(atom1, pBond->atomToName) == 0))
    {
      return TRUE;
    }
  }
  return FALSE;
}

BOOL ResidueIntraBond13Check(char* atom1, char* atom2, BondSet* pBondSet)
{
  for (int i = 0; i < pBondSet->count; i++)
  {
    Bond* pBond = pBondSet->bonds + i;
    if (strcmp(atom1, pBond->atomFromName) == 0)
    {
      if (ResidueIntraBond12Check(pBond->atomToName, atom2, pBondSet))
      {
        return TRUE;
      }
    }
    else if (strcmp(atom1, pBond->atomToName) == 0)
    {
      if (ResidueIntraBond12Check(pBond->atomFromName, atom2, pBondSet))
      {
        return TRUE;
      }
    }
  }
  return FALSE;
}

BOOL ResidueIntraBond14Check(char* atom1, char* atom2, BondSet* pBondSet)
{
  for (int i = 0; i < pBondSet->count; i++)
  {
    Bond* pBond = pBondSet->bonds + i;
    if (strcmp(atom1, pBond->atomFromName) == 0)
    {
      if (ResidueIntraBond13Check(pBond->atomToName, atom2, pBondSet))
      {
        return TRUE;
      }
    }
    else if (strcmp(atom1, pBond->atomToName) == 0)
    {
      if (ResidueIntraBond13Check(pBond->atomFromName, atom2, pBondSet))
      {
        return TRUE;
      }
    }
  }
  return FALSE;
}

int ResidueIntraBondConnectionCheck(char* atom1, char* atom2, BondSet* pBondSet)
{
  if (ResidueIntraBond12Check(atom1, atom2, pBondSet))
  {
    return 12;
  }
  else if (ResidueIntraBond13Check(atom1, atom2, pBondSet))
  {
    return 13;
  }
  else if (ResidueIntraBond14Check(atom1, atom2, pBondSet))
  {
    return 14;
  }
  else
  {
    return 15;
  }
}

int ResidueAndNextResidueInterBondConnectionCheck_charmm22(char* atomOnPreResi, char* atomOnNextResi, Residue* pPreResi, Residue* pNextResi)
{
  if (strcmp(atomOnPreResi, "C") == 0)
  {
    if (strcmp(atomOnNextResi, "N") == 0)
    {
      return 12;
    }
    else if (strcmp(atomOnNextResi, "CA") == 0
      || strcmp(atomOnNextResi, "HN") == 0
      || (strcmp(atomOnNextResi, "CD") == 0 && strcmp(pNextResi->name, "PRO") == 0))
    {
      return 13;
    }
    else if (strcmp(atomOnNextResi, "CB") == 0
      || strcmp(atomOnNextResi, "C") == 0
      || strcmp(atomOnNextResi, "HA") == 0
      || strcmp(atomOnNextResi, "HA1") == 0
      || strcmp(atomOnNextResi, "HA2") == 0
      || (strcmp(atomOnNextResi, "HD1") == 0 && strcmp(pNextResi->name, "PRO") == 0)
      || (strcmp(atomOnNextResi, "HD2") == 0 && strcmp(pNextResi->name, "PRO") == 0)
      || (strcmp(atomOnNextResi, "CG") == 0 && strcmp(pNextResi->name, "PRO") == 0))
    {
      return 14;
    }
  }
  else if (strcmp(atomOnPreResi, "O") == 0
    || strcmp(atomOnPreResi, "CA") == 0)
  {
    if (strcmp(atomOnNextResi, "N") == 0)
    {
      return 13;
    }
    else if (strcmp(atomOnNextResi, "CA") == 0
      || strcmp(atomOnNextResi, "HN") == 0
      || (strcmp(atomOnNextResi, "CD") == 0 && strcmp(pNextResi->name, "PRO") == 0))
    {
      return 14;
    }
  }
  else if (strcmp(atomOnPreResi, "CB") == 0
    || strcmp(atomOnPreResi, "HA") == 0
    || strcmp(atomOnPreResi, "HA1") == 0
    || strcmp(atomOnPreResi, "HA2") == 0
    || strcmp(atomOnPreResi, "N") == 0)
  {
    if (strcmp(atomOnNextResi, "N") == 0)
    {
      return 14;
    }
  }
  return 15;
}


int ResidueAndNextResidueInterBondConnectionCheck_charmm19(char* atomOnPreResi, char* atomOnNextResi, char* nextResiName)
{
  if (strcmp(atomOnPreResi, "C") == 0)
  {
    if (strcmp(atomOnNextResi, "N") == 0)
    {
      return 12;
    }
    else if (strcmp(atomOnNextResi, "CA") == 0
      || strcmp(atomOnNextResi, "H") == 0
      || (strcmp(atomOnNextResi, "CD") == 0 && strcmp(nextResiName, "PRO") == 0))
    {
      return 13;
    }
    else if (strcmp(atomOnNextResi, "CB") == 0
      || strcmp(atomOnNextResi, "C") == 0
      || (strcmp(atomOnNextResi, "CG") == 0 && strcmp(nextResiName, "PRO") == 0))
    {
      return 14;
    }
  }
  else if (strcmp(atomOnPreResi, "O") == 0
    || strcmp(atomOnPreResi, "CA") == 0)
  {
    if (strcmp(atomOnNextResi, "N") == 0)
    {
      return 13;
    }
    else if (strcmp(atomOnNextResi, "CA") == 0
      || strcmp(atomOnNextResi, "H") == 0
      || (strcmp(atomOnNextResi, "CD") == 0 && strcmp(nextResiName, "PRO") == 0))
    {
      return 14;
    }
  }
  else if (strcmp(atomOnPreResi, "CB") == 0
    || strcmp(atomOnPreResi, "N") == 0)
  {
    if (strcmp(atomOnNextResi, "N") == 0)
    {
      return 14;
    }
  }

  return 15;
}


#define DEAL_WITH_BASIC_ENERGY_TERMS
int VdwAttEnergyAtomAndAtom(Atom* pAtom1, Atom* pAtom2, double distance, int bondType, double* vdwAtt)
{
  if (distance >= ENERGY_DISTANCE_CUTOFF)
  {
    return Success;
  }
  if (bondType == 12 || bondType == 13)
  {
    return Success;
  }
  //if (AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2))
  //{ // do not consider H-H vdw attractive energy
  //  return Success;
  //}
  double rsum = RADIUS_SCALE_FOR_VDW * (pAtom1->vdw_radius + pAtom2->vdw_radius);
  double ratio = distance / rsum;
  double energy = 0.0;
  double scale = 0.0;
  if (bondType == 14)
  {
    scale = ENERGY_SCALE_FACTOR_BOND_14;
  }
  else if (bondType == 15)
  {
    scale = ENERGY_SCALE_FACTOR_BOND_15;
  }

  if (ratio < 0.8909)
  {
    energy = 0.0;
  }
  else if (distance <= 5.0)
  {
    double epsilon = sqrt(pAtom1->vdw_epsilon * pAtom2->vdw_epsilon);
    double B6 = pow(1 / ratio, 6.0);
    double A12 = B6 * B6;
    energy = epsilon * (A12 - 2.0 * B6);
  }
  else if (distance > 5.0 && distance < ENERGY_DISTANCE_CUTOFF)
  {
    double epsilon = sqrt(pAtom1->vdw_epsilon * pAtom2->vdw_epsilon);
    double B6 = pow((double)rsum / 5.0, 6.0);
    double A12 = B6 * B6;
    double M = epsilon * (A12 - 2.0 * B6);
    double N = 2.4 * epsilon * (B6 - A12);
    double a = 2 * M + N;
    double b = -33 * M - 17 * N;
    double c = 180 * M + 96 * N;
    double d = -324 * M - 180 * N;
    energy = a * distance * distance * distance + b * distance * distance + c * distance + d;
  }
  energy *= scale;
  *vdwAtt = energy;
  if (ENERGY_DEBUG_MODE_VDW_ATT)
  {
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %f, ratio: %f, vdwAtt: %f\n",
      AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1),
      AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
      bondType, distance, ratio, energy);
  }
  return Success;
}


int VdwRepEnergyAtomAndAtom(Atom* pAtom1, Atom* pAtom2, double distance, int bondType, double* vdwRep)
{
  if (bondType == 12 || bondType == 13)
  {
    return Success;
  }
  double r1 = 0, r2 = 0;
  if (AtomIsHydrogen(pAtom1) == TRUE && (pAtom1->polarity == Type_AtomPolarity_P || pAtom1->polarity == Type_AtomPolarity_C))
  {
    r1 = 0.2 * pAtom1->vdw_radius;
  }
  else
  {
    r1 = RADIUS_SCALE_FOR_VDW * pAtom1->vdw_radius;
  }
  if (AtomIsHydrogen(pAtom2) == TRUE && (pAtom2->polarity == Type_AtomPolarity_P || pAtom2->polarity == Type_AtomPolarity_C))
  {
    r2 = 0.2 * pAtom2->vdw_radius;
  }
  else
  {
    r2 = RADIUS_SCALE_FOR_VDW * pAtom2->vdw_radius;
  }
  // double rsum = RADIUS_SCALE_FOR_VDW * (pAtom1->vdw_radius + pAtom2->vdw_radius);
  double rsum = r1 + r2;
  double ratio = distance / rsum;
  double epsilon = sqrt(pAtom1->vdw_epsilon * pAtom2->vdw_epsilon);
  double RATIO_CUTOFF = 0.70; // can be adjusted
  double energy = 0.0;
  double scale = 0.0;
  if (bondType == 14)
  {
    scale = ENERGY_SCALE_FACTOR_BOND_14;
  }
  else if (bondType == 15)
  {
    scale = ENERGY_SCALE_FACTOR_BOND_15;
  }

  if (ratio > 0.8909)
  {
    energy = 0.0;
  }
  else if (ratio >= RATIO_CUTOFF)
  {
    double B6 = pow(1 / ratio, 6.0);
    double A12 = B6 * B6;
    energy = epsilon * (A12 - 2.0 * B6);
  }
  else
  {
    double B6_0 = pow(1 / RATIO_CUTOFF, 6.0);
    double a = epsilon * (B6_0 * B6_0 - 2.0 * B6_0);
    double b = epsilon * 12.0 * (B6_0 / RATIO_CUTOFF - B6_0 * B6_0 / RATIO_CUTOFF);
    double y0 = a * epsilon;
    double k = b * epsilon;
    energy = k * (ratio - RATIO_CUTOFF) + y0;
  }
  energy *= scale;
  //set a cutoff for maximum clash
  //double MAX_CLASH=5.0*epsilon;
  //if(energy>MAX_CLASH)
  //{
  //  energy = MAX_CLASH;
  //}
  *vdwRep = energy;

  if (ENERGY_DEBUG_MODE_VDW_REP)
  {
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %f, ratio: %f, vdwRep: %f\n",
      AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1),
      AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
      bondType, distance, ratio, energy);
  }
  return Success;
}


int HBondEnergyAtomAndAtom(Atom* atomH, Atom* atomA, Atom* atomD, Atom* atomB, double distHA, int bondType, double* etot, double* edist, double* eDHA, double* eHAB)
{
  if (bondType == 12 || bondType == 13)
  {
    return Success;
  }

  if (atomB == NULL && strcmp(AtomGetHbHorA(atomA), "D") == 0 && !(atomD == NULL && strcmp(AtomGetHbHorA(atomH), "D") == 0))
  { // atomA is in fact the water oxygen atom
    if (distHA > HB_HA_DIST_CUTOFF_MAX || distHA < HB_HA_DIST_CUTOFF_MIN)
    {
      return Success;
    }
    XYZ xyzDH = XYZDifference(&atomD->xyz, &atomH->xyz);
    XYZ xyzHA = XYZDifference(&atomH->xyz, &atomA->xyz);
    double angleDHA = PI - XYZAngle(&xyzDH, &xyzHA);
    if (RadToDeg(angleDHA) < ANGLE_DHA_CUTOFF_MIN)
    {
      return Success;
    }

    double energyDist = 0.0;
    if (distHA < HB_HA_OPT_DIST)
    {
      energyDist = -1.0 * HB_HA_WELL_DEPTH * cos((distHA - HB_HA_OPT_DIST) * PI);
    }
    else
    {
      energyDist = -0.5 * cos(PI / (HB_HA_DIST_CUTOFF_MAX - HB_HA_OPT_DIST) * (distHA - HB_HA_OPT_DIST)) - 0.5;
    }
    if (energyDist > 0.0)
    {
      energyDist = 0.0;
    }

    double energyDHA = -1.0 * cos(angleDHA) * cos(angleDHA) * cos(angleDHA) * cos(angleDHA);

    double energy = 0.0;
    if (RadToDeg(angleDHA) >= ANGLE_DHA_CUTOFF_MIN && distHA <= HB_HA_DIST_CUTOFF_MAX && distHA >= HB_HA_DIST_CUTOFF_MIN)
    {
      energy = energyDist + energyDHA;
      if (energy > 0.0)
      {
        energy = 0.0;
      }
      *etot = energy;
      *edist = energyDist;
      *eDHA = energyDHA;
      *eHAB = 0;
    }

    if (ENERGY_DEBUG_MODE_HBOND)
    {
      printf("AtomH: %1s %4d %4s, AtomA: %1s %4d %4s, dist: %5.2f, theta: %5.1f, edist: %5.2f, ethe: %5.2f, \n",
        AtomGetChainName(atomH), AtomGetPosInChain(atomH), AtomGetName(atomH),
        AtomGetChainName(atomA), AtomGetPosInChain(atomA), AtomGetName(atomA),
        distHA, RadToDeg(angleDHA), energyDist, energyDHA);
    }
  }
  else if (atomD == NULL && strcmp(AtomGetHbHorA(atomH), "D") == 0 && !(atomB == NULL && strcmp(AtomGetHbHorA(atomA), "D") == 0))
  { // atomH is in fact the water oxygen atom
    if (distHA > HB_DA_DIST_CUTOFF_MAX || distHA < HB_DA_DIST_CUTOFF_MIN)
    {
      return Success;
    }
    XYZ xyzHA = XYZDifference(&atomH->xyz, &atomA->xyz);
    XYZ xyzAAB = XYZDifference(&atomA->xyz, &atomB->xyz);
    double angleHAB = PI - XYZAngle(&xyzHA, &xyzAAB);
    if (RadToDeg(angleHAB) < ANGLE_DAB_CUTOFF_MIN)
    {
      return Success;
    }

    double energyDist = 0.0;
    if (distHA < HB_DA_OPT_DIST)
    {
      energyDist = -1.0 * HB_DA_WELL_DEPTH * cos((distHA - HB_DA_OPT_DIST) * PI);
    }
    else
    {
      energyDist = -0.5 * cos(PI / (HB_DA_DIST_CUTOFF_MAX - HB_DA_OPT_DIST) * (distHA - HB_DA_OPT_DIST)) - 0.5;
    }
    if (energyDist > 0.0)
    {
      energyDist = 0.0;
    }

    double energyHAB = 0.0;
    if (atomH->isBBAtom == TRUE && atomA->isBBAtom == TRUE)
    {
      energyHAB = -1.0 * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150));
    }
    else
    {
      if (atomA->hybridType == Type_AtomHybridType_SP3)
      {
        energyHAB = -1.0 * cos(angleHAB - DegToRad(135)) * cos(angleHAB - DegToRad(135)) * cos(angleHAB - DegToRad(135)) * cos(angleHAB - DegToRad(135));
      }
      else if (atomA->hybridType == Type_AtomHybridType_SP2)
      {
        energyHAB = -1.0 * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150));
      }
    }

    double energy = 0.0;
    if (RadToDeg(angleHAB) >= ANGLE_DAB_CUTOFF_MIN && distHA <= HB_DA_DIST_CUTOFF_MAX && distHA >= HB_DA_DIST_CUTOFF_MIN)
    {
      energy = energyDist + energyHAB;
      if (energy > 0.0)
      {
        energy = 0.0;
      }
      *etot = energy;
      *edist = energyDist;
      *eDHA = 0;
      *eHAB = energyHAB;
    }

    if (ENERGY_DEBUG_MODE_HBOND)
    {
      printf("AtomH: %1s %4d %4s, AtomA: %1s %4d %4s, dist: %5.2f, phi: %5.1f, edist: %5.2f, ephi: %5.2f\n",
        AtomGetChainName(atomH), AtomGetPosInChain(atomH), AtomGetName(atomH),
        AtomGetChainName(atomA), AtomGetPosInChain(atomA), AtomGetName(atomA),
        distHA, RadToDeg(angleHAB), energyDist, energyHAB);
    }
  }
  else if ((atomB == NULL && strcmp(AtomGetHbHorA(atomA), "D") == 0) && (atomD == NULL && strcmp(AtomGetHbHorA(atomH), "D") == 0))
  { // both atomA and atomH are water oxygen atoms
    if (distHA > HB_DA_DIST_CUTOFF_MAX)
    {
      return Success;
    }
    double energyDist = 0.0;
    if (distHA < HB_DA_OPT_DIST)
    {
      energyDist = -1.0 * HB_DA_WELL_DEPTH * cos((distHA - HB_DA_OPT_DIST) * PI);
    }
    else
    {
      energyDist = -0.5 * cos(PI / (HB_DA_DIST_CUTOFF_MAX - HB_DA_OPT_DIST) * (distHA - HB_DA_OPT_DIST)) - 0.5;
    }
    if (energyDist > 0.0)
    {
      energyDist = 0.0;
    }
    *etot = energyDist;
    *edist = energyDist;
    *eDHA = 0;
    *eHAB = 0;

    if (ENERGY_DEBUG_MODE_HBOND)
    {
      printf("AtomH: %1s %4d %4s, AtomA: %1s %4d %4s, dist: %5.2f, edist: %5.2f\n",
        AtomGetChainName(atomH), AtomGetPosInChain(atomH), AtomGetName(atomH),
        AtomGetChainName(atomA), AtomGetPosInChain(atomA), AtomGetName(atomA),
        distHA, energyDist);
    }
  }
  else
  { // regular hydrogen bond
    if (distHA > HB_HA_DIST_CUTOFF_MAX || distHA < HB_HA_DIST_CUTOFF_MIN)
    {
      return Success;
    }
    XYZ xyzDH = XYZDifference(&atomD->xyz, &atomH->xyz);
    XYZ xyzHA = XYZDifference(&atomH->xyz, &atomA->xyz);
    double angleDHA = PI - XYZAngle(&xyzDH, &xyzHA);
    if (RadToDeg(angleDHA) < ANGLE_DHA_CUTOFF_MIN)
    {
      return Success;
    }
    XYZ xyzAAB = XYZDifference(&atomA->xyz, &atomB->xyz);
    double angleHAB = PI - XYZAngle(&xyzHA, &xyzAAB);
    if (RadToDeg(angleHAB) < ANGLE_HAB_CUTOFF_MIN)
    {
      return Success;
    }

    // calculate HA distance-based energy term
    double energyDist = 0.0;
    if (distHA < HB_HA_OPT_DIST)
    {
      energyDist = -1.0 * HB_HA_WELL_DEPTH * cos((distHA - HB_HA_OPT_DIST) * PI);
    }
    else
    {
      energyDist = -0.5 * cos(PI / (HB_HA_DIST_CUTOFF_MAX - HB_HA_OPT_DIST) * (distHA - HB_HA_OPT_DIST)) - 0.5;
    }
    if (energyDist > 0.0)
    {
      energyDist = 0.0;
    }

    // calculate D-H-A angle-based energy term
    double energyDHA = -1.0 * cos(angleDHA) * cos(angleDHA) * cos(angleDHA) * cos(angleDHA);
    double energyHAB = 0.0;
    if (atomH->isBBAtom == TRUE && atomA->isBBAtom == TRUE)
    {
      energyHAB = -1.0 * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150));
    }
    else
    {
      if (atomA->hybridType == Type_AtomHybridType_SP3)
      {
        energyHAB = -1.0 * cos(angleHAB - DegToRad(135)) * cos(angleHAB - DegToRad(135)) * cos(angleHAB - DegToRad(135)) * cos(angleHAB - DegToRad(135));
      }
      else if (atomA->hybridType == Type_AtomHybridType_SP2)
      {
        energyHAB = -1.0 * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150)) * cos(angleHAB - DegToRad(150));
      }
    }

    // calcualte H-A-AB angle-based energy term
    // original function is error, we should add angle restriction to calculate energy
    double energy = 0.0;
    if (RadToDeg(angleDHA) >= ANGLE_DHA_CUTOFF_MIN && RadToDeg(angleHAB) >= ANGLE_HAB_CUTOFF_MIN && distHA <= HB_HA_DIST_CUTOFF_MAX && distHA >= HB_HA_DIST_CUTOFF_MIN)
    {
      energy = energyDist + energyDHA + energyHAB;
      if (energy > 0.0)
      {
        energy = 0.0;
      }
      *etot = energy;
      *edist = energyDist;
      *eDHA = energyDHA;
      *eHAB = energyHAB;
    }

    if (ENERGY_DEBUG_MODE_HBOND)
    {
      printf("AtomH: %1s %4d %4s, AtomA: %1s %4d %4s, dist: %5.2f, theta: %5.1f, phi: %5.1f, edist: %5.2f, ethe: %5.2f, ephi: %5.2f\n",
        AtomGetChainName(atomH), AtomGetPosInChain(atomH), AtomGetName(atomH),
        AtomGetChainName(atomA), AtomGetPosInChain(atomA), AtomGetName(atomA),
        distHA, RadToDeg(angleDHA), RadToDeg(angleHAB), energyDist, energyDHA, energyHAB);
    }
  }

  return Success;
}





int ElecEnergyAtomAndAtom(Atom* pAtom1, Atom* pAtom2, double dist12, int bondType, double* elec)
{
  if (bondType == 12 || bondType == 13)
  {
    return Success;
  }
  if (dist12 > ELEC_DISTANCE_CUTOFF)
  {
    return Success;
  }
  if (fabs(pAtom1->charge) < 1e-2 || fabs(pAtom2->charge) < 1e-2)
  {
    return Success;
  }
  else if (dist12 < 0.8 * (pAtom1->vdw_radius + pAtom2->vdw_radius))
  {
    dist12 = 0.8 * (pAtom1->vdw_radius + pAtom2->vdw_radius);
  }

  double energy = COULOMB_CONSTANT * pAtom1->charge * pAtom2->charge / dist12 / dist12 / 40.0;
  double scale = 0.0;
  if (bondType == 14)
  {
    scale = ENERGY_SCALE_FACTOR_BOND_14;
  }
  else if (bondType == 15)
  {
    scale = ENERGY_SCALE_FACTOR_BOND_15;
  }
  energy *= scale;
  *elec = energy;

  if (ENERGY_DEBUG_MODE_ELEC)
  {
    if (fabs(energy) > 0.0 && dist12 < 4.0)
    {
      printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %5.2f, elec: %5.2f\n",
        AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1),
        AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
        bondType, dist12, energy);
    }
  }
  return Success;
}

int LKDesolvationEnergyAtomAndAtom(Atom* pAtom1, Atom* pAtom2, double distance, int bondType, double* energyP, double* energyH)
{
  if (bondType == 12 || bondType == 13)
  {
    return Success;
  }
  if (AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2))
  {
    return Success;
  }
  if (distance > ENERGY_DISTANCE_CUTOFF)
  {
    return Success;
  }
  double volume1 = pAtom1->EEF1_volume;
  double volume2 = pAtom2->EEF1_volume;
  double dGFreeAtom1 = pAtom1->EEF1_freeDG;
  double dGFreeAtom2 = pAtom2->EEF1_freeDG;
  double coefficient = -0.089793561062582974; // 0.5/(pi^1.5)
  double r1 = pAtom1->vdw_radius * RADIUS_SCALE_FOR_DESOLV;
  double r2 = pAtom2->vdw_radius * RADIUS_SCALE_FOR_DESOLV;
  double r12 = r1 + r2;

  distance = distance < r12 ? r12 : distance;
  double lamda1 = pAtom1->EEF1_lamda_ * distance * distance;
  double lamda2 = pAtom2->EEF1_lamda_ * distance * distance;
  double x1 = (distance - r1) / pAtom1->EEF1_lamda_;
  double x2 = (distance - r2) / pAtom2->EEF1_lamda_;

  double desolv12 = coefficient * volume2 * dGFreeAtom1 / lamda1;
  desolv12 *= exp(-1.0 * x1 * x1);
  double desolv21 = coefficient * volume1 * dGFreeAtom2 / lamda2;
  desolv21 *= exp(-1.0 * x2 * x2);
  if (pAtom1->polarity == Type_AtomPolarity_P || pAtom1->polarity == Type_AtomPolarity_C)
  {
    *energyP += desolv12;
  }
  else
  {
    *energyH += desolv12;
  }
  if (pAtom2->polarity == Type_AtomPolarity_P || pAtom2->polarity == Type_AtomPolarity_C)
  {
    *energyP += desolv21;
  }
  else
  {
    *energyH += desolv21;
  }

  if (ENERGY_DEBUG_MODE_DESOLV)
  {
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

int SSbondEnergyAtomAndAtom(Atom* pAtomS1, Atom* pAtomS2, Atom* pAtomCB1, Atom* pAtomCB2, Atom* pAtomCA1, Atom* pAtomCA2, double* sse)
{
  XYZ xyzC1S1 = XYZDifference(&pAtomCB1->xyz, &pAtomS1->xyz);
  XYZ xyzS1S2 = XYZDifference(&pAtomS1->xyz, &pAtomS2->xyz);
  XYZ xyzS2C2 = XYZDifference(&pAtomS2->xyz, &pAtomCB2->xyz);
  double dSS = XYZDistance(&pAtomS1->xyz, &pAtomS2->xyz);
  double aC1S1S2 = RadToDeg(PI - XYZAngle(&xyzS1S2, &xyzS2C2));
  double aC2S2S1 = RadToDeg(PI - XYZAngle(&xyzC1S1, &xyzS1S2));
  double xC1S1S2C2 = GetTorsionAngle(&pAtomCB1->xyz, &pAtomS1->xyz, &pAtomS2->xyz, &pAtomCB2->xyz);
  double xCA1CB1SG1SG2 = GetTorsionAngle(&pAtomCA1->xyz, &pAtomCB1->xyz, &pAtomS1->xyz, &pAtomS2->xyz);
  double xCA2CB2SG2SG1 = GetTorsionAngle(&pAtomCA2->xyz, &pAtomCB2->xyz, &pAtomS2->xyz, &pAtomS1->xyz);
  *sse = (0.8 * (1 - exp(-10.0 * (dSS - SSBOND_DISTANCE))) * (1 - exp(-10.0 * (dSS - SSBOND_DISTANCE))))
    + 0.005 * (aC1S1S2 - SSBOND_ANGLE) * (aC1S1S2 - SSBOND_ANGLE)
    + 0.005 * (aC2S2S1 - SSBOND_ANGLE) * (aC2S2S1 - SSBOND_ANGLE)
    + cos(2.0 * xC1S1S2C2) + 1.0
    + 1.25 * sin(xCA1CB1SG1SG2 + 2.0 * PI / 3.0) - 1.75
    + 1.25 * sin(xCA2CB2SG2SG1 + 2.0 * PI / 3.0) - 1.75;
  if (*sse > 0.0)
  {
    *sse = 0.0;
  }
  return Success;
}


#define RESIDUE_PAIRWISE_ENERGY_COMPUTATION
int AminoAcidReferenceEnergy(char* AAname, double energyTerm[MAX_ENERGY_TERM])
{
  if (strcmp(AAname, "ALA") == 0)
  {
    energyTerm[1] += 1.0;
  }
  else if (strcmp(AAname, "CYS") == 0)
  {
    energyTerm[2] += 1.0;
  }
  else if (strcmp(AAname, "ASP") == 0)
  {
    energyTerm[3] += 1.0;
  }
  else if (strcmp(AAname, "GLU") == 0)
  {
    energyTerm[4] += 1.0;
  }
  else if (strcmp(AAname, "PHE") == 0)
  {
    energyTerm[5] += 1.0;
  }
  else if (strcmp(AAname, "GLY") == 0)
  {
    energyTerm[6] += 1.0;
  }
  else if (strcmp(AAname, "HIS") == 0)
  {
    energyTerm[7] += 1.0;
  }
  else if (strcmp(AAname, "HSE") == 0)
  {
    energyTerm[7] += 1.0;
  }
  else if (strcmp(AAname, "HSD") == 0)
  {
    energyTerm[7] += 1.0;
  }
  else if (strcmp(AAname, "ILE") == 0)
  {
    energyTerm[8] += 1.0;
  }
  else if (strcmp(AAname, "LYS") == 0)
  {
    energyTerm[9] += 1.0;
  }
  else if (strcmp(AAname, "LEU") == 0)
  {
    energyTerm[10] += 1.0;
  }
  else if (strcmp(AAname, "MET") == 0)
  {
    energyTerm[11] += 1.0;
  }
  else if (strcmp(AAname, "ASN") == 0)
  {
    energyTerm[12] += 1.0;
  }
  else if (strcmp(AAname, "PRO") == 0)
  {
    energyTerm[13] += 1.0;
  }
  else if (strcmp(AAname, "GLN") == 0)
  {
    energyTerm[14] += 1.0;
  }
  else if (strcmp(AAname, "ARG") == 0)
  {
    energyTerm[15] += 1.0;
  }
  else if (strcmp(AAname, "SER") == 0)
  {
    energyTerm[16] += 1.0;
  }
  else if (strcmp(AAname, "THR") == 0)
  {
    energyTerm[17] += 1.0;
  }
  else if (strcmp(AAname, "VAL") == 0)
  {
    energyTerm[18] += 1.0;
  }
  else if (strcmp(AAname, "TRP") == 0)
  {
    energyTerm[19] += 1.0;
  }
  else if (strcmp(AAname, "TYR") == 0)
  {
    energyTerm[20] += 1.0;
  }
  return Success;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// the following are energy between residue and residue
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int EnergyIntraResidue(Residue* pThis, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < ResidueGetAtomCount(pThis); i++)
  {
    Atom* pAtom1 = ResidueGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
    {
      continue;
    }
    for (int j = i + 1; j < ResidueGetAtomCount(pThis); j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pThis, j);
      if (pAtom2->isXyzValid == FALSE)
      {
        continue;
      }
      if (pAtom2->isBBAtom == TRUE && pAtom1->isBBAtom == TRUE)
      {
        continue;
      }
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
      {
        continue;
      }
      if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
      {
        if (strcmp(ResidueGetName(pThis), "ILE") == 0
          || strcmp(ResidueGetName(pThis), "MET") == 0
          || strcmp(ResidueGetName(pThis), "GLN") == 0
          || strcmp(ResidueGetName(pThis), "GLU") == 0
          || strcmp(ResidueGetName(pThis), "LYS") == 0
          || strcmp(ResidueGetName(pThis), "ARG") == 0)
        {
          int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom1), AtomGetName(pAtom2), ResidueGetBonds(pThis));
          if (bondType == 12 || bondType == 13)
          {
            continue;
          }
          double vdwAtt = 0, vdwRep = 0, desolvP = 0, desolvH = 0;
          VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwAtt);
          VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwRep);
          LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desolvP, &desolvH);
          energyTerms[21] += vdwAtt;
          energyTerms[22] += vdwRep;
          energyTerms[24] += desolvP;
          energyTerms[25] += desolvH;
        }
      }
      else if ((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE))
      {
        if (strcmp(AtomGetName(pAtom1), "CB") == 0 || strcmp(AtomGetName(pAtom2), "CB") == 0)
        {
          continue;
        }
        int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom1), AtomGetName(pAtom2), ResidueGetBonds(pThis));
        if (bondType == 12 || bondType == 13)
        {
          continue;
        }

        double vdwAtt = 0, vdwRep = 0, desolvP = 0, desolvH = 0, ele = 0;
        VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwAtt);
        VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwRep);
        LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desolvP, &desolvH);
        ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
        energyTerms[21] += vdwAtt;
        energyTerms[22] += vdwRep;
        energyTerms[23] += ele;
        energyTerms[24] += desolvP;
        energyTerms[25] += desolvH;
        if (distance < HB_HA_DIST_CUTOFF_MAX)
        {
          double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
          if (pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE)
          {
            HBondEnergyAtomAndAtom(pAtom1, pAtom2, ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom2)),
              distance, bondType, &hbtot, &hbd, &hbt, &hbp);
          }
          else if (pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE)
          {
            HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom2)), ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
              distance, bondType, &hbtot, &hbd, &hbt, &hbp);
          }
          energyTerms[26] += hbd;
          energyTerms[27] += hbt;
          energyTerms[28] += hbp;
        }
      }
    }

  }
  return Success;
}


int EnergyResidueAndNextResidue(Residue* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
  {
    Atom* pAtom1 = ResidueGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
    {
      continue;
    }
    for (int j = 0;j < ResidueGetAtomCount(pOther);j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)
      {
        continue;
      }
      if (pAtom1->isBBAtom && pAtom2->isBBAtom)
      {
        continue;
      }
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
      {
        continue;
      }
      if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
      {
        int bondType = 15;
        double att = 0, rep = 0, desP = 0, desH = 0, ele = 0;
        VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
        VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
        LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
        ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
        energyTerms[31] += att;
        energyTerms[32] += rep;
        energyTerms[33] += ele;
        energyTerms[34] += desP;
        energyTerms[35] += desH;
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        if (pThis->posInChain - pOther->posInChain <= 2 && pThis->posInChain - pOther->posInChain >= -2)
        {
          hbd *= HB_LOCAL_SCALE;
          hbt *= HB_LOCAL_SCALE;
          hbp *= HB_LOCAL_SCALE;
        }
        if (pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE)
        {
          energyTerms[41] += hbd;
          energyTerms[42] += hbt;
          energyTerms[43] += hbp;
        }
        else if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
        {
          energyTerms[47] += hbd;
          energyTerms[48] += hbt;
          energyTerms[49] += hbp;
        }
        else
        {
          energyTerms[44] += hbd;
          energyTerms[45] += hbt;
          energyTerms[46] += hbp;
        }
      }
      else if ((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE)
        || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE))
      {
        int bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1), AtomGetName(pAtom2), ResidueGetName(pOther));
        if (bondType == 12 || bondType == 13)
        {
          continue;
        }
        double vdwAtt = 0, vdwRep = 0, desolvP = 0, desolvH = 0, ele = 0;
        VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwAtt);
        VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwRep);
        LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desolvP, &desolvH);
        ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
        energyTerms[31] += vdwAtt;
        energyTerms[32] += vdwRep;
        energyTerms[33] += ele;
        energyTerms[34] += desolvP;
        energyTerms[35] += desolvH;
        if (distance < HB_HA_DIST_CUTOFF_MAX)
        {
          double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
          if (pAtom1->isHBatomH && pAtom2->isHBatomA)
          {
            HBondEnergyAtomAndAtom(pAtom1, pAtom2, ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
              distance, bondType, &hbtot, &hbd, &hbt, &hbp);
          }
          else if (pAtom2->isHBatomH && pAtom1->isHBatomA)
          {
            HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
              distance, bondType, &hbtot, &hbd, &hbt, &hbp);
          }
          if (pThis->posInChain - pOther->posInChain <= 2 && pThis->posInChain - pOther->posInChain >= -2)
          {
            hbd *= HB_LOCAL_SCALE;
            hbt *= HB_LOCAL_SCALE;
            hbp *= HB_LOCAL_SCALE;
          }
          if (pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE)
          {
            energyTerms[41] += hbd;
            energyTerms[42] += hbt;
            energyTerms[43] += hbp;
          }
          else if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
          {
            energyTerms[47] += hbd;
            energyTerms[48] += hbt;
            energyTerms[49] += hbp;
          }
          else
          {
            energyTerms[44] += hbd;
            energyTerms[45] += hbt;
            energyTerms[46] += hbp;
          }
        }
      }
    }
  }
  return Success;
}

int EnergyResidueAndOtherResidueSameChain(Residue* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
  {
    Atom* pAtom1 = ResidueGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)continue;
    for (int j = 0;j < ResidueGetAtomCount(pOther);j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      int bondType = 15;
      double vdwAtt = 0, vdwRep = 0, ele = 0, desolvP = 0, desolvH = 0;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwAtt);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwRep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desolvP, &desolvH);
      energyTerms[31] += vdwAtt;
      energyTerms[32] += vdwRep;
      energyTerms[33] += ele;
      energyTerms[34] += desolvP;
      energyTerms[35] += desolvH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom2->isHBatomH && pAtom1->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        if (pThis->posInChain - pOther->posInChain <= 2 && pThis->posInChain - pOther->posInChain >= -2)
        {
          hbd *= HB_LOCAL_SCALE;
          hbt *= HB_LOCAL_SCALE;
          hbp *= HB_LOCAL_SCALE;
        }
        if (pAtom1->isBBAtom && pAtom2->isBBAtom)
        {
          energyTerms[41] += hbd;
          energyTerms[42] += hbt;
          energyTerms[43] += hbp;
        }
        else if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
        {
          energyTerms[47] += hbd;
          energyTerms[48] += hbt;
          energyTerms[49] += hbp;
        }
        else
        {
          energyTerms[44] += hbd;
          energyTerms[45] += hbt;
          energyTerms[46] += hbp;
        }
      }
    }
  }

  // consider the SSBOND
  if (strcmp("CYS", ResidueGetName(pThis)) == 0 && strcmp("CYS", ResidueGetName(pOther)) == 0)
  {
    Atom* pSG1 = NULL, * pSG2 = NULL, * pCB1 = NULL, * pCB2 = NULL, * pCA1 = NULL, * pCA2 = NULL;
    for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
    {
      Atom* pAtom = ResidueGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA1 = pAtom;
    }
    for (int i = 0;i < ResidueGetAtomCount(pOther);i++)
    {
      Atom* pAtom = ResidueGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDistance(&pSG1->xyz, &pSG2->xyz);
      if (dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[36] += ssbond;
      }
    }
  }

  return Success;
}


int EnergyResidueAndOtherResidueDiffChain(Residue* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < ResidueGetAtomCount(pThis); i++)
  {
    Atom* pAtom1 = ResidueGetAtom(pThis, i);
    if (!pAtom1->isXyzValid)continue;
    for (int j = 0; j < ResidueGetAtomCount(pOther); j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if (!pAtom2->isXyzValid)continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      int bondType = 15;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[51] += att;
      energyTerms[52] += rep;
      energyTerms[53] += ele;
      energyTerms[54] += desP;
      energyTerms[55] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        if (pAtom1->isBBAtom && pAtom2->isBBAtom)
        {
          energyTerms[61] += hbd;
          energyTerms[62] += hbt;
          energyTerms[63] += hbp;
        }
        else if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
        {
          energyTerms[67] += hbd;
          energyTerms[68] += hbt;
          energyTerms[69] += hbp;
        }
        else
        {
          energyTerms[64] += hbd;
          energyTerms[65] += hbt;
          energyTerms[66] += hbp;
        }
      }
    }
  }

  //consider the SSBOND
  if (strcmp("CYS", ResidueGetName(pThis)) == 0 && strcmp("CYS", ResidueGetName(pOther)) == 0)
  {
    Atom* pSG1 = NULL, * pSG2 = NULL, * pCB1 = NULL, * pCB2 = NULL, * pCA1 = NULL, * pCA2 = NULL;
    for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
    {
      Atom* pAtom = ResidueGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA1 = pAtom;
    }
    for (int i = 0;i < ResidueGetAtomCount(pOther);i++)
    {
      Atom* pAtom = ResidueGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDistance(&pSG1->xyz, &pSG2->xyz);
      if (dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[56] += ssbond;
      }
    }
  }

  return Success;
}


int EnergyResidueAndLigandResidue(Residue* pProtein, Residue* pLigand, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < ResidueGetAtomCount(pProtein); i++)
  {
    Atom* pAtom1 = ResidueGetAtom(pProtein, i);
    if (pAtom1->isXyzValid == FALSE)continue;
    for (int j = 0; j < ResidueGetAtomCount(pLigand); j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pLigand, j);
      if (pAtom2->isXyzValid == FALSE)continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      int bondType = 15;
      double vdwAtt = 0, vdwRep = 0, ele = 0, desolvP = 0, desolvH = 0;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwAtt);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwRep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desolvP, &desolvH);
      energyTerms[71] += vdwAtt;
      energyTerms[72] += vdwRep;
      energyTerms[73] += ele;
      energyTerms[74] += desolvP;
      energyTerms[75] += desolvH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, ResidueGetAtomByName(pProtein, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pLigand, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pLigand, AtomGetHbDorB(pAtom2)), ResidueGetAtomByName(pProtein, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
        {
          energyTerms[84] += hbd;
          energyTerms[85] += hbt;
          energyTerms[86] += hbp;
        }
        else
        {
          energyTerms[81] += hbd;
          energyTerms[82] += hbt;
          energyTerms[83] += hbp;
        }
      }
    }
  }
  return Success;
}


#define ROTAMER_PAIRWISE_ENERGY_COMPUTATION
////////////////////////////////////////////////////////////////////////////////////////////
//these functions are used to calculate energy for rots
////////////////////////////////////////////////////////////////////////////////////////////

int EnergyIntraRotamer(Rotamer* pThis, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < RotamerGetAtomCount(pThis); ++i)
  {
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE) continue;
    for (int j = i + 1; j < RotamerGetAtomCount(pThis);++j)
    {
      Atom* pAtom2 = RotamerGetAtom(pThis, j);
      if (pAtom2->isXyzValid == FALSE) continue;
      if (pAtom2->isBBAtom && pAtom1->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
      {
        if (strcmp(RotamerGetType(pThis), "ILE") == 0 || strcmp(RotamerGetType(pThis), "MET") == 0 ||
          strcmp(RotamerGetType(pThis), "GLN") == 0 || strcmp(RotamerGetType(pThis), "GLU") == 0 ||
          strcmp(RotamerGetType(pThis), "LYS") == 0 || strcmp(RotamerGetType(pThis), "ARG") == 0)
        {
          int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom1), AtomGetName(pAtom2), RotamerGetBonds(pThis));
          if (bondType == 12 || bondType == 13) continue;
          double vdwAtt = 0, vdwRep = 0, desolvP = 0, desolvH = 0;
          VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwAtt);
          VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwRep);
          LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desolvP, &desolvH);
          energyTerms[21] += vdwAtt;
          energyTerms[22] += vdwRep;
          energyTerms[24] += desolvP;
          energyTerms[25] += desolvH;
        }
      }
      else if ((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE))
      {
        if (strcmp(AtomGetName(pAtom1), "CB") == 0 || strcmp(AtomGetName(pAtom2), "CB") == 0) continue;
        int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom1), AtomGetName(pAtom2), RotamerGetBonds(pThis));
        if (bondType == 12 || bondType == 13) continue;
        double vdwAtt = 0, vdwRep = 0, desolvP = 0, desolvH = 0, ele = 0;
        VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwAtt);
        VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &vdwRep);
        LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desolvP, &desolvH);
        ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
        energyTerms[21] += vdwAtt;
        energyTerms[22] += vdwRep;
        energyTerms[23] += ele;
        energyTerms[24] += desolvP;
        energyTerms[25] += desolvH;
        if (distance < HB_HA_DIST_CUTOFF_MAX)
        {
          double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
          if (pAtom1->isHBatomH == TRUE && pAtom2->isHBatomA == TRUE)
          {
            HBondEnergyAtomAndAtom(pAtom1, pAtom2, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom2)),
              distance, bondType, &hbtot, &hbd, &hbt, &hbp);
          }
          else if (pAtom2->isHBatomH == TRUE && pAtom1->isHBatomA == TRUE)
          {
            HBondEnergyAtomAndAtom(pAtom2, pAtom1, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom2)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
              distance, bondType, &hbtot, &hbd, &hbt, &hbp);
          }
          energyTerms[26] += hbd;
          energyTerms[27] += hbt;
          energyTerms[28] += hbp;
        }
      }
    }
  }
  return Success;
}


int EnergyRotamerAndRotamerSameChain(Rotamer* pThis, Rotamer* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0;i < RotamerGetAtomCount(pThis);i++)
  {
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE) continue;
    if (pAtom1->isBBAtom) continue;
    for (int j = 0;j < RotamerGetAtomCount(pOther);j++)
    {
      Atom* pAtom2 = RotamerGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE) continue;
      if (pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      int bondType = 15;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[31] += att;
      energyTerms[32] += rep;
      energyTerms[33] += ele;
      energyTerms[34] += desP;
      energyTerms[35] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), RotamerGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
          if (pThis->posInChain - pOther->posInChain <= 2 && pThis->posInChain - pOther->posInChain >= -2)
          {
            hbd *= HB_LOCAL_SCALE;
            hbt *= HB_LOCAL_SCALE;
            hbp *= HB_LOCAL_SCALE;
          }
          // only hbscsc
          energyTerms[47] += hbd;
          energyTerms[48] += hbt;
          energyTerms[49] += hbp;
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, RotamerGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
          if (pThis->posInChain - pOther->posInChain <= 2 && pThis->posInChain - pOther->posInChain >= -2)
          {
            hbd *= HB_LOCAL_SCALE;
            hbt *= HB_LOCAL_SCALE;
            hbp *= HB_LOCAL_SCALE;
          }
          // only hbscsc
          energyTerms[47] += hbd;
          energyTerms[48] += hbt;
          energyTerms[49] += hbp;
        }
      }
    }
  }

  //consider the SSBOND
  if (strcmp("CYS", RotamerGetType(pThis)) == 0 && strcmp("CYS", RotamerGetType(pOther)) == 0)
  {
    Atom* pSG1 = NULL, * pSG2 = NULL, * pCB1 = NULL, * pCB2 = NULL, * pCA1 = NULL, * pCA2 = NULL;
    for (int i = 0;i < RotamerGetAtomCount(pThis);i++)
    {
      Atom* pAtom = RotamerGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA1 = pAtom;
    }
    for (int i = 0;i < RotamerGetAtomCount(pOther);i++)
    {
      Atom* pAtom = RotamerGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDistance(&pSG1->xyz, &pSG2->xyz);
      if (dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[36] += ssbond;
      }
    }
  }

  return Success;
}

int EnergyRotamerAndRotamerDiffChain(Rotamer* pThis, Rotamer* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < RotamerGetAtomCount(pThis); i++)
  {
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE) continue;
    if (pAtom1->isBBAtom) continue;
    for (int j = 0; j < RotamerGetAtomCount(pOther); j++)
    {
      Atom* pAtom2 = RotamerGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE) continue;
      if (pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      int bondType = 15;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[51] += att;
      energyTerms[52] += rep;
      energyTerms[53] += ele;
      energyTerms[54] += desP;
      energyTerms[55] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), RotamerGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, RotamerGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        /*if(pAtom1->isBBAtom && pAtom2->isBBAtom){
          energyTerms[61]+=hbd;
          energyTerms[62]+=hbt;
          energyTerms[63]+=hbp;
        }
        else if(pAtom1->isBBAtom==FALSE && pAtom2->isBBAtom==FALSE){*/
        energyTerms[67] += hbd;
        energyTerms[68] += hbt;
        energyTerms[69] += hbp;
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
  if (strcmp("CYS", RotamerGetType(pThis)) == 0 && strcmp("CYS", RotamerGetType(pOther)) == 0)
  {
    Atom* pSG1 = NULL, * pSG2 = NULL, * pCB1 = NULL, * pCB2 = NULL, * pCA1 = NULL, * pCA2 = NULL;
    for (int i = 0;i < RotamerGetAtomCount(pThis);i++)
    {
      Atom* pAtom = RotamerGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA1 = pAtom;
    }
    for (int i = 0;i < RotamerGetAtomCount(pOther);i++)
    {
      Atom* pAtom = RotamerGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDistance(&pSG1->xyz, &pSG2->xyz);
      if (dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[56] += ssbond;
      }
    }
  }
  return Success;
}



int EnergyRotamerAndLigandRotamer(Rotamer* pThis, Rotamer* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < RotamerGetAtomCount(pThis); i++)
  {
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE) continue;
    if (pAtom1->isBBAtom) continue;
    for (int j = 0; j < RotamerGetAtomCount(pOther); j++)
    {
      Atom* pAtom2 = RotamerGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE) continue;
      if (pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      int bondType = 15;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), RotamerGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, RotamerGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        energyTerms[84] += hbd;
        energyTerms[85] += hbt;
        energyTerms[86] += hbp;
      }
    }
  }

  return Success;
}



int EnergyRotamerAndFixedResidueSameChain(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  int neighborCheck = 0;
  if (RotamerGetPosInChain(pThis) + 1 == ResidueGetPosInChain(pOther)) neighborCheck = 12;
  else if (RotamerGetPosInChain(pThis) - 1 == ResidueGetPosInChain(pOther)) neighborCheck = 21;

  for (int i = 0; i < RotamerGetAtomCount(pThis); i++)
  {
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE) continue;
    if (pAtom1->isBBAtom) continue;
    for (int j = 0; j < ResidueGetAtomCount(pOther); j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      int bondType = 15;
      if (neighborCheck == 12) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1), AtomGetName(pAtom2), ResidueGetName(pOther));
      else if (neighborCheck == 21) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom2), AtomGetName(pAtom1), RotamerGetType(pThis));
      if (bondType == 12 || bondType == 13) continue;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[31] += att;
      energyTerms[32] += rep;
      energyTerms[33] += ele;
      energyTerms[34] += desP;
      energyTerms[35] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        if (pThis->posInChain - pOther->posInChain <= 2 && pThis->posInChain - pOther->posInChain >= -2)
        {
          hbd *= HB_LOCAL_SCALE;
          hbt *= HB_LOCAL_SCALE;
          hbp *= HB_LOCAL_SCALE;
        }
        /*if(pAtom1->isBBAtom && pAtom2->isBBAtom){
          energyTerms[41]+=hbd;
          energyTerms[42]+=hbt;
          energyTerms[43]+=hbp;
        }
        else */if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
        {
          energyTerms[47] += hbd;
          energyTerms[48] += hbt;
          energyTerms[49] += hbp;
        }
        else
        {
          energyTerms[44] += hbd;
          energyTerms[45] += hbt;
          energyTerms[46] += hbp;
        }

      }
    }
  }

  //consider the SSBOND
  if (strcmp("CYS", RotamerGetType(pThis)) == 0 && strcmp("CYS", ResidueGetName(pOther)) == 0)
  {
    Atom* pSG1 = NULL, * pSG2 = NULL, * pCB1 = NULL, * pCB2 = NULL, * pCA1 = NULL, * pCA2 = NULL;
    for (int i = 0;i < RotamerGetAtomCount(pThis);i++)
    {
      Atom* pAtom = RotamerGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA1 = pAtom;
    }
    for (int i = 0;i < ResidueGetAtomCount(pOther);i++)
    {
      Atom* pAtom = ResidueGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDistance(&pSG1->xyz, &pSG2->xyz);
      if (dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[36] += ssbond;
      }
    }
  }


  return Success;
}

int EnergyRotamerAndFixedResidueDiffChain(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < RotamerGetAtomCount(pThis); i++)
  {
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if (!pAtom1->isXyzValid) continue;
    if (pAtom1->isBBAtom) continue;
    for (int j = 0; j < ResidueGetAtomCount(pOther); j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if (!pAtom2->isXyzValid) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      int bondType = 15;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[51] += att;
      energyTerms[52] += rep;
      energyTerms[53] += ele;
      energyTerms[54] += desP;
      energyTerms[55] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        /*if(pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE){
          energyTerms[61]+=hbd;
          energyTerms[62]+=hbt;
          energyTerms[63]+=hbp;
        }
        else */if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
        {
          energyTerms[67] += hbd;
          energyTerms[68] += hbt;
          energyTerms[69] += hbp;
        }
        else
        {
          energyTerms[64] += hbd;
          energyTerms[65] += hbt;
          energyTerms[66] += hbp;
        }
      }
    }
  }

  //consider the SSBOND
  if (strcmp("CYS", RotamerGetType(pThis)) == 0 && strcmp("CYS", ResidueGetName(pOther)) == 0)
  {
    Atom* pSG1 = NULL, * pSG2 = NULL, * pCB1 = NULL, * pCB2 = NULL, * pCA1 = NULL, * pCA2 = NULL;
    for (int i = 0;i < RotamerGetAtomCount(pThis);i++)
    {
      Atom* pAtom = RotamerGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA1 = pAtom;
    }
    for (int i = 0;i < ResidueGetAtomCount(pOther);i++)
    {
      Atom* pAtom = ResidueGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDistance(&pSG1->xyz, &pSG2->xyz);
      if (dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[56] += ssbond;
      }
    }
  }


  return Success;
}


int EnergyRotamerAndFixedLigResidue(Rotamer* pThis, Residue* pLigand, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < RotamerGetAtomCount(pThis); i++)
  {
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE) continue;
    if (pAtom1->isBBAtom) continue;
    for (int j = 0; j < ResidueGetAtomCount(pLigand); j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pLigand, j);
      if (pAtom2->isXyzValid == FALSE) continue;
      //if(pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      int bondType = 15;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pLigand, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pLigand, AtomGetHbDorB(pAtom2)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        //if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
        energyTerms[84] += hbd;
        energyTerms[85] += hbt;
        energyTerms[86] += hbp;
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



int EnergyLigRotamerAndFixedResidue(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < RotamerGetAtomCount(pThis); i++)
  {
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE) continue;
    for (int j = 0; j < ResidueGetAtomCount(pOther); j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      int bondType = 15;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
        {
          energyTerms[84] += hbd;
          energyTerms[85] += hbt;
          energyTerms[86] += hbp;
        }
        else
        {
          energyTerms[81] += hbd;
          energyTerms[82] += hbt;
          energyTerms[83] += hbp;
        }
      }
    }
  }
  return Success;
}



int EnergyRotamerAndDesignResidueSameChain(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  int neighborCheck = 0;
  if (RotamerGetPosInChain(pThis) + 1 == ResidueGetPosInChain(pOther)) neighborCheck = 12;
  else if (RotamerGetPosInChain(pThis) - 1 == ResidueGetPosInChain(pOther)) neighborCheck = 21;

  for (int i = 0;i < RotamerGetAtomCount(pThis);i++)
  {
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE) continue;
    if (pAtom1->isBBAtom) continue;
    for (int j = 0;j < ResidueGetAtomCount(pOther);j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE) continue;
      if (pAtom2->isBBAtom == FALSE) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      int bondType = 15;
      if (neighborCheck == 12) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1), AtomGetName(pAtom2), ResidueGetName(pOther));
      else if (neighborCheck == 21) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom2), AtomGetName(pAtom1), RotamerGetType(pThis));
      if (bondType == 12 || bondType == 13) continue;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[31] += att;
      energyTerms[32] += rep;
      energyTerms[33] += ele;
      energyTerms[34] += desP;
      energyTerms[35] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        if (pThis->posInChain - pOther->posInChain <= 2 && pThis->posInChain - pOther->posInChain >= -2)
        {
          hbd *= HB_LOCAL_SCALE;
          hbt *= HB_LOCAL_SCALE;
          hbp *= HB_LOCAL_SCALE;
        }
        energyTerms[44] += hbd;
        energyTerms[45] += hbt;
        energyTerms[46] += hbp;
      }
    }
  }

  return Success;
}



int EnergyRotamerAndDesignResidueDiffChain(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < RotamerGetAtomCount(pThis); i++)
  {
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if (!pAtom1->isXyzValid) continue;
    if (pAtom1->isBBAtom) continue;
    for (int j = 0; j < ResidueGetAtomCount(pOther); j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if (!pAtom2->isXyzValid) continue;
      if (!pAtom2->isBBAtom) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      int bondType = 15;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[51] += att;
      energyTerms[52] += rep;
      energyTerms[53] += ele;
      energyTerms[54] += desP;
      energyTerms[55] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        energyTerms[64] += hbd;
        energyTerms[65] += hbt;
        energyTerms[66] += hbp;
      }
    }
  }
  return Success;
}





int EnergyLigRotamerAndDesignResidue(Rotamer* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < RotamerGetAtomCount(pThis); i++)
  {
    Atom* pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE) continue;
    if (pAtom1->isBBAtom) continue;
    for (int j = 0; j < ResidueGetAtomCount(pOther); j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE) continue;
      if (pAtom2->isBBAtom == FALSE) continue;
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      int bondType = 15;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[71] += att;
      energyTerms[72] += rep;
      energyTerms[73] += ele;
      energyTerms[74] += desP;
      energyTerms[75] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), RotamerGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
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
int RamaTableReadFromFile(RamaTable* pRama, char* ramafile)
{
  FILE* pFile = fopen(ramafile, "r");
  if (!pFile) return IOError;
  char buffer[MAX_LEN_ONE_LINE_CONTENT + 1];
  while (fgets(buffer, MAX_LEN_ONE_LINE_CONTENT, pFile))
  {
    if (buffer[0] == ' ' || buffer[0] == '#') continue;
    int phipsi[2];
    double aae[20];
    sscanf(buffer, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
      &phipsi[0], &phipsi[1], &aae[0], &aae[1], &aae[2], &aae[3], &aae[4], &aae[5], &aae[6], &aae[7], &aae[8], &aae[9],
      &aae[10], &aae[11], &aae[12], &aae[13], &aae[14], &aae[15], &aae[16], &aae[17], &aae[18], &aae[19]);
    int phiindex = (phipsi[0] + 180) / 10;
    int psiindex = (phipsi[1] + 180) / 10;
    for (int j = 0;j < 20;j++)
    {
      pRama->ramatable[phiindex][psiindex][j] = aae[j];
    }
  }
  fclose(pFile);
  return Success;
}

int AApropensityTableReadFromFile(AAppTable* pAAppTable, char* aappfile)
{
  FILE* pFile = fopen(aappfile, "r");
  if (!pFile) return IOError;
  char buffer[MAX_LEN_ONE_LINE_CONTENT + 1];
  while (fgets(buffer, MAX_LEN_ONE_LINE_CONTENT, pFile))
  {
    if (buffer[0] == ' ' || buffer[0] == '#') continue;
    int phipsi[2];
    double aae[20];
    sscanf(buffer, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
      &phipsi[0], &phipsi[1], &aae[0], &aae[1], &aae[2], &aae[3], &aae[4], &aae[5], &aae[6], &aae[7], &aae[8], &aae[9],
      &aae[10], &aae[11], &aae[12], &aae[13], &aae[14], &aae[15], &aae[16], &aae[17], &aae[18], &aae[19]);
    int phiindex = (phipsi[0] + 180) / 10;
    int psiindex = (phipsi[1] + 180) / 10;
    for (int j = 0;j < 20;j++)
    {
      pAAppTable->aapptable[phiindex][psiindex][j] = aae[j];
    }
  }
  fclose(pFile);
  return Success;
}

int AminoAcidPropensityAndRamachandranEnergy(Residue* pThis, AAppTable* pAAppTable, RamaTable* pRama)
{
  char aa1 = AA3ToAA1(ResidueGetName(pThis));
  int aaindex = AA1GetIndex(aa1);
  int phiindex = ((int)pThis->phipsi[0] + 180) / 10;
  int psiindex = ((int)pThis->phipsi[1] + 180) / 10;
  pThis->aapp += pAAppTable->aapptable[phiindex][psiindex][aaindex];
  pThis->rama += pRama->ramatable[phiindex][psiindex][aaindex];
  return Success;
}


int RotamerPropensityAndRamachandranEnergy(Rotamer* pThis, Residue* pResidue, AAppTable* pAAppTable, RamaTable* pRama, double energyTerms[MAX_ENERGY_TERM])
{
  char aa1 = AA3ToAA1(RotamerGetType(pThis));
  int aaindex = AA1GetIndex(aa1);
  int phiindex = ((int)pResidue->phipsi[0] + 180) / 10;
  int psiindex = ((int)pResidue->phipsi[1] + 180) / 10;
  energyTerms[91] += pAAppTable->aapptable[phiindex][psiindex][aaindex];
  energyTerms[92] += pRama->ramatable[phiindex][psiindex][aaindex];
  return Success;
}


int AminoAcidDunbrackEnergy(Residue* pThis, BBdepRotamerLib* pBBdepRotLib)
{
  if (strcmp(ResidueGetName(pThis), "ALA") == 0
    || strcmp(ResidueGetName(pThis), "GLY") == 0
    || AA3GetIndex(ResidueGetName(pThis)) < 0 // indicates a non-canonical residue
    )
  {
    return Success;
  }
  if (pThis->isSCIntact)
  {
    int binIdx = ((int)(pThis->phipsi[0] + 180) / 10) * 36 + (int)(pThis->phipsi[1] + 180) / 10;
    RotLibPhiPsi* pRotLibPhiPsi = &pBBdepRotLib->rotlibphipsis[binIdx];
    int rotTypeIdx = -1;
    StringArrayFind(&pRotLibPhiPsi->rotTypes, ResidueGetName(pThis), &rotTypeIdx);
    DoubleArray* pTorsionsArray = pRotLibPhiPsi->torsions[rotTypeIdx];
    DoubleArray* pDeviationsArray = pRotLibPhiPsi->deviations[rotTypeIdx];
    int matchIdx = -1;
    double pMatch;
    double pMin = 10;
    for (int i = 0;i < IntArrayGet(&pRotLibPhiPsi->rotamerCounts, rotTypeIdx);i++)
    {
      double p = DoubleArrayGet(&pRotLibPhiPsi->probability[rotTypeIdx], i);
      if (p < CUT_EXCL_LOW_PROB_ROT) break;
      if (p < pMin) pMin = p;
      DoubleArray* pTorsions = &pTorsionsArray[i];
      DoubleArray* pDeviations = &pDeviationsArray[i];

      BOOL match = TRUE;
      for (int j = 0;j < DoubleArrayGetLength(&pThis->Xs);j++)
      {
        //use a strict criteria
        //double min=DoubleArrayGet(pTorsions,j)-DegToRad(5.0);
        //double max=DoubleArrayGet(pTorsions,j)+DegToRad(5.0);
        double min = DoubleArrayGet(pTorsions, j) - DoubleArrayGet(pDeviations, j);
        double max = DoubleArrayGet(pTorsions, j) + DoubleArrayGet(pDeviations, j);
        double torsion = DoubleArrayGet(&pThis->Xs, j);
        double torsionm2pi = torsion - 2 * PI;
        double torsionp2pi = torsion + 2 * PI;
        double torsion2 = torsion;
        if ((strcmp(ResidueGetName(pThis), "PHE") == 0 && j == 1) || (strcmp(ResidueGetName(pThis), "TYR") == 0 && j == 1) ||
          (strcmp(ResidueGetName(pThis), "ASP") == 0 && j == 1) || strcmp(ResidueGetName(pThis), "GLU") == 0 && j == 2)
        {
          torsion2 = torsion + PI;
          torsion2 = torsion > 0 ? torsion - PI : torsion2;
        }
        double torsion2m2pi = torsion2 - 2 * PI;
        double torsion2p2pi = torsion2 + 2 * PI;
        if (!((torsion <= max && torsion >= min) || (torsionm2pi <= max && torsionm2pi >= min) || (torsionp2pi <= max && torsionp2pi >= min) ||
          (torsion2 <= max && torsion2 >= min) || (torsion2m2pi <= max && torsion2m2pi >= min) || (torsion2p2pi <= max && torsion2p2pi >= min)))
        {
          match = FALSE;
          break;
        }
      }
      if (match)
      {
        matchIdx = i;
        pMatch = p;
        break;
      }
    }
    double pDelta = 1e-7;
    if (matchIdx != -1) ResidueSetDunbrack(pThis, -log(pMatch + pDelta));
    else ResidueSetDunbrack(pThis, -log(pMin + pDelta));
  }

  return Success;
}


int RotamerDunbrackEnergy(Rotamer* pThis, double energyTerms[MAX_ENERGY_TERM])
{
  energyTerms[93] += pThis->dunbrack;
  return Success;
}

int EnergyResidueAndResidueSameChain(Residue* pThis, Residue* pOther, double energyTerms[MAX_ENERGY_TERM])
{
  int neighborCheck = 0;
  if (ResidueGetPosInChain(pThis) + 1 == ResidueGetPosInChain(pOther)) neighborCheck = 12;
  else if (ResidueGetPosInChain(pThis) - 1 == ResidueGetPosInChain(pOther)) neighborCheck = 21;

  for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
  {
    Atom* pAtom1 = ResidueGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE) continue;
    for (int j = 0;j < ResidueGetAtomCount(pOther);j++)
    {
      Atom* pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE) continue;
      if ((neighborCheck == 12 || neighborCheck == 21) && (pAtom1->isBBAtom && pAtom2->isBBAtom))
      {
        continue;
      }
      double distance = XYZDistance(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF) continue;
      int bondType = 15;
      if (neighborCheck == 12) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1), AtomGetName(pAtom2), ResidueGetName(pOther));
      else if (neighborCheck == 21) bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom2), AtomGetName(pAtom1), ResidueGetName(pThis));
      if (bondType == 12 || bondType == 13) continue;
      double att = 0, rep = 0, ele = 0, desP = 0, desH = 0;
      VdwAttEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &att);
      VdwRepEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &rep);
      ElecEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &ele);
      LKDesolvationEnergyAtomAndAtom(pAtom1, pAtom2, distance, bondType, &desP, &desH);
      energyTerms[31] += att;
      energyTerms[32] += rep;
      energyTerms[33] += ele;
      energyTerms[34] += desP;
      energyTerms[35] += desH;
      if (distance < HB_HA_DIST_CUTOFF_MAX)
      {
        double hbtot = 0, hbd = 0, hbt = 0, hbp = 0;
        if (pAtom1->isHBatomH && pAtom2->isHBatomA)
        {
          HBondEnergyAtomAndAtom(pAtom1, pAtom2, ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)), ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        else if (pAtom1->isHBatomA && pAtom2->isHBatomH)
        {
          HBondEnergyAtomAndAtom(pAtom2, pAtom1, ResidueGetAtomByName(pOther, AtomGetHbDorB(pAtom2)), ResidueGetAtomByName(pThis, AtomGetHbDorB(pAtom1)),
            distance, bondType, &hbtot, &hbd, &hbt, &hbp);
        }
        if (pThis->posInChain - pOther->posInChain <= 2 && pThis->posInChain - pOther->posInChain >= -2)
        {
          hbd *= HB_LOCAL_SCALE;
          hbt *= HB_LOCAL_SCALE;
          hbp *= HB_LOCAL_SCALE;
        }
        if (pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == TRUE)
        {
          energyTerms[41] += hbd;
          energyTerms[42] += hbt;
          energyTerms[43] += hbp;
        }
        else if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
        {
          energyTerms[47] += hbd;
          energyTerms[48] += hbt;
          energyTerms[49] += hbp;
        }
        else
        {
          energyTerms[44] += hbd;
          energyTerms[45] += hbt;
          energyTerms[46] += hbp;
        }
      }
    }
  }

  //consider the SSBOND
  if (strcmp("CYS", ResidueGetName(pThis)) == 0 && strcmp("CYS", ResidueGetName(pOther)) == 0)
  {
    Atom* pSG1 = NULL, * pSG2 = NULL, * pCB1 = NULL, * pCB2 = NULL, * pCA1 = NULL, * pCA2 = NULL;
    for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
    {
      Atom* pAtom = ResidueGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA1 = pAtom;
    }
    for (int i = 0;i < ResidueGetAtomCount(pOther);i++)
    {
      Atom* pAtom = ResidueGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0) pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0) pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0) pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDistance(&pSG1->xyz, &pSG2->xyz);
      if (dist<SSBOND_CUTOFF_MAX && dist>SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[36] += ssbond;
      }
    }
  }

  return Success;
}
