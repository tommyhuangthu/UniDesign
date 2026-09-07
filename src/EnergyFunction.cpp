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
//#include <ctype.h>
#include <time.h>

extern double ROTLIB_LOW_PROB_CUTOFF;

double WEIGHTS[MAX_ENERGY_TERM];

int EnergyTermInitialize(double *energyTerms)
{
  memset(energyTerms, 0, sizeof(double) * MAX_ENERGY_TERM);
  return Success;
}

int EnergyTermWeighting(double *energyTerms)
{
  for (int i = 1; i < MAX_ENERGY_TERM; i++)
  {
    energyTerms[i] *= WEIGHTS[i];
    energyTerms[0] += energyTerms[i];
  }
  return Success;
}

int EnergyTermShowMonomer(double *energyTerms)
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
  printf("Total                 =            %12.6f\n", energyTerms[0]);

  return Success;
}

int EnergyTermShowComplex(double *energyTerms)
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
  printf("Total                 =            %12.6f\n", energyTerms[0]);
  return Success;
}

int EnergyWeightRead(char *file)
{
  for (int i = 0; i < MAX_ENERGY_TERM; i++)
  {
    WEIGHTS[i] = 1.0;
  }
  FILE *pIn = fopen(file, "r");
  if (pIn != NULL)
  {
    char line[MAX_LEN_ONE_LINE_CONTENT + 1];
    while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
    {
      char term[MAX_LEN_ONE_LINE_CONTENT + 1];
      double val = 0.0;
      sscanf(line, "%s %lf", term, &val);
      if (strcmp(term, "reference_ALA") == 0)
      {
        WEIGHTS[1] = val;
      }
      else if (strcmp(term, "reference_CYS") == 0)
      {
        WEIGHTS[2] = val;
      }
      else if (strcmp(term, "reference_ASP") == 0)
      {
        WEIGHTS[3] = val;
      }
      else if (strcmp(term, "reference_GLU") == 0)
      {
        WEIGHTS[4] = val;
      }
      else if (strcmp(term, "reference_PHE") == 0)
      {
        WEIGHTS[5] = val;
      }
      else if (strcmp(term, "reference_GLY") == 0)
      {
        WEIGHTS[6] = val;
      }
      else if (strcmp(term, "reference_HIS") == 0)
      {
        WEIGHTS[7] = val;
      }
      else if (strcmp(term, "reference_ILE") == 0)
      {
        WEIGHTS[8] = val;
      }
      else if (strcmp(term, "reference_LYS") == 0)
      {
        WEIGHTS[9] = val;
      }
      else if (strcmp(term, "reference_LEU") == 0)
      {
        WEIGHTS[10] = val;
      }
      else if (strcmp(term, "reference_MET") == 0)
      {
        WEIGHTS[11] = val;
      }
      else if (strcmp(term, "reference_ASN") == 0)
      {
        WEIGHTS[12] = val;
      }
      else if (strcmp(term, "reference_PRO") == 0)
      {
        WEIGHTS[13] = val;
      }
      else if (strcmp(term, "reference_GLN") == 0)
      {
        WEIGHTS[14] = val;
      }
      else if (strcmp(term, "reference_ARG") == 0)
      {
        WEIGHTS[15] = val;
      }
      else if (strcmp(term, "reference_SER") == 0)
      {
        WEIGHTS[16] = val;
      }
      else if (strcmp(term, "reference_THR") == 0)
      {
        WEIGHTS[17] = val;
      }
      else if (strcmp(term, "reference_VAL") == 0)
      {
        WEIGHTS[18] = val;
      }
      else if (strcmp(term, "reference_TRP") == 0)
      {
        WEIGHTS[19] = val;
      }
      else if (strcmp(term, "reference_TYR") == 0)
      {
        WEIGHTS[20] = val;
      }

      else if (strcmp(term, "intraR_vdwatt") == 0)
      {
        WEIGHTS[21] = val;
      }
      else if (strcmp(term, "intraR_vdwrep") == 0)
      {
        WEIGHTS[22] = val;
      }
      else if (strcmp(term, "intraR_electr") == 0)
      {
        WEIGHTS[23] = val;
      }
      else if (strcmp(term, "intraR_deslvP") == 0)
      {
        WEIGHTS[24] = val;
      }
      else if (strcmp(term, "intraR_deslvH") == 0)
      {
        WEIGHTS[25] = val;
      }
      else if (strcmp(term, "intraR_hbscbb_dis") == 0)
      {
        WEIGHTS[26] = val;
      }
      else if (strcmp(term, "intraR_hbscbb_the") == 0)
      {
        WEIGHTS[27] = val;
      }
      else if (strcmp(term, "intraR_hbscbb_phi") == 0)
      {
        WEIGHTS[28] = val;
      }
      else if (strcmp(term, "aapropensity") == 0)
      {
        WEIGHTS[91] = val;
      }
      else if (strcmp(term, "ramachandran") == 0)
      {
        WEIGHTS[92] = val;
      }
      else if (strcmp(term, "dunbrack") == 0)
      {
        WEIGHTS[93] = val;
      }
      else if (strcmp(term, "interS_vdwatt") == 0)
      {
        WEIGHTS[31] = val;
      }
      else if (strcmp(term, "interS_vdwrep") == 0)
      {
        WEIGHTS[32] = val;
      }
      else if (strcmp(term, "interS_electr") == 0)
      {
        WEIGHTS[33] = val;
      }
      else if (strcmp(term, "interS_deslvP") == 0)
      {
        WEIGHTS[34] = val;
      }
      else if (strcmp(term, "interS_deslvH") == 0)
      {
        WEIGHTS[35] = val;
      }
      else if (strcmp(term, "interS_ssbond") == 0)
      {
        WEIGHTS[36] = val;
      }
      else if (strcmp(term, "interS_hbbbbb_dis") == 0)
      {
        WEIGHTS[41] = val;
      }
      else if (strcmp(term, "interS_hbbbbb_the") == 0)
      {
        WEIGHTS[42] = val;
      }
      else if (strcmp(term, "interS_hbbbbb_phi") == 0)
      {
        WEIGHTS[43] = val;
      }
      else if (strcmp(term, "interS_hbscbb_dis") == 0)
      {
        WEIGHTS[44] = val;
      }
      else if (strcmp(term, "interS_hbscbb_the") == 0)
      {
        WEIGHTS[45] = val;
      }
      else if (strcmp(term, "interS_hbscbb_phi") == 0)
      {
        WEIGHTS[46] = val;
      }
      else if (strcmp(term, "interS_hbscsc_dis") == 0)
      {
        WEIGHTS[47] = val;
      }
      else if (strcmp(term, "interS_hbscsc_the") == 0)
      {
        WEIGHTS[48] = val;
      }
      else if (strcmp(term, "interS_hbscsc_phi") == 0)
      {
        WEIGHTS[49] = val;
      }
      else if (strcmp(term, "interD_vdwatt") == 0)
      {
        WEIGHTS[51] = val;
      }
      else if (strcmp(term, "interD_vdwrep") == 0)
      {
        WEIGHTS[52] = val;
      }
      else if (strcmp(term, "interD_electr") == 0)
      {
        WEIGHTS[53] = val;
      }
      else if (strcmp(term, "interD_deslvP") == 0)
      {
        WEIGHTS[54] = val;
      }
      else if (strcmp(term, "interD_deslvH") == 0)
      {
        WEIGHTS[55] = val;
      }
      else if (strcmp(term, "interD_ssbond") == 0)
      {
        WEIGHTS[56] = val;
      }
      else if (strcmp(term, "interD_hbbbbb_dis") == 0)
      {
        WEIGHTS[61] = val;
      }
      else if (strcmp(term, "interD_hbbbbb_the") == 0)
      {
        WEIGHTS[62] = val;
      }
      else if (strcmp(term, "interD_hbbbbb_phi") == 0)
      {
        WEIGHTS[63] = val;
      }
      else if (strcmp(term, "interD_hbscbb_dis") == 0)
      {
        WEIGHTS[64] = val;
      }
      else if (strcmp(term, "interD_hbscbb_the") == 0)
      {
        WEIGHTS[65] = val;
      }
      else if (strcmp(term, "interD_hbscbb_phi") == 0)
      {
        WEIGHTS[66] = val;
      }
      else if (strcmp(term, "interD_hbscsc_dis") == 0)
      {
        WEIGHTS[67] = val;
      }
      else if (strcmp(term, "interD_hbscsc_the") == 0)
      {
        WEIGHTS[68] = val;
      }
      else if (strcmp(term, "interD_hbscsc_phi") == 0)
      {
        WEIGHTS[69] = val;
      }
      else if (strcmp(term, "ligand_vdwatt") == 0)
      {
        WEIGHTS[71] = val;
      }
      else if (strcmp(term, "ligand_vdwrep") == 0)
      {
        WEIGHTS[72] = val;
      }
      else if (strcmp(term, "ligand_electr") == 0)
      {
        WEIGHTS[73] = val;
      }
      else if (strcmp(term, "ligand_deslvP") == 0)
      {
        WEIGHTS[74] = val;
      }
      else if (strcmp(term, "ligand_deslvH") == 0)
      {
        WEIGHTS[75] = val;
      }
      else if (strcmp(term, "ligand_hbscbb_dis") == 0)
      {
        WEIGHTS[84] = val;
      }
      else if (strcmp(term, "ligand_hbscbb_the") == 0)
      {
        WEIGHTS[85] = val;
      }
      else if (strcmp(term, "ligand_hbscbb_phi") == 0)
      {
        WEIGHTS[86] = val;
      }
      else if (strcmp(term, "ligand_hbscsc_dis") == 0)
      {
        WEIGHTS[87] = val;
      }
      else if (strcmp(term, "ligand_hbscsc_the") == 0)
      {
        WEIGHTS[88] = val;
      }
      else if (strcmp(term, "ligand_hbscsc_phi") == 0)
      {
        WEIGHTS[89] = val;
      }
    }
    fclose(pIn);
  }
  else
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];
    sprintf(errMsg, "in file %s line %d, cannot read energy weights from file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, IOError);
    return IOError;
  }

  return Success;
}

int EnergyWeightWrite(char *file)
{
  FILE *pOut = fopen(file, "w");
  if (pOut == NULL)
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];
    sprintf(errMsg, "in file %s line %d, cannot write energy weights to file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, IOError);
    return IOError;
  }
  fprintf(pOut, "reference_ALA      %8.3f\n", WEIGHTS[1]);
  fprintf(pOut, "reference_CYS      %8.3f\n", WEIGHTS[2]);
  fprintf(pOut, "reference_ASP      %8.3f\n", WEIGHTS[3]);
  fprintf(pOut, "reference_GLU      %8.3f\n", WEIGHTS[4]);
  fprintf(pOut, "reference_PHE      %8.3f\n", WEIGHTS[5]);
  fprintf(pOut, "reference_GLY      %8.3f\n", WEIGHTS[6]);
  fprintf(pOut, "reference_HIS      %8.3f\n", WEIGHTS[7]);
  fprintf(pOut, "reference_ILE      %8.3f\n", WEIGHTS[8]);
  fprintf(pOut, "reference_LYS      %8.3f\n", WEIGHTS[9]);
  fprintf(pOut, "reference_LEU      %8.3f\n", WEIGHTS[10]);
  fprintf(pOut, "reference_MET      %8.3f\n", WEIGHTS[11]);
  fprintf(pOut, "reference_ASN      %8.3f\n", WEIGHTS[12]);
  fprintf(pOut, "reference_PRO      %8.3f\n", WEIGHTS[13]);
  fprintf(pOut, "reference_GLN      %8.3f\n", WEIGHTS[14]);
  fprintf(pOut, "reference_ARG      %8.3f\n", WEIGHTS[15]);
  fprintf(pOut, "reference_SER      %8.3f\n", WEIGHTS[16]);
  fprintf(pOut, "reference_THR      %8.3f\n", WEIGHTS[17]);
  fprintf(pOut, "reference_VAL      %8.3f\n", WEIGHTS[18]);
  fprintf(pOut, "reference_TRP      %8.3f\n", WEIGHTS[19]);
  fprintf(pOut, "reference_TYR      %8.3f\n", WEIGHTS[20]);

  fprintf(pOut, "intraR_vdwatt      %8.3f\n", WEIGHTS[21]);
  fprintf(pOut, "intraR_vdwrep      %8.3f\n", WEIGHTS[22]);
  fprintf(pOut, "intraR_electr      %8.3f\n", WEIGHTS[23]);
  fprintf(pOut, "intraR_deslvP      %8.3f\n", WEIGHTS[24]);
  fprintf(pOut, "intraR_deslvH      %8.3f\n", WEIGHTS[25]);
  fprintf(pOut, "intraR_hbscbb_dis  %8.3f\n", WEIGHTS[26]);
  fprintf(pOut, "intraR_hbscbb_the  %8.3f\n", WEIGHTS[27]);
  fprintf(pOut, "intraR_hbscbb_phi  %8.3f\n", WEIGHTS[28]);
  fprintf(pOut, "aapropensity       %8.3f\n", WEIGHTS[91]);
  fprintf(pOut, "ramachandran       %8.3f\n", WEIGHTS[92]);
  fprintf(pOut, "dunbrack           %8.3f\n", WEIGHTS[93]);

  fprintf(pOut, "interS_vdwatt      %8.3f\n", WEIGHTS[31]);
  fprintf(pOut, "interS_vdwrep      %8.3f\n", WEIGHTS[32]);
  fprintf(pOut, "interS_electr      %8.3f\n", WEIGHTS[33]);
  fprintf(pOut, "interS_deslvP      %8.3f\n", WEIGHTS[34]);
  fprintf(pOut, "interS_deslvH      %8.3f\n", WEIGHTS[35]);
  fprintf(pOut, "interS_ssbond      %8.3f\n", WEIGHTS[36]);
  fprintf(pOut, "interS_hbbbbb_dis  %8.3f\n", WEIGHTS[41]);
  fprintf(pOut, "interS_hbbbbb_the  %8.3f\n", WEIGHTS[42]);
  fprintf(pOut, "interS_hbbbbb_phi  %8.3f\n", WEIGHTS[43]);
  fprintf(pOut, "interS_hbscbb_dis  %8.3f\n", WEIGHTS[44]);
  fprintf(pOut, "interS_hbscbb_the  %8.3f\n", WEIGHTS[45]);
  fprintf(pOut, "interS_hbscbb_phi  %8.3f\n", WEIGHTS[46]);
  fprintf(pOut, "interS_hbscsc_dis  %8.3f\n", WEIGHTS[47]);
  fprintf(pOut, "interS_hbscsc_the  %8.3f\n", WEIGHTS[48]);
  fprintf(pOut, "interS_hbscsc_phi  %8.3f\n", WEIGHTS[49]);

  fprintf(pOut, "interD_vdwatt      %8.3f\n", WEIGHTS[51]);
  fprintf(pOut, "interD_vdwrep      %8.3f\n", WEIGHTS[52]);
  fprintf(pOut, "interD_electr      %8.3f\n", WEIGHTS[53]);
  fprintf(pOut, "interD_deslvP      %8.3f\n", WEIGHTS[54]);
  fprintf(pOut, "interD_deslvH      %8.3f\n", WEIGHTS[55]);
  fprintf(pOut, "interD_ssbond      %8.3f\n", WEIGHTS[56]);
  fprintf(pOut, "interD_hbbbbb_dis  %8.3f\n", WEIGHTS[61]);
  fprintf(pOut, "interD_hbbbbb_the  %8.3f\n", WEIGHTS[62]);
  fprintf(pOut, "intraD_hbbbbb_phi  %8.3f\n", WEIGHTS[63]);
  fprintf(pOut, "interD_hbscbb_dis  %8.3f\n", WEIGHTS[64]);
  fprintf(pOut, "interD_hbscbb_the  %8.3f\n", WEIGHTS[65]);
  fprintf(pOut, "interD_hbscbb_phi  %8.3f\n", WEIGHTS[66]);
  fprintf(pOut, "interD_hbscsc_dis  %8.3f\n", WEIGHTS[67]);
  fprintf(pOut, "interD_hbscsc_the  %8.3f\n", WEIGHTS[68]);
  fprintf(pOut, "interD_hbscsc_phi  %8.3f\n", WEIGHTS[69]);

  fprintf(pOut, "ligand_vdwatt      %8.3f\n", WEIGHTS[71]);
  fprintf(pOut, "ligand_vdwrep      %8.3f\n", WEIGHTS[72]);
  fprintf(pOut, "ligand_electr      %8.3f\n", WEIGHTS[73]);
  fprintf(pOut, "ligand_deslvP      %8.3f\n", WEIGHTS[74]);
  fprintf(pOut, "ligand_deslvH      %8.3f\n", WEIGHTS[75]);
  fprintf(pOut, "ligand_hbscbb_dis  %8.3f\n", WEIGHTS[84]);
  fprintf(pOut, "ligand_hbscbb_the  %8.3f\n", WEIGHTS[85]);
  fprintf(pOut, "ligand_hbscbb_phi  %8.3f\n", WEIGHTS[86]);
  fprintf(pOut, "ligand_hbscsc_dis  %8.3f\n", WEIGHTS[87]);
  fprintf(pOut, "ligand_hbscsc_the  %8.3f\n", WEIGHTS[88]);
  fprintf(pOut, "ligand_hbscsc_phi  %8.3f\n", WEIGHTS[89]);
  fclose(pOut);
  return Success;
}

BOOL ResidueIntraBond12Check(char *atom1, char *atom2, Bonds *pBondSet)
{
  for (int i = 0; i < BondsGetNum(pBondSet); i++)
  {
    Bond *pBond = BondsGet(pBondSet, i);
    if ((strcmp(atom1, BondGetFromName(pBond)) == 0 && strcmp(atom2, BondGetToName(pBond)) == 0) ||
        (strcmp(atom2, BondGetFromName(pBond)) == 0 && strcmp(atom1, BondGetToName(pBond)) == 0))
    {
      return TRUE;
    }
  }
  return FALSE;
}

BOOL ResidueIntraBond13Check(char *atom1, char *atom2, Bonds *pBondSet)
{
  for (int i = 0; i < BondsGetNum(pBondSet); i++)
  {
    Bond *pBond = BondsGet(pBondSet, i);
    if (strcmp(atom1, BondGetFromName(pBond)) == 0)
    {
      if (ResidueIntraBond12Check(BondGetToName(pBond), atom2, pBondSet))
      {
        return TRUE;
      }
    }
    else if (strcmp(atom1, BondGetToName(pBond)) == 0)
    {
      if (ResidueIntraBond12Check(BondGetFromName(pBond), atom2, pBondSet))
      {
        return TRUE;
      }
    }
  }
  return FALSE;
}

BOOL ResidueIntraBond14Check(char *atom1, char *atom2, Bonds *pBondSet)
{
  for (int i = 0; i < BondsGetNum(pBondSet); i++)
  {
    Bond *pBond = BondsGet(pBondSet, i);
    if (strcmp(atom1, BondGetFromName(pBond)) == 0)
    {
      if (ResidueIntraBond13Check(BondGetToName(pBond), atom2, pBondSet))
      {
        return TRUE;
      }
    }
    else if (strcmp(atom1, BondGetToName(pBond)) == 0)
    {
      if (ResidueIntraBond13Check(BondGetFromName(pBond), atom2, pBondSet))
      {
        return TRUE;
      }
    }
  }
  return FALSE;
}

int ResidueIntraBondConnectionCheck(char *atom1, char *atom2, Bonds *pBondSet)
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
  return 15;
}

int ResidueAndNextResidueInterBondConnectionCheck_charmm22(char *atomPreRes, char *atomNextRes, Residue *pPreRes, Residue *pNextRes)
{
  if (strcmp(atomPreRes, "C") == 0)
  {
    if (strcmp(atomNextRes, "N") == 0)
    {
      return 12;
    }
    else if (strcmp(atomNextRes, "CA") == 0 ||
             strcmp(atomNextRes, "HN") == 0 ||
             (strcmp(atomNextRes, "CD") == 0 && strcmp(pNextRes->name, "PRO") == 0))
    {
      return 13;
    }
    else if (strcmp(atomNextRes, "CB") == 0 ||
             strcmp(atomNextRes, "C") == 0 ||
             strcmp(atomNextRes, "HA") == 0 ||
             strcmp(atomNextRes, "HA1") == 0 ||
             strcmp(atomNextRes, "HA2") == 0 ||
             (strcmp(atomNextRes, "HD1") == 0 && strcmp(pNextRes->name, "PRO") == 0) ||
             (strcmp(atomNextRes, "HD2") == 0 && strcmp(pNextRes->name, "PRO") == 0) ||
             (strcmp(atomNextRes, "CG") == 0 && strcmp(pNextRes->name, "PRO") == 0))
    {
      return 14;
    }
  }
  else if (strcmp(atomPreRes, "O") == 0 || strcmp(atomPreRes, "CA") == 0)
  {
    if (strcmp(atomNextRes, "N") == 0)
    {
      return 13;
    }
    else if (strcmp(atomNextRes, "CA") == 0 ||
             strcmp(atomNextRes, "HN") == 0 ||
             (strcmp(atomNextRes, "CD") == 0 && strcmp(pNextRes->name, "PRO") == 0))
    {
      return 14;
    }
  }
  else if (strcmp(atomPreRes, "CB") == 0 || strcmp(atomPreRes, "HA") == 0 || strcmp(atomPreRes, "HA1") == 0 || strcmp(atomPreRes, "HA2") == 0 || strcmp(atomPreRes, "N") == 0)
  {
    if (strcmp(atomNextRes, "N") == 0)
    {
      return 14;
    }
  }
  return 15;
}

int ResidueAndNextResidueInterBondConnectionCheck_charmm19(char *atomPreRes, char *atomNextRes, char *nextResName)
{
  if (strcmp(atomPreRes, "C") == 0)
  {
    if (strcmp(atomNextRes, "N") == 0)
    {
      return 12;
    }
    else if (strcmp(atomNextRes, "CA") == 0 ||
             strcmp(atomNextRes, "H") == 0 ||
             (strcmp(atomNextRes, "CD") == 0 && strcmp(nextResName, "PRO") == 0))
    {
      return 13;
    }
    else if (strcmp(atomNextRes, "CB") == 0 ||
             strcmp(atomNextRes, "C") == 0 ||
             (strcmp(atomNextRes, "CG") == 0 && strcmp(nextResName, "PRO") == 0))
    {
      return 14;
    }
  }
  else if (strcmp(atomPreRes, "O") == 0 || strcmp(atomPreRes, "CA") == 0)
  {
    if (strcmp(atomNextRes, "N") == 0)
    {
      return 13;
    }
    else if (strcmp(atomNextRes, "CA") == 0 ||
             strcmp(atomNextRes, "H") == 0 ||
             (strcmp(atomNextRes, "CD") == 0 && strcmp(nextResName, "PRO") == 0))
    {
      return 14;
    }
  }
  else if (strcmp(atomPreRes, "CB") == 0 || strcmp(atomPreRes, "N") == 0)
  {
    if (strcmp(atomNextRes, "N") == 0)
    {
      return 14;
    }
  }

  return 15;
}

int VdwAttEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance, int bondType, double *vdwAtt)
{
  if (distance >= ENERGY_DISTANCE_CUTOFF)
  {
    return Success;
  }
  if (bondType == 12 || bondType == 13)
  {
    return Success;
  }
  // if (AtomIsHydrogen(pAtom1) || AtomIsHydrogen(pAtom2))
  //{ // do not consider H-H vdw attractive energy
  //   return Success;
  // }
  double rsum = RADIUS_SCALE_FOR_VDW * (pAtom1->vdw_rad + pAtom2->vdw_rad);
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
    double epsilon = sqrt(pAtom1->vdw_eps * pAtom2->vdw_eps);
    double B6 = pow(1 / ratio, 6.0);
    double A12 = B6 * B6;
    energy = epsilon * (A12 - 2.0 * B6);
  }
  else if (distance > 5.0 && distance < ENERGY_DISTANCE_CUTOFF)
  {
    double epsilon = sqrt(pAtom1->vdw_eps * pAtom2->vdw_eps);
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
  if (FLAG_DEBUGGING_VDW_ATT == TRUE)
  {
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %f, ratio: %f, vdwAtt: %f\n",
           AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1),
           AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
           bondType, distance, ratio, energy);
  }
  return Success;
}

int VdwRepEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance, int bondType, double *vdwRep)
{
  if (bondType == 12 || bondType == 13)
  {
    return Success;
  }
  double r1 = 0, r2 = 0;
  if (AtomIsHydrogen(pAtom1) == TRUE && (pAtom1->polarity == Type_AtomPolarity_P || pAtom1->polarity == Type_AtomPolarity_C))
  {
    r1 = 0.2 * pAtom1->vdw_rad;
  }
  else
  {
    r1 = RADIUS_SCALE_FOR_VDW * pAtom1->vdw_rad;
  }
  if (AtomIsHydrogen(pAtom2) == TRUE && (pAtom2->polarity == Type_AtomPolarity_P || pAtom2->polarity == Type_AtomPolarity_C))
  {
    r2 = 0.2 * pAtom2->vdw_rad;
  }
  else
  {
    r2 = RADIUS_SCALE_FOR_VDW * pAtom2->vdw_rad;
  }
  // double rsum = RADIUS_SCALE_FOR_VDW * (pAtom1->vdw_radius + pAtom2->vdw_radius);
  double rsum = r1 + r2;
  double ratio = distance / rsum;
  double epsilon = sqrt(pAtom1->vdw_eps * pAtom2->vdw_eps);
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
  // set a cutoff for maximum clash
  // double MAX_CLASH=5.0*epsilon;
  // if(energy>MAX_CLASH)
  //{
  //   energy = MAX_CLASH;
  // }
  *vdwRep = energy;

  if (FLAG_DEBUGGING_VDW_REP)
  {
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, bondType: %2d, dist: %f, ratio: %f, vdwRep: %f\n",
           AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1),
           AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
           bondType, distance, ratio, energy);
  }
  return Success;
}

int HBondEnergyAtomAndAtom(Atom *atomH, Atom *atomA, Atom *atomD, Atom *atomB, double distHA, int bondType, double *etot, double *edist, double *eDHA, double *eHAB)
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
    XYZ xyzDH = XYZDiff(&atomD->xyz, &atomH->xyz);
    XYZ xyzHA = XYZDiff(&atomH->xyz, &atomA->xyz);
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

    if (FLAG_DEBUGGING_HBOND)
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
    XYZ xyzHA = XYZDiff(&atomH->xyz, &atomA->xyz);
    XYZ xyzAAB = XYZDiff(&atomA->xyz, &atomB->xyz);
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

    if (FLAG_DEBUGGING_HBOND)
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

    if (FLAG_DEBUGGING_HBOND)
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
    XYZ xyzDH = XYZDiff(&atomD->xyz, &atomH->xyz);
    XYZ xyzHA = XYZDiff(&atomH->xyz, &atomA->xyz);
    double angleDHA = PI - XYZAngle(&xyzDH, &xyzHA);
    if (RadToDeg(angleDHA) < ANGLE_DHA_CUTOFF_MIN)
    {
      return Success;
    }
    XYZ xyzAAB = XYZDiff(&atomA->xyz, &atomB->xyz);
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

    if (FLAG_DEBUGGING_HBOND)
    {
      printf("AtomH: %1s %4d %4s, AtomA: %1s %4d %4s, dist: %5.2f, theta: %5.1f, phi: %5.1f, edist: %5.2f, ethe: %5.2f, ephi: %5.2f\n",
             AtomGetChainName(atomH), AtomGetPosInChain(atomH), AtomGetName(atomH),
             AtomGetChainName(atomA), AtomGetPosInChain(atomA), AtomGetName(atomA),
             distHA, RadToDeg(angleDHA), RadToDeg(angleHAB), energyDist, energyDHA, energyHAB);
    }
  }

  return Success;
}

int ElecEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double dist12, int bondType, double *elec)
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
  else if (dist12 < 0.8 * (pAtom1->vdw_rad + pAtom2->vdw_rad))
  {
    dist12 = 0.8 * (pAtom1->vdw_rad + pAtom2->vdw_rad);
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

  if (FLAG_DEBUGGING_ELEC)
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

int LKDesolvationEnergyAtomAndAtom(Atom *pAtom1, Atom *pAtom2, double distance, int bondType, double *energyP, double *energyH)
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
  double r1 = pAtom1->vdw_rad * RADIUS_SCALE_FOR_DESOLV;
  double r2 = pAtom2->vdw_rad * RADIUS_SCALE_FOR_DESOLV;
  double r12 = r1 + r2;

  distance = distance < r12 ? r12 : distance;
  double lamda1 = pAtom1->EEF1_lamda * distance * distance;
  double lamda2 = pAtom2->EEF1_lamda * distance * distance;
  double x1 = (distance - r1) / pAtom1->EEF1_lamda;
  double x2 = (distance - r2) / pAtom2->EEF1_lamda;

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

  if (FLAG_DEBUGGING_DESOLV)
  {
    printf("Atom1: %1s %4d %4s, Atom2: %1s %4d %4s, dist: %5.2f, desolv12: %7.3f, desolv21: %7.3f\n",
           AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), AtomGetName(pAtom1),
           AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), AtomGetName(pAtom2),
           distance, desolv12, desolv21);
  }

  return Success;
}

///////////////////////////////////////////////////////////////////
// User can define new energy term here
///////////////////////////////////////////////////////////////////

int SSbondEnergyAtomAndAtom(Atom *pAtomS1, Atom *pAtomS2, Atom *pAtomCB1, Atom *pAtomCB2, Atom *pAtomCA1, Atom *pAtomCA2, double *sse)
{
  XYZ xyzC1S1 = XYZDiff(&pAtomCB1->xyz, &pAtomS1->xyz);
  XYZ xyzS1S2 = XYZDiff(&pAtomS1->xyz, &pAtomS2->xyz);
  XYZ xyzS2C2 = XYZDiff(&pAtomS2->xyz, &pAtomCB2->xyz);
  double dSS = XYZDist(&pAtomS1->xyz, &pAtomS2->xyz);
  double aC1S1S2 = RadToDeg(PI - XYZAngle(&xyzS1S2, &xyzS2C2));
  double aC2S2S1 = RadToDeg(PI - XYZAngle(&xyzC1S1, &xyzS1S2));
  double xC1S1S2C2 = GetTorsionAngle(&pAtomCB1->xyz, &pAtomS1->xyz, &pAtomS2->xyz, &pAtomCB2->xyz);
  double xCA1CB1SG1SG2 = GetTorsionAngle(&pAtomCA1->xyz, &pAtomCB1->xyz, &pAtomS1->xyz, &pAtomS2->xyz);
  double xCA2CB2SG2SG1 = GetTorsionAngle(&pAtomCA2->xyz, &pAtomCB2->xyz, &pAtomS2->xyz, &pAtomS1->xyz);
  *sse = (0.8 * (1 - exp(-10.0 * (dSS - SSBOND_DISTANCE))) * (1 - exp(-10.0 * (dSS - SSBOND_DISTANCE)))) + 0.005 * (aC1S1S2 - SSBOND_ANGLE) * (aC1S1S2 - SSBOND_ANGLE) + 0.005 * (aC2S2S1 - SSBOND_ANGLE) * (aC2S2S1 - SSBOND_ANGLE) + cos(2.0 * xC1S1S2C2) + 1.0 + 1.25 * sin(xCA1CB1SG1SG2 + 2.0 * PI / 3.0) - 1.75 + 1.25 * sin(xCA2CB2SG2SG1 + 2.0 * PI / 3.0) - 1.75;
  if (*sse > 0.0)
  {
    *sse = 0.0;
  }
  return Success;
}

#define RESIDUE_PAIRWISE_ENERGY_COMPUTATION
int AminoAcidReferenceEnergy(char *AAname, double *energyTerms)
{
  if (strcmp(AAname, "ALA") == 0)
  {
    energyTerms[1] += 1.0;
  }
  else if (strcmp(AAname, "CYS") == 0)
  {
    energyTerms[2] += 1.0;
  }
  else if (strcmp(AAname, "ASP") == 0)
  {
    energyTerms[3] += 1.0;
  }
  else if (strcmp(AAname, "GLU") == 0)
  {
    energyTerms[4] += 1.0;
  }
  else if (strcmp(AAname, "PHE") == 0)
  {
    energyTerms[5] += 1.0;
  }
  else if (strcmp(AAname, "GLY") == 0)
  {
    energyTerms[6] += 1.0;
  }
  else if (strcmp(AAname, "HIS") == 0)
  {
    energyTerms[7] += 1.0;
  }
  else if (strcmp(AAname, "HSE") == 0)
  {
    energyTerms[7] += 1.0;
  }
  else if (strcmp(AAname, "HSD") == 0)
  {
    energyTerms[7] += 1.0;
  }
  else if (strcmp(AAname, "ILE") == 0)
  {
    energyTerms[8] += 1.0;
  }
  else if (strcmp(AAname, "LYS") == 0)
  {
    energyTerms[9] += 1.0;
  }
  else if (strcmp(AAname, "LEU") == 0)
  {
    energyTerms[10] += 1.0;
  }
  else if (strcmp(AAname, "MET") == 0)
  {
    energyTerms[11] += 1.0;
  }
  else if (strcmp(AAname, "ASN") == 0)
  {
    energyTerms[12] += 1.0;
  }
  else if (strcmp(AAname, "PRO") == 0)
  {
    energyTerms[13] += 1.0;
  }
  else if (strcmp(AAname, "GLN") == 0)
  {
    energyTerms[14] += 1.0;
  }
  else if (strcmp(AAname, "ARG") == 0)
  {
    energyTerms[15] += 1.0;
  }
  else if (strcmp(AAname, "SER") == 0)
  {
    energyTerms[16] += 1.0;
  }
  else if (strcmp(AAname, "THR") == 0)
  {
    energyTerms[17] += 1.0;
  }
  else if (strcmp(AAname, "VAL") == 0)
  {
    energyTerms[18] += 1.0;
  }
  else if (strcmp(AAname, "TRP") == 0)
  {
    energyTerms[19] += 1.0;
  }
  else if (strcmp(AAname, "TYR") == 0)
  {
    energyTerms[20] += 1.0;
  }
  return Success;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// the following are energy between residue and residue
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int EnergyIntraResidue(Residue *pThis, double *energyTerms)
{
  for (int i = 0; i < ResidueGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = ResidueGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
    {
      continue;
    }
    for (int j = i + 1; j < ResidueGetAtomNum(pThis); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pThis, j);
      if (pAtom2->isXyzValid == FALSE)
      {
        continue;
      }
      if (pAtom2->isBBAtom == TRUE && pAtom1->isBBAtom == TRUE)
      {
        continue;
      }
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
      {
        continue;
      }
      if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
      {
        if (strcmp(ResidueGetName(pThis), "ILE") == 0 || strcmp(ResidueGetName(pThis), "MET") == 0 || strcmp(ResidueGetName(pThis), "GLN") == 0 || strcmp(ResidueGetName(pThis), "GLU") == 0 || strcmp(ResidueGetName(pThis), "LYS") == 0 || strcmp(ResidueGetName(pThis), "ARG") == 0)
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

int EnergyResidueAndNextResidue(Residue *pThis, Residue *pOther, double *energyTerms)
{
  for (int i = 0; i < ResidueGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = ResidueGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
    {
      continue;
    }
    for (int j = 0; j < ResidueGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)
      {
        continue;
      }
      if (pAtom1->isBBAtom && pAtom2->isBBAtom)
      {
        continue;
      }
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
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
      else if ((pAtom1->isBBAtom == TRUE && pAtom2->isBBAtom == FALSE) || (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == TRUE))
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

int EnergyResidueAndOtherResidueSameChain(Residue *pThis, Residue *pOther, double *energyTerms)
{
  for (int i = 0; i < ResidueGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = ResidueGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    for (int j = 0; j < ResidueGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
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
    Atom *pSG1 = NULL, *pSG2 = NULL, *pCB1 = NULL, *pCB2 = NULL, *pCA1 = NULL, *pCA2 = NULL;
    for (int i = 0; i < ResidueGetAtomNum(pThis); i++)
    {
      Atom *pAtom = ResidueGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA1 = pAtom;
    }
    for (int i = 0; i < ResidueGetAtomNum(pOther); i++)
    {
      Atom *pAtom = ResidueGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDist(&pSG1->xyz, &pSG2->xyz);
      if (dist < SSBOND_CUTOFF_MAX && dist > SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[36] += ssbond;
      }
    }
  }

  return Success;
}

int EnergyResidueAndOtherResidueDiffChain(Residue *pThis, Residue *pOther, double *energyTerms)
{
  for (int i = 0; i < ResidueGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = ResidueGetAtom(pThis, i);
    if (!pAtom1->isXyzValid)
      continue;
    for (int j = 0; j < ResidueGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pOther, j);
      if (!pAtom2->isXyzValid)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
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

  // consider the SSBOND
  if (strcmp("CYS", ResidueGetName(pThis)) == 0 && strcmp("CYS", ResidueGetName(pOther)) == 0)
  {
    Atom *pSG1 = NULL, *pSG2 = NULL, *pCB1 = NULL, *pCB2 = NULL, *pCA1 = NULL, *pCA2 = NULL;
    for (int i = 0; i < ResidueGetAtomNum(pThis); i++)
    {
      Atom *pAtom = ResidueGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA1 = pAtom;
    }
    for (int i = 0; i < ResidueGetAtomNum(pOther); i++)
    {
      Atom *pAtom = ResidueGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDist(&pSG1->xyz, &pSG2->xyz);
      if (dist < SSBOND_CUTOFF_MAX && dist > SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[56] += ssbond;
      }
    }
  }

  return Success;
}

int EnergyResidueAndLigandResidue(Residue *pProtein, Residue *pLigand, double *energyTerms)
{
  for (int i = 0; i < ResidueGetAtomNum(pProtein); i++)
  {
    Atom *pAtom1 = ResidueGetAtom(pProtein, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    for (int j = 0; j < ResidueGetAtomNum(pLigand); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pLigand, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
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

int EnergyIntraRotamer(Rotamer *pThis, double *energyTerms)
{
  for (int i = 0; i < RotamerGetAtomNum(pThis); ++i)
  {
    Atom *pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    for (int j = i + 1; j < RotamerGetAtomNum(pThis); ++j)
    {
      Atom *pAtom2 = RotamerGetAtom(pThis, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      if (pAtom2->isBBAtom && pAtom1->isBBAtom)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
      if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
      {
        if (strcmp(RotamerGetType(pThis), "ILE") == 0 || strcmp(RotamerGetType(pThis), "MET") == 0 ||
            strcmp(RotamerGetType(pThis), "GLN") == 0 || strcmp(RotamerGetType(pThis), "GLU") == 0 ||
            strcmp(RotamerGetType(pThis), "LYS") == 0 || strcmp(RotamerGetType(pThis), "ARG") == 0)
        {
          int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom1), AtomGetName(pAtom2), RotamerGetBonds(pThis));
          if (bondType == 12 || bondType == 13)
            continue;
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
          continue;
        int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom1), AtomGetName(pAtom2), RotamerGetBonds(pThis));
        if (bondType == 12 || bondType == 13)
          continue;
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

int EnergyRotamerAndRotamerSameChain(Rotamer *pThis, Rotamer *pOther, double *energyTerms)
{
  for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    if (pAtom1->isBBAtom)
      continue;
    for (int j = 0; j < RotamerGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = RotamerGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      if (pAtom2->isBBAtom)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
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

  // consider the SSBOND
  if (strcmp("CYS", RotamerGetType(pThis)) == 0 && strcmp("CYS", RotamerGetType(pOther)) == 0)
  {
    Atom *pSG1 = NULL, *pSG2 = NULL, *pCB1 = NULL, *pCB2 = NULL, *pCA1 = NULL, *pCA2 = NULL;
    for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
    {
      Atom *pAtom = RotamerGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA1 = pAtom;
    }
    for (int i = 0; i < RotamerGetAtomNum(pOther); i++)
    {
      Atom *pAtom = RotamerGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDist(&pSG1->xyz, &pSG2->xyz);
      if (dist < SSBOND_CUTOFF_MAX && dist > SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[36] += ssbond;
      }
    }
  }

  return Success;
}

int EnergyRotamerAndRotamerDiffChain(Rotamer *pThis, Rotamer *pOther, double *energyTerms)
{
  for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    if (pAtom1->isBBAtom)
      continue;
    for (int j = 0; j < RotamerGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = RotamerGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      if (pAtom2->isBBAtom)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
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

  // consider the SSBOND
  if (strcmp("CYS", RotamerGetType(pThis)) == 0 && strcmp("CYS", RotamerGetType(pOther)) == 0)
  {
    Atom *pSG1 = NULL, *pSG2 = NULL, *pCB1 = NULL, *pCB2 = NULL, *pCA1 = NULL, *pCA2 = NULL;
    for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
    {
      Atom *pAtom = RotamerGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA1 = pAtom;
    }
    for (int i = 0; i < RotamerGetAtomNum(pOther); i++)
    {
      Atom *pAtom = RotamerGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDist(&pSG1->xyz, &pSG2->xyz);
      if (dist < SSBOND_CUTOFF_MAX && dist > SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[56] += ssbond;
      }
    }
  }
  return Success;
}

int EnergyRotamerAndLigRotamer(Rotamer *pThis, Rotamer *pOther, double *energyTerms)
{
  for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    if (pAtom1->isBBAtom)
      continue;
    for (int j = 0; j < RotamerGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = RotamerGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      if (pAtom2->isBBAtom)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
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

int EnergyRotamerAndFixedResidueSameChain(Rotamer *pThis, Residue *pOther, double *energyTerms)
{
  int neighborCheck = 0;
  if (RotamerGetPosInChain(pThis) + 1 == ResidueGetPosInChain(pOther))
    neighborCheck = 12;
  else if (RotamerGetPosInChain(pThis) - 1 == ResidueGetPosInChain(pOther))
    neighborCheck = 21;

  for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    if (pAtom1->isBBAtom)
      continue;
    for (int j = 0; j < ResidueGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
      int bondType = 15;
      if (neighborCheck == 12)
        bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1), AtomGetName(pAtom2), ResidueGetName(pOther));
      else if (neighborCheck == 21)
        bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom2), AtomGetName(pAtom1), RotamerGetType(pThis));
      if (bondType == 12 || bondType == 13)
        continue;
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
        else */
        if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
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
  if (strcmp("CYS", RotamerGetType(pThis)) == 0 && strcmp("CYS", ResidueGetName(pOther)) == 0)
  {
    Atom *pSG1 = NULL, *pSG2 = NULL, *pCB1 = NULL, *pCB2 = NULL, *pCA1 = NULL, *pCA2 = NULL;
    for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
    {
      Atom *pAtom = RotamerGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA1 = pAtom;
    }
    for (int i = 0; i < ResidueGetAtomNum(pOther); i++)
    {
      Atom *pAtom = ResidueGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDist(&pSG1->xyz, &pSG2->xyz);
      if (dist < SSBOND_CUTOFF_MAX && dist > SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[36] += ssbond;
      }
    }
  }

  return Success;
}

int EnergyRotamerAndFixedResidueDiffChain(Rotamer *pThis, Residue *pOther, double *energyTerms)
{
  for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = RotamerGetAtom(pThis, i);
    if (!pAtom1->isXyzValid)
      continue;
    if (pAtom1->isBBAtom)
      continue;
    for (int j = 0; j < ResidueGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pOther, j);
      if (!pAtom2->isXyzValid)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
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
        else */
        if (pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE)
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

  // consider the SSBOND
  if (strcmp("CYS", RotamerGetType(pThis)) == 0 && strcmp("CYS", ResidueGetName(pOther)) == 0)
  {
    Atom *pSG1 = NULL, *pSG2 = NULL, *pCB1 = NULL, *pCB2 = NULL, *pCA1 = NULL, *pCA2 = NULL;
    for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
    {
      Atom *pAtom = RotamerGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA1 = pAtom;
    }
    for (int i = 0; i < ResidueGetAtomNum(pOther); i++)
    {
      Atom *pAtom = ResidueGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDist(&pSG1->xyz, &pSG2->xyz);
      if (dist < SSBOND_CUTOFF_MAX && dist > SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[56] += ssbond;
      }
    }
  }

  return Success;
}

int EnergyRotamerAndFixedLigResidue(Rotamer *pThis, Residue *pLigand, double *energyTerms)
{
  for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    if (pAtom1->isBBAtom)
      continue;
    for (int j = 0; j < ResidueGetAtomNum(pLigand); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pLigand, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      // if(pAtom2->isBBAtom) continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
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
        // if(pAtom1->isBBAtom == FALSE && pAtom2->isBBAtom == FALSE){
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

int EnergyLigRotamerAndFixedResidue(Rotamer *pThis, Residue *pOther, double *energyTerms)
{
  for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    for (int j = 0; j < ResidueGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
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

int EnergyRotamerAndDesignResidueSameChain(Rotamer *pThis, Residue *pOther, double *energyTerms)
{
  int neighborCheck = 0;
  if (RotamerGetPosInChain(pThis) + 1 == ResidueGetPosInChain(pOther))
    neighborCheck = 12;
  else if (RotamerGetPosInChain(pThis) - 1 == ResidueGetPosInChain(pOther))
    neighborCheck = 21;

  for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    if (pAtom1->isBBAtom)
      continue;
    for (int j = 0; j < ResidueGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      if (pAtom2->isBBAtom == FALSE)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
      int bondType = 15;
      if (neighborCheck == 12)
        bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1), AtomGetName(pAtom2), ResidueGetName(pOther));
      else if (neighborCheck == 21)
        bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom2), AtomGetName(pAtom1), RotamerGetType(pThis));
      if (bondType == 12 || bondType == 13)
        continue;
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

int EnergyRotamerAndDesignResidueDiffChain(Rotamer *pThis, Residue *pOther, double *energyTerms)
{
  for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = RotamerGetAtom(pThis, i);
    if (!pAtom1->isXyzValid)
      continue;
    if (pAtom1->isBBAtom)
      continue;
    for (int j = 0; j < ResidueGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pOther, j);
      if (!pAtom2->isXyzValid)
        continue;
      if (!pAtom2->isBBAtom)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
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

int EnergyLigRotamerAndDesignResidue(Rotamer *pThis, Residue *pOther, double *energyTerms)
{
  for (int i = 0; i < RotamerGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = RotamerGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    if (pAtom1->isBBAtom)
      continue;
    for (int j = 0; j < ResidueGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      if (pAtom2->isBBAtom == FALSE)
        continue;
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
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

int RamaTableReadFromFile(RamaTable *pRama, char *filename)
{
  FILE *pFile = fopen(filename, "r");
  if (!pFile)
    return IOError;
  char buffer[MAX_LEN_ONE_LINE_CONTENT + 1];
  while (fgets(buffer, MAX_LEN_ONE_LINE_CONTENT, pFile))
  {
    if (buffer[0] == ' ' || buffer[0] == '#')
      continue;
    int phipsi[2];
    double aae[20];
    sscanf(buffer, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
           &phipsi[0], &phipsi[1], &aae[0], &aae[1], &aae[2], &aae[3], &aae[4], &aae[5], &aae[6], &aae[7], &aae[8], &aae[9],
           &aae[10], &aae[11], &aae[12], &aae[13], &aae[14], &aae[15], &aae[16], &aae[17], &aae[18], &aae[19]);
    int phiindex = (phipsi[0] + 180) / 10;
    int psiindex = (phipsi[1] + 180) / 10;
    for (int j = 0; j < 20; j++)
    {
      pRama->ramatable[phiindex][psiindex][j] = aae[j];
    }
  }
  fclose(pFile);
  return Success;
}

int AApropensityTableReadFromFile(AAppTable *pAAPP, char *filename)
{
  FILE *pFile = fopen(filename, "r");
  if (pFile == NULL)
  {
    return IOError;
  }
  char buffer[MAX_LEN_ONE_LINE_CONTENT + 1];
  while (fgets(buffer, MAX_LEN_ONE_LINE_CONTENT, pFile))
  {
    if (buffer[0] == ' ' || buffer[0] == '#')
    {
      continue;
    }
    int phipsi[2];
    double aae[20];
    sscanf(buffer, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
           &phipsi[0], &phipsi[1], &aae[0], &aae[1], &aae[2], &aae[3], &aae[4], &aae[5], &aae[6], &aae[7], &aae[8], &aae[9],
           &aae[10], &aae[11], &aae[12], &aae[13], &aae[14], &aae[15], &aae[16], &aae[17], &aae[18], &aae[19]);
    int phiindex = (phipsi[0] + 180) / 10;
    int psiindex = (phipsi[1] + 180) / 10;
    for (int j = 0; j < 20; j++)
    {
      pAAPP->aapptable[phiindex][psiindex][j] = aae[j];
    }
  }
  fclose(pFile);
  return Success;
}

int AminoAcidPropensityAndRamachandranEnergy(Residue *pRes, AAppTable *pAAPP, RamaTable *pRama)
{
  char aa1 = AA3ToAA1(ResidueGetName(pRes));
  int aaNdx = AA1GetIndex(aa1);
  int phiNdx = ((int)pRes->phipsi[0] + 180) / 10;
  int psiNdx = ((int)pRes->phipsi[1] + 180) / 10;
  pRes->aapp += pAAPP->aapptable[phiNdx][psiNdx][aaNdx];
  pRes->rama += pRama->ramatable[phiNdx][psiNdx][aaNdx];
  return Success;
}

int RotamerPropensityAndRamachandranEnergy(Rotamer *pRot, Residue *pRes, AAppTable *pAAPP, RamaTable *pRama, double *energyTerms)
{
  char aa1 = AA3ToAA1(RotamerGetType(pRot));
  int aaNdx = AA1GetIndex(aa1);
  int phiNdx = ((int)pRes->phipsi[0] + 180) / 10;
  int psiNdx = ((int)pRes->phipsi[1] + 180) / 10;
  energyTerms[91] += pAAPP->aapptable[phiNdx][psiNdx][aaNdx];
  energyTerms[92] += pRama->ramatable[phiNdx][psiNdx][aaNdx];
  return Success;
}

int AminoAcidDunbrackEnergy(Residue *pRes, BBdepRotLib *rotlib)
{
  if (strcmp(ResidueGetName(pRes), "ALA") == 0 ||
      strcmp(ResidueGetName(pRes), "GLY") == 0 ||
      AA3GetIndex(ResidueGetName(pRes)) < 0)
  {
    return Success;
  }

  if (pRes->isSCIntact == TRUE)
  {
    int binIndex = ((int)(pRes->phipsi[0] + 180) / 10) * 36 + (int)(pRes->phipsi[1] + 180) / 10;
    RotLibPhiPsi *pRotLibPhiPsi = &rotlib->pps[binIndex];
    int rotIndex = -1;
    StringArrayFind(&pRotLibPhiPsi->rotTypes, ResidueGetName(pRes), &rotIndex);
    DoubleArray *pTorsionsArray = pRotLibPhiPsi->tors[rotIndex];
    DoubleArray *pDeviationsArray = pRotLibPhiPsi->devs[rotIndex];
    int matchedIndex = -1;
    double matchedProbability;
    double accumulativeProbability = 0;
    for (int i = 0; i < IntArrayGet(&pRotLibPhiPsi->rotNums, rotIndex); i++)
    {
      double prob = DoubleArrayGet(&pRotLibPhiPsi->prob[rotIndex], i);
      accumulativeProbability += prob;

      DoubleArray *pTorsions = &pTorsionsArray[i];
      DoubleArray *pDeviations = &pDeviationsArray[i];
      BOOL rotamerMatched = TRUE;
      for (int j = 0; j < DoubleArrayGetLength(&pRes->Xs); j++)
      {
        double Xmin = DoubleArrayGet(pTorsions, j) - DoubleArrayGet(pDeviations, j);
        double Xmax = DoubleArrayGet(pTorsions, j) + DoubleArrayGet(pDeviations, j);
        double X = DoubleArrayGet(&pRes->Xs, j);
        double Xm2pi = X - 2 * PI;
        double Xp2pi = X + 2 * PI;
        double X2 = X;
        if ((strcmp(ResidueGetName(pRes), "PHE") == 0 && j == 1) ||
            (strcmp(ResidueGetName(pRes), "TYR") == 0 && j == 1) ||
            (strcmp(ResidueGetName(pRes), "ASP") == 0 && j == 1) ||
            (strcmp(ResidueGetName(pRes), "GLU") == 0 && j == 2))
        {
          X2 = X + PI;
          X2 = X > 0 ? X - PI : X2;
        }
        double X2m2pi = X2 - 2 * PI;
        double X2p2pi = X2 + 2 * PI;
        if (!((X <= Xmax && X >= Xmin) ||
              (Xm2pi <= Xmax && Xm2pi >= Xmin) ||
              (Xp2pi <= Xmax && Xp2pi >= Xmin) ||
              (X2 <= Xmax && X2 >= Xmin) ||
              (X2m2pi <= Xmax && X2m2pi >= Xmin) ||
              (X2p2pi <= Xmax && X2p2pi >= Xmin)))
        {
          rotamerMatched = FALSE;
          break;
        }
      }
      if (rotamerMatched == TRUE)
      {
        matchedIndex = i;
        matchedProbability = prob;
        break;
      }
    }

    if (matchedIndex != -1)
    {
      ResidueSetDunbrack(pRes, -log(matchedProbability));
    }
    else
    {
      double unmatchedProbability = 1 - accumulativeProbability;
      if (unmatchedProbability < 1e-3)
        unmatchedProbability = 1e-3;
      else if (unmatchedProbability > ROTLIB_LOW_PROB_CUTOFF)
        unmatchedProbability = ROTLIB_LOW_PROB_CUTOFF;
      printf("sidechain of native residue %s%d%s has no matched rotamer conformation in library; "
        "its rotamer frequency was set to an unmatchedProbability of %f\n",
        ResidueGetChainName(pRes), ResidueGetPosInChain(pRes), ResidueGetName(pRes), unmatchedProbability);
      ResidueSetDunbrack(pRes, -log(unmatchedProbability));
    }
  }

  return Success;
}

int RotamerDunbrackEnergy(Rotamer *pThis, double *energyTerms)
{
  energyTerms[93] += pThis->dunbrack;
  return Success;
}

int EnergyResidueAndResidueSameChain(Residue *pThis, Residue *pOther, double *energyTerms)
{
  int neighborCheck = 0;
  if (ResidueGetPosInChain(pThis) + 1 == ResidueGetPosInChain(pOther))
    neighborCheck = 12;
  else if (ResidueGetPosInChain(pThis) - 1 == ResidueGetPosInChain(pOther))
    neighborCheck = 21;

  for (int i = 0; i < ResidueGetAtomNum(pThis); i++)
  {
    Atom *pAtom1 = ResidueGetAtom(pThis, i);
    if (pAtom1->isXyzValid == FALSE)
      continue;
    for (int j = 0; j < ResidueGetAtomNum(pOther); j++)
    {
      Atom *pAtom2 = ResidueGetAtom(pOther, j);
      if (pAtom2->isXyzValid == FALSE)
        continue;
      if ((neighborCheck == 12 || neighborCheck == 21) && (pAtom1->isBBAtom && pAtom2->isBBAtom))
      {
        continue;
      }
      double distance = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
      if (distance > ENERGY_DISTANCE_CUTOFF)
        continue;
      int bondType = 15;
      if (neighborCheck == 12)
        bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom1), AtomGetName(pAtom2), ResidueGetName(pOther));
      else if (neighborCheck == 21)
        bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(AtomGetName(pAtom2), AtomGetName(pAtom1), ResidueGetName(pThis));
      if (bondType == 12 || bondType == 13)
        continue;
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

  // consider the SSBOND
  if (strcmp("CYS", ResidueGetName(pThis)) == 0 && strcmp("CYS", ResidueGetName(pOther)) == 0)
  {
    Atom *pSG1 = NULL, *pSG2 = NULL, *pCB1 = NULL, *pCB2 = NULL, *pCA1 = NULL, *pCA2 = NULL;
    for (int i = 0; i < ResidueGetAtomNum(pThis); i++)
    {
      Atom *pAtom = ResidueGetAtom(pThis, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB1 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA1 = pAtom;
    }
    for (int i = 0; i < ResidueGetAtomNum(pOther); i++)
    {
      Atom *pAtom = ResidueGetAtom(pOther, i);
      if (strcmp(AtomGetName(pAtom), "SG") == 0)
        pSG2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CB") == 0)
        pCB2 = pAtom;
      else if (strcmp(AtomGetName(pAtom), "CA") == 0)
        pCA2 = pAtom;
    }
    if (pSG1 != NULL && pCB1 != NULL && pSG2 != NULL && pCB2 != NULL && pCA1 != NULL && pCA2 != NULL)
    {
      double dist = XYZDist(&pSG1->xyz, &pSG2->xyz);
      if (dist < SSBOND_CUTOFF_MAX && dist > SSBOND_CUTOFF_MIN)
      {
        double ssbond = 0;
        SSbondEnergyAtomAndAtom(pSG1, pSG2, pCB1, pCB2, pCA1, pCA2, &ssbond);
        energyTerms[36] += ssbond;
      }
    }
  }

  return Success;
}
