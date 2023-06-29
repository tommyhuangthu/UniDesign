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
#pragma warning(disable:6031)

#include "Utility.h"
#include "Atom.h"
#include "Residue.h"
#include "SmallMolEEF1.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "SmallMolParAndTopo.h"


int FindAtomC(Atom* pAtomD, AtomArray* pAtomArray, IntArray* icFlags, IntArray* bondsFromToAndType, AtomArray* pAtomCArray)
{
  AtomArrayDestroy(pAtomCArray);
  AtomArrayCreate(pAtomCArray);
  for (int i = 0; i < IntArrayGetLength(bondsFromToAndType); i += 2)
  {
    int atomIndex = IntArrayGet(bondsFromToAndType, i);
    Atom* pAtom = AtomArrayGet(pAtomArray, atomIndex);
    if (strcmp(AtomGetName(pAtomD), AtomGetName(pAtom)) == 0)
    {
      int atomCIndex = IntArrayGet(bondsFromToAndType, i + 1);
      if (IntArrayGet(icFlags, atomCIndex) == 0) continue;
      AtomArrayAppend(pAtomCArray, AtomArrayGet(pAtomArray, atomCIndex));
    }

    atomIndex = IntArrayGet(bondsFromToAndType, i + 1);
    pAtom = AtomArrayGet(pAtomArray, atomIndex);
    if (strcmp(AtomGetName(pAtomD), AtomGetName(pAtom)) == 0)
    {
      int atomCIndex = IntArrayGet(bondsFromToAndType, i);
      if (IntArrayGet(icFlags, atomCIndex) == 0) continue;
      AtomArrayAppend(pAtomCArray, AtomArrayGet(pAtomArray, atomCIndex));
    }
  }

  return Success;
}


int FindAtomB(Atom* pAtomC, Atom* pAtomD, AtomArray* pAtoms, IntArray* icFlags, IntArray* bondsFromToAndType, AtomArray* pAtomBArray)
{
  AtomArrayDestroy(pAtomBArray);
  AtomArrayCreate(pAtomBArray);
  for (int i = 0; i < IntArrayGetLength(bondsFromToAndType); i += 2)
  {
    int atomIndex = IntArrayGet(bondsFromToAndType, i);
    Atom* pAtom = AtomArrayGet(pAtoms, atomIndex);
    if (strcmp(AtomGetName(pAtomC), AtomGetName(pAtom)) == 0)
    {
      int atomCIndex = IntArrayGet(bondsFromToAndType, i + 1);
      if (IntArrayGet(icFlags, atomCIndex) == 0) continue;
      if (strcmp(AtomGetName(pAtomD), AtomGetName(AtomArrayGet(pAtoms, atomCIndex))) == 0) continue;
      AtomArrayAppend(pAtomBArray, AtomArrayGet(pAtoms, atomCIndex));
    }

    atomIndex = IntArrayGet(bondsFromToAndType, i + 1);
    pAtom = AtomArrayGet(pAtoms, atomIndex);
    if (strcmp(AtomGetName(pAtomC), AtomGetName(pAtom)) == 0)
    {
      int atomCIndex = IntArrayGet(bondsFromToAndType, i);
      if (IntArrayGet(icFlags, atomCIndex) == 0) continue;
      if (strcmp(AtomGetName(pAtomD), AtomGetName(AtomArrayGet(pAtoms, atomCIndex))) == 0) continue;
      AtomArrayAppend(pAtomBArray, AtomArrayGet(pAtoms, atomCIndex));
    }
  }

  return Success;
}


int GetAttachedAtomTypeNum(AtomArray* pAtoms, int bonds[][2], int bondNum, int atomIndex, int attachAtomTypeNum[10], int* attachCount, int attachIndex[10])
{
  // attach[ 0]:C attached
  // attach[ 1]:H
  // attach[ 2]:O
  // attach[ 3]:N
  // attach[ 4]:P
  // attach[ 5]:S
  // attach[ 6]:F
  // attach[ 7]:Cl
  // attach[ 8]:Br
  // attach[ 9]:I
  // attach[10]:B
  // attach[11]:Fe/Zn/X,metal
  for (int i = 0;i < bondNum;i++)
  {
    if (atomIndex == bonds[i][0] - 1)
    {
      attachIndex[*attachCount] = bonds[i][1] - 1;
      (*attachCount)++;
      Atom* pAtom = AtomArrayGet(pAtoms, bonds[i][1] - 1);
      if (pAtom->name[0] == 'C') attachAtomTypeNum[0]++;
      else if (pAtom->name[0] == 'H')attachAtomTypeNum[1]++;
      else if (pAtom->name[0] == 'O')attachAtomTypeNum[2]++;
      else if (pAtom->name[0] == 'N')attachAtomTypeNum[3]++;
      else attachAtomTypeNum[4]++;
    }
    else if (atomIndex == bonds[i][1] - 1)
    {
      attachIndex[*attachCount] = bonds[i][0] - 1;
      (*attachCount)++;
      Atom* pAtom = AtomArrayGet(pAtoms, bonds[i][0] - 1);
      if (pAtom->name[0] == 'C') attachAtomTypeNum[0]++;
      else if (pAtom->name[0] == 'H')attachAtomTypeNum[1]++;
      else if (pAtom->name[0] == 'O')attachAtomTypeNum[2]++;
      else if (pAtom->name[0] == 'N')attachAtomTypeNum[3]++;
      else attachAtomTypeNum[4]++;
    }
  }

  return Success;
}


int GenerateSmallMolParameterAndTopologyFromMol2(char* mol2file, char* parfile, char* topfile, char* iniAtom1, char* iniAtom2, char* iniAtom3)
{
  int result = Success;
  char errMsg[MAX_LEN_ERR_MSG + 1];
  char line[MAX_LEN_ONE_LINE_CONTENT + 1];
  BOOL readingAtom = FALSE;
  BOOL readingBond = FALSE;
  char resiName[MAX_LEN_RES_NAME + 1] = "";
  int bondNum = 0;
  int bonds[500][2];
  IntArray bondsFromToAndType;
  IntArrayCreate(&bondsFromToAndType, 0);

  AtomArray atoms;
  Atom atom;
  AtomArrayCreate(&atoms);
  AtomCreate(&atom);
  FILE* pMol2 = fopen(mol2file, "r");
  BOOL atomConflict = FALSE;
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pMol2) != NULL)
  {
    char keyword[MAX_LEN_ONE_LINE_CONTENT + 1] = "";
    sscanf(line, "%s", keyword);
    if (strcmp(keyword, "@<TRIPOS>ATOM") == 0)
    {
      readingAtom = TRUE;
      readingBond = FALSE;
      continue;
    }
    else if (strcmp(keyword, "@<TRIPOS>BOND") == 0)
    {
      readingAtom = FALSE;
      readingBond = TRUE;
      continue;
    }
    else if (keyword[0] == '@')
    {
      readingAtom = readingBond = FALSE;
      continue;
    }
    if (readingAtom == TRUE)
    {
      int atomId;
      char subst_name[MAX_LEN_RES_NAME + 1];
      // the subst_id in field 7 of TRIPOS mol2 is used as the posInChain of the ligand
      int count = sscanf(line, "%d %s %lf %lf %lf %s %d %s %lf", &atomId, atom.name, &atom.xyz.X, &atom.xyz.Y, &atom.xyz.Z, atom.type, &atom.posInChain, subst_name, &atom.charge);
      if (count < 9)
      {
        char errMsg[MAX_LEN_ERR_MSG + 1];
        sprintf(errMsg, "in file %s line %d, probable format error when reading '%s' in file %s", __FILE__, __LINE__, line, mol2file);
        TraceError(errMsg, FormatError);
        exit(FormatError);
      }
      if (strcmp(resiName, "") == 0)
      {
        strcpy(resiName, subst_name);
        if (strlen(resiName) > 3) resiName[3] = '\0';
        if (strcmp(resiName, "***") == 0) strcpy(resiName, "LIG");
        if (LigandResidueNameConflictWithAminoAcid(resiName)) strcpy(resiName, "LIG");
      }
      AtomArrayAppend(&atoms, &atom);
      // check if the newly added atom has an identical name with existing ones
      if (atomConflict == FALSE)
      {
        for (int ndx = 0; ndx < AtomArrayGetCount(&atoms) - 1; ndx++)
        {
          Atom* pAtom = AtomArrayGet(&atoms, ndx);
          if (strcmp(AtomGetName(pAtom), AtomGetName(&atom)) == 0)
          {
            atomConflict = TRUE;
            break;
          }
        }
      }
    }
    if (readingBond == TRUE)
    {
      int id;
      char type[5];
      sscanf(line, "%d %d %d %s", &id, &bonds[bondNum][0], &bonds[bondNum][1], type);
      IntArrayAppend(&bondsFromToAndType, bonds[bondNum][0] - 1);
      IntArrayAppend(&bondsFromToAndType, bonds[bondNum][1] - 1);
      bondNum++;
    }
  }
  AtomDestroy(&atom);

  // rename atom name if atomConflict == TRUE
  if (atomConflict == TRUE)
  {
    for (int ndx = 0; ndx < AtomArrayGetCount(&atoms); ndx++)
    {
      char newName[MAX_LEN_ATOM_NAME + 1];
      sprintf(newName, "%s%d", AtomGetName(AtomArrayGet(&atoms, ndx)), ndx + 1);
      AtomSetName(AtomArrayGet(&atoms, ndx), newName);
    }
  }

  // create atom parameter for ligand
  for (int i = 0; i < AtomArrayGetCount(&atoms); i++)
  {
    int attachAtomTypeNum[10];
    int attachCount = 0;
    int attachIndex[10];
    for (int j = 0;j < 10;j++) attachAtomTypeNum[j] = 0;
    for (int j = 0;j < 10;j++) attachIndex[j] = -1;
    Atom* pAtom = AtomArrayGet(&atoms, i);
    pAtom->isBBAtom = FALSE;
    strcpy(pAtom->hbB2, "-");
    strcpy(pAtom->hbDorB, "-");
    strcpy(pAtom->hbHorA, "-");
    pAtom->hybridType = Type_AtomHybridType_None;
    pAtom->polarity = Type_AtomPolarity_NPAliphatic;

    if (pAtom->name[0] == 'C' && strstr(pAtom->type, "C.") != NULL)
    {
      GetAttachedAtomTypeNum(&atoms, bonds, bondNum, i, attachAtomTypeNum, &attachCount, attachIndex);
      if (attachAtomTypeNum[1] >= 3)
      { // atom i has >= 3 hydrogen atoms attached
        pAtom->EEF1_atType = EEF1_AtomType_CH3E;
      }
      else if (attachAtomTypeNum[1] == 2)
      { // atom i has 2 hydrogen atoms attached
        pAtom->EEF1_atType = EEF1_AtomType_CH2E;
      }
      else if (attachAtomTypeNum[1] == 1)
      { // atom i has 1 hydrogen atoms attached
        pAtom->EEF1_atType = EEF1_AtomType_CH1E;
        if (strcmp(pAtom->type, "C.ar") == 0)
        {
          pAtom->EEF1_atType = EEF1_AtomType_CR1E;
        }
      }
      else
      {
        if (strcmp(pAtom->type, "C.ar") == 0)
        {
          pAtom->EEF1_atType = EEF1_AtomType_CR;
        }
        else if (strcmp(pAtom->type, "C.cat") == 0)
        {
          pAtom->EEF1_atType = EEF1_AtomType_Carg;
        }
        else if (strcmp(pAtom->type, "C.2") == 0)
        {
          pAtom->EEF1_atType = EEF1_AtomType_CR;
          if (attachAtomTypeNum[2] > 0)
          {
            pAtom->EEF1_atType = EEF1_AtomType_CO;
          }
          for (int k = 0;k < attachCount;k++)
          {
            Atom* pAttachAtom = AtomArrayGet(&atoms, attachIndex[k]);
            if (strcmp(pAttachAtom->type, "O.co2") == 0 || pAttachAtom->EEF1_atType == EEF1_AtomType_OOC)
            {
              pAtom->EEF1_atType = EEF1_AtomType_COO;
              break;
            }
          }
        }
        else if (strcmp(pAtom->type, "C.3") == 0)
        {
          if (attachCount == 4) pAtom->EEF1_atType = EEF1_AtomType_CH1E;
          else if (attachCount == 3) pAtom->EEF1_atType = EEF1_AtomType_CH1E;
          else if (attachCount == 2) pAtom->EEF1_atType = EEF1_AtomType_CH2E;
          else if (attachCount == 1) pAtom->EEF1_atType = EEF1_AtomType_CH3E;
        }
        else pAtom->EEF1_atType = EEF1_AtomType_CR;
      }
    }
    else if (pAtom->name[0] == 'H')
    {
      pAtom->EEF1_atType = EEF1_AtomType_Hnpl;
      GetAttachedAtomTypeNum(&atoms, bonds, bondNum, i, attachAtomTypeNum, &attachCount, attachIndex);
      if (attachAtomTypeNum[2] > 0 || attachAtomTypeNum[3] > 0)
      {
        pAtom->EEF1_atType = EEF1_AtomType_Hpol;
        pAtom->polarity = Type_AtomPolarity_P;
        strcpy(pAtom->hbHorA, "H");
        strcpy(pAtom->hbDorB, AtomGetName(AtomArrayGet(&atoms, attachIndex[0])));
      }
    }
    else if (pAtom->name[0] == 'O')
    {
      GetAttachedAtomTypeNum(&atoms, bonds, bondNum, i, attachAtomTypeNum, &attachCount, attachIndex);
      pAtom->polarity = Type_AtomPolarity_P;
      if (strcmp(pAtom->type, "O.co2") == 0)
      {
        pAtom->EEF1_atType = EEF1_AtomType_OOC;
        pAtom->polarity = Type_AtomPolarity_C;
        strcpy(pAtom->hbHorA, "A");
        pAtom->hybridType = Type_AtomHybridType_SP2;
        for (int k = 0;k < attachCount;k++)
        {
          Atom* pAttachAtom = AtomArrayGet(&atoms, attachIndex[k]);
          if (pAttachAtom->name[0] != 'H')
          {
            strcpy(pAtom->hbDorB, AtomGetName(AtomArrayGet(&atoms, attachIndex[0])));
            break;
          }
        }
      }//COO-
      else if (strcmp(pAtom->type, "O.2") == 0)
      {
        pAtom->EEF1_atType = EEF1_AtomType_OC;
        strcpy(pAtom->hbHorA, "A");
        pAtom->hybridType = Type_AtomHybridType_SP2;
        if (attachAtomTypeNum[1] == 1)
        {
          pAtom->EEF1_atType = EEF1_AtomType_OH1;
        }
        for (int k = 0;k < attachCount;k++)
        {
          Atom* pAttachAtom = AtomArrayGet(&atoms, attachIndex[k]);
          if (pAttachAtom->name[0] != 'H')
          {
            strcpy(pAtom->hbDorB, AtomGetName(AtomArrayGet(&atoms, attachIndex[0])));
            break;
          }
        }
      }
      else if (strcmp(pAtom->type, "O.3") == 0)
      {
        pAtom->EEF1_atType = EEF1_AtomType_Oest;
        if (attachAtomTypeNum[1] == 1)
        {
          pAtom->EEF1_atType = EEF1_AtomType_OH1;
          strcpy(pAtom->hbHorA, "A");
          pAtom->hybridType = Type_AtomHybridType_SP3;
          for (int k = 0;k < attachCount;k++)
          {
            Atom* pAttachAtom = AtomArrayGet(&atoms, attachIndex[k]);
            if (pAttachAtom->name[0] != 'H')
            {
              strcpy(pAtom->hbDorB, AtomGetName(AtomArrayGet(&atoms, attachIndex[0])));
              break;
            }
          }
        }
        else if (attachAtomTypeNum[1] == 2)
        {
          pAtom->EEF1_atType = EEF1_AtomType_OH2;
          strcpy(pAtom->hbHorA, "A");
          pAtom->hybridType = Type_AtomHybridType_SP3;
          for (int k = 0;k < attachCount;k++)
          {
            Atom* pAttachAtom = AtomArrayGet(&atoms, attachIndex[k]);
            if (pAttachAtom->name[0] != 'H')
            {
              strcpy(pAtom->hbDorB, AtomGetName(AtomArrayGet(&atoms, attachIndex[0])));
              break;
            }
          }
        }
      }
      else
      {
        pAtom->EEF1_atType = EEF1_AtomType_Oest;
      }
    }
    else if (pAtom->name[0] == 'N' && strstr(pAtom->type, "N.") != NULL)
    {
      GetAttachedAtomTypeNum(&atoms, bonds, bondNum, i, attachAtomTypeNum, &attachCount, attachIndex);
      pAtom->polarity = Type_AtomPolarity_P;
      if (strcmp(pAtom->type, "N.4") == 0) pAtom->polarity = Type_AtomPolarity_C;
      if (attachAtomTypeNum[1] == 3)
      {
        pAtom->EEF1_atType = EEF1_AtomType_NH3;
        pAtom->hybridType = Type_AtomHybridType_SP3;
      }
      else if (attachAtomTypeNum[1] == 2)
      {
        pAtom->EEF1_atType = EEF1_AtomType_NH2;
        pAtom->hybridType = Type_AtomHybridType_SP3;
        if (strcmp(pAtom->type, "N.2") == 0) pAtom->hybridType = Type_AtomHybridType_SP2;
        for (int k = 0;k < attachCount;k++)
        {
          Atom* pAttachAtom = AtomArrayGet(&atoms, attachIndex[k]);
          if (strcmp(pAttachAtom->type, "C.cat") == 0 || pAttachAtom->EEF1_atType == EEF1_AtomType_Carg)
          {
            pAtom->EEF1_atType = EEF1_AtomType_Narg;
            break;
          }
        }
      }
      else if (attachAtomTypeNum[1] == 1)
      {
        pAtom->EEF1_atType = EEF1_AtomType_NH1;
        pAtom->hybridType = Type_AtomHybridType_SP3;
        if (strcmp(pAtom->type, "N.2") == 0) pAtom->hybridType = Type_AtomHybridType_SP2;
        for (int k = 0;k < attachCount;k++)
        {
          Atom* pAttachAtom = AtomArrayGet(&atoms, attachIndex[k]);
          if (strcmp(pAttachAtom->type, "C.cat") == 0 || pAttachAtom->EEF1_atType == EEF1_AtomType_Carg)
          {
            pAtom->EEF1_atType = EEF1_AtomType_Narg;
            break;
          }
        }
      }
      else
      {
        if (strcmp(pAtom->type, "N.ar") == 0 || strcmp(pAtom->type, "N.2") == 0)
        {
          pAtom->EEF1_atType = EEF1_AtomType_NR;
          if (attachCount < 3)
          {
            strcpy(pAtom->hbHorA, "A");
            pAtom->hybridType = Type_AtomHybridType_SP2;
            for (int k = 0;k < attachCount;k++)
            {
              Atom* pAttachAtom = AtomArrayGet(&atoms, attachIndex[k]);
              if (pAttachAtom->name[0] != 'H')
              {
                strcpy(pAtom->hbDorB, AtomGetName(AtomArrayGet(&atoms, attachIndex[0])));
                break;
              }
            }
          }
          else pAtom->EEF1_atType = EEF1_AtomType_Npro;
        }
        else if (strcmp(pAtom->type, "N.am") == 0) pAtom->EEF1_atType = EEF1_AtomType_Npro;
        else if (strcmp(pAtom->type, "N.3") == 0 || strcmp(pAtom->type, "N.pl3") == 0) pAtom->EEF1_atType = EEF1_AtomType_Npro;
        else if (attachCount == 3) pAtom->EEF1_atType = EEF1_AtomType_Npro;
        else if (strcmp(pAtom->type, "N.1") == 0 && attachCount == 1) pAtom->EEF1_atType = EEF1_AtomType_NR;
        else pAtom->EEF1_atType = EEF1_AtomType_NR;
      }
    }
    else if (pAtom->name[0] == 'P')
    {
      GetAttachedAtomTypeNum(&atoms, bonds, bondNum, i, attachAtomTypeNum, &attachCount, attachIndex);
      pAtom->EEF1_atType = EEF1_AtomType_P;
    }
    else if (pAtom->name[0] == 'S')
    {
      GetAttachedAtomTypeNum(&atoms, bonds, bondNum, i, attachAtomTypeNum, &attachCount, attachIndex);
      if (attachAtomTypeNum[1] == 1) pAtom->EEF1_atType = EEF1_AtomType_SH1E;
      else pAtom->EEF1_atType = EEF1_AtomType_S;
    }
    else if (pAtom->name[0] == 'B' && strcmp(pAtom->type, "B") == 0) pAtom->EEF1_atType = EEF1_AtomType_B;
    else if (pAtom->name[0] == 'F' && strcmp(pAtom->type, "F") == 0) pAtom->EEF1_atType = EEF1_AtomType_F;
    else if (pAtom->name[0] == 'C' && strcmp(pAtom->type, "Cl") == 0) pAtom->EEF1_atType = EEF1_AtomType_Cl;
    else if (pAtom->name[0] == 'B' && strcmp(pAtom->type, "Br") == 0) pAtom->EEF1_atType = EEF1_AtomType_Br;
    else if (pAtom->name[0] == 'I' && strcmp(pAtom->type, "I") == 0) pAtom->EEF1_atType = EEF1_AtomType_I;
    else if (pAtom->name[0] == 'F' && strcmp(pAtom->type, "Fe") == 0) pAtom->EEF1_atType = EEF1_AtomType_Fe;
    else if (pAtom->name[0] == 'Z' && strcmp(pAtom->type, "F") == 0) pAtom->EEF1_atType = EEF1_AtomType_Zn;
    else if (pAtom->name[0] == 'A' && strcmp(pAtom->type, "Al") == 0) pAtom->EEF1_atType = EEF1_AtomType_Al;
    else if (pAtom->name[0] == 'N' && strcmp(pAtom->type, "Na") == 0) pAtom->EEF1_atType = EEF1_AtomType_Na;
    else if (pAtom->name[0] == 'K' && strcmp(pAtom->type, "K") == 0) pAtom->EEF1_atType = EEF1_AtomType_K;
    else if (pAtom->name[0] == 'C' && strcmp(pAtom->type, "Ca") == 0) pAtom->EEF1_atType = EEF1_AtomType_Ca;
    else if (pAtom->name[0] == 'M' && strcmp(pAtom->type, "Mg") == 0) pAtom->EEF1_atType = EEF1_AtomType_Mg;
    else pAtom->EEF1_atType = EEF1_AtomType_Other;

    AssignAtomParameterByEEF1Type(pAtom, pAtom->EEF1_atType);
  }

  // output the parameter file
  FILE* pFilePara = fopen(parfile, "w");
  fprintf(pFilePara, "!this is an automatically generated atom parameter file for residue %s\n", resiName);
  fprintf(pFilePara, "!Resi   Name   Type   Backbone   Polar   epsilon   Rmin     charge  HbH/A  HbD/B  HbB2  Hybrid  DG_free  Volume  Lambda\n");
  for (int i = 0; i < AtomArrayGetCount(&atoms); i++)
  {
    Atom* pAtom = AtomArrayGet(&atoms, i);
    char polarity[10];
    char hybrid[5];
    if (pAtom->polarity == Type_AtomPolarity_P) strcpy(polarity, "P");
    else if (pAtom->polarity == Type_AtomPolarity_C) strcpy(polarity, "C");
    else if (pAtom->polarity == Type_AtomPolarity_NPAliphatic) strcpy(polarity, "NP1");
    else strcpy(polarity, "NP2");
    if (pAtom->hybridType == Type_AtomHybridType_SP) strcpy(hybrid, "SP1");
    else if (pAtom->hybridType == Type_AtomHybridType_SP2) strcpy(hybrid, "SP2");
    else if (pAtom->hybridType == Type_AtomHybridType_SP3) strcpy(hybrid, "SP3");
    else strcpy(hybrid, "-");
    fprintf(pFilePara, "%-6s  %-6s %-6s %-10s %-7s %-8.4f %7.4f  %7.4f  %-5s  %-5s  %-5s %-4s   %7.2f     %4.1f    %3.1f\n",
      resiName, pAtom->name, pAtom->type, "N", polarity,
      pAtom->vdw_epsilon, pAtom->vdw_radius, pAtom->charge,
      pAtom->hbHorA, pAtom->hbDorB, pAtom->hbB2, hybrid,
      pAtom->EEF1_freeDG, pAtom->EEF1_volume, pAtom->EEF1_lamda_);
  }
  fclose(pFilePara);

  //////////////////////////////////////////////////////////////////////////////////////
  // output the topology file
  //////////////////////////////////////////////////////////////////////////////////////
  // write atoms and bonds
  FILE* pFileTopo = fopen(topfile, "w");
  if (pFileTopo == NULL)
  {
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, topfile);
    result = IOError;
    TraceError(errMsg, result);
    return result;
  }
  fprintf(pFileTopo, "!this is an automatically generated topology file for residue %s\n", resiName);
  fprintf(pFileTopo, "RESI %-3.3s\n", resiName);
  for (int i = 0; i < AtomArrayGetCount(&atoms); i++)
  {
    fprintf(pFileTopo, "ATOM  %s\n", AtomGetName(AtomArrayGet(&atoms, i)));
  }
  for (int i = 0; i < IntArrayGetLength(&bondsFromToAndType); i += 2)
  {
    fprintf(pFileTopo, "BOND  %s %s\n", AtomGetName(AtomArrayGet(&atoms, IntArrayGet(&bondsFromToAndType, i))), AtomGetName(AtomArrayGet(&atoms, IntArrayGet(&bondsFromToAndType, i + 1))));
  }
  // write ICs
  fprintf(pFileTopo, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  fprintf(pFileTopo, "! Internal Coordinates (ICs) are associated with ligand placement for enzyme design\n");
  fprintf(pFileTopo, "! The IC orders depend on three initial atoms on the ligand. If the three atoms\n");
  fprintf(pFileTopo, "! were not provided by user, they will be chosen automatically based on bond connection\n");
  fprintf(pFileTopo, "! Anyway, it is very important to check the ICs manually!\n");
  fprintf(pFileTopo, "! The format of an IC entry is: \n");
  fprintf(pFileTopo, "! IC atomA atomB  atomC atomD distAB angleABC angleABCD angleBCD distCD    or:\n");
  fprintf(pFileTopo, "! IC atomA atomB *atomC atomD distAC angleACB angleABCD angleBCD distCD\n");
  fprintf(pFileTopo, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

  Atom* pIniAtomA = AtomArrayGetByName(&atoms, iniAtom1);
  Atom* pIniAtomB = AtomArrayGetByName(&atoms, iniAtom2);
  Atom* pIniAtomC = AtomArrayGetByName(&atoms, iniAtom3);
  if (pIniAtomA == NULL || pIniAtomB == NULL || pIniAtomC == NULL)
  {
    pIniAtomA = AtomArrayGet(&atoms, IntArrayGet(&bondsFromToAndType, 0));
    pIniAtomB = AtomArrayGet(&atoms, IntArrayGet(&bondsFromToAndType, 1));
    //the third atom should form a bond with either A or B
    for (int i = 2; i < IntArrayGetLength(&bondsFromToAndType);i += 2)
    {
      int from = IntArrayGet(&bondsFromToAndType, i);
      int to = IntArrayGet(&bondsFromToAndType, i + 1);
      if (from == IntArrayGet(&bondsFromToAndType, 0) && to != IntArrayGet(&bondsFromToAndType, 1))
      {
        pIniAtomC = AtomArrayGet(&atoms, to);
        break;
      }
      else if (from == IntArrayGet(&bondsFromToAndType, 1) && to != IntArrayGet(&bondsFromToAndType, 0))
      {
        pIniAtomC = AtomArrayGet(&atoms, to);
        break;
      }
      else if (from != IntArrayGet(&bondsFromToAndType, 0) && to == IntArrayGet(&bondsFromToAndType, 1))
      {
        pIniAtomC = AtomArrayGet(&atoms, from);
        break;
      }
      else if (from != IntArrayGet(&bondsFromToAndType, 1) && to == IntArrayGet(&bondsFromToAndType, 0))
      {
        pIniAtomC = AtomArrayGet(&atoms, from);
        break;
      }
    }
    printf("The three initial atoms used for generating ligand topology are: %s, %s, and %s\n", AtomGetName(pIniAtomA), AtomGetName(pIniAtomB), AtomGetName(pIniAtomC));
    fprintf(pFileTopo, "The three initial atoms used for generating ligand topology are: %s, %s, and %s\n", AtomGetName(pIniAtomA), AtomGetName(pIniAtomB), AtomGetName(pIniAtomC));
    fprintf(pFileTopo, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  }


  IntArray icFlags;
  IntArrayCreate(&icFlags, AtomArrayGetCount(&atoms));
  for (int i = 0; i < IntArrayGetLength(&icFlags); i++)
  {
    Atom* pAtom = AtomArrayGet(&atoms, i);
    if (pAtom == pIniAtomA || pAtom == pIniAtomB || pAtom == pIniAtomC) IntArraySet(&icFlags, i, 1);
    else IntArraySet(&icFlags, i, 0);
  }

  BOOL allICsCalculated = FALSE;
  int maxIter = 0;
  while (allICsCalculated == FALSE && maxIter < AtomArrayGetCount(&atoms) + 1)
  {
    allICsCalculated = TRUE;
    for (int i = 0; i < AtomArrayGetCount(&atoms); i++)
    {
      if (IntArrayGet(&icFlags, i) == 1)
      {
        continue;
      }
      BOOL icFoundFlag = FALSE;
      Atom* pAtomD = AtomArrayGet(&atoms, i);
      AtomArray atomCArray;
      AtomArrayCreate(&atomCArray);
      // find atomCs bonded with atomD, where atomC's IC should have been calculated
      FindAtomC(pAtomD, &atoms, &icFlags, &bondsFromToAndType, &atomCArray);
      for (int j = 0; j < AtomArrayGetCount(&atomCArray); j++)
      {
        Atom* pAtomC = AtomArrayGet(&atomCArray, j);
        AtomArray atomBArray;
        AtomArrayCreate(&atomBArray);
        FindAtomB(pAtomC, pAtomD, &atoms, &icFlags, &bondsFromToAndType, &atomBArray);
        if (AtomArrayGetCount(&atomBArray) >= 2)
        { // improper IC
          Atom* pAtomA = AtomArrayGet(&atomBArray, 0);
          Atom* pAtomB = AtomArrayGet(&atomBArray, 1);
          Atom* pAtomX = NULL;
          if (strcmp(AtomGetName(pAtomB), AtomGetName(pIniAtomA)) == 0 ||
              strcmp(AtomGetName(pAtomB), AtomGetName(pIniAtomB)) == 0 ||
              strcmp(AtomGetName(pAtomB), AtomGetName(pIniAtomC)) == 0)
          {
            pAtomX = pAtomA;
            pAtomA = pAtomB;
            pAtomB = pAtomX;
          }
          XYZ AB = XYZDifference(&pAtomA->xyz, &pAtomB->xyz);
          XYZ BC = XYZDifference(&pAtomB->xyz, &pAtomC->xyz);
          XYZ CD = XYZDifference(&pAtomC->xyz, &pAtomD->xyz);
          XYZ AC = XYZDifference(&pAtomA->xyz, &pAtomC->xyz);
          double dAC = XYZDistance(&pAtomA->xyz, &pAtomC->xyz);
          double aACB = XYZAngle(&AC, &BC);
          double xABCD = GetTorsionAngle(&pAtomA->xyz, &pAtomB->xyz, &pAtomC->xyz, &pAtomD->xyz);
          double aBCD = PI - XYZAngle(&BC, &CD);
          double dCD = XYZDistance(&pAtomC->xyz, &pAtomD->xyz);
          fprintf(pFileTopo, "IC   %-7.7s %-7.7s *%-7.7s %-7.7s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
            AtomGetName(pAtomA), AtomGetName(pAtomB), AtomGetName(pAtomC), AtomGetName(pAtomD), dAC, RadToDeg(aACB), RadToDeg(xABCD), RadToDeg(aBCD), dCD);
          IntArraySet(&icFlags, i, 1);
          icFoundFlag = TRUE;
        }
        else
        { // proper IC
          for (int k = 0; k < AtomArrayGetCount(&atomBArray); k++)
          {
            Atom* pAtomB = AtomArrayGet(&atomBArray, k);
            AtomArray atomAArray;
            AtomArrayCreate(&atomAArray);
            // the method for finding aomA is identical to that of finding atomB
            FindAtomB(pAtomB, pAtomC, &atoms, &icFlags, &bondsFromToAndType, &atomAArray);
            if (AtomArrayGetCount(&atomAArray) == 0)
            {
              fprintf(pFileTopo, "! Atom %s's IC cannot be calculated because atomA cannot be found; its IC is ignored\n", AtomGetName(pAtomD));
              IntArraySet(&icFlags, i, 1);
              icFoundFlag = TRUE;
            }
            else
            {
              Atom* pAtomA = AtomArrayGet(&atomAArray, 0);
              XYZ AB = XYZDifference(&pAtomA->xyz, &pAtomB->xyz);
              XYZ BC = XYZDifference(&pAtomB->xyz, &pAtomC->xyz);
              XYZ CD = XYZDifference(&pAtomC->xyz, &pAtomD->xyz);
              XYZ AC = XYZDifference(&pAtomA->xyz, &pAtomC->xyz);
              double dAB = XYZDistance(&pAtomA->xyz, &pAtomB->xyz);
              double aABC = PI - XYZAngle(&AB, &BC);
              double xABCD = GetTorsionAngle(&pAtomA->xyz, &pAtomB->xyz, &pAtomC->xyz, &pAtomD->xyz);
              double aBCD = PI - XYZAngle(&BC, &CD);
              double dCD = XYZDistance(&pAtomC->xyz, &pAtomD->xyz);
              fprintf(pFileTopo, "IC   %-7.7s %-7.7s  %-7.7s %-7.7s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
                AtomGetName(pAtomA), AtomGetName(pAtomB), AtomGetName(pAtomC), AtomGetName(pAtomD), dAB, RadToDeg(aABC), RadToDeg(xABCD), RadToDeg(aBCD), dCD);
              IntArraySet(&icFlags, i, 1);
              icFoundFlag = TRUE;
            }
            AtomArrayDestroy(&atomAArray);
            if (icFoundFlag == TRUE) break;
          } // for k
        }
        AtomArrayDestroy(&atomBArray);
        if (icFoundFlag == TRUE) break;
      } // for j
      AtomArrayDestroy(&atomCArray);
      if (icFoundFlag == TRUE) break;
    } // for i

    for (int i = 0; i < IntArrayGetLength(&icFlags); i++)
    {
      if (IntArrayGet(&icFlags, i) == 0)
      {
        allICsCalculated = FALSE;
        break;
      }
    }
    maxIter++;
  } // while

  for (int i = 0; i < AtomArrayGetCount(&atoms); i++)
  {
    Atom* pAtom = AtomArrayGet(&atoms, i);
    if (IntArrayGet(&icFlags, i) == 0)
    {
      fprintf(pFileTopo, "! Atom %s's IC cannot be calculated because it is not connected to other atoms; its IC is ignored\n", AtomGetName(pAtom));
    }
  }
  fclose(pFileTopo);

  IntArrayDestroy(&icFlags);
  IntArrayDestroy(&bondsFromToAndType);
  AtomArrayDestroy(&atoms);

  return Success;
}
