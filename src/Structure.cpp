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
#pragma warning(disable:6301)
#pragma warning(disable:28182)

#include "Structure.h"
#include <string.h>

//#define DEBUGGING_STRUCTURE


int StructureCreate(Structure* pThis)
{
  strcpy(pThis->name, "");
  pThis->chainNum = 0;
  pThis->chains = NULL;
  pThis->desSiteCount = 0;
  pThis->designSites = NULL;
  return Success;
}
int StructureDestroy(Structure* pThis)
{
  //first free the memory for design sites
  for (int i = 0; i < pThis->desSiteCount; i++)
  {
    DesignSiteDestroy(&pThis->designSites[i]);
  }
  pThis->designSites = NULL;
  pThis->desSiteCount = 0;
  for (int i = 0;i < pThis->chainNum;i++)
  {
    ChainDestroy(&pThis->chains[i]);
  }
  free(pThis->chains);
  strcpy(pThis->name, "");
  pThis->chainNum = 0;
  pThis->chains = NULL;
  return Success;
}

char* StructureGetName(Structure* pThis)
{
  return pThis->name;
}

int StructureSetName(Structure* pThis, char* newName)
{
  if (strlen(newName) > MAX_LEN_STRUCTURE_NAME)
  {
    char errMsg[MAX_LEN_ERR_MSG + 1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s line %d, name is too long", __FILE__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->name, newName);
  return Success;
}

int StructureGetChainCount(Structure* pThis)
{
  return pThis->chainNum;
}

Chain* StructureGetChain(Structure* pThis, int index)
{
  if (index < 0 || index >= StructureGetChainCount(pThis))
  {
    return NULL;
  }
  return &pThis->chains[index];
}

Chain* StructureFindChainByName(Structure* pThis, char* chainName)
{
  int index = -1;
  int result = StructureFindChainIndex(pThis, chainName, &index);
  if (FAILED(result))
  {
    return NULL;
  }
  else
  {
    return StructureGetChain(pThis, index);
  }
}

int StructureFindChainIndex(Structure* pThis, char* chainName, int* index)
{
  int i;
  for (i = 0;i < pThis->chainNum;i++)
  {
    if (strcmp(ChainGetName(&pThis->chains[i]), chainName) == 0)
    {
      *index = i;
      return Success;
    }
  }
  return DataNotExistError;
}


int StructureAddChain(Structure* pThis, Chain* newChain)
{
  int index = -1;
  int result = StructureFindChainIndex(pThis, ChainGetName(newChain), &index);
  if (FAILED(result))
  {
    (pThis->chainNum)++;
    pThis->chains = (Chain*)realloc(pThis->chains, sizeof(Chain) * pThis->chainNum);
    ChainCreate(&pThis->chains[pThis->chainNum - 1]);
    return ChainCopy(&pThis->chains[pThis->chainNum - 1], newChain);
  }
  else
  {
    return ChainCopy(&pThis->chains[index], newChain);
  }
}

int StructureDeleteChain(Structure* pThis, char* chainName)
{
  int index;
  int result = StructureFindChainIndex(pThis, chainName, &index);
  if (FAILED(result))
    return result;
  for (int i = index;i < pThis->chainNum - 1;i++)
  {
    ChainCopy(&pThis->chains[i], &pThis->chains[i + 1]);
  }
  ChainDestroy(&pThis->chains[pThis->chainNum - 1]);
  (pThis->chainNum)--;
  return Success;
}

int StructureShowInPDBFormat(Structure* pThis, FILE* pFile)
{
  int atomIndex = 1;
  for (int i = 0;i < StructureGetChainCount(pThis);i++)
  {
    Chain* pChain = StructureGetChain(pThis, i);
    for (int j = 0; j < ChainGetResidueCount(pChain); j++)
    {
      Residue* pResi = ChainGetResidue(pChain, j);
      ResidueShowInPDBFormat(pResi, "ATOM", ResidueGetChainName(pResi), atomIndex, ResidueGetPosInChain(pResi), pFile);
      atomIndex += ResidueGetAtomCount(pResi);
    }
  }
  return Success;
}

int StructureReadPDB(Structure* pStructure, char* pdbFile, AtomParamsSet* pAtomParams, ResiTopoSet* pTopos)
{
  char initChainID[MAX_LEN_CHAIN_NAME + 1];
  char initResPos[6];
  strcpy(initChainID, "UNK");
  initChainID[strlen(initChainID)] = '\0';
  strcpy(initResPos, "UNK");
  initResPos[strlen(initResPos)] = '\0';
  BOOL firstResidueInChain = TRUE;

  FileReader file;
  if (FAILED(FileReaderCreate(&file, pdbFile)))
  {
    char errMsg[MAX_LEN_ERR_MSG + 1];
    sprintf(errMsg, "in file %s line %d, cannot read file %s", __FILE__, __LINE__, pdbFile);
    TraceError(errMsg, IOError);
    return IOError;
  }
  char line[MAX_LEN_ONE_LINE_CONTENT + 1];
  while (!FAILED(FileReaderGetNextLine(&file, line)))
  {
    char keyword[6] = "";
    char strAtomName[MAX_LEN_ATOM_NAME + 1] = "";
    char strResName[MAX_LEN_RES_NAME + 1] = "";
    char strChainID[MAX_LEN_CHAIN_NAME + 1] = "";
    char strResPos[6] = "";

    ExtractTargetStringFromSourceString(keyword, line, 0, 4);
    if (strcmp(keyword, "ATOM") != 0 && strcmp(keyword, "HETA") != 0) continue;

    ExtractTargetStringFromSourceString(strAtomName, line, 12, 4);
    ExtractTargetStringFromSourceString(strResName, line, 17, 4);
    ExtractTargetStringFromSourceString(strChainID, line, 21, 1);
    ExtractTargetStringFromSourceString(strResPos, line, 22, 5);
    if (strcmp(strChainID, "") == 0 || strcmp(strChainID, " ") == 0)
    {
      strcpy(strChainID, "A");
    }

    if (strcmp(initChainID, strChainID))
    { // new chain
      Chain* pChain = StructureGetChain(pStructure, StructureGetChainCount(pStructure) - 1);
      if (pChain)
      {
        if (ChainGetType(pChain) == Type_Chain_Protein)
        {
          Residue* pFirsResidueInChain = ChainGetResidue(pChain, 0);
          Residue* pLastResidueInChain = ChainGetResidue(pChain, ChainGetResidueCount(pChain) - 1);
          if (ResidueGetAtomByName(pFirsResidueInChain, "HT1") || ResidueGetAtomByName(pFirsResidueInChain, "HN1"))
          {
            ResiduePatchCTER(pLastResidueInChain, "CTER", pAtomParams, pTopos);
            pLastResidueInChain->terminalType = Type_ResIsCter;
          }
        }
        for (int i = 0;i < ChainGetResidueCount(pChain);i++)
        {
          Residue* pResi = ChainGetResidue(pChain, i);
          if (pResi->isSCIntact)
          {
            ResidueCalcAllAtomXYZ(pResi, pTopos, ChainGetResidue(pChain, i - 1), ChainGetResidue(pChain, i + 1));
            ResidueCalcSidechainTorsion(pResi, pTopos);
          }
          else ResidueCalcAllBackboneAtomXYZ(pResi, pTopos, ChainGetResidue(pChain, i - 1), ChainGetResidue(pChain, i + 1));
        }
      }

      //deal with new chains
      Type_Chain chainType = ChainTypeIdentifiedFromResidueName(strResName);
      if (chainType == Type_Chain_Unknown) continue;
      strcpy(initChainID, strChainID);
      Chain newChain;
      ChainCreate(&newChain);
      ChainSetType(&newChain, chainType);
      ChainSetName(&newChain, strChainID);
      StructureAddChain(pStructure, &newChain);
      ChainDestroy(&newChain);
      FileReaderSetCurrentPos(&file, FileReaderGetCurrentPos(&file) - 1);
      firstResidueInChain = TRUE;
    }
    else if (strcmp(initResPos, strResPos) != 0)
    { // new residue
      Residue newResi;
      ResidueCreate(&newResi);
      if (strcmp(strResName, "HIS") == 0)
      {
        FileReaderSetCurrentPos(&file, FileReaderGetCurrentPos(&file) - 1);
        int hisState = 1;
        PDBReaderCheckHisState(&file, &hisState);
        if (hisState == 1) strcpy(strResName, "HSD");
        else if (hisState == 2) strcpy(strResName, "HSE");
        else strcpy(strResName, "HSP");
      }

      ResidueSetName(&newResi, strResName);
      ResidueSetChainName(&newResi, strChainID);
      ResidueSetPosInChain(&newResi, atoi(strResPos));
      ResidueAddAtomsFromAtomParams(&newResi, pAtomParams);
      ResidueAddBondsFromResiTopos(&newResi, pTopos);
      if (firstResidueInChain && StructureGetChain(pStructure, StructureGetChainCount(pStructure) - 1)->type == Type_Chain_Protein)
      {
        ResiduePatchNTERorCTER(&newResi, "NTER", pAtomParams, pTopos);
        newResi.terminalType = Type_ResIsNter;
        firstResidueInChain = FALSE;
      }
      FileReaderSetCurrentPos(&file, FileReaderGetCurrentPos(&file) - 1);
      ResidueReadXYZFromPDB(&newResi, &file);
      ResidueCheckAtomCoordinateValidity(&newResi);
      if (StructureGetChain(pStructure, StructureGetChainCount(pStructure) - 1) == NULL)
      {
        char errMsg[MAX_LEN_ERR_MSG + 1];
        sprintf(errMsg, "in file %s line %d, cannot find the target chain in structure", __FILE__, __LINE__);
        TraceError(errMsg, DataNotExistError);
        return DataNotExistError;
      }
      else ChainAppendResidue(StructureGetChain(pStructure, StructureGetChainCount(pStructure) - 1), &newResi);
      ResidueDestroy(&newResi);
      strcpy(initResPos, strResPos);
    }
  }

  // handle the last chain
  Chain* pChain = StructureGetChain(pStructure, StructureGetChainCount(pStructure) - 1);
  if (pChain != NULL)
  {
    if (ChainGetType(pChain) == Type_Chain_Protein)
    {
      Residue* pFirsResidueInChain = ChainGetResidue(pChain, 0);
      Residue* pLastResidueInChain = ChainGetResidue(pChain, ChainGetResidueCount(pChain) - 1);
      if (ResidueGetAtomByName(pFirsResidueInChain, "HT1") != NULL || ResidueGetAtomByName(pFirsResidueInChain, "HN1") != NULL)
      {
        ResiduePatchCTER(pLastResidueInChain, "CTER", pAtomParams, pTopos);
        pLastResidueInChain->terminalType = Type_ResIsCter;
      }
    }
    for (int i = 0;i < ChainGetResidueCount(pChain);i++)
    {
      Residue* pResi = ChainGetResidue(pChain, i);
      if (pResi->isSCIntact)
      {
        ResidueCalcAllAtomXYZ(pResi, pTopos, ChainGetResidue(pChain, i - 1), ChainGetResidue(pChain, i + 1));
        ResidueCalcSidechainTorsion(pResi, pTopos);
      }
      else ResidueCalcAllBackboneAtomXYZ(pResi, pTopos, ChainGetResidue(pChain, i - 1), ChainGetResidue(pChain, i + 1));
    }
  }

  FileReaderDestroy(&file);
  return Success;
}


int StructureReadMol2(Structure* pStructure, char* mol2file, AtomParamsSet* pAtomParams, ResiTopoSet* pTopos)
{
  FILE* pMol2 = fopen(mol2file, "r");

  char line[MAX_LEN_ONE_LINE_CONTENT + 1];
  char keyword[MAX_LEN_ONE_LINE_CONTENT + 1];
  BOOL readingAtom = FALSE;
  BOOL readingBond = FALSE;

  char resiName[MAX_LEN_RES_NAME + 1] = "#";

  //add a new chain
  Chain newChain;
  Type_Chain chainType;
  ChainCreate(&newChain);
  chainType = Type_Chain_SmallMol;
  ChainSetType(&newChain, chainType);

  Residue newResi;
  ResidueCreate(&newResi);
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pMol2) != NULL)
  {
    strcpy(keyword, "");
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

    if (readingAtom)
    {
      char strAtomName[MAX_LEN_ATOM_NAME + 1] = "";
      char strResName[MAX_LEN_RES_NAME + 1] = "";
      char strChainID[MAX_LEN_CHAIN_NAME + 1] = "";
      char strResPos[MAX_LEN_ONE_LINE_CONTENT + 1] = "";
      double charge = 0.0;
      double x = 0, y = 0, z = 0;
      char type[MAX_LEN_ATOM_TYPE + 1] = "";

      int atomNdx;
      sscanf(line, "%d %s %lf %lf %lf %s %s %s %lf", &atomNdx, strAtomName, &x, &y, &z, type, strChainID, strResName, &charge);
      int posInChain = atoi(strChainID);
      if (resiName[0] == '#')
      {
        strcpy(resiName, strResName);
        if (strlen(resiName) > 3) resiName[3] = '\0';
        if (strcmp(resiName, "***") == 0) strcpy(resiName, "LIG");
        if (LigandResidueNameConflictWithAminoAcid(resiName)) strcpy(resiName, "LIG");
      }
      //deploy the residue first if it's not ready
      if (ResidueGetAtomCount(&newResi) == 0)
      {
        ResidueSetName(&newResi, resiName);
        ResidueSetPosInChain(&newResi, posInChain);
        ResidueAddAtomsFromAtomParams(&newResi, pAtomParams);
        ResidueAddBondsFromResiTopos(&newResi, pTopos);
        //set chain name
        //ChainSetName(&newChain, strChainID);
        ChainSetName(&newChain, "X");
      }
      //Atom* pAtom=ResidueGetAtomByName(&newResi,strAtomName);
      Atom* pAtom = ResidueGetAtom(&newResi, atomNdx - 1);
      pAtom->xyz.X = x;
      pAtom->xyz.Y = y;
      pAtom->xyz.Z = z;
      pAtom->isXyzValid = TRUE;
    }
  }
  ChainAppendResidue(&newChain, &newResi);
  StructureAddChain(pStructure, &newChain);
  ChainDestroy(&newChain);
  ResidueDestroy(&newResi);
  fclose(pMol2);

  return Success;
}



int ChainComputeResiduePosition(Structure* pStructure, int chainIndex)
{
  Chain* pChainI = StructureGetChain(pStructure, chainIndex);
  for (int ir = 0; ir < ChainGetResidueCount(pChainI); ir++)
  {
    Residue* pResiIR = ChainGetResidue(pChainI, ir);
    pResiIR->nCbIn10A = 0;
  }
  for (int ir = 0; ir < ChainGetResidueCount(pChainI); ir++)
  {
    Residue* pResiIR = ChainGetResidue(pChainI, ir);
    Atom* pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CB");
    if (pAtomCAorCB1 == NULL) pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CA");
    for (int is = ir + 1; is < ChainGetResidueCount(pChainI); is++)
    {
      Residue* pResiIS = ChainGetResidue(pChainI, is);
      Atom* pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CB");
      if (pAtomCAorCB2 == NULL) pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CA");
      if (XYZDistance(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 10.0)
      {
        pResiIR->nCbIn10A++;
        pResiIS->nCbIn10A++;
      }
    }
    //printf("Residue: %s %d %s, NumCBwithin8AtoCurResi: %d\n",ResidueGetChainName(pResiIR), ResidueGetPosInChain(pResiIR), ResidueGetName(pResiIR), pResiIR->numCBwithin8AtoCurResidue);
  }

  return Success;
}

int StructureComputeResiduePosition(Structure* pStructure)
{
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    for (int ir = 0; ir < ChainGetResidueCount(pChainI); ir++)
    {
      Residue* pResiIR = ChainGetResidue(pChainI, ir);
      pResiIR->nCbIn10A = 0;
    }
  }
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    for (int ir = 0; ir < ChainGetResidueCount(pChainI); ir++)
    {
      Residue* pResiIR = ChainGetResidue(pChainI, ir);
      Atom* pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CB");
      if (pAtomCAorCB1 == NULL) pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CA");
      for (int is = ir + 1; is < ChainGetResidueCount(pChainI); is++)
      {
        Residue* pResiIS = ChainGetResidue(pChainI, is);
        Atom* pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CB");
        if (pAtomCAorCB2 == NULL) pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CA");
        if (XYZDistance(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 10.0)
        {
          pResiIR->nCbIn10A++;
          pResiIS->nCbIn10A++;
        }
      }
      for (int k = i + 1; k < StructureGetChainCount(pStructure); k++)
      {
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) != Type_Chain_Protein) continue;
        for (int ks = 0; ks < ChainGetResidueCount(pChainK); ks++)
        {
          Residue* pResiKS = ChainGetResidue(pChainK, ks);
          Atom* pAtomCAorCB2 = pAtomCAorCB2 = ResidueGetAtomByName(pResiKS, "CB");
          if (pAtomCAorCB2 == NULL) pAtomCAorCB2 = ResidueGetAtomByName(pResiKS, "CA");
          if (XYZDistance(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 10.0)
          {
            pResiIR->nCbIn10A++;
            pResiKS->nCbIn10A++;
          }
        }
      }
      //printf("Residue: %s %d %s, NumCBwithin8AtoCurResi: %d\n",ResidueGetChainName(pResiIR), ResidueGetPosInChain(pResiIR), ResidueGetName(pResiIR), pResiIR->numCBwithin8AtoCurResidue);
    }
  }

  return Success;
}

int StructureShowAtomParameter(Structure* pStructure)
{
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    ChainShowAtomParameter(StructureGetChain(pStructure, i));
  }
  return Success;
}

int StructureShowBondInformation(Structure* pStructure)
{
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    ChainShowBondInformation(StructureGetChain(pStructure, i));
  }
  return Success;
}


int StructureCheckIntraBondType(Structure* pStructure)
{
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChain = StructureGetChain(pStructure, i);
    for (int j = 0; j < ChainGetResidueCount(pChain); j++)
    {
      Residue* pResidue = ChainGetResidue(pChain, j);
      for (int k = 0; k < ResidueGetAtomCount(pResidue); k++)
      {
        Atom* pAtom = ResidueGetAtom(pResidue, k);
        for (int m = k + 1; m < ResidueGetAtomCount(pResidue); m++)
        {
          Atom* pAtom2 = ResidueGetAtom(pResidue, m);
          //int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom), AtomGetName(pAtom2), pRes);
          int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom), AtomGetName(pAtom2), ResidueGetBonds(pResidue));
          //printf("Residue: %s %d %s, Atom1: %s, Atom2: %s, BondType: %d\n",ResidueGetName(pRes), ResidueGetPosInChain(pRes), ResidueGetChainName(pRes),AtomGetName(pAtom),AtomGetName(pAtom2),bondType);
        }
      }
    }
  }
  return Success;
}

int StructureCheckNeighbouringBondType(Structure* pStructure)
{
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChain = StructureGetChain(pStructure, i);
    for (int j = 0; j < ChainGetResidueCount(pChain) - 1; j++)
    {
      Residue* pResi1 = ChainGetResidue(pChain, j);
      Residue* pResi2 = ChainGetResidue(pChain, j + 1);
      for (int p = 0; p < pResi1->atoms.atomNum; p++)
      {
        Atom* pAtom1 = ResidueGetAtom(pResi1, p);
        for (int q = 0; q < pResi2->atoms.atomNum; q++)
        {
          Atom* pAtom2 = ResidueGetAtom(pResi2, q);
          int bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(pAtom1->name, pAtom2->name, pResi2->name);
          printf("Residue: %s %d %s, Atom1: %s, Residue: %s %d %s, Atom2: %s, BondType: %d\n", ResidueGetName(pResi1), ResidueGetPosInChain(pResi1), ResidueGetChainName(pResi1), AtomGetName(pAtom1), ResidueGetName(pResi2), ResidueGetPosInChain(pResi2), ResidueGetChainName(pResi2), AtomGetName(pAtom2), bondType);
        }
      }
    }
  }

  return Success;
}

///////////////////////////////////////////////
//functions for dealing with design sites
//////////////////////////////////////////////
int StructureGetDesignSiteCount(Structure* pThis)
{
  return pThis->desSiteCount;
}

DesignSite* StructureGetDesignSite(Structure* pThis, int index)
{
  if (index < 0 || index >= pThis->desSiteCount) return NULL;
  else return &pThis->designSites[index];
}

DesignSite* StructureFindDesignSite(Structure* pThis, int chainIndex, int resiIndex)
{
  if (chainIndex < 0 || chainIndex >= StructureGetChainCount(pThis)) return NULL;
  Chain* pChain = StructureGetChain(pThis, chainIndex);
  if (resiIndex < 0 || resiIndex >= ChainGetResidueCount(pChain)) return NULL;
  for (int i = 0; i < StructureGetDesignSiteCount(pThis); i++)
  {
    DesignSite* pDesignSite = StructureGetDesignSite(pThis, i);
    if (pDesignSite->chnNdx == chainIndex && pDesignSite->resNdx == resiIndex)
    {
      return pDesignSite;
    }
  }
  return NULL;
}


DesignSite* StructureFindDesignSiteByChainName(Structure* pStructure, char* chainName, int posInChain)
{
  for (int i = 0; i < StructureGetDesignSiteCount(pStructure); i++)
  {
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure, i);
    if (strcmp(chainName, DesignSiteGetChainName(pDesignSite)) == 0 &&
      DesignSiteGetPosInChain(pDesignSite) == posInChain)
    {
      return pDesignSite;
    }
  }
  return NULL;
}


int StructureFindDesignSiteIndexByChainNameAndPosInChain(Structure* pStructure, char* chainName, int posInChain)
{
  for (int i = 0; i < StructureGetDesignSiteCount(pStructure); i++)
  {
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure, i);
    if (strcmp(DesignSiteGetChainName(pDesignSite), chainName) == 0 &&
      DesignSiteGetPosInChain(pDesignSite) == posInChain)
    {
      return i;
    }
  }
  return -1;
}


int StructureGetDesignSiteIndex(Structure* pStructure, int chainIndex, int resiIndex)
{
  for (int i = 0; i < StructureGetDesignSiteCount(pStructure); i++)
  {
    DesignSite* pDesignSite = StructureGetDesignSite(pStructure, i);
    if (pDesignSite->chnNdx == chainIndex && pDesignSite->resNdx == resiIndex)
    {
      return i;
    }
  }
  return -1;
}


int StructureShowDesignSites(Structure* pThis)
{
  double confSpace = 0.0;
  int totRotCount = 0;
  for (int i = 0; i < pThis->desSiteCount; i++)
  {
    DesignSite* pSite = StructureGetDesignSite(pThis, i);
    int rotCount = RotamerSetGetCount(DesignSiteGetRotamers(pSite));
    printf("site %3d : %3s %s %4d, %6d rotamers, ",
      i, ResidueGetName(pSite->pRes), ResidueGetChainName(pSite->pRes), ResidueGetPosInChain(pSite->pRes), rotCount);
    totRotCount += rotCount;
    if (rotCount > 0)
    {
      confSpace += log((double)rotCount) / log(10.0);
    }
    switch (pSite->pRes->desType)
    {
    case Type_DesType_Catalytic:
      printf("catalytic\n"); break;
    case Type_DesType_Fixed:
      printf("fixed\n"); break;
    case Type_DesType_Mutable:
      printf("mutated\n"); break;
    case Type_DesType_SmallMol:
      printf("smallmol\n"); break;
    case Type_DesType_Repackable:
      printf("rotameric\n"); break;
    case Type_DesType_NatRot:
      printf("natrot\n"); break;
    default:
      break;
    }
  }
  printf("total rotamer count: %d, conformation space: 1e^%.1f\n", totRotCount, confSpace);
  return Success;
}

int ProteinSiteAddDesignSite(Structure* pThis, int chainIndex, int resiIndex)
{
  //this function add a new empty design site
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if (pCurrentDesignSite != NULL)
  {
    return Success;
  }
  (pThis->desSiteCount)++;
  pThis->designSites = (DesignSite*)realloc(pThis->designSites, sizeof(DesignSite) * pThis->desSiteCount);
  DesignSiteCreate(&pThis->designSites[pThis->desSiteCount - 1]);
  pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->desSiteCount - 1);
  pCurrentDesignSite->pRes = ChainGetResidue(StructureGetChain(pThis, chainIndex), resiIndex);
  pCurrentDesignSite->chnNdx = chainIndex;
  pCurrentDesignSite->resNdx = resiIndex;
  return result;
}

int StructureRemoveAllDesignSites(Structure* pThis)
{
  for (int i = 0; i < pThis->desSiteCount; i++)
  {
    DesignSite* pDesignSite = &pThis->designSites[i];
    DesignSiteDestroy(pDesignSite);
    pDesignSite = NULL;
  }
  free(pThis->designSites);
  pThis->designSites = NULL;
  pThis->desSiteCount = 0;
  return Success;
}

int ProteinSiteRemoveDesignSite(Structure* pThis, int chainIndex, int resiIndex)
{
  DesignSite* pDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if (pDesignSite != NULL)
  {
    DesignSiteDestroy(pDesignSite);
    pDesignSite = NULL;
  }
  (pThis->desSiteCount)--;
  return Success;
}



int StructureFindSmallMol(Structure* pThis, Residue** ppSmallMol)
{
  int result = DataNotExistError;
  for (int i = 0;i < pThis->chainNum;i++)
  {
    if (pThis->chains[i].type == Type_Chain_SmallMol)
    {
      *ppSmallMol = ChainGetResidue(&pThis->chains[i], 0);
      if (*ppSmallMol != NULL)
      {
        result = Success;
        break;
      }
    }
  }
  return result;
}


int StructureCalcProteinResidueSidechainTorsion(Structure* pThis, ResiTopoSet* pTopos)
{
  for (int i = 0;i < StructureGetChainCount(pThis);i++)
  {
    Chain* pChain = StructureGetChain(pThis, i);
    if (ChainGetType(pChain) == Type_Chain_Protein)
    {
      for (int j = 0;j < ChainGetResidueCount(pChain);j++)
      {
        Residue* pResidue = ChainGetResidue(pChain, j);
        if (strcmp(ResidueGetName(pResidue), "ALA") == 0 || strcmp(ResidueGetName(pResidue), "GLY") == 0) continue;
        ResidueCalcSidechainTorsion(pResidue, pTopos);
      }
    }
  }
  return Success;
}


int StructureCalcPhiPsi(Structure* pStructure)
{
  for (int i = 0;i < StructureGetChainCount(pStructure);i++)
  {
    Chain* pChain = StructureGetChain(pStructure, i);
    if (ChainGetType(pChain) != Type_Chain_Protein) continue;
    for (int j = 0;j < ChainGetResidueCount(pChain);j++)
    {
      Residue* pResi0 = ChainGetResidue(pChain, j - 1);
      Residue* pResi1 = ChainGetResidue(pChain, j);
      Residue* pResi2 = ChainGetResidue(pChain, j + 1);

      if (pResi0 != NULL && AA3GetIndex(ResidueGetName(pResi0)) < 0)  continue;
      if (pResi1 != NULL && AA3GetIndex(ResidueGetName(pResi1)) < 0)  continue;
      if (pResi2 != NULL && AA3GetIndex(ResidueGetName(pResi2)) < 0)  continue;


      Atom* pC0 = NULL;
      Atom* pN1 = ResidueGetAtomByName(pResi1, "N");
      Atom* pCA1 = ResidueGetAtomByName(pResi1, "CA");
      Atom* pC1 = ResidueGetAtomByName(pResi1, "C");
      Atom* pN2 = NULL;

      double phi = 0.0;
      double psi = 0.0;
      if (pResi0 == NULL)
      {
        phi = -60.0;
        pC0 = NULL;
      }
      else
      {
        pC0 = ResidueGetAtomByName(pResi0, "C");
      }

      if (pResi2 == NULL)
      {
        psi = 60.0;
        pN2 = NULL;
      }
      else
      {
        pN2 = ResidueGetAtomByName(pResi2, "N");
      }

      if (pC0 != NULL)
      {
        double dist_C0N1 = XYZDistance(&pC0->xyz, &pN1->xyz);
        if (dist_C0N1 < 1.45 && dist_C0N1>1.25)
        {// indicative of an amide bond
          phi = RadToDeg(GetTorsionAngle(&pC0->xyz, &pN1->xyz, &pCA1->xyz, &pC1->xyz));
        }
        else
        {
          phi = -60.0;
        }
      }

      if (pN2 != NULL)
      {
        double dist_C1N2 = XYZDistance(&pC1->xyz, &pN2->xyz);
        if (dist_C1N2 < 1.45 && dist_C1N2>1.25)
        {// indicative of an amide bond
          psi = RadToDeg(GetTorsionAngle(&pN1->xyz, &pCA1->xyz, &pC1->xyz, &pN2->xyz));
        }
        else
        {
          psi = 60.0;
        }
      }

      //assign phi and psi values to residue1
      pResi1->phipsi[0] = phi;
      pResi1->phipsi[1] = psi;
    }
  }

  return Success;
}


int StructureCopy(Structure* pThis, Structure* pOther)
{
  for (int i = 0;i < StructureGetChainCount(pOther);i++)
  {
    //Chain tempChain;
    //ChainCreate(&tempChain);
    //ChainCopy(&tempChain,&pOther->chains[i]);
    StructureAddChain(pThis, &pOther->chains[i]);
    //ChainDestroy(&tempChain);
  }
  pThis->chainNum = pOther->chainNum;
  pThis->desSiteCount = pOther->desSiteCount;
  strcpy(pThis->name, pOther->name);
  for (int i = 0;i < pOther->desSiteCount;i++)
  {
    DesignSiteCopy(&pThis->designSites[i], &pOther->designSites[i]);
  }
  return Success;
}

//deal with nucleic acid
Residue* StructureFindNucleotidePair(Structure* pStructure, Residue* pQuery)
{
  Chain* pChainI = StructureFindChainByName(pStructure, ResidueGetChainName(pQuery));
  if (ChainGetType(pChainI) == Type_Chain_DNA)
  {
    if (!strcmp(ResidueGetName(pQuery), "DA"))
    {
      //paired residue is DT or U;
      for (int k = 0;k < StructureGetChainCount(pStructure);k++)
      {
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            //paired residue must be DT;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pQuery) - ResidueGetPosInChain(pResiS)) <= 1) continue;
              if (!strcmp(ResidueGetName(pResiS), "DT"))
              {
                //check the distance
                //atom N1 and N6 from nucleotide DA
                //atom N3 and O4 from nucleotide DT
                double distN1N3 = XYZDistance(&ResidueGetAtomByName(pQuery, "N1")->xyz, &ResidueGetAtomByName(pResiS, "N3")->xyz);
                double distN6O4 = XYZDistance(&ResidueGetAtomByName(pQuery, "N6")->xyz, &ResidueGetAtomByName(pResiS, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
          else
          {
            //paired residue must be DT or U;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResiS), "DT") == 0 || strcmp(ResidueGetName(pResiS), "U") == 0)
              {
                //check the distance
                //atom N1 and N6 from nucleotide DA
                //atom N3 and O4 from nucleotide DT or U
                double distN1N3 = XYZDistance(&ResidueGetAtomByName(pQuery, "N1")->xyz, &ResidueGetAtomByName(pResiS, "N3")->xyz);
                double distN6O4 = XYZDistance(&ResidueGetAtomByName(pQuery, "N6")->xyz, &ResidueGetAtomByName(pResiS, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
        }
      }
    }
    else if (!strcmp(ResidueGetName(pQuery), "DC"))
    {
      //paired residue is DG or G;
      for (int k = 0;k < StructureGetChainCount(pStructure);k++)
      {
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            //paired residue must be DG;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pQuery) - ResidueGetPosInChain(pResiS)) <= 1) continue;
              if (!strcmp(ResidueGetName(pResiS), "DG"))
              {
                //check the distance
                //atom O2, N3 and N4 from nucleotide DC
                //atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDistance(&ResidueGetAtomByName(pQuery, "O2")->xyz, &ResidueGetAtomByName(pResiS, "N2")->xyz);
                double distN3N1 = XYZDistance(&ResidueGetAtomByName(pQuery, "N3")->xyz, &ResidueGetAtomByName(pResiS, "N1")->xyz);
                double distN4O6 = XYZDistance(&ResidueGetAtomByName(pQuery, "N4")->xyz, &ResidueGetAtomByName(pResiS, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
          else
          {
            //paired residue must be DG or G;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResiS), "DG") == 0 || strcmp(ResidueGetName(pResiS), "G") == 0)
              {
                //check the distance
                //atom O2, N3 and N4 from nucleotide DC
                //atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDistance(&ResidueGetAtomByName(pQuery, "O2")->xyz, &ResidueGetAtomByName(pResiS, "N2")->xyz);
                double distN3N1 = XYZDistance(&ResidueGetAtomByName(pQuery, "N3")->xyz, &ResidueGetAtomByName(pResiS, "N1")->xyz);
                double distN4O6 = XYZDistance(&ResidueGetAtomByName(pQuery, "N4")->xyz, &ResidueGetAtomByName(pResiS, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
        }
      }
    }
    else if (!strcmp(ResidueGetName(pQuery), "DG"))
    {
      //paired residue is DC or C;
      for (int k = 0;k < StructureGetChainCount(pStructure);k++)
      {
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            //paired residue must be DC;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pQuery) - ResidueGetPosInChain(pResiS)) <= 1) continue;
              if (!strcmp(ResidueGetName(pResiS), "DC"))
              {
                //check the distance
                //atom O2, N3 and N4 from nucleotide DC
                //atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDistance(&ResidueGetAtomByName(pResiS, "O2")->xyz, &ResidueGetAtomByName(pQuery, "N2")->xyz);
                double distN3N1 = XYZDistance(&ResidueGetAtomByName(pResiS, "N3")->xyz, &ResidueGetAtomByName(pQuery, "N1")->xyz);
                double distN4O6 = XYZDistance(&ResidueGetAtomByName(pResiS, "N4")->xyz, &ResidueGetAtomByName(pQuery, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
          else
          {
            //paired residue must be DC or C;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResiS), "DC") == 0 || strcmp(ResidueGetName(pResiS), "C") == 0)
              {
                //check the distance
                //atom O2, N3 and N4 from nucleotide DC
                //atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDistance(&ResidueGetAtomByName(pResiS, "O2")->xyz, &ResidueGetAtomByName(pQuery, "N2")->xyz);
                double distN3N1 = XYZDistance(&ResidueGetAtomByName(pResiS, "N3")->xyz, &ResidueGetAtomByName(pQuery, "N1")->xyz);
                double distN4O6 = XYZDistance(&ResidueGetAtomByName(pResiS, "N4")->xyz, &ResidueGetAtomByName(pQuery, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
        }
      }
    }
    else
    {//!strcmp(ResidueGetName(pQuery),"DT")
  //paired residue is DA or A;
      for (int k = 0;k < StructureGetChainCount(pStructure);k++)
      {
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            //paired residue must be DA;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pQuery) - ResidueGetPosInChain(pResiS)) <= 1) continue;
              if (!strcmp(ResidueGetName(pResiS), "DA"))
              {
                //check the distance
                //atom N1 and N6 from nucleotide DA or A
                //atom N3 and O4 from nucleotide DT
                double distN1N3 = XYZDistance(&ResidueGetAtomByName(pResiS, "N1")->xyz, &ResidueGetAtomByName(pQuery, "N3")->xyz);
                double distN6O4 = XYZDistance(&ResidueGetAtomByName(pResiS, "N6")->xyz, &ResidueGetAtomByName(pQuery, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
          else
          {
            //paired residue must be DA or A;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResiS), "DA") == 0 || strcmp(ResidueGetName(pResiS), "A") == 0)
              {
                //check the distance
                //atom N1 and N6 from nucleotide DA or A
                //atom N3 and O4 from nucleotide DT
                double distN1N3 = XYZDistance(&ResidueGetAtomByName(pResiS, "N1")->xyz, &ResidueGetAtomByName(pQuery, "N3")->xyz);
                double distN6O4 = XYZDistance(&ResidueGetAtomByName(pResiS, "N6")->xyz, &ResidueGetAtomByName(pQuery, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
        }
      }
    }
  }
  else if (ChainGetType(pChainI) == Type_Chain_RNA)
  {
    if (!strcmp(ResidueGetName(pQuery), "A"))
    {
      //paired residue is DT or U;
      for (int k = 0;k < StructureGetChainCount(pStructure);k++)
      {
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            //paired residue must be U;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pQuery) - ResidueGetPosInChain(pResiS)) <= 1) continue;
              if (!strcmp(ResidueGetName(pResiS), "U"))
              {
                //check the distance
                //atom N1 and N6 from nucleotide A
                //atom N3 and O4 from nucleotide U
                double distN1N3 = XYZDistance(&ResidueGetAtomByName(pQuery, "N1")->xyz, &ResidueGetAtomByName(pResiS, "N3")->xyz);
                double distN6O4 = XYZDistance(&ResidueGetAtomByName(pQuery, "N6")->xyz, &ResidueGetAtomByName(pResiS, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
          else
          {
            //paired residue must be DT or U;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResiS), "DT") == 0 || strcmp(ResidueGetName(pResiS), "U") == 0)
              {
                //check the distance
                //atom N1 and N6 from nucleotide A
                //atom N3 and O4 from nucleotide DT or U
                double distN1N3 = XYZDistance(&ResidueGetAtomByName(pQuery, "N1")->xyz, &ResidueGetAtomByName(pResiS, "N3")->xyz);
                double distN6O4 = XYZDistance(&ResidueGetAtomByName(pQuery, "N6")->xyz, &ResidueGetAtomByName(pResiS, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
        }
      }
    }
    else if (!strcmp(ResidueGetName(pQuery), "C"))
    {
      //paired residue is DG or G;
      for (int k = 0;k < StructureGetChainCount(pStructure);k++)
      {
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            //paired residue must be G;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pQuery) - ResidueGetPosInChain(pResiS)) <= 1) continue;
              if (!strcmp(ResidueGetName(pResiS), "G"))
              {
                //check the distance
                //atom O2, N3 and N4 from nucleotide C
                //atom N2, N1 and O6 from nucleotide G
                double distO2N2 = XYZDistance(&ResidueGetAtomByName(pQuery, "O2")->xyz, &ResidueGetAtomByName(pResiS, "N2")->xyz);
                double distN3N1 = XYZDistance(&ResidueGetAtomByName(pQuery, "N3")->xyz, &ResidueGetAtomByName(pResiS, "N1")->xyz);
                double distN4O6 = XYZDistance(&ResidueGetAtomByName(pQuery, "N4")->xyz, &ResidueGetAtomByName(pResiS, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
          else
          {
            //paired residue must be DG or G;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResiS), "DG") == 0 || strcmp(ResidueGetName(pResiS), "G") == 0)
              {
                //check the distance
                //atom O2, N3 and N4 from nucleotide DC
                //atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDistance(&ResidueGetAtomByName(pQuery, "O2")->xyz, &ResidueGetAtomByName(pResiS, "N2")->xyz);
                double distN3N1 = XYZDistance(&ResidueGetAtomByName(pQuery, "N3")->xyz, &ResidueGetAtomByName(pResiS, "N1")->xyz);
                double distN4O6 = XYZDistance(&ResidueGetAtomByName(pQuery, "N4")->xyz, &ResidueGetAtomByName(pResiS, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
        }
      }
    }
    else if (!strcmp(ResidueGetName(pQuery), "G"))
    {
      //paired residue is DC or C;
      for (int k = 0;k < StructureGetChainCount(pStructure);k++)
      {
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            //paired residue must be C;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pQuery) - ResidueGetPosInChain(pResiS)) <= 1) continue;
              if (!strcmp(ResidueGetName(pResiS), "C"))
              {
                //check the distance
                //atom O2, N3 and N4 from nucleotide C
                //atom N2, N1 and O6 from nucleotide G
                double distO2N2 = XYZDistance(&ResidueGetAtomByName(pResiS, "O2")->xyz, &ResidueGetAtomByName(pQuery, "N2")->xyz);
                double distN3N1 = XYZDistance(&ResidueGetAtomByName(pResiS, "N3")->xyz, &ResidueGetAtomByName(pQuery, "N1")->xyz);
                double distN4O6 = XYZDistance(&ResidueGetAtomByName(pResiS, "N4")->xyz, &ResidueGetAtomByName(pQuery, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
          else
          {
            //paired residue must be DC or C;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResiS), "DC") == 0 || strcmp(ResidueGetName(pResiS), "C") == 0)
              {
                //check the distance
                //atom O2, N3 and N4 from nucleotide DC
                //atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDistance(&ResidueGetAtomByName(pResiS, "O2")->xyz, &ResidueGetAtomByName(pQuery, "N2")->xyz);
                double distN3N1 = XYZDistance(&ResidueGetAtomByName(pResiS, "N3")->xyz, &ResidueGetAtomByName(pQuery, "N1")->xyz);
                double distN4O6 = XYZDistance(&ResidueGetAtomByName(pResiS, "N4")->xyz, &ResidueGetAtomByName(pQuery, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
        }
      }
    }
    else
    {//!strcmp(ResidueGetName(pQuery),"U")
  //paired residue is DA or A;
      for (int k = 0;k < StructureGetChainCount(pStructure);k++)
      {
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            //paired residue must be A;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pQuery) - ResidueGetPosInChain(pResiS)) <= 1) continue;
              if (!strcmp(ResidueGetName(pResiS), "A"))
              {
                //check the distance
                //atom N1 and N6 from nucleotide A
                //atom N3 and O4 from nucleotide U
                double distN1N3 = XYZDistance(&ResidueGetAtomByName(pResiS, "N1")->xyz, &ResidueGetAtomByName(pQuery, "N3")->xyz);
                double distN6O4 = XYZDistance(&ResidueGetAtomByName(pResiS, "N6")->xyz, &ResidueGetAtomByName(pQuery, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
          else
          {
            //paired residue must be DA or A;
            for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
            {
              Residue* pResiS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResiS), "DA") == 0 || strcmp(ResidueGetName(pResiS), "A") == 0)
              {
                //check the distance
                //atom N1 and N6 from nucleotide DA or A
                //atom N3 and O4 from nucleotide U
                double distN1N3 = XYZDistance(&ResidueGetAtomByName(pResiS, "N1")->xyz, &ResidueGetAtomByName(pQuery, "N3")->xyz);
                double distN6O4 = XYZDistance(&ResidueGetAtomByName(pResiS, "N6")->xyz, &ResidueGetAtomByName(pQuery, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResiS;
                }
              }
            }
          }
        }
      }
    }
  }

  return NULL;
}

int StructureMutateNucleotide(Residue* pQuery, char oneLetterName, Structure* pStructure, AtomParamsSet* pAtomParaSet, ResiTopoSet* pResiTopos)
{
  char newCleotideName[MAX_LEN_RES_NAME + 1];
  Chain* pChainI = StructureFindChainByName(pStructure, ResidueGetChainName(pQuery));
  if (ChainGetType(pChainI) == Type_Chain_DNA)
  {
    if (oneLetterName == 'a') { strcpy(newCleotideName, "DA"); }
    else if (oneLetterName == 'c') { strcpy(newCleotideName, "DC"); }
    else if (oneLetterName == 'g') { strcpy(newCleotideName, "DG"); }
    else if (oneLetterName == 't') { strcpy(newCleotideName, "DT"); }
    else if (oneLetterName == 'n') { strcpy(newCleotideName, "DN"); }
    else
    {
      char errMsg[MAX_LEN_ERR_MSG + 1];
      sprintf(errMsg, "in file %s line %d, DNA does not have nucleotide type '%c'", __FILE__, __LINE__, oneLetterName);
      TraceError(errMsg, ValueError);
      return ValueError;
    }
  }
  else if (ChainGetType(pChainI) == Type_Chain_RNA)
  {
    if (oneLetterName == 'a') { strcpy(newCleotideName, "A"); }
    else if (oneLetterName == 'c') { strcpy(newCleotideName, "C"); }
    else if (oneLetterName == 'g') { strcpy(newCleotideName, "G"); }
    else if (oneLetterName == 'u') { strcpy(newCleotideName, "U"); }
    else if (oneLetterName == 'n') { strcpy(newCleotideName, "N"); }
    else
    {
      char errMsg[MAX_LEN_ERR_MSG + 1];
      sprintf(errMsg, "in file %s line %d, RNA does not have nucleotide type '%c'", __FILE__, __LINE__, oneLetterName);
      TraceError(errMsg, ValueError);
      return ValueError;
    }
  }

  Residue newResi;
  ResidueCreate(&newResi);
  ResidueCopy(&newResi, pQuery);
  AtomArrayDestroy(&newResi.atoms);
  BondSetDestroy(&newResi.bonds);
  //modify the new residue
  ResidueSetName(&newResi, newCleotideName);
  ResidueAddAtomsFromAtomParams(&newResi, pAtomParaSet);
  ResidueAddBondsFromResiTopos(&newResi, pResiTopos);
  //copy coordinates for main-chain atoms
  for (int i = 0;i < ResidueGetAtomCount(&newResi);i++)
  {
    Atom* pAtomNew = ResidueGetAtom(&newResi, i);
    if (pAtomNew->isBBAtom)
    {
      Atom* pAtomOld = ResidueGetAtomByName(pQuery, AtomGetName(pAtomNew));
      if (pAtomOld != NULL)
      {
        pAtomNew->xyz = pAtomOld->xyz;
        pAtomNew->isXyzValid = pAtomOld->isXyzValid;
      }
      else
      {
        pAtomNew->xyz = { 0, 0, 0 };
        pAtomNew->isXyzValid = FALSE;
      }
    }
  }
  // calculate coordinates for side-chain atoms
  if (strcmp(ResidueGetName(pQuery), "DA") == 0 || strcmp(ResidueGetName(pQuery), "DG") == 0 || strcmp(ResidueGetName(pQuery), "A") == 0 || strcmp(ResidueGetName(pQuery), "G") == 0)
  {
    if (strcmp(ResidueGetName(&newResi), "DA") == 0 || strcmp(ResidueGetName(&newResi), "DG") == 0 || strcmp(ResidueGetName(&newResi), "A") == 0 || strcmp(ResidueGetName(&newResi), "G") == 0)
    {
      Atom* pAtomNew = ResidueGetAtomByName(&newResi, "N9");
      Atom* pAtomOld = ResidueGetAtomByName(pQuery, "N9");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;
      pAtomNew = ResidueGetAtomByName(&newResi, "C4");
      pAtomOld = ResidueGetAtomByName(pQuery, "C4");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;
    }
    else if (strcmp(ResidueGetName(&newResi), "DC") == 0 || strcmp(ResidueGetName(&newResi), "DT") == 0 || strcmp(ResidueGetName(&newResi), "C") == 0 || strcmp(ResidueGetName(&newResi), "T") == 0 || strcmp(ResidueGetName(&newResi), "U") == 0)
    {
      Atom* pAtomNew = ResidueGetAtomByName(&newResi, "N1");
      Atom* pAtomOld = ResidueGetAtomByName(pQuery, "N9");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;

      pAtomNew = ResidueGetAtomByName(&newResi, "C2");
      pAtomOld = ResidueGetAtomByName(pQuery, "C4");
      ResidueTopology topNew;
      CharmmIC icNew;
      ResidueTopologyCreate(&topNew);
      CharmmICCreate(&icNew);
      ResiTopoSetGet(pResiTopos, ResidueGetName(&newResi), &topNew);
      ResidueTopologyFindCharmmIC(&topNew, AtomGetName(pAtomNew), &icNew);
      icNew.icParam[2] = GetTorsionAngle(&ResidueGetAtomByName(pQuery, "O4'")->xyz,
        &ResidueGetAtomByName(pQuery, "C1'")->xyz,
        &ResidueGetAtomByName(pQuery, "N9")->xyz,
        &pAtomOld->xyz);
      GetFourthAtom(&ResidueGetAtomByName(&newResi, "O4'")->xyz,
        &ResidueGetAtomByName(&newResi, "C1'")->xyz,
        &ResidueGetAtomByName(&newResi, "N1")->xyz,
        icNew.icParam,
        &pAtomNew->xyz);
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;
    }
    else
    {
      Atom* pAtomNew = ResidueGetAtomByName(&newResi, "N10");
      Atom* pAtomOld = ResidueGetAtomByName(pQuery, "N9");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;

      pAtomNew = ResidueGetAtomByName(&newResi, "C4");
      pAtomOld = ResidueGetAtomByName(pQuery, "C4");
      /*pAtomNew->xyz=pAtomOld->xyz;*/
      ResidueTopology topNew;
      CharmmIC icNew;
      ResidueTopologyCreate(&topNew);
      CharmmICCreate(&icNew);
      ResiTopoSetGet(pResiTopos, ResidueGetName(&newResi), &topNew);
      ResidueTopologyFindCharmmIC(&topNew, AtomGetName(pAtomNew), &icNew);
      icNew.icParam[2] = GetTorsionAngle(&ResidueGetAtomByName(pQuery, "O4'")->xyz,
        &ResidueGetAtomByName(pQuery, "C1'")->xyz,
        &ResidueGetAtomByName(pQuery, "N9")->xyz,
        &pAtomOld->xyz);
      GetFourthAtom(&ResidueGetAtomByName(&newResi, "O4'")->xyz,
        &ResidueGetAtomByName(&newResi, "C1'")->xyz,
        &ResidueGetAtomByName(&newResi, "N10")->xyz,
        icNew.icParam,
        &pAtomNew->xyz);
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;
    }
  }

  else
  {
    if (strcmp(ResidueGetName(&newResi), "DA") == 0 || strcmp(ResidueGetName(&newResi), "DG") == 0 || strcmp(ResidueGetName(&newResi), "A") == 0 || strcmp(ResidueGetName(&newResi), "G") == 0)
    {
      Atom* pAtomNew = ResidueGetAtomByName(&newResi, "N9");
      Atom* pAtomOld = ResidueGetAtomByName(pQuery, "N1");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;

      pAtomNew = ResidueGetAtomByName(&newResi, "C4");
      pAtomOld = ResidueGetAtomByName(pQuery, "C2");
      ResidueTopology topNew;
      CharmmIC icNew;
      ResidueTopologyCreate(&topNew);
      CharmmICCreate(&icNew);
      ResiTopoSetGet(pResiTopos, ResidueGetName(&newResi), &topNew);
      ResidueTopologyFindCharmmIC(&topNew, AtomGetName(pAtomNew), &icNew);
      icNew.icParam[2] = GetTorsionAngle(&ResidueGetAtomByName(pQuery, "O4'")->xyz,
        &ResidueGetAtomByName(pQuery, "C1'")->xyz,
        &ResidueGetAtomByName(pQuery, "N1")->xyz,
        &pAtomOld->xyz);
      GetFourthAtom(&ResidueGetAtomByName(&newResi, "O4'")->xyz,
        &ResidueGetAtomByName(&newResi, "C1'")->xyz,
        &ResidueGetAtomByName(&newResi, "N9")->xyz,
        icNew.icParam,
        &pAtomNew->xyz);
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;
    }
    else if (strcmp(ResidueGetName(&newResi), "DC") == 0 || strcmp(ResidueGetName(&newResi), "DT") == 0 || strcmp(ResidueGetName(&newResi), "C") == 0 || strcmp(ResidueGetName(&newResi), "T") == 0 || strcmp(ResidueGetName(&newResi), "U") == 0)
    {
      Atom* pAtomNew = ResidueGetAtomByName(&newResi, "N1");
      Atom* pAtomOld = ResidueGetAtomByName(pQuery, "N1");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;

      pAtomNew = ResidueGetAtomByName(&newResi, "C2");
      pAtomOld = ResidueGetAtomByName(pQuery, "C2");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;
    }
    else
    {
      Atom* pAtomNew = ResidueGetAtomByName(&newResi, "N10");
      Atom* pAtomOld = ResidueGetAtomByName(pQuery, "N1");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;

      pAtomNew = ResidueGetAtomByName(&newResi, "C4");
      pAtomOld = ResidueGetAtomByName(pQuery, "C2");
      ResidueTopology topNew;
      CharmmIC icNew;
      ResidueTopologyCreate(&topNew);
      CharmmICCreate(&icNew);
      ResiTopoSetGet(pResiTopos, ResidueGetName(&newResi), &topNew);
      ResidueTopologyFindCharmmIC(&topNew, AtomGetName(pAtomNew), &icNew);
      icNew.icParam[2] = GetTorsionAngle(&ResidueGetAtomByName(pQuery, "O4'")->xyz,
        &ResidueGetAtomByName(pQuery, "C1'")->xyz,
        &ResidueGetAtomByName(pQuery, "N1")->xyz,
        &pAtomOld->xyz);
      GetFourthAtom(&ResidueGetAtomByName(&newResi, "O4'")->xyz,
        &ResidueGetAtomByName(&newResi, "C1'")->xyz,
        &ResidueGetAtomByName(&newResi, "N10")->xyz,
        icNew.icParam,
        &pAtomNew->xyz);
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;
    }
  }
  // calculate coordinates for other side-chain atoms
  ResidueCalcAllSidechainAtomXYZ(&newResi, pResiTopos);

  ResidueCopy(pQuery, &newResi);
  ResidueDestroy(&newResi);

  return Success;
}

int StructureCalcAminoAcidPropensityAndRamaEnergy(Structure* pStructure, AAppTable* pAAPP, RamaTable* pRama)
{
  for (int i = 0;i < StructureGetChainCount(pStructure);i++)
  {
    Chain* pChain = StructureGetChain(pStructure, i);
    if (ChainGetType(pChain) == Type_Chain_Protein)
    {
      for (int j = 0;j < ChainGetResidueCount(pChain);j++)
      {
        Residue* pResi = ChainGetResidue(pChain, j);
        AminoAcidPropensityAndRamachandranEnergy(pResi, pAAPP, pRama);
      }
    }
  }
  return Success;
}


int StructureCalcAminoAcidDunbrackEnergy(Structure* pStructure, BBdepRotamerLib* pBBdep)
{
  for (int i = 0;i < StructureGetChainCount(pStructure);i++)
  {
    Chain* pChain = StructureGetChain(pStructure, i);
    if (ChainGetType(pChain) == Type_Chain_Protein)
    {
      for (int j = 0;j < ChainGetResidueCount(pChain);j++)
      {
        Residue* pResi = ChainGetResidue(pChain, j);
        AminoAcidDunbrackEnergy(pResi, pBBdep);
      }
    }
  }
  return Success;
}