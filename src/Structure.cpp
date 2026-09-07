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
#pragma warning(disable : 6301)
#pragma warning(disable : 28182)

#include "Structure.h"
#include <string.h>

int StructureCreate(Structure *pStruct)
{
  strcpy(pStruct->name, "");
  pStruct->chnNum = 0;
  pStruct->chains = NULL;
  pStruct->desSiteNum = 0;
  pStruct->desSites = NULL;
  return Success;
}
int StructureDestroy(Structure *pStruct)
{
  for (int i = 0; i < pStruct->desSiteNum; i++)
  {
    DesignSiteDestroy(&pStruct->desSites[i]);
  }
  pStruct->desSites = NULL;
  pStruct->desSiteNum = 0;
  for (int i = 0; i < pStruct->chnNum; i++)
  {
    ChainDestroy(&pStruct->chains[i]);
  }
  free(pStruct->chains);
  strcpy(pStruct->name, "");
  pStruct->chnNum = 0;
  pStruct->chains = NULL;
  return Success;
}

char *StructureGetName(Structure *pStruct)
{
  return pStruct->name;
}

int StructureSetName(Structure *pStruct, char *name)
{
  if (strlen(name) > MAX_LEN_STRUCTURE_NAME)
  {
    char errMsg[MAX_LEN_ERR_MSG + 1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s line %d, name is too long", __FILE__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pStruct->name, name);
  return Success;
}

int StructureGetChainNum(Structure *pStruct)
{
  return pStruct->chnNum;
}

Chain *StructureGetChain(Structure *pStruct, int index)
{
  if (index < 0 || index >= StructureGetChainNum(pStruct))
  {
    return NULL;
  }
  return &pStruct->chains[index];
}

Chain *StructureFindChainByName(Structure *pStruct, char *chainName)
{
  int index = -1;
  int result = StructureFindChainIndex(pStruct, chainName, &index);
  if (FAILED(result) == TRUE || index == -1)
  {
    return NULL;
  }
  return StructureGetChain(pStruct, index);
}

int StructureFindChainIndex(Structure *pStruct, char *chnName, int *ndx)
{
  int i;
  for (i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    if (strcmp(ChainGetName(&pStruct->chains[i]), chnName) == 0)
    {
      *ndx = i;
      return Success;
    }
  }
  return DataNotExistError;
}

int StructureAddChain(Structure *pStruct, Chain *pChain)
{
  int ndx = -1;
  int result = StructureFindChainIndex(pStruct, ChainGetName(pChain), &ndx);
  if (FAILED(result))
  {
    (pStruct->chnNum)++;
    pStruct->chains = (Chain *)realloc(pStruct->chains, sizeof(Chain) * pStruct->chnNum);
    ChainCreate(&pStruct->chains[pStruct->chnNum - 1]);
    return ChainCopy(&pStruct->chains[pStruct->chnNum - 1], pChain);
  }
  return ChainCopy(&pStruct->chains[ndx], pChain);
}

int StructureDeleteChain(Structure *pStruct, char *chnName)
{
  int ndx;
  int result = StructureFindChainIndex(pStruct, chnName, &ndx);
  if (FAILED(result))
  {
    return result;
  }
  for (int i = ndx; i < pStruct->chnNum - 1; i++)
  {
    ChainCopy(&pStruct->chains[i], &pStruct->chains[i + 1]);
  }
  ChainDestroy(&pStruct->chains[pStruct->chnNum - 1]);
  (pStruct->chnNum)--;
  return Success;
}

int StructureShowInPDBFormat(Structure *pStruct, FILE *pFile)
{
  int atomNdx = 1;
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChn = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChn); j++)
    {
      Residue *pRes = ChainGetResidue(pChn, j);
      ResidueShowInPDBFormat(pRes, "ATOM", ResidueGetChainName(pRes), atomNdx, ResidueGetPosInChain(pRes), pFile);
      atomNdx += ResidueGetAtomNum(pRes);
    }
  }
  return Success;
}

int StructureReadPDB(Structure *pStruct, char *pdbfile, AtomParams *pParams, ResTopos *pTops)
{
  char iniChainID[MAX_LEN_CHAIN_NAME + 1] = "UNK";
  char iniResPos[6] = "UNK";
  BOOL firstResInChn = TRUE;

  FileReader fr;
  if (FAILED(FileReaderCreate(&fr, pdbfile)))
  {
    char errMsg[MAX_LEN_ERR_MSG + 1];
    sprintf(errMsg, "in file %s line %d, cannot read file %s", __FILE__, __LINE__, pdbfile);
    TraceError(errMsg, IOError);
    return IOError;
  }

  char line[MAX_LEN_ONE_LINE_CONTENT + 1];
  while (!FAILED(FileReaderGetNextLine(&fr, line)))
  {
    char keyword[6] = "";
    ExtractTargetStringFromSourceString(keyword, line, 0, 4);

    // Only read ATOM records.
    if (strcmp(keyword, "ATOM") != 0)
    {
      continue;
    }

    char strResName[MAX_LEN_RES_NAME + 1] = "";
    ExtractTargetStringFromSourceString(strResName, line, 17, 4);
    Type_Chain chainType = ChainTypeIdentifiedFromResidueName(strResName);
    if (chainType == Type_Chain_Unknown)
    {
      printf("residue name %s cannot be regonized, ignored it.\n", strResName);
      continue;
    }
    char strAtomName[MAX_LEN_ATOM_NAME + 1] = "";
    ExtractTargetStringFromSourceString(strAtomName, line, 12, 4);

    char strChainID[MAX_LEN_CHAIN_NAME + 1] = "";
    ExtractTargetStringFromSourceString(strChainID, line, 21, 1);
    if (strcmp(strChainID, "") == 0 || strcmp(strChainID, " ") == 0)
    {
      strcpy(strChainID, "A");
    }

    char strResPos[6] = "";
    ExtractTargetStringFromSourceString(strResPos, line, 22, 5);
    // Start a new chain.
    if (strcmp(iniChainID, strChainID) != 0)
    {
      // Finish processing the previous chain.
      Chain *pPreviousChain = StructureGetChain(pStruct, StructureGetChainNum(pStruct) - 1);
      if (pPreviousChain != NULL && ChainGetResidueNum(pPreviousChain) > 0)
      {
        // Patch the C-terminus of the previous protein chain.
        if (ChainGetType(pPreviousChain) == Type_Chain_Protein)
        {
          Residue *pFirstResInChn = ChainGetResidue(pPreviousChain, 0);
          Residue *pLastResInChn = ChainGetResidue(pPreviousChain, ChainGetResidueNum(pPreviousChain) - 1);
          if (pFirstResInChn != NULL &&
              pLastResInChn != NULL &&
              (ResidueGetAtomByName(pFirstResInChn, "HT1") 
              || ResidueGetAtomByName(pFirstResInChn, "HN1")))
          {
            ResiduePatchCTER(pLastResInChn, "CTER", pParams, pTops);
            pLastResInChn->terType = Type_ResIsCter;
          }
        }
        for (int i = 0; i < ChainGetResidueNum(pPreviousChain); i++)
        {
          Residue *pRes = ChainGetResidue(pPreviousChain, i);
          if (pRes->isSCIntact == TRUE)
          {
            ResidueCalcAllAtomXYZ(pRes, pTops, ChainGetResidue(pPreviousChain, i - 1), ChainGetResidue(pPreviousChain, i + 1));
            ResidueCalcSidechainTorsion(pRes, pTops);
          }
          else
          {
            ResidueCalcAllBackboneAtomXYZ(pRes, pTops, ChainGetResidue(pPreviousChain, i - 1), ChainGetResidue(pPreviousChain, i + 1));
          }
        }
      }
      
      //Initialize the new chain.
      strcpy(iniChainID, strChainID);

      // Important: residue state belongs to the previous chain,
      // so it must be reset when a new chain begins.
      strcpy(iniResPos, "UNK");
      firstResInChn = TRUE;

      Chain chain;
      ChainCreate(&chain);
      ChainSetType(&chain, chainType);
      ChainSetName(&chain, strChainID);
      StructureAddChain(pStruct, &chain);
      ChainDestroy(&chain);

      // Rewind so the current line is processed again as the
      // first residue of the newly created chain.
      FileReaderSetCurrentPos(&fr, FileReaderGetCurrentPos(&fr) - 1);
      continue;
    }

    // Start of a new residue in the current chain.
    if (strcmp(iniResPos, strResPos) != 0)
    {
      Residue res;
      ResidueCreate(&res);

      // Determine the protonation state of histidine.
      if (strcmp(strResName, "HIS") == 0)
      {
        FileReaderSetCurrentPos(&fr, FileReaderGetCurrentPos(&fr) - 1);
        int hisState = 1;
        PDBReaderCheckHisState(&fr, &hisState);
        if (hisState == 1)
        {
          strcpy(strResName, "HSD");
        }
        else if (hisState == 2)
        {
          strcpy(strResName, "HSE");
        }
        else
        {
          strcpy(strResName, "HSP");
        }
      }

      ResidueSetName(&res, strResName);
      ResidueSetChainName(&res, strChainID);

      // atoi() reads the numeric part and ignores a possible
      // insertion code. This preserves the original behavior.
      ResidueSetPosInChain(&res, atoi(strResPos));
      ResidueAddAtomsFromAtomParams(&res, pParams);
      ResidueAddBondsFromResiTopos(&res, pTops);

      Chain *pCurrentChain = StructureFindChainByName(pStruct, strChainID);
      if (pCurrentChain == NULL)
      {
        char errMsg[MAX_LEN_ERR_MSG + 1];
        sprintf(errMsg, "in file %s line %d, cannot find the target chain in structure", __FILE__, __LINE__);
        TraceError(errMsg, DataNotExistError);
        ResidueDestroy(&res);
        FileReaderDestroy(&fr);
        return DataNotExistError;
      }

      // Patch the N-terminus only for the first residue of
      // a protein chain.
      if (firstResInChn == TRUE && 
          ChainGetType(pCurrentChain) == Type_Chain_Protein)
      {
        ResiduePatchNTERorCTER(&res, "NTER", pParams, pTops);
        res.terType = Type_ResIsNter;
      }
      // firstResInChn must be cleared for every chain type,
      // not only protein chains.
      firstResInChn = FALSE;

      // Rewind so ResidueReadXYZFromPDB() starts from the
      // first atom record of this residue.
      FileReaderSetCurrentPos(&fr, FileReaderGetCurrentPos(&fr) - 1);
      ResidueReadXYZFromPDB(&res, &fr);
      ResidueCheckAtomCoordinateValidity(&res);
      ChainAppendResidue(pCurrentChain, &res);
      ResidueDestroy(&res);
      strcpy(iniResPos, strResPos);
    }
  }

  // Finish processing the last polymer chain.
  Chain *pLastChain = StructureGetChain(pStruct, StructureGetChainNum(pStruct) - 1);
  if (pLastChain != NULL && ChainGetResidueNum(pLastChain) > 0)
  {
    if (ChainGetType(pLastChain) == Type_Chain_Protein)
    {
      Residue *pFirstResInChn = ChainGetResidue(pLastChain, 0);
      Residue *pLastResInChn = ChainGetResidue(pLastChain, ChainGetResidueNum(pLastChain) - 1);
      if (pFirstResInChn != NULL &&
          pLastResInChn != NULL &&
          (ResidueGetAtomByName(pFirstResInChn, "HT1") != NULL 
          || ResidueGetAtomByName(pFirstResInChn, "HN1") != NULL))
      {
        ResiduePatchCTER(pLastResInChn, "CTER", pParams, pTops);
        pLastResInChn->terType = Type_ResIsCter;
      }
    }

    for (int i = 0; i < ChainGetResidueNum(pLastChain); i++)
    {
      Residue *pRes = ChainGetResidue(pLastChain, i);
      if (pRes->isSCIntact == TRUE)
      {
        ResidueCalcAllAtomXYZ(pRes, pTops, ChainGetResidue(pLastChain, i - 1), ChainGetResidue(pLastChain, i + 1));
        ResidueCalcSidechainTorsion(pRes, pTops);
      }
      else
      {
        ResidueCalcAllBackboneAtomXYZ(pRes, pTops, ChainGetResidue(pLastChain, i - 1), ChainGetResidue(pLastChain, i + 1));
      }
    }
  }

  FileReaderDestroy(&fr);
  return Success;
}

int StructureReadMol2(Structure *pStruct, char *mol2file, AtomParams *pParams, ResTopos *pTops)
{
  FILE *pMol2 = fopen(mol2file, "r");

  char line[MAX_LEN_ONE_LINE_CONTENT + 1];
  char keyword[MAX_LEN_ONE_LINE_CONTENT + 1];
  BOOL readingAtom = FALSE;
  BOOL readingBond = FALSE;

  char resiName[MAX_LEN_RES_NAME + 1] = "#";

  // add a new chain
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
        if (strlen(resiName) > 3)
          resiName[3] = '\0';
        if (strcmp(resiName, "***") == 0)
          strcpy(resiName, "LIG");
        if (LigandResidueNameConflictWithAminoAcid(resiName))
          strcpy(resiName, "LIG");
      }
      // deploy the residue first if it's not ready
      if (ResidueGetAtomNum(&newResi) == 0)
      {
        ResidueSetName(&newResi, resiName);
        ResidueSetPosInChain(&newResi, posInChain);
        ResidueAddAtomsFromAtomParams(&newResi, pParams);
        ResidueAddBondsFromResiTopos(&newResi, pTops);
        // set chain name
        // ChainSetName(&newChain, strChainID);
        ChainSetName(&newChain, "X");
      }
      // Atom* pAtom=ResidueGetAtomByName(&newResi,strAtomName);
      Atom *pAtom = ResidueGetAtom(&newResi, atomNdx - 1);
      pAtom->xyz.X = x;
      pAtom->xyz.Y = y;
      pAtom->xyz.Z = z;
      pAtom->isXyzValid = TRUE;
    }
  }
  ChainAppendResidue(&newChain, &newResi);
  StructureAddChain(pStruct, &newChain);
  ChainDestroy(&newChain);
  ResidueDestroy(&newResi);
  fclose(pMol2);

  return Success;
}

int ChainComputeResiduePosition(Structure *pStructure, int chainIndex)
{
  Chain *pChainI = StructureGetChain(pStructure, chainIndex);
  for (int ir = 0; ir < ChainGetResidueNum(pChainI); ir++)
  {
    Residue *pResiIR = ChainGetResidue(pChainI, ir);
    pResiIR->nCbIn10A = 0;
  }
  for (int ir = 0; ir < ChainGetResidueNum(pChainI); ir++)
  {
    Residue *pResiIR = ChainGetResidue(pChainI, ir);
    Atom *pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CB");
    if (pAtomCAorCB1 == NULL)
      pAtomCAorCB1 = ResidueGetAtomByName(pResiIR, "CA");
    for (int is = ir + 1; is < ChainGetResidueNum(pChainI); is++)
    {
      Residue *pResiIS = ChainGetResidue(pChainI, is);
      Atom *pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CB");
      if (pAtomCAorCB2 == NULL)
        pAtomCAorCB2 = ResidueGetAtomByName(pResiIS, "CA");
      if (XYZDist(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 10.0)
      {
        pResiIR->nCbIn10A++;
        pResiIS->nCbIn10A++;
      }
    }
    // printf("Residue: %s %d %s, NumCBwithin8AtoCurResi: %d\n",ResidueGetChainName(pResiIR), ResidueGetPosInChain(pResiIR), ResidueGetName(pResiIR), pResiIR->numCBwithin8AtoCurResidue);
  }

  return Success;
}

int StructureComputeResiduePosition(Structure *pStructure)
{
  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChainI = StructureGetChain(pStructure, i);
    for (int ir = 0; ir < ChainGetResidueNum(pChainI); ir++)
    {
      Residue *pResiIR = ChainGetResidue(pChainI, ir);
      pResiIR->nCbIn10A = 0;
    }
  }
  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein)
      continue;
    for (int ir = 0; ir < ChainGetResidueNum(pChainI); ir++)
    {
      Residue *pResIR = ChainGetResidue(pChainI, ir);
      Atom *pAtomCAorCB1 = ResidueGetAtomByName(pResIR, "CB");
      if (pAtomCAorCB1 == NULL)
        pAtomCAorCB1 = ResidueGetAtomByName(pResIR, "CA");
      for (int is = ir + 1; is < ChainGetResidueNum(pChainI); is++)
      {
        Residue *pResIS = ChainGetResidue(pChainI, is);
        Atom *pAtomCAorCB2 = ResidueGetAtomByName(pResIS, "CB");
        if (pAtomCAorCB2 == NULL)
          pAtomCAorCB2 = ResidueGetAtomByName(pResIS, "CA");
        if (XYZDist(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 10.0)
        {
          pResIR->nCbIn10A++;
          pResIS->nCbIn10A++;
        }
      }
      for (int k = i + 1; k < StructureGetChainNum(pStructure); k++)
      {
        Chain *pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) != Type_Chain_Protein)
          continue;
        for (int ks = 0; ks < ChainGetResidueNum(pChainK); ks++)
        {
          Residue *pResKS = ChainGetResidue(pChainK, ks);
          Atom *pAtomCAorCB2 = pAtomCAorCB2 = ResidueGetAtomByName(pResKS, "CB");
          if (pAtomCAorCB2 == NULL)
            pAtomCAorCB2 = ResidueGetAtomByName(pResKS, "CA");
          if (XYZDist(&pAtomCAorCB1->xyz, &pAtomCAorCB2->xyz) < 10.0)
          {
            pResIR->nCbIn10A++;
            pResKS->nCbIn10A++;
          }
        }
      }
      // printf("Residue: %s %d %s, NumCBwithin8AtoCurResi: %d\n",ResidueGetChainName(pResiIR), ResidueGetPosInChain(pResiIR), ResidueGetName(pResiIR), pResiIR->numCBwithin8AtoCurResidue);
    }
  }

  return Success;
}

int StructureShowAtomParameter(Structure *pStruct)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    ChainShowAtomParameter(StructureGetChain(pStruct, i));
  }
  return Success;
}

int StructureShowBondInformation(Structure *pStruct)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    ChainShowBondInformation(StructureGetChain(pStruct, i));
  }
  return Success;
}

int StructureCheckIntraBondType(Structure *pStruct)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChain = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChain); j++)
    {
      Residue *pRes = ChainGetResidue(pChain, j);
      for (int k = 0; k < ResidueGetAtomNum(pRes); k++)
      {
        Atom *pAtom = ResidueGetAtom(pRes, k);
        for (int m = k + 1; m < ResidueGetAtomNum(pRes); m++)
        {
          Atom *pAtom2 = ResidueGetAtom(pRes, m);
          int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtom), AtomGetName(pAtom2), ResidueGetBonds(pRes));
          printf("Residue: %s %d %s, Atom1: %s, Atom2: %s, BondType: %d\n",
            ResidueGetName(pRes), ResidueGetPosInChain(pRes), ResidueGetChainName(pRes),AtomGetName(pAtom),AtomGetName(pAtom2),bondType);
        }
      }
    }
  }
  return Success;
}

int StructureCheckNeighbouringBondType(Structure *pStruct)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChain = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChain) - 1; j++)
    {
      Residue *pRes1 = ChainGetResidue(pChain, j);
      Residue *pRes2 = ChainGetResidue(pChain, j + 1);
      for (int p = 0; p < pRes1->atoms.atomNum; p++)
      {
        Atom *pAtom1 = ResidueGetAtom(pRes1, p);
        for (int q = 0; q < pRes2->atoms.atomNum; q++)
        {
          Atom *pAtom2 = ResidueGetAtom(pRes2, q);
          int bondType = ResidueAndNextResidueInterBondConnectionCheck_charmm19(pAtom1->name, pAtom2->name, pRes2->name);
          printf("Residue: %s %d %s, Atom1: %s, Residue: %s %d %s, Atom2: %s, BondType: %d\n", ResidueGetName(pRes1), ResidueGetPosInChain(pRes1), ResidueGetChainName(pRes1), AtomGetName(pAtom1), ResidueGetName(pRes2), ResidueGetPosInChain(pRes2), ResidueGetChainName(pRes2), AtomGetName(pAtom2), bondType);
        }
      }
    }
  }

  return Success;
}


int StructureGetDesignSiteNum(Structure *pStruct)
{
  return pStruct->desSiteNum;
}

DesignSite *StructureGetDesignSite(Structure *pStruct, int index)
{
  if (index < 0 || index >= pStruct->desSiteNum)
    return NULL;
  else
    return &pStruct->desSites[index];
}

DesignSite *StructureFindDesignSite(Structure *pStruct, int chainIndex, int resIndex)
{
  if (chainIndex < 0 || chainIndex >= StructureGetChainNum(pStruct))
    return NULL;
  Chain *pChain = StructureGetChain(pStruct, chainIndex);
  if (resIndex < 0 || resIndex >= ChainGetResidueNum(pChain))
    return NULL;
  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++)
  {
    DesignSite *pSite = StructureGetDesignSite(pStruct, i);
    if (pSite->chnNdx == chainIndex && pSite->resNdx == resIndex)
    {
      return pSite;
    }
  }
  return NULL;
}

DesignSite *StructureFindDesignSiteByChainName(Structure *pStruct, char *chainName, int posInChain)
{
  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++)
  {
    DesignSite *pSite = StructureGetDesignSite(pStruct, i);
    if (strcmp(chainName, DesignSiteGetChainName(pSite)) == 0 &&
        DesignSiteGetPosInChain(pSite) == posInChain)
    {
      return pSite;
    }
  }
  return NULL;
}

int StructureFindDesignSiteIndexByChainNameAndPosInChain(Structure *pStruct, char *chainName, int posInChain)
{
  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++)
  {
    DesignSite *pDesignSite = StructureGetDesignSite(pStruct, i);
    if (strcmp(DesignSiteGetChainName(pDesignSite), chainName) == 0 &&
        DesignSiteGetPosInChain(pDesignSite) == posInChain)
    {
      return i;
    }
  }
  return -1;
}

int StructureGetDesignSiteIndex(Structure *pStruct, int chainIndex, int resIndex)
{
  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++)
  {
    DesignSite *pDesignSite = StructureGetDesignSite(pStruct, i);
    if (pDesignSite->chnNdx == chainIndex && pDesignSite->resNdx == resIndex)
    {
      return i;
    }
  }
  return -1;
}

int StructureShowDesignSites(Structure *pStruct)
{
  double confSpace = 0.0;
  int totRotNum = 0;
  for (int i = 0; i < pStruct->desSiteNum; i++)
  {
    DesignSite *pSite = StructureGetDesignSite(pStruct, i);
    int rotNum = RotamersGetNum(DesignSiteGetRotamers(pSite));
    printf("site %3d: %s %4d %3s, %6d rotamers, ",
           i, ResidueGetChainName(pSite->pRes), ResidueGetPosInChain(pSite->pRes), ResidueGetName(pSite->pRes), rotNum);
    totRotNum += rotNum;
    if (rotNum > 0)
      confSpace += log((double)rotNum) / log(10.0);
    switch (pSite->pRes->desType)
    {
    case Type_DesType_Catalytic:
      printf("Catalytic\n");
      break;
    case Type_DesType_Fixed:
      printf("Fixed\n");
      break;
    case Type_DesType_Mutable:
      printf("Mutable\n");
      break;
    case Type_DesType_SmallMol:
      printf("Smallmol\n");
      break;
    case Type_DesType_Repackable:
      printf("Repackable\n");
      break;
    case Type_DesType_NatRot:
      printf("Fixed_Natrot\n");
      break;
    default:
      break;
    }
  }
  printf("number of total rotamers: %d, conformation space: 1e^%.1f\n", totRotNum, confSpace);
  return Success;
}

int ProteinSiteAddDesignSite(Structure *pStruct, int chainIndex, int resIndex)
{
  int result = Success;
  DesignSite *pSite = StructureFindDesignSite(pStruct, chainIndex, resIndex);
  if (pSite != NULL)
  {
    return Success;
  }
  (pStruct->desSiteNum)++;
  pStruct->desSites = (DesignSite *)realloc(pStruct->desSites, sizeof(DesignSite) * pStruct->desSiteNum);
  DesignSiteCreate(&pStruct->desSites[pStruct->desSiteNum - 1]);
  pSite = StructureGetDesignSite(pStruct, pStruct->desSiteNum - 1);
  pSite->pRes = ChainGetResidue(StructureGetChain(pStruct, chainIndex), resIndex);
  pSite->chnNdx = chainIndex;
  pSite->resNdx = resIndex;
  return result;
}

int StructureRemoveAllDesignSites(Structure *pStruct)
{
  for (int i = 0; i < pStruct->desSiteNum; i++)
  {
    DesignSite *pSite = &pStruct->desSites[i];
    DesignSiteDestroy(pSite);
    pSite = NULL;
  }
  free(pStruct->desSites);
  pStruct->desSites = NULL;
  pStruct->desSiteNum = 0;
  return Success;
}

int ProteinSiteRemoveDesignSite(Structure *pStruct, int chainIndex, int resIndex)
{
  DesignSite *pSite = StructureFindDesignSite(pStruct, chainIndex, resIndex);
  if (pSite != NULL)
  {
    DesignSiteDestroy(pSite);
    pSite = NULL;
  }
  (pStruct->desSiteNum)--;
  return Success;
}

int StructureFindSmallMol(Structure *pStruct, Residue **ppSmallMol)
{
  int result = DataNotExistError;
  for (int i = 0; i < pStruct->chnNum; i++)
  {
    if (pStruct->chains[i].type == Type_Chain_SmallMol)
    {
      *ppSmallMol = ChainGetResidue(&pStruct->chains[i], 0);
      if (*ppSmallMol != NULL)
      {
        result = Success;
        break;
      }
    }
  }
  return result;
}

int StructureCalcProteinResidueSidechainTorsion(Structure *pStruct, ResTopos *pTops)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChain = StructureGetChain(pStruct, i);
    if (ChainGetType(pChain) == Type_Chain_Protein)
    {
      for (int j = 0; j < ChainGetResidueNum(pChain); j++)
      {
        Residue *pResidue = ChainGetResidue(pChain, j);
        if (strcmp(ResidueGetName(pResidue), "ALA") == 0 || strcmp(ResidueGetName(pResidue), "GLY") == 0)
          continue;
        ResidueCalcSidechainTorsion(pResidue, pTops);
      }
    }
  }
  return Success;
}

int StructureCalcPhiPsi(Structure *pStruct)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChain = StructureGetChain(pStruct, i);
    if (ChainGetType(pChain) != Type_Chain_Protein)
      continue;
    for (int j = 0; j < ChainGetResidueNum(pChain); j++)
    {
      Residue *pRes0 = ChainGetResidue(pChain, j - 1);
      Residue *pRes1 = ChainGetResidue(pChain, j);
      Residue *pRes2 = ChainGetResidue(pChain, j + 1);

      if (pRes0 != NULL && AA3GetIndex(ResidueGetName(pRes0)) < 0)
      {
        pRes0 = NULL;
      }
      if (pRes1 != NULL && AA3GetIndex(ResidueGetName(pRes1)) < 0)
      {
        pRes1 = NULL;
        continue;
      }
      if (pRes2 != NULL && AA3GetIndex(ResidueGetName(pRes2)) < 0)
      {
        pRes2 = NULL;
      }

      Atom *pC0 = NULL;
      Atom *pN1 = ResidueGetAtomByName(pRes1, "N");
      Atom *pCA1 = ResidueGetAtomByName(pRes1, "CA");
      Atom *pC1 = ResidueGetAtomByName(pRes1, "C");
      Atom *pN2 = NULL;

      double phi = 0.0;
      double psi = 0.0;
      if (pRes0 == NULL)
      {
        phi = -60.0;
        pC0 = NULL;
      }
      else
      {
        pC0 = ResidueGetAtomByName(pRes0, "C");
      }

      if (pRes2 == NULL)
      {
        psi = 60.0;
        pN2 = NULL;
      }
      else
      {
        pN2 = ResidueGetAtomByName(pRes2, "N");
      }

      if (pC0 != NULL)
      {
        double dist_C0N1 = XYZDist(&pC0->xyz, &pN1->xyz);
        if (dist_C0N1 < 1.45 && dist_C0N1 > 1.25)
        { // indicative of an amide bond
          phi = RadToDeg(GetTorsionAngle(&pC0->xyz, &pN1->xyz, &pCA1->xyz, &pC1->xyz));
        }
        else
        {
          phi = -60.0;
        }
      }

      if (pN2 != NULL)
      {
        double dist_C1N2 = XYZDist(&pC1->xyz, &pN2->xyz);
        if (dist_C1N2 < 1.45 && dist_C1N2 > 1.25)
        { // indicative of an amide bond
          psi = RadToDeg(GetTorsionAngle(&pN1->xyz, &pCA1->xyz, &pC1->xyz, &pN2->xyz));
        }
        else
        {
          psi = 60.0;
        }
      }

      // assign phi and psi values to residue1
      pRes1->phipsi[0] = phi;
      pRes1->phipsi[1] = psi;
    }
  }

  return Success;
}

int StructureCopy(Structure *pStruct1, Structure *pStruct2)
{
  for (int i = 0; i < StructureGetChainNum(pStruct2); i++)
  {
    StructureAddChain(pStruct1, &pStruct2->chains[i]);
  }
  pStruct1->chnNum = pStruct2->chnNum;
  pStruct1->desSiteNum = pStruct2->desSiteNum;
  strcpy(pStruct1->name, pStruct2->name);
  for (int i = 0; i < pStruct2->desSiteNum; i++)
  {
    DesignSiteCopy(&pStruct1->desSites[i], &pStruct2->desSites[i]);
  }
  return Success;
}


Residue *StructureFindNucleotidePair(Structure *pStruct, Residue *pRes)
{
  Chain *pChainI = StructureFindChainByName(pStruct, ResidueGetChainName(pRes));
  if (ChainGetType(pChainI) == Type_Chain_DNA)
  {
    if (!strcmp(ResidueGetName(pRes), "DA"))
    {
      // paired residue is DT or U;
      for (int k = 0; k < StructureGetChainNum(pStruct); k++)
      {
        Chain *pChainK = StructureGetChain(pStruct, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            // paired residue must be DT;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pRes) - ResidueGetPosInChain(pResS)) <= 1)
                continue;
              if (!strcmp(ResidueGetName(pResS), "DT"))
              {
                // check the distance
                // atom N1 and N6 from nucleotide DA
                // atom N3 and O4 from nucleotide DT
                double distN1N3 = XYZDist(&ResidueGetAtomByName(pRes, "N1")->xyz, &ResidueGetAtomByName(pResS, "N3")->xyz);
                double distN6O4 = XYZDist(&ResidueGetAtomByName(pRes, "N6")->xyz, &ResidueGetAtomByName(pResS, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
          else
          {
            // paired residue must be DT or U;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResS), "DT") == 0 || strcmp(ResidueGetName(pResS), "U") == 0)
              {
                // check the distance
                // atom N1 and N6 from nucleotide DA
                // atom N3 and O4 from nucleotide DT or U
                double distN1N3 = XYZDist(&ResidueGetAtomByName(pRes, "N1")->xyz, &ResidueGetAtomByName(pResS, "N3")->xyz);
                double distN6O4 = XYZDist(&ResidueGetAtomByName(pRes, "N6")->xyz, &ResidueGetAtomByName(pResS, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
        }
      }
    }
    else if (!strcmp(ResidueGetName(pRes), "DC"))
    {
      // paired residue is DG or G;
      for (int k = 0; k < StructureGetChainNum(pStruct); k++)
      {
        Chain *pChainK = StructureGetChain(pStruct, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            // paired residue must be DG;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pRes) - ResidueGetPosInChain(pResS)) <= 1)
                continue;
              if (!strcmp(ResidueGetName(pResS), "DG"))
              {
                // check the distance
                // atom O2, N3 and N4 from nucleotide DC
                // atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDist(&ResidueGetAtomByName(pRes, "O2")->xyz, &ResidueGetAtomByName(pResS, "N2")->xyz);
                double distN3N1 = XYZDist(&ResidueGetAtomByName(pRes, "N3")->xyz, &ResidueGetAtomByName(pResS, "N1")->xyz);
                double distN4O6 = XYZDist(&ResidueGetAtomByName(pRes, "N4")->xyz, &ResidueGetAtomByName(pResS, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
          else
          {
            // paired residue must be DG or G;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResS), "DG") == 0 || strcmp(ResidueGetName(pResS), "G") == 0)
              {
                // check the distance
                // atom O2, N3 and N4 from nucleotide DC
                // atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDist(&ResidueGetAtomByName(pRes, "O2")->xyz, &ResidueGetAtomByName(pResS, "N2")->xyz);
                double distN3N1 = XYZDist(&ResidueGetAtomByName(pRes, "N3")->xyz, &ResidueGetAtomByName(pResS, "N1")->xyz);
                double distN4O6 = XYZDist(&ResidueGetAtomByName(pRes, "N4")->xyz, &ResidueGetAtomByName(pResS, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
        }
      }
    }
    else if (!strcmp(ResidueGetName(pRes), "DG"))
    {
      // paired residue is DC or C;
      for (int k = 0; k < StructureGetChainNum(pStruct); k++)
      {
        Chain *pChainK = StructureGetChain(pStruct, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            // paired residue must be DC;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pRes) - ResidueGetPosInChain(pResS)) <= 1)
                continue;
              if (!strcmp(ResidueGetName(pResS), "DC"))
              {
                // check the distance
                // atom O2, N3 and N4 from nucleotide DC
                // atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDist(&ResidueGetAtomByName(pResS, "O2")->xyz, &ResidueGetAtomByName(pRes, "N2")->xyz);
                double distN3N1 = XYZDist(&ResidueGetAtomByName(pResS, "N3")->xyz, &ResidueGetAtomByName(pRes, "N1")->xyz);
                double distN4O6 = XYZDist(&ResidueGetAtomByName(pResS, "N4")->xyz, &ResidueGetAtomByName(pRes, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
          else
          {
            // paired residue must be DC or C;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResS), "DC") == 0 || strcmp(ResidueGetName(pResS), "C") == 0)
              {
                // check the distance
                // atom O2, N3 and N4 from nucleotide DC
                // atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDist(&ResidueGetAtomByName(pResS, "O2")->xyz, &ResidueGetAtomByName(pRes, "N2")->xyz);
                double distN3N1 = XYZDist(&ResidueGetAtomByName(pResS, "N3")->xyz, &ResidueGetAtomByName(pRes, "N1")->xyz);
                double distN4O6 = XYZDist(&ResidueGetAtomByName(pResS, "N4")->xyz, &ResidueGetAtomByName(pRes, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
        }
      }
    }
    else
    { //! strcmp(ResidueGetName(pQuery),"DT")
      // paired residue is DA or A;
      for (int k = 0; k < StructureGetChainNum(pStruct); k++)
      {
        Chain *pChainK = StructureGetChain(pStruct, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            // paired residue must be DA;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pRes) - ResidueGetPosInChain(pResS)) <= 1)
                continue;
              if (!strcmp(ResidueGetName(pResS), "DA"))
              {
                // check the distance
                // atom N1 and N6 from nucleotide DA or A
                // atom N3 and O4 from nucleotide DT
                double distN1N3 = XYZDist(&ResidueGetAtomByName(pResS, "N1")->xyz, &ResidueGetAtomByName(pRes, "N3")->xyz);
                double distN6O4 = XYZDist(&ResidueGetAtomByName(pResS, "N6")->xyz, &ResidueGetAtomByName(pRes, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
          else
          {
            // paired residue must be DA or A;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResS), "DA") == 0 || strcmp(ResidueGetName(pResS), "A") == 0)
              {
                // check the distance
                // atom N1 and N6 from nucleotide DA or A
                // atom N3 and O4 from nucleotide DT
                double distN1N3 = XYZDist(&ResidueGetAtomByName(pResS, "N1")->xyz, &ResidueGetAtomByName(pRes, "N3")->xyz);
                double distN6O4 = XYZDist(&ResidueGetAtomByName(pResS, "N6")->xyz, &ResidueGetAtomByName(pRes, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResS;
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
    if (!strcmp(ResidueGetName(pRes), "A"))
    {
      // paired residue is DT or U;
      for (int k = 0; k < StructureGetChainNum(pStruct); k++)
      {
        Chain *pChainK = StructureGetChain(pStruct, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            // paired residue must be U;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pRes) - ResidueGetPosInChain(pResS)) <= 1)
                continue;
              if (!strcmp(ResidueGetName(pResS), "U"))
              {
                // check the distance
                // atom N1 and N6 from nucleotide A
                // atom N3 and O4 from nucleotide U
                double distN1N3 = XYZDist(&ResidueGetAtomByName(pRes, "N1")->xyz, &ResidueGetAtomByName(pResS, "N3")->xyz);
                double distN6O4 = XYZDist(&ResidueGetAtomByName(pRes, "N6")->xyz, &ResidueGetAtomByName(pResS, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
          else
          {
            // paired residue must be DT or U;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResS), "DT") == 0 || strcmp(ResidueGetName(pResS), "U") == 0)
              {
                // check the distance
                // atom N1 and N6 from nucleotide A
                // atom N3 and O4 from nucleotide DT or U
                double distN1N3 = XYZDist(&ResidueGetAtomByName(pRes, "N1")->xyz, &ResidueGetAtomByName(pResS, "N3")->xyz);
                double distN6O4 = XYZDist(&ResidueGetAtomByName(pRes, "N6")->xyz, &ResidueGetAtomByName(pResS, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
        }
      }
    }
    else if (!strcmp(ResidueGetName(pRes), "C"))
    {
      // paired residue is DG or G;
      for (int k = 0; k < StructureGetChainNum(pStruct); k++)
      {
        Chain *pChainK = StructureGetChain(pStruct, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            // paired residue must be G;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pRes) - ResidueGetPosInChain(pResS)) <= 1)
                continue;
              if (!strcmp(ResidueGetName(pResS), "G"))
              {
                // check the distance
                // atom O2, N3 and N4 from nucleotide C
                // atom N2, N1 and O6 from nucleotide G
                double distO2N2 = XYZDist(&ResidueGetAtomByName(pRes, "O2")->xyz, &ResidueGetAtomByName(pResS, "N2")->xyz);
                double distN3N1 = XYZDist(&ResidueGetAtomByName(pRes, "N3")->xyz, &ResidueGetAtomByName(pResS, "N1")->xyz);
                double distN4O6 = XYZDist(&ResidueGetAtomByName(pRes, "N4")->xyz, &ResidueGetAtomByName(pResS, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
          else
          {
            // paired residue must be DG or G;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResS), "DG") == 0 || strcmp(ResidueGetName(pResS), "G") == 0)
              {
                // check the distance
                // atom O2, N3 and N4 from nucleotide DC
                // atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDist(&ResidueGetAtomByName(pRes, "O2")->xyz, &ResidueGetAtomByName(pResS, "N2")->xyz);
                double distN3N1 = XYZDist(&ResidueGetAtomByName(pRes, "N3")->xyz, &ResidueGetAtomByName(pResS, "N1")->xyz);
                double distN4O6 = XYZDist(&ResidueGetAtomByName(pRes, "N4")->xyz, &ResidueGetAtomByName(pResS, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
        }
      }
    }
    else if (!strcmp(ResidueGetName(pRes), "G"))
    {
      // paired residue is DC or C;
      for (int k = 0; k < StructureGetChainNum(pStruct); k++)
      {
        Chain *pChainK = StructureGetChain(pStruct, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            // paired residue must be C;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pRes) - ResidueGetPosInChain(pResS)) <= 1)
                continue;
              if (!strcmp(ResidueGetName(pResS), "C"))
              {
                // check the distance
                // atom O2, N3 and N4 from nucleotide C
                // atom N2, N1 and O6 from nucleotide G
                double distO2N2 = XYZDist(&ResidueGetAtomByName(pResS, "O2")->xyz, &ResidueGetAtomByName(pRes, "N2")->xyz);
                double distN3N1 = XYZDist(&ResidueGetAtomByName(pResS, "N3")->xyz, &ResidueGetAtomByName(pRes, "N1")->xyz);
                double distN4O6 = XYZDist(&ResidueGetAtomByName(pResS, "N4")->xyz, &ResidueGetAtomByName(pRes, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
          else
          {
            // paired residue must be DC or C;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResS), "DC") == 0 || strcmp(ResidueGetName(pResS), "C") == 0)
              {
                // check the distance
                // atom O2, N3 and N4 from nucleotide DC
                // atom N2, N1 and O6 from nucleotide DG
                double distO2N2 = XYZDist(&ResidueGetAtomByName(pResS, "O2")->xyz, &ResidueGetAtomByName(pRes, "N2")->xyz);
                double distN3N1 = XYZDist(&ResidueGetAtomByName(pResS, "N3")->xyz, &ResidueGetAtomByName(pRes, "N1")->xyz);
                double distN4O6 = XYZDist(&ResidueGetAtomByName(pResS, "N4")->xyz, &ResidueGetAtomByName(pRes, "O6")->xyz);
                if (distO2N2 >= 2.0 && distO2N2 <= 4.0 && distN3N1 >= 2.0 && distN3N1 <= 4.0 && distN4O6 >= 2.0 && distN4O6 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
        }
      }
    }
    else
    { //! strcmp(ResidueGetName(pQuery),"U")
      // paired residue is DA or A;
      for (int k = 0; k < StructureGetChainNum(pStruct); k++)
      {
        Chain *pChainK = StructureGetChain(pStruct, k);
        if (ChainGetType(pChainK) == Type_Chain_DNA || ChainGetType(pChainK) == Type_Chain_RNA)
        {
          if (!strcmp(ChainGetName(pChainK), ChainGetName(pChainI)))
          {
            // paired residue must be A;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (abs(ResidueGetPosInChain(pRes) - ResidueGetPosInChain(pResS)) <= 1)
                continue;
              if (!strcmp(ResidueGetName(pResS), "A"))
              {
                // check the distance
                // atom N1 and N6 from nucleotide A
                // atom N3 and O4 from nucleotide U
                double distN1N3 = XYZDist(&ResidueGetAtomByName(pResS, "N1")->xyz, &ResidueGetAtomByName(pRes, "N3")->xyz);
                double distN6O4 = XYZDist(&ResidueGetAtomByName(pResS, "N6")->xyz, &ResidueGetAtomByName(pRes, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResS;
                }
              }
            }
          }
          else
          {
            // paired residue must be DA or A;
            for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
            {
              Residue *pResS = ChainGetResidue(pChainK, s);
              if (strcmp(ResidueGetName(pResS), "DA") == 0 || strcmp(ResidueGetName(pResS), "A") == 0)
              {
                // check the distance
                // atom N1 and N6 from nucleotide DA or A
                // atom N3 and O4 from nucleotide U
                double distN1N3 = XYZDist(&ResidueGetAtomByName(pResS, "N1")->xyz, &ResidueGetAtomByName(pRes, "N3")->xyz);
                double distN6O4 = XYZDist(&ResidueGetAtomByName(pResS, "N6")->xyz, &ResidueGetAtomByName(pRes, "O4")->xyz);
                if (distN1N3 >= 2.0 && distN1N3 <= 4.0 && distN6O4 >= 2.0 && distN6O4 <= 4.0)
                {
                  return pResS;
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

int StructureMutateNucleotide(Residue *pQuery, char oneLetterName, Structure *pStructure, AtomParams *pAtomParaSet, ResTopos *pResiTopos)
{
  char newCleotideName[MAX_LEN_RES_NAME + 1];
  Chain *pChainI = StructureFindChainByName(pStructure, ResidueGetChainName(pQuery));
  if (ChainGetType(pChainI) == Type_Chain_DNA)
  {
    if (oneLetterName == 'a')
    {
      strcpy(newCleotideName, "DA");
    }
    else if (oneLetterName == 'c')
    {
      strcpy(newCleotideName, "DC");
    }
    else if (oneLetterName == 'g')
    {
      strcpy(newCleotideName, "DG");
    }
    else if (oneLetterName == 't')
    {
      strcpy(newCleotideName, "DT");
    }
    else if (oneLetterName == 'n')
    {
      strcpy(newCleotideName, "DN");
    }
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
    if (oneLetterName == 'a')
    {
      strcpy(newCleotideName, "A");
    }
    else if (oneLetterName == 'c')
    {
      strcpy(newCleotideName, "C");
    }
    else if (oneLetterName == 'g')
    {
      strcpy(newCleotideName, "G");
    }
    else if (oneLetterName == 'u')
    {
      strcpy(newCleotideName, "U");
    }
    else if (oneLetterName == 'n')
    {
      strcpy(newCleotideName, "N");
    }
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
  BondsDestroy(&newResi.bonds);
  // modify the new residue
  ResidueSetName(&newResi, newCleotideName);
  ResidueAddAtomsFromAtomParams(&newResi, pAtomParaSet);
  ResidueAddBondsFromResiTopos(&newResi, pResiTopos);
  // copy coordinates for main-chain atoms
  for (int i = 0; i < ResidueGetAtomNum(&newResi); i++)
  {
    Atom *pAtomNew = ResidueGetAtom(&newResi, i);
    if (pAtomNew->isBBAtom)
    {
      Atom *pAtomOld = ResidueGetAtomByName(pQuery, AtomGetName(pAtomNew));
      if (pAtomOld != NULL)
      {
        pAtomNew->xyz = pAtomOld->xyz;
        pAtomNew->isXyzValid = pAtomOld->isXyzValid;
      }
      else
      {
        pAtomNew->xyz.X = pAtomNew->xyz.Y = pAtomNew->xyz.Z = 0;
        pAtomNew->isXyzValid = FALSE;
      }
    }
  }
  // calculate coordinates for side-chain atoms
  if (strcmp(ResidueGetName(pQuery), "DA") == 0 || strcmp(ResidueGetName(pQuery), "DG") == 0 || strcmp(ResidueGetName(pQuery), "A") == 0 || strcmp(ResidueGetName(pQuery), "G") == 0)
  {
    if (strcmp(ResidueGetName(&newResi), "DA") == 0 || strcmp(ResidueGetName(&newResi), "DG") == 0 || strcmp(ResidueGetName(&newResi), "A") == 0 || strcmp(ResidueGetName(&newResi), "G") == 0)
    {
      Atom *pAtomNew = ResidueGetAtomByName(&newResi, "N9");
      Atom *pAtomOld = ResidueGetAtomByName(pQuery, "N9");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;
      pAtomNew = ResidueGetAtomByName(&newResi, "C4");
      pAtomOld = ResidueGetAtomByName(pQuery, "C4");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;
    }
    else if (strcmp(ResidueGetName(&newResi), "DC") == 0 || strcmp(ResidueGetName(&newResi), "DT") == 0 || strcmp(ResidueGetName(&newResi), "C") == 0 || strcmp(ResidueGetName(&newResi), "T") == 0 || strcmp(ResidueGetName(&newResi), "U") == 0)
    {
      Atom *pAtomNew = ResidueGetAtomByName(&newResi, "N1");
      Atom *pAtomOld = ResidueGetAtomByName(pQuery, "N9");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;

      pAtomNew = ResidueGetAtomByName(&newResi, "C2");
      pAtomOld = ResidueGetAtomByName(pQuery, "C4");
      ResTopology topNew;
      CharmmIC icNew;
      ResTopologyCreate(&topNew);
      CharmmICCreate(&icNew);
      ResToposGet(pResiTopos, ResidueGetName(&newResi), &topNew);
      ResTopologyFindCharmmIC(&topNew, AtomGetName(pAtomNew), &icNew);
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
      Atom *pAtomNew = ResidueGetAtomByName(&newResi, "N10");
      Atom *pAtomOld = ResidueGetAtomByName(pQuery, "N9");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;

      pAtomNew = ResidueGetAtomByName(&newResi, "C4");
      pAtomOld = ResidueGetAtomByName(pQuery, "C4");
      /*pAtomNew->xyz=pAtomOld->xyz;*/
      ResTopology topNew;
      CharmmIC icNew;
      ResTopologyCreate(&topNew);
      CharmmICCreate(&icNew);
      ResToposGet(pResiTopos, ResidueGetName(&newResi), &topNew);
      ResTopologyFindCharmmIC(&topNew, AtomGetName(pAtomNew), &icNew);
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
      Atom *pAtomNew = ResidueGetAtomByName(&newResi, "N9");
      Atom *pAtomOld = ResidueGetAtomByName(pQuery, "N1");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;

      pAtomNew = ResidueGetAtomByName(&newResi, "C4");
      pAtomOld = ResidueGetAtomByName(pQuery, "C2");
      ResTopology topNew;
      CharmmIC icNew;
      ResTopologyCreate(&topNew);
      CharmmICCreate(&icNew);
      ResToposGet(pResiTopos, ResidueGetName(&newResi), &topNew);
      ResTopologyFindCharmmIC(&topNew, AtomGetName(pAtomNew), &icNew);
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
      Atom *pAtomNew = ResidueGetAtomByName(&newResi, "N1");
      Atom *pAtomOld = ResidueGetAtomByName(pQuery, "N1");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;

      pAtomNew = ResidueGetAtomByName(&newResi, "C2");
      pAtomOld = ResidueGetAtomByName(pQuery, "C2");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;
    }
    else
    {
      Atom *pAtomNew = ResidueGetAtomByName(&newResi, "N10");
      Atom *pAtomOld = ResidueGetAtomByName(pQuery, "N1");
      pAtomNew->xyz = pAtomOld->xyz;
      pAtomNew->isXyzValid = pAtomOld->isXyzValid;

      pAtomNew = ResidueGetAtomByName(&newResi, "C4");
      pAtomOld = ResidueGetAtomByName(pQuery, "C2");
      ResTopology topNew;
      CharmmIC icNew;
      ResTopologyCreate(&topNew);
      CharmmICCreate(&icNew);
      ResToposGet(pResiTopos, ResidueGetName(&newResi), &topNew);
      ResTopologyFindCharmmIC(&topNew, AtomGetName(pAtomNew), &icNew);
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

int StructureCalcAminoAcidPropensityAndRamaEnergy(Structure *pStruct, AAppTable *pAAPP, RamaTable *pRama)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChain = StructureGetChain(pStruct, i);
    if (ChainGetType(pChain) == Type_Chain_Protein)
    {
      for (int j = 0; j < ChainGetResidueNum(pChain); j++)
      {
        Residue *pResi = ChainGetResidue(pChain, j);
        AminoAcidPropensityAndRamachandranEnergy(pResi, pAAPP, pRama);
      }
    }
  }
  return Success;
}

int StructureCalcAminoAcidDunbrackEnergy(Structure *pStructure, BBdepRotLib *pBBdep)
{
  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChain = StructureGetChain(pStructure, i);
    if (ChainGetType(pChain) == Type_Chain_Protein)
    {
      for (int j = 0; j < ChainGetResidueNum(pChain); j++)
      {
        Residue *pResi = ChainGetResidue(pChain, j);
        AminoAcidDunbrackEnergy(pResi, pBBdep);
      }
    }
  }
  return Success;
}
