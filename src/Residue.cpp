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

#include "Residue.h"
#include <string.h>
#include <ctype.h>

extern BOOL FLAG_READ_HYDROGEN;


int ResidueCreate(Residue* pThis)
{
  strcpy(pThis->name, "");
  strcpy(pThis->chainName, "");
  pThis->posInChain = -10000;
  pThis->desType = Type_DesType_Fixed;
  AtomArrayCreate(&pThis->atoms);
  StringArrayCreate(&pThis->patches);
  BondSetCreate(&pThis->bonds);
  pThis->terminalType = Type_ResIsNotTer;
  pThis->nCbIn10A = 0;
  pThis->internalEnergy = 0;
  pThis->backboneEnergy = 0;
  pThis->phipsi[0] = -60;
  pThis->phipsi[1] = 60;
  DoubleArrayCreate(&pThis->Xs, 0);
  pThis->dunbrack = 0;
  pThis->aapp = 0;
  pThis->rama = 0;
  pThis->isBBIntact = TRUE;
  pThis->isSCIntact = TRUE;
  strcpy(pThis->designAATypes, "");
  return Success;
}

int ResidueDestroy(Residue* pThis)
{
  BondSetDestroy(&pThis->bonds);
  AtomArrayDestroy(&pThis->atoms);
  StringArrayDestroy(&pThis->patches);
  DoubleArrayDestroy(&pThis->Xs);
  return Success;
}

int ResidueCopy(Residue* pThis, Residue* pOther)
{
  ResidueDestroy(pThis);
  strcpy(pThis->name, pOther->name);
  strcpy(pThis->chainName, pOther->chainName);
  pThis->posInChain = pOther->posInChain;
  pThis->desType = pOther->desType;
  pThis->nCbIn10A = pOther->nCbIn10A;
  AtomArrayCreate(&pThis->atoms);
  AtomArrayCopy(&pThis->atoms, &pOther->atoms);
  StringArrayCreate(&pThis->patches);
  StringArrayCopy(&pThis->patches, &pOther->patches);
  BondSetCopy(&pThis->bonds, &pOther->bonds);
  pThis->terminalType = pOther->terminalType;
  pThis->phipsi[0] = pOther->phipsi[0];
  pThis->phipsi[1] = pOther->phipsi[1];
  DoubleArrayCopy(&pThis->Xs, &pOther->Xs);
  pThis->dunbrack = pOther->dunbrack;
  pThis->aapp = pOther->aapp;
  pThis->rama = pOther->rama;
  pThis->isBBIntact = pOther->isBBIntact;
  pThis->isSCIntact = pOther->isSCIntact;
  return Success;
}

char* ResidueGetName(Residue* pThis)
{
  return pThis->name;
}

int ResidueSetName(Residue* pThis, char* newName)
{
  if (strlen(newName) > MAX_LEN_RES_NAME)
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

char* ResidueGetChainName(Residue* pThis)
{
  return pThis->chainName;
}

int ResidueSetChainName(Residue* pThis, char* newChainName)
{
  int i;
  if (strlen(newChainName) > MAX_LEN_CHAIN_NAME)
  {
    char errMsg[MAX_LEN_ERR_MSG + 1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s line %d, name is too long", __FILE__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->chainName, newChainName);
  for (i = 0;i < ResidueGetAtomCount(pThis);i++)
  {
    Atom* pCurAtom = ResidueGetAtom(pThis, i);
    AtomSetChainName(pCurAtom, newChainName);
  }
  return Success;
}

int ResidueGetPosInChain(Residue* pThis)
{
  return pThis->posInChain;
}

int ResidueSetPosInChain(Residue* pThis, int newPosInChain)
{
  int i;
  pThis->posInChain = newPosInChain;
  for (i = 0;i < ResidueGetAtomCount(pThis);i++)
  {
    Atom* pCurAtom = ResidueGetAtom(pThis, i);
    AtomSetPosInChain(pCurAtom, newPosInChain);
  }

  return Success;
}

int ResidueGetDesignType(Residue* pThis)
{
  return pThis->desType;
}


int ResidueSetDesignType(Residue* pThis, Type_ResidueDesignType newFlag)
{
  pThis->desType = newFlag;
  return Success;
}

double ResidueGetCharge(Residue* pThis)
{
  return AtomArrayCalcTotalCharge(&pThis->atoms);
}

int ResidueGetPolarity(Residue* pThis, Type_ResiduePolarity* pPolarity)
{
  // if the total charge is not zero, the type is Type_ResiduePolarity_Charged;
  // if the total charge is zero, but any non-root atom has non-zero charge, the type is Type_ResiduePolarity_Polar;
  // if the total charge is zero and all non-root atoms have zero charge, the type is Type_ResiduePolarity_NonPolar;
  BOOL nonRootAtomCharged = FALSE;
  for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
  {
    Atom* pAtom;
    BOOL nonRoot;
    BOOL charged;
    pAtom = ResidueGetAtom(pThis, i);
    nonRoot = !(pAtom->isBBAtom);
    charged = fabs(pAtom->charge) > MIN_ZERO_TOLERANCE;
    if (nonRoot && charged)
      nonRootAtomCharged = TRUE;
  }
  BOOL residueCharged = fabs(ResidueGetCharge(pThis)) > MIN_ZERO_TOLERANCE;
  if (residueCharged)
  {
    *pPolarity = Type_ResiduePolarity_Charged;
  }
  else if (nonRootAtomCharged)
  {
    *pPolarity = Type_ResiduePolarity_Polar;
  }
  else
  {
    *pPolarity = Type_ResiduePolarity_NonPolar;
  }
  return Success;
}

int ResidueGetAtomCount(Residue* pThis)
{
  return AtomArrayGetCount(&pThis->atoms);
}

Atom* ResidueGetAtom(Residue* pThis, int index)
{
  return AtomArrayGet(&pThis->atoms, index);
}

Atom* ResidueGetAtomByName(Residue* pThis, char* atomName)
{
  int index;
  int result = ResidueFindAtom(pThis, atomName, &index);
  if (FAILED(result))
  {
    return NULL;
  }
  else
  {
    return ResidueGetAtom(pThis, index);
  }
}

int ResidueFindAtom(Residue* pThis, char* atomName, int* pIndex)
{
  return AtomArrayFind(&pThis->atoms, atomName, pIndex);
}

int ResidueGetAtomXYZ(Residue* pThis, char* atomName, XYZ* pXYZ)
{
  Atom* pAtom = ResidueGetAtomByName(pThis, atomName);
  if (pAtom == NULL || pAtom->isXyzValid == FALSE)
  {
    return DataNotExistError;
  }
  *pXYZ = pAtom->xyz;
  return Success;
}

AtomArray* ResidueGetAllAtoms(Residue* pThis)
{
  return &pThis->atoms;
}

int ResidueInsertAtom(Residue* pThis, int newIndex, Atom* pNewAtom)
{
  int index;
  if (newIndex < 0 || newIndex > AtomArrayGetCount(&pThis->atoms))
  {
    char errMsg[MAX_LEN_ERR_MSG + 1];
    int errorCode = IndexError;
    sprintf(errMsg, "in file %s line %d, index is invalid", __FILE__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  else if (FAILED(ResidueFindAtom(pThis, AtomGetName(pNewAtom), &index)))
  {
    return AtomArrayInsert(&pThis->atoms, newIndex, pNewAtom);
  }
  else
  {
    return AtomCopy(ResidueGetAtom(pThis, index), pNewAtom);
  }
}

int ResidueAddAtom(Residue* pThis, Atom* pNewAtom)
{
  int index;
  if (FAILED(ResidueFindAtom(pThis, AtomGetName(pNewAtom), &index)))
  {
    return AtomArrayAppend(&pThis->atoms, pNewAtom);
  }
  else
  {
    return AtomCopy(ResidueGetAtom(pThis, index), pNewAtom);
  }
}

int ResidueAddAtomsFromAtomParams(Residue* pThis, AtomParamsSet* pAtomParams)
{
  int count;
  int result = AtomParamsSetGetAtomCount(pAtomParams, ResidueGetName(pThis), &count);
  if (FAILED(result))
  {
    return result;
  }

  for (int i = 0;i < count;i++)
  {
    BOOL atomExist = FALSE;
    Atom newAtom;
    AtomCreate(&newAtom);
    AtomParamsSetGetAtomParam(pAtomParams, ResidueGetName(pThis), i, &newAtom);
    for (int j = 0; j < ResidueGetAtomCount(pThis); j++)
    {
      Atom* pAtom1 = ResidueGetAtom(pThis, j);
      if (strcmp(AtomGetName(pAtom1), AtomGetName(&newAtom)) == 0)
      {
        Atom tempAtom;
        AtomCreate(&tempAtom);
        AtomCopy(&tempAtom, pAtom1);
        AtomCopy(pAtom1, &newAtom);
        AtomSetChainName(pAtom1, AtomGetChainName(&tempAtom));
        AtomSetPosInChain(pAtom1, AtomGetPosInChain(&tempAtom));
        pAtom1->xyz = tempAtom.xyz;
        pAtom1->isXyzValid = tempAtom.isXyzValid;
        pAtom1->bfactor = tempAtom.bfactor;
        AtomDestroy(&tempAtom);
        atomExist = TRUE;
        break;
      }
    }
    if (atomExist == FALSE)
    {
      AtomSetChainName(&newAtom, ResidueGetChainName(pThis));
      AtomSetPosInChain(&newAtom, ResidueGetPosInChain(pThis));
      ResidueAddAtom(pThis, &newAtom);
    }
    AtomDestroy(&newAtom);
  }
  return Success;
}

int ResidueDeleteAtom(Residue* pThis, char* atomName)
{
  int index;
  if (FAILED(ResidueFindAtom(pThis, atomName, &index)))
  {
    return DataNotExistError;
  }
  else
  {
    return AtomArrayRemove(&pThis->atoms, index);
  }
}


int PDBReaderCheckHisState(FileReader* pPDBFileReader, int* state)
{
  int curPosInReader = FileReaderGetCurrentPos(pPDBFileReader);
  // the following three lines show the format of a pdb file
  // type, serial, name, altLoc, resName, chainID, resPos, ..., X,   Y,   Z
  // 0,    6,      12,   16,     17,      21,      22,     ..., 30,  38,  46
  // 6,    5,       4,    1,      4,       1,       5,   ...,    8,   8,  8
  BOOL isHD1 = FALSE;
  BOOL isHE2 = FALSE;
  BOOL firstLine = TRUE;
  char iniSeqPos[MAX_LEN_ONE_LINE_CONTENT + 1] = "UNKNOWN";
  char line[MAX_LEN_ONE_LINE_CONTENT + 1];
  while (!FAILED(FileReaderGetNextLine(pPDBFileReader, line)))
  {
    char strType[MAX_LEN_ONE_LINE_CONTENT + 1];
    char strAtomName[MAX_LEN_ONE_LINE_CONTENT + 1];
    char strResName[MAX_LEN_ONE_LINE_CONTENT + 1];
    char strResPos[MAX_LEN_ONE_LINE_CONTENT + 1];
    ExtractTargetStringFromSourceString(strType, line, 0, 4);

    if (strcmp(strType, "ENDM") == 0) break;
    if (strcmp(strType, "ATOM") != 0 && strcmp(strType, "HETA") != 0) continue;
    ExtractTargetStringFromSourceString(strAtomName, line, 12, 4);
    ExtractTargetStringFromSourceString(strResName, line, 17, 4);
    ExtractTargetStringFromSourceString(strResPos, line, 22, 5);

    if (firstLine)
    {
      strcpy(iniSeqPos, strResPos);
      firstLine = FALSE;
    }
    else
    {
      // a new residue encountered, break the loop
      if (strcmp(iniSeqPos, strResPos) != 0) break;
    }

    // Only a histidine residue is allowed for checking His status
    if (strcmp(strResName, "HIS") != 0 && strcmp(strResName, "HSD") != 0 && strcmp(strResName, "HSE") != 0 && strcmp(strResName, "HSP") != 0)
    {
      char errMsg[MAX_LEN_ERR_MSG + 1];
      sprintf(errMsg, "in file %s line %d, a histidine residue (with name HIS, HSD, HSE, or HSP) is expected", __FILE__, __LINE__);
      TraceError(errMsg, FormatError);
    }

    // read the hydrogen atoms
    if (strAtomName[0] == 'H')
    {
      //continue;
    }
    else if (isdigit(strAtomName[0]) && strAtomName[1] == 'H' && (int)strlen(strAtomName) == 4)
    {
      char tempName[MAX_LEN_ATOM_NAME + 1];
      tempName[0] = strAtomName[1];
      tempName[1] = strAtomName[2];
      tempName[2] = strAtomName[3];
      tempName[3] = strAtomName[0];
      tempName[4] = '\0';
      strcpy(strAtomName, tempName);
      //continue;
    }

    if (!strcmp(strAtomName, "HD1")) isHD1 = TRUE;
    else if (!strcmp(strAtomName, "HE2")) isHE2 = TRUE;

  }

  //remember to reset the pointer of FileReader
  FileReaderSetCurrentPos(pPDBFileReader, curPosInReader + 1);

  if (isHD1 && isHE2 == FALSE) *state = 1;      // HIS -> HSD
  else if (isHD1 == FALSE && isHE2) *state = 2; // HIS -> HSE
  else if (isHD1 && isHE2) *state = 3;  // HIS -> HSP
  else *state = 1;                                 // HIS -> HSD
  return Success;
}


int ResidueReadXYZFromPDB(Residue* pThis, FileReader* pPDBFileReader)
{
  // the following three lines show the format of a pdb file
  // type, serial, name, altLoc, resName, chainID, resPos, ..., X,   Y,   Z
  // 0,    6,      12,   16,     17,      21,      22,     ..., 30,  38,  46
  // 6,    5,       4,    1,      4,       1,       5,   ...,    8,   8,  8
  BOOL firstLine = TRUE;
  char iniSeqPos[MAX_LEN_ONE_LINE_CONTENT + 1] = "UNKNOWN";
  char line[MAX_LEN_ONE_LINE_CONTENT + 1];
  while (!FAILED(FileReaderGetNextLine(pPDBFileReader, line)))
  {
    char strType[MAX_LEN_ONE_LINE_CONTENT + 1];
    char strAtomName[MAX_LEN_ONE_LINE_CONTENT + 1];
    char strResName[MAX_LEN_ONE_LINE_CONTENT + 1];
    char strResPos[MAX_LEN_ONE_LINE_CONTENT + 1];
    char strX[MAX_LEN_ONE_LINE_CONTENT + 1];
    char strY[MAX_LEN_ONE_LINE_CONTENT + 1];
    char strZ[MAX_LEN_ONE_LINE_CONTENT + 1];
    char bFactor[MAX_LEN_ONE_LINE_CONTENT + 1] = "";
    ExtractTargetStringFromSourceString(strType, line, 0, 4);

    //for smallmol residue, read its vdw energies if possible
    if (strcmp(strType, "ENER") == 0)
    {
      StringArray strings;
      StringArrayCreate(&strings);
      StringArraySplitString(&strings, line, ' ');
      pThis->internalEnergy = atof(StringArrayGet(&strings, 2));
      pThis->backboneEnergy = atof(StringArrayGet(&strings, 4));
      StringArrayDestroy(&strings);
    }

    if (strcmp(strType, "ENDM") == 0) return Success;
    if (strcmp(strType, "ATOM") && strcmp(strType, "HETA")) continue;

    ExtractTargetStringFromSourceString(strAtomName, line, 12, 4);
    ExtractTargetStringFromSourceString(strResName, line, 17, 4);
    ExtractTargetStringFromSourceString(strResPos, line, 22, 5);
    ExtractTargetStringFromSourceString(strX, line, 30, 8);
    ExtractTargetStringFromSourceString(strY, line, 38, 8);
    ExtractTargetStringFromSourceString(strZ, line, 46, 8);
    ExtractTargetStringFromSourceString(bFactor, line, 60, 6);
    // if the residue is HIS, use HSD or HSE
    if (strcmp(strResName, "HIS") == 0 && strcmp(pThis->name, "HSD") == 0) strcpy(strResName, "HSD");
    else if (strcmp(strResName, "HIS") == 0 && strcmp(pThis->name, "HSE") == 0) strcpy(strResName, "HSE");
    // the CD1 atom of ILE in pdb is altered into CD
    if (strcmp(strResName, "ILE") == 0 && strcmp(strAtomName, "CD1") == 0) strcpy(strAtomName, "CD");

    // use flag to determine if the hydrogen atoms' coordinates are read
    if (isdigit(strAtomName[0]) && isalpha(strAtomName[1]))
    {
      char tempName[MAX_LEN_ATOM_NAME + 1];
      strcpy(tempName, strAtomName + 1);
      tempName[strlen(strAtomName) - 1] = strAtomName[0];
      tempName[strlen(strAtomName)] = '\0';
      strcpy(strAtomName, tempName);
    }
    if (FLAG_READ_HYDROGEN == FALSE && strAtomName[0] == 'H')
    {
      continue;
    }

    if (firstLine)
    {
      strcpy(iniSeqPos, strResPos);
      firstLine = FALSE;
    }
    else
    {
      // read the line of a new residue
      if (strcmp(iniSeqPos, strResPos) != 0)
      {
        FileReaderSetCurrentPos(pPDBFileReader, FileReaderGetCurrentPos(pPDBFileReader) - 1);
        return Success;
      }
    }
    // read OXT coordinate from PDB file instead of recalculation
    if (strcmp(strAtomName, "OXT") == 0)
    {
      Atom atomOXT;
      AtomCreate(&atomOXT);
      AtomCopy(&atomOXT, ResidueGetAtomByName(pThis, "O"));
      atomOXT.xyz.X = atof(strX); atomOXT.xyz.Y = atof(strY); atomOXT.xyz.Z = atof(strZ);
      atomOXT.isXyzValid = TRUE;
      atomOXT.charge = -0.55;
      strcpy(atomOXT.name, "OXT");
      ResidueGetAtomByName(pThis, "O")->charge = -0.55;
      ResidueGetAtomByName(pThis, "C")->charge = 0.1;
      if (strcmp(bFactor, "") != 0 && strcmp(bFactor, "0.00") != 0)
      {
        atomOXT.bfactor = atof(bFactor);
      }
      AtomArrayAppend(&pThis->atoms, &atomOXT);
      AtomDestroy(&atomOXT);
    }
    else
    {
      Atom* pAtom = ResidueGetAtomByName(pThis, strAtomName);
      if (pAtom == NULL || pAtom->isXyzValid) continue;
      pAtom->xyz.X = atof(strX); pAtom->xyz.Y = atof(strY); pAtom->xyz.Z = atof(strZ);
      pAtom->isXyzValid = TRUE;
      if (strcmp(bFactor, "") != 0 && strcmp(bFactor, "0.00") != 0)
      {
        pAtom->bfactor = atof(bFactor);
      }
    }
  }
  return Success;
}


BondSet* ResidueGetBonds(Residue* pThis)
{
  return &pThis->bonds;
}


int ResidueAddBondsFromResiTopos(Residue* pThis, ResiTopoSet* pResiTopoCollection)
{
  ResidueTopology resiTopo;
  ResidueTopologyCreate(&resiTopo);
  if (FAILED(ResiTopoSetGet(pResiTopoCollection, pThis->name, &resiTopo)))
  {
    return DataNotExistError;
  }
  BondSetCopy(&pThis->bonds, ResidueTopologyGetBonds(&resiTopo));
  ResidueTopologyDestroy(&resiTopo);
  return Success;
}


int ResidueShowInPDBFormat(Residue* pThis, char* header, char* chainName, int atomIndex, int resiIndex, FILE* pFile)
{
  char resiName[MAX_LEN_RES_NAME + 1];
  strcpy(resiName, ResidueGetName(pThis));
  if (strcmp(resiName, "HSD") == 0 || strcmp(resiName, "HSE") == 0 || strcmp(resiName, "HSP") == 0) strcpy(resiName, "HIS");
  AtomArrayShowInPDBFormat(&pThis->atoms, header, resiName, chainName, atomIndex, resiIndex, pFile);
  return Success;
}


int ResiduePatch(Residue* pThis, char* patchName, AtomParamsSet* pAtomParam, ResiTopoSet* pTopos)
{
  // delete old atoms
  ResidueTopology patchResiTopo;
  ResidueTopologyCreate(&patchResiTopo);
  int result = ResiTopoSetGet(pTopos, patchName, &patchResiTopo);
  if (FAILED(result)) return result;
  StringArray* deleteAtoms = ResidueTopologyGetDeletes(&patchResiTopo);
  for (int i = 0;i < StringArrayGetCount(deleteAtoms);i++)
  {
    ResidueDeleteAtom(pThis, StringArrayGet(deleteAtoms, i));
  }

  // Temporarily reset its name to that of the patch residue, to find topology in the topology collection
  char resiName[MAX_LEN_RES_NAME];
  strcpy(resiName, ResidueGetName(pThis));
  ResidueSetName(pThis, patchName);
  result = ResidueAddAtomsFromAtomParams(pThis, pAtomParam);
  if (FAILED(result)) return result;
  ResidueSetName(pThis, resiName);

  // new patches are stored at the head
  StringArrayInsert(&pThis->patches, 0, patchName);

  // deal with the bonds
  BondSet* pPatchBonds = ResidueTopologyGetBonds(&patchResiTopo);
  for (int i = 0;i < BondSetGetCount(pPatchBonds);i++)
  {
    Bond* pCurBond = BondSetGet(pPatchBonds, i);
    BondSetAdd(&pThis->bonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSet newBonds;
  BondSetCreate(&newBonds);
  for (int i = 0;i < BondSetGetCount(&pThis->bonds);i++)
  {
    Bond* pCurBond = BondSetGet(&pThis->bonds, i);
    char* atomName = BondGetFromName(pCurBond);
    if (atomName[0] != '+' && atomName[0] != '-' && ResidueGetAtomByName(pThis, atomName) == NULL)
    {
      continue;
    }
    atomName = BondGetToName(pCurBond);
    if (atomName[0] != '+' && atomName[0] != '-' && ResidueGetAtomByName(pThis, atomName) == NULL)
    {
      continue;
    }
    BondSetAdd(&newBonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSetCopy(&pThis->bonds, &newBonds);
  BondSetDestroy(&newBonds);
  ResidueTopologyDestroy(&patchResiTopo);
  return Success;
}


// there is a bug in previous ResiduePatch for CTER patching, because it will recalculate the coordinate of atom O and OXT
// but most time the coordinates of O and OXT will be determined, so it should be directly read from file instead of recalculating
int ResiduePatchCTER(Residue* pThis, char* patchName, AtomParamsSet* pAtomParam, ResiTopoSet* pTopos)
{
  // delete old atoms
  ResidueTopology patchResiTopo;
  ResidueTopologyCreate(&patchResiTopo);
  int result = ResiTopoSetGet(pTopos, "CTER", &patchResiTopo);
  if (FAILED(result)) return result;
  // do not delete atom O;
  //StringArray* deleteAtoms = ResidueTopologyGetDeletes(&patchResiTopo);
  //for(int i=0;i<StringArrayGetCount(deleteAtoms);i++){
  //  ResidueDeleteAtom(pThis, StringArrayGet(deleteAtoms, i));
  //}

  // Temporarily reset its name to that of the patch residue, to find topology in the topology collection
  char resiName[MAX_LEN_RES_NAME];
  strcpy(resiName, ResidueGetName(pThis));
  ResidueSetName(pThis, patchName);
  result = ResidueAddAtomsFromAtomParams(pThis, pAtomParam);
  // after patching, please check the chain name of the atom
  //for(int i=1; i<ResidueGetAtomCount(pThis); i++){
  //  Atom* pAtom = ResidueGetAtom(pThis,i);
  //  if(strcmp(AtomGetChainName(pAtom),"")==0){
  //    AtomSetChainName(pAtom,AtomGetChainName(ResidueGetAtom(pThis,0)));
  //  }
  //}
  if (FAILED(result)) return result;
  ResidueSetName(pThis, resiName);

  // new patches are stored at the head
  StringArrayInsert(&pThis->patches, 0, patchName);
  // deal with the bonds
  BondSet* pPatchBonds = ResidueTopologyGetBonds(&patchResiTopo);
  for (int i = 0;i < BondSetGetCount(pPatchBonds);i++)
  {
    Bond* pCurBond = BondSetGet(pPatchBonds, i);
    BondSetAdd(&pThis->bonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSet newBonds;
  BondSetCreate(&newBonds);
  for (int i = 0;i < BondSetGetCount(&pThis->bonds);i++)
  {
    Bond* pCurBond = BondSetGet(&pThis->bonds, i);
    char* atomName = BondGetFromName(pCurBond);
    if (atomName[0] != '+' && atomName[0] != '-' && ResidueGetAtomByName(pThis, atomName) == NULL)
    {
      continue;
    }
    atomName = BondGetToName(pCurBond);
    if (atomName[0] != '+' && atomName[0] != '-' && ResidueGetAtomByName(pThis, atomName) == NULL)
    {
      continue;
    }
    BondSetAdd(&newBonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSetCopy(&pThis->bonds, &newBonds);
  BondSetDestroy(&newBonds);
  ResidueTopologyDestroy(&patchResiTopo);

  // remove the useless bonds
  char deleteBondWithPreviousOrNextResidue = '+';
  for (int i = 0;i < BondSetGetCount(&pThis->bonds);i++)
  {
    Bond* pCurBond = BondSetGet(&pThis->bonds, i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    if (fromAtomName[0] == deleteBondWithPreviousOrNextResidue || toAtomName[0] == deleteBondWithPreviousOrNextResidue)
    {
      BondSetRemove(&pThis->bonds, fromAtomName, toAtomName);
      break;
    }
  }
  return Success;
}

int ResiduePatchNTERorCTER(Residue* pThis, char* NTERorCTER, AtomParamsSet* pAtomParam, ResiTopoSet* pTopos)
{
  int result;
  char deleteBondWithPreviousOrNextResidue = ' ';
  if (strcmp(NTERorCTER, "NTER") == 0)
  {
    if (strcmp(ResidueGetName(pThis), "GLY") == 0)
    {
      result = ResiduePatch(pThis, "GLYP", pAtomParam, pTopos);
    }
    else if (strcmp(ResidueGetName(pThis), "PRO") == 0)
    {
      result = ResiduePatch(pThis, "PROP", pAtomParam, pTopos);
    }
    else
    {
      result = ResiduePatch(pThis, "NTER", pAtomParam, pTopos);
    }
    deleteBondWithPreviousOrNextResidue = '-';
  }
  else if (strcmp(NTERorCTER, "CTER") == 0)
  {
    result = ResiduePatch(pThis, "CTER", pAtomParam, pTopos);
    deleteBondWithPreviousOrNextResidue = '+';
  }
  else
  {
    result = ValueError;
  }

  if (FAILED(result))
  {
    return result;
  }
  else
  {
    Bond* pCurBond;
    char* fromAtomName;
    char* toAtomName;
    int i;
    for (i = 0;i < BondSetGetCount(&pThis->bonds);i++)
    {
      pCurBond = BondSetGet(&pThis->bonds, i);
      fromAtomName = BondGetFromName(pCurBond);
      toAtomName = BondGetToName(pCurBond);
      if (fromAtomName[0] == deleteBondWithPreviousOrNextResidue || toAtomName[0] == deleteBondWithPreviousOrNextResidue)
      {
        BondSetRemove(&pThis->bonds, fromAtomName, toAtomName);
        break;
      }
    }
  }
  return Success;
}

StringArray* ResidueGetPatchingHistory(Residue* pThis)
{
  return &pThis->patches;
}


int ResidueCalcAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos, Residue* pPrevResi, Residue* pNextResi, char* atomName, XYZ* pDestXYZ)
{
  // find IC from patch topology orderly;
  BOOL icFound = FALSE;
  CharmmIC ic;
  ResidueTopology topo;
  CharmmICCreate(&ic);
  ResidueTopologyCreate(&topo);
  for (int i = 0;i < StringArrayGetCount(&pThis->patches);i++)
  {
    ResiTopoSetGet(pResiTopos, StringArrayGet(&pThis->patches, i), &topo);
    if (!FAILED(ResidueTopologyFindCharmmIC(&topo, atomName, &ic)))
    {
      icFound = TRUE;
      break;
    }
  }
  // find IC in the residue topology;
  if (icFound == FALSE)
  {
    ResiTopoSetGet(pResiTopos, ResidueGetName(pThis), &topo);
    if (FAILED(ResidueTopologyFindCharmmIC(&topo, atomName, &ic)))
    {
      int result = DataNotExistError;
      return result;
    }
  }

  char* namesOfAtomABC[3];
  XYZ xyzsOfAtomABC[3];
  namesOfAtomABC[0] = CharmmICGetAtomA(&ic);
  namesOfAtomABC[1] = CharmmICGetAtomB(&ic);
  namesOfAtomABC[2] = CharmmICGetAtomC(&ic);
  // find atom coordinate from the current residue, preceding residue and next residue;
  for (int i = 0;i < 3;i++)
  {
    if (namesOfAtomABC[i][0] == '-' && pPrevResi != NULL && !FAILED(ResidueGetAtomXYZ(pPrevResi, namesOfAtomABC[i] + 1, &xyzsOfAtomABC[i])))
    {
      continue;
    }
    else if (namesOfAtomABC[i][0] == '+' && pNextResi != NULL && !FAILED(ResidueGetAtomXYZ(pNextResi, namesOfAtomABC[i] + 1, &xyzsOfAtomABC[i])))
    {
      continue;
    }
    else if (namesOfAtomABC[i][0] != '+' && namesOfAtomABC[i][0] != '-' && !FAILED(ResidueGetAtomXYZ(pThis, namesOfAtomABC[i], &xyzsOfAtomABC[i])))
    {
      continue;
    }
    else
    {
      ResidueTopologyDestroy(&topo);
      CharmmICDestroy(&ic);
      return DataNotExistError;
    }
  }

  GetFourthAtom(&xyzsOfAtomABC[0], &xyzsOfAtomABC[1], &xyzsOfAtomABC[2], CharmmICGetICParams(&ic), pDestXYZ);
  ResidueTopologyDestroy(&topo);
  CharmmICDestroy(&ic);
  return Success;
}


int ResidueCalcAllAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos, Residue* pPrevResi, Residue* pNextResi)
{
  BOOL allAtomsXYZAreValid = FALSE;
  BOOL done = FALSE;
  while (!done)
  {
    done = TRUE;
    allAtomsXYZAreValid = TRUE;
    for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
    {
      XYZ newXYZ;
      Atom* pCurAtom = ResidueGetAtom(pThis, i);
      if (pCurAtom->isXyzValid) continue;
      int result = ResidueCalcAtomXYZ(pThis, pResiTopos, pPrevResi, pNextResi, AtomGetName(pCurAtom), &newXYZ);
      if (FAILED(result))
      {
        allAtomsXYZAreValid = FALSE;
        continue;
      }
      else
      {
        pCurAtom->xyz = newXYZ;
        pCurAtom->isXyzValid = TRUE;
        done = FALSE; // New atom XYZ has been calculated in this 'while' loop, go on and try to find more
      }
    }
  }

  if (!allAtomsXYZAreValid) return DataNotExistError;

  return Success;
}


int ResidueCalcAllBackboneAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos, Residue* pPrevResi, Residue* pNextResi)
{
  BOOL allBackboneAtomsXYZAreValid = FALSE;
  BOOL done = FALSE;
  while (!done)
  {
    done = TRUE;
    allBackboneAtomsXYZAreValid = TRUE;
    for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
    {
      XYZ newXYZ;
      Atom* pCurAtom = ResidueGetAtom(pThis, i);
      if (!pCurAtom->isBBAtom) continue;
      if (pCurAtom->isXyzValid) continue;
      int result = ResidueCalcAtomXYZ(pThis, pResiTopos, pPrevResi, pNextResi, AtomGetName(pCurAtom), &newXYZ);
      if (FAILED(result))
      {
        allBackboneAtomsXYZAreValid = FALSE;
        continue;
      }
      else
      {
        pCurAtom->xyz = newXYZ;
        pCurAtom->isXyzValid = TRUE;
        done = FALSE; // New atom XYZ has been calculated in this 'while' loop, go on and try to find more
      }
    }
  }

  if (!allBackboneAtomsXYZAreValid) return DataNotExistError;

  return Success;
}


int ResidueCalcAllSidechainAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos)
{
  BOOL allSidechainAtomsXYZAreValid = FALSE;
  BOOL done = FALSE;
  while (!done)
  {
    done = TRUE;
    allSidechainAtomsXYZAreValid = TRUE;
    for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
    {
      XYZ newXYZ;
      Atom* pCurAtom = ResidueGetAtom(pThis, i);
      if (pCurAtom->isBBAtom) continue;
      if (pCurAtom->isXyzValid) continue;
      int result = ResidueCalcAtomXYZ(pThis, pResiTopos, NULL, NULL, AtomGetName(pCurAtom), &newXYZ);
      if (FAILED(result))
      {
        allSidechainAtomsXYZAreValid = FALSE;
        continue;
      }
      else
      {
        pCurAtom->xyz = newXYZ;
        pCurAtom->isXyzValid = TRUE;
        done = FALSE; // New atom XYZ has been calculated in this 'while' loop, go on and try to find more
      }
    }
  }

  if (!allSidechainAtomsXYZAreValid) return DataNotExistError;

  return Success;
}

int ResidueShowAtomParameter(Residue* pThis)
{
  for (int i = 0; i < ResidueGetAtomCount(pThis); ++i)
  {
    Atom* pAtom = ResidueGetAtom(pThis, i);
    AtomShowAtomParameter(pAtom);
  }
  return Success;
}

int ResidueShowBondInformation(Residue* pThis)
{
  BondSetShow(ResidueGetBonds(pThis));
  return Success;
}

double ResidueAndResidueSidechainRMSD(Residue* pThis, Residue* pOther)
{
  double rmsd = 0.0;
  int count = 0;
  if (strcmp(ResidueGetName(pThis), ResidueGetName(pOther)) == 0 ||
    (strcmp(ResidueGetName(pThis), "HSD") == 0 && strcmp(ResidueGetName(pOther), "HSE") == 0) ||
    (strcmp(ResidueGetName(pThis), "HSE") == 0 && strcmp(ResidueGetName(pOther), "HSD") == 0))
  {
    for (int i = 0; i < ResidueGetAtomCount(pThis); ++i)
    {
      Atom* pAtom1 = ResidueGetAtom(pThis, i);
      if (AtomIsHydrogen(pAtom1) || pAtom1->isBBAtom || strcmp(pAtom1->name, "CB") == 0) continue;
      Atom* pAtom2 = ResidueGetAtomByName(pOther, AtomGetName(pAtom1));
      rmsd += (pAtom1->xyz.X - pAtom2->xyz.X) * (pAtom1->xyz.X - pAtom2->xyz.X) + (pAtom1->xyz.Y - pAtom2->xyz.Y) * (pAtom1->xyz.Y - pAtom2->xyz.Y) + (pAtom1->xyz.Z - pAtom2->xyz.Z) * (pAtom1->xyz.Z - pAtom2->xyz.Z);
      count++;
    }
    if (count > 0) return sqrt(rmsd / count);
    else return 0.0;
  }
  //when residues are different, return a big value
  else
  {
    return 1e8;
  }
}

BOOL LigandResidueNameConflictWithAminoAcid(char* ligname)
{
  BOOL result = FALSE;
  char resiname[][MAX_LEN_RES_NAME + 1] = {
    "ALA","CYS","ASP","GLU","PHE","GLY","HIS","HSE","HSD","ILE","LYS","LEU",
    "MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR" };
  for (int i = 0; i < 22; i++)
  {
    if (strcmp(ligname, resiname[i]) == 0)
    {
      result = TRUE;
      break;
    }
  }
  return result;
}


int ResidueCheckAtomCoordinateValidity(Residue* pThis)
{
  char errMsg[MAX_LEN_ERR_MSG + 1];
  for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
  {
    Atom* pAtom = ResidueGetAtom(pThis, i);
    if (AtomIsHydrogen(pAtom)) continue;
    if (pAtom->isBBAtom && pAtom->isXyzValid == FALSE)
    {
      if (strcmp(AtomGetName(pAtom), "N") == 0 || strcmp(AtomGetName(pAtom), "CA") == 0 || strcmp(AtomGetName(pAtom), "C") == 0)
      {
        sprintf(errMsg, "in file %s line %d, the coordinate of backbone atom %s on residue %s%d%s is missing, please check", __FILE__, __LINE__, AtomGetName(pAtom), ResidueGetChainName(pThis), ResidueGetPosInChain(pThis), ResidueGetName(pThis));
        TraceError(errMsg, ValueError);
        pThis->isBBIntact = FALSE;
      }
    }
    else if (pAtom->isBBAtom == FALSE && pAtom->isXyzValid == FALSE)
    {
      pThis->isSCIntact = FALSE;
    }
  }
  return Success;
}

BOOL ResidueIsSymmetricalCheck(Residue* pResidue)
{
  if (strcmp(ResidueGetName(pResidue), "PHE") == 0 ||
    strcmp(ResidueGetName(pResidue), "TYR") == 0 ||
    strcmp(ResidueGetName(pResidue), "ASP") == 0 ||
    strcmp(ResidueGetName(pResidue), "GLU") == 0 ||
    strcmp(ResidueGetName(pResidue), "ARG") == 0)
  {
    return TRUE;
  }
  return FALSE;
}


int SymmetricalResidueGenerate(Residue* pSymmetrical, Residue* pResidue)
{
  ResidueCopy(pSymmetrical, pResidue);
  if (strcmp(ResidueGetName(pResidue), "PHE") == 0 || strcmp(ResidueGetName(pResidue), "TYR") == 0)
  {
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "CD1"), "TMP");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "CD2"), "CD1");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "TMP"), "CD2");

    AtomSetName(ResidueGetAtomByName(pSymmetrical, "CE1"), "TMP");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "CE2"), "CE1");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "TMP"), "CE2");
  }
  else if (strcmp(ResidueGetName(pResidue), "ASP") == 0)
  {
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "OD1"), "TMP");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "OD2"), "OD1");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "TMP"), "OD2");
  }
  else if (strcmp(ResidueGetName(pResidue), "GLU") == 0)
  {
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "OE1"), "TMP");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "OE2"), "OE1");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "TMP"), "OE2");
  }
  else if (strcmp(ResidueGetName(pResidue), "ARG") == 0)
  {
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "NH1"), "TMP");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "NH2"), "NH1");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "TMP"), "NH2");
  }
  else if (strcmp(ResidueGetName(pResidue), "HSD") == 0 || strcmp(ResidueGetName(pResidue), "HSE") == 0)
  {
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "ND1"), "TMP");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "CD2"), "ND1");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "TMP"), "CD2");

    AtomSetName(ResidueGetAtomByName(pSymmetrical, "CE1"), "TMP");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "ND2"), "CE1");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "TMP"), "NE2");
  }
  else if (strcmp(ResidueGetName(pResidue), "ASN") == 0)
  {
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "OD1"), "TMP");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "ND2"), "OD1");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "TMP"), "ND2");
  }
  else if (strcmp(ResidueGetName(pResidue), "GLN") == 0)
  {
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "OE1"), "TMP");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "NE2"), "OE1");
    AtomSetName(ResidueGetAtomByName(pSymmetrical, "TMP"), "NE2");
  }

  return Success;
}

BOOL ResidueAndResidueAllTorsionsAreSimilar(Residue* pThis, Residue* pOther, double cutoff)
{
  if (ResidueAndResidueInSameType(pThis, pOther))
  {
    BOOL match = TRUE;
    for (int j = 0;j < DoubleArrayGetLength(&pOther->Xs);j++)
    {
      double min = DoubleArrayGet(&pOther->Xs, j) - DegToRad(cutoff);
      double max = DoubleArrayGet(&pOther->Xs, j) + DegToRad(cutoff);
      double torsion = DoubleArrayGet(&pThis->Xs, j);
      double torsionm2pi = torsion - 2 * PI;
      double torsionp2pi = torsion + 2 * PI;
      double torsion2 = torsion;
      if ((strcmp(ResidueGetName(pThis), "PHE") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "TYR") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "ASP") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "GLU") == 0 && j == 2) ||
        (strcmp(ResidueGetName(pThis), "GLN") == 0 && j == 2) ||
        (strcmp(ResidueGetName(pThis), "ASN") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "HSD") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "HSE") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "ARG") == 0 && j == 3))
      {
        torsion2 = torsion + PI;
        torsion2 = torsion > 0 ? torsion - PI : torsion2;
      }
      double torsion2m2pi = torsion2 - 2 * PI;
      double torsion2p2pi = torsion2 + 2 * PI;
      if (!(
        (torsion <= max && torsion >= min) ||
        (torsionm2pi <= max && torsionm2pi >= min) ||
        (torsionp2pi <= max && torsionp2pi >= min) ||
        (torsion2 <= max && torsion2 >= min) ||
        (torsion2m2pi <= max && torsion2m2pi >= min) ||
        (torsion2p2pi <= max && torsion2p2pi >= min)
        ))
      {
        match = FALSE;
        break;
      }
    }
    return match;
  }
  return FALSE;
}

BOOL ResidueAndResidueCheckTorsionSimilarity(Residue* pThis, Residue* pOther, double cutoff, IntArray* pSimArray)
{
  BOOL allmatch = TRUE;
  if (ResidueAndResidueInSameType(pThis, pOther))
  {
    for (int j = 0;j < DoubleArrayGetLength(&pOther->Xs);j++)
    {
      //if(j>1) break;
      //int match=1;
      double min = DoubleArrayGet(&pOther->Xs, j) - DegToRad(cutoff);
      double max = DoubleArrayGet(&pOther->Xs, j) + DegToRad(cutoff);
      double torsion = DoubleArrayGet(&pThis->Xs, j);
      double torsionm2pi = torsion - 2 * PI;
      double torsionp2pi = torsion + 2 * PI;
      double torsion2 = torsion;
      if ((strcmp(ResidueGetName(pThis), "PHE") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "TYR") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "ASP") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "GLU") == 0 && j == 2) ||
        (strcmp(ResidueGetName(pThis), "GLN") == 0 && j == 2) ||
        (strcmp(ResidueGetName(pThis), "ASN") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "HSD") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "HSE") == 0 && j == 1) ||
        (strcmp(ResidueGetName(pThis), "ARG") == 0 && j == 3))
      {
        torsion2 = torsion + PI;
        torsion2 = torsion > 0 ? torsion - PI : torsion2;
      }
      double torsion2m2pi = torsion2 - 2 * PI;
      double torsion2p2pi = torsion2 + 2 * PI;
      if ((torsion <= max && torsion >= min) ||
        (torsionm2pi <= max && torsionm2pi >= min) ||
        (torsionp2pi <= max && torsionp2pi >= min) ||
        (torsion2 <= max && torsion2 >= min) ||
        (torsion2m2pi <= max && torsion2m2pi >= min) ||
        (torsion2p2pi <= max && torsion2p2pi >= min))
      {
        IntArrayAppend(pSimArray, 1);
        //IntArraySet(pSimArray,j,1);//1 stands for similarity
      }
      else
      {
        IntArrayAppend(pSimArray, 2);
        //IntArraySet(pSimArray,j,2);//2 stands for dissimilarity
        allmatch = FALSE;
      }

    }
  }
  else
  {
    allmatch = FALSE;
  }
  return allmatch;
}




BOOL ResidueAndResidueInSameType(Residue* pThis, Residue* pOther)
{
  if (strcmp(ResidueGetName(pThis), ResidueGetName(pOther)) == 0 ||
    (strcmp(ResidueGetName(pThis), "HSD") == 0 && strcmp(ResidueGetName(pOther), "HSE") == 0) ||
    (strcmp(ResidueGetName(pThis), "HSE") == 0 && strcmp(ResidueGetName(pOther), "HSD") == 0))
  {
    return TRUE;
  }
  return FALSE;
}


double ResidueGetAverageBfactor(Residue* pThis)
{
  double bfac = 0.0;
  int nhcount = 0;
  for (int i = 0;i < ResidueGetAtomCount(pThis);i++)
  {
    Atom* pAtom = ResidueGetAtom(pThis, i);
    if (AtomIsHydrogen(pAtom)) continue;
    if (strcmp(AtomGetName(pAtom), "OXT") == 0) continue;
    bfac += pAtom->bfactor;
    nhcount++;
  }
  return bfac / nhcount;
}


double ResidueGetDunbrack(Residue* pThis)
{
  return pThis->dunbrack;
}


int ResidueSetDunbrack(Residue* pThis, double dun)
{
  pThis->dunbrack = dun;
  return Success;
}


BOOL IsResidueHistidine(char* resiName)
{
  if (!strcmp(resiName, "HIS") || !strcmp(resiName, "HSD") || !strcmp(resiName, "HSE") || !strcmp(resiName, "HSP"))
  {
    return TRUE;
  }
  return FALSE;
}


