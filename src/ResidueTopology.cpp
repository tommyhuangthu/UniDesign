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

#include "ResidueTopology.h"
#include <string.h>

int BondCreate(Bond* pThis)
{
  strcpy(pThis->atomFromName, "");
  strcpy(pThis->atomToName, "");
  pThis->type = Type_Bond_None;
  return Success;
}

int BondDestroy(Bond* pThis)
{
  strcpy(pThis->atomFromName, "");
  strcpy(pThis->atomToName, "");
  pThis->type = Type_Bond_None;
  return Success;
}

int BondCopy(Bond* pThis, Bond* pOther)
{
  *pThis = *pOther;
  return Success;
}

char* BondGetFromName(Bond* pThis)
{
  return pThis->atomFromName;
}

int BondSetFromName(Bond* pThis, char* from)
{
  if (strlen(from) > MAX_LEN_ATOM_NAME)
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s line %d, name is too long", __FILE__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->atomFromName, from);
  return Success;
}

char* BondGetToName(Bond* pThis)
{
  return pThis->atomToName;
}

int BondSetToName(Bond* pThis, char* to)
{
  if (strlen(to) > MAX_LEN_ATOM_NAME)
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s line %d, name is too long", __FILE__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->atomToName, to);
  return Success;
}

Type_Bond BondGetType(Bond* pThis)
{
  return pThis->type;
}

int BondSetType(Bond* pThis, Type_Bond newType)
{
  if (newType<Type_Bond_Single || newType>Type_Bond_None)
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s line %d, name is too long", __FILE__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  pThis->type = newType;
  return Success;
}

int BondShow(Bond* pThis)
{
  char bondType[MAX_LEN_ONE_LINE_CONTENT+1];
  switch (BondGetType(pThis))
  {
  case Type_Bond_None:
    strcpy(bondType, "NONE  ");break;
  case Type_Bond_Single:
    strcpy(bondType, "BOND  ");break;
  case Type_Bond_double:
    strcpy(bondType, "DOUBLE");break;
  case Type_Bond_Triple:
    strcpy(bondType, "TRIPLE");break;
  default:
    strcpy(bondType, "");break;
  }
  printf("%s  %s  %s", bondType, BondGetFromName(pThis), BondGetToName(pThis));
  return Success;
}


int BondsCreate(Bonds* pThis)
{
  pThis->num = 0;
  pThis->bonds = NULL;
  return Success;
}

int BondsDestroy(Bonds* pThis)
{
  for (int i = 0;i < pThis->num;i++)
  {
    BondDestroy(&pThis->bonds[i]);
  }
  free(pThis->bonds);
  pThis->bonds = NULL;
  pThis->num = 0;
  return Success;
}

int BondsCopy(Bonds* pThis, Bonds* pOther)
{
  BondsDestroy(pThis);
  BondsCreate(pThis);
  pThis->num = pOther->num;
  pThis->bonds = (Bond*)malloc(sizeof(Bond) * pOther->num);
  for (int i = 0;i < pThis->num;i++)
  {
    BondCreate(&pThis->bonds[i]);
    BondCopy(&pThis->bonds[i], &pOther->bonds[i]);
  }
  return Success;
}

int BondsGetNum(Bonds* pThis)
{
  return pThis->num;
}

int BondsAdd(Bonds* pThis, char* atom1, char* atom2, Type_Bond bondType)
{
  if (BondsFind(pThis, atom1, atom2) == Type_Bond_None)
  {
    Bond* pNewlyAddedBond;
    (pThis->num)++;
    pThis->bonds = (Bond*)realloc(pThis->bonds, sizeof(Bond) * pThis->num);
    pNewlyAddedBond = &pThis->bonds[pThis->num - 1];
    BondCreate(pNewlyAddedBond);
    if (FAILED(BondSetFromName(pNewlyAddedBond, atom1)) || 
      FAILED(BondSetToName(pNewlyAddedBond, atom2)) || 
      FAILED(BondSetType(pNewlyAddedBond, bondType)))
    {
      return ValueError;
    }
  }
  return Success;
}

int BondsRemove(Bonds* pThis, char* atom1, char* atom2)
{
  for (int i = 0;i < pThis->num;i++)
  {
    Bond* pCurBond = &pThis->bonds[i];
    if ((strcmp(BondGetFromName(pCurBond), atom1) == 0 && strcmp(BondGetToName(pCurBond), atom2) == 0) || 
      (strcmp(BondGetFromName(pCurBond), atom2) == 0 && strcmp(BondGetToName(pCurBond), atom1) == 0))
    {
      BondCopy(pCurBond, &pThis->bonds[pThis->num - 1]);
      BondDestroy(&pThis->bonds[pThis->num - 1]);
      (pThis->num)--;
      return Success;
    }
  }
  return DataNotExistError;
}

Type_Bond BondsFind(Bonds* pThis, char* atom1, char* atom2)
{
  for (int i = 0;i < pThis->num;i++)
  {
    Bond* pCurBond = &pThis->bonds[i];
    if ((strcmp(BondGetFromName(pCurBond), atom1) == 0 && strcmp(BondGetToName(pCurBond), atom2) == 0) || 
      (strcmp(BondGetFromName(pCurBond), atom2) == 0 && strcmp(BondGetToName(pCurBond), atom1) == 0))
    {
      return BondGetType(pCurBond);
    }
  }
  return Type_Bond_None;
}

Bond* BondsGet(Bonds* pThis, int index)
{
  if (index < 0 || index >= pThis->num)
  {
    return NULL;
  }
  return &pThis->bonds[index];
}

int BondsShow(Bonds* pThis)
{
  for (int i = 0;i < pThis->num;i++)
  {
    BondShow(&pThis->bonds[i]);
    printf("\n");
  }
  return Success;
}


//------------Internal Coordinates-----------------------

int CharmmICCreate(CharmmIC* pThis)
{
  for (int i = 0;i < 4;i++)
    strcpy(pThis->atomNames[i], "");
  return Success;
}

int CharmmICCreateFromStringArray(CharmmIC* pThis, StringArray* params)
{
  CharmmICCreate(pThis);
  if (StringArrayGetNum(params) != 10)
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s line %d, error format", __FILE__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }

  for (int i = 0;i < 4;i++)
  {
    if (strlen(StringArrayGet(params, i + 1)) > MAX_LEN_ATOM_NAME)
    {
      char errMsg[MAX_LEN_ERR_MSG+1];
      int errorCode = ValueError;
      sprintf(errMsg, "in file %s line %d, name is too long", __FILE__, __LINE__);
      TraceError(errMsg, errorCode);
      return errorCode;
    }
    strcpy(pThis->atomNames[i], StringArrayGet(params, i + 1));
  }

  for (int i = 0;i < 5;i++)
  {
    pThis->icParam[i] = atof(StringArrayGet(params, i + 5));
    if (i >= 1 && i <= 3)
    {
      // covert degree to arc
      pThis->icParam[i] = DegToRad(pThis->icParam[i]);
    }
  }

  char buffer[MAX_LEN_ONE_LINE_CONTENT+1];
  strcpy(buffer, CharmmICGetAtomC(pThis));
  // deal with improper torsion angle
  if (buffer[0] == '*')
  {
    pThis->torFlag = FALSE;
    strcpy(pThis->atomNames[2], buffer + 1);
  }
  else
  {
    pThis->torFlag = TRUE;
  }
  return Success;
}

int CharmmICDestroy(CharmmIC* pThis)
{
  for (int i = 0;i < 4;i++)
    strcpy(pThis->atomNames[i], "");
  return Success;
}

int CharmmICCopy(CharmmIC* pThis, CharmmIC* pOther)
{
  *pThis = *pOther;
  return Success;
}

char* CharmmICGetAtomA(CharmmIC* pThis)
{
  return pThis->atomNames[0];
}

char* CharmmICGetAtomB(CharmmIC* pThis)
{
  return pThis->atomNames[1];
}

char* CharmmICGetAtomC(CharmmIC* pThis)
{
  return pThis->atomNames[2];
}

char* CharmmICGetAtomD(CharmmIC* pThis)
{
  return pThis->atomNames[3];
}

double* CharmmICGetICParams(CharmmIC* pThis)
{
  return pThis->icParam;
}

int CharmmICSetTorsionAngle(CharmmIC* pThis, double newTorsionAngle)
{
  pThis->icParam[2] = newTorsionAngle;
  return Success;
}

BOOL CharmmICGetTorsionProperFlag(CharmmIC* pThis)
{
  return pThis->torFlag;
}

int CharmmICCalcXYZ(CharmmIC* pThis, AtomArray* atomArray, XYZ* pDestXYZ)
{
  XYZ* pAtomsABC[3];
  for (int i = 0;i < 3;i++)
  {
    Atom* pAtom = AtomArrayGetByName(atomArray, pThis->atomNames[i]);
    if (pAtom == NULL)
    {
      return DataNotExistError;
    }

    if (pAtom->isXyzValid == FALSE)
    {
      return DataNotExistError;
    }

    pAtomsABC[i] = &pAtom->xyz;
  }
  return GetFourthAtom(pAtomsABC[0], pAtomsABC[1], pAtomsABC[2], pThis->icParam, pDestXYZ);
}

int CharmmICShow(CharmmIC* pThis)
{
  printf("IC  %-7.7s %-7.7s %c%-7.7s %-7.7s  ", 
    CharmmICGetAtomA(pThis), 
    CharmmICGetAtomB(pThis), 
    CharmmICGetTorsionProperFlag(pThis) ? ' ' : '*', 
    CharmmICGetAtomC(pThis), 
    CharmmICGetAtomD(pThis));
  printf("%8.3f ", pThis->icParam[0]);
  for (int i = 1;i < 4;i++)
  {
    printf("%8.3f ", RadToDeg(pThis->icParam[i]));
  }
  printf("%8.3f", pThis->icParam[4]);
  return Success;
}


int ResTopologyCreate(ResTopology* pThis)
{
  StringArrayCreate(&pThis->atoms);
  StringArrayCreate(&pThis->dels);
  BondsCreate(&pThis->bonds);
  pThis->icNum = 0;
  pThis->ics = NULL;
  return Success;
}

int ResTopologyCreateFromFileReader(ResTopology* pThis, FileReader* fr)
{
  BOOL resNameAlreadySet = FALSE;
  ResTopologyCreate(pThis);
  while (FileReaderEndOfFile(fr) == FALSE)
  {
    char line[MAX_LEN_ONE_LINE_CONTENT+1];
    char* keyword;
    StringArray strs;
    StringArrayCreate(&strs);

    FileReaderGetNextLine(fr, line);
    StringArraySplitString(&strs, line, ' ');
    keyword = StringArrayGet(&strs, 0);

    if (strcmp(keyword, "RESI") == 0 || strcmp(keyword, "PRES") == 0)
    {
      if (resNameAlreadySet == TRUE)
      {
        FileReaderSetCurrentPos(fr, FileReaderGetCurrentPos(fr) - 1);
        StringArrayDestroy(&strs);
        return Success;
      }
      else
      {
        char* resName = StringArrayGet(&strs, 1);
        if (resName == NULL || FAILED(ResTopologySetName(pThis, resName)))
        {
          char errMsg[MAX_LEN_ERR_MSG+1];
          int errCode = ValueError;
          sprintf(errMsg, "in file %s line %d, failed to set residue name", __FILE__, __LINE__);
          TraceError(errMsg, errCode);
          return errCode;
        }
        resNameAlreadySet = TRUE;
      }
    }
    else if (strcmp(keyword, "ATOM") == 0)
    {
      char* atomName = StringArrayGet(&strs, 1);
      if (atomName == NULL || FAILED(StringArrayAppend(&pThis->atoms, atomName)))
      {
        char errMsg[MAX_LEN_ERR_MSG+1];
        int errCode = ValueError;
        sprintf(errMsg, "in file %s line %d, failed to add atom name", __FILE__, __LINE__);
        TraceError(errMsg, errCode);
        return errCode;
      }
    }
    else if (strcmp(keyword, "DELETE") == 0)
    {
      char* delAtomName = StringArrayGet(&strs, 2);
      if (delAtomName == NULL || FAILED(StringArrayAppend(&pThis->dels, delAtomName)))
      {
        char errMsg[MAX_LEN_ERR_MSG+1];
        int errCode = ValueError;
        sprintf(errMsg, "in file %s line %d, failed to add atom name", __FILE__, __LINE__);
        TraceError(errMsg, errCode);
        return errCode;
      }
    }
    else if (strcmp(keyword, "BOND") == 0 || strcmp(keyword, "DOUBLE") == 0 || strcmp(keyword, "TRIPLE") == 0)
    {
      for (int i = 1;i < StringArrayGetNum(&strs);i += 2)
      {
        char* atomFromName = StringArrayGet(&strs, i);
        char* atomToName = StringArrayGet(&strs, i + 1);
        Type_Bond bondType;

        if (atomFromName == NULL || atomToName == NULL)
        {
          char errMsg[MAX_LEN_ERR_MSG+1];
          int errCode = ValueError;
          sprintf(errMsg, "in file %s line %d, atom %s or %s cannot be found", __FILE__, __LINE__, atomFromName, atomToName);
          TraceError(errMsg, errCode);
          return errCode;
        }
        if (strcmp(keyword, "BOND") == 0)
        {
          bondType = Type_Bond_Single;
        }
        else if (strcmp(keyword, "DOUBLE") == 0)
        {
          bondType = Type_Bond_double;
        }
        else if (strcmp(keyword, "TRIPLE") == 0)
        {
          bondType = Type_Bond_Triple;
        }

        if (FAILED(BondsAdd(&pThis->bonds, atomFromName, atomToName, bondType)))
        {
          char errMsg[MAX_LEN_ERR_MSG+1];
          int errCode = ValueError;
          sprintf(errMsg, "in file %s line %d, cannot add bond to set", __FILE__, __LINE__);
          TraceError(errMsg, errCode);
          return errCode;
        }
      }
    }
    else if (strcmp(keyword, "IC") == 0)
    {
      CharmmIC ic;
      if (FAILED(CharmmICCreateFromStringArray(&ic, &strs)) || FAILED(ResTopologyAddCharmmIC(pThis, &ic)))
      {
        char errMsg[MAX_LEN_ERR_MSG+1];
        int errCode = ValueError;
        sprintf(errMsg, "in file %s line %d, IC error", __FILE__, __LINE__);
        TraceError(errMsg, errCode);
        return errCode;
      }
      CharmmICDestroy(&ic);
    }
    else
    {
      ;
    }
    StringArrayDestroy(&strs);
  }
  return Success;
}

int ResTopologyDestroy(ResTopology* pThis)
{
  StringArrayDestroy(&pThis->atoms);
  StringArrayDestroy(&pThis->dels);
  BondsDestroy(&pThis->bonds);

  for (int i = 0;i < pThis->icNum;i++)
  {
    CharmmICDestroy(&pThis->ics[i]);
  }
  free(pThis->ics);
  pThis->ics = NULL;
  pThis->icNum = 0;
  return Success;
}

int ResTopologyCopy(ResTopology* pThis, ResTopology* pOther)
{
  ResTopologyDestroy(pThis);
  ResTopologyCreate(pThis);
  ResTopologySetName(pThis, ResTopologyGetName(pOther));
  StringArrayCopy(&pThis->atoms, &pOther->atoms);
  StringArrayCopy(&pThis->dels, &pOther->dels);
  BondsCopy(&pThis->bonds, &pOther->bonds);
  pThis->icNum = pOther->icNum;
  pThis->ics = (CharmmIC*)malloc(sizeof(CharmmIC) * pThis->icNum);
  for (int i = 0;i < pThis->icNum;i++)
  {
    CharmmICCreate(&pThis->ics[i]);
    CharmmICCopy(&pThis->ics[i], &pOther->ics[i]);
  }
  return Success;
}

char* ResTopologyGetName(ResTopology* pThis)
{
  return pThis->resName;
}

int ResTopologySetName(ResTopology* pThis, char* newName)
{
  if (strlen(newName) > MAX_LEN_RES_NAME)
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(errMsg, "in file %s line %d, name is too long", __FILE__, __LINE__);
    TraceError(errMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->resName, newName);
  return Success;
}

StringArray* ResTopologyGetAtoms(ResTopology* pThis)
{
  return &pThis->atoms;
}

StringArray* ResTopologyGetDeletes(ResTopology* pThis)
{
  return &pThis->dels;
}

Bonds* ResTopologyGetBonds(ResTopology* pThis)
{
  return &pThis->bonds;
}

int ResTopologyGetCharmmICNum(ResTopology* pThis)
{
  return pThis->icNum;
}

int ResTopologyGetCharmmIC(ResTopology* pThis, int ndx, CharmmIC* ic)
{
  if (ndx < 0 || ndx >= pThis->icNum)
    return IndexError;
  CharmmICCopy(ic, &pThis->ics[ndx]);
  return Success;
}

int ResTopologyFindCharmmIC(ResTopology* pThis, char* atomDName, CharmmIC* ic)
{
  for (int i = 0;i < pThis->icNum;i++)
  {
    if (strcmp(CharmmICGetAtomD(&pThis->ics[i]), atomDName) == 0)
    {
      CharmmICCopy(ic, &pThis->ics[i]);
      return Success;
    }
  }
  return DataNotExistError;
}

int ResTopologyFindCharmmICIndex(ResTopology* pThis, char* atomName, int* ndx)
{
  for (int i = 0;i < pThis->icNum;i++)
  {
    if (strcmp(CharmmICGetAtomD(&pThis->ics[i]), atomName) == 0)
    {
      *ndx = i;
      return Success;
    }
  }
  return DataNotExistError;
}

int ResTopologyAddCharmmIC(ResTopology* pThis, CharmmIC* ic)
{
  int num = pThis->icNum + 1;
  pThis->ics = (CharmmIC*)realloc(pThis->ics, sizeof(CharmmIC) * num);
  CharmmICCreate(&pThis->ics[num - 1]);
  CharmmICCopy(&pThis->ics[num - 1], ic);
  pThis->icNum = num;
  return Success;
}

int ResTopologyShow(ResTopology* pThis)
{
  printf("RESI  %s\n", ResTopologyGetName(pThis));
  for (int i = 0;i < StringArrayGetNum(&pThis->atoms);i++)
  {
    printf("ATOM   %s\n", StringArrayGet(&pThis->atoms, i));
  }
  BondsShow(&pThis->bonds);
  for (int i = 0;i < pThis->icNum;i++)
  {
    CharmmICShow(&pThis->ics[i]);
    printf("\n");
  }
  return Success;
}


int ResToposCreate(ResTopos* pThis)
{
  pThis->num = 0;
  pThis->topos = NULL;
  return Success;
}

int ResToposDestroy(ResTopos* pThis)
{
  for (int i = 0;i < pThis->num;i++)
  {
    ResTopologyDestroy(&pThis->topos[i]);
  }
  free(pThis->topos);
  pThis->topos = NULL;
  pThis->num = 0;
  return Success;
}

int ResToposCopy(ResTopos* pThis, ResTopos* pOther)
{
  ResToposDestroy(pThis);
  pThis->num = pOther->num;
  pThis->topos = (ResTopology*)malloc(sizeof(ResTopology) * pThis->num);
  for (int i = 0;i < pThis->num;i++)
  {
    ResTopologyCreate(&pThis->topos[i]);
    ResTopologyCopy(&pThis->topos[i], &pOther->topos[i]);
  }
  return Success;
}

int ResToposGet(ResTopos* pThis, char* resiName, ResTopology* pDestTopo)
{
  for (int i = 0;i < pThis->num;i++)
  {
    if (strcmp(ResTopologyGetName(&pThis->topos[i]), resiName) == 0)
    {
      ResTopologyCopy(pDestTopo, &pThis->topos[i]);
      return Success;
    }
  }
  return DataNotExistError;
}

int ResToposFindIndexByName(ResTopos* pThis, char* resiName, int* index)
{
  for (int i = 0;i < pThis->num;i++)
  {
    if (strcmp(ResTopologyGetName(&pThis->topos[i]), resiName) == 0)
    {
      *index = i;
      return Success;
    }
  }
  return DataNotExistError;
}

int ResToposAdd(ResTopos* pThis, ResTopology* pNewTopo)
{
  for (int i = 0;i < pThis->num;i++)
  {
    if (strcmp(ResTopologyGetName(&pThis->topos[i]), ResTopologyGetName(pNewTopo)) == 0)
    {
      ResTopologyCopy(&pThis->topos[i], pNewTopo);
      return Success;
    }
  }

  int num = pThis->num + 1;
  pThis->topos = (ResTopology*)realloc(pThis->topos, sizeof(ResTopology) * num);
  ResTopologyCreate(&pThis->topos[num - 1]);
  ResTopologyCopy(&pThis->topos[num - 1], pNewTopo);
  pThis->num = num;
  return Success;
}

int ResToposAddFromFile(ResTopos* pThis, char* filepath)
{
  FileReader file;
  int result = FileReaderCreate(&file, filepath);
  if (FAILED(result))
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    sprintf(errMsg, "in file %s line %d, cannot read file %s", __FILE__, __LINE__, filepath);
    TraceError(errMsg, result);
    return result;
  }

  while (!FileReaderEndOfFile(&file))
  {
    ResTopology newTopo;
    result = ResTopologyCreateFromFileReader(&newTopo, &file);
    if (FAILED(result))
    {
      char errMsg[MAX_LEN_ERR_MSG+1];
      sprintf(errMsg, "in file %s line %d, cannot read file %s", __FILE__, __LINE__, filepath);
      TraceError(errMsg, result);
      return result;
    }
    ResToposAdd(pThis, &newTopo);
    ResTopologyDestroy(&newTopo);
  }
  FileReaderDestroy(&file);
  return Success;
}


int ResToposRead(ResTopos* pResiTopo, char* filePath)
{
  if (filePath == NULL || FAILED(ResToposAddFromFile(pResiTopo, filePath)))
  {
    return FormatError;
  }
  return Success;
}


int ResToposShow(ResTopos* pThis)
{
  for (int i = 0;i < pThis->num;i++)
  {
    printf("Residue No. %d\n", i);
    ResTopologyShow(&pThis->topos[i]);
  }
  return Success;
}
