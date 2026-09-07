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

#include "SmallMol.h"
#include <string.h>

#include <ctype.h>
#include <time.h>

//#define DEBUGGING_SMALLMOL

int CataConsItemShow(CataConsItem* pThis)
{
  char type[MAX_LEN_ONE_LINE_CONTENT+1];
  char location[MAX_LEN_ONE_LINE_CONTENT+1];

  switch (pThis->type)
  {
  case Type_CataCons_NewPseudoAtom:
    strcpy(type, "PSEUDO_ATOM");break;
  case Type_CataCons_Distance:
    strcpy(type, "DISTANCE   ");break;
  case Type_CataCons_Angle:
    strcpy(type, "ANGLE      ");break;
  case Type_CataCons_Torsion:
    strcpy(type, "TORSION    ");break;
  }
  printf("%s ", type);
  for (int i = 0;i < 4;i++)
  {
    printf("%4.4s", pThis->atomName[i]);
    switch (pThis->atomLoc[i])
    {
    case Type_CataConsAtomLocation_Pseudo:
      strcpy(location, "(PSU)");break;
    case Type_CataConsAtomLocation_FirstSite:
      strcpy(location, "(1st)");break;
    case Type_CataConsAtomLocation_SecondSite:
      strcpy(location, "(2nd)");break;
    default:
      strcpy(location, "(UKN)");break;
    }
    printf("%s ", location);
    printf("%2d ", pThis->atomIndex[i]);
  }
  printf("%.3f %.3f\n", pThis->min, pThis->max);
  return Success;
}


int CataConsGroupDeploy(CataConsGroup* pThis, Rotamer* pFirstSiteRotamer, Rotamer* pSecondSiteRotamer)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  int  result = Success;
  StringArray pseudoAtomNames;
  StringArrayCreate(&pseudoAtomNames);

  for (int itemIndex = 0; itemIndex < pThis->nconsItems; itemIndex++)
  {
    CataConsItem* pCurItem = &pThis->consItems[itemIndex];
    if (pCurItem->type == Type_CataCons_NewPseudoAtom)
    {
      StringArrayAppend(&pseudoAtomNames, pCurItem->atomName[0]);
    }
  }

  for (int itemIndex = 0; itemIndex < pThis->nconsItems; itemIndex++)
  {
    int atomNum;
    CataConsItem* pCurItem = &pThis->consItems[itemIndex];
    switch (pCurItem->type)
    {
    case Type_CataCons_NewPseudoAtom:
      atomNum = 3;
      break;
    case Type_CataCons_Distance:
      atomNum = 2;
      break;
    case Type_CataCons_Angle:
      atomNum = 3;
      break;
    case Type_CataCons_Torsion:
      atomNum = 4;
      break;
    default:
      atomNum = 4;
    }

    for (int i = 0;i < atomNum;i++)
    {
      switch (pCurItem->atomLoc[i])
      {
      case Type_CataConsAtomLocation_FirstSite:
        result = RotamerFindAtom(pFirstSiteRotamer, pCurItem->atomName[i], &pCurItem->atomIndex[i]);
        if (FAILED(result))
        {
          sprintf(errMsg, "in file %s line %d, cannot find atom %s on rotamer %s of the first site", __FILE__, __LINE__, pCurItem->atomName[i], pThis->site1RotType);
          TraceError(errMsg, result);
          return result;
        }
        break;

      case Type_CataConsAtomLocation_SecondSite:
        result = RotamerFindAtom(pSecondSiteRotamer, pCurItem->atomName[i], &pCurItem->atomIndex[i]);
        if (FAILED(result))
        {
          sprintf(errMsg, "in file %s line %d, cannot find atom %s on rotamer %s of the second site", __FILE__, __LINE__, pCurItem->atomName[i], pThis->site2RotType);
          TraceError(errMsg, result);
          return result;
        }
        break;

      case Type_CataConsAtomLocation_Pseudo:
        result = StringArrayFind(&pseudoAtomNames, pCurItem->atomName[i], &pCurItem->atomIndex[i]);
        if (FAILED(result))
        {
          sprintf(errMsg, "in file %s line %d, cannot find atom %s in the pseudo atom set", __FILE__, __LINE__, pCurItem->atomName[i]);
          TraceError(errMsg, result);
          return result;
        }
        break;

      default:
        sprintf(errMsg, "in file %s line %d, unspecified location for atom %s", __FILE__, __LINE__, pCurItem->atomName[i]);
        TraceError(errMsg, result);
        return result;
      }
    }
  }

  StringArrayDestroy(&pseudoAtomNames);

  return Success;
}


int CataConsSitePairCreate(CataConsSitePair* pThis, FileReader* pFile)
{
  char line[MAX_LEN_ONE_LINE_CONTENT+1];
  char keyword[MAX_LEN_ONE_LINE_CONTENT+1];
  char errMsg[MAX_LEN_ERR_MSG+1];
  char atomName[MAX_LEN_ONE_LINE_CONTENT+1];
  int  result;

  CataConsItem* pCurItem;
  CataConsGroup* pCurGroup;

  BOOL done;

  //the deployed flad is initialized as FALSE;
  pThis->deployedFlag = FALSE;
  pThis->groups = NULL;
  pThis->ngroups = 0;
  //initialize the constraint type to be HBbond
  pThis->pairConsType = Type_SitePairCons_Hbond;

  done = FALSE;
  pCurGroup = NULL;
  while (!FAILED(FileReaderGetNextLine(pFile, line)))
  {
    ExtractFirstStringFromSourceString(keyword, line);

    // when the keyword 'SITEPAIR' is first met, the value of flag is false;
    // and the information following the 'SITEPAIR' should be read, and the flag is set as true.
    // when the keyword 'SITEPAIR' is met again, the value of flag is true;
    // the position of filereader should be reset and exit this function;
    if (strcmp(keyword, "SITEPAIR") == 0 && done == FALSE)
    {
      sscanf(line, "%s %s %d %s %s %d", pThis->chnName1, pThis->resName1, &pThis->pos1, pThis->chnName2, pThis->resName2, &pThis->pos2);
      done = TRUE;
    }
    else if (strcmp(keyword, "SITEPAIR") == 0 && done)
    {
      FileReaderSetCurrentPos(pFile, FileReaderGetCurrentPos(pFile) - 1);
      return Success;
    }
    else if (strcmp(keyword, "TYPE") == 0)
    {
      ExtractFirstStringFromSourceString(keyword, line);
      if (strcmp(keyword, "COVALENT_BOND") == 0) pThis->pairConsType = Type_SitePairCons_Covalent;
      else if (strcmp(keyword, "HBOND") == 0) pThis->pairConsType = Type_SitePairCons_Hbond;
      else if (strcmp(keyword, "SALT_BRIDGE") == 0) pThis->pairConsType = Type_SitePairCons_SaltBridge;
      else if (strcmp(keyword, "VDW") == 0) pThis->pairConsType = Type_SitePairCons_VDW;
      else if (strcmp(keyword, "PPSTACK") == 0) pThis->pairConsType = Type_SitePairCons_PiPiStack;
      else if (strcmp(keyword, "TSTACK") == 0) pThis->pairConsType = Type_SitePairCons_TStack;
      else pThis->pairConsType = Type_SitePairCons_Other;
    }
    else if (strcmp(keyword, "GROUP") == 0)
    {
      //Create New Group
      (pThis->ngroups)++;
      pThis->groups = (CataConsGroup*)realloc(pThis->groups, sizeof(CataConsGroup) * pThis->ngroups);
      pCurGroup = &pThis->groups[pThis->ngroups - 1];
      sscanf(line, "%s %s", pCurGroup->site1RotType, pCurGroup->site2RotType);
      pCurGroup->nconsItems = 0;
      pCurGroup->consItems = NULL;
      pCurGroup->npseudoAtoms = 0;
      pCurGroup->pseudoAtoms = NULL;
    }
    else if (strcmp(keyword, "ITEM") == 0)
    {
      //Create New Item
      int atomIndex;
      char constraintType[MAX_LEN_ONE_LINE_CONTENT+1];

      // allocate memories for new catalytic constraint items;
      (pCurGroup->nconsItems)++;
      pCurGroup->consItems = (CataConsItem*)realloc(pCurGroup->consItems,
        sizeof(CataConsItem) * pCurGroup->nconsItems);
      pCurItem = &pCurGroup->consItems[pCurGroup->nconsItems - 1];

      // read the catalytic constraint type;
      if (FAILED(ExtractFirstStringFromSourceString(constraintType, line)))
      {
        result = FormatError;
        sprintf(errMsg, "in file %s line %d, unrecognized format in line:\n%s", __FILE__, __LINE__, line);
        TraceError(errMsg, result);
        return result;
      }

      if (strcmp(constraintType, "PSEUDO_ATOM") == 0)
      {
        pCurItem->type = Type_CataCons_NewPseudoAtom;
        (pCurGroup->npseudoAtoms)++;
        pCurGroup->pseudoAtoms = (XYZ*)realloc(pCurGroup->pseudoAtoms,
          sizeof(XYZ) * pCurGroup->npseudoAtoms);
      }
      else if (strcmp(constraintType, "DISTANCE") == 0)
      {
        pCurItem->type = Type_CataCons_Distance;
      }
      else if (strcmp(constraintType, "ANGLE") == 0)
      {
        pCurItem->type = Type_CataCons_Angle;
      }
      else if (strcmp(constraintType, "TORSION") == 0)
      {
        pCurItem->type = Type_CataCons_Torsion;
      }
      else
      {
        result = FormatError;
        sprintf(errMsg, "in file %s line %d, unrecognized keyword %s", __FILE__, __LINE__, constraintType);
        TraceError(errMsg, result);
        return result;
      }

      // read the four atoms for each atom;
      for (atomIndex = 0; atomIndex < 4; atomIndex++)
      {
        if (FAILED(ExtractFirstStringFromSourceString(atomName, line)))
        {
          result = FormatError;
          sprintf(errMsg, "in file %s line %d, unrecognized format in line:\n%s", __FILE__, __LINE__, line);
          TraceError(errMsg, result);
          return result;
        }
        if (atomName[0] == SECOND_CATACON_SITE_SYMBOL)
        {
          pCurItem->atomLoc[atomIndex] = Type_CataConsAtomLocation_SecondSite;
          // remove the symbol '*' before the atom name;
          strcpy(pCurItem->atomName[atomIndex], atomName+1);
        }
        else if (atomName[0] == PSEUDO_ATOM_SYMBOL)
        {
          pCurItem->atomLoc[atomIndex] = Type_CataConsAtomLocation_Pseudo;
          // remove the symbol '+' before the atom name;
          strcpy(pCurItem->atomName[atomIndex], atomName+1);
        }
        else
        {
          pCurItem->atomLoc[atomIndex] = Type_CataConsAtomLocation_FirstSite;
          strcpy(pCurItem->atomName[atomIndex], atomName);
        }
      }

      // read the range for constraint;
      sscanf(line, "%lf %lf", &pCurItem->min, &pCurItem->max);

      if (pCurItem->type == Type_CataCons_Angle ||
        pCurItem->type == Type_CataCons_Torsion)
      {
        pCurItem->min = DegToRad(pCurItem->min);
        pCurItem->max = DegToRad(pCurItem->max);
      }

      if (pCurItem->type != Type_CataCons_NewPseudoAtom)
      {
        // set a disturbance for the upper and lower bound;
        pCurItem->min -= fabs(DISTURBANCE_IN_RANGE_CHECK);
        pCurItem->max += fabs(DISTURBANCE_IN_RANGE_CHECK);
      }

      // initialize the index for the related four atoms, and set the index as -1;
      for (atomIndex = 0; atomIndex < 4; atomIndex++)
      {
        pCurItem->atomIndex[atomIndex] = -1;
      }
    }

    else
    {
      ;
    }

  }
  if (done)
  {
    return Success;
  }
  else
  {
    return ValueError;
  }

}


void CataConsSitePairDestroy(CataConsSitePair* pThis)
{
  for (int i = 0;i < pThis->ngroups;i++)
  {
    free(pThis->groups[i].consItems);
    free(pThis->groups[i].pseudoAtoms);
  }
  free(pThis->groups);
}


BOOL CataConsSitePairGetDeployedFlag(CataConsSitePair* pThis)
{
  return pThis->deployedFlag;
}


int CataConsSitePairDeploy(CataConsSitePair* pThis, Rotamers* pFirstSiteRotSet, Rotamers* pSecondSiteRotSet)
{
  int result;
  char errMsg[MAX_LEN_ERR_MSG+1];

  if (CataConsSitePairGetDeployedFlag(pThis))
  {
    result = ValueError;
    sprintf(errMsg, "in file %s line %d, catalytic constraint already deployed; if this pair of constraint is duplicated, please merge it with the existing one", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }


  for (int i = 0;i < pThis->ngroups;i++)
  {
    Rotamer* pFirstSiteRotamer = NULL;
    Rotamer* pSecondSiteRotamer = NULL;
    char* firstSiteRotamerType = pThis->groups[i].site1RotType;
    char* secondSiteRotamerType = pThis->groups[i].site2RotType;

    pFirstSiteRotamer = RotamersFindRep(pFirstSiteRotSet, firstSiteRotamerType);
    pSecondSiteRotamer = RotamersFindRep(pSecondSiteRotSet, secondSiteRotamerType);

    if (pFirstSiteRotamer == NULL)
    {
      result = DataNotExistError;
      sprintf(errMsg, "in file %s line %d, no %s rotamer exists on the 1st site, catalytic constraint may not apply", __FILE__, __LINE__, firstSiteRotamerType);
      TraceError(errMsg, result);
      return result;

    }
    if (pSecondSiteRotamer == NULL)
    {
      result = DataNotExistError;
      sprintf(errMsg, "in file %s line %d, no %s rotamer exists on the 2nd site, catalytic constraint may not apply", __FILE__, __LINE__, secondSiteRotamerType);
      TraceError(errMsg, result);
      return result;
    }

    result = CataConsGroupDeploy(&pThis->groups[i], pFirstSiteRotamer, pSecondSiteRotamer);
    if (FAILED(result))
    {
      sprintf(errMsg, "in file %s line %d, failed deploy", __FILE__, __LINE__);
      TraceError(errMsg, result);
      return result;
    }
  }

  pThis->deployedFlag = TRUE;

  return Success;
}


int CataConsSitePairShow(CataConsSitePair* pThis)
{
  CataConsGroup* pCurGroup = NULL;
  printf("SITE    %s   %s   %d    %s   %s   %d\n",
    pThis->chnName1, pThis->resName1, pThis->pos1,
    pThis->chnName2, pThis->resName2, pThis->pos2);

  for (int groupIndex = 0; groupIndex < pThis->ngroups; groupIndex++)
  {
    pCurGroup = &pThis->groups[groupIndex];
    printf("GROUP    %s    %s\n", pCurGroup->site1RotType, pCurGroup->site2RotType);
    for (int itemIndex = 0; itemIndex < pCurGroup->nconsItems; itemIndex++)
    {
      CataConsItemShow(&pCurGroup->consItems[itemIndex]);
    }
  }
  return Success;
}


int CataConsSitePairArrayCreate(CataConsSitePairArray* pThis, char* cataConsFile)
{
  FileReader file;
  int result = FileReaderCreate(&file, cataConsFile);
  if (FAILED(result))
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    sprintf(errMsg, "in file %s line %d, cannot open file", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }
  pThis->num = 0;
  pThis->sites = NULL;

  while (TRUE)
  {
    int newNum = pThis->num + 1;
    CataConsSitePair* pCurSite;
    pThis->sites = (CataConsSitePair*)realloc(pThis->sites, sizeof(CataConsSitePair) * newNum);
    pCurSite = &pThis->sites[newNum - 1];
    result = CataConsSitePairCreate(pCurSite, &file);
    if (FAILED(result))
    {
      break;
    }
    else
    {
      pThis->num = newNum;
    }
  }
  FileReaderDestroy(&file);
  return Success;
}


void CataConsSitePairArrayDestroy(CataConsSitePairArray* pThis)
{
  for (int i = 0;i < pThis->num;i++)
  {
    CataConsSitePairDestroy(&pThis->sites[i]);
  }
  free(pThis->sites);
}


int CataConsSitePairArrayGetNum(CataConsSitePairArray* pThis)
{
  return pThis->num;
}


CataConsSitePair* CataConsSitePairArrayGet(CataConsSitePairArray* pThis, int index)
{
  if (index<0 || index>pThis->num)
  {
    return NULL;
  }
  return &pThis->sites[index];
}


CataConsSitePair* CataConsSitePairArrayFind(CataConsSitePairArray* pThis, char* chnName1, int pos1, char* resiName1, char* chnName2, int pos2, char* resiName2)
{
  for (int i = 0;i < pThis->num;i++)
  {
    if (pos1 != pThis->sites[i].pos1) continue;
    if (pos2 != pThis->sites[i].pos2) continue;
    if (strcmp(chnName1, pThis->sites[i].chnName1) != 0) continue;
    if (strcmp(chnName2, pThis->sites[i].chnName2) != 0) continue;
    return &pThis->sites[i];
  }
  return NULL;
}


int CataConsSitePairArrayShow(CataConsSitePairArray* pThis)
{
  for (int i = 0;i < CataConsSitePairArrayGetNum(pThis);i++)
  {
    CataConsSitePairShow(CataConsSitePairArrayGet(pThis, i));
    printf("\n");
  }
  return Success;
}


int CataConsSitePairArrayTester(char* file)
{
  CataConsSitePairArray cataCons;
  CataConsSitePairArrayCreate(&cataCons, file);
  CataConsSitePairArrayShow(&cataCons);
  for (int i = 0;i < 100;i++)
  {
    CataConsSitePairArrayDestroy(&cataCons);
    CataConsSitePairArrayCreate(&cataCons, file);
  }

  CataConsSitePairArrayShow(&cataCons);
  CataConsSitePairArrayDestroy(&cataCons);
  return Success;
}


BOOL CataConsItemCheck(CataConsItem* pThis, XYZArray* pOnFirstSite, XYZArray* pOnSecondSite, XYZ* pseudoAtoms)
{
  double value;
  XYZ* pAtoms[4];
  Atom* pAtomOnProtein;
  Atom* pAtomOnRotamer;
  XYZ* pPseudoAtom;
  XYZ vectorI, vectorJ;

  pAtomOnProtein = NULL;
  pAtomOnRotamer = NULL;
  pPseudoAtom = NULL;

  for (int i = 0;i < 4;i++)
  {
    int atomIndex = pThis->atomIndex[i];
    switch (pThis->atomLoc[i])
    {
    case Type_CataConsAtomLocation_FirstSite:
      pAtoms[i] = XYZArrayGet(pOnFirstSite, atomIndex);
      break;
    case Type_CataConsAtomLocation_SecondSite:
      pAtoms[i] = XYZArrayGet(pOnSecondSite, atomIndex);
      break;
    case Type_CataConsAtomLocation_Pseudo:
      pPseudoAtom = &pseudoAtoms[atomIndex];
      pAtoms[i] = pPseudoAtom;
      break;
    case Type_CataConsAtomLocation_Undefined:
      break;
    }
  }

  switch (pThis->type)
  {
  case Type_CataCons_NewPseudoAtom:
    pPseudoAtom->X = pThis->min * (pAtoms[1]->X) + pThis->max * (pAtoms[2]->X);
    pPseudoAtom->Y = pThis->min * (pAtoms[1]->Y) + pThis->max * (pAtoms[2]->Y);
    pPseudoAtom->Z = pThis->min * (pAtoms[1]->Z) + pThis->max * (pAtoms[2]->Z);
    return TRUE;

  case Type_CataCons_Distance:
    value = XYZDist(pAtoms[1], pAtoms[0]);
    //printf("dist = %f\n", value);
    if (value > pThis->max || value < pThis->min)
    {
      return FALSE;
    }
    return TRUE;

  case Type_CataCons_Angle:
    vectorI = XYZDiff(pAtoms[0], pAtoms[1]);
    vectorJ = XYZDiff(pAtoms[2], pAtoms[1]);
    value = XYZAngle(&vectorI, &vectorJ);
    //printf("angle = %f\n", RadToDeg(value));
    if (value > pThis->max || value < pThis->min)
    {
      return FALSE;
    }
    return TRUE;

  case Type_CataCons_Torsion:
    value = GetTorsionAngle(pAtoms[0], pAtoms[1], pAtoms[2], pAtoms[3]);
    //printf("torsion = %f\n", RadToDeg(value));
    return RadInRange(value, pThis->min, pThis->max);
  }

  return FALSE;
}


BOOL CataConsGroupCheck(CataConsGroup* pThis, Rotamer* pOnFirstSite, Rotamer* pOnSecondSite)
{
  if (strcmp(pOnFirstSite->type, pThis->site1RotType) != 0 || strcmp(pOnSecondSite->type, pThis->site2RotType) != 0)
  {
    return FALSE;
  }
  for (int i = 0;i < pThis->nconsItems;i++)
  {
    if (!CataConsItemCheck(&pThis->consItems[i], &pOnFirstSite->xyzs, &pOnSecondSite->xyzs, pThis->pseudoAtoms))
    {
      return FALSE;
    }
  }
  return TRUE;
}


BOOL CataConsSitePairCheck(CataConsSitePair* pThis, Rotamer* pOnFirstSite, Rotamer* pOnSecondSite)
{
  if (pThis->deployedFlag == FALSE)
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    sprintf(errMsg, "in file %s line %d, catalytic constraint pair between %s%d%s and %s%d%s has not been deployed yet", __FILE__, __LINE__,
      pThis->chnName1, pThis->pos1, pThis->resName1,
      pThis->chnName2, pThis->pos2, pThis->resName2);
    TraceError(errMsg, ValueError);
    return FALSE;
  }
  for (int i = 0;i < pThis->ngroups;i++)
  {
    if (CataConsGroupCheck(&pThis->groups[i], pOnFirstSite, pOnSecondSite))
    {
      return TRUE;
    }
  }
  return FALSE;
}


//////////////////////////////////////////////////////////////////////////
//Methods of Placing Rule                                             
//////////////////////////////////////////////////////////////////////////
//Various Actions In Placing SmallMol
int PlacingActionCreate(PlacingAction* pThis, Type_PlacingAction type)
{
  pThis->actionType = type;
  switch (type)
  {
  case Type_PlacingAction_Load:
    IntArrayCreate(&pThis->load_atoms, 0);
    break;
  case Type_PlacingAction_Evaluate:
    IntArrayCreate(&pThis->evaluate_atoms, 0);
    break;
  case Type_PlacingAction_Calc:
    IntArrayCreate(&pThis->calc_atoms, 0);
    IntArrayCreate(&pThis->calc_params, 0);
    break;
  case Type_PlacingAction_CheckVDW_Backbone:
    IntArrayCreate(&pThis->checkVDW_backbone_smallmolAtomHasXyz, 0);
    break;
  case Type_PlacingAction_CheckVDW_Internal:
    pThis->checkVDW_internal_smallMolAtomNum = 0;
    pThis->checkVDW_internal_smallMolAtom13bondedMatrix = NULL;
    IntArrayCreate(&pThis->checkVDW_internal_smallmolAtomHasXyz, 0);
    break;
  case Type_PlacingAction_CheckRMSD:
    IntArrayCreate(&pThis->checkRMSD_smallmolAtomHasXyz, 0);
    break;
  case Type_PlacingAction_CheckMultiCons:
    pThis->checkMultiCons_cataConsNum = 0;
    StringArrayCreate(&pThis->checkMultiCons_firstSiteChainNames);
    StringArrayCreate(&pThis->checkMultiCons_firstSiteResidueNames);
    IntArrayCreate(&pThis->checkMultiCons_firstSitePosInChains, 0);
    StringArrayCreate(&pThis->checkMultiCons_secondSiteChainNames);
    StringArrayCreate(&pThis->checkMultiCons_secondSiteResidueNames);
    IntArrayCreate(&pThis->checkMultiCons_secondSitePosInChains, 0);

    for (int i = 0;i < MAX_NUM_CHECK_MULTI_CONS;i++)
    {
      pThis->checkMultiCons_rotNumOnEachSite[i][0] =
        pThis->checkMultiCons_rotNumOnEachSite[i][1] = 0;
      pThis->checkMultiCons_predeterminedConsRelations[i] = NULL;
    }

    break;
  default:
    break;
  }
  return Success;
}


int PlacingActionDestroy(PlacingAction* pThis)
{
  switch (pThis->actionType)
  {
  case Type_PlacingAction_Load:
    IntArrayDestroy(&pThis->load_atoms);
    break;
  case Type_PlacingAction_Evaluate:
    IntArrayDestroy(&pThis->evaluate_atoms);
    break;
  case Type_PlacingAction_Calc:
    IntArrayDestroy(&pThis->calc_atoms);
    IntArrayDestroy(&pThis->calc_params);
    break;
  case Type_PlacingAction_CheckVDW_Backbone:
    IntArrayDestroy(&pThis->checkVDW_backbone_smallmolAtomHasXyz);
    break;
  case Type_PlacingAction_CheckVDW_Internal:
    for (int i = 0;i < pThis->checkVDW_internal_smallMolAtomNum;i++)
    {
      free(pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i]);
    }
    free(pThis->checkVDW_internal_smallMolAtom13bondedMatrix);
    pThis->checkVDW_internal_smallMolAtom13bondedMatrix = NULL;
    pThis->checkVDW_internal_smallMolAtomNum = 0;
    IntArrayDestroy(&pThis->checkVDW_internal_smallmolAtomHasXyz);
    break;
  case Type_PlacingAction_CheckRMSD:
    IntArrayDestroy(&pThis->checkRMSD_smallmolAtomHasXyz);
    break;
  case Type_PlacingAction_CheckMultiCons:
    pThis->checkMultiCons_cataConsNum = 0;
    StringArrayDestroy(&pThis->checkMultiCons_firstSiteChainNames);
    StringArrayDestroy(&pThis->checkMultiCons_firstSiteResidueNames);
    IntArrayDestroy(&pThis->checkMultiCons_firstSitePosInChains);
    StringArrayDestroy(&pThis->checkMultiCons_secondSiteChainNames);
    StringArrayDestroy(&pThis->checkMultiCons_secondSiteResidueNames);
    IntArrayDestroy(&pThis->checkMultiCons_secondSitePosInChains);
    for (int i = 0;i < MAX_NUM_CHECK_MULTI_CONS;i++)
    {
      if (pThis->checkMultiCons_predeterminedConsRelations[i] != NULL)
      {
        for (int j = 0;j < pThis->checkMultiCons_rotNumOnEachSite[i][0];j++)
        {
          free(pThis->checkMultiCons_predeterminedConsRelations[i][j]);
        }
        free(pThis->checkMultiCons_predeterminedConsRelations[i]);
        pThis->checkMultiCons_predeterminedConsRelations[i] = NULL;
      }
    }
    break;
  default:
    break;
  }
  return Success;
}


BOOL PlacingActionValidParamName(char* paramName)
{
  if (isalpha(paramName[0]) || paramName[0] == '_') return TRUE;
  else
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    sprintf(errMsg, "in file %s line %d, parameter name %s is invalid, it can only begin with a letter or an underscore, e.g. D1, D2 or D3 for distance, A1, A2 or A3 for angles, T1, T2 or T3 for torsions", __FILE__, __LINE__, paramName);
    TraceError(errMsg, ValueError);
    return FALSE;
  }
}


int PlacingActionRead_CALC(PlacingAction* pThis, StringArray* pAtomNames, XYZArray* pAtomXYZs, StringArray* pParamNames, DoubleArray* pParams, StringArray* pContent)
{
  int index;
  char* atomName;
  char* paramName;
  char errMsg[MAX_LEN_ERR_MSG+1];
  int result;

  if (StringArrayGetNum(pContent) != 9)
  {
    return FormatError;
  }
  for (int i = 0;i < 3;i++)
  {
    atomName = StringArrayGet(pContent, i);
    // skip the '*' which denotes an improper dihedral
    if (i == 2 && atomName[0] == '*') atomName++;
    result = StringArrayFind(pAtomNames, atomName, &index);
    if (FAILED(result))
    {
      result = ValueError;
      sprintf(errMsg, "in file %s line %d, atom %s not exist or calculated, please use $%s rather than %s for ligand atom", __FILE__, __LINE__, atomName, atomName, atomName);
      TraceError(errMsg, result);
      return result;
    }
    IntArrayAppend(&pThis->calc_atoms, index);
  }
  atomName = StringArrayGet(pContent, 3);
  result = StringArrayFind(pAtomNames, atomName, &index);
  if (FAILED(result))
  {
    StringArrayAppend(pAtomNames, atomName);
    index = StringArrayGetNum(pAtomNames) - 1;
    XYZArrayResize(pAtomXYZs, XYZArrayGetNum(pAtomXYZs) + 1);
  }
  IntArrayAppend(&pThis->calc_atoms, index);

  for (int i = 0;i < 5;i++)
  {
    paramName = StringArrayGet(pContent, 4 + i);
    if (paramName[0] == '+' || paramName[0] == '-' || isdigit(paramName[0]))
    {
      StringArrayAppend(pParamNames, "");
      if (1 <= i && i <= 3)
      {
        DoubleArrayAppend(pParams, DegToRad(atof(paramName)));
      }
      else
      {
        DoubleArrayAppend(pParams, atof(paramName));
      }
      index = StringArrayGetNum(pParamNames) - 1;
      IntArrayAppend(&pThis->calc_params, index);
    }
    else
    {
      result = StringArrayFind(pParamNames, paramName, &index);
      if (FAILED(result))
      {
        result = ValueError;
        sprintf(errMsg, "in file %s line %d, parameter '%s' was not defined", __FILE__, __LINE__, paramName);
        TraceError(errMsg, result);
        return result;
      }
      IntArrayAppend(&pThis->calc_params, index);
    }
  }
  return Success;
}


int PlacingActionRead_CHECK_CATA_CONS(PlacingAction* pThis, StringArray* pContent)
{
  if (StringArrayGetNum(pContent) != 6) return FormatError;

  strcpy(pThis->checkCataCons_firstSiteChainName, StringArrayGet(pContent, 0));
  strcpy(pThis->checkCataCons_firstSiteResidueName, StringArrayGet(pContent, 1));
  pThis->checkCataCons_firstSitePosInChain = atoi(StringArrayGet(pContent, 2));

  strcpy(pThis->checkCataCons_secondSiteChainName, StringArrayGet(pContent, 3));
  strcpy(pThis->checkCataCons_secondSiteResidueName, StringArrayGet(pContent, 4));
  pThis->checkCataCons_secondSitePosInChain = atoi(StringArrayGet(pContent, 5));

  return Success;
}


int PlacingActionRead_CHECK_MULTI_CONS(PlacingAction* pThis, StringArray* pContent)
{
  pThis->checkMultiCons_cataConsNum = StringArrayGetNum(pContent) / 6;
  for (int i = 0;i < pThis->checkMultiCons_cataConsNum;i++)
  {
    StringArrayAppend(&pThis->checkMultiCons_firstSiteChainNames, StringArrayGet(pContent, 6 * i));
    StringArrayAppend(&pThis->checkMultiCons_firstSiteResidueNames, StringArrayGet(pContent, 6 * i + 1));
    IntArrayAppend(&pThis->checkMultiCons_firstSitePosInChains, atoi(StringArrayGet(pContent, 6 * i + 2)));
    StringArrayAppend(&pThis->checkMultiCons_secondSiteChainNames, StringArrayGet(pContent, 6 * i + 3));
    StringArrayAppend(&pThis->checkMultiCons_secondSiteResidueNames, StringArrayGet(pContent, 6 * i + 4));
    IntArrayAppend(&pThis->checkMultiCons_secondSitePosInChains, atoi(StringArrayGet(pContent, 6 * i + 5)));
  }
  return Success;
}


int PlacingActionRead_CHECK_RMSD(PlacingAction* pThis, StringArray* pContent)
{
  char* withHydrogen;
  char errMsg[MAX_LEN_ERR_MSG+1];
  int result;

  if (StringArrayGetNum(pContent) != 2) return FormatError;

  withHydrogen = StringArrayGet(pContent, 0);
  if (strcmp(withHydrogen, "WITH_HYDROGEN") == 0)
  {
    pThis->checkRMSD_withHydrogen = TRUE;
  }
  else if (strcmp(withHydrogen, "WITHOUT_HYDROGEN") == 0)
  {
    pThis->checkRMSD_withHydrogen = FALSE;
  }
  else
  {
    result = FormatError;
    sprintf(errMsg, "in file %s line %d, unrecognized keyword %s", __FILE__, __LINE__, withHydrogen);
    TraceError(errMsg, result);
    return result;
  }

  pThis->checkRMSD_minDifference = atof(StringArrayGet(pContent, 1));
  return Success;
}


int    PlacingActionRead_CHECK_VDW_BACKBONE(PlacingAction* pThis, StringArray* pContent)
{
  char* withHydrogen;
  char errMsg[MAX_LEN_ERR_MSG+1];
  int result;

  if (StringArrayGetNum(pContent) != 2) return FormatError;

  withHydrogen = StringArrayGet(pContent, 0);
  if (strcmp(withHydrogen, "WITH_HYDROGEN") == 0)
  {
    pThis->checkVDW_backbone_withHydrogen = TRUE;
  }
  else if (strcmp(withHydrogen, "WITHOUT_HYDROGEN") == 0)
  {
    pThis->checkVDW_backbone_withHydrogen = FALSE;
  }
  else
  {
    result = FormatError;
    sprintf(errMsg, "in file %s line %d, unrecognized keyword %s", __FILE__, __LINE__, withHydrogen);
    TraceError(errMsg, result);
    return result;
  }

  pThis->checkVDW_backbone_maxAllowed = atof(StringArrayGet(pContent, 1));

  return Success;
}


int    PlacingActionRead_CHECK_VDW_INTERNAL(PlacingAction* pThis, StringArray* pContent)
{
  char* withHydrogen;
  char errMsg[MAX_LEN_ERR_MSG+1];
  int result;

  if (StringArrayGetNum(pContent) != 2) return FormatError;

  withHydrogen = StringArrayGet(pContent, 0);
  if (strcmp(withHydrogen, "WITH_HYDROGEN") == 0)
  {
    pThis->checkVDW_internal_withHydrogen = TRUE;
  }
  else if (strcmp(withHydrogen, "WITHOUT_HYDROGEN") == 0)
  {
    pThis->checkVDW_internal_withHydrogen = FALSE;
  }
  else
  {
    result = FormatError;
    sprintf(errMsg, "in file %s line %d, unrecognized keyword %s", __FILE__, __LINE__, withHydrogen);
    TraceError(errMsg, result);
    return result;
  }
  pThis->checkVDW_internal_maxAllowed = atof(StringArrayGet(pContent, 1));

  return Success;
}


int    PlacingActionRead_EVALUATE(PlacingAction* pThis, StringArray* pAtomNames, StringArray* pParamNames, DoubleArray* pParams, StringArray* pContent)
{
  char* paramName;
  char* atomName;
  char* typeName;
  int    result;
  char errMsg[MAX_LEN_ERR_MSG+1];
  int index;

  int wordNum = StringArrayGetNum(pContent);

  if (wordNum < 4)
    return FormatError;

  paramName = StringArrayGet(pContent, 0);
  if (!PlacingActionValidParamName(paramName))
  {
    return FormatError;
  }

  if (FAILED(StringArrayFind(pParamNames, paramName, &index)))
  {
    StringArrayAppend(pParamNames, paramName);
    index = StringArrayGetNum(pParamNames) - 1;
    DoubleArrayAppend(pParams, 0.0);
    pThis->evaluate_param = index;
  }
  else
  {
    result = FormatError;
    sprintf(errMsg, "in file %s line %d, parameter %s already defined", __FILE__, __LINE__, paramName);
    TraceError(errMsg, result);
    return result;
  }

  typeName = StringArrayGet(pContent, 1);

  int atomNum = 0;
  if (strcmp(typeName, "HALF_DISTANCE") == 0)
  {
    pThis->evaluate_type = Type_PlacingActionEvaluate_HalfDistance;
    atomNum = 2;
  }
  else if (strcmp(typeName, "DISTANCE") == 0)
  {
    pThis->evaluate_type = Type_PlacingActionEvaluate_Distance;
    atomNum = 2;
  }
  else if (strcmp(typeName, "ANGLE") == 0)
  {
    pThis->evaluate_type = Type_PlacingActionEvaluate_Angle;
    atomNum = 3;
  }
  else if (strcmp(typeName, "TORSION") == 0)
  {
    pThis->evaluate_type = Type_PlacingActionEvaluate_Torsion;
    atomNum = 4;
  }
  else
  {
    result = FormatError;
    sprintf(errMsg, "in file %s line %d, keyword '%s' not defined", __FILE__, __LINE__, typeName);
    TraceError(errMsg, result);
    return result;
  }

  if (atomNum + 2 != wordNum)
  {
    return FormatError;
  }

  for (int i = 0;i < atomNum;i++)
  {
    atomName = StringArrayGet(pContent, i + 2);
    if (FAILED(StringArrayFind(pAtomNames, atomName, &index)))
    {
      sprintf(errMsg, "in file %s line %d, atom %s not defined or calculated yet", __FILE__, __LINE__, atomName);
      result = FormatError;
      TraceError(errMsg, result);
      return result;
    }
    IntArrayAppend(&pThis->evaluate_atoms, index);
  }

  return Success;
}


int    PlacingActionRead_LOAD(PlacingAction* pThis, StringArray* pAtomNames, XYZArray* pAtomXYZs, StringArray* pContent)
{
  int index;
  char* newAtomName;
  for (int i = 0; i < StringArrayGetNum(pContent); i++)
  {
    newAtomName = StringArrayGet(pContent, i);
    if (strlen(newAtomName) > MAX_LEN_ATOM_NAME)
    {
      int result = NameError;
      char errMsg[MAX_LEN_ERR_MSG+1];
      sprintf(errMsg, "in file %s line %d, atom name %s is too long", __FILE__, __LINE__, newAtomName);
      TraceError(errMsg, result);
      return result;
    }
    if (FAILED(StringArrayFind(pAtomNames, newAtomName, &index)))
    {
      StringArrayAppend(pAtomNames, newAtomName);
      index = StringArrayGetNum(pAtomNames) - 1;
      IntArrayAppend(&pThis->load_atoms, index);
      XYZArrayResize(pAtomXYZs, XYZArrayGetNum(pAtomXYZs) + 1);
    }
  }
  return Success;
}


int PlacingActionRead_VARIATE(PlacingAction* pThis, StringArray* pParamNames, DoubleArray* pParams, StringArray* pContent)
{
  char* name;
  char* degreeOrAnsgrom;
  int index;
  char errMsg[MAX_LEN_ERR_MSG+1];
  int    result;

  if (StringArrayGetNum(pContent) != 5) return FormatError;
  name = StringArrayGet(pContent, 0);
  if (!PlacingActionValidParamName(name)) return FormatError;

  degreeOrAnsgrom = StringArrayGet(pContent, 1);
  if (FAILED(StringArrayFind(pParamNames, name, &index)))
  {
    StringArrayAppend(pParamNames, name);
    DoubleArrayAppend(pParams, 0.0);

    index = StringArrayGetNum(pParamNames) - 1;
    pThis->variate_param = index;

    pThis->variate_from = atof(StringArrayGet(pContent, 2));
    pThis->variate_to = atof(StringArrayGet(pContent, 3));
    pThis->variate_increment = atof(StringArrayGet(pContent, 4));

    pThis->variate_to += DISTURBANCE_IN_RANGE_CHECK;
    if (strcmp(degreeOrAnsgrom, "DEGREE") == 0)
    {
      pThis->variate_from = DegToRad(pThis->variate_from);
      pThis->variate_to = DegToRad(pThis->variate_to);
      pThis->variate_increment = DegToRad(pThis->variate_increment);
      if (pThis->variate_to > pThis->variate_from + 2 * PI - DISTURBANCE_IN_RANGE_CHECK)
        pThis->variate_to = pThis->variate_from + 2 * PI - DISTURBANCE_IN_RANGE_CHECK;
    }
    else if (strcmp(degreeOrAnsgrom, "ANGSTROM") == 0)
    {
      ;
    }
    else
    {
      result = FormatError;
      sprintf(errMsg, "in file %s line %d, invalid unit '%s'; only DEGREE or ANGSTROM is allowed", __FILE__, __LINE__, degreeOrAnsgrom);
      TraceError(errMsg, result);
      return result;
    }
  }
  else
  {
    result = NameError;
    sprintf(errMsg, "in file %s line %d, parameter %s has already been defined", __FILE__, __LINE__, name);
    TraceError(errMsg, result);
    return result;
  }
  return Success;
}


int    PlacingRuleReadFile(PlacingRule* pThis, FileReader* pFileReader)
{
  BOOL done;
  char line[MAX_LEN_ONE_LINE_CONTENT+1];
  char errMsg[MAX_LEN_ERR_MSG+1];
  int    result;

  done = FALSE;
  while (!FAILED(FileReaderGetNextLine(pFileReader, line)))
  {
    char keyword[MAX_LEN_ONE_LINE_CONTENT+1];
    char buffer[MAX_LEN_ONE_LINE_CONTENT+1];
    PlacingAction* pCurAction;
    StringArray wordsInLine;
    StringArrayCreate(&wordsInLine);
    ExtractFirstStringFromSourceString(keyword, line);
    if (strcmp(keyword, "PLACING") == 0 && done == FALSE)
    {
      ExtractFirstStringFromSourceString(pThis->chainName, line);
      ExtractFirstStringFromSourceString(pThis->residueName, line);
      ExtractFirstStringFromSourceString(buffer, line);
      pThis->posInChain = atoi(buffer);
      ExtractFirstStringFromSourceString(pThis->rotamerType, line);
      done = TRUE;
      StringArrayDestroy(&wordsInLine);
      continue;
    }
    else if (strcmp(keyword, "PLACING") == 0 && done)
    {
      FileReaderSetCurrentPos(pFileReader, FileReaderGetCurrentPos(pFileReader) - 1);
      StringArrayDestroy(&wordsInLine);
      return Success;
    }


    (pThis->actionNum)++;
    pThis->actions = (PlacingAction*)realloc(pThis->actions,
      sizeof(PlacingAction) * pThis->actionNum);
    pCurAction = &pThis->actions[pThis->actionNum - 1];

    StringArraySplitString(&wordsInLine, line, ' ');

    if (strcmp(keyword, "LOAD") == 0)
    {
      PlacingActionCreate(pCurAction, Type_PlacingAction_Load);
      result = PlacingActionRead_LOAD(pCurAction, &pThis->atomNames, &pThis->atomXYZs, &wordsInLine);
    }
    else if (strcmp(keyword, "VARIATE") == 0)
    {
      PlacingActionCreate(pCurAction, Type_PlacingAction_Variate);
      result = PlacingActionRead_VARIATE(pCurAction, &pThis->paramNames, &pThis->params, &wordsInLine);
    }
    else if (strcmp(keyword, "CALC") == 0)
    {
      PlacingActionCreate(pCurAction, Type_PlacingAction_Calc);
      result = PlacingActionRead_CALC(pCurAction, &pThis->atomNames, &pThis->atomXYZs,
        &pThis->paramNames, &pThis->params, &wordsInLine);
    }
    else if (strcmp(keyword, "EVALUATE") == 0)
    {
      PlacingActionCreate(pCurAction, Type_PlacingAction_Evaluate);
      result = PlacingActionRead_EVALUATE(pCurAction, &pThis->atomNames, &pThis->paramNames,
        &pThis->params, &wordsInLine);
    }
    else if (strcmp(keyword, "CHECK_CATA_CONS") == 0)
    {
      PlacingActionCreate(pCurAction, Type_PlacingAction_CheckCataCons);
      result = PlacingActionRead_CHECK_CATA_CONS(pCurAction, &wordsInLine);
    }
    else if (strcmp(keyword, "CHECK_MULTI_CONS_BEGIN") == 0)
    {
      StringArray allContent;
      StringArrayCreate(&allContent);
      PlacingActionCreate(pCurAction, Type_PlacingAction_CheckMultiCons);
      while (FileReaderGetNextLine(pFileReader, line) == Success)
      {
        StringArraySplitString(&wordsInLine, line, ' ');
        if (StringArrayGetNum(&wordsInLine) != 6)
        {
          break;
        }
        else
        {
          int i;
          for (i = 0;i < 6;i++)
          {
            StringArrayAppend(&allContent, StringArrayGet(&wordsInLine, i));
          }
        }
      }
      if (strcmp(line, "CHECK_MULTI_CONS_END") == 0)
      {
        int num = StringArrayGetNum(&allContent) / 6;
        if (num > MAX_NUM_CHECK_MULTI_CONS)
        {
          sprintf(errMsg, "in file %s line %d, no more than %d constraints can be checked simultaneously", __FILE__, __LINE__, MAX_NUM_CHECK_MULTI_CONS);
          result = IndexError;
          TraceError(errMsg, result);
          return result;
        }
        result = PlacingActionRead_CHECK_MULTI_CONS(pCurAction, &allContent);
      }
      else
      {
        result = FormatError;
      }
      StringArrayDestroy(&allContent);
    }
    else if (strcmp(keyword, "CHECK_VDW_BACKBONE") == 0)
    {
      PlacingActionCreate(pCurAction, Type_PlacingAction_CheckVDW_Backbone);
      result = PlacingActionRead_CHECK_VDW_BACKBONE(pCurAction, &wordsInLine);
    }
    else if (strcmp(keyword, "CHECK_VDW_INTERNAL") == 0)
    {
      PlacingActionCreate(pCurAction, Type_PlacingAction_CheckVDW_Internal);
      result = PlacingActionRead_CHECK_VDW_INTERNAL(pCurAction, &wordsInLine);
    }
    else if (strcmp(keyword, "CHECK_RMSD") == 0)
    {
      PlacingActionCreate(pCurAction, Type_PlacingAction_CheckRMSD);
      result = PlacingActionRead_CHECK_RMSD(pCurAction, &wordsInLine);
    }
    else
    {
      sprintf(errMsg, "in file %s line %d, unrecognized keyword %s", __FILE__, __LINE__, keyword);
      result = FormatError;
      TraceError(errMsg, result);
      return result;
    }

    if (FAILED(result))
    {
      sprintf(errMsg, "in file %s line %d, when reading : \n%s  %s", __FILE__, __LINE__, keyword, line);
      TraceError(errMsg, result);
      StringArrayDestroy(&wordsInLine);
      return result;
    }
    StringArrayDestroy(&wordsInLine);
  }

  return Success;
}

int PlacingRuleCreate(PlacingRule* pThis, char* file)
{
  FileReader fr;
  char errMsg[MAX_LEN_ERR_MSG+1];
  int result = FileReaderCreate(&fr, file);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, cannot read file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, result);
    return result;
  }

  pThis->vdwInternal = pThis->vdwBackbone = 0.0;

  StringArrayCreate(&pThis->atomNames);
  XYZArrayCreate(&pThis->atomXYZs, 0);

  IntArrayCreate(&pThis->atomPosOnSmallMol, 0);
  XYZArrayCreate(&pThis->smallMolAtomXYZs, 0);

  pThis->xyzValidArray = NULL;

  StringArrayCreate(&pThis->paramNames);
  DoubleArrayCreate(&pThis->params, 0);

  pThis->actions = NULL;
  pThis->actionNum = 0;

  pThis->deployedFlag = FALSE;

  result = PlacingRuleReadFile(pThis, &fr);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, failed to read ligand placement file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, result);
    return result;
  }

  FileReaderDestroy(&fr);

  return Success;
}


void PlacingRuleDestroy(PlacingRule* pThis)
{
  int i;
  for (i = 0;i < pThis->actionNum;i++)
  {
    PlacingActionDestroy(&pThis->actions[i]);
  }

  StringArrayDestroy(&pThis->atomNames);
  XYZArrayDestroy(&pThis->atomXYZs);

  IntArrayDestroy(&pThis->atomPosOnSmallMol);
  XYZArrayDestroy(&pThis->smallMolAtomXYZs);

  if (pThis->xyzValidArray != NULL)
  {
    free(pThis->xyzValidArray);
    pThis->xyzValidArray = NULL;
  }

  StringArrayDestroy(&pThis->paramNames);
  DoubleArrayDestroy(&pThis->params);

  free(pThis->actions);
  pThis->actions = NULL;
  pThis->actionNum = 0;
}


BOOL PlacingRuleGetDeployedFlag(PlacingRule* pThis)
{
  return pThis->deployedFlag;
}


char* PlacingRuleGetResiName(PlacingRule* pThis)
{
  return pThis->residueName;
}


int PlacingRuleGetPosInChain(PlacingRule* pThis)
{
  return pThis->posInChain;
}


char* PlacingRuleGetChainName(PlacingRule* pThis)
{
  return pThis->chainName;
}


char* PlacingRuleGetRotamerType(PlacingRule* pThis)
{
  return pThis->rotamerType;
}


//double PlacingRuleGetTruncatedBackboneRange(PlacingRule* pThis)
//{
//  for (int i = 0;i < pThis->actionNum;i++)
//  {
//    if (pThis->actions[i].actionType == Type_PlacingAction_CheckVDW_Backbone)
//    {
//      return pThis->actions[i].checkVDW_backbone_activeRange;
//    }
//  }
//  return 0.0;
//}


int PlacingActionDeploy_CheckCataCons(PlacingAction* pThis,
  Residue* pSmallMol,
  CataConsSitePairArray* pCataConsArray,
  char* startingSiteChainName,
  int startingSitePosInChain,
  int relatedProteinSiteNum,
  DesignSite** relatedProteinSites)
{
  pThis->checkCataCons_pCataCon = CataConsSitePairArrayFind(
    pCataConsArray,
    pThis->checkCataCons_firstSiteChainName,
    pThis->checkCataCons_firstSitePosInChain,
    pThis->checkCataCons_firstSiteResidueName,
    pThis->checkCataCons_secondSiteChainName,
    pThis->checkCataCons_secondSitePosInChain,
    pThis->checkCataCons_secondSiteResidueName);

  if (pThis->checkCataCons_pCataCon == NULL)
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    sprintf(errMsg, "in file %s line %d, cannot find the catalytic constraint between site %s%d%s and site %s%d%s", __FILE__, __LINE__,
      pThis->checkCataCons_firstSiteChainName,
      pThis->checkCataCons_firstSitePosInChain,
      pThis->checkCataCons_firstSiteResidueName,
      pThis->checkCataCons_secondSiteChainName,
      pThis->checkCataCons_secondSitePosInChain,
      pThis->checkCataCons_secondSiteResidueName);
    TraceError(errMsg, DataNotExistError);
    return DataNotExistError;
  }

  //Determine the identity of the first and second site of the catalytic constraint
  for (int i = 0;i < 2;i++)
  {
    int posInChain;
    char* chainName;
    if (i == 0)
    {
      posInChain = pThis->checkCataCons_pCataCon->pos1;
      chainName = pThis->checkCataCons_pCataCon->chnName1;
    }
    else
    {
      posInChain = pThis->checkCataCons_pCataCon->pos2;
      chainName = pThis->checkCataCons_pCataCon->chnName2;
    }

    pThis->checkCataCons_pSite[i] = -1;

    if (posInChain == ResidueGetPosInChain(pSmallMol) &&
      strcmp(chainName, ResidueGetChainName(pSmallMol)) == 0)
    {
      pThis->checkCataCons_pSite[i] = relatedProteinSiteNum;  //the small molecule
    }
    else if (
      posInChain == startingSitePosInChain &&
      strcmp(chainName, startingSiteChainName) == 0)
    {
      pThis->checkCataCons_pSite[i] = relatedProteinSiteNum + 1;  //the starting site
    }
    else
    {
      for (int j = 0;j < relatedProteinSiteNum;j++)
      {
        //this is probably a smallmol design site, not finished, skip it
        if (relatedProteinSites[j] == NULL) continue;
        if (posInChain == DesignSiteGetPosInChain(relatedProteinSites[j]) &&
          strcmp(chainName, DesignSiteGetChainName(relatedProteinSites[j])) == 0)
        {
          pThis->checkCataCons_pSite[i] = j;
          break;
        }
      }
    }

    if (pThis->checkCataCons_pSite[i] == -1)
    {
      char errMsg[MAX_LEN_ERR_MSG+1];
      sprintf(errMsg, "in file %s line %d, cannot find design site %s%d needed by checking of catalytic constraint", __FILE__, __LINE__, chainName, posInChain);
      TraceError(errMsg, DataNotExistError);
      return DataNotExistError;
    }
  }

  return Success;
}


int PlacingActionDeploy_CheckMultiCons(PlacingAction* pThis,
  Residue* pSmallMol,
  CataConsSitePairArray* pCataConsArray,
  char* startingSiteChainName,
  char* startingSiteResidueName,
  int startingSitePosInChain,
  int relatedProteinSiteNum,
  DesignSite** relatedProteinSites)
{
  //Find Identity of every CataConsPair
  for (int i = 0;i < pThis->checkMultiCons_cataConsNum;i++)
  {
    pThis->checkMultiCons_pCataCons[i] = CataConsSitePairArrayFind(
      pCataConsArray,
      StringArrayGet(&pThis->checkMultiCons_firstSiteChainNames, i),
      IntArrayGet(&pThis->checkMultiCons_firstSitePosInChains, i),
      StringArrayGet(&pThis->checkMultiCons_firstSiteResidueNames, i),
      StringArrayGet(&pThis->checkMultiCons_secondSiteChainNames, i),
      IntArrayGet(&pThis->checkMultiCons_secondSitePosInChains, i),
      StringArrayGet(&pThis->checkMultiCons_secondSiteResidueNames, i));

    if (pThis->checkMultiCons_pCataCons[i] == NULL)
    {
      char errMsg[MAX_LEN_ERR_MSG+1];
      sprintf(errMsg, "in file %s line %d, cannot find the catalytic constraint between site %s%d%s and site%s%d%s", __FILE__, __LINE__,
        StringArrayGet(&pThis->checkMultiCons_firstSiteChainNames, i),
        IntArrayGet(&pThis->checkMultiCons_firstSitePosInChains, i),
        StringArrayGet(&pThis->checkMultiCons_firstSiteResidueNames, i),
        StringArrayGet(&pThis->checkMultiCons_secondSiteChainNames, i),
        IntArrayGet(&pThis->checkMultiCons_secondSitePosInChains, i),
        StringArrayGet(&pThis->checkMultiCons_secondSiteResidueNames, i));
      TraceError(errMsg, DataNotExistError);
      return DataNotExistError;
    }
  }
  //Determine the identity of the first and second site of every catalytic constraint
  for (int i = 0;i < pThis->checkMultiCons_cataConsNum;i++)
  {
    for (int j = 0;j < 2;j++)
    {
      int posInChain;
      char* chainName;
      if (j == 0)
      {
        posInChain = pThis->checkMultiCons_pCataCons[i]->pos1;
        chainName = pThis->checkMultiCons_pCataCons[i]->chnName1;
      }
      else
      {
        posInChain = pThis->checkMultiCons_pCataCons[i]->pos2;
        chainName = pThis->checkMultiCons_pCataCons[i]->chnName2;
      }

      pThis->checkMultiCons_siteIndexes[i][j] = -1;

      if (posInChain == ResidueGetPosInChain(pSmallMol) &&
        strcmp(chainName, ResidueGetChainName(pSmallMol)) == 0)
      {
        pThis->checkMultiCons_siteIndexes[i][j] = relatedProteinSiteNum;
      }
      else if (
        posInChain == startingSitePosInChain &&
        strcmp(chainName, startingSiteChainName) == 0)
      {
        pThis->checkMultiCons_siteIndexes[i][j] = relatedProteinSiteNum + 1;
      }
      else
      {
        int s;
        for (s = 0;s < relatedProteinSiteNum;s++)
        {
          //this is probably a smallmol design site, not finished, skip it
          if (relatedProteinSites[s] == NULL) continue;
          if (posInChain == DesignSiteGetPosInChain(relatedProteinSites[s]) &&
            strcmp(chainName, DesignSiteGetChainName(relatedProteinSites[s])) == 0)
          {
            pThis->checkMultiCons_siteIndexes[i][j] = s;
            break;
          }
        }
      }

      if (pThis->checkMultiCons_siteIndexes[i][j] == -1)
      {
        char errMsg[MAX_LEN_ERR_MSG+1];
        sprintf(errMsg, "in file %s line %d, cannot find design site %s%d needed by checking of catalytic constraint", __FILE__, __LINE__, chainName, posInChain);
        TraceError(errMsg, DataNotExistError);
        return DataNotExistError;
      }
    }
  }

  //Determine whether the catalytic constraints between only protein sites are satisfied
  for (int i = 0;i < pThis->checkMultiCons_cataConsNum;i++)
  {
    int indexOfSite1 = pThis->checkMultiCons_siteIndexes[i][0];
    int indexOfSite2 = pThis->checkMultiCons_siteIndexes[i][1];
    if (indexOfSite1 < relatedProteinSiteNum && indexOfSite2 < relatedProteinSiteNum)
    {
      int rotOnSite1;
      int rotOnSite2;
      int rotNumOnSite1;
      int rotNumOnSite2;
      pThis->checkMultiCons_rotNumOnEachSite[i][0] = rotNumOnSite1 =
        RotamersGetNum(DesignSiteGetRotamers(relatedProteinSites[indexOfSite1]));
      pThis->checkMultiCons_rotNumOnEachSite[i][1] = rotNumOnSite2 =
        RotamersGetNum(DesignSiteGetRotamers(relatedProteinSites[indexOfSite2]));
      pThis->checkMultiCons_predeterminedConsRelations[i] = (BOOL**)malloc(sizeof(BOOL*) * rotNumOnSite1);
      for (rotOnSite1 = 0;rotOnSite1 < rotNumOnSite1;rotOnSite1++)
      {
        pThis->checkMultiCons_predeterminedConsRelations[i][rotOnSite1] =
          (BOOL*)malloc(sizeof(BOOL) * rotNumOnSite2);
        for (rotOnSite2 = 0; rotOnSite2 < rotNumOnSite2; rotOnSite2++)
        {
          BOOL checkResult;
          checkResult = CataConsSitePairCheck(pThis->checkMultiCons_pCataCons[i],
            RotamersGet(DesignSiteGetRotamers(relatedProteinSites[indexOfSite1]), rotOnSite1),
            RotamersGet(DesignSiteGetRotamers(relatedProteinSites[indexOfSite2]), rotOnSite2)
          );
          pThis->checkMultiCons_predeterminedConsRelations[i][rotOnSite1][rotOnSite2] = checkResult;
        }
      }
    }
  }

  //Construct the data structure needed by the recurrence algorithm in PlacingRuleProcess_CHECK_MULTI_CONS()
  pThis->checkMultiCons_stepNum = 0;
  for (int i = 0;i < pThis->checkMultiCons_cataConsNum;i++)
  {
    for (int j = 0;j < 2;j++)
    {
      int siteToFindInPreviousSteps = pThis->checkMultiCons_siteIndexes[i][j];
      BOOL addNewStep = TRUE;
      for (int k = 0;k < pThis->checkMultiCons_stepNum;k++)
      {
        if (pThis->checkMultiCons_steps[k].designSite == siteToFindInPreviousSteps)
        {
          addNewStep = FALSE;
          break;
        }
      }
      if (addNewStep)
      {
        pThis->checkMultiCons_steps[pThis->checkMultiCons_stepNum].designSite = siteToFindInPreviousSteps;
        pThis->checkMultiCons_steps[pThis->checkMultiCons_stepNum].nConsToCheckAtThisStep = 0;
        pThis->checkMultiCons_stepNum++;
      }
    }
  }

  for (int i = 0;i < pThis->checkMultiCons_cataConsNum;i++)
  {
    BOOL site1Found = FALSE;
    BOOL site2Found = FALSE;
    int stepIndex = 0;
    CheckMultiConsStep* pCurStep = NULL;
    for (stepIndex = 0; stepIndex < pThis->checkMultiCons_stepNum; stepIndex++)
    {
      pCurStep = &pThis->checkMultiCons_steps[stepIndex];
      if (pCurStep->designSite == pThis->checkMultiCons_siteIndexes[i][0])
      {
        site1Found = TRUE;
      }
      if (pCurStep->designSite == pThis->checkMultiCons_siteIndexes[i][1])
      {
        site2Found = TRUE;
      }
      if (site1Found && site2Found)
      {
        break;
      }
    }
    pCurStep->consToCheckAtThisStep[pCurStep->nConsToCheckAtThisStep] = i;
    pCurStep->nConsToCheckAtThisStep++;
  }

  return Success;
}


int PlacingActionDeploy_CheckRMSD(PlacingAction* pThis, IntArray* pSmallmolAtomGetXyzByThisStep)
{
  IntArrayCopy(&pThis->checkRMSD_smallmolAtomHasXyz, pSmallmolAtomGetXyzByThisStep);
  return Success;
}


int PlacingActionDeploy_CheckVDWBackbone(PlacingAction* pThis, IntArray* pSmallmolAtomGetXyzByThisStep)
{
  IntArrayCopy(&pThis->checkVDW_backbone_smallmolAtomHasXyz, pSmallmolAtomGetXyzByThisStep);
  return Success;
}


int PlacingActionDeploy_CheckVDWInternal(PlacingAction* pThis, Residue* pSmallMol, IntArray* pSmallmolAtomGetXyzByThisStep)
{
  pThis->checkVDW_internal_smallMolAtomNum = ResidueGetAtomNum(pSmallMol);
  pThis->checkVDW_internal_smallMolAtom13bondedMatrix =
    (BOOL**)malloc(sizeof(BOOL*) * ResidueGetAtomNum(pSmallMol));
  for (int i = 0;i < ResidueGetAtomNum(pSmallMol);i++)
  {
    pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i] =
      (BOOL*)malloc(sizeof(BOOL) * ResidueGetAtomNum(pSmallMol));
    memset(pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i], FALSE,
      sizeof(BOOL) * ResidueGetAtomNum(pSmallMol));
  }
  //Find 1-2 bond relation
  for (int i = 0;i < ResidueGetAtomNum(pSmallMol);i++)
  {
    for (int j = 0;j < ResidueGetAtomNum(pSmallMol);j++)
    {
      Atom* pAtomI = ResidueGetAtom(pSmallMol, i);
      Atom* pAtomJ = ResidueGetAtom(pSmallMol, j);
      if (BondsFind(&pSmallMol->bonds, AtomGetName(pAtomI), AtomGetName(pAtomJ)) == Type_Bond_None)
      {
        pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i][j] = FALSE;
      }
      else
      {
        pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i][j] = TRUE;
      }
    }
  }
  //Find 1-3 bond relation
  for (int i = 0;i < ResidueGetAtomNum(pSmallMol);i++)
  {
    for (int j = 0;j < ResidueGetAtomNum(pSmallMol);j++)
    {
      for (int k = 0;k < ResidueGetAtomNum(pSmallMol);k++)
      {
        Atom* pAtomI = ResidueGetAtom(pSmallMol, i);
        Atom* pAtomJ = ResidueGetAtom(pSmallMol, j);
        Atom* pAtomK = ResidueGetAtom(pSmallMol, k);
        if (
          (BondsFind(&pSmallMol->bonds, AtomGetName(pAtomI), AtomGetName(pAtomK)) != Type_Bond_None) &&
          (BondsFind(&pSmallMol->bonds, AtomGetName(pAtomK), AtomGetName(pAtomJ)) != Type_Bond_None)
          )
        {
          pThis->checkVDW_internal_smallMolAtom13bondedMatrix[i][j] = TRUE;
        }
      }
    }
  }

  IntArrayCopy(&pThis->checkVDW_internal_smallmolAtomHasXyz, pSmallmolAtomGetXyzByThisStep);
  return Success;
}


int PlacingActionDeploy_Calc(PlacingAction* pThis,
  IntArray* pAtomPosOnSmallMol,
  IntArray* pSmallmolAtomGetXyzByThisStep)
{
  int atomIndexInPlacingRule = IntArrayGet(&pThis->calc_atoms, 3);
  int atomIndexOnSmallmol = IntArrayGet(pAtomPosOnSmallMol, atomIndexInPlacingRule);
  if (atomIndexOnSmallmol != -1)
  {
    IntArraySet(pSmallmolAtomGetXyzByThisStep, atomIndexOnSmallmol, 1);
  }
  return Success;
}


int PlacingActionDeploy_Load(PlacingAction* pThis,
  IntArray* pAtomPosOnSmallMol,
  IntArray* pSmallmolAtomGetXyzByThisStep)
{
  int i;
  for (i = 0;i < IntArrayGetLength(&pThis->load_atoms);i++)
  {
    int atomIndexInPlacingRule = IntArrayGet(&pThis->load_atoms, i);
    int atomIndexOnSmallmol = IntArrayGet(pAtomPosOnSmallMol, atomIndexInPlacingRule);
    if (atomIndexOnSmallmol != -1)
    {
      IntArraySet(pSmallmolAtomGetXyzByThisStep, atomIndexOnSmallmol, 1);
    }
  }
  return Success;
}


int PlacingRuleDeploy(PlacingRule* pThis, Residue* pSmallMol, CataConsSitePairArray* pCataConsArray, int relatedProteinSiteNum, DesignSite** relatedProteinSites, AtomArray* pTruncatedBackbone)
{
  int result;
  char errMsg[MAX_LEN_ERR_MSG+1];
  IntArray smallmolAtomGetXyzByThisStep;

  if (PlacingRuleGetDeployedFlag(pThis))
  {
    result = ValueError;
    sprintf(errMsg, "in file %s line %d, placing rule %s %s %d has already been deployed", __FILE__, __LINE__, pThis->chainName, pThis->residueName, pThis->posInChain);
    TraceError(errMsg, result);
    return result;
  }
  else
  {
    pThis->deployedFlag = TRUE;
    pThis->pSmallMol = pSmallMol;
    pThis->pTruncBone = pTruncatedBackbone;
  }

  // maintain the set of small molecule atom xyzs stored in this placing rule
  XYZArrayResize(&pThis->smallMolAtomXYZs, ResidueGetAtomNum(pSmallMol));
  pThis->xyzValidArray = (BOOL*)malloc(sizeof(BOOL) * ResidueGetAtomNum(pSmallMol));
  for (int i = 0; i < ResidueGetAtomNum(pSmallMol); i++)
  {
    pThis->xyzValidArray[i] = FALSE;
  }

  // find the exact locations of each Atom, special operations are needed if it is on the small molecule
  for (int i = 0;i < StringArrayGetNum(&pThis->atomNames);i++)
  {
    int index;
    char atomName[MAX_LEN_ONE_LINE_CONTENT];
    strcpy(atomName, StringArrayGet(&pThis->atomNames, i));
    if (atomName[0] == SMALLMOL_ATOM_SYMBOL)
    {
      int result = ResidueFindAtom(pSmallMol, atomName+1, &index);
      // use atomName+1 to skip the '$' character used to represent small molecular atoms
      if (FAILED(result))
      {
        sprintf(errMsg, "in file %s line %d,cannot find atom %s on ligand", __FILE__, __LINE__, atomName);
        TraceError(errMsg, result);
        return result;
      }
      IntArrayAppend(&pThis->atomPosOnSmallMol, index);
    }
    else
    {
      IntArrayAppend(&pThis->atomPosOnSmallMol, -1);
    }
  }


  IntArrayCreate(&smallmolAtomGetXyzByThisStep, XYZArrayGetNum(&pThis->smallMolAtomXYZs));
  for (int j = 0; j < IntArrayGetLength(&smallmolAtomGetXyzByThisStep); j++)
  {
    IntArraySet(&smallmolAtomGetXyzByThisStep, j, 0);
  }
  for (int i = 0;i < pThis->actionNum;i++)
  {
    PlacingAction* pCurAction = &pThis->actions[i];
    //printf("deploy placing action %d\n", i);
    switch (pCurAction->actionType)
    {
    case Type_PlacingAction_Calc:
      result = PlacingActionDeploy_Calc(pCurAction, &pThis->atomPosOnSmallMol, &smallmolAtomGetXyzByThisStep);
      break;
    case Type_PlacingAction_CheckCataCons:
      result = PlacingActionDeploy_CheckCataCons(pCurAction, pSmallMol, pCataConsArray,
        pThis->chainName, pThis->posInChain, relatedProteinSiteNum, relatedProteinSites);
      break;
    case Type_PlacingAction_CheckMultiCons:
      result = PlacingActionDeploy_CheckMultiCons(pCurAction, pSmallMol, pCataConsArray, pThis->chainName,
        pThis->residueName, pThis->posInChain, relatedProteinSiteNum, relatedProteinSites);
      break;
    case Type_PlacingAction_CheckRMSD:
      result = PlacingActionDeploy_CheckRMSD(pCurAction, &smallmolAtomGetXyzByThisStep);
      break;
    case Type_PlacingAction_CheckVDW_Internal:
      result = PlacingActionDeploy_CheckVDWInternal(pCurAction, pSmallMol, &smallmolAtomGetXyzByThisStep);
      break;
    case Type_PlacingAction_CheckVDW_Backbone:
      result = PlacingActionDeploy_CheckVDWBackbone(pCurAction, &smallmolAtomGetXyzByThisStep);
      break;
    case Type_PlacingAction_Load:
      result = PlacingActionDeploy_Load(pCurAction, &pThis->atomPosOnSmallMol, &smallmolAtomGetXyzByThisStep);
      break;
    default:
      result = Success;
      break;
    }
    if (FAILED(result))
    {
      sprintf(errMsg, "in file %s line %d, error occurred when deploy the following placing action", __FILE__, __LINE__);
      PlacingRuleShowAction(pThis, i);
      TraceError(errMsg, result);
      return result;
    }
  }

  IntArrayDestroy(&smallmolAtomGetXyzByThisStep);

  return Success;
}


int PlacingRuleProcess_CALC(PlacingRule* pThis, PlacingAction* pAction)
{
  XYZ atomXYZs[4];
  double icParam[5];
  int atomIndex;
  for (int i = 0;i < 3;i++)
  {
    atomIndex = IntArrayGet(&pAction->calc_atoms, i);
    atomXYZs[i] = *XYZArrayGet(&pThis->atomXYZs, atomIndex);
  }
  for (int i = 0;i < 5;i++)
  {
    int paramIndex = IntArrayGet(&pAction->calc_params, i);
    icParam[i] = DoubleArrayGet(&pThis->params, paramIndex);
  }

  GetFourthAtom(&atomXYZs[0], &atomXYZs[1], &atomXYZs[2], icParam, &atomXYZs[3]);
  atomIndex = IntArrayGet(&pAction->calc_atoms, 3);
  XYZArraySet(&pThis->atomXYZs, atomIndex, &atomXYZs[3]);

  // if the atom is on small mol, additional operations must be taken to maintain consistency
  int atomPosOnSmallMol = IntArrayGet(&pThis->atomPosOnSmallMol, atomIndex);
  if (atomPosOnSmallMol != -1)
  {
    XYZArraySet(&pThis->smallMolAtomXYZs, atomPosOnSmallMol, &atomXYZs[3]);
    pThis->xyzValidArray[atomPosOnSmallMol] = TRUE;
  }
  return Success;
}


BOOL   PlacingRuleProcess_CHECK_CATA_CONS(PlacingRule* pThis, Rotamer* pStartRot, int siteNum, DesignSite** ppSites, PlacingAction* pAction)
{
  BOOL checkResult = FALSE;
  Rotamers rotSetOnlySmallMol;
  Rotamers rotSetOnlyStartingProteinRotamer;
  Rotamers* pRotSetOnFirstAndSecondSite[2] = { NULL,NULL };
  RotamersCreate(&rotSetOnlySmallMol);
  RotamersCreate(&rotSetOnlyStartingProteinRotamer);

  RotamersAdd(&rotSetOnlyStartingProteinRotamer, pStartRot);
  if (pThis->pSmallMol != NULL)
  {
    Rotamer tempRotamer;
    RotamerCreate(&tempRotamer);
    RotamerSetType(&tempRotamer, ResidueGetName(pThis->pSmallMol));
    RotamerAddAtoms(&tempRotamer, ResidueGetAllAtoms(pThis->pSmallMol));
    XYZArrayCopy(&tempRotamer.xyzs, &pThis->smallMolAtomXYZs);
    RotamersAdd(&rotSetOnlySmallMol, &tempRotamer);
    RotamerDestroy(&tempRotamer);
  }

  // determine the identity of the first and second site of the catalytic constraint
  for (int i = 0;i < 2;i++)
  {
    int siteIndex = pAction->checkCataCons_pSite[i];
    if (siteIndex == siteNum)
    {
      pRotSetOnFirstAndSecondSite[i] = &rotSetOnlySmallMol;
    }
    else if (siteIndex == siteNum + 1)
    {
      pRotSetOnFirstAndSecondSite[i] = &rotSetOnlyStartingProteinRotamer;
    }
    else
    {
      pRotSetOnFirstAndSecondSite[i] = DesignSiteGetRotamers(ppSites[siteIndex]);
    }
  }

  // check the constraints
  checkResult = FALSE;
  for (int i = 0;i < RotamersGetNum(pRotSetOnFirstAndSecondSite[0]);i++)
  {
    if (checkResult)
    {
      break;
    }
    for (int j = 0;j < RotamersGetNum(pRotSetOnFirstAndSecondSite[1]);j++)
    {
      checkResult = CataConsSitePairCheck(pAction->checkCataCons_pCataCon,
        RotamersGet(pRotSetOnFirstAndSecondSite[0], i),
        RotamersGet(pRotSetOnFirstAndSecondSite[1], j));
#ifdef DEBUGGING_SMALLMOL
      if (checkResult == TRUE)
      {
        printf("(%s, %d) and (%s, %d) cons is true\n",
          RotamerGetType(RotamerSetGet(pRotSetOnFirstAndSecondSite[0], i)),
          i,
          RotamerGetType(RotamerSetGet(pRotSetOnFirstAndSecondSite[1], j)),
          j);
        break;
      }
#endif
    }
  }

  RotamersDestroy(&rotSetOnlySmallMol);
  RotamersDestroy(&rotSetOnlyStartingProteinRotamer);
  return checkResult;
}


BOOL   PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence(PlacingAction* pAction,
  int relatedProteinSiteNum,
  DesignSite** relatedProteinSites,
  int* rotamersSelectedOnPreviousSteps,
  Rotamer* pSmallMolRotamer,
  Rotamer* pStartingProteinRotamer,
  int currentStepIndex)
{
  int rotNum;
  CheckMultiConsStep* pCurStep = &pAction->checkMultiCons_steps[currentStepIndex];
  BOOL satisfied = FALSE;

  if (currentStepIndex == pAction->checkMultiCons_stepNum)
  {
    return TRUE;
  }

  if (pCurStep->designSite < relatedProteinSiteNum)
  {
    rotNum = RotamersGetNum(DesignSiteGetRotamers(relatedProteinSites[pCurStep->designSite]));
  }
  else
  {
    rotNum = 1;
  }

  for (int i = 0;i < rotNum;i++)
  {
    BOOL checkResult = TRUE;
    rotamersSelectedOnPreviousSteps[pCurStep->designSite] = i;

    for (int j = 0;j < pCurStep->nConsToCheckAtThisStep;j++)
    {
      int consIndex = pCurStep->consToCheckAtThisStep[j];
      CataConsSitePair* pCurCons = pAction->checkMultiCons_pCataCons[consIndex];
      int site1Index = pAction->checkMultiCons_siteIndexes[consIndex][0];
      int site2Index = pAction->checkMultiCons_siteIndexes[consIndex][1];
      int rotOnSite1 = rotamersSelectedOnPreviousSteps[site1Index];
      int rotOnSite2 = rotamersSelectedOnPreviousSteps[site2Index];
      Rotamer* pRotOnSite1;
      Rotamer* pRotOnSite2;

      if (site1Index < relatedProteinSiteNum && site2Index < relatedProteinSiteNum)
      {
        checkResult = pAction->checkMultiCons_predeterminedConsRelations[consIndex][rotOnSite1][rotOnSite2];
#ifdef DEBUGGING_SMALLMOL
        if (checkResult)
        {
          printf("(%d %d, %d %d) satisfy constraint %d\n", site1Index, rotOnSite1, site2Index, rotOnSite2, consIndex);
        }
#endif
      }
      else
      {
        if (site1Index == relatedProteinSiteNum)
        {
          pRotOnSite1 = pSmallMolRotamer;
        }
        else if (site1Index == relatedProteinSiteNum + 1)
        {
          pRotOnSite1 = pStartingProteinRotamer;
        }
        else
        {
          pRotOnSite1 = RotamersGet(DesignSiteGetRotamers(relatedProteinSites[site1Index]), rotOnSite1);
        }

        if (site2Index == relatedProteinSiteNum)
        {
          pRotOnSite2 = pSmallMolRotamer;
        }
        else if (site2Index == relatedProteinSiteNum + 1)
        {
          pRotOnSite2 = pStartingProteinRotamer;
        }
        else
        {
          pRotOnSite2 = RotamersGet(DesignSiteGetRotamers(relatedProteinSites[site2Index]), rotOnSite2);
        }

        checkResult = CataConsSitePairCheck(pCurCons, pRotOnSite1, pRotOnSite2);
#ifdef DEBUGGING_SMALLMOL
        if (checkResult)
        {
          if (site1Index == relatedProteinSiteNum + 1)
          {
            printf("(%d start, %d %d) satisfy constraint %d\n", site1Index, site2Index, rotOnSite2, consIndex);
          }
          else if (site2Index == relatedProteinSiteNum + 1)
          {
            printf("(%d %d, %d start) satisfy constraint %d\n", site1Index, rotOnSite1, site2Index, consIndex);
          }
          if (site1Index == relatedProteinSiteNum)
          {
            printf("(%d ligand, %d %d) satisfy constraint %d\n", site1Index, site2Index, rotOnSite2, consIndex);
          }
          else if (site2Index == relatedProteinSiteNum)
          {
            printf("(%d %d, %d ligand) satisfy constraint %d\n", site1Index, rotOnSite1, site2Index, consIndex);
          }
        }

#endif
      }

      if (checkResult == FALSE)
      {
        break;
      }
    }

    if (checkResult == FALSE)
    {
      continue;
    }

    checkResult = PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence(pAction, relatedProteinSiteNum,
      relatedProteinSites, rotamersSelectedOnPreviousSteps, pSmallMolRotamer, pStartingProteinRotamer,
      currentStepIndex + 1);
    if (checkResult)
    {
      satisfied = TRUE;
#ifdef DEBGGING_PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence
      continue;
#else
      break;
#endif
    }
  }
  return satisfied;
#undef DEBGGING_PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence
}


BOOL   PlacingRuleProcess_CHECK_MULTI_CONS(PlacingRule* pThis,
  Rotamer* pProteinRotamerOnStartingSite,
  int relatedProteinSiteNum,
  DesignSite** relatedProteinSites,
  PlacingAction* pAction)
{
  int checkResult;
  int* rotamersSelectedOnPreviousSteps = (int*)malloc(sizeof(int) * (relatedProteinSiteNum + 2));
  Rotamer smallMolRotamer;
  RotamerCreate(&smallMolRotamer);
  RotamerSetType(&smallMolRotamer, ResidueGetName(pThis->pSmallMol));
  RotamerAddAtoms(&smallMolRotamer, ResidueGetAllAtoms(pThis->pSmallMol));
  XYZArrayCopy(&smallMolRotamer.xyzs, &pThis->smallMolAtomXYZs);

  for (int i = 0;i < relatedProteinSiteNum + 2;i++)
  {
    rotamersSelectedOnPreviousSteps[i] = -1;
  }

  checkResult = PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence(pAction, relatedProteinSiteNum,
    relatedProteinSites, rotamersSelectedOnPreviousSteps, &smallMolRotamer, pProteinRotamerOnStartingSite, 0);

  RotamerDestroy(&smallMolRotamer);
  free(rotamersSelectedOnPreviousSteps);

  return checkResult;

}


BOOL   PlacingRuleProcess_CHECK_RMSD(PlacingRule* pThis, Rotamers* pSmallMolRotSet, PlacingAction* pAction)
{
  int rotNdx;
  int atomNdx;
  int rotNum = RotamersGetNum(pSmallMolRotSet);
  int atomNum = ResidueGetAtomNum(pThis->pSmallMol);
  BOOL* omitThisAtom;

  if (rotNum == 0)
  {
    return TRUE;
  }

  omitThisAtom = (BOOL*)calloc(atomNum, sizeof(BOOL));
  for (atomNdx = 0; atomNdx < atomNum; atomNdx++)
  {
    if (AtomIsHydrogen(ResidueGetAtom(pThis->pSmallMol, atomNdx)) &&
      (!pAction->checkRMSD_withHydrogen))
    {
      omitThisAtom[atomNdx] = TRUE;
    }
    else if (IntArrayGet(&pAction->checkRMSD_smallmolAtomHasXyz, atomNdx) == FALSE)
    {
      omitThisAtom[atomNdx] = TRUE;
    }
    else
    {
      omitThisAtom[atomNdx] = FALSE;
    }
  }

  for (rotNdx = 0; rotNdx < rotNum; rotNdx++)
  {
    double totRMSD = 0.0;
    int totAtomNum = 0;
    Rotamer* pRot = RotamersGet(pSmallMolRotSet, rotNdx);
    for (atomNdx = 0; atomNdx < atomNum; atomNdx++)
    {
      XYZ* atomInRotSet;
      XYZ* atomInNewRot;
      if (omitThisAtom[atomNdx])
      {
        continue;
      }
      atomInRotSet = XYZArrayGet(&pRot->xyzs, atomNdx);
      atomInNewRot = XYZArrayGet(&pThis->smallMolAtomXYZs, atomNdx);
      totRMSD += pow(XYZDist(atomInRotSet, atomInNewRot), 2);
      totAtomNum++;
    }
    totRMSD = sqrt(totRMSD / totAtomNum);
    //printf("%f ",totalRMSD);
    if (totRMSD < pAction->checkRMSD_minDifference)
    {
      free(omitThisAtom);
      return FALSE;
    }
  }

  free(omitThisAtom);
  return TRUE;
}


BOOL   PlacingRuleProcess_CHECK_VDW_BACKBONE(PlacingRule* pThis, PlacingAction* pAction)
{
  Atom* pAtomOnSmallMol;
  Atom* pAtomOnBackbone;
  pAction->checkVDW_backbone_totalEnergy = 0.0;

  for (int i = 0;i < AtomArrayGetNum(pThis->pTruncBone);i++)
  {
    pAtomOnBackbone = AtomArrayGet(pThis->pTruncBone, i);
    if (pAction->checkVDW_backbone_withHydrogen == FALSE && AtomIsHydrogen(pAtomOnBackbone))
    {
      continue;
    }
    for (int j = 0;j < ResidueGetAtomNum(pThis->pSmallMol);j++)
    {
      double dist;
      double rminSum;
      double ratio;
      pAtomOnSmallMol = ResidueGetAtom(pThis->pSmallMol, j);
      if (pAction->checkVDW_backbone_withHydrogen == FALSE && AtomIsHydrogen(pAtomOnSmallMol))
      {
        continue;
      }
      if (IntArrayGet(&pAction->checkVDW_backbone_smallmolAtomHasXyz, j) == 0)
      {
        continue;
      }
      dist = XYZDist(&(pAtomOnBackbone->xyz), XYZArrayGet(&pThis->smallMolAtomXYZs, j));
      rminSum = 0.95 * (pAtomOnBackbone->vdw_rad + pAtomOnSmallMol->vdw_rad);

      ratio = dist / rminSum;
      if (ratio < 0.8909)
      {
        double B6 = pow(1 / ratio, 6.0);
        double A12 = B6 * B6;
        double epsilon = sqrt(pAtomOnBackbone->vdw_eps * pAtomOnSmallMol->vdw_eps);
        double energy = epsilon * (A12 - 2.0 * B6);
        pAction->checkVDW_backbone_totalEnergy += energy;
        //pAction->checkVDW_backbone_totalPotential += 10.0 - 11.225*ratio ;
      }
      pThis->vdwBackbone = pAction->checkVDW_backbone_totalEnergy;

    }

    if (pAction->checkVDW_backbone_totalEnergy > pAction->checkVDW_backbone_maxAllowed)
    {
      return FALSE;
    }
  }

  return TRUE;
}


BOOL   PlacingRuleProcess_CHECK_VDW_INTERNAL(PlacingRule* pThis, PlacingAction* pAction)
{
  Atom* pAtomI;
  Atom* pAtomJ;
  pAction->checkVDW_internal_totalEnergy = 0.0;

  for (int i = 0;i < pAction->checkVDW_internal_smallMolAtomNum;i++)
  {
    pAtomI = ResidueGetAtom(pThis->pSmallMol, i);
    if (pAction->checkVDW_internal_withHydrogen == FALSE && AtomIsHydrogen(pAtomI))
    {
      continue;
    }
    if (IntArrayGet(&pAction->checkVDW_internal_smallmolAtomHasXyz, i) == 0)
    {
      //This atom's Xyz has not been calculated yet
      continue;
    }

    for (int j = i + 1;j < pAction->checkVDW_internal_smallMolAtomNum;j++)
    {
      double dist;
      double rminSum;
      double ratio;
      pAtomJ = ResidueGetAtom(pThis->pSmallMol, j);
      if (pAction->checkVDW_internal_withHydrogen == FALSE && AtomIsHydrogen(pAtomJ))
      {
        continue;
      }
      if (IntArrayGet(&pAction->checkVDW_internal_smallmolAtomHasXyz, j) == 0)
      {
        //This atom's Xyz has not been calculated yet
        continue;
      }

      //If atomI and atomJ is bonded or 1-3 bonded, continue
      if (pAction->checkVDW_internal_smallMolAtom13bondedMatrix[i][j] == TRUE)
      {
        continue;
      }
      //A rough estimation of VDW potential
      dist = XYZDist(XYZArrayGet(&pThis->smallMolAtomXYZs, i), XYZArrayGet(&pThis->smallMolAtomXYZs, j));
      rminSum = 0.95 * (pAtomI->vdw_rad + pAtomJ->vdw_rad);
      ratio = dist / rminSum;

      if (ratio < 0.8909)
      {
        double B6 = pow(1 / ratio, 6.0);
        double A12 = B6 * B6;
        double epsilon = sqrt(pAtomI->vdw_eps * pAtomJ->vdw_eps);
        double energy = epsilon * (A12 - 2.0 * B6);
        pAction->checkVDW_internal_totalEnergy += energy;
        //pAction->checkVDW_internal_totalPotential += 10.0 - 11.225*ratio;
      }
      //Debug
      pThis->vdwInternal = pAction->checkVDW_internal_totalEnergy;

      if (pAction->checkVDW_internal_totalEnergy > pAction->checkVDW_internal_maxAllowed)
      {
        return FALSE;
      }
    }
  }
  return TRUE;
}


int    PlacingRuleProcess_EVALUATE(PlacingRule* pThis, PlacingAction* pAction)
{
  XYZ atomXYZs[4];
  double newValue = 0.0;
  XYZ vectorI;
  XYZ vectorJ;

  for (int i = 0; i < IntArrayGetLength(&pAction->evaluate_atoms); i++)
  {
    atomXYZs[i] = *XYZArrayGet(&pThis->atomXYZs, IntArrayGet(&pAction->evaluate_atoms, i));
  }

  switch (pAction->evaluate_type)
  {
  case Type_PlacingActionEvaluate_HalfDistance:
    newValue = XYZDist(&atomXYZs[0], &atomXYZs[1]) * 0.5;
    break;
  case Type_PlacingActionEvaluate_Distance:
    newValue = XYZDist(&atomXYZs[0], &atomXYZs[1]);
    break;
  case Type_PlacingActionEvaluate_Angle:
    vectorI = XYZDiff(&atomXYZs[1], &atomXYZs[0]);
    vectorJ = XYZDiff(&atomXYZs[1], &atomXYZs[2]);
    newValue = XYZAngle(&vectorI, &vectorJ);
    break;
  case Type_PlacingActionEvaluate_Torsion:
    newValue = GetTorsionAngle(
      &atomXYZs[0], &atomXYZs[1], &atomXYZs[2], &atomXYZs[3]);
    break;
  default:
    //Should not get here
    break;
  }
  DoubleArraySet(&pThis->params, pAction->evaluate_param, newValue);
  return Success;
}


int    PlacingRuleProcess_LOAD(PlacingRule* pThis,
  Rotamer* pProteinRotamerOnStartingSite,
  PlacingAction* pAction)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  int    result;
  for (int i = 0;i < IntArrayGetLength(&pAction->load_atoms);i++)
  {
    char* atomName;
    Atom* pAtomOnProteinOrSmallMol;
    int index = IntArrayGet(&pAction->load_atoms, i);
    atomName = StringArrayGet(&pThis->atomNames, index);

    if (atomName[0] == SMALLMOL_ATOM_SYMBOL)
    {
      //Find Atom On SmallMol
      pAtomOnProteinOrSmallMol = ResidueGetAtomByName(pThis->pSmallMol, atomName+1);
      //use 'atomName+1' instead of 'atomName' to skip the '$' symbol
      if (pAtomOnProteinOrSmallMol != NULL)
      {
        int atomPosOnSmallMol = IntArrayGet(&pThis->atomPosOnSmallMol, index);
        XYZArraySet(&pThis->atomXYZs, index, &pAtomOnProteinOrSmallMol->xyz);
        //The atom is on small mol, additional operations must be taken to maintain consistency
        XYZArraySet(&pThis->smallMolAtomXYZs, atomPosOnSmallMol, &pAtomOnProteinOrSmallMol->xyz);
        continue; //continue 'for'
      }
    }
    else
    {
      //Find Atom On Protein Rotamer
      Atom* pTempAtom = RotamerGetAtomByName(pProteinRotamerOnStartingSite, atomName);
      if (pTempAtom != NULL)
      {
        XYZArraySet(&pThis->atomXYZs, index, &pTempAtom->xyz);
        continue; //continue 'for'
      }
    }

    //Cannot find the atom, report an error
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d, cannot find atom %s anywhere, when processing the 'LOAD' step", __FILE__, __LINE__, atomName);
    TraceError(errMsg, result);
    PlacingRuleShowAction(pThis, (int)(pAction - pThis->actions));
    return result;
  }

  return Success;
}


int PlacingRuleProcess(PlacingRule* pRule, Rotamer* pStartRot, CataConsSitePairArray* pConsArray, int siteNum, DesignSite** ppSites, Rotamers* pSmallSet, int step)
{
  int result = Success;
  for (; step < pRule->actionNum; step++)
  {
    PlacingAction* pCurAction = &pRule->actions[step];
    BOOL cataConsSatisfied = FALSE;
    BOOL rmsdSatisfied = FALSE;
    BOOL vdwSatisfied = FALSE;

    if (pCurAction->actionType == Type_PlacingAction_Variate) break;
    switch (pCurAction->actionType)
    {
    case Type_PlacingAction_Calc:
      result = PlacingRuleProcess_CALC(pRule, pCurAction);
      break;
    case Type_PlacingAction_CheckCataCons:
      cataConsSatisfied = PlacingRuleProcess_CHECK_CATA_CONS(pRule, pStartRot, siteNum, ppSites, pCurAction);
      result = Success;
      break;
    case Type_PlacingAction_CheckMultiCons:
      cataConsSatisfied = PlacingRuleProcess_CHECK_MULTI_CONS(pRule, pStartRot, siteNum, ppSites, pCurAction);
      result = Success;
      break;
    case Type_PlacingAction_CheckRMSD:
      rmsdSatisfied = PlacingRuleProcess_CHECK_RMSD(pRule, pSmallSet, pCurAction);
      result = Success;
      break;
    case Type_PlacingAction_CheckVDW_Backbone:
      vdwSatisfied = PlacingRuleProcess_CHECK_VDW_BACKBONE(pRule, pCurAction);
      result = Success;
      break;
    case Type_PlacingAction_CheckVDW_Internal:
      vdwSatisfied = PlacingRuleProcess_CHECK_VDW_INTERNAL(pRule, pCurAction);
      result = Success;
      break;
    case Type_PlacingAction_Evaluate:
      result = PlacingRuleProcess_EVALUATE(pRule, pCurAction);
      break;
    case Type_PlacingAction_Load:
      result = PlacingRuleProcess_LOAD(pRule, pStartRot, pCurAction);
      break;
    default:
      printf("small molecule process type %d in placing rule can not be handled\n", pCurAction->actionType);
      result = FormatError;
      break;
    }

    if (FAILED(result))
    {
      return result;
    }
    if (pCurAction->actionType == Type_PlacingAction_CheckCataCons && !cataConsSatisfied)
    {
      //printf("single-constraint does not satisfy at step %d\n", step+1);
      return Success;
    }
    if (pCurAction->actionType == Type_PlacingAction_CheckMultiCons && !cataConsSatisfied)
    {
      //printf("multi-constraints do not satisfy at step %d\n", step+1);
      return Success;
    }
    if (pCurAction->actionType == Type_PlacingAction_CheckRMSD && !rmsdSatisfied)
    {
      return Success;
    }
    if (pCurAction->actionType == Type_PlacingAction_CheckVDW_Backbone && !vdwSatisfied)
    {
      return Success;
    }
    if (pCurAction->actionType == Type_PlacingAction_CheckVDW_Internal && !vdwSatisfied)
    {
      return Success;
    }

  }

  if (step == pRule->actionNum)
  {
    Rotamer newRot;
    RotamerCreate(&newRot);
    RotamerAddAtoms(&newRot, ResidueGetAllAtoms(pRule->pSmallMol));
    strcpy(newRot.type, ResidueGetName(pRule->pSmallMol));
    RotamerCopyAtomXYZ(&newRot, &pRule->smallMolAtomXYZs);
    RotamerAtomsCopyXYZFromArray(&newRot, &pRule->smallMolAtomXYZs);
    for (int i = 0; i < RotamerGetAtomNum(&newRot); i++)
    {
      RotamerGetAtom(&newRot, i)->isXyzValid = pRule->xyzValidArray[i];
    }
    BondsCopy(&newRot.bonds, &pRule->pSmallMol->bonds);
    newRot.vdwBackbone = pRule->vdwBackbone;
    newRot.vdwInternal = pRule->vdwInternal;
    RotamersAdd(pSmallSet, &newRot);
    RotamerDestroy(&newRot);
    return Success;
  }

  if (pRule->actions[step].actionType == Type_PlacingAction_Variate)
  {
    for (double curValue = pRule->actions[step].variate_from; curValue < pRule->actions[step].variate_to; curValue += pRule->actions[step].variate_increment)
    {
      DoubleArraySet(&pRule->params, pRule->actions[step].variate_param, curValue);
      result = PlacingRuleProcess(pRule,
        pStartRot,
        pConsArray,
        siteNum,
        ppSites,
        pSmallSet, step + 1);
      if (FAILED(result)) return result;
    }
  }

  return Success;
}


int PlacingRulePlaceSmallMol(PlacingRule* pThis, CataConsSitePairArray* pConsArray, DesignSite* pStartSite, int siteNum, DesignSite** ppSites, Rotamers* pSmallSet)
{
  int result = Success;
  char errMsg[MAX_LEN_ERR_MSG+1];

  if (PlacingRuleGetDeployedFlag(pThis) == FALSE)
  {
    sprintf(errMsg, "in file %s line %d, placing rule on %s%d%s was not deployed", __FILE__, __LINE__, pThis->chainName, pThis->posInChain, pThis->residueName);
    result = ValueError;
    TraceError(errMsg, result);
    return result;
  }

  time_t begin = time(NULL);
  int progressBarWidth = 100;
  ShowProgress(progressBarWidth, 0.0);
  for (int i = 0; i < RotamersGetNum(DesignSiteGetRotamers(pStartSite)); i++)
  {
    Rotamer* pCurRot = RotamersGet(DesignSiteGetRotamers(pStartSite), i);
    RotamerRestore(pCurRot, DesignSiteGetRotamers(pStartSite));

    if (strcmp(pThis->rotamerType, RotamerGetType(pCurRot)) == 0)
    {
      result = PlacingRuleProcess(pThis, pCurRot, pConsArray, siteNum, ppSites, pSmallSet, 0);
      if (FAILED(result))
      {
        sprintf(errMsg, "in file %s line %d, when placing ligand from site %s%d%s with rotamer type %s", __FILE__, __LINE__, pThis->chainName, pThis->posInChain, pThis->residueName, pThis->rotamerType);
        TraceError(errMsg, result);
        return result;
      }
    }
    RotamerExtract(pCurRot);

    time_t now = time(NULL);
    int elapsed = (int)(now - begin);
    double percentage = (double)(i + 1) / (double)(RotamersGetNum(DesignSiteGetRotamers(pStartSite))) * 100.0;
    ShowProgress(progressBarWidth, percentage);
    printf("time elapsed: %2d h %2d m %2d s; %d rots created\n", elapsed / 3600, (elapsed % 3600) / 60, (elapsed % 3600) % 60, RotamersGetNum(pSmallSet));
  }
  return Success;
}


//Methods used for debugging
int    PlacingRuleShowAtom(PlacingRule* pThis, int index)
{
  int value;

  if (index < 0 || index >= StringArrayGetNum(&pThis->atomNames))
    return IndexError;

  printf("ATOM  %2d  %5.5s  ", index, StringArrayGet(&pThis->atomNames, index));

  XYZShow(XYZArrayGet(&pThis->atomXYZs, index));
  if (pThis->deployedFlag)
  {
    value = IntArrayGet(&pThis->atomPosOnSmallMol, index);
    printf("  %2d  ", value);
    if (value != -1)
    {
      printf("%5.5s  ", AtomGetName(ResidueGetAtom(pThis->pSmallMol, value)));
      XYZShow(XYZArrayGet(&pThis->smallMolAtomXYZs, value));
    }
  }
  printf("\n");
  return Success;
}


int    PlacingRuleShowAtomDistances(PlacingRule* pThis)
{
  int num = XYZArrayGetNum(&pThis->atomXYZs);
  for (int i = -1;i < num;i++)
  {
    if (i == -1)
    {
      printf("%5.5s|", "");
    }
    else
    {
      printf("%5.5s|", StringArrayGet(&pThis->atomNames, i));
    }
    for (int j = 0;j < num;j++)
    {
      if (i == -1)
      {
        printf("%5.5s|", StringArrayGet(&pThis->atomNames, j));
        continue;
      }

      if (i != j)
      {
        printf("%5.2f|", XYZDist(
          XYZArrayGet(&pThis->atomXYZs, i),
          XYZArrayGet(&pThis->atomXYZs, j))
        );
      }
      else
      {
        printf("  \\  |");
      }
    }
    printf("\n");
  }
  return Success;
}


int    PlacingRuleShowParam(PlacingRule* pThis, int index)
{
  if (index < 0 || index >= StringArrayGetNum(&pThis->paramNames))
    return IndexError;
  printf("PARAM  %2d  %8.8s  %f\n", index,
    StringArrayGet(&pThis->paramNames, index),
    DoubleArrayGet(&pThis->params, index));
  return Success;
}


int    PlacingRuleShowAction(PlacingRule* pThis, int index)
{
  PlacingAction* pCurAction;

  if (index < 0 || index >= pThis->actionNum)
    return IndexError;

  pCurAction = &pThis->actions[index];
  switch (pCurAction->actionType)
  {
  case Type_PlacingAction_Load:
    for (int j = 0;j < IntArrayGetLength(&pCurAction->load_atoms);j++)
    {
      index = IntArrayGet(&pCurAction->load_atoms, j);
      printf("LOAD  %s ", StringArrayGet(&pThis->atomNames, index));
      XYZShow(XYZArrayGet(&pThis->atomXYZs, index));
      printf("\n");
    }
    break;
  case Type_PlacingAction_Variate:
    printf("VARIATE  %s  %f  %f  %f\n",
      StringArrayGet(&pThis->paramNames, pCurAction->variate_param),
      pCurAction->variate_from,
      pCurAction->variate_to, pCurAction->variate_increment);
    break;
  case Type_PlacingAction_Calc:
    printf("CALC   ");
    for (int j = 0;j < IntArrayGetLength(&pCurAction->calc_atoms);j++)
    {
      index = IntArrayGet(&pCurAction->calc_atoms, j);
      printf("%4.4s ", StringArrayGet(&pThis->atomNames, index));
    }
    for (int j = 0;j < IntArrayGetLength(&pCurAction->calc_params);j++)
    {
      double paramValue;
      index = IntArrayGet(&pCurAction->calc_params, j);
      paramValue = DoubleArrayGet(&pThis->params, index);
      if (strcmp(StringArrayGet(&pThis->paramNames, index), "") != 0)
      {
        printf("%8.8s", StringArrayGet(&pThis->paramNames, index));
      }
      else
      {
        if (1 <= j && j <= 3)
        {
          printf("%8.2f", RadToDeg(paramValue));
        }
        else
        {
          printf("%8.2f", paramValue);
        }
      }
      printf("  ");
    }
    printf("\n");
    break;
  case Type_PlacingAction_Evaluate:
    printf("EVALUATE  ");
    printf("%s  %d  ",
      StringArrayGet(&pThis->paramNames, pCurAction->evaluate_param),
      pCurAction->evaluate_type);
    for (int j = 0;j < IntArrayGetLength(&pCurAction->evaluate_atoms);j++)
    {
      index = IntArrayGet(&pCurAction->evaluate_atoms, j);
      printf("%4.4s ", StringArrayGet(&pThis->atomNames, index));
    }
    printf("\n");
    break;
  case Type_PlacingAction_CheckCataCons:
    printf("CHECK_CATA_CONS  ");
    printf("%s  %s  %d\t%s	%s	%d\n",
      pCurAction->checkCataCons_firstSiteChainName,
      pCurAction->checkCataCons_firstSiteResidueName,
      pCurAction->checkCataCons_firstSitePosInChain,
      pCurAction->checkCataCons_secondSiteChainName,
      pCurAction->checkCataCons_secondSiteResidueName,
      pCurAction->checkCataCons_secondSitePosInChain);
    break;
  case Type_PlacingAction_CheckMultiCons:
    printf("CHECK_MULTI_CONS_BEGIN\n");
    for (int j = 0;j < pCurAction->checkMultiCons_cataConsNum;j++)
    {
      printf("%s  %s  %d\t%s	%s	%d\n",
        StringArrayGet(&pCurAction->checkMultiCons_firstSiteChainNames, j),
        StringArrayGet(&pCurAction->checkMultiCons_firstSiteResidueNames, j),
        IntArrayGet(&pCurAction->checkMultiCons_firstSitePosInChains, j),
        StringArrayGet(&pCurAction->checkMultiCons_secondSiteChainNames, j),
        StringArrayGet(&pCurAction->checkMultiCons_secondSiteResidueNames, j),
        IntArrayGet(&pCurAction->checkMultiCons_secondSitePosInChains, j));
    }
    printf("CHECK_MULTI_CONS_END\n");
    break;
  case Type_PlacingAction_CheckVDW_Backbone:
    printf("CHECK_VDW_BACKBONE   %s  %f\n",
      pCurAction->checkVDW_backbone_withHydrogen ? "WITH_HYDROGEN" : "WITHOUT_HYDROGEN",
      pCurAction->checkVDW_backbone_maxAllowed);
    for (int j = 0;j < IntArrayGetLength(&pCurAction->checkVDW_backbone_smallmolAtomHasXyz);j++)
    {
      Atom* pAtom = ResidueGetAtom(pThis->pSmallMol, j);
      BOOL include = IntArrayGet(&pCurAction->checkVDW_backbone_smallmolAtomHasXyz, j);
      if (include)
      {
        printf("%s ", AtomGetName(pAtom));
      }
    }printf("\n");
    break;
  case Type_PlacingAction_CheckVDW_Internal:
    printf("CHECK_VDW_INTERNAL   %s  %f\n",
      pCurAction->checkVDW_internal_withHydrogen ? "WITH_HYDROGEN" : "WITHOUT_HYDROGEN",
      pCurAction->checkVDW_internal_maxAllowed);
    for (int j = 0;j < IntArrayGetLength(&pCurAction->checkVDW_internal_smallmolAtomHasXyz);j++)
    {
      Atom* pAtom = ResidueGetAtom(pThis->pSmallMol, j);
      BOOL include = IntArrayGet(&pCurAction->checkVDW_internal_smallmolAtomHasXyz, j);
      if (include)
      {
        printf("%s ", AtomGetName(pAtom));
      }
    }printf("\n");
    break;
  case Type_PlacingAction_CheckRMSD:
    printf("CHECK_RMSD   %s  %f\n",
      pCurAction->checkRMSD_withHydrogen ? "WITH_HYDROGEN" : "WITHOUT_HYDROGEN",
      pCurAction->checkRMSD_minDifference);
    break;
  default:
    printf("small molecule placing action cannot be handled\n");
  }
  return Success;
}


int PlacingRuleShow(PlacingRule* pThis)
{
  for (int i = 0;i < StringArrayGetNum(&pThis->atomNames);i++)
  {
    PlacingRuleShowAtom(pThis, i);
  }
  for (int i = 0;i < StringArrayGetNum(&pThis->paramNames);i++)
  {
    PlacingRuleShowParam(pThis, i);
  }
  for (int i = 0;i < pThis->actionNum;i++)
  {
    PlacingRuleShowAction(pThis, i);
  }
  return Success;
}


int PlacingRuleTester(char* file)
{
  PlacingRule rule;
  int result = PlacingRuleCreate(&rule, file);
  if (FAILED(result))
  {
    return result;
  }

  PlacingRuleDestroy(&rule);
  return Success;
}


int ScreenSmallmolRotamersByRMSD(char* oldFile, char* newFile, double rmsdcut)
{
  FILE* pIn = fopen(oldFile, "r");
  FILE* pOut = fopen(newFile, "w");
  char errMsg[MAX_LEN_ERR_MSG+1];

  if (pIn == NULL || pOut == NULL)
  {
    sprintf(errMsg, "in file %s line %d, cannot open %s or %s", __FILE__, __LINE__, oldFile, newFile);
    TraceError(errMsg, IOError);
    return IOError;
  }

  int atomNum = 0;
  char line[MAX_LEN_ONE_LINE_CONTENT];
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    char keyword[10];
    ExtractFirstStringFromSourceString(keyword, line);
    if (strcmp(keyword, "ATOM") == 0) atomNum++;
    else if (strcmp(keyword, "ENDMDL") == 0) break;
  }
  Rotamer currentRotamer;
  RotamerCreate(&currentRotamer);
  XYZArrayResize(&currentRotamer.xyzs, atomNum);

  fseek(pIn, 0, SEEK_SET);

  int totRotNum = 0;
  BOOL lastAccepted = FALSE;
  Rotamers acceptedRotamers;
  RotamersCreate(&acceptedRotamers);
  StringArray buffer;
  StringArrayCreate(&buffer);
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    char keyword[10] = "";
    ExtractTargetStringFromSourceString(keyword, line, 0, 4);
    if (strcmp(keyword, "ENDM") == 0)
    {
      lastAccepted = TRUE;
      for (int i = 0;i < RotamersGetNum(&acceptedRotamers);i++)
      {
        double rmsd = XYZArrayRMSD(&RotamersGet(&acceptedRotamers, i)->xyzs, &currentRotamer.xyzs);
        if (rmsd < rmsdcut)
        {
          lastAccepted = FALSE;
          break;
        }
      }
      if (lastAccepted)
      {
        fprintf(pOut, "MODEL     %d\n", RotamersGetNum(&acceptedRotamers) + 1);
        for (int i = 0;i < StringArrayGetNum(&buffer);i++)
        {
          fprintf(pOut, "%s", StringArrayGet(&buffer, i));
        }
        fprintf(pOut, "ENDMDL\n");
        RotamersAdd(&acceptedRotamers, &currentRotamer);
      }

      totRotNum++;
      printf("\r %d / %d rotamers processed", RotamersGetNum(&acceptedRotamers), totRotNum);
      //fflush(stdout);
    }
    else if (strcmp(keyword, "MODE") == 0)
    {
      atomNum = 0;
      StringArrayDestroy(&buffer);
      StringArrayCreate(&buffer);
    }
    else if (strcmp(keyword, "ATOM") == 0)
    {
      XYZ* pXYZ = XYZArrayGet(&currentRotamer.xyzs, atomNum);
      ExtractTargetStringFromSourceString(keyword, line, 31, 7);
      pXYZ->X = atof(keyword);
      ExtractTargetStringFromSourceString(keyword, line, 39, 7);
      pXYZ->Y = atof(keyword);
      ExtractTargetStringFromSourceString(keyword, line, 47, 7);
      pXYZ->Z = atof(keyword);
      StringArrayAppend(&buffer, line);
      atomNum++;
    }
    else if (strcmp(keyword, "ENER") == 0)
    {
      StringArrayAppend(&buffer, line);
    }
  }

  StringArrayDestroy(&buffer);
  RotamerDestroy(&currentRotamer);
  RotamersDestroy(&acceptedRotamers);
  fclose(pIn);
  fclose(pOut);

  return Success;
}


int AnalyzeSmallMolRotamers(char* rotFile, Residue* pSmallMol)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  char line[MAX_LEN_ONE_LINE_CONTENT];

  FILE* pIn = fopen(rotFile, "r");
  if (pIn == NULL)
  {
    sprintf(errMsg, "in file %s line %d, cannot read file %s", __FILE__, __LINE__, rotFile);
    TraceError(errMsg, IOError);
    return IOError;
  }

  int atomNum = 0;
  XYZArray currentRotamerXYZ;
  IntArray atomPosOnNativeSmallMol;
  IntArrayCreate(&atomPosOnNativeSmallMol, 0);
  XYZArrayCreate(&currentRotamerXYZ, 0);
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    char keyword[10];
    char atomName[10];
    ExtractTargetStringFromSourceString(keyword, line, 0, 6);
    if (strcmp(keyword, "ATOM") == 0)
    {
      atomNum++;
      ExtractTargetStringFromSourceString(atomName, line, 12, 4);
      // 05/08/2023: add this if(){} to reformat the atom name
      if (isdigit(atomName[0]) && isalpha(atomName[1]))
      {
        char tempName[MAX_LEN_ATOM_NAME+1];
        strcpy(tempName, atomName+1);
        tempName[strlen(atomName)-1] = atomName[0];
        tempName[strlen(atomName)] = '\0';
        strcpy(atomName, tempName);
      }
      int posOnNativeSmallMol = -1;
      ResidueFindAtom(pSmallMol, atomName, &posOnNativeSmallMol);
      if (posOnNativeSmallMol == -1)
      {
        sprintf(errMsg, "in file %s line %d, cannot find atom %s on smallmol", __FILE__, __LINE__, atomName);
        TraceError(errMsg, DataNotExistError);
        return DataNotExistError;
      }
      IntArrayAppend(&atomPosOnNativeSmallMol, posOnNativeSmallMol);
    }
    else if (strcmp(keyword, "ENDMDL") == 0) break;
  }
  XYZArrayResize(&currentRotamerXYZ, atomNum);
  fseek(pIn, 0, SEEK_SET);

  int poseNum = 0;
  // minRMSD pose
  double minRMSD = 1e8;
  int indexOfMinRMSD = -1;
  double internalVDWOfMinRMSD = 1e8;
  double backboneVDWOfMinRMSD = 1e8;
  // minInternalVDW pose
  double minInternalVDW = 1e8;
  int indexOfMinInternalVDW = -1;
  double backboneVDWOfMinInternalVDW = 1e8;
  double rmsdOfMinInternalVDW = 1e8;

  //minBackboneVDW pose
  double minBackboneVDW = 1e8;
  int indexOfMinBackboneVDW = -1;
  double internalVDWOfMinBackboneVDW = 1e8;
  double rmsdOfMinBackboneVDW = 1e8;

  //number of poses below 0.5, 1, 1.5, 2.0, 2.5, and 3.0 angstroms;
  int numOfPoseBelow0_5A = 0;
  int numOfPoseBelow1_0A = 0;
  int numOfPoseBelow1_5A = 0;
  int numOfPoseBelow2_0A = 0;
  int numOfPoseBelow2_5A = 0;
  int numOfPoseBelow3_0A = 0;

  double internalVDW = 1e8;
  double backboneVDW = 1e8;
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    char keyword[MAX_LEN_ONE_LINE_CONTENT] = "";
    ExtractTargetStringFromSourceString(keyword, line, 0, 4);
    if (strcmp(keyword, "ENDM") == 0)
    {
      double rmsd = 0.0;
      for (int i = 0;i < atomNum;i++)
      {
        Atom* pNativeAtom = ResidueGetAtom(pSmallMol, IntArrayGet(&atomPosOnNativeSmallMol, i));
        double dist = XYZDist(XYZArrayGet(&currentRotamerXYZ, i), &pNativeAtom->xyz);
        rmsd += dist * dist;
      }
      rmsd = sqrt(rmsd / atomNum);

      if (rmsd < minRMSD)
      {
        minRMSD = rmsd;
        indexOfMinRMSD = poseNum;
        internalVDWOfMinRMSD = internalVDW;
        backboneVDWOfMinRMSD = backboneVDW;
      }

      if (rmsd < 0.5)
      {
        numOfPoseBelow0_5A++;
      }
      else if (rmsd < 1.0)
      {
        numOfPoseBelow1_0A++;
      }
      else if (rmsd < 1.5)
      {
        numOfPoseBelow1_5A++;
      }
      else if (rmsd < 2.0)
      {
        numOfPoseBelow2_0A++;
      }
      else if (rmsd < 2.5)
      {
        numOfPoseBelow2_5A++;
      }
      else if (rmsd < 3.0)
      {
        numOfPoseBelow3_0A++;
      }
    }
    else if (strcmp(keyword, "MODE") == 0)
    {
      atomNum = 0;
      poseNum++;
    }
    else if (strcmp(keyword, "ATOM") == 0)
    {
      XYZ* pXYZ = XYZArrayGet(&currentRotamerXYZ, atomNum);
      ExtractTargetStringFromSourceString(keyword, line, 31, 7);
      pXYZ->X = atof(keyword);
      ExtractTargetStringFromSourceString(keyword, line, 39, 7);
      pXYZ->Y = atof(keyword);
      ExtractTargetStringFromSourceString(keyword, line, 47, 7);
      pXYZ->Z = atof(keyword);
      atomNum++;
    }
    else if (strcmp(keyword, "ENER") == 0)
    {
      sscanf(line, "%s %s %lf %s %lf", keyword, keyword, &internalVDW, keyword, &backboneVDW);
      if (internalVDW < minInternalVDW)
      {
        minInternalVDW = internalVDW;
        indexOfMinInternalVDW = poseNum;
        backboneVDWOfMinInternalVDW = backboneVDW;
      }
      if (backboneVDW < minBackboneVDW)
      {
        minBackboneVDW = backboneVDW;
        indexOfMinBackboneVDW = poseNum;
        internalVDWOfMinBackboneVDW = internalVDW;
      }
    }
  }

  printf("Summary of ligand pose analysis:\n");
  printf("----------------------------------------------\n");
  printf("No. of total poses                : %d\n", poseNum);
  printf("No. of poses < 0.5 A              : %d\n", numOfPoseBelow0_5A);
  printf("No. of poses < 1.0 A              : %d\n", numOfPoseBelow1_0A);
  printf("No. of poses < 1.5 A              : %d\n", numOfPoseBelow1_5A);
  printf("No. of poses < 2.0 A              : %d\n", numOfPoseBelow2_0A);
  printf("No. of poses < 2.5 A              : %d\n", numOfPoseBelow2_5A);
  printf("No. of poses < 3.0 A              : %d\n", numOfPoseBelow3_0A);
  printf("----------------------------------------------\n");
  printf("minRMSD (angstroms)               : %f\n", minRMSD);
  printf("Index of minRMSD pose (from 1)    : %d\n", indexOfMinRMSD);
  printf("InternalVDW of minRMSD pose       : %f\n", internalVDWOfMinRMSD);
  printf("BackboneVDW of minRMSD pose       : %f\n", backboneVDWOfMinRMSD);
  printf("----------------------------------------------\n");
  printf("minInternalVDW                    : %f\n", minInternalVDW);
  printf("Index of MinInternalVDW pose      : %d\n", indexOfMinInternalVDW);
  printf("BackboneVDW of MinInternalVDW pose: %f\n", backboneVDWOfMinInternalVDW);
  printf("----------------------------------------------\n");
  printf("minBackboneVDW                    : %f\n", minBackboneVDW);
  printf("Index of MinBackboneVDW pose      : %d\n", indexOfMinBackboneVDW);
  printf("InternalVDW of MinBackboneVDW pose: %f\n", internalVDWOfMinBackboneVDW);

  fclose(pIn);
  XYZArrayDestroy(&currentRotamerXYZ);
  IntArrayDestroy(&atomPosOnNativeSmallMol);

  return Success;
}


int AnalyzeSmallMolRotamersForSpecifiedAtoms(char* poseFile, Residue* pSmall, char* specificAtomFile)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  char line[MAX_LEN_ONE_LINE_CONTENT];
  double minInternalVDW = 1e8;
  double minBackboneVDW = 1e8;
  double internalVDWOfMinRMSD = 1e4;
  double backboneVDWOfMinRMSD = 1e4;
  double minRMSD = 1e8;
  int minRMSDndx = -1;

  FileReader fr;
  StringArray specificAtomNames;

  XYZArray curRotXYZs;
  IntArray atomNdxArray;
  FILE* pIn = fopen(poseFile, "r");
  if (pIn == NULL)
  {
    sprintf(errMsg, "in file %s line %d, cannot open file %s", __FILE__, __LINE__, poseFile);
    TraceError(errMsg, IOError);
    return IOError;
  }
  IntArrayCreate(&atomNdxArray, 0);
  XYZArrayCreate(&curRotXYZs, 0);

  StringArrayCreate(&specificAtomNames);
  FileReaderCreate(&fr, specificAtomFile);
  while (!FAILED(FileReaderGetNextLine(&fr, line)))
  {
    StringArray wordsInLine;
    StringArrayCreate(&wordsInLine);
    StringArraySplitString(&wordsInLine, line, ' ');
    for (int i = 0; i < StringArrayGetNum(&wordsInLine); i++)
    {
      int posOnNativeSmallMol = -1;
      ResidueFindAtom(pSmall, StringArrayGet(&wordsInLine, i), &posOnNativeSmallMol);
      if (posOnNativeSmallMol == -1)
      {
        sprintf(errMsg, "in file %s line %d, cannot find atom %s on ligand", __FILE__, __LINE__, StringArrayGet(&wordsInLine, i));
        TraceError(errMsg, DataNotExistError);
        return DataNotExistError;
      }
      StringArrayAppend(&specificAtomNames, StringArrayGet(&wordsInLine, i));
    }
    StringArrayDestroy(&wordsInLine);
  }

  int atomNum = 0;
  int poseNum = 0;
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    char keyword[10];
    char atomName[10];
    ExtractTargetStringFromSourceString(keyword, line, 0, 6);
    if (strcmp(keyword, "ATOM") == 0)
    {
      ExtractTargetStringFromSourceString(atomName, line, 12, 4);
      if (isdigit(atomName[0]) && isalpha(atomName[1]))
      {
        char temp[MAX_LEN_ATOM_NAME+1];
        strcpy(temp, atomName+1);
        temp[strlen(atomName)-1] = atomName[0];
        temp[strlen(atomName)] = '\0';
        strcpy(atomName, temp);
      }
      for (int i = 0; i < StringArrayGetNum(&specificAtomNames); i++)
      {
        if (strcmp(atomName, StringArrayGet(&specificAtomNames, i)) == 0)
        {
          atomNum++;
          int pos = -1;
          ResidueFindAtom(pSmall, atomName, &pos);
          if (pos == -1)
          {
            sprintf(errMsg, "in file %s line %d, cannot find atom %s on ligand", __FILE__, __LINE__, atomName);
            TraceError(errMsg, DataNotExistError);
            return DataNotExistError;
          }
          IntArrayAppend(&atomNdxArray, pos);
        }
      }
    }
    else if (strcmp(keyword, "ENDMDL") == 0)
    {
      break;
    }
  }

  XYZArrayResize(&curRotXYZs, atomNum);

  fseek(pIn, 0, SEEK_SET);

  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    char keyword[10];
    strcpy(keyword, "");
    ExtractTargetStringFromSourceString(keyword, line, 0, 4);
    if (strcmp(keyword, "ENDM") == 0)
    {
      double rmsd = 0.0;
      for (int i = 0;i < atomNum;i++)
      {
        Atom* pAtom = ResidueGetAtom(pSmall, IntArrayGet(&atomNdxArray, i));
        double dist = XYZDist(XYZArrayGet(&curRotXYZs, i), &pAtom->xyz);
        rmsd += dist * dist;
      }
      rmsd = sqrt(rmsd / atomNum);

      if (rmsd < minRMSD)
      {
        minRMSD = rmsd;
        minRMSDndx = poseNum - 1;
      }
    }
    else if (strcmp(keyword, "MODE") == 0)
    {
      atomNum = 0;
      poseNum++;
    }
    else if (strcmp(keyword, "ATOM") == 0)
    {
      char atomName[MAX_LEN_ATOM_NAME+1];

      ExtractTargetStringFromSourceString(atomName, line, 12, 4);
      for (int i = 0; i < StringArrayGetNum(&specificAtomNames); i++)
      {
        if (strcmp(atomName, StringArrayGet(&specificAtomNames, i)) == 0)
        {
          XYZ* pXYZ = XYZArrayGet(&curRotXYZs, atomNum);
          ExtractTargetStringFromSourceString(keyword, line, 31, 7);
          pXYZ->X = atof(keyword);
          ExtractTargetStringFromSourceString(keyword, line, 39, 7);
          pXYZ->Y = atof(keyword);
          ExtractTargetStringFromSourceString(keyword, line, 47, 7);
          pXYZ->Z = atof(keyword);
          atomNum++;
        }
      }
    }
    else if (strcmp(keyword, "ENER") == 0)
    {
      char internalVDW[20];
      char backboneVDW[20];
      ExtractFirstStringFromSourceString(keyword, line);
      ExtractFirstStringFromSourceString(internalVDW, line);
      ExtractFirstStringFromSourceString(internalVDW, line);
      ExtractFirstStringFromSourceString(backboneVDW, line);
      ExtractFirstStringFromSourceString(backboneVDW, line);
      if (atof(internalVDW) < minInternalVDW)
      {
        minInternalVDW = atof(internalVDW);
      }
      if (atof(backboneVDW) < minBackboneVDW)
      {
        minBackboneVDW = atof(backboneVDW);
      }
      if (poseNum == minRMSDndx + 1)
      {
        internalVDWOfMinRMSD = atof(internalVDW);
        backboneVDWOfMinRMSD = atof(backboneVDW);
      }
    }
    else
    {
      continue;
    }
  }

  printf("Num of Total Poses     : %d\n", poseNum);
  printf("Min RMSD               : %f\n", minRMSD);
  printf("Min RMSD Rotamer Index : %d\n", minRMSDndx);
  printf("InternalVDW of minRMSD : %f\n", internalVDWOfMinRMSD);
  printf("BackboneVDW of minRMSD : %f\n", backboneVDWOfMinRMSD);
  printf("Min InternalVDW        : %f\n", minInternalVDW);
  printf("Min BackboneVDW        : %f\n", minBackboneVDW);

  fclose(pIn);
  XYZArrayDestroy(&curRotXYZs);
  IntArrayDestroy(&atomNdxArray);
  FileReaderDestroy(&fr);
  StringArrayDestroy(&specificAtomNames);

  return Success;
}


typedef struct
{
  double vdwInternal;
  double vdwBackbone;
  int ndx;
} RotamerEnergy;


int CompareByInternalEnergy(const void* a, const void* b)
{
  double e1 = (*(RotamerEnergy*)a).vdwInternal;
  double e2 = (*(RotamerEnergy*)b).vdwInternal;
  return e1 < e2 ? -1 : 1;
}


int CompareByBackboneEnergy(const void* a, const void* b)
{
  double e1 = (*(RotamerEnergy*)a).vdwBackbone;
  double e2 = (*(RotamerEnergy*)b).vdwBackbone;
  return e1 < e2 ? -1 : 1;
}


int SmallMolRotamersGetBothHighRankOfBackboneVdwAndInternalVdw(char* oldRotamersFile, char* newRotamersFile, double bbVdwRankPerc, double internalVdwRankPerc)
{
  int result = Success;
  char errMsg[MAX_LEN_ERR_MSG+1];
  printf("select smallmol poses that internal and backbone VDW energies ranked in top %.0f%% and %.0f%% respectively\n", internalVdwRankPerc*100, bbVdwRankPerc*100);
  FILE* pIn = fopen(oldRotamersFile, "r");
  if (pIn == NULL)
  {
    sprintf(errMsg, "in file %s line %d, cannot open file %s", __FILE__, __LINE__, oldRotamersFile);
    result = IOError;
    TraceError(errMsg, result);
    return result;
  }

  int poseNum = 0;
  char line[MAX_LEN_ONE_LINE_CONTENT+1];
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    char keyword[MAX_LEN_ONE_LINE_CONTENT+1];
    ExtractFirstStringFromSourceString(keyword, line);
    if (strcmp(keyword, "MODEL") == 0)
    {
      poseNum++;
    }
  }
  printf("total smallmol pose number: %d\n", poseNum);

  RotamerEnergy* pRotamerEnergy = (RotamerEnergy*)malloc(sizeof(RotamerEnergy) * poseNum);
  BOOL* flagRotamerWithinRank = (BOOL*)malloc(sizeof(BOOL) * poseNum);
  for (int i=0; i<poseNum; i++)
  {
    flagRotamerWithinRank[i] = FALSE;
  }

  fseek(pIn, 0, SEEK_SET);
  int rotNdx = 0;
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    char keyword[MAX_LEN_ONE_LINE_CONTENT+1];
    ExtractFirstStringFromSourceString(keyword, line);
    if (strcmp(keyword, "ENERGY") != 0) continue;
    ExtractFirstStringFromSourceString(keyword, line);
    ExtractFirstStringFromSourceString(keyword, line);
    double vdwInternal = atof(keyword);
    ExtractFirstStringFromSourceString(keyword, line);
    ExtractFirstStringFromSourceString(keyword, line);
    double vdwBackbone = atof(keyword);
    pRotamerEnergy[rotNdx].vdwInternal = vdwInternal;
    pRotamerEnergy[rotNdx].vdwBackbone = vdwBackbone;
    pRotamerEnergy[rotNdx].ndx = rotNdx;
    rotNdx++;
  }

  printf("rank smallmol rotamers by vdwInternal and vdwBackbone energy\n");
  int bbVdwRankThreshold = (int)(poseNum * bbVdwRankPerc);
  int internalVdwRankThreshold = (int)(poseNum * internalVdwRankPerc);
  // rank by internalVDW;
  IntArray highRankIndexByInternalVDW;
  IntArrayCreate(&highRankIndexByInternalVDW, internalVdwRankThreshold);
  qsort(pRotamerEnergy, poseNum, sizeof(RotamerEnergy), CompareByInternalEnergy);
  rotNdx = 0;
  while (rotNdx < internalVdwRankThreshold)
  {
    IntArraySet(&highRankIndexByInternalVDW, rotNdx, pRotamerEnergy[rotNdx].ndx);
    rotNdx++;
  }

  // rank by backboneVDW;
  IntArray highRankIndexByBBVDW;
  IntArrayCreate(&highRankIndexByBBVDW, bbVdwRankThreshold);
  qsort(pRotamerEnergy, poseNum, sizeof(RotamerEnergy), CompareByBackboneEnergy);
  rotNdx = 0;
  while (rotNdx < bbVdwRankThreshold)
  {
    IntArraySet(&highRankIndexByBBVDW, rotNdx, pRotamerEnergy[rotNdx].ndx);
    rotNdx++;
  }

  // find out the rots that are ranked in top 'highRank' by internalVDW as well as backboneVDW;
  for (int i = 0; i < IntArrayGetLength(&highRankIndexByInternalVDW); i++)
  {
    int indexI = IntArrayGet(&highRankIndexByInternalVDW, i);
    for (int j = 0; j < IntArrayGetLength(&highRankIndexByBBVDW); j++)
    {
      int indexJ = IntArrayGet(&highRankIndexByBBVDW, j);
      if (indexJ == indexI)
      {
        printf("RotamerIndex: %6d, InternalRank: %6d, BackboneRank: %6d\n", indexJ, i, j);
        flagRotamerWithinRank[indexI] = TRUE;
      }
    }
  }

  //write top ranked smallmol rotamers to a new file
  FILE* pOut = fopen(newRotamersFile, "w");
  if (pOut == NULL)
  {
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, newRotamersFile);
    result = IOError;
    TraceError(errMsg, result);
    return result;
  }

  fseek(pIn, 0, SEEK_SET);
  rotNdx = 0;
  int accRotNdx = 0;
  StringArray buffer;
  StringArrayCreate(&buffer);
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    char keyword[MAX_LEN_ONE_LINE_CONTENT+1];
    ExtractTargetStringFromSourceString(keyword, line, 0, 4);
    if (strcmp(keyword, "MODE") == 0)
    {
      rotNdx++;
      StringArrayDestroy(&buffer);
      StringArrayCreate(&buffer);
    }
    else if (strcmp(keyword, "ATOM") == 0)
    {
      StringArrayAppend(&buffer, line);
    }
    else if(strcmp(keyword, "ENER") == 0)
    {
      StringArrayAppend(&buffer, line);
    }
    else if (strcmp(keyword, "ENDM") == 0)
    {
      if (flagRotamerWithinRank[rotNdx-1] == TRUE)
      {
        fprintf(pOut, "MODEL     %d\n", accRotNdx+1);
        for (int j = 0;j < StringArrayGetNum(&buffer);j++)
        {
          fprintf(pOut, "%s", StringArrayGet(&buffer, j));
        }
        fprintf(pOut, "ENDMDL\n");
        accRotNdx++;
      }
    }
  }
  fclose(pOut);
  StringArrayDestroy(&buffer);

  IntArrayDestroy(&highRankIndexByBBVDW);
  IntArrayDestroy(&highRankIndexByInternalVDW);
  free(pRotamerEnergy);
  free(flagRotamerWithinRank);

  return Success;
}

