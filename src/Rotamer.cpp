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

#include "Rotamer.h"
#include <string.h>
#include <ctype.h>

extern double ROTLIB_LOW_PROB_CUTOFF;

int Type_ProteinAtomOrder_ToInt(Type_ProteinAtomOrder order)
{
  // 0 for alpha, 1 for Beta, and etc
  return (int)order;
}

Type_ProteinAtomOrder Type_ProteinAtomOrder_FromInt(int order)
{
  // 0 for alpha, 1 for Beta, and etc
  return (Type_ProteinAtomOrder)order;
}

Type_ProteinAtomOrder Type_ProteinAtomOrder_JudgedByAtomName(char* atName)
{
  int i;
  char sequence[] = "ABGDEZ";
  char canBeJudged = atName[strlen(atName) - 1];

  if (atName[0] == 'H')
    return Type_ProteinAtomOrder_Other;

  if (isdigit(canBeJudged))
  {
    if (canBeJudged == '1')
    {
      canBeJudged = atName[strlen(atName) - 2];
    }
    else
    {
      return Type_ProteinAtomOrder_Other;
    }
  }

  for (i = 0;i < sizeof(sequence) / sizeof(char); i++)
  {
    if (canBeJudged == sequence[i])
    {
      return Type_ProteinAtomOrder_FromInt(i);
    }
  }
  return Type_ProteinAtomOrder_Other;
}

int RotamerCreate(Rotamer* pRot)
{
  strcpy(pRot->type, "");
  AtomArrayCreate(&pRot->atoms);
  XYZArrayCreate(&pRot->xyzs, 0);
  BondsCreate(&pRot->bonds);
  strcpy(pRot->chainName, "");
  pRot->posInChain = -1;
  pRot->vdwInternal = 0;
  pRot->vdwBackbone = 0;
  pRot->selfEnergy = 0;
  pRot->selfEnergyBin = 0;
  pRot->dunbrack = 0;
  DoubleArrayCreate(&pRot->Xs, 0);
  return Success;
}

int RotamerDestroy(Rotamer* pRot)
{
  strcpy(pRot->type, "");
  AtomArrayDestroy(&pRot->atoms);
  XYZArrayDestroy(&pRot->xyzs);
  BondsDestroy(&pRot->bonds);
  DoubleArrayDestroy(&pRot->Xs);
  return Success;
}

int RotamerCopy(Rotamer* pDest, Rotamer* pSrc)
{
  RotamerDestroy(pDest);
  RotamerCreate(pDest);
  strcpy(pDest->type, pSrc->type);
  strcpy(pDest->chainName, pSrc->chainName);
  pDest->posInChain = pSrc->posInChain;
  AtomArrayCopy(&pDest->atoms, &pSrc->atoms);
  XYZArrayCopy(&pDest->xyzs, &pSrc->xyzs);
  BondsCopy(&pDest->bonds, &pSrc->bonds);
  pDest->vdwInternal = pSrc->vdwInternal;
  pDest->vdwBackbone = pSrc->vdwBackbone;
  pDest->selfEnergy = pSrc->selfEnergy;
  pDest->dunbrack = pSrc->dunbrack;
  DoubleArrayCopy(&pDest->Xs, &pSrc->Xs);
  return Success;
}

char* RotamerGetType(Rotamer* pRotlib)
{
  return pRotlib->type;
}

int RotamerSetType(Rotamer* pRotlib, char* type)
{
  if (strlen(type) > MAX_LEN_RES_NAME)
  {
    return NameError;
  }
  strcpy(pRotlib->type, type);
  return Success;
}

int RotamerCopyAtomXYZ(Rotamer* pThis, XYZArray* pNewXYZ)
{
  if (XYZArrayGetNum(&pThis->xyzs) != XYZArrayGetNum(pNewXYZ))
  {
    XYZArrayResize(&pThis->xyzs, XYZArrayGetNum(pNewXYZ));
  }
  XYZArrayCopy(&pThis->xyzs, pNewXYZ);
  return Success;
}

int RotamerAtomsCopyXYZFromArray(Rotamer* pThis, XYZArray* pNewXYZ)
{
  if (AtomArrayGetNum(&pThis->atoms) != XYZArrayGetNum(pNewXYZ))
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    sprintf(errMsg, "in file %s line %d, rotamer's atom number is not equal to the length of XYZArray", __FILE__, __LINE__);
    TraceError(errMsg, ValueError);
    exit(ValueError);
  }
  for (int i = 0; i < AtomArrayGetNum(&pThis->atoms); i++)
  {
    Atom* pAtom = AtomArrayGet(&pThis->atoms, i);
    pAtom->xyz = *XYZArrayGet(pNewXYZ, i);
    pAtom->isXyzValid = TRUE;
  }
  return Success;
}

char* RotamerGetChainName(Rotamer* pThis)
{
  return pThis->chainName;
}

int RotamerSetChainName(Rotamer* pThis, char* newChainName)
{
  if (strlen(newChainName) > MAX_LEN_CHAIN_NAME)
  {
    return NameError;
  }
  else
  {
    strcpy(pThis->chainName, newChainName);
    return Success;
  }
}

int RotamerGetPosInChain(Rotamer* pThis)
{
  return pThis->posInChain;
}

int RotamerSetPosInChain(Rotamer* pThis, int newPosInChain)
{
  pThis->posInChain = newPosInChain;
  return Success;
}

double RotamerGetDunbrack(Rotamer* pThis)
{
  return pThis->dunbrack;
}

int RotamerSetDunbrack(Rotamer* pThis, double dun)
{
  pThis->dunbrack = dun;
  return Success;
}

int RotamerShow(Rotamer* pThis)
{
  printf("Rotamer : %s for %s%d \n", pThis->type, pThis->chainName, pThis->posInChain);
  printf("%d atoms\n", AtomArrayGetNum(&pThis->atoms));
  printf("%d bonds\n", BondsGetNum(&pThis->bonds));
  printf("%d xyzs\n", XYZArrayGetNum(&pThis->xyzs));
  if (AtomArrayGetNum(&pThis->atoms) != 0)
  {
    AtomArrayShowInPDBFormat(&pThis->atoms, "ATOM", pThis->type, " ", 0, 0, NULL);
  }
  BondsShow(&pThis->bonds);
  return Success;
}

int RotamerGetAtomNum(Rotamer* pThis)
{
  if (AtomArrayGetNum(&pThis->atoms) == 0)
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
    sprintf(errMsg, "in file %s line %d, the rotamer %s%d%s may need to be restored", __FILE__, __LINE__,
      RotamerGetChainName(pThis), RotamerGetPosInChain(pThis), RotamerGetType(pThis));
    TraceError(errMsg, Warning);
  }
  return AtomArrayGetNum(&pThis->atoms);
}

Atom* RotamerGetAtom(Rotamer* pThis, int index)
{
  if (AtomArrayGetNum(&pThis->atoms) == 0)
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
    sprintf(errMsg, "in file %s line %d, the rotamer %s%d%s may need to be restored", __FILE__, __LINE__,
      RotamerGetChainName(pThis), RotamerGetPosInChain(pThis), RotamerGetType(pThis));
    TraceError(errMsg, Warning);
  }
  return AtomArrayGet(&pThis->atoms, index);
}

Atom* RotamerGetAtomByName(Rotamer* pThis, char* atomName)
{
  if (AtomArrayGetNum(&pThis->atoms) == 0)
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
    sprintf(errMsg, "in file %s line %d, the rotamer %s%d%s may need to be restored", __FILE__, __LINE__,
      RotamerGetChainName(pThis), RotamerGetPosInChain(pThis), RotamerGetType(pThis));
    TraceError(errMsg, Warning);
  }
  int index;
  int result = RotamerFindAtom(pThis, atomName, &index);
  if (FAILED(result))
  {
    return NULL;
  }
  else
  {
    return RotamerGetAtom(pThis, index);
  }
}

int RotamerFindAtom(Rotamer* pThis, char* atomName, int* index)
{
  if (AtomArrayGetNum(&pThis->atoms) == 0)
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
    sprintf(errMsg, "in file %s line %d, the rotamer %s%d%s may need to be restored", __FILE__, __LINE__,
      RotamerGetChainName(pThis), RotamerGetPosInChain(pThis), RotamerGetType(pThis));
    TraceError(errMsg, Warning);
  }
  return AtomArrayFind(&pThis->atoms, atomName, index);
}

int RotamerAddAtoms(Rotamer* pThis, AtomArray* pNewAtoms)
{
  XYZArrayResize(&pThis->xyzs, AtomArrayGetNum(pNewAtoms));
  for (int i = 0;i < AtomArrayGetNum(pNewAtoms);i++)
  {
    XYZArraySet(&pThis->xyzs, i, &AtomArrayGet(pNewAtoms, i)->xyz);
  }
  AtomArrayCopy(&pThis->atoms, pNewAtoms);
  return Success;
}

Bonds* RotamerGetBonds(Rotamer* pThis)
{
  if (AtomArrayGetNum(&pThis->atoms) == 0)
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
    sprintf(errMsg, "in file %s line %d, the rotamer %s%d%s may need to be restored", __FILE__, __LINE__,
      RotamerGetChainName(pThis), RotamerGetPosInChain(pThis), RotamerGetType(pThis));
    TraceError(errMsg, Warning);
  }
  return &pThis->bonds;
}

double RotamerAndRotamerSidechainRMSD(Rotamer* pThis, Rotamer* pOther)
{
  double rmsd = 0.0;
  int num = 0;
  if (RotamerAndRotamerInSameType(pThis, pOther))
  {
    for (int i = 0; i < RotamerGetAtomNum(pThis); ++i)
    {
      Atom* pAtom1 = RotamerGetAtom(pThis, i);
      if (AtomIsHydrogen(pAtom1) || pAtom1->isBBAtom || strcmp(pAtom1->name, "CB") == 0) continue;
      Atom* pAtom2 = RotamerGetAtomByName(pOther, AtomGetName(pAtom1));
      rmsd += (pAtom1->xyz.X - pAtom2->xyz.X) * (pAtom1->xyz.X - pAtom2->xyz.X) 
            + (pAtom1->xyz.Y - pAtom2->xyz.Y) * (pAtom1->xyz.Y - pAtom2->xyz.Y) 
            + (pAtom1->xyz.Z - pAtom2->xyz.Z) * (pAtom1->xyz.Z - pAtom2->xyz.Z);
      num++;
    }
    if (num > 0) return sqrt(rmsd / num);
    else return 0.0;
  }
  else return 1e8;
}

double RotamerAndResidueSidechainRMSD(Rotamer* pThis, Residue* pOther)
{
  double rmsd = 0.0;
  int num = 0;
  if (RotamerAndResidueInSameType(pThis, pOther))
  {
    for (int i = 0; i < RotamerGetAtomNum(pThis); ++i)
    {
      Atom* pAtom1 = RotamerGetAtom(pThis, i);
      if (AtomIsHydrogen(pAtom1) || pAtom1->isBBAtom || strcmp(pAtom1->name, "CB") == 0) continue;
      Atom* pAtom2 = ResidueGetAtomByName(pOther, AtomGetName(pAtom1));
      rmsd += (pAtom1->xyz.X - pAtom2->xyz.X) * (pAtom1->xyz.X - pAtom2->xyz.X) 
            + (pAtom1->xyz.Y - pAtom2->xyz.Y) * (pAtom1->xyz.Y - pAtom2->xyz.Y) 
            + (pAtom1->xyz.Z - pAtom2->xyz.Z) * (pAtom1->xyz.Z - pAtom2->xyz.Z);
      num++;
    }
    if (num > 0) return sqrt(rmsd / num);
    else return 0.0;
  }
  else return 1e8;
}

int RotamerExtract(Rotamer* pThis)
{
  AtomArrayDestroy(&pThis->atoms);
  BondsDestroy(&pThis->bonds);
  return Success;
}

int RotamerRestore(Rotamer* pThis, Rotamers* pRotamerSet)
{
  int result;
  char errMsg[MAX_LEN_ERR_MSG+1];
  if (AtomArrayGetNum(&pThis->atoms) == XYZArrayGetNum(&pThis->xyzs))
  {
    return Success;
  }
  Rotamer* pRepresentative = NULL;
  for (int i = 0;i < pRotamerSet->repNum;i++)
  {
    if (!strcmp(RotamerGetType(pThis), RotamerGetType(&pRotamerSet->reps[i])))
    {
      pRepresentative = &pRotamerSet->reps[i];
      break;
    }
  }
  if (pRepresentative == NULL)
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d, cannot find representative rotamer for %s",
      __FILE__, __LINE__, RotamerGetType(pThis));
    TraceError(errMsg, result);
    return result;
  }
  // Restore atoms and bonds
  AtomArrayCopy(&pThis->atoms, &pRepresentative->atoms);
  BondsCopy(&pThis->bonds, &pRepresentative->bonds);
  if (AtomArrayGetNum(&pThis->atoms) != XYZArrayGetNum(&pThis->xyzs))
  {
    result = ValueError;
    sprintf(errMsg, "In file %s line %d, rotamer's atom number (%d) is not equal to representative rotamer's atom number (%d)",
      __FILE__, __LINE__, XYZArrayGetNum(&pThis->xyzs), AtomArrayGetNum(&pThis->atoms));
    TraceError(errMsg, result);
    return result;
  }
  // Copy atom XYZ from Xyzs
  for (int i = 0; i < AtomArrayGetNum(&pThis->atoms); i++)
  {
    AtomArrayGet(&pThis->atoms, i)->xyz = *XYZArrayGet(&pThis->xyzs, i);
  }

  return Success;
}

int RotamerOfProteinInitAtomsAndBonds_Charmm22(Rotamer* pThis, Residue* pResi, AtomParams* atomParams, ResTopos* resiTopos)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  for (int i = 0;i < ResidueGetAtomNum(pResi);i++)
  {
    Atom* pResiAtom = ResidueGetAtom(pResi, i);
    if (pResiAtom->isBBAtom)
    {
      AtomArrayAppend(&pThis->atoms, pResiAtom);
    }
  }

  if (strcmp(ResidueGetName(pResi), "GLY") == 0)
  {
    AtomArrayRemoveByName(&pThis->atoms, "HA1");
    AtomArrayRemoveByName(&pThis->atoms, "HA2");
  }
  else
  {
    AtomArrayRemoveByName(&pThis->atoms, "HA");
    AtomArrayRemoveByName(&pThis->atoms, "CB");
  }

  int num;
  int result = AtomParamsGetAtomNumOfResidue(atomParams, RotamerGetType(pThis), &num);
  if (FAILED(result))
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d,cannot find atom parameters for %s", __FILE__, __LINE__, RotamerGetType(pThis));
    TraceError(errMsg, result);
    return result;
  }

  Atom rotamerAtom;
  AtomCreate(&rotamerAtom);
  for (int i = 0;i < num;i++)
  {
    AtomParamsGetAtomParamOfResidue(atomParams, RotamerGetType(pThis), i, &rotamerAtom);
    if (rotamerAtom.isBBAtom == FALSE)
    {
      rotamerAtom.isXyzValid = FALSE;
      AtomArrayAppend(&pThis->atoms, &rotamerAtom);
    }
  }

  if (strcmp(ResidueGetName(pResi), "PRO") == 0 && strcmp(RotamerGetType(pThis), "PRO") != 0)
  {
    XYZ xyzCD, xyzN, xyzHN, diff;
    AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "HN", &rotamerAtom);
    ResidueGetAtomXYZ(pResi, "CD", &xyzCD);
    ResidueGetAtomXYZ(pResi, "N", &xyzN);
    diff = XYZDiff(&xyzN, &xyzCD);
    XYZScale(&diff, 1.0 / XYZDist(&xyzN, &xyzCD));
    xyzHN = XYZSum(&xyzN, &diff);
    rotamerAtom.xyz = xyzHN;
    rotamerAtom.isXyzValid = TRUE;
    AtomArrayAppend(&pThis->atoms, &rotamerAtom);
  }
  else if (strcmp(ResidueGetName(pResi), "PRO") == 0 && strcmp(RotamerGetType(pThis), "PRO") == 0)
  {
    // do nothing;
  }

  if (strcmp(RotamerGetType(pThis), "GLY") == 0)
  {
    AtomParamsGetAtomParamByResidueName(atomParams, "GLY", "HA1", &rotamerAtom);
    rotamerAtom.isXyzValid = FALSE;
    AtomArrayAppend(&pThis->atoms, &rotamerAtom);

    AtomParamsGetAtomParamByResidueName(atomParams, "GLY", "HA2", &rotamerAtom);
    rotamerAtom.isXyzValid = FALSE;
    AtomArrayAppend(&pThis->atoms, &rotamerAtom);
  }
  else
  {
    XYZ xyz;

    AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "CB", &rotamerAtom);
    // calculate the coordinate for CB;
    // if the original residue is glycine, calculate CB according to the chirality;
    // else get the coordinate from the crystal structure;
    if (strcmp(pResi->name, "GLY") == 0)
    {
      ResTopology rotamerTopo;
      CharmmIC ic;
      ResTopologyCreate(&rotamerTopo);
      CharmmICCreate(&ic);
      ResToposGet(resiTopos, RotamerGetType(pThis), &rotamerTopo);
      ResTopologyFindCharmmIC(&rotamerTopo, "CB", &ic);
      GetFourthAtom(&RotamerGetAtomByName(pThis, "N")->xyz, &RotamerGetAtomByName(pThis, "C")->xyz, &RotamerGetAtomByName(pThis, "CA")->xyz, ic.icParam, &xyz);
      ResTopologyDestroy(&rotamerTopo);
      CharmmICDestroy(&ic);
    }
    else
    {
      result = ResidueGetAtomXYZ(pResi, "CB", &xyz);
    }

    rotamerAtom.xyz = xyz;
    rotamerAtom.isXyzValid = TRUE;
    AtomArrayAppend(&pThis->atoms, &rotamerAtom);
    AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "HA", &rotamerAtom);
    rotamerAtom.isXyzValid = FALSE;
    AtomArrayAppend(&pThis->atoms, &rotamerAtom);
  }
  AtomDestroy(&rotamerAtom);

  // Deal with the bonds
  Bonds newBonds;
  BondsCreate(&newBonds);
  // Add the bonds between mainchain atoms in the original residue
  Bonds* pOriginalBondsInResidue = ResidueGetBonds(pResi);
  for (int i = 0;i < BondsGetNum(pOriginalBondsInResidue);i++)
  {
    Bond* pCurBond = BondsGet(pOriginalBondsInResidue, i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms, fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms, toAtomName);
    if (fromAtomName[0] == '+' || fromAtomName[0] == '-' || toAtomName[0] == '+' || toAtomName[0] == '-' || (fromAtom != NULL && fromAtom->isBBAtom && toAtom != NULL && toAtom->isBBAtom))
    {
      BondsAdd(&newBonds, fromAtomName, toAtomName, BondGetType(pCurBond));
    }
  }

  // Add the bonds between sidechain atoms in the rotamer
  ResTopology rotTopo;
  ResTopologyCreate(&rotTopo);
  if (FAILED(ResToposGet(resiTopos, pThis->type, &rotTopo)))
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d,cannot find residue topology for %s", __FILE__, __LINE__, RotamerGetType(pThis));
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0;i < BondsGetNum(ResTopologyGetBonds(&rotTopo));i++)
  {
    Bond* pCurBond = BondsGet(ResTopologyGetBonds(&rotTopo), i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms, fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms, toAtomName);
    if (fromAtom == NULL || fromAtom->isBBAtom)
    {
      continue;
    }
    if (toAtom == NULL || toAtom->isBBAtom)
    {
      continue;
    }
    BondsAdd(&newBonds, fromAtomName, toAtomName, BondGetType(pCurBond));
  }

  // Add the bonds between CB and other atoms
  for (int i = 0;i < BondsGetNum(ResTopologyGetBonds(&rotTopo));i++)
  {
    Bond* pCurBond = BondsGet(ResTopologyGetBonds(&rotTopo), i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms, fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms, toAtomName);
    if ((strcmp(fromAtomName, "CB") == 0 && toAtom != NULL) || (strcmp(toAtomName, "CB") == 0 && fromAtom != NULL))
    {
      BondsAdd(&newBonds, fromAtomName, toAtomName, BondGetType(pCurBond));
    }
  }
  // deal with special cases
  if (strcmp(ResidueGetName(pResi), "GLY") == 0 && strcmp(RotamerGetType(pThis), "GLY") != 0)
  {
    BondsAdd(&newBonds, "CA", "HA", Type_Bond_Single);
  }
  if (strcmp(ResidueGetName(pResi), "GLY") != 0 && strcmp(RotamerGetType(pThis), "GLY") == 0)
  {
    BondsAdd(&newBonds, "CA", "HA1", Type_Bond_Single);
    BondsAdd(&newBonds, "CA", "HA2", Type_Bond_Single);
  }
  BondsCopy(&pThis->bonds, &newBonds);
  BondsDestroy(&newBonds);

  // set the chain name and posInchain for every atom in the rotamer
  for (int i = 0;i < AtomArrayGetNum(&pThis->atoms);i++)
  {
    AtomSetChainName(AtomArrayGet(&pThis->atoms, i), ResidueGetChainName(pResi));
    AtomSetPosInChain(AtomArrayGet(&pThis->atoms, i), ResidueGetPosInChain(pResi));
  }

  XYZArrayResize(&pThis->xyzs, AtomArrayGetNum(&pThis->atoms));
  ResTopologyDestroy(&rotTopo);
  return Success;
}

int RotamerOfProteinInitAtomsAndBonds_Charmm19(Rotamer* pThis, Residue* pResi, AtomParams* atomParams, ResTopos* resiTopos)
{
  char errMsg[MAX_LEN_ERR_MSG+1];

  // copy backbone N/CA/C/O atoms;
  for (int i = 0;i < ResidueGetAtomNum(pResi);i++)
  {
    Atom* pResiAtom = ResidueGetAtom(pResi, i);
    if (pResiAtom->isBBAtom && AtomIsHydrogen(pResiAtom) == FALSE)
    {
      AtomArrayAppend(&pThis->atoms, pResiAtom);
    }
  }
  // copy sidechain atoms;
  int num;
  int result = AtomParamsGetAtomNumOfResidue(atomParams, RotamerGetType(pThis), &num);
  if (FAILED(result))
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d,cannot find atom parameters for %s", __FILE__, __LINE__, RotamerGetType(pThis));
    TraceError(errMsg, result);
    return result;
  }
  Atom rotamerAtom;
  AtomCreate(&rotamerAtom);
  for (int i = 0;i < num;i++)
  {
    AtomParamsGetAtomParamOfResidue(atomParams, RotamerGetType(pThis), i, &rotamerAtom);
    if (rotamerAtom.isBBAtom == FALSE)
    {
      rotamerAtom.isXyzValid = FALSE;
      rotamerAtom.posInChain = pResi->posInChain;
      AtomArrayAppend(&pThis->atoms, &rotamerAtom);
    }
  }

  //handle backbone hydrogen atoms
  if (pResi->terType == Type_ResIsNter)
  {//Residue is at N-ter
    if (strcmp(ResidueGetName(pResi), "GLY") == 0)
    {//Residue is Glycine
      if (strcmp(RotamerGetType(pThis), "GLY") == 0)
      {
        AtomArrayAppend(&pThis->atoms, ResidueGetAtomByName(pResi, "HT1"));
        AtomArrayAppend(&pThis->atoms, ResidueGetAtomByName(pResi, "HT2"));
        AtomArrayAppend(&pThis->atoms, ResidueGetAtomByName(pResi, "HT3"));
      }
      else if (strcmp(RotamerGetType(pThis), "PRO") == 0)
      {//Rotamer is proline, copy HT1 and HT2 from residue, but recalculate the coordinates
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "HT1"));
        rotamerAtom.isXyzValid = FALSE;
        AtomArrayAppend(&pThis->atoms, &rotamerAtom);
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "HT2"));
        rotamerAtom.isXyzValid = FALSE;
        AtomArrayAppend(&pThis->atoms, &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, "PROP", "N", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "N"), &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, "PROP", "CA", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CA"), &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, "PROP", "CD", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CD"), &rotamerAtom);
      }
      else
      {
        ResTopology nterTopo;
        CharmmIC charmm;
        StringArray hydrogens;
        StringArrayCreate(&hydrogens);
        StringArrayAppend(&hydrogens, "HT1");
        StringArrayAppend(&hydrogens, "HT2");
        StringArrayAppend(&hydrogens, "HT3");
        ResTopologyCreate(&nterTopo);
        ResToposGet(resiTopos, "NTER", &nterTopo);
        CharmmICCreate(&charmm);
        for (int i = 0; i < StringArrayGetNum(&hydrogens); i++)
        {
          char* atomName = StringArrayGet(&hydrogens, i);
          AtomParamsGetAtomParamByResidueName(atomParams, "NTER", atomName, &rotamerAtom);
          ResTopologyFindCharmmIC(&nterTopo, atomName, &charmm);
          GetFourthAtom(&RotamerGetAtomByName(pThis, charmm.atomNames[0])->xyz,
            &RotamerGetAtomByName(pThis, charmm.atomNames[1])->xyz,
            &RotamerGetAtomByName(pThis, charmm.atomNames[2])->xyz,
            charmm.icParam,
            &rotamerAtom.xyz);
          rotamerAtom.isXyzValid = TRUE;
          AtomArrayAppend(&pThis->atoms, &rotamerAtom);
        }
        StringArrayDestroy(&hydrogens);
        CharmmICDestroy(&charmm);
        ResTopologyDestroy(&nterTopo);
        AtomParamsGetAtomParamByResidueName(atomParams, "NTER", "N", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "N"), &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, "NTER", "CA", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CA"), &rotamerAtom);
      }
    }
    else if (strcmp(ResidueGetName(pResi), "PRO") == 0)
    {//Residue is Proline
      if (strcmp(RotamerGetType(pThis), "GLY") == 0)
      {//Rotamer is Glycine
        ResTopology nterTopo;
        CharmmIC charmm;
        StringArray hydrogens;
        StringArrayCreate(&hydrogens);
        StringArrayAppend(&hydrogens, "HT1");
        StringArrayAppend(&hydrogens, "HT2");
        StringArrayAppend(&hydrogens, "HT3");
        ResTopologyCreate(&nterTopo);
        ResToposGet(resiTopos, "GLYP", &nterTopo);
        CharmmICCreate(&charmm);
        for (int i = 0; i < StringArrayGetNum(&hydrogens); i++)
        {
          char* atomName = StringArrayGet(&hydrogens, i);
          AtomParamsGetAtomParamByResidueName(atomParams, "GLYP", atomName, &rotamerAtom);
          ResTopologyFindCharmmIC(&nterTopo, atomName, &charmm);
          GetFourthAtom(&RotamerGetAtomByName(pThis, charmm.atomNames[0])->xyz,
            &RotamerGetAtomByName(pThis, charmm.atomNames[1])->xyz,
            &RotamerGetAtomByName(pThis, charmm.atomNames[2])->xyz,
            charmm.icParam,
            &rotamerAtom.xyz);
          rotamerAtom.isXyzValid = TRUE;
          AtomArrayAppend(&pThis->atoms, &rotamerAtom);
        }
        StringArrayDestroy(&hydrogens);
        CharmmICDestroy(&charmm);
        ResTopologyDestroy(&nterTopo);
        AtomParamsGetAtomParamByResidueName(atomParams, "GLYP", "N", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "N"), &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, "GLYP", "CA", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CA"), &rotamerAtom);
      }
      else if (strcmp(RotamerGetType(pThis), "PRO") == 0)
      {//Rotamer is Proline; cannot use PROP because the coordinate of CD has not been determined
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "HT1"));
        rotamerAtom.isXyzValid = FALSE;
        AtomArrayAppend(&pThis->atoms, &rotamerAtom);
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "HT2"));
        rotamerAtom.isXyzValid = FALSE;
        AtomArrayAppend(&pThis->atoms, &rotamerAtom);
      }
      else
      {//Rotamer is in any other type
        ResTopology nterTopo;
        CharmmIC charmm;
        StringArray hydrogens;
        StringArrayCreate(&hydrogens);
        StringArrayAppend(&hydrogens, "HT1");
        StringArrayAppend(&hydrogens, "HT2");
        StringArrayAppend(&hydrogens, "HT3");
        ResTopologyCreate(&nterTopo);
        ResToposGet(resiTopos, "NTER", &nterTopo);
        CharmmICCreate(&charmm);
        for (int i = 0; i < StringArrayGetNum(&hydrogens); i++)
        {
          char* atomName = StringArrayGet(&hydrogens, i);
          AtomParamsGetAtomParamByResidueName(atomParams, "NTER", atomName, &rotamerAtom);
          ResTopologyFindCharmmIC(&nterTopo, atomName, &charmm);
          GetFourthAtom(&RotamerGetAtomByName(pThis, charmm.atomNames[0])->xyz,
            &RotamerGetAtomByName(pThis, charmm.atomNames[1])->xyz,
            &RotamerGetAtomByName(pThis, charmm.atomNames[2])->xyz,
            charmm.icParam,
            &rotamerAtom.xyz);
          rotamerAtom.isXyzValid = TRUE;
          AtomArrayAppend(&pThis->atoms, &rotamerAtom);
        }
        StringArrayDestroy(&hydrogens);
        CharmmICDestroy(&charmm);
        ResTopologyDestroy(&nterTopo);
        AtomParamsGetAtomParamByResidueName(atomParams, "NTER", "N", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "N"), &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, "NTER", "CA", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CA"), &rotamerAtom);
      }
    }
    else
    {//Residue is in other type (not Gly and Pro)
      if (strcmp(RotamerGetType(pThis), "GLY") == 0)
      {
        ResTopology nterTopo;
        CharmmIC charmm;
        StringArray hydrogens;
        StringArrayCreate(&hydrogens);
        StringArrayAppend(&hydrogens, "HT1");
        StringArrayAppend(&hydrogens, "HT2");
        StringArrayAppend(&hydrogens, "HT3");
        ResTopologyCreate(&nterTopo);
        ResToposGet(resiTopos, "GLYP", &nterTopo);
        CharmmICCreate(&charmm);
        for (int i = 0; i < StringArrayGetNum(&hydrogens); i++)
        {
          char* atomName = StringArrayGet(&hydrogens, i);
          AtomParamsGetAtomParamByResidueName(atomParams, "GLYP", atomName, &rotamerAtom);
          ResTopologyFindCharmmIC(&nterTopo, atomName, &charmm);
          GetFourthAtom(&RotamerGetAtomByName(pThis, charmm.atomNames[0])->xyz,
            &RotamerGetAtomByName(pThis, charmm.atomNames[1])->xyz,
            &RotamerGetAtomByName(pThis, charmm.atomNames[2])->xyz,
            charmm.icParam,
            &rotamerAtom.xyz);
          rotamerAtom.isXyzValid = TRUE;
          AtomArrayAppend(&pThis->atoms, &rotamerAtom);
        }
        StringArrayDestroy(&hydrogens);
        CharmmICDestroy(&charmm);
        ResTopologyDestroy(&nterTopo);
        AtomParamsGetAtomParamByResidueName(atomParams, "GLYP", "N", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "N"), &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, "GLYP", "CA", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CA"), &rotamerAtom);
      }
      else if (strcmp(RotamerGetType(pThis), "PRO") == 0)
      {
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "HT1"));
        rotamerAtom.isXyzValid = FALSE;
        AtomArrayAppend(&pThis->atoms, &rotamerAtom);
        AtomCopy(&rotamerAtom, ResidueGetAtomByName(pResi, "HT2"));
        rotamerAtom.isXyzValid = FALSE;
        AtomArrayAppend(&pThis->atoms, &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, "PROP", "N", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "N"), &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, "PROP", "CA", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CA"), &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, "PROP", "CD", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CD"), &rotamerAtom);
      }
      else
      {
        AtomArrayAppend(&pThis->atoms, ResidueGetAtomByName(pResi, "HT1"));
        AtomArrayAppend(&pThis->atoms, ResidueGetAtomByName(pResi, "HT2"));
        AtomArrayAppend(&pThis->atoms, ResidueGetAtomByName(pResi, "HT3"));
        /*ResidueTopology nterTopo;
        CharmmIC charmm;
        StringArray hydrogens;
        StringArrayCreate(&hydrogens);
        StringArrayAppend(&hydrogens, "HT1");
        StringArrayAppend(&hydrogens, "HT2");
        StringArrayAppend(&hydrogens, "HT3");
        ResidueTopologyCreate(&nterTopo);
        ResiTopoSetGet(resiTopos, "NTER", &nterTopo);
        CharmmICCreate(&charmm);
        for(int i = 0; i < StringArrayGetNum(&hydrogens); i++){
          char *atName = StringArrayGet(&hydrogens, i);
          AtomParamsSetGetAtomParamByName(atomParams, "NTER", atName, &rotamerAtom);
          ResidueTopologyFindCharmmIC(&nterTopo, atName, &charmm);
          GetFourthAtom(&RotamerGetAtomByName(pRot, charmm.atomNames[0])->xyz,
            &RotamerGetAtomByName(pRot, charmm.atomNames[1])->xyz,
            &RotamerGetAtomByName(pRot, charmm.atomNames[2])->xyz,
            charmm.icParam,
            &rotamerAtom.xyz);
          rotamerAtom.isXyzValid = TRUE;
          AtomArrayAppend(&pRot->atoms,&rotamerAtom);
        }
        StringArrayDestroy(&hydrogens);
        CharmmICDestroy(&charmm);
        ResidueTopologyDestroy(&nterTopo);*/
      }
    }
  }
  else if (pResi->terType != Type_ResIsNter)
  {//Residue is not at N-ter
    if (strcmp(ResidueGetName(pResi), "GLY") == 0)
    {//Residue is Glycine
      if (strcmp(RotamerGetType(pThis), "GLY") == 0)
      {
        AtomArrayAppend(&pThis->atoms, ResidueGetAtomByName(pResi, "H"));
      }
      else if (strcmp(RotamerGetType(pThis), "PRO") == 0)
      {
        AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "N", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "N"), &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "CA", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CA"), &rotamerAtom);
      }
      else
      {
        AtomArrayAppend(&pThis->atoms, ResidueGetAtomByName(pResi, "H"));
        AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "CA", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CA"), &rotamerAtom);
      }
    }
    else if (strcmp(ResidueGetName(pResi), "PRO") == 0)
    {//Residue is Proline
      if (strcmp(RotamerGetType(pThis), "PRO") != 0)
      {
        XYZ xyzCD, xyzN, xyzHN, diff;
        AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "H", &rotamerAtom);
        ResidueGetAtomXYZ(pResi, "CD", &xyzCD);
        ResidueGetAtomXYZ(pResi, "N", &xyzN);
        diff = XYZDiff(&xyzN, &xyzCD);
        XYZScale(&diff, 1.0 / XYZDist(&xyzN, &xyzCD));
        xyzHN = XYZSum(&xyzN, &diff);
        rotamerAtom.xyz = xyzHN;
        rotamerAtom.isXyzValid = TRUE;
        AtomArrayAppend(&pThis->atoms, &rotamerAtom);
        //change the parameters of Pro N
        AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "N", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "N"), &rotamerAtom);
        if (strcmp(RotamerGetType(pThis), "GLY") == 0)
        {
          AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "CA", &rotamerAtom);
          AtomCopyParameter(RotamerGetAtomByName(pThis, "CA"), &rotamerAtom);
        }
      }
    }
    else
    {
      if (strcmp(RotamerGetType(pThis), "GLY") == 0)
      {
        AtomArrayAppend(&pThis->atoms, ResidueGetAtomByName(pResi, "H"));
        AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "CA", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CA"), &rotamerAtom);
      }
      else if (strcmp(RotamerGetType(pThis), "PRO") == 0)
      {
        AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "N", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "N"), &rotamerAtom);
        AtomParamsGetAtomParamByResidueName(atomParams, RotamerGetType(pThis), "CA", &rotamerAtom);
        AtomCopyParameter(RotamerGetAtomByName(pThis, "CA"), &rotamerAtom);
      }
      else
      {
        AtomArrayAppend(&pThis->atoms, ResidueGetAtomByName(pResi, "H"));
      }
    }
  }
  AtomDestroy(&rotamerAtom);

  //deal with the CB atoms
  if (strcmp(ResidueGetName(pResi), "GLY") == 0)
  {
    if (strcmp(RotamerGetType(pThis), "GLY") != 0)
    {
      ResTopology rotamerTopo;
      CharmmIC ic;

      ResTopologyCreate(&rotamerTopo);
      CharmmICCreate(&ic);
      ResToposGet(resiTopos, RotamerGetType(pThis), &rotamerTopo);
      ResTopologyFindCharmmIC(&rotamerTopo, "CB", &ic);
      GetFourthAtom(&RotamerGetAtomByName(pThis, "N")->xyz, &RotamerGetAtomByName(pThis, "C")->xyz, &RotamerGetAtomByName(pThis, "CA")->xyz, ic.icParam, &AtomArrayGetByName(&pThis->atoms, "CB")->xyz);
      AtomArrayGetByName(&pThis->atoms, "CB")->isXyzValid = TRUE;
      ResTopologyDestroy(&rotamerTopo);
      CharmmICDestroy(&ic);
    }
    else
    {
      // do nothing;
    }
  }
  else
  {//both residue and rotamer are not glycine
    if (strcmp(RotamerGetType(pThis), "GLY") != 0)
    {
      Atom* pAtomCB = AtomArrayGetByName(&pThis->atoms, "CB");
      pAtomCB->xyz = ResidueGetAtomByName(pResi, "CB")->xyz;
      pAtomCB->isXyzValid = TRUE;
    }
    //else{//residue is not glycine, rotamer is glycine
    //  AtomParamsSetGetAtomParamByName(atomParams,RotamerGetType(pRot),"CA",&rotamerAtom);
    //  AtomCopyParameter(RotamerGetAtomByName(pRot,"CA"),&rotamerAtom);
    //}
  }

  // Deal with the bonds
  Bonds newBonds;
  BondsCreate(&newBonds);
  // Add the bonds between mainchain atoms in the original residue
  Bonds* pOriginalBondsInResidue = ResidueGetBonds(pResi);
  for (int i = 0;i < BondsGetNum(pOriginalBondsInResidue);i++)
  {
    Bond* pCurBond = BondsGet(pOriginalBondsInResidue, i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms, fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms, toAtomName);
    if (fromAtomName[0] == '+' || fromAtomName[0] == '-' || toAtomName[0] == '+' || toAtomName[0] == '-' || (fromAtom != NULL && fromAtom->isBBAtom && toAtom != NULL && toAtom->isBBAtom))
    {
      BondsAdd(&newBonds, fromAtomName, toAtomName, BondGetType(pCurBond));
    }
  }

  // Add the bonds between sidechain atoms in the rotamer residue
  ResTopology rotTopo;
  ResTopologyCreate(&rotTopo);
  if (FAILED(ResToposGet(resiTopos, pThis->type, &rotTopo)))
  {
    // Cannot Find this type of rotamer in Residue Topologies
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d,cannot find residue topology for %s", __FILE__, __LINE__, RotamerGetType(pThis));
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0;i < BondsGetNum(ResTopologyGetBonds(&rotTopo));i++)
  {
    Bond* pCurBond = BondsGet(ResTopologyGetBonds(&rotTopo), i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms, fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms, toAtomName);
    if (fromAtom == NULL || fromAtom->isBBAtom)
    {
      continue;
    }
    if (toAtom == NULL || toAtom->isBBAtom)
    {
      continue;
    }
    BondsAdd(&newBonds, fromAtomName, toAtomName, BondGetType(pCurBond));
  }

  // Add the bonds between CB and other atoms
  for (int i = 0;i < BondsGetNum(ResTopologyGetBonds(&rotTopo));i++)
  {
    Bond* pCurBond = BondsGet(ResTopologyGetBonds(&rotTopo), i);
    char* fromAtomName = BondGetFromName(pCurBond);
    char* toAtomName = BondGetToName(pCurBond);
    Atom* fromAtom = AtomArrayGetByName(&pThis->atoms, fromAtomName);
    Atom* toAtom = AtomArrayGetByName(&pThis->atoms, toAtomName);
    if ((strcmp(fromAtomName, "CB") == 0 && toAtom != NULL) || (strcmp(toAtomName, "CB") == 0 && fromAtom != NULL))
    {
      BondsAdd(&newBonds, fromAtomName, toAtomName, BondGetType(pCurBond));
    }
  }

  // deal with special cases for prolines
  if (strcmp(ResidueGetName(pResi), "PRO") == 0)
  {
    if (pResi->terType == Type_ResIsNter)
    {
      if (strcmp(RotamerGetType(pThis), "PRO") == 0)
      {
        BondsAdd(&newBonds, "N", "CD", Type_Bond_Single);
      }
      else
      {
        BondsAdd(&newBonds, "HT3", "N", Type_Bond_Single);
      }
    }
    else
    {
      if (strcmp(RotamerGetType(pThis), "PRO") == 0)
      {
        BondsAdd(&newBonds, "N", "CD", Type_Bond_Single);
      }
      else
      {
        BondsAdd(&newBonds, "N", "H", Type_Bond_Single);
      }
    }
  }
  else
  {
    if (strcmp(RotamerGetType(pThis), "PRO") == 0)
    {
      BondsAdd(&newBonds, "N", "CD", Type_Bond_Single);
    }
  }

  // copy bonds back to rotamer
  BondsCopy(&pThis->bonds, &newBonds);
  BondsDestroy(&newBonds);


  // set the chain name and posInChain for every atom in the residue
  for (int i = 0;i < AtomArrayGetNum(&pThis->atoms);i++)
  {
    AtomSetChainName(AtomArrayGet(&pThis->atoms, i), ResidueGetChainName(pResi));
    AtomSetPosInChain(AtomArrayGet(&pThis->atoms, i), ResidueGetPosInChain(pResi));
  }

  XYZArrayResize(&pThis->xyzs, AtomArrayGetNum(&pThis->atoms));
  ResTopologyDestroy(&rotTopo);
  return Success;
}

int RotamerOfProteinPatch(Rotamer* pThis, char* patchType, AtomParams* atomParams, ResTopos* resiTopo)
{
  // Copy atoms and bonds to a temporary residue
  Residue tempResi;
  ResidueCreate(&tempResi);
  ResidueSetName(&tempResi, pThis->type);
  for (int i = 0;i < AtomArrayGetNum(&pThis->atoms);i++)
  {
    ResidueAddAtom(&tempResi, AtomArrayGet(&pThis->atoms, i));
  }
  BondsCopy(ResidueGetBonds(&tempResi), &pThis->bonds);

  // Call the ResiduePatch to patch the temporary residue
  int result = ResiduePatch(&tempResi, patchType, atomParams, resiTopo);
  if (FAILED(result))
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    sprintf(errMsg, "in file %s line %d, patching %s on rotamer %s failed", __FILE__, __LINE__, patchType, pThis->type);
    TraceError(errMsg, result);
    return result;
  }

  // Copy back atoms and bonds
  AtomArrayCopy(&pThis->atoms, ResidueGetAllAtoms(&tempResi));
  BondsCopy(&pThis->bonds, ResidueGetBonds(&tempResi));

  XYZArrayResize(&pThis->xyzs, AtomArrayGetNum(&pThis->atoms));
  for (int i = 0;i < AtomArrayGetNum(&pThis->atoms);i++)
  {
    XYZArraySet(&pThis->xyzs, i, &(AtomArrayGet(&pThis->atoms, i)->xyz));
  }

  ResidueDestroy(&tempResi);
  return Success;
}

int RotamerOfProteinCalcXYZ(Rotamer* pThis, Residue* pResi, char* patchName, DoubleArray* torsions, ResTopos* resiTopos)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  ResTopology rotamerTopology;
  ResTopologyCreate(&rotamerTopology);
  int result = ResToposGet(resiTopos, pThis->type, &rotamerTopology);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, cannot find the Topology for rotamer %s", __FILE__, __LINE__, RotamerGetType(pThis));
    TraceError(errMsg, result);
    return result;
  }

  int torNum = DoubleArrayGetLength(torsions);
  torNum = torNum > 5 ? 5 : torNum;
  // Calculate the XYZ of atoms which are directly determined by the rotamer's torsion angles
  for (int torsionIndex = 0;torsionIndex < torNum;torsionIndex++)
  {
    Type_ProteinAtomOrder desiredAtomBOrder = Type_ProteinAtomOrder_FromInt(torsionIndex);
    Type_ProteinAtomOrder desiredAtomCOrder = Type_ProteinAtomOrder_FromInt(torsionIndex + 1);
    CharmmIC icOfCurrentTorsion;
    CharmmICCreate(&icOfCurrentTorsion);
    BOOL icFound = FALSE;
    for (int icIndex = 0; icIndex < ResTopologyGetCharmmICNum(&rotamerTopology); icIndex++)
    {
      Type_ProteinAtomOrder atomBOrder;
      Type_ProteinAtomOrder atomCOrder;
      ResTopologyGetCharmmIC(&rotamerTopology, icIndex, &icOfCurrentTorsion);
      atomBOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomB(&icOfCurrentTorsion));
      atomCOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomC(&icOfCurrentTorsion));
      if (desiredAtomBOrder == atomBOrder && desiredAtomCOrder == atomCOrder)
      {
        icFound = TRUE;
        break;
      }
    }

    if (!icFound)
    {
      result = DataNotExistError;
      sprintf(errMsg, "In file %s line %d, cannot find IC of the %dth torsion for rotamer %s", __FILE__, __LINE__, torsionIndex + 1, pThis->type);
      TraceError(errMsg, result);
      CharmmICDestroy(&icOfCurrentTorsion);
      return result;
    }

    double torsionValue = DoubleArrayGet(torsions, torsionIndex);
    // correct the torsion X2 for trp rotamer;
    //if(strcmp(RotamerGetType(pRot), "TRP") == 0 && torsionIndex == 1){
    //  torsionValue += PI;
    //  if(torsionValue > PI) torsionValue -= 2.0*PI;
    //}

    CharmmICSetTorsionAngle(&icOfCurrentTorsion, torsionValue);
    Atom* pAtomD = AtomArrayGetByName(&pThis->atoms, CharmmICGetAtomD(&icOfCurrentTorsion));
    XYZ newXYZofAtomD;
    if (pAtomD == NULL || FAILED(CharmmICCalcXYZ(&icOfCurrentTorsion, &pThis->atoms, &newXYZofAtomD)))
    {
      result = DataNotExistError;
      sprintf(errMsg, "In file %s line %d, cannot calculate coordinate for atom %s", __FILE__, __LINE__, AtomGetName(pAtomD));
      TraceError(errMsg, result);
      CharmmICShow(&icOfCurrentTorsion);printf("\n");
      CharmmICDestroy(&icOfCurrentTorsion);
      return result;
    }
    pAtomD->xyz = newXYZofAtomD;
    pAtomD->isXyzValid = TRUE;
    CharmmICDestroy(&icOfCurrentTorsion);
  }

  // calculate backbone atoms with invalid coordinates;
  if (strcmp(RotamerGetType(pThis), "PRO") == 0)
  {
    if (pResi->terType == Type_ResIsNter)
    {
      double tempIC[5];
      Atom* pAtomCA = RotamerGetAtomByName(pThis, "CA");
      Atom* pAtomCD = RotamerGetAtomByName(pThis, "CD");
      Atom* pAtomN = RotamerGetAtomByName(pThis, "N");
      Atom* pAtomHT1 = RotamerGetAtomByName(pThis, "HT1");
      Atom* pAtomHT2 = RotamerGetAtomByName(pThis, "HT2");
      tempIC[0] = 0.0; tempIC[1] = 0.0; tempIC[2] = DegToRad(120.0); tempIC[3] = DegToRad(109.5); tempIC[4] = 1.04;
      GetFourthAtom(&pAtomCA->xyz, &pAtomCD->xyz, &pAtomN->xyz, tempIC, &pAtomHT1->xyz);
      pAtomHT1->isXyzValid = TRUE;
      tempIC[2] = DegToRad(-120.0);
      GetFourthAtom(&pAtomCA->xyz, &pAtomCD->xyz, &pAtomN->xyz, tempIC, &pAtomHT2->xyz);
      pAtomHT2->isXyzValid = TRUE;
    }
  }

  // Calculate the XYZ of other side chain atoms
  if (!AtomArrayAllAtomXYZAreValid(&pThis->atoms))
  {
    Residue tempResi;
    ResidueCreate(&tempResi);
    ResidueSetName(&tempResi, pThis->type);
    if (patchName != NULL && strcmp(patchName, "") != 0)
    {
      StringArrayAppend(&tempResi.patches, patchName);
    }
    for (int i = 0;i < AtomArrayGetNum(&pThis->atoms);i++)
    {
      ResidueAddAtom(&tempResi, AtomArrayGet(&pThis->atoms, i));
    }
    result = ResidueCalcAllAtomXYZ(&tempResi, resiTopos, NULL, NULL);
    if (FAILED(result))
    {
      result = DataNotExistError;
      sprintf(errMsg, "in file %s line %d, some atoms' XYZ cannot be calculated", __FILE__, __LINE__);
      TraceError(errMsg, result);
      AtomArrayShowInPDBFormat(&pThis->atoms, "ATOM", pThis->type, " ", 0, 0, NULL);
      return result;
    }
    for (int i = 0;i < AtomArrayGetNum(&pThis->atoms);i++)
    {
      AtomCopy(AtomArrayGet(&pThis->atoms, i), ResidueGetAtom(&tempResi, i));
    }
    ResidueDestroy(&tempResi);
  }

  for (int i = 0;i < AtomArrayGetNum(&pThis->atoms);i++)
  {
    XYZArraySet(&pThis->xyzs, i, &(AtomArrayGet(&pThis->atoms, i)->xyz));
  }

  ResTopologyDestroy(&rotamerTopology);
  return Success;
}


int RotamerOfProteinGenerate(Rotamer* pThis, Residue* pResi, char* rotamerType, char* patchType, DoubleArray* torsions, AtomParams* atomParams, ResTopos* resiTopos)
{
  if (strlen(rotamerType) > MAX_LEN_RES_NAME) return NameError;
  RotamerDestroy(pThis);
  RotamerCreate(pThis);
  strcpy(pThis->type, rotamerType);

  //First step, initialize atoms and bonds
  int result = RotamerOfProteinInitAtomsAndBonds_Charmm19(pThis, pResi, atomParams, resiTopos);
  if (FAILED(result))
  {
    return result;
  }

  //Second step, patch the rotamer if necessary
  if (patchType != NULL && strcmp(patchType, "") != 0)
  {
    result = RotamerOfProteinPatch(pThis, patchType, atomParams, resiTopos);
    if (FAILED(result))
    {
      return result;
    }
  }

  //Third step, calculate coordinates of all atoms
  result = RotamerOfProteinCalcXYZ(pThis, pResi, patchType, torsions, resiTopos);
  if (FAILED(result))
  {
    return result;
  }

  //Fourth step, set the chain name and posInChain to be consistent with the residue
  RotamerSetChainName(pThis, ResidueGetChainName(pResi));
  RotamerSetPosInChain(pThis, ResidueGetPosInChain(pResi));
  return Success;
}


int RotamerShowInPDBFormat(Rotamer* pThis, char* header, char* chainName, int atomIndex, int resiIndex, FILE* pFile)
{
  char rotType[MAX_LEN_RES_NAME+1];
  strcpy(rotType, RotamerGetType(pThis));
  if (strcmp(rotType, "HSD") == 0 || strcmp(rotType, "HSE") == 0 || strcmp(rotType, "HSP") == 0)
  {
    strcpy(rotType, "HIS");
  }
  return AtomArrayShowInPDBFormat(&pThis->atoms, header, rotType, chainName, atomIndex, resiIndex, pFile);
}

int RotamerShowAtomParameter(Rotamer* pThis)
{
  for (int i = 0; i < RotamerGetAtomNum(pThis); ++i)
  {
    Atom* pAtom = RotamerGetAtom(pThis, i);
    AtomShowAtomParameter(pAtom);
  }
  return Success;
}

int RotamerShowBondInformation(Rotamer* pThis)
{
  BondsShow(RotamerGetBonds(pThis));
  return Success;
}

int RotamersCreate(Rotamers* pThis)
{
  pThis->capacity = 2;
  pThis->num = 0;
  pThis->rots = (Rotamer*)malloc(pThis->capacity * sizeof(Rotamer));
  for (int i = 0;i < pThis->capacity;i++)
  {
    RotamerCreate(&pThis->rots[i]);
  }
  pThis->repNum = 0;
  pThis->reps = NULL;
  return Success;
}


int RotamersDestroy(Rotamers* pThis)
{
  for (int i = 0;i < pThis->capacity;i++)
  {
    RotamerDestroy(&pThis->rots[i]);
  }
  free(pThis->rots);
  pThis->rots = NULL;
  pThis->capacity = pThis->num = 0;
  for (int i = 0;i < pThis->repNum;i++)
  {
    RotamerDestroy(&pThis->reps[i]);
  }
  free(pThis->reps);
  pThis->reps = NULL;
  pThis->repNum = 0;
  return Success;
}


int RotamersCopy(Rotamers* pThis, Rotamers* pOther)
{
  RotamersDestroy(pThis);
  RotamersCreate(pThis);
  for (int i = 0;i < pThis->capacity;i++)
  {
    RotamerDestroy(&pThis->rots[i]);
  }
  free(pThis->rots);
  pThis->rots = (Rotamer*)malloc(sizeof(Rotamer) * pOther->capacity);
  for (int i = 0; i < pOther->capacity; i++)
  {
    RotamerCreate(&pThis->rots[i]);
  }
  for (int i = 0;i < pOther->num;i++)
  {
    RotamerCopy(&pThis->rots[i], &pOther->rots[i]);
  }
  pThis->capacity = pOther->capacity;
  pThis->num = pOther->num;
  pThis->reps = (Rotamer*)malloc(sizeof(Rotamer) * pOther->repNum);
  for (int i = 0;i < pOther->repNum;i++)
  {
    RotamerCreate(&pThis->reps[i]);
    RotamerCopy(&pThis->reps[i], &pOther->reps[i]);
  }
  pThis->repNum = pOther->repNum;
  return Success;
}


int RotamersGetNum(Rotamers* pThis)
{
  return pThis->num;
}


Rotamer* RotamersGet(Rotamers* pThis, int index)
{
  if (index < 0 || index >= pThis->num)
  {
    return NULL;
  }
  return &pThis->rots[index];
}


Rotamer* RotamersFindRep(Rotamers* pThis, char* type)
{
  for (int i = 0;i < pThis->repNum;i++)
  {
    if (strcmp(pThis->reps[i].type, type) == 0)
    {
      return &(pThis->reps[i]);
    }
  }
  return NULL;
}


int RotamersGetRepNum(Rotamers* pThis)
{
  return pThis->repNum;
}


Rotamer* RotamersGetRepByIndex(Rotamers* pThis, int ndx)
{
  if (ndx<0 || ndx>RotamersGetRepNum(pThis))
  {
    return NULL;
  }
  return &(pThis->reps[ndx]);
}


int RotamersAdd(Rotamers* pThis, Rotamer* pNewRotamer)
{
  // For the newly added rotamer, record only the atom coordinates
  //Rotamer* pNewlyAddedRotInTheSet = &pThis->rotamers[pThis->num];
  //RotamerCopy(pNewlyAddedRotInTheSet,pNewRotamer);
  //AtomArrayDestroy(&pNewlyAddedRotInTheSet->atoms);
  //AtomArrayCreate(&pNewlyAddedRotInTheSet->atoms);
  //BondSetDestroy(&pNewlyAddedRotInTheSet->bonds);
  //BondSetCreate(&pNewlyAddedRotInTheSet->bonds);

  //another method the create extracted rotamer
  Rotamer* pNewlyAddedRotInTheSet = &pThis->rots[pThis->num];
  strcpy(pNewlyAddedRotInTheSet->type, pNewRotamer->type);
  strcpy(pNewlyAddedRotInTheSet->chainName, pNewRotamer->chainName);
  pNewlyAddedRotInTheSet->posInChain = pNewRotamer->posInChain;
  XYZArrayCopy(&pNewlyAddedRotInTheSet->xyzs, &pNewRotamer->xyzs);
  pNewlyAddedRotInTheSet->vdwInternal = pNewRotamer->vdwInternal;
  pNewlyAddedRotInTheSet->vdwBackbone = pNewRotamer->vdwBackbone;
  pNewlyAddedRotInTheSet->selfEnergy = pNewRotamer->selfEnergy;
  pNewlyAddedRotInTheSet->dunbrack = pNewRotamer->dunbrack;
  DoubleArrayCopy(&pNewlyAddedRotInTheSet->Xs, &pNewRotamer->Xs);

  (pThis->num)++;
  if (pThis->num == pThis->capacity)
  {
    pThis->rots = (Rotamer*)realloc(pThis->rots, sizeof(Rotamer) * pThis->capacity * 2);
    for (int i = pThis->capacity; i < pThis->capacity * 2; i++)
    {
      RotamerCreate(&pThis->rots[i]);
    }
    pThis->capacity *= 2;
  }
  RotamerCopy(&pThis->rots[pThis->num - 1], pNewRotamer);

  if (RotamersFindRep(pThis, pNewRotamer->type) == NULL)
  {
    (pThis->repNum)++;
    pThis->reps = (Rotamer*)realloc(pThis->reps, sizeof(Rotamer) * pThis->repNum);
    RotamerCreate(&pThis->reps[pThis->repNum - 1]);
    RotamerCopy(&pThis->reps[pThis->repNum - 1], pNewRotamer);
  }

  return Success;
}


int RotamersShow(Rotamers* pThis, FILE* pFile)
{
  for (int i = 0;i < pThis->repNum;i++)
  {
    RotamerShow(&pThis->reps[i]);
  }
  return Success;
}


//////////////////////////////////////////////////////////////////////////
//the followings are backbone-dependent rotamer libs
//////////////////////////////////////////////////////////////////////////
int RotLibPhiPsiCreate(RotLibPhiPsi* pThis)
{
  pThis->phi = 0;
  pThis->psi = 0;
  pThis->prob = NULL;
  StringArrayCreate(&pThis->rotTypes);
  IntArrayCreate(&pThis->rotNums, 0);
  pThis->tors = NULL;
  pThis->devs = NULL;
  return Success;
}


int RotLibPhiPsiGetNum(RotLibPhiPsi* pThis, char* typeName)
{
  int resiIndex;
  if (FAILED(StringArrayFind(&pThis->rotTypes, typeName, &resiIndex)))
  {
    return Success;
  }
  return IntArrayGet(&pThis->rotNums, resiIndex);
}


int RotLibPhiPsiGet(RotLibPhiPsi* pThis, char* typeName, int rotIndex, DoubleArray* pDestTorsion, double* probability)
{
  int rotTypeIndex;
  if (FAILED(StringArrayFind(&pThis->rotTypes, typeName, &rotTypeIndex)))
  {
    return DataNotExistError;
  }
  int num = IntArrayGet(&pThis->rotNums, rotTypeIndex);
  if (rotIndex < 0 || rotIndex >= num)
  {
    return IndexError;
  }

  DoubleArrayCopy(pDestTorsion, &pThis->tors[rotTypeIndex][rotIndex]);
  *probability = DoubleArrayGet(&pThis->prob[rotTypeIndex], rotIndex);
  return Success;
}


int BBdepRotamerLibCreate(BBdepRotLib* pRotLib, char* rotlibfile)
{
  // set all the possible rotamer types
  pRotLib->num = 1296;
  pRotLib->pps = (RotLibPhiPsi*)malloc(sizeof(RotLibPhiPsi) * pRotLib->num);
  for (int i = 0;i < pRotLib->num;i++)
  {
    RotLibPhiPsiCreate(&pRotLib->pps[i]);
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "ALA");//index: 0
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "ARG");//index: 1
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "ASN");//index: 2
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "ASP");//index: 3
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "CYS");//index: 4
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "GLN");//index: 5
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "GLU");//index: 6
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "GLY");//index: 7
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "HSD");//index: 8
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "HSE");//index: 9
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "ILE");//index:10
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "LEU");//index:11
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "LYS");//index:12
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "MET");//index:13
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "PHE");//index:14
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "PRO");//index:15
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "SER");//index:16
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "THR");//index:17
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "TRP");//index:18
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "TYR");//index:19
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "VAL");//index:20
    IntArrayResize(&pRotLib->pps[i].rotNums, StringArrayGetNum(&pRotLib->pps[i].rotTypes));
    pRotLib->pps[i].tors = (DoubleArray**)malloc(sizeof(DoubleArray*) * StringArrayGetNum(&pRotLib->pps[i].rotTypes));
    pRotLib->pps[i].devs = (DoubleArray**)malloc(sizeof(DoubleArray*) * StringArrayGetNum(&pRotLib->pps[i].rotTypes));
    pRotLib->pps[i].prob = (DoubleArray*)malloc(sizeof(DoubleArray) * StringArrayGetNum(&pRotLib->pps[i].rotTypes));
    for (int j = 0;j < StringArrayGetNum(&pRotLib->pps[i].rotTypes);j++)
    {
      pRotLib->pps[i].tors[j] = NULL;
      pRotLib->pps[i].devs[j] = NULL;
      DoubleArrayCreate(&pRotLib->pps[i].prob[j], 0);
    }
  }

  FILE* pIn = fopen(rotlibfile, "r");
  char buffer[MAX_LEN_ONE_LINE_CONTENT+1];
  while (fgets(buffer, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    if (buffer[0] == ' ' || buffer[0] == '#') continue;
    int phi, psi;
    char resName[MAX_LEN_RES_NAME+1];
    double prob;
    int num, r[4];
    double x[4];
    double s[4];
    sscanf(buffer, "%s %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
      resName, &phi, &psi, &num, &r[0], &r[1], &r[2], &r[3], &prob, &x[0], &x[1], &x[2], &x[3], &s[0], &s[1], &s[2], &s[3]);
    if (phi == 180 || psi == 180) continue;
    int binNdx = ((phi + 180) / 10) * 36 + (psi + 180) / 10;
    if (binNdx >= pRotLib->num) continue;
    //if(prob<ROT_PROB_CUT_MIN) continue;
    pRotLib->pps[binNdx].phi = phi;
    pRotLib->pps[binNdx].psi = psi;
    int rtNdx = -1;
    int nx = 0;
    if (!strcmp(resName, "ALA")) { rtNdx = 0; nx = 0; }
    else if (!strcmp(resName, "ARG")) { rtNdx = 1; nx = 4; }
    else if (!strcmp(resName, "ASN")) { rtNdx = 2; nx = 2; }
    else if (!strcmp(resName, "ASP")) { rtNdx = 3; nx = 2; }
    else if (!strcmp(resName, "CYS")) { rtNdx = 4; nx = 1; }
    else if (!strcmp(resName, "GLN")) { rtNdx = 5; nx = 3; }
    else if (!strcmp(resName, "GLU")) { rtNdx = 6; nx = 3; }
    else if (!strcmp(resName, "GLY")) { rtNdx = 7; nx = 0; }
    else if (!strcmp(resName, "HSD")) { rtNdx = 8; nx = 2; }
    else if (!strcmp(resName, "HSE")) { rtNdx = 9; nx = 2; }
    else if (!strcmp(resName, "ILE")) { rtNdx = 10;nx = 2; }
    else if (!strcmp(resName, "LEU")) { rtNdx = 11;nx = 2; }
    else if (!strcmp(resName, "LYS")) { rtNdx = 12;nx = 4; }
    else if (!strcmp(resName, "MET")) { rtNdx = 13;nx = 3; }
    else if (!strcmp(resName, "PHE")) { rtNdx = 14;nx = 2; }
    else if (!strcmp(resName, "PRO")) { rtNdx = 15;nx = 2; }
    else if (!strcmp(resName, "SER")) { rtNdx = 16;nx = 1; }
    else if (!strcmp(resName, "THR")) { rtNdx = 17;nx = 1; }
    else if (!strcmp(resName, "TRP")) { rtNdx = 18;nx = 2; }
    else if (!strcmp(resName, "TYR")) { rtNdx = 19;nx = 2; }
    else if (!strcmp(resName, "VAL")) { rtNdx = 20;nx = 1; }

    //printf("rotamer type %s cannot be identified, continue\n",resname);
    if (rtNdx < 0) continue;
    //the number of current rotamer type + 1
    IntArraySet(&pRotLib->pps[binNdx].rotNums, rtNdx, IntArrayGet(&pRotLib->pps[binNdx].rotNums, rtNdx) + 1);
    pRotLib->pps[binNdx].tors[rtNdx] = (DoubleArray*)realloc(pRotLib->pps[binNdx].tors[rtNdx],
      sizeof(DoubleArray) * IntArrayGet(&pRotLib->pps[binNdx].rotNums, rtNdx));
    pRotLib->pps[binNdx].devs[rtNdx] = (DoubleArray*)realloc(pRotLib->pps[binNdx].devs[rtNdx],
      sizeof(DoubleArray) * IntArrayGet(&pRotLib->pps[binNdx].rotNums, rtNdx));
    DoubleArrayCreate(&pRotLib->pps[binNdx].tors[rtNdx][IntArrayGet(&pRotLib->pps[binNdx].rotNums, rtNdx) - 1], nx);
    DoubleArrayCreate(&pRotLib->pps[binNdx].devs[rtNdx][IntArrayGet(&pRotLib->pps[binNdx].rotNums, rtNdx) - 1], nx);
    for (int j = 0;j < nx;j++)
    {
      //remember to convert degree to rad
      DoubleArraySet(&pRotLib->pps[binNdx].tors[rtNdx][IntArrayGet(&pRotLib->pps[binNdx].rotNums, rtNdx) - 1], j, DegToRad(x[j]));
      DoubleArraySet(&pRotLib->pps[binNdx].devs[rtNdx][IntArrayGet(&pRotLib->pps[binNdx].rotNums, rtNdx) - 1], j, DegToRad(s[j]));
    }
    DoubleArrayResize(&pRotLib->pps[binNdx].prob[rtNdx], IntArrayGet(&pRotLib->pps[binNdx].rotNums, rtNdx));
    DoubleArraySet(&pRotLib->pps[binNdx].prob[rtNdx], IntArrayGet(&pRotLib->pps[binNdx].rotNums, rtNdx) - 1, prob);
    //printf("%s",buffer);

  }
  fclose(pIn);

  return Success;
}


int BBdepRotamerLibDestroy(BBdepRotLib* pThis)
{
  if (pThis->pps != NULL)
  {
    for (int i = 0;i < pThis->num;i++)
    {
      RotLibPhiPsi* pRotLibPhiPsi = &pThis->pps[i];
      for (int j = 0;j < StringArrayGetNum(&pRotLibPhiPsi->rotTypes);j++)
      {
        //free the memory for torsion angles
        DoubleArray* pTorsionsArrayForTypeJ = pRotLibPhiPsi->tors[j];
        for (int k = 0;k < IntArrayGet(&pRotLibPhiPsi->rotNums, j);k++)
        {
          DoubleArrayDestroy(&pTorsionsArrayForTypeJ[k]);
        }
        free(pTorsionsArrayForTypeJ);
        pTorsionsArrayForTypeJ = NULL;
        DoubleArray* pDeviationsArrayForTypeJ = pRotLibPhiPsi->devs[j];
        for (int k = 0;k < IntArrayGet(&pRotLibPhiPsi->rotNums, j);k++)
        {
          DoubleArrayDestroy(&pDeviationsArrayForTypeJ[k]);
        }
        free(pDeviationsArrayForTypeJ);
        pDeviationsArrayForTypeJ = NULL;
        //free the memory for probabilities
        DoubleArrayDestroy(&pRotLibPhiPsi->prob[j]);
      }
      free(pRotLibPhiPsi->prob);
      pRotLibPhiPsi->prob = NULL;
      free(pRotLibPhiPsi->tors);
      pRotLibPhiPsi->tors = NULL;
      free(pRotLibPhiPsi->devs);
      pRotLibPhiPsi->devs = NULL;
      IntArrayDestroy(&pRotLibPhiPsi->rotNums);
      StringArrayDestroy(&pRotLibPhiPsi->rotTypes);
    }
    free(pThis->pps);
    pThis->pps = NULL;
    pThis->num = 0;
  }
  return Success;
}


int RotamerOfProteinGenerateByBBdepRot(Rotamer* pThis, Residue* pResi, char* rotamerType, char* patchType, DoubleArray* torsions, AtomParams* atomParams, ResTopos* resiTopos)
{
  if (strlen(rotamerType) > MAX_LEN_RES_NAME) return NameError;
  RotamerDestroy(pThis);
  RotamerCreate(pThis);
  strcpy(pThis->type, rotamerType);

  //First step, initialize atoms and bonds
  int result = RotamerOfProteinInitAtomsAndBonds_Charmm19(pThis, pResi, atomParams, resiTopos);
  if (FAILED(result))
  {
    return result;
  }

  //Second step, patch the rotamer if necessary
  if (patchType != NULL && strcmp(patchType, "") != 0)
  {
    result = RotamerOfProteinPatch(pThis, patchType, atomParams, resiTopos);
    if (FAILED(result))
    {
      return result;
    }
  }

  //Third step, calculate coordinates of all atoms
  result = RotamerOfProteinCalcXYZ(pThis, pResi, patchType, torsions, resiTopos);
  if (FAILED(result))
  {
    return result;
  }

  //Fourth step, set the chain name and posInChain to be consistent with the residue
  RotamerSetChainName(pThis, ResidueGetChainName(pResi));
  RotamerSetPosInChain(pThis, ResidueGetPosInChain(pResi));
  return Success;
}


int RotamerSetOfProteinGenerateByBBdepRotLib(Rotamers* pThis, Residue* pRes, StringArray* designTypes, StringArray* patchTypes, BBdepRotLib* bbrotlib, AtomParams* atomParams, ResTopos* resiTopo)
{
  int binNdx = -1;
  if (pRes->phipsi[0] >= -180 && pRes->phipsi[0] <= 180 && pRes->phipsi[1] >= -180 && pRes->phipsi[1] <= 180)
  {
    int phiindex = (int)(pRes->phipsi[0] + 180) / 10;
    int psiindex = (int)(pRes->phipsi[1] + 180) / 10;
    phiindex = phiindex < 36 ? phiindex : 35;
    phiindex = phiindex >= 0 ? phiindex : 0;
    psiindex = psiindex < 36 ? psiindex : 35;
    psiindex = psiindex >= 0 ? psiindex : 0;
    binNdx = (phiindex * 36 + psiindex);
  }

  RotLibPhiPsi* pRotLibPhiPsi = &bbrotlib->pps[binNdx];
  int typeNum = StringArrayGetNum(designTypes);
  for (int typeIndex = 0; typeIndex < typeNum; typeIndex++) {
    char* typeName = StringArrayGet(designTypes, typeIndex);
    char* patchName = StringArrayGet(patchTypes, typeIndex);
    int rotNum = RotLibPhiPsiGetNum(pRotLibPhiPsi, typeName);

    // new code, faster than below method
    Rotamer newRot;
    double probability;
    DoubleArray torsions;
    RotamerCreate(&newRot);
    DoubleArrayCreate(&torsions, 0);
    // for each rotamer type, calculate the first rotamer coordinates
    RotLibPhiPsiGet(pRotLibPhiPsi, typeName, 0, &torsions, &probability);
    int result = RotamerOfProteinGenerateByBBdepRot(&newRot, pRes, typeName, patchName, &torsions, atomParams, resiTopo);
    if (FAILED(result)) {
      char errMsg[MAX_LEN_ERR_MSG+1];
      sprintf(errMsg, "in file %s line %d, failed to create rotamer %s on residue %s", __FILE__, __LINE__, typeName, ResidueGetName(pRes));
      TraceError(errMsg, result);
      return result;
    }
    DoubleArrayCopy(&newRot.Xs, &torsions);
    RotamerCalcDunbrackEnergy(&newRot, probability);
    RotamersAdd(pThis, &newRot);

    // for the other rots, just calculate the coordinates, don't have to deal with atoms and bonds again
    for (int rotNdx = 1; rotNdx < rotNum; ++rotNdx) {
      RotLibPhiPsiGet(pRotLibPhiPsi, typeName, rotNdx, &torsions, &probability);
      if (probability < ROTLIB_LOW_PROB_CUTOFF) break;
      // set the coordinates of side-chain atoms to be false
      for (int i = 0; i < RotamerGetAtomNum(&newRot); i++) {
        Atom* pAtom = RotamerGetAtom(&newRot, i);
        if (pAtom->isBBAtom == FALSE && strcmp(pAtom->name, "CB") != 0) pAtom->isXyzValid = FALSE;
      }
      RotamerOfProteinCalcXYZ(&newRot, pRes, patchName, &torsions, resiTopo);
      DoubleArrayCopy(&newRot.Xs, &torsions);
      RotamerCalcDunbrackEnergy(&newRot, probability);
      RotamersAdd(pThis, &newRot);
    }
    RotamerDestroy(&newRot);
    DoubleArrayDestroy(&torsions);
  }

  return Success;
}


int RotamerCalcDunbrackEnergy(Rotamer* pThis, double probability)
{
  if (strcmp(RotamerGetType(pThis), "ALA") == 0 || strcmp(RotamerGetType(pThis), "GLY") == 0) pThis->dunbrack = 0;
  else
  {
    double delta = 1e-7;
    pThis->dunbrack = -1.0 * log(probability + delta);
  }
  return Success;
}


BOOL RotamerAndRotamerInSameType(Rotamer* pThis, Rotamer* pOther)
{
  if (strcmp(RotamerGetType(pThis), RotamerGetType(pOther)) == 0 ||
    (strcmp(RotamerGetType(pThis), "HSD") == 0 && strcmp(RotamerGetType(pOther), "HSE") == 0) ||
    (strcmp(RotamerGetType(pThis), "HSE") == 0 && strcmp(RotamerGetType(pOther), "HSD") == 0))
  {
    return TRUE;
  }
  return FALSE;
}


BOOL RotamerAndResidueInSameType(Rotamer* pThis, Residue* pOther)
{
  if (strcmp(RotamerGetType(pThis), ResidueGetName(pOther)) == 0 ||
    (strcmp(RotamerGetType(pThis), "HSD") == 0 && strcmp(ResidueGetName(pOther), "HSE") == 0) ||
    (strcmp(RotamerGetType(pThis), "HSE") == 0 && strcmp(ResidueGetName(pOther), "HSD") == 0))
  {
    return TRUE;
  }
  return FALSE;
}


BOOL RotamerAndResidueWithSimilarTorsions(Rotamer* pThis, Residue* pOther, double cutoff)
{
  if (RotamerAndResidueInSameType(pThis, pOther))
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
      if ((strcmp(RotamerGetType(pThis), "PHE") == 0 && j == 1) ||
        (strcmp(RotamerGetType(pThis), "TYR") == 0 && j == 1) ||
        (strcmp(RotamerGetType(pThis), "ASP") == 0 && j == 1) ||
        strcmp(RotamerGetType(pThis), "GLU") == 0 && j == 2)
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


int ResidueCalcSidechainTorsion(Residue* pThis, ResTopos* pResiTopos)
{
  char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
  int torNum = 0;
  char resName[MAX_LEN_RES_NAME+1];
  strcpy(resName, ResidueGetName(pThis));
  if (!strcmp(resName, "ALA")) { torNum = 0; }
  else if (!strcmp(resName, "ARG")) { torNum = 4; }
  else if (!strcmp(resName, "ASN")) { torNum = 2; }
  else if (!strcmp(resName, "ASP")) { torNum = 2; }
  else if (!strcmp(resName, "CYS")) { torNum = 1; }
  else if (!strcmp(resName, "GLN")) { torNum = 3; }
  else if (!strcmp(resName, "GLU")) { torNum = 3; }
  else if (!strcmp(resName, "GLY")) { torNum = 0; }
  else if (!strcmp(resName, "HSD")) { torNum = 2; }
  else if (!strcmp(resName, "HSE")) { torNum = 2; }
  else if (!strcmp(resName, "ILE")) { torNum = 2; }
  else if (!strcmp(resName, "LEU")) { torNum = 2; }
  else if (!strcmp(resName, "LYS")) { torNum = 4; }
  else if (!strcmp(resName, "MET")) { torNum = 3; }
  else if (!strcmp(resName, "PHE")) { torNum = 2; }
  else if (!strcmp(resName, "PRO")) { torNum = 2; }
  else if (!strcmp(resName, "SER")) { torNum = 1; }
  else if (!strcmp(resName, "THR")) { torNum = 1; }
  else if (!strcmp(resName, "TRP")) { torNum = 2; }
  else if (!strcmp(resName, "TYR")) { torNum = 2; }
  else if (!strcmp(resName, "VAL")) { torNum = 1; }

  ResTopology resiTop;
  ResTopologyCreate(&resiTop);
  ResToposGet(pResiTopos, resName, &resiTop);
  for (int torsionIndex = 0;torsionIndex < torNum;torsionIndex++)
  {
    Type_ProteinAtomOrder desiredAtomBOrder = Type_ProteinAtomOrder_FromInt(torsionIndex);
    Type_ProteinAtomOrder desiredAtomCOrder = Type_ProteinAtomOrder_FromInt(torsionIndex + 1);
    CharmmIC icOfCurrentTorsion;
    CharmmICCreate(&icOfCurrentTorsion);
    BOOL icFound = FALSE;
    for (int icIndex = 0; icIndex < ResTopologyGetCharmmICNum(&resiTop); icIndex++)
    {
      Type_ProteinAtomOrder atomBOrder;
      Type_ProteinAtomOrder atomCOrder;
      ResTopologyGetCharmmIC(&resiTop, icIndex, &icOfCurrentTorsion);
      atomBOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomB(&icOfCurrentTorsion));
      atomCOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomC(&icOfCurrentTorsion));
      if (desiredAtomBOrder == atomBOrder && desiredAtomCOrder == atomCOrder)
      {
        icFound = TRUE;
        break;
      }
    }

    if (!icFound)
    {
      sprintf(errMsg, "in file %s line %d, cannot find IC of the %d-th torsion for residue %s", __FILE__, __LINE__, torsionIndex + 1, resName);
      TraceError(errMsg, DataNotExistError);
      CharmmICDestroy(&icOfCurrentTorsion);
      return DataNotExistError;
    }

    Atom* pAtomA = ResidueGetAtomByName(pThis, CharmmICGetAtomA(&icOfCurrentTorsion));
    Atom* pAtomB = ResidueGetAtomByName(pThis, CharmmICGetAtomB(&icOfCurrentTorsion));
    Atom* pAtomC = ResidueGetAtomByName(pThis, CharmmICGetAtomC(&icOfCurrentTorsion));
    Atom* pAtomD = ResidueGetAtomByName(pThis, CharmmICGetAtomD(&icOfCurrentTorsion));
    double torsion = GetTorsionAngle(&pAtomA->xyz, &pAtomB->xyz, &pAtomC->xyz, &pAtomD->xyz);
    DoubleArrayAppend(&pThis->Xs, torsion);

    CharmmICDestroy(&icOfCurrentTorsion);
  }
  ResTopologyDestroy(&resiTop);
  return Success;
}


Rotamer* RotamerSetFindFirstGivenTypeRotamer(Rotamers* pThis, char* type)
{
  for (int i = 0;i < RotamersGetNum(pThis);i++)
  {
    Rotamer* pRotamer = RotamersGet(pThis, i);
    if (strcmp(RotamerGetType(pRotamer), type) == 0)
    {
      return pRotamer;
    }
  }
  return NULL;
}


BOOL RotamerIsSymmetricalCheck(Rotamer* pRotamer)
{
  if (strcmp(RotamerGetType(pRotamer), "PHE") == 0 ||
    strcmp(RotamerGetType(pRotamer), "TYR") == 0 ||
    strcmp(RotamerGetType(pRotamer), "ASP") == 0 ||
    strcmp(RotamerGetType(pRotamer), "GLU") == 0 ||
    strcmp(RotamerGetType(pRotamer), "ARG") == 0)
  {
    return TRUE;
  }
  return FALSE;
}


int SymmetricalRotamerGenerate(Rotamer* pSymmetrical, Rotamer* pRotamer)
{
  RotamerCopy(pSymmetrical, pRotamer);
  if (strcmp(RotamerGetType(pRotamer), "PHE") == 0 || strcmp(RotamerGetType(pRotamer), "TYR") == 0)
  {
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "CD1"), "TMP");
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "CD2"), "CD1");
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "TMP"), "CD2");

    AtomSetName(RotamerGetAtomByName(pSymmetrical, "CE1"), "TMP");
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "CE2"), "CE1");
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "TMP"), "CE2");
  }
  else if (strcmp(RotamerGetType(pRotamer), "ASP") == 0)
  {
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "OD1"), "TMP");
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "OD2"), "OD1");
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "TMP"), "OD2");
  }
  else if (strcmp(RotamerGetType(pRotamer), "GLU") == 0)
  {
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "OE1"), "TMP");
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "OE2"), "OE1");
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "TMP"), "OE2");
  }
  else if (strcmp(RotamerGetType(pRotamer), "ARG") == 0)
  {
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "NH1"), "TMP");
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "NH2"), "NH1");
    AtomSetName(RotamerGetAtomByName(pSymmetrical, "TMP"), "NH2");
  }

  return Success;
}


int BBdepRotamerLibFromText2Binary(char* textlibfile, char* binlibfile)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  FILE* FOUT = fopen(binlibfile, "wb");
  FILE* FIN = fopen(textlibfile, "r");
  if (FIN == NULL)
  {
    sprintf(errMsg, "in file %s line %d, cannot read file %s", __FILE__, __LINE__, textlibfile);
    TraceError(errMsg, IOError);
    exit(IOError);
  }
  char line[MAX_LEN_ONE_LINE_CONTENT+1];
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, FIN))
  {
    if (line[0] == '#') continue;
    char aaname[4];
    int phi, psi;
    int a, b, c, d, e;
    float prob, x1, x2, x3, x4, v1, v2, v3, v4;
    sscanf(line, "%s %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f\n",
      aaname, &phi, &psi, &a, &b, &c, &d, &e, &prob, &x1, &x2, &x3, &x4, &v1, &v2, &v3, &v4);
    if (phi == 180 || psi == 180) continue;
    fwrite((char*)&prob, sizeof(char), 4, FOUT);
    /*short s1=(short)(10*x1);
    short s2=(short)(10*x2);
    short s3=(short)(10*x3);
    short s4=(short)(10*x4);
    short s5=(short)(10*v1);
    short s6=(short)(10*v2);
    short s7=(short)(10*v3);
    short s8=(short)(10*v4);
    fwrite((char*)&s1,sizeof(char),2,FOUT);
    fwrite((char*)&s2,sizeof(char),2,FOUT);
    fwrite((char*)&s3,sizeof(char),2,FOUT);
    fwrite((char*)&s4,sizeof(char),2,FOUT);
    fwrite((char*)&s5,sizeof(char),2,FOUT);
    fwrite((char*)&s6,sizeof(char),2,FOUT);
    fwrite((char*)&s7,sizeof(char),2,FOUT);
    fwrite((char*)&s8,sizeof(char),2,FOUT);*/
    fwrite((char*)&x1, sizeof(char), 4, FOUT);
    fwrite((char*)&x2, sizeof(char), 4, FOUT);
    fwrite((char*)&x3, sizeof(char), 4, FOUT);
    fwrite((char*)&x4, sizeof(char), 4, FOUT);
    fwrite((char*)&v1, sizeof(char), 4, FOUT);
    fwrite((char*)&v2, sizeof(char), 4, FOUT);
    fwrite((char*)&v3, sizeof(char), 4, FOUT);
    fwrite((char*)&v4, sizeof(char), 4, FOUT);
  }
  fclose(FIN);
  fclose(FOUT);
  return Success;
}


int BBdepRotamerLibFromBinary2Text(char* binlibfile, char* textlibfile)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  FILE* FOUT = fopen(textlibfile, "w");
  FILE* FIN = fopen(binlibfile, "rb");
  if (FIN == NULL)
  {
    sprintf(errMsg, "in file %s line %d, cannot read binary file %s", __FILE__, __LINE__, binlibfile);
    TraceError(errMsg, IOError);
    exit(IOError);
  }
  while (!feof(FIN))
  {
    float p;
    if (fread((char*)&p, sizeof(char), 4, FIN) < 4)
    {
      break;
    }
    fprintf(FOUT, "%8.6f", p);
    //short val;
    float val;
    for (int i = 0;i < 8;++i)
    {
      //fread((char*)&val,sizeof(char),2,FIN);
      //fprintf(FOUT,"%7.1f",(float)val/10.0);
      fread((char*)&val, sizeof(char), 4, FIN);
      fprintf(FOUT, "%7.1f", val);
    }
    fprintf(FOUT, "\n");
  }
  fclose(FIN);
  fclose(FOUT);
  return Success;
}


int BBdepRotamerLibCreate2(BBdepRotLib* pRotLib, char* binlibfile)
{
  //ACDEFGHIKLMNPQRSTVWY, only for regular amino acid
  int nx[20] = { 0,1,2,3,2,0,2,2,4,2,3,2,2,3,4,1,1,1,2,2 };
  int nrot[20] = { 0,3,18,54,18,0,36,9,73,9,27,36,2,108,75,3,3,3,36,18 };
  int lrot[20] = { 0,129,111,240,448,0,294,330,348,339,421,75,466,132,0,468,471,528,474,510 };
  // set all the possible rotamer types
  pRotLib->num = 1296;
  pRotLib->pps = (RotLibPhiPsi*)malloc(sizeof(RotLibPhiPsi) * pRotLib->num);
  for (int i = 0;i < pRotLib->num;i++)
  {
    RotLibPhiPsiCreate(&pRotLib->pps[i]);
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "ALA");//index: 0
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "ARG");//index: 1
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "ASN");//index: 2
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "ASP");//index: 3
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "CYS");//index: 4
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "GLN");//index: 5
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "GLU");//index: 6
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "GLY");//index: 7
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "HSD");//index: 8
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "HSE");//index: 9
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "ILE");//index:10
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "LEU");//index:11
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "LYS");//index:12
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "MET");//index:13
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "PHE");//index:14
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "PRO");//index:15
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "SER");//index:16
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "THR");//index:17
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "TRP");//index:18
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "TYR");//index:19
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "VAL");//index:20
    StringArrayAppend(&pRotLib->pps[i].rotTypes, "HSP");//index:21
    IntArrayResize(&pRotLib->pps[i].rotNums, StringArrayGetNum(&pRotLib->pps[i].rotTypes));
    pRotLib->pps[i].tors = (DoubleArray**)malloc(sizeof(DoubleArray*) * StringArrayGetNum(&pRotLib->pps[i].rotTypes));
    pRotLib->pps[i].devs = (DoubleArray**)malloc(sizeof(DoubleArray*) * StringArrayGetNum(&pRotLib->pps[i].rotTypes));
    pRotLib->pps[i].prob = (DoubleArray*)malloc(sizeof(DoubleArray) * StringArrayGetNum(&pRotLib->pps[i].rotTypes));
    for (int rotTypeIndex = 0;rotTypeIndex < StringArrayGetNum(&pRotLib->pps[i].rotTypes);rotTypeIndex++)
    {
      int aaIndex = AA3GetIndex(StringArrayGet(&pRotLib->pps[i].rotTypes, rotTypeIndex));
      IntArraySet(&pRotLib->pps[i].rotNums, rotTypeIndex, nrot[aaIndex]);
      pRotLib->pps[i].tors[rotTypeIndex] = (DoubleArray*)malloc(sizeof(DoubleArray) * nrot[aaIndex]);
      pRotLib->pps[i].devs[rotTypeIndex] = (DoubleArray*)malloc(sizeof(DoubleArray) * nrot[aaIndex]);
      for (int rotIdx = 0;rotIdx < nrot[aaIndex];rotIdx++)
      {
        DoubleArrayCreate(&pRotLib->pps[i].tors[rotTypeIndex][rotIdx], nx[aaIndex]);
        DoubleArrayCreate(&pRotLib->pps[i].devs[rotTypeIndex][rotIdx], nx[aaIndex]);
      }
      DoubleArrayCreate(&pRotLib->pps[i].prob[rotTypeIndex], nrot[aaIndex]);
    }
  }

  // parse the binary rotamer library file
  FILE* FIN = fopen(binlibfile, "rb");
  for (int x = 0;x < 36;x++)
  {
    for (int y = 0;y < 36;y++)
    {
      int binIdx = 36 * x + y;
      for (int rotTypeIdx = 0;rotTypeIdx < StringArrayGetNum(&pRotLib->pps[binIdx].rotTypes);rotTypeIdx++)
      {
        int aaIndex = AA3GetIndex(StringArrayGet(&pRotLib->pps[binIdx].rotTypes, rotTypeIdx));
        if (aaIndex == 0 || aaIndex == 5) continue;
        int nchi = nx[aaIndex];
        //the last '36' means 36 bytes
        fseek(FIN, (1296 * lrot[aaIndex] + binIdx * nrot[aaIndex]) * 36, SEEK_SET);
        for (int rotIdx = 0;rotIdx < nrot[aaIndex];rotIdx++)
        {
          float p = 0.;
          fread((char*)&p, sizeof(char), 4, FIN);
          DoubleArraySet(&pRotLib->pps[binIdx].prob[rotTypeIdx], rotIdx, p);
          float val;
          for (int k = 0;k < nchi;k++)
          {
            fread((char*)&val, sizeof(char), 4, FIN);
            DoubleArraySet(&pRotLib->pps[binIdx].tors[rotTypeIdx][rotIdx], k, DegToRad((double)val));
          }
          fseek(FIN, (4 - nchi) * 4, SEEK_CUR);//skip the other chi angles
          for (int k = 0;k < nchi;k++)
          {
            fread((char*)&val, sizeof(char), 4, FIN);
            DoubleArraySet(&pRotLib->pps[binIdx].devs[rotTypeIdx][rotIdx], k, DegToRad((double)val));
          }
          fseek(FIN, (4 - nchi) * 4, SEEK_CUR);//skip the other chi angles
        }
      }
    }
  }
  fclose(FIN);

  return Success;
}


