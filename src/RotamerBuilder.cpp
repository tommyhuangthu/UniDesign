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
#pragma warning(disable:6308)
#pragma warning(disable:28182)
#pragma warning(disable:26812)
#pragma warning(disable:6387)

#include "RotamerBuilder.h"
#include "EnergyFunction.h"
#include <string.h>
#include <ctype.h>

extern double CUT_PPI_DIST_SHELL1;
extern double CUT_PPI_DIST_SHELL2;
extern double CUT_PLI_DIST_SHELL1;
extern double CUT_PLI_DIST_SHELL2;
extern BOOL FLAG_DESIGN_MON;
extern BOOL FLAG_DESIGN_PPI;
extern BOOL FLAG_DESIGN_PLI;
extern BOOL FLAG_DESIGN_ENZ;
extern BOOL FLAG_USE_INPUT_SC;
extern BOOL FLAG_ROTATE_HYDROXYL;
extern BOOL FLAG_WILDTYPE_ONLY;
extern char DES_CHAINS[MAX_LEN_ONE_LINE_CONTENT+1];
extern BOOL FLAG_EXCL_CYS_ROTS;

int StructureGetTruncatedBackbone(Structure* pStruct, Residue* pSmallMol, double range, BOOL withH, AtomArray* pBackbone)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain* pChain = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChain); j++)
    {
      Residue* pRes = ChainGetResidue(pChain, j);
      if (strcmp(ResidueGetChainName(pRes), ResidueGetChainName(pSmallMol)) == 0 && ResidueGetPosInChain(pRes) == ResidueGetPosInChain(pSmallMol))
      {
        continue;
      }
      if (AtomArrayCalcMinDistance(ResidueGetAllAtoms(pRes), ResidueGetAllAtoms(pSmallMol)) > range)
      {
        continue;
      }
      for (int k = 0; k < ResidueGetAtomNum(pRes); k++)
      {
        Atom* pAtom = ResidueGetAtom(pRes, k);
        if (AtomIsHydrogen(pAtom) == TRUE && withH == FALSE)
        {
          continue;
        }
        if (pAtom->isBBAtom == FALSE && strcmp(AtomGetName(pAtom), "CB") != 0)
        {
          continue;
        }
        AtomArrayAppend(pBackbone, pAtom);
      }
    }
  }

  return Success;
}


int StructureGetTruncatedBackbone2(Structure* pStruct, Residue* pSmallMol, BOOL withH, AtomArray* pBackbone)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain* pChain = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChain); j++)
    {
      Residue* pRes = ChainGetResidue(pChain, j);
      if (strcmp(ResidueGetChainName(pRes), ResidueGetChainName(pSmallMol)) == 0 && ResidueGetPosInChain(pRes) == ResidueGetPosInChain(pSmallMol))
      {
        continue;
      }
      for (int k = 0; k < ResidueGetAtomNum(pRes); k++)
      {
        Atom* pAtom = ResidueGetAtom(pRes, k);
        if (AtomIsHydrogen(pAtom) == TRUE && withH == FALSE)
        {
          continue;
        }
        if (pAtom->isBBAtom == FALSE && strcmp(AtomGetName(pAtom), "CB") != 0)
        {
          continue;
        }
        AtomArrayAppend(pBackbone, pAtom);
      }
    }
  }

  return Success;
}


int StructureDeployCataConsSitePair(Structure* pStruct, CataConsSitePair* pCataConsSitePair)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  Residue* pSmallMol = NULL;
  Rotamers* pRotSetOfFirstAndSecondSite[2] = {NULL, NULL};

  int result = StructureFindSmallMol(pStruct, &pSmallMol);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, cannot find small molecule", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }

  Rotamers tempSmallRotSet;
  RotamersCreate(&tempSmallRotSet);
  Rotamer tempSmallRot;
  RotamerCreate(&tempSmallRot);
  RotamerSetType(&tempSmallRot, ResidueGetName(pSmallMol));
  RotamerAddAtoms(&tempSmallRot, ResidueGetAllAtoms(pSmallMol));
  RotamersAdd(&tempSmallRotSet, &tempSmallRot);
  RotamerDestroy(&tempSmallRot);
  for (int i = 0; i < 2; i++)
  {
    int posInChain;
    char* chainName;
    char* residueName;
    if (i == 0)
    {
      posInChain = pCataConsSitePair->pos1;
      chainName = pCataConsSitePair->chnName1;
      residueName = pCataConsSitePair->resName1;
    }
    else
    {
      posInChain = pCataConsSitePair->pos2;
      chainName = pCataConsSitePair->chnName2;
      residueName = pCataConsSitePair->resName2;
    }

    pRotSetOfFirstAndSecondSite[i] = NULL;
    if (posInChain == ResidueGetPosInChain(pSmallMol) && strcmp(chainName, ResidueGetChainName(pSmallMol)) == 0)
    {
      pRotSetOfFirstAndSecondSite[i] = &tempSmallRotSet;
    }
    else
    {
      int chainIndex = -1, resiIndex = -1;
      StructureFindChainIndex(pStruct, chainName, &chainIndex);
      ChainFindResidueIndexByPosInChain(StructureGetChain(pStruct, chainIndex), posInChain, &resiIndex);
      DesignSite* pDesignSite = StructureFindDesignSite(pStruct, chainIndex, resiIndex);
      pRotSetOfFirstAndSecondSite[i] = DesignSiteGetRotamers(pDesignSite);
    }
    if (pRotSetOfFirstAndSecondSite[i] == NULL)
    {
      result = DataNotExistError;
      sprintf(errMsg, "in file %s line %d, cannot find site %s %s %d", __FILE__, __LINE__, chainName, residueName, posInChain);
      TraceError(errMsg, result);
      return result;
    }
  }
  result = CataConsSitePairDeploy(pCataConsSitePair, pRotSetOfFirstAndSecondSite[0], pRotSetOfFirstAndSecondSite[1]);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, cannot deploy catalytic constraints between %s%d%s and %s%d%s", __FILE__, __LINE__,
      pCataConsSitePair->chnName1, pCataConsSitePair->pos1, pCataConsSitePair->resName1,
      pCataConsSitePair->chnName2, pCataConsSitePair->pos2, pCataConsSitePair->resName2);
    TraceError(errMsg, result);
    return result;
  }

  RotamersDestroy(&tempSmallRotSet);
  return result;
}


int StructurePlaceSmallMol(Structure* pStruct, PlacingRule* pRule, CataConsSitePairArray* pConsArray, int siteNum, DesignSite** ppSites, Rotamers* pSmallSet)
{
  int result;
  char errMsg[MAX_LEN_ERR_MSG+1];
  int chnNdx = -1, resNdx = -1;
  StructureFindChainIndex(pStruct, PlacingRuleGetChainName(pRule), &chnNdx);
  ChainFindResidueIndexByPosInChain(StructureGetChain(pStruct, chnNdx), PlacingRuleGetPosInChain(pRule), &resNdx);
  DesignSite* pStartSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pStartSite == NULL)
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d, starting site %s%d%s is required for ligand placement", __FILE__, __LINE__,
      PlacingRuleGetChainName(pRule), PlacingRuleGetPosInChain(pRule), PlacingRuleGetResiName(pRule));
    TraceError(errMsg, result);
    return result;
  }

  result = PlacingRulePlaceSmallMol(pRule, pConsArray, pStartSite, siteNum, ppSites, pSmallSet);
  if (FAILED(result))
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d, failed to place ligand poses", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }

  return Success;
}


int StructureGenerateSmallMolRotamers(Structure* pStruct, char* cataconsFile, char* ligplaceFile)
{
  int result = Success;
  char errMsg[MAX_LEN_ERR_MSG+1];

  printf("make ligand poses for ligand pose sampling in enzyme design\n");

  printf("1. read and deploy catalytic constraints ...\n");
  CataConsSitePairArray cataCons;
  result = CataConsSitePairArrayCreate(&cataCons, cataconsFile);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, failed to read catalytic-constraint file %s", __FILE__, __LINE__, cataconsFile);
    result = ValueError;
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0;i < CataConsSitePairArrayGetNum(&cataCons);i++)
  {
    result = StructureDeployCataConsSitePair(pStruct, CataConsSitePairArrayGet(&cataCons, i));
    if (FAILED(result))
    {
      sprintf(errMsg, "in file %s line %d, failed to deploy catalytic constraints", __FILE__, __LINE__);
      result = ValueError;
      TraceError(errMsg, result);
      return result;
    }
  }

  printf("2. read and deploy ligand placing rules ... \n");
  PlacingRule placingRule;
  result = PlacingRuleCreate(&placingRule, ligplaceFile);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, failed to read ligand-placing file %s", __FILE__, __LINE__, ligplaceFile);
    result = ValueError;
    return result;
  }

  Residue* pSmallMol = NULL;
  result = StructureFindSmallMol(pStruct, &pSmallMol);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, cannot find small molecule", __FILE__, __LINE__);
    result = DataNotExistError;
    return result;
  }
  AtomArray truncBackbone;
  AtomArrayCreate(&truncBackbone);
  //StructureGetTruncatedBackbone(pStruct, pSmallMol, PlacingRuleGetTruncatedBackboneRange(&placingRule), TRUE, &truncBackbone);
  StructureGetTruncatedBackbone2(pStruct, pSmallMol, TRUE, &truncBackbone);

  int nRelatedSites = 0;
  DesignSite** ppRelatedSites = NULL;
  for (int i = 0; i < CataConsSitePairArrayGetNum(&cataCons); i++)
  {
    CataConsSitePair* pSitePair = CataConsSitePairArrayGet(&cataCons, i);
    for (int k = 0; k < 2; k++)
    {
      DesignSite* pSite = NULL;
      if (k == 0)
      {
        if (ResidueGetPosInChain(pSmallMol) == pSitePair->pos1 && strcmp(ResidueGetChainName(pSmallMol), pSitePair->chnName1) == 0)
        {
          continue;
        }
        int chainIndex = -1, resiIndex = -1;
        StructureFindChainIndex(pStruct, pSitePair->chnName1, &chainIndex);
        ChainFindResidueIndexByPosInChain(StructureGetChain(pStruct, chainIndex), pSitePair->pos1, &resiIndex);
        pSite = StructureFindDesignSite(pStruct, chainIndex, resiIndex);
      }
      else
      {
        if (ResidueGetPosInChain(pSmallMol) == pSitePair->pos2 && strcmp(ResidueGetChainName(pSmallMol), pSitePair->chnName2) == 0)
        {
          continue;
        }
        int chainIndex = -1, resiIndex = -1;
        StructureFindChainIndex(pStruct, pSitePair->chnName2, &chainIndex);
        ChainFindResidueIndexByPosInChain(StructureGetChain(pStruct, chainIndex), pSitePair->pos2, &resiIndex);
        pSite = StructureFindDesignSite(pStruct, chainIndex, resiIndex);
      }
      // small molecule site is excluded
      BOOL designsiteExist = FALSE;
      for (int j = 0; j < nRelatedSites; j++)
      {
        DesignSite* pSite1 = ppRelatedSites[j];
        // skip identical sites
        if (pSite1 == pSite)
        {
          designsiteExist = TRUE;
          break;
        }
      }
      if (designsiteExist == FALSE)
      {
        nRelatedSites++;
        ppRelatedSites = (DesignSite**)realloc(ppRelatedSites, sizeof(DesignSite*) * nRelatedSites);
        ppRelatedSites[nRelatedSites - 1] = pSite;
      }
    }
  }
  result = PlacingRuleDeploy(&placingRule, pSmallMol, &cataCons, nRelatedSites, ppRelatedSites, &truncBackbone);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, failed to deploy placing rule", __FILE__, __LINE__);
    result = ValueError;
    TraceError(errMsg, result);
    return result;
  }

  printf("3. start placing ligand poses ...\n");
  Rotamers ligRots;
  RotamersCreate(&ligRots);
  result = StructurePlaceSmallMol(pStruct, &placingRule, &cataCons, nRelatedSites, ppRelatedSites, &ligRots);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, failed to place ligand poses", __FILE__, __LINE__);
    result = ValueError;
    return result;
  }
  if (RotamersGetNum(&ligRots) > 0)
  {
    BOOL smallmolRotSetGenerated = FALSE;
    for (int i = 0;i < nRelatedSites;i++)
    {
      if (ppRelatedSites[i]->pRes == pSmallMol)
      {
        for (int j = 0;j < RotamersGetNum(&ligRots);j++)
        {
          RotamersAdd(DesignSiteGetRotamers(ppRelatedSites[i]), RotamersGet(&ligRots, j));
        }
        smallmolRotSetGenerated = TRUE;
        break;
      }
    }
    if (smallmolRotSetGenerated == FALSE)
    {
      int chnNdx = -1, resNdx = -1;
      StructureFindChainIndex(pStruct, ResidueGetChainName(pSmallMol), &chnNdx);
      ChainFindResidueIndexByPosInChain(StructureGetChain(pStruct, chnNdx), ResidueGetPosInChain(pSmallMol), &resNdx);
      ProteinSiteAddDesignSite(pStruct, chnNdx, resNdx);
      RotamersCopy(DesignSiteGetRotamers(StructureGetDesignSite(pStruct, pStruct->desSiteNum - 1)), &ligRots);
    }
  }

  RotamersDestroy(&ligRots);
  PlacingRuleDestroy(&placingRule);
  CataConsSitePairArrayDestroy(&cataCons);
  AtomArrayDestroy(&truncBackbone);
  free(ppRelatedSites);
  ppRelatedSites = NULL;

  return Success;
}


int StructureWriteSmallMolRotamers(Structure* pStruct, char* smallmolFile)
{
  int result;
  char errMsg[MAX_LEN_ERR_MSG+1];
  Residue* pSmallMol = NULL;
  StructureFindSmallMol(pStruct, &pSmallMol);
  if (pSmallMol == NULL)
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d, cannot find the small molecule", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }
  int chainIndex = -1, resiIndex = -1;
  StructureFindChainIndex(pStruct, ResidueGetChainName(pSmallMol), &chainIndex);
  ChainFindResidueIndexByPosInChain(StructureGetChain(pStruct, chainIndex), ResidueGetPosInChain(pSmallMol), &resiIndex);
  DesignSite* pSmallSite = StructureFindDesignSite(pStruct, chainIndex, resiIndex);
  Rotamers* pSmallSet = DesignSiteGetRotamers(pSmallSite);
  if (pSmallSet == NULL)
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d, cannot find the small-molecule design site", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }

  FILE* outputFile = NULL;
  if (smallmolFile == NULL || strcmp(smallmolFile, "") == 0) outputFile = NULL;
  else
  {
    outputFile = fopen(smallmolFile, "w");
    if (outputFile == NULL)
    {
      sprintf(errMsg, "In file %s line %d, cannot write to file %s", __FILE__, __LINE__, smallmolFile);
      result = IOError;
      TraceError(errMsg, result);
      return result;
    }
  }

  for (int i = 0; i < RotamersGetNum(pSmallSet); i++)
  {
    Rotamer newRot;
    RotamerCreate(&newRot);
    Model(i + 1, outputFile);
    RotamerCopy(&newRot, RotamersGet(pSmallSet, i));
    RotamerRestore(&newRot, pSmallSet);
    RotamerShowInPDBFormat(&newRot, "ATOM", ResidueGetChainName(pSmallMol), 1, 1, outputFile);
    RotamerExtract(&newRot);
    // save energy to the ligand pose file
    fprintf(outputFile, "ENERGY INTERNAL: %f BACKBONE: %f\n", newRot.vdwInternal, newRot.vdwBackbone);
    EndModel(outputFile);
    RotamerDestroy(&newRot);
  }
  if (outputFile != NULL) fclose(outputFile);

  return Success;
}


int StructureReadSmallMolRotamers(Structure* pStruct, ResTopos* resTopos, char* smallmolFile)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  Residue* pSmallMol = NULL;
  int result = StructureFindSmallMol(pStruct, &pSmallMol);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, cannot find small molecule", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }

  int chnNdx = -1, resNdx = -1;
  StructureFindChainIndex(pStruct, ResidueGetChainName(pSmallMol), &chnNdx);
  ChainFindResidueIndexByPosInChain(StructureGetChain(pStruct, chnNdx), ResidueGetPosInChain(pSmallMol), &resNdx);
  DesignSite* pSmallMolSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pSmallMolSite == NULL)
  {
    ProteinSiteAddDesignSite(pStruct, chnNdx, resNdx);
    pSmallMolSite = StructureGetDesignSite(pStruct, pStruct->desSiteNum - 1);
  }
  ResidueSetDesignType(pSmallMol, Type_DesType_SmallMol);

  FileReader fr;
  result = FileReaderCreate(&fr, smallmolFile);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, cannot read file %s", __FILE__, __LINE__, smallmolFile);
    TraceError(errMsg, result);
    FileReaderDestroy(&fr);
    return result;
  }
  while (!FileReaderEndOfFile(&fr))
  {
    char line[MAX_LEN_ONE_LINE_CONTENT+1];
    char keyword[MAX_LEN_ONE_LINE_CONTENT+1];
    FileReaderGetNextLine(&fr, line);
    ExtractFirstStringFromSourceString(keyword, line);
    if (strcmp(keyword, "MODEL") != 0) continue;

    Residue tmpRes;
    ResidueCreate(&tmpRes);
    ResidueCopy(&tmpRes, pSmallMol);
    for (int i = 0;i < ResidueGetAtomNum(&tmpRes);i++)
    {
      ResidueGetAtom(&tmpRes, i)->isXyzValid = FALSE;
    }
    if (FAILED(ResidueReadXYZFromPDB(&tmpRes, &fr)))
    {
      sprintf(errMsg, "in file %s line %d, error reading the ligand pose file", __FILE__, __LINE__);
      result = FormatError;
      TraceError(errMsg, result);
      return result;
    }

    ResidueCheckAtomCoordinateValidity(&tmpRes);
    result = ResidueCalcAllAtomXYZ(&tmpRes, resTopos, NULL, NULL);
    if (FAILED(result))
    {
      sprintf(errMsg, "in file %s line %d, not all coordinates can be calculated", __FILE__, __LINE__);
      result = ValueError;
      TraceError(errMsg, result);
      return result;
    }
    else
    {
      Rotamer newRot;
      RotamerCreate(&newRot);
      strcpy(newRot.type, ResidueGetName(&tmpRes));
      RotamerAddAtoms(&newRot, ResidueGetAllAtoms(&tmpRes));
      BondsCopy(RotamerGetBonds(&newRot), ResidueGetBonds(&tmpRes));
      RotamerSetChainName(&newRot, ResidueGetChainName(&tmpRes));
      //read the energy from smallmol rotamer file
      newRot.vdwInternal = tmpRes.internalEnergy;
      newRot.vdwBackbone = tmpRes.backboneEnergy;
      //set smallmol internal energy as self energy and recalculate backbone energy later
      //newRot.selfEnergy = tmpRes.internalEnergy;
      RotamerSetPosInChain(&newRot, ResidueGetPosInChain(&tmpRes));
      RotamersAdd(DesignSiteGetRotamers(pSmallMolSite), &newRot);
      RotamerDestroy(&newRot);
    }
    ResidueDestroy(&tmpRes);
  }
  FileReaderDestroy(&fr);
  return Success;
}


////////////////////////////////////////////////////////////////////////////////////////
// This function is used to screen out ligand poses with bad (meaningless) orientations
////////////////////////////////////////////////////////////////////////////////////////
// oldConfsFile records the ligand poses to be screened
// newConfsFile records the ligand poses to be saved
// orntRuleFile records the ligand pose screening rule
////////////////////////////////////////////////////////////////////////////////////////
// the following information is needed in the screening rule file:
// two atoms on small molecule: a catalytic atom and a binding atom;
// a set of residues in protein scaffold;
// the information is organized in the following format:
// CATA_ATOM  C15
// BIND_ATOM  C11
// BIND_SITE_GROUP
// BIND_SITE  CHAINB VAL 55
// BIND_SITE_GROUP
// BIND_SITE  CHAINB PHE 56
// BIND_SITE_GROUP
// BIND_SITE  CHAINB SER 66
// BIND_SITE_GROUP
// BIND_SITE  CHAINB TRP 153
////////////////////////////////////////////////////////////////////////////////////////
// CRITERIA:
// (1) in one BIND_SITE_GROUP, all BIND_SITEs must satisfy dist(CATA_ATOM, CA_ATOM_OF_BIND_SITE) <= dist(BIND_ATOM, CA_ATOM_OF_BIND_SITE)
// (2) if multi BIND_SITE_GROUPs, it is okay to have any one satisfy (1)
//////////////////////////////////////////////////////////////////////////////////////////
int StructureSmallmolOrientationScreen(Structure* pStruct, ResTopos* resTopos, char* oldposeFile, char* newposeFile, char* orientRuleFile)
{
  typedef struct _BindSite
  {
    char chainName[MAX_LEN_CHAIN_NAME+1];
    int  posInChain;
    char resiName[MAX_LEN_RES_NAME+1];
  } BindSite;

  typedef struct _BindSiteGroup
  {
    int       siteNum;
    BindSite* pSites;
  } BindSiteGroup;

  char           cataAtomName[MAX_LEN_ATOM_NAME+1];
  char           bindAtomName[MAX_LEN_ATOM_NAME+1];
  int            groupNum;
  BindSiteGroup* pGroups = NULL;
  double         MIN_DELTA_DIST = 0;

  FileReader     file;
  FileReaderCreate(&file, orientRuleFile);
  groupNum = 0;
  while (!FileReaderEndOfFile(&file))
  {
    BOOL doneInThisGroup = FALSE;
    int siteNum;
    char           line[MAX_LEN_ONE_LINE_CONTENT+1];
    while (!FAILED(FileReaderGetNextLine(&file, line)))
    {
      StringArray wordsInLine;
      StringArrayCreate(&wordsInLine);
      StringArraySplitString(&wordsInLine, line, ' ');

      if (strcmp(StringArrayGet(&wordsInLine, 0), "CATA_ATOM") == 0)
      {
        strcpy(cataAtomName, StringArrayGet(&wordsInLine, 1));
      }
      else if (strcmp(StringArrayGet(&wordsInLine, 0), "BIND_ATOM") == 0)
      {
        strcpy(bindAtomName, StringArrayGet(&wordsInLine, 1));
      }
      else if (strcmp(StringArrayGet(&wordsInLine, 0), "BIND_SITE_GROUP") == 0 && doneInThisGroup == TRUE)
      {
        FileReaderSetCurrentPos(&file, FileReaderGetCurrentPos(&file) - 1);
        break;
      }
      else if (strcmp(StringArrayGet(&wordsInLine, 0), "BIND_SITE_GROUP") == 0 && doneInThisGroup == FALSE)
      {
        doneInThisGroup = TRUE;
        groupNum++;
        pGroups = (BindSiteGroup*)realloc(pGroups, sizeof(BindSiteGroup) * groupNum);
        siteNum = 0;
        pGroups[groupNum - 1].siteNum = 0;
        pGroups[groupNum - 1].pSites = NULL;
        continue;
      }
      else if (strcmp(StringArrayGet(&wordsInLine, 0), "BIND_SITE") == 0)
      {
        siteNum++;
        pGroups[groupNum - 1].pSites = (BindSite*)realloc(pGroups[groupNum - 1].pSites, sizeof(BindSite) * siteNum);
        strcpy(pGroups[groupNum - 1].pSites[siteNum - 1].chainName, StringArrayGet(&wordsInLine, 1));
        strcpy(pGroups[groupNum - 1].pSites[siteNum - 1].resiName, StringArrayGet(&wordsInLine, 2));
        pGroups[groupNum - 1].pSites[siteNum - 1].posInChain = atoi(StringArrayGet(&wordsInLine, 3));
        pGroups[groupNum - 1].siteNum = siteNum;
      }
      else if (strcmp(StringArrayGet(&wordsInLine, 0), "MIN_DELTA_DIST") == 0)
      {
        MIN_DELTA_DIST = atof(StringArrayGet(&wordsInLine, 1));
      }
    }
  }
  FileReaderDestroy(&file);

  char errMsg[MAX_LEN_ERR_MSG+1];
  printf("reading ligand poses from file %s\n", oldposeFile);
  FILE* pIn = fopen(oldposeFile, "r");
  if (pIn == NULL)
  {
    sprintf(errMsg, "in file %s line %d, can not open file %s for reading", __FILE__, __LINE__, oldposeFile);
    TraceError(errMsg, IOError);
    return IOError;
  }
  FILE* pOut = fopen(newposeFile, "w");
  if (pOut == NULL)
  {
    sprintf(errMsg, "in file %s line %d, can not open file %s for reading", __FILE__, __LINE__, newposeFile);
    TraceError(errMsg, IOError);
    return IOError;
  }

  char badOrntPoseFile[MAX_LEN_FILE_NAME+1];
  sprintf(badOrntPoseFile, "%s_BAD_ORNT.pdb",newposeFile);
  FILE* pOutBad = fopen(badOrntPoseFile, "w");
  if (pOutBad == NULL)
  {
    sprintf(errMsg, "in file %s line %d, can not open file %s for reading", __FILE__, __LINE__, badOrntPoseFile);
    TraceError(errMsg, IOError);
    return IOError;
  }

  int totRotNum = 0;
  int accRotNdx = 0;
  int rejectedRotamerNdx = 0;
  XYZ xyzCataAtom;
  XYZ xyzBindAtom;
  BOOL cataAtomExist = FALSE;
  BOOL bindAtomExist = FALSE;
  StringArray buffer;
  StringArrayCreate(&buffer);
  char           line[MAX_LEN_ONE_LINE_CONTENT+1];
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    char keyword[10] = "";
    ExtractTargetStringFromSourceString(keyword, line, 0, 4);
    if (strcmp(keyword, "ENDM") == 0)
    {
      if (cataAtomExist == FALSE || bindAtomExist == FALSE)
      {
        sprintf(errMsg, "in file %s line %d, CataAtom or BindAtom does not exist in file %s", __FILE__, __LINE__, oldposeFile);
        TraceError(errMsg, FormatError);
        return FormatError;
      }

      BOOL curRotamerAccepted = FALSE;
      for (int i = 0; i < groupNum; i++)
      {
        curRotamerAccepted = TRUE;
        for (int j = 0; j < pGroups[i].siteNum; j++)
        {
          //05/09/2023, fix the bug for finding the binding residues
          int resNdx = -1;
          Chain* pChain = StructureFindChainByName(pStruct, pGroups[i].pSites[j].chainName);
          ChainFindResidueIndexByPosInChain(pChain, pGroups[i].pSites[j].posInChain, &resNdx);
          if (resNdx == -1)
          {
            sprintf(errMsg, "in file %s line %d, can not find binding site pos %d on chain %s", __FILE__, __LINE__, pGroups[i].pSites[j].posInChain, pGroups[i].pSites[j].chainName);
            TraceError(errMsg, ValueError);
            exit(ValueError);
          }
          Residue* pResidue = ChainGetResidue(pChain, resNdx);
          XYZ* pXyzCA = &ResidueGetAtomByName(pResidue, "CA")->xyz;
          // check the distance criterion
          if ( XYZDist(&xyzCataAtom, pXyzCA) - XYZDist(&xyzBindAtom, pXyzCA) < MIN_DELTA_DIST)
          {
            curRotamerAccepted = FALSE;
            break;
          }
        }
        if (curRotamerAccepted == TRUE)
        {
          fprintf(pOut, "MODEL     %d\n", accRotNdx + 1);
          for (int j = 0; j < StringArrayGetNum(&buffer); j++)
          {
            fprintf(pOut, "%s", StringArrayGet(&buffer, j));
          }
          fprintf(pOut, "ENDMDL\n");
          accRotNdx++;
          break;
        }
      }

      if (curRotamerAccepted == FALSE)
      {
        fprintf(pOutBad, "MODEL     %d\n", rejectedRotamerNdx + 1);
        for (int j = 0; j < StringArrayGetNum(&buffer); j++)
        {
          fprintf(pOutBad, "%s", StringArrayGet(&buffer, j));
        }
        fprintf(pOutBad, "ENDMDL\n");
        rejectedRotamerNdx++;
      }

      printf("%d / %d rotamers have been processed      \r", accRotNdx, totRotNum);
      //fflush(stdout);
    }
    else if (strcmp(keyword, "MODE") == 0)
    {
      StringArrayDestroy(&buffer);
      StringArrayCreate(&buffer);
      totRotNum++;
    }
    else if (strcmp(keyword, "ATOM") == 0)
    {
      char atomName[MAX_LEN_ATOM_NAME+1];
      char strX[MAX_LEN_ONE_LINE_CONTENT+1];
      char strY[MAX_LEN_ONE_LINE_CONTENT+1];
      char strZ[MAX_LEN_ONE_LINE_CONTENT+1];
      ExtractTargetStringFromSourceString(atomName, line, 12, 4);
      // 05/08/2023: add this if(){} to reformat the atom name
      if (isdigit(atomName[0]) && isalpha(atomName[1]))
      {
        char tempName[MAX_LEN_ATOM_NAME+1];
        strcpy(tempName, atomName);
        tempName[strlen(atomName)-1] = atomName[0];
        tempName[strlen(atomName)] = '\0';
        strcpy(atomName, tempName);
      }
      if (strcmp(atomName, cataAtomName) == 0)
      {
        cataAtomExist = TRUE;
        ExtractTargetStringFromSourceString(strX, line, 30, 8);
        xyzCataAtom.X = atof(strX);
        ExtractTargetStringFromSourceString(strY, line, 38, 8);
        xyzCataAtom.Y = atof(strY);
        ExtractTargetStringFromSourceString(strZ, line, 46, 8);
        xyzCataAtom.Z = atof(strZ);
      }

      if (strcmp(atomName, bindAtomName) == 0)
      {
        bindAtomExist = TRUE;
        ExtractTargetStringFromSourceString(strX, line, 30, 8);
        xyzBindAtom.X = atof(strX);
        ExtractTargetStringFromSourceString(strY, line, 38, 8);
        xyzBindAtom.Y = atof(strY);
        ExtractTargetStringFromSourceString(strZ, line, 46, 8);
        xyzBindAtom.Z = atof(strZ);
      }
      StringArrayAppend(&buffer, line);
    }
    else if (strcmp(keyword, "ENER") == 0)
    {
      StringArrayAppend(&buffer, line);
    }
  }
  printf("\n");
  fclose(pIn);
  fclose(pOut);
  fclose(pOutBad);

  StringArrayDestroy(&buffer);
  for (int j = 0; j < groupNum; j++)
  {
    free(pGroups[j].pSites);
    pGroups[j].siteNum = 0;
  }
  free(pGroups);

  return Success;
}


int ProteinSiteBuildAllRotamersByBBdepRotLib(Structure* pStruct, int chnNdx, int resNdx, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos)
{
  int result = Success;
  DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pSite != NULL)
  {
    DesignSiteRemoveRotamers(pSite);
  }
  else
  {
    (pStruct->desSiteNum)++;
    pStruct->desSites = (DesignSite*)realloc(pStruct->desSites, sizeof(DesignSite) * pStruct->desSiteNum);
    DesignSiteCreate(&pStruct->desSites[pStruct->desSiteNum - 1]);
    pSite = StructureGetDesignSite(pStruct, pStruct->desSiteNum - 1);
    Chain* pDestChain = StructureGetChain(pStruct, chnNdx);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resNdx);
    pSite->pRes = pDestResidue;
    pSite->chnNdx = chnNdx;
    pSite->resNdx = resNdx;
  }

  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  StringArrayAppend(&designTypes, "ALA"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ARG"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ASN"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ASP"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "CYS"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLN"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLU"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "GLY"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "HSD"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "HSE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "ILE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "LEU"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "LYS"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "MET"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "PHE"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "PRO"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "SER"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "THR"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "TRP"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "TYR"); StringArrayAppend(&patchTypes, "");
  StringArrayAppend(&designTypes, "VAL"); StringArrayAppend(&patchTypes, "");
  result = RotamerSetOfProteinGenerateByBBdepRotLib(DesignSiteGetRotamers(pSite), pSite->pRes, &designTypes, &patchTypes, rotlib, atomParams, resiTopos);
  ResidueSetDesignType(pSite->pRes, Type_DesType_Mutable);
  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}


int ProteinSiteBuildMutatedRotamersByBBdepRotLib(Structure* pStruct, int chnNdx, int resNdx, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos, StringArray* pDesTypes, StringArray* pPatchTypes)
{
  int result = Success;
  DesignSite* pCurSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pCurSite != NULL)
  {
    DesignSiteRemoveRotamers(pCurSite);
  }
  else
  {
    (pStruct->desSiteNum)++;
    pStruct->desSites = (DesignSite*)realloc(pStruct->desSites, sizeof(DesignSite) * pStruct->desSiteNum);
    DesignSiteCreate(&pStruct->desSites[pStruct->desSiteNum - 1]);
    pCurSite = StructureGetDesignSite(pStruct, pStruct->desSiteNum - 1);
    Chain* pDestChain = StructureGetChain(pStruct, chnNdx);
    Residue* pDestRes = ChainGetResidue(pDestChain, resNdx);
    pCurSite->pRes = pDestRes;
    pCurSite->chnNdx = chnNdx;
    pCurSite->resNdx = resNdx;
  }

  result = RotamerSetOfProteinGenerateByBBdepRotLib(DesignSiteGetRotamers(pCurSite), pCurSite->pRes, pDesTypes, pPatchTypes, rotlib, atomParams, resiTopos);
  ResidueSetDesignType(pCurSite->pRes, Type_DesType_Mutable);
  return result;
}


int ProteinSiteBuildWildtypeRotamersByBBdepRotLib(Structure* pStruct, int chnNdx, int resNdx, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos)
{
  int result = Success;
  DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pSite != NULL)
  {
    DesignSiteRemoveRotamers(pSite);
  }
  else
  {
    (pStruct->desSiteNum)++;
    pStruct->desSites = (DesignSite*)realloc(pStruct->desSites, sizeof(DesignSite) * pStruct->desSiteNum);
    DesignSiteCreate(&pStruct->desSites[pStruct->desSiteNum - 1]);
    pSite = StructureGetDesignSite(pStruct, pStruct->desSiteNum - 1);
    Chain* pDestChain = StructureGetChain(pStruct, chnNdx);
    Residue* pResidue = ChainGetResidue(pDestChain, resNdx);
    pSite->pRes = pResidue;
    pSite->chnNdx = chnNdx;
    pSite->resNdx = resNdx;
  }

  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  StringArrayAppend(&designTypes, ResidueGetName(pSite->pRes));
  StringArrayAppend(&patchTypes, "");
  if (!strcmp(ResidueGetName(pSite->pRes), "HSD"))
  {
    StringArrayAppend(&designTypes, "HSE");
    StringArrayAppend(&patchTypes, "");
  }
  else if (!strcmp(ResidueGetName(pSite->pRes), "HSE"))
  {
    StringArrayAppend(&designTypes, "HSD");
    StringArrayAppend(&patchTypes, "");
  }

  result = RotamerSetOfProteinGenerateByBBdepRotLib(DesignSiteGetRotamers(pSite), pSite->pRes, &designTypes, &patchTypes, rotlib, atomParams, resiTopos);
  ResidueSetDesignType(pSite->pRes, Type_DesType_Repackable);
  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}


int ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(Structure* pStruct, int chnNdx, int resNdx, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos)
{
  DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pSite != NULL) DesignSiteRemoveRotamers(pSite);
  else {
    (pStruct->desSiteNum)++;
    pStruct->desSites = (DesignSite*)realloc(pStruct->desSites, sizeof(DesignSite) * pStruct->desSiteNum);
    DesignSiteCreate(&pStruct->desSites[pStruct->desSiteNum - 1]);
    pSite = StructureGetDesignSite(pStruct, pStruct->desSiteNum - 1);
    Chain* pDestChn = StructureGetChain(pStruct, chnNdx);
    Residue* pDestRes = ChainGetResidue(pDestChn, resNdx);
    pSite->pRes = pDestRes;
    pSite->chnNdx = chnNdx;
    pSite->resNdx = resNdx;
  }

  if (ResidueGetDesignType(pSite->pRes) == Type_DesType_Fixed || ResidueGetDesignType(pSite->pRes) == Type_DesType_NatRot){
    ;
  }
  else {
    StringArray designTypes;
    StringArray patchTypes;
    StringArrayCreate(&designTypes);
    StringArrayCreate(&patchTypes);
    if (ResidueGetDesignType(pSite->pRes) == Type_DesType_Mutable) {
      int len = strlen(pSite->pRes->desAATypes);
      char* desStr = "ACDEFGHIKLMNPQRSTVWY";
      if (len == 0) {
        for (int i = 0;i < 20;i++) {
          char aa = desStr[i];
          char rotType[MAX_LEN_RES_NAME+1];
          AA1ToAA3(aa, rotType);
          if (aa == 'C' && FLAG_EXCL_CYS_ROTS) continue;
          else if (aa == 'H') {
            StringArrayAppend(&designTypes, "HSD"); 
            StringArrayAppend(&patchTypes, "");
            StringArrayAppend(&designTypes, "HSE"); 
            StringArrayAppend(&patchTypes, "");
          }
          else if (aa != 'X') { 
            StringArrayAppend(&designTypes, rotType); 
            StringArrayAppend(&patchTypes, ""); 
          }
        }
      }
      else {
        for (int i = 0;i < len;i++) {
          char aa = pSite->pRes->desAATypes[i];
          char rotType[MAX_LEN_RES_NAME+1];
          AA1ToAA3(aa, rotType);
          if (aa == 'C' && FLAG_EXCL_CYS_ROTS) continue;
          else if (aa == 'H') {
            StringArrayAppend(&designTypes, "HSD"); 
            StringArrayAppend(&patchTypes, "");
            StringArrayAppend(&designTypes, "HSE"); 
            StringArrayAppend(&patchTypes, "");
          }
          else if (aa != 'X') { 
            StringArrayAppend(&designTypes, rotType); 
            StringArrayAppend(&patchTypes, ""); 
          }
        }
      }
    }
    else if (ResidueGetDesignType(pSite->pRes) == Type_DesType_Catalytic)
    {
      int len = strlen(pSite->pRes->desAATypes);
      if (len == 0) {
        StringArrayAppend(&designTypes, ResidueGetName(pSite->pRes)); 
        StringArrayAppend(&patchTypes, "");
      }
      else {
        for (int i = 0;i < len;i++) {
          char aa = pSite->pRes->desAATypes[i];
          char rotType[MAX_LEN_RES_NAME+1];
          AA1ToAA3(aa, rotType);
          if (aa == 'H') {
            StringArrayAppend(&designTypes, "HSD"); 
            StringArrayAppend(&patchTypes, "");
            StringArrayAppend(&designTypes, "HSE"); 
            StringArrayAppend(&patchTypes, "");
          }
          else if (aa != 'X') { 
            StringArrayAppend(&designTypes, rotType); 
            StringArrayAppend(&patchTypes, ""); 
          }
        }
      }
    }
    else if (ResidueGetDesignType(pSite->pRes) == Type_DesType_Repackable) {
      if (IsResidueHistidine(ResidueGetName(pSite->pRes))) {
        StringArrayAppend(&designTypes, "HSD"); 
        StringArrayAppend(&patchTypes, "");
        StringArrayAppend(&designTypes, "HSE"); 
        StringArrayAppend(&patchTypes, "");
      }
      else {
        StringArrayAppend(&designTypes, ResidueGetName(pSite->pRes)); 
        StringArrayAppend(&patchTypes, "");
      }
    }
    RotamerSetOfProteinGenerateByBBdepRotLib(DesignSiteGetRotamers(pSite), pSite->pRes, &designTypes, &patchTypes, rotlib, atomParams, resiTopos);
    StringArrayDestroy(&patchTypes);
    StringArrayDestroy(&designTypes);
  }

  return Success;
}


int StructureBuildAllRotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
    Chain* pChnI = StructureGetChain(pStruct, i);
    if (ChainGetType(pChnI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChnI)) == NULL) continue;
    for (int j = 0; j < ChainGetResidueNum(pChnI); j++) {
      Residue* pRes = ChainGetResidue(pChnI, j);
      if (pRes->desType != Type_DesType_Fixed) continue;
      if (FLAG_WILDTYPE_ONLY == TRUE) ResidueSetDesignType(pRes, Type_DesType_Repackable);
      else ResidueSetDesignType(pRes, Type_DesType_Mutable);
      ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resiTopos);
      if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE) ProteinSiteBuildNativeRotamer(pStruct, i, j, resiTopos);
      if (FLAG_ROTATE_HYDROXYL == TRUE) ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resiTopos);
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
    Chain* pChnI = StructureGetChain(pStruct, i);
    if (ChainGetType(pChnI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChnI)) != NULL) continue;
    for (int j = 0; j < ChainGetResidueNum(pChnI); j++) {
      Residue* pRes = ChainGetResidue(pChnI, j);
      if (pRes->desType != Type_DesType_Fixed) continue;
      BOOL isInteract = FALSE;
      for (int k = 0; k < StructureGetChainNum(pStruct); k++) {
        if (k == i) continue;
        Chain* pChainK = StructureGetChain(pStruct, k);
        if (ChainGetType(pChainK) != Type_Chain_Protein) continue;
        if (strstr(DES_CHAINS, ChainGetName(pChainK)) == 0) continue;
        for (int s = 0;s < ChainGetResidueNum(pChainK);s++) {
          Residue* pResidueKS = ChainGetResidue(pChainK, s);
          if (AtomArrayCalcMinDistance(&pRes->atoms, &pResidueKS->atoms) < CUT_PPI_DIST_SHELL1) {
            isInteract = TRUE;
            break;
          }
        }
        if (isInteract == TRUE) break;
      }
      if (isInteract == TRUE) {
        ResidueSetDesignType(pRes, Type_DesType_Repackable);
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resiTopos);
        if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE) ProteinSiteBuildNativeRotamer(pStruct, i, j, resiTopos);
        if (FLAG_ROTATE_HYDROXYL == TRUE) ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resiTopos);
      }
    }
  }

  return Success;
}


int StructureBuildPPIRotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resTopos)
{
  if (StructureGetChainNum(pStruct) < 2) {
    printf("there is only one chain in the whole structure, no protein-protein interface found\n");
    exit(ValueError);
  }

  IntArray* flagMutable = (IntArray*)malloc(sizeof(IntArray) * StructureGetChainNum(pStruct));
  IntArray* flagRepackable = (IntArray*)malloc(sizeof(IntArray) * StructureGetChainNum(pStruct));
  for (int i = 0;i < StructureGetChainNum(pStruct);i++) {
    Chain* pChnI = StructureGetChain(pStruct, i);
    IntArrayCreate(&flagMutable[i], pChnI->resNum);
    IntArrayCreate(&flagRepackable[i], pChnI->resNum);
    for (int j = 0; j < ChainGetResidueNum(pChnI); j++) {
      IntArraySet(&flagMutable[i], j, 0);
      IntArraySet(&flagRepackable[i], j, 0);
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
    Chain* pChnI = StructureGetChain(pStruct, i);
    for (int k = i + 1; k < StructureGetChainNum(pStruct); k++) {
      Chain* pChnK = StructureGetChain(pStruct, k);
      for (int j = 0; j < ChainGetResidueNum(pChnI); j++) {
        Residue* pResIJ = ChainGetResidue(pChnI, j);
        for (int s = 0; s < ChainGetResidueNum(pChnK); s++) {
          Residue* pResKS = ChainGetResidue(pChnK, s);
          if (AtomArrayCalcMinDistance(&pResIJ->atoms, &pResKS->atoms) < CUT_PPI_DIST_SHELL1) {
            if (strstr(DES_CHAINS, ChainGetName(pChnI)) != NULL) IntArraySet(&flagMutable[i], j, 1);
            else IntArraySet(&flagRepackable[i], j, 1);
            if (strstr(DES_CHAINS, ChainGetName(pChnK)) != NULL) IntArraySet(&flagMutable[k], s, 1);
            else IntArraySet(&flagRepackable[k], s, 1);
          }
          else if (AtomArrayCalcMinDistance(&pResIJ->atoms, &pResKS->atoms) < CUT_PPI_DIST_SHELL2) {
            if (strstr(DES_CHAINS, ChainGetName(pChnI)) != NULL) IntArraySet(&flagRepackable[i], j, 1);
            if (strstr(DES_CHAINS, ChainGetName(pChnK)) != NULL) IntArraySet(&flagRepackable[k], s, 1);
          }
        }
      }
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
    Chain* pChnI = StructureGetChain(pStruct, i);
    if (ChainGetType(pChnI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChnI)) != NULL) {
      for (int j = 0; j < ChainGetResidueNum(pChnI); j++) {
        Residue* pRes = ChainGetResidue(pChnI, j);
        if (pRes->desType != Type_DesType_Fixed) continue;
        if (IntArrayGet(&flagMutable[i], j) == 1) {
          if (FLAG_WILDTYPE_ONLY == TRUE) ResidueSetDesignType(pRes, Type_DesType_Repackable);
          else ResidueSetDesignType(pRes, Type_DesType_Mutable);
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resTopos);
          if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE) ProteinSiteBuildNativeRotamer(pStruct, i, j, resTopos);
          if (FLAG_ROTATE_HYDROXYL == TRUE) ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resTopos);
        }
        else if (IntArrayGet(&flagRepackable[i], j) == 1) {
          ResidueSetDesignType(pRes, Type_DesType_Repackable);
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resTopos);
          if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE) ProteinSiteBuildNativeRotamer(pStruct, i, j, resTopos);
          if (FLAG_ROTATE_HYDROXYL == TRUE) ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resTopos);
        }
      }
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
    Chain* pChnI = StructureGetChain(pStruct, i);
    if (ChainGetType(pChnI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChnI)) == NULL) {
      for (int j = 0; j < ChainGetResidueNum(pChnI); j++) {
        Residue* pRes = ChainGetResidue(pChnI, j);
        if (pRes->desType != Type_DesType_Fixed) continue;
        if (IntArrayGet(&flagRepackable[i], j) == 1) {
          ResidueSetDesignType(pRes, Type_DesType_Repackable);
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resTopos);
          if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE) ProteinSiteBuildNativeRotamer(pStruct, i, j, resTopos);
          if (FLAG_ROTATE_HYDROXYL == TRUE) ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resTopos);
        }
      }
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
    IntArrayDestroy(&flagMutable[i]);
    IntArrayDestroy(&flagRepackable[i]);
  }
  free(flagMutable);
  flagMutable = NULL;
  free(flagRepackable);
  flagRepackable = NULL;

  return Success;
}


int StructureBuildResfileRotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* pLib, AtomParams* pParams, ResTopos* pTopos, char* resfile)
{
  FileReader fr;
  StringArray strings;
  StringArrayCreate(&strings);
  if (!FAILED(FileReaderCreate(&fr, resfile))) {
    char buffer[MAX_LEN_ONE_LINE_CONTENT+1];
    while (!FAILED(FileReaderGetNextLine(&fr, buffer))) {
      char keyword[MAX_LEN_ONE_LINE_CONTENT+1];
      ExtractFirstStringFromSourceString(keyword, buffer);
      if (strcmp(keyword, "SITES_CATALYTIC_START") == 0) {
        while (!FAILED(FileReaderGetNextLine(&fr, buffer))) {
          StringArraySplitString(&strings, buffer, ' ');
          if (strcmp(StringArrayGet(&strings, 0), "SITES_CATALYTIC_END") == 0) break;
          char* chnName = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStruct, chnName);
          int resNdx = -1;
          ChainFindResidueIndexByPosInChain(pChain, posInChain, &resNdx);
          if (resNdx != -1) {
            if (StringArrayGetNum(&strings) == 3) {
              if (strcmp(StringArrayGet(&strings, 2), "NATROT") == 0) 
                ResidueSetDesignType(ChainGetResidue(pChain, resNdx), Type_DesType_NatRot);
              else {
                ResidueSetDesignType(ChainGetResidue(pChain, resNdx), Type_DesType_Catalytic);
                strcpy(ChainGetResidue(pChain, resNdx)->desAATypes, StringArrayGet(&strings, 2));
              }
            }
          }
        }
      }
      else if (strcmp(keyword, "SITES_DESIGN_START") == 0) {
        while (!FAILED(FileReaderGetNextLine(&fr, buffer))) {
          StringArraySplitString(&strings, buffer, ' ');
          if (strcmp(StringArrayGet(&strings, 0), "SITES_DESIGN_END") == 0) break;
          char* chnName = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStruct, chnName);
          int resNdx = -1;
          ChainFindResidueIndexByPosInChain(pChain, posInChain, &resNdx);
          if (resNdx != -1) {
            ResidueSetDesignType(ChainGetResidue(pChain, resNdx), Type_DesType_Mutable);
            if (StringArrayGetNum(&strings) == 3) strcpy(ChainGetResidue(pChain, resNdx)->desAATypes, StringArrayGet(&strings, 2));
          }
        }
      }
      else if (strcmp(keyword, "SITES_REPACK_START") == 0) {
        while (!FAILED(FileReaderGetNextLine(&fr, buffer))) {
          StringArraySplitString(&strings, buffer, ' ');
          if (strcmp(StringArrayGet(&strings, 0), "SITES_REPACK_END") == 0) break;
          char* chnName = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStruct, chnName);
          int resNdx = -1;
          ChainFindResidueIndexByPosInChain(pChain, posInChain, &resNdx);
          if (resNdx != -1) {
            ResidueSetDesignType(ChainGetResidue(pChain, resNdx), Type_DesType_Repackable);
            if (StringArrayGetNum(&strings) == 3) strcpy(ChainGetResidue(pChain, resNdx)->desAATypes, StringArrayGet(&strings, 2));
          }
        }
      }
      else if (strcmp(keyword, "SITES_FIX_START") == 0) {
        while (!FAILED(FileReaderGetNextLine(&fr, buffer))) {
          StringArraySplitString(&strings, buffer, ' ');
          if (strcmp(StringArrayGet(&strings, 0), "SITES_FIX_END") == 0) break;
          char* chnName = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStruct, chnName);
          int resNdx = -1;
          ChainFindResidueIndexByPosInChain(pChain, posInChain, &resNdx);
          if (resNdx != -1) ResidueSetDesignType(ChainGetResidue(pChain, resNdx), Type_DesType_Fixed);
        }
      }
    }
  }
  FileReaderDestroy(&fr);
  StringArrayDestroy(&strings);

  for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
    Chain* pChain = StructureGetChain(pStruct, i);
    if (ChainGetType(pChain) != Type_Chain_Protein) continue;
    for (int j = 0; j < ChainGetResidueNum(pChain);j++) {
      Residue* pRes = ChainGetResidue(pChain, j);
      if (ResidueGetDesignType(pRes) == Type_DesType_Mutable && FLAG_WILDTYPE_ONLY == TRUE) 
        ResidueSetDesignType(pRes, Type_DesType_Repackable);
      if (ResidueGetDesignType(pRes) == Type_DesType_Repackable && (strcmp(ResidueGetName(pRes), "ALA") == 0 || strcmp(ResidueGetName(pRes), "GLY") == 0)) 
        ResidueSetDesignType(pRes, Type_DesType_Fixed);
      if (ResidueGetDesignType(pRes) == Type_DesType_Fixed) 
        continue; // skip if current site is fixed
      ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, pLib, pParams, pTopos);
      if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE) 
        ProteinSiteBuildNativeRotamer(pStruct, i, j, pTopos);
      if (FLAG_ROTATE_HYDROXYL == TRUE) 
        ProteinSiteExpandHydroxylRotamers(pStruct, i, j, pTopos);
    }
  }

  return Success;
}

int StructureBuildCatalyticRotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resTopos, char* resfile)
{
  FileReader fr;
  if (!FAILED(FileReaderCreate(&fr, resfile))) {
    BOOL flagCataSite = FALSE;
    char buffer[MAX_LEN_ONE_LINE_CONTENT+1];
    while (!FAILED(FileReaderGetNextLine(&fr, buffer))) {
      char keyword[MAX_LEN_ONE_LINE_CONTENT+1];
      ExtractFirstStringFromSourceString(keyword, buffer);
      if (strcmp(keyword, "SITES_CATALYTIC_START") == 0) {
        flagCataSite = TRUE;
        StringArray strings;
        StringArrayCreate(&strings);
        while (!FAILED(FileReaderGetNextLine(&fr, buffer))) {
          StringArraySplitString(&strings, buffer, ' ');
          if (strcmp(StringArrayGet(&strings, 0), "SITES_CATALYTIC_END") == 0) 
            break;
          char* chnName = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStruct, chnName);
          int resNdx = -1;
          ChainFindResidueIndexByPosInChain(pChain, posInChain, &resNdx);
          if (resNdx != -1) {
            if (StringArrayGetNum(&strings) == 3) {
              if (strcmp(StringArrayGet(&strings, 2), "NATROT") == 0) 
                ResidueSetDesignType(ChainGetResidue(pChain, resNdx), Type_DesType_NatRot);
              else {
                ResidueSetDesignType(ChainGetResidue(pChain, resNdx), Type_DesType_Catalytic);
                strcpy(ChainGetResidue(pChain, resNdx)->desAATypes, StringArrayGet(&strings, 2));
              }
            }
          }
        }
        StringArrayDestroy(&strings);
      }
    }
    FileReaderDestroy(&fr);

    if (flagCataSite == FALSE) {
      char errMsg[MAX_LEN_ERR_MSG+1];
      sprintf(errMsg, "in file %s line %d, catalytic sites were not designated. Note: user should specify catalytic sites using keyword 'SITES_CATALYTIC_START' "
        "and 'SITES_CATALYTIC_END' in the resfile; in most cases, a catalytic site is not designable", __FILE__, __LINE__);
      TraceError(errMsg, FormatError);
      return FormatError;
    }

    for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
      Chain* pChain = StructureGetChain(pStruct, i);
      if (ChainGetType(pChain) != Type_Chain_Protein) continue;
      for (int j = 0; j < ChainGetResidueNum(pChain); j++) {
        Residue* pRes = ChainGetResidue(pChain, j);
        if (ResidueGetDesignType(pRes) == Type_DesType_NatRot || ResidueGetDesignType(pRes) == Type_DesType_Catalytic) {
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resTopos);
          if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE) 
            ProteinSiteBuildNativeRotamer(pStruct, i, j, resTopos);
          if (FLAG_ROTATE_HYDROXYL == TRUE) 
            ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resTopos);
        }        
      }
    }
  }
  else {
    char errMsg[MAX_LEN_ERR_MSG+1];
    sprintf(errMsg, "in file %s line %d, catalytic sites were not designated. Note: user should specify catalytic sites using keyword 'SITES_CATALYTIC_START' "
      "and 'SITES_CATALYTIC_END' in the resfile; in most cases, a catalytic site is not designable", __FILE__, __LINE__);
    TraceError(errMsg, FormatError);
    return FormatError;
  }

  return Success;
}


int StructureBuildPLIShell1RotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos, char* resfile)
{
  FileReader fr;
  int mutSiteNum = 0;
  if (!FAILED(FileReaderCreate(&fr, resfile))) {
    char buffer[MAX_LEN_ONE_LINE_CONTENT+1];
    while (!FAILED(FileReaderGetNextLine(&fr, buffer))) {
      char keyword[MAX_LEN_ONE_LINE_CONTENT+1];
      ExtractFirstStringFromSourceString(keyword, buffer);
      if (strcmp(keyword, "SITES_DESIGN_START") == 0) {
        StringArray strings;
        StringArrayCreate(&strings);
        while (!FAILED(FileReaderGetNextLine(&fr, buffer))) {
          StringArraySplitString(&strings, buffer, ' ');
          if (strcmp(StringArrayGet(&strings, 0), "SITES_DESIGN_END") == 0) break;
          char* chain = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStruct, chain);
          int resNdx = -1;
          ChainFindResidueIndexByPosInChain(pChain, posInChain, &resNdx);
          if (resNdx != -1) {
            Residue* pRes = ChainGetResidue(pChain, resNdx);
            ResidueSetDesignType(pRes, Type_DesType_Mutable);
            if (StringArrayGetNum(&strings) == 3) {
              strcpy(pRes->desAATypes, StringArrayGet(&strings, 2));
              if ((strcmp(ResidueGetName(pRes), "ALA") == 0 && strlen(pRes->desAATypes) == 1 && pRes->desAATypes[0] == 'A')
                || (strcmp(ResidueGetName(pRes), "GLY") == 0 && strlen(pRes->desAATypes) == 1 && pRes->desAATypes[0] == 'G')) {
                ResidueSetDesignType(pRes, Type_DesType_Fixed);
                continue;
              }
            }
          }
          mutSiteNum++;
        }
        StringArrayDestroy(&strings);
      }
    }
    FileReaderDestroy(&fr);
  }

  if (mutSiteNum > 0) { // Resfile has specified at least one mutable site; build rotamers for each of these mutable sites
    for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
      Chain* pChain = StructureGetChain(pStruct, i);
      if (ChainGetType(pChain) != Type_Chain_Protein) continue;
      for (int j = 0; j < ChainGetResidueNum(pChain); j++) {
        Residue* pRes = ChainGetResidue(pChain, j);
        if (ResidueGetDesignType(pRes) == Type_DesType_Mutable) {
          if (FLAG_WILDTYPE_ONLY == TRUE) {
            if (strcmp(ResidueGetName(pRes), "ALA") == 0 || strcmp(ResidueGetName(pRes), "GLY") == 0) {
              ResidueSetDesignType(pRes, Type_DesType_Fixed);
              continue;
            }
            else ResidueSetDesignType(pRes, Type_DesType_Repackable);
          }
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resiTopos);
          if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE) ProteinSiteBuildNativeRotamer(pStruct, i, j, resiTopos);
          if (FLAG_ROTATE_HYDROXYL == TRUE) ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resiTopos);
        }
      }
    }
  }
  else { // Resfile did not specify any mutable site; thus find mutable sites based on PLI distance and build rotamers for such
    char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
    Residue* pSmallMol = NULL;
    StructureFindSmallMol(pStruct, &pSmallMol);
    if (pSmallMol == NULL) {
      int result = DataNotExistError;
      sprintf(errMsg, "in file %s line %d, cannot find small molecule", __FILE__, __LINE__);
      TraceError(errMsg, result);
      return result;
    }
    for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
      Chain* pChain = StructureGetChain(pStruct, i);
      if (strcmp(ResidueGetChainName(pSmallMol), ChainGetName(pChain)) == 0) continue;
      if (ChainGetType(pChain) != Type_Chain_Protein) continue;
      for (int j = 0; j < ChainGetResidueNum(pChain); j++) {
        Residue* pRes = ChainGetResidue(pChain, j);
        double minDist = AtomArrayCalcMinDistance(&pSmallMol->atoms, &pRes->atoms);
        if (minDist > CUT_PLI_DIST_SHELL1) continue;
        if (ResidueGetDesignType(pRes) == Type_DesType_Fixed) {
          ResidueSetDesignType(pRes, Type_DesType_Mutable);
          if (FLAG_WILDTYPE_ONLY == TRUE) {
            if (strcmp(ResidueGetName(pRes), "ALA") == 0 || strcmp(ResidueGetName(pRes), "GLY") == 0) {
              ResidueSetDesignType(pRes, Type_DesType_Fixed);
              continue;
            }
            else ResidueSetDesignType(pRes, Type_DesType_Repackable);
          }
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resiTopos);
          if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE) ProteinSiteBuildNativeRotamer(pStruct, i, j, resiTopos);
          if (FLAG_ROTATE_HYDROXYL == TRUE) ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resiTopos);
        }
      }
    }
  }

  return Success;
}


int StructureBuildPLIShell2RotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos, char* resfile)
{
  FileReader fr;
  int repackSiteNum = 0;
  if (!FAILED(FileReaderCreate(&fr, resfile))) {
    char buffer[MAX_LEN_ONE_LINE_CONTENT+1];
    while (!FAILED(FileReaderGetNextLine(&fr, buffer))) {
      char keyword[MAX_LEN_ONE_LINE_CONTENT+1];
      ExtractFirstStringFromSourceString(keyword, buffer);
      if (strcmp(keyword, "SITES_REPACK_START") == 0) {
        StringArray strs;
        StringArrayCreate(&strs);
        while (!FAILED(FileReaderGetNextLine(&fr, buffer))) {
          StringArraySplitString(&strs, buffer, ' ');
          if (strcmp(StringArrayGet(&strs, 0), "SITES_REPACK_END") == 0) break;
          char* chnname = StringArrayGet(&strs, 0);
          int posInChain = atoi(StringArrayGet(&strs, 1));
          Chain* pChain = StructureFindChainByName(pStruct, chnname);
          int resNdx = -1;
          ChainFindResidueIndexByPosInChain(pChain, posInChain, &resNdx);
          if (resNdx != -1) ChainGetResidue(pChain, resNdx)->desType = Type_DesType_Repackable;
          repackSiteNum++;
        }
        StringArrayDestroy(&strs);
      }
    }
    FileReaderDestroy(&fr);
  }

  if (repackSiteNum > 0) { // Resfile has specified at least one mutable site; build rotamers for each of these mutable sites
    for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
      Chain* pChain = StructureGetChain(pStruct, i);
      if (ChainGetType(pChain) != Type_Chain_Protein) continue;
      for (int j = 0; j < ChainGetResidueNum(pChain); j++) {
        Residue* pRes = ChainGetResidue(pChain, j);
        if (ResidueGetDesignType(pRes) == Type_DesType_Repackable) {
          if (strcmp(ResidueGetName(pRes), "ALA") == 0 || strcmp(ResidueGetName(pRes), "GLY") == 0) {
            ResidueSetDesignType(pRes, Type_DesType_Fixed);
            continue;
          }
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resiTopos);
          if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE) ProteinSiteBuildNativeRotamer(pStruct, i, j, resiTopos);
          if (FLAG_ROTATE_HYDROXYL == TRUE) ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resiTopos);
        }
      }
    }
  }
  else {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
    Residue* pSmallMol = NULL;
    StructureFindSmallMol(pStruct, &pSmallMol);
    if (pSmallMol == NULL) {
      int result = DataNotExistError;
      sprintf(errMsg, "in file %s line %d, cannot find small molecule", __FILE__, __LINE__);
      TraceError(errMsg, result);
      return result;
    }
    for (int i = 0; i < StructureGetChainNum(pStruct); i++) {
      Chain* pChainI = StructureGetChain(pStruct, i);
      if (strcmp(ResidueGetChainName(pSmallMol), ChainGetName(pChainI)) == 0) continue;
      if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
      for (int j = 0; j < ChainGetResidueNum(pChainI); j++) {
        Residue* pRes = ChainGetResidue(pChainI, j);
        if (AtomArrayCalcMinDistance(&pSmallMol->atoms, &pRes->atoms) > CUT_PLI_DIST_SHELL2) continue;
        if (ResidueGetDesignType(pRes) == Type_DesType_Fixed) {
          ResidueSetDesignType(pRes, Type_DesType_Repackable);
          if (strcmp(ResidueGetName(pRes), "ALA") == 0 || strcmp(ResidueGetName(pRes), "GLY") == 0) {
            ResidueSetDesignType(pRes, Type_DesType_Fixed);
            continue;
          }
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resiTopos);
          if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE) ProteinSiteBuildNativeRotamer(pStruct, i, j, resiTopos);
          if (FLAG_ROTATE_HYDROXYL == TRUE) ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resiTopos);
        }
      }
    }
  }

  return Success;
}


BOOL IsNativeRotamerInBBdepRotLib(Structure* pStruct, int chnNdx, int resNdx, ResTopos* resiTopos, BBdepRotLib* rotlib, double cut)
{
  char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
  Chain* pChn = StructureGetChain(pStruct, chnNdx);
  Residue* pRes = ChainGetResidue(pChn, resNdx);
  if (pChn->type == Type_Chain_Protein)
  {
    if (strcmp(ResidueGetName(pRes), "ALA") == 0 || strcmp(ResidueGetName(pRes), "GLY") == 0)
    {
      return TRUE;
    }

    // 1. calculate all sidechain torsions for current residue
    double DUNBRACK_ENERGY = 0.0;
    int torNum = 0;
    char resName[MAX_LEN_RES_NAME+1];
    strcpy(resName, ResidueGetName(pRes));
    if (      strcmp(resName, "ALA") == 0 ) { torNum = 0; }
    else if ( strcmp(resName, "ARG") == 0 ) { torNum = 4; }
    else if ( strcmp(resName, "ASN") == 0 ) { torNum = 2; }
    else if ( strcmp(resName, "ASP") == 0 ) { torNum = 2; }
    else if ( strcmp(resName, "CYS") == 0 ) { torNum = 1; }
    else if ( strcmp(resName, "GLN") == 0 ) { torNum = 3; }
    else if ( strcmp(resName, "GLU") == 0 ) { torNum = 3; }
    else if ( strcmp(resName, "GLY") == 0 ) { torNum = 0; }
    else if ( strcmp(resName, "HSD") == 0 ) { torNum = 2; }
    else if ( strcmp(resName, "HSE") == 0 ) { torNum = 2; }
    else if ( strcmp(resName, "ILE") == 0 ) { torNum = 2; }
    else if ( strcmp(resName, "LEU") == 0 ) { torNum = 2; }
    else if ( strcmp(resName, "LYS") == 0 ) { torNum = 4; }
    else if ( strcmp(resName, "MET") == 0 ) { torNum = 3; }
    else if ( strcmp(resName, "PHE") == 0 ) { torNum = 2; }
    else if ( strcmp(resName, "PRO") == 0 ) { torNum = 2; }
    else if ( strcmp(resName, "SER") == 0 ) { torNum = 1; }
    else if ( strcmp(resName, "THR") == 0 ) { torNum = 1; }
    else if ( strcmp(resName, "TRP") == 0 ) { torNum = 2; }
    else if ( strcmp(resName, "TYR") == 0 ) { torNum = 2; }
    else if ( strcmp(resName, "VAL") == 0 ) { torNum = 1; }

    DoubleArray Xs;
    DoubleArrayCreate(&Xs, 0);
    ResTopology topo;
    ResTopologyCreate(&topo);
    ResToposGet(resiTopos, resName, &topo);
    for (int ndx = 0; ndx < torNum; ndx++)
    {
      Type_ProteinAtomOrder desiredAtomBOrder = Type_ProteinAtomOrder_FromInt(ndx);
      Type_ProteinAtomOrder desiredAtomCOrder = Type_ProteinAtomOrder_FromInt(ndx + 1);
      CharmmIC curTorIC;
      CharmmICCreate(&curTorIC);
      BOOL icFound = FALSE;
      for (int icNdx = 0; icNdx < ResTopologyGetCharmmICNum(&topo); icNdx++)
      {
        Type_ProteinAtomOrder atomBOrder;
        Type_ProteinAtomOrder atomCOrder;
        ResTopologyGetCharmmIC(&topo, icNdx, &curTorIC);
        atomBOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomB(&curTorIC));
        atomCOrder = Type_ProteinAtomOrder_JudgedByAtomName(CharmmICGetAtomC(&curTorIC));
        if (desiredAtomBOrder == atomBOrder && desiredAtomCOrder == atomCOrder)
        {
          icFound = TRUE;
          break;
        }
      }

      if (icFound == FALSE)
      {
        sprintf(errMsg, "in file %s line %d, cannot find IC of the %dth torsion for residue %s", __FILE__, __LINE__, ndx + 1, resName);
        TraceError(errMsg, DataNotExistError);
        CharmmICDestroy(&curTorIC);
        return DataNotExistError;
      }

      Atom* pAtomA = ResidueGetAtomByName(pRes, CharmmICGetAtomA(&curTorIC));
      Atom* pAtomB = ResidueGetAtomByName(pRes, CharmmICGetAtomB(&curTorIC));
      Atom* pAtomC = ResidueGetAtomByName(pRes, CharmmICGetAtomC(&curTorIC));
      Atom* pAtomD = ResidueGetAtomByName(pRes, CharmmICGetAtomD(&curTorIC));
      DoubleArrayAppend(&Xs, GetTorsionAngle(&pAtomA->xyz, &pAtomB->xyz, &pAtomC->xyz, &pAtomD->xyz));

      CharmmICDestroy(&curTorIC);
    }
    ResTopologyDestroy(&topo);

    // 2. get the phi, psi angles for current residue
    int binNdx = ((int)(pRes->phipsi[0] + 180) / 10) * 36 + (int)(pRes->phipsi[1] + 180) / 10;
    RotLibPhiPsi* pRotLibPhiPsi = &rotlib->pps[binNdx];
    int rtNdx = -1, matchNdx = -1;
    StringArrayFind(&pRotLibPhiPsi->rotTypes, ResidueGetName(pRes), &rtNdx);
    DoubleArray* pTorsionsArrayForTypeI = pRotLibPhiPsi->tors[rtNdx];
    DoubleArray* pDeviationsArrayForTypeI = pRotLibPhiPsi->devs[rtNdx];
    for (int i = 0;i < IntArrayGet(&pRotLibPhiPsi->rotNums, rtNdx);i++)
    {
      DoubleArray* pTors = &pTorsionsArrayForTypeI[i];
      DoubleArray* pDevs = &pDeviationsArrayForTypeI[i];
      BOOL isMatched = TRUE;
      for (int j = 0; j < DoubleArrayGetLength(&Xs); j++)
      {
        double min = DoubleArrayGet(pTors, j) - DegToRad(cut);
        double max = DoubleArrayGet(pTors, j) + DegToRad(cut);
        double tor = DoubleArrayGet(&Xs, j);
        double m2pi = tor - 2 * PI;
        double p2pi = tor + 2 * PI;
        double tor2 = tor;
        if ( (strcmp(resName, "PHE") == 0 && j == 1) ||
             (strcmp(resName, "TYR") == 0 && j == 1) ||
             (strcmp(resName, "ASP") == 0 && j == 1) ||
             (strcmp(resName, "GLU") == 0 && j == 2) )
        {
          tor2 = tor + PI;
          tor2 = tor > 0 ? tor - PI : tor2;
        }
        double tor2m2pi = tor2 - 2 * PI;
        double tor2p2pi = tor2 + 2 * PI;
        if ( !( (tor <= max && tor >= min) ||
                (m2pi <= max && m2pi >= min) ||
                (p2pi <= max && p2pi >= min) ||
                (tor2 <= max && tor2 >= min) ||
                (tor2m2pi <= max && tor2m2pi >= min) ||
                (tor2p2pi <= max && tor2p2pi >= min)
             ) )
        {
          isMatched = FALSE;
          break;
        }
      }
      if (isMatched == TRUE)
      {
        matchNdx = i;
        break;
      }
    }
    DoubleArrayDestroy(&Xs);

    if (matchNdx != -1)
    {
      return TRUE;
    }
    return FALSE;
  }

  return FALSE;
}

int ProteinSiteWriteRotamers(Structure* pStruct, int chnNdx, int resNdx, char* file)
{
  DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  Rotamers* pSet = DesignSiteGetRotamers(pSite);
  FILE* pOut = fopen(file, "w");
  if (pOut != NULL)
  {
    for (int i = 0; i < RotamersGetNum(pSet); i++)
    {
      Rotamer* pRotamer = RotamersGet(pSet, i);
      RotamerRestore(pRotamer, pSet);
      Model(i, pOut);
      RotamerShowInPDBFormat(pRotamer, "ATOM", RotamerGetChainName(pRotamer), 1, i, pOut);
      EndModel(pOut);
      RotamerExtract(pRotamer);
    }
    fclose(pOut);
  }

  return Success;
}


int ProteinSiteBuildNativeRotamer(Structure* pStruct, int chnNdx, int resNdx, ResTopos* resiTopos)
{
  Chain* pChn = StructureGetChain(pStruct, chnNdx);
  Residue* pRes = ChainGetResidue(pChn, resNdx);

  if (ChainGetType(pChn) == Type_Chain_Protein) {
    DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
    if (pSite == NULL) {
      (pStruct->desSiteNum)++;
      pStruct->desSites = (DesignSite*)realloc(pStruct->desSites, sizeof(DesignSite) * pStruct->desSiteNum);
      DesignSiteCreate(&pStruct->desSites[pStruct->desSiteNum - 1]);
      pSite = StructureGetDesignSite(pStruct, pStruct->desSiteNum - 1);
      pSite->pRes = pRes;
      pSite->chnNdx = chnNdx;
      pSite->resNdx = resNdx;
    }

    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    Rotamer rot;
    RotamerCreate(&rot);
    RotamerSetType(&rot, ResidueGetName(pSite->pRes));
    RotamerSetChainName(&rot, ResidueGetChainName(pSite->pRes));
    RotamerSetPosInChain(&rot, ResidueGetPosInChain(pSite->pRes));
    Rotamer* pRep = RotamersFindRep(pSet, RotamerGetType(&rot));
    if (pRep != NULL) {
      AtomArrayCopy(&rot.atoms, &pRep->atoms);
      for (int i = 0; i < RotamerGetAtomNum(&rot); i++) {
        Atom* pAtom = RotamerGetAtom(&rot, i);
        pAtom->xyz = ResidueGetAtomByName(pSite->pRes, AtomGetName(pAtom))->xyz;
      }
      BondsCopy(&rot.bonds, &pRep->bonds);
      XYZArrayResize(&rot.xyzs, RotamerGetAtomNum(pRep));
      for (int i = 0; i < XYZArrayGetNum(&rot.xyzs); i++) XYZArraySet(&rot.xyzs, i, &AtomArrayGet(&rot.atoms, i)->xyz);
      RotamerSetDunbrack(&rot, ResidueGetDunbrack(pRes));
      DoubleArrayCopy(&rot.Xs, &pRes->Xs);
      RotamersAdd(&pSite->rots, &rot);
    }
    else
    { // cannot find a representative rotamer
      if (ResidueGetDesignType(pRes) == Type_DesType_NatRot) {
        AtomArrayCopy(&rot.atoms, &pRes->atoms);
        BondsCopy(&rot.bonds, &pRes->bonds);
        XYZArrayResize(&rot.xyzs, AtomArrayGetNum(&rot.atoms));
        for (int i = 0; i < XYZArrayGetNum(&rot.xyzs); i++) XYZArraySet(&rot.xyzs, i, &AtomArrayGet(&rot.atoms, i)->xyz);
        RotamerSetDunbrack(&rot, ResidueGetDunbrack(pRes));
        DoubleArrayCopy(&rot.Xs, &pRes->Xs);
        RotamersAdd(&pSite->rots, &rot);
      }
    }
    RotamerDestroy(&rot);

    if (strcmp(ResidueGetName(pSite->pRes), "HSD") == 0) {
      if ((pRep = RotamersFindRep(pSet, "HSE")) != NULL) {
        Residue resi;
        ResidueCreate(&resi);
        ResidueSetName(&resi, "HSE");
        AtomArrayCopy(&resi.atoms, &pRep->atoms);
        for (int i = 0; i < ResidueGetAtomNum(&resi); i++) {
          Atom* pAtom = ResidueGetAtom(&resi, i);
          if (pAtom->isBBAtom == FALSE && AtomIsHydrogen(pAtom) == TRUE) {
            pAtom->isXyzValid = FALSE;
            continue;
          }
          pAtom->xyz = ResidueGetAtomByName(pSite->pRes, AtomGetName(pAtom))->xyz;
        }
        ResidueCalcAllAtomXYZ(&resi, resiTopos, NULL, NULL);

        RotamerCreate(&rot);
        RotamerSetType(&rot, "HSE");
        RotamerSetChainName(&rot, ResidueGetChainName(pSite->pRes));
        RotamerSetPosInChain(&rot, ResidueGetPosInChain(pSite->pRes));
        AtomArrayCopy(&rot.atoms, &pRep->atoms);
        for (int i = 0; i < RotamerGetAtomNum(&rot); i++) {
          Atom* pAtom = RotamerGetAtom(&rot, i);
          pAtom->xyz = ResidueGetAtomByName(&resi, AtomGetName(pAtom))->xyz;
        }
        BondsCopy(&rot.bonds, &pRep->bonds);

        XYZArrayResize(&rot.xyzs, RotamerGetAtomNum(pRep));
        for (int i = 0; i < XYZArrayGetNum(&rot.xyzs); i++) XYZArraySet(&rot.xyzs, i, &AtomArrayGet(&rot.atoms, i)->xyz);
        RotamerSetDunbrack(&rot, ResidueGetDunbrack(pRes));
        DoubleArrayCopy(&rot.Xs, &pRes->Xs);
        RotamersAdd(&pSite->rots, &rot);

        RotamerDestroy(&rot);
        ResidueDestroy(&resi);
      }
    }
    else if (strcmp(ResidueGetName(pSite->pRes), "HSE") == 0) {
      if ((pRep = RotamersFindRep(pSet, "HSD")) != NULL) {
        Residue resi;
        ResidueCreate(&resi);
        ResidueSetName(&resi, "HSD");
        AtomArrayCopy(&resi.atoms, &pRep->atoms);
        for (int i = 0; i < ResidueGetAtomNum(&resi); i++) {
          Atom* pAtom = ResidueGetAtom(&resi, i);
          if (pAtom->isBBAtom == FALSE && AtomIsHydrogen(pAtom) == TRUE) {
            pAtom->isXyzValid = FALSE;
            continue;
          }
          pAtom->xyz = ResidueGetAtomByName(pSite->pRes, AtomGetName(pAtom))->xyz;
        }
        ResidueCalcAllAtomXYZ(&resi, resiTopos, NULL, NULL);

        RotamerCreate(&rot);
        RotamerSetType(&rot, "HSD");
        RotamerSetChainName(&rot, ResidueGetChainName(pSite->pRes));
        RotamerSetPosInChain(&rot, ResidueGetPosInChain(pSite->pRes));

        AtomArrayCopy(&rot.atoms, &pRep->atoms);
        for (int i = 0; i < RotamerGetAtomNum(&rot); i++) {
          Atom* pAtom = RotamerGetAtom(&rot, i);
          pAtom->xyz = ResidueGetAtomByName(&resi, AtomGetName(pAtom))->xyz;
        }
        BondsCopy(&rot.bonds, &pRep->bonds);

        XYZArrayResize(&rot.xyzs, RotamerGetAtomNum(pRep));
        for (int i = 0; i < XYZArrayGetNum(&rot.xyzs); i++) XYZArraySet(&rot.xyzs, i, &AtomArrayGet(&rot.atoms, i)->xyz);
        RotamerSetDunbrack(&rot, ResidueGetDunbrack(pRes));
        DoubleArrayCopy(&rot.Xs, &pRes->Xs);
        RotamersAdd(&pSite->rots, &rot);
        RotamerDestroy(&rot);
        ResidueDestroy(&resi);
      }
    }
  }

  return Success;
}


int ProteinSiteBuildFlippedNativeRotamer(Structure* pStruct, int chnNdx, int resNdx, ResTopos* resiTopos)
{
  Chain* pChn = StructureGetChain(pStruct, chnNdx);
  Residue* pRes = ChainGetResidue(pChn, resNdx);
  if (pChn->type == Type_Chain_Protein)
  {
    if ( strcmp(ResidueGetName(pRes), "ASN") != 0 && 
         strcmp(ResidueGetName(pRes), "GLN") != 0 && 
         strcmp(ResidueGetName(pRes), "HSD") != 0 && 
         strcmp(ResidueGetName(pRes), "HSE") != 0 )
    {
      return Success;
    }

    DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    int rotNum = RotamersGetNum(pSet);
    for (int i = 0; i < rotNum; i++)
    {
      Rotamer* pRot = RotamersGet(pSet, i);
      Rotamer* pRepresentative = RotamersFindRep(pSet, RotamerGetType(pRot));
      Rotamer rot;
      RotamerCreate(&rot);
      RotamerCopy(&rot, pRepresentative);
      if (strcmp(RotamerGetType(pRot), "ASN") == 0)
      {
        int ndx1, ndx2;
        RotamerFindAtom(&rot, "ND2", &ndx1);
        RotamerFindAtom(&rot, "OD1", &ndx2);
        XYZ xyz1 = RotamerGetAtom(&rot, ndx1)->xyz;
        XYZ xyz2 = RotamerGetAtom(&rot, ndx2)->xyz;
        RotamerGetAtom(&rot, ndx1)->xyz = xyz2;
        RotamerGetAtom(&rot, ndx2)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz2);
        XYZArraySet(&rot.xyzs, ndx2, &xyz1);
        ResTopology topo;
        CharmmIC ic;
        ResTopologyCreate(&topo);
        CharmmICCreate(&ic);
        ResToposGet(resiTopos, "ASN", &topo);
        ResTopologyFindCharmmIC(&topo, "HD21", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&rot, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &xyz1);
        RotamerFindAtom(&rot, "HD21", &ndx1);
        RotamerGetAtom(&rot, ndx1)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz1);
        ResTopologyFindCharmmIC(&topo, "HD22", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&rot, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &xyz1);
        RotamerFindAtom(&rot, "HD22", &ndx1);
        RotamerGetAtom(&rot, ndx1)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz1);
        ResTopologyDestroy(&topo);
        CharmmICDestroy(&ic);
      }
      else if (strcmp(RotamerGetType(pRot), "GLN") == 0)
      {
        int ndx1, ndx2;
        RotamerFindAtom(&rot, "NE2", &ndx1);
        RotamerFindAtom(&rot, "OE1", &ndx2);
        XYZ xyz1 = RotamerGetAtom(&rot, ndx1)->xyz;
        XYZ xyz2 = RotamerGetAtom(&rot, ndx2)->xyz;
        RotamerGetAtom(&rot, ndx1)->xyz = xyz2;
        RotamerGetAtom(&rot, ndx2)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz2);
        XYZArraySet(&rot.xyzs, ndx2, &xyz1);
        ResTopology topo;
        CharmmIC ic;
        ResTopologyCreate(&topo);
        CharmmICCreate(&ic);
        ResToposGet(resiTopos, "GLN", &topo);
        ResTopologyFindCharmmIC(&topo, "HE21", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&rot, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &xyz1);
        RotamerFindAtom(&rot, "HE21", &ndx1);
        RotamerGetAtom(&rot, ndx1)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz1);
        ResTopologyFindCharmmIC(&topo, "HE22", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&rot, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &xyz1);
        RotamerFindAtom(&rot, "HE22", &ndx1);
        RotamerGetAtom(&rot, ndx1)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz1);
        ResTopologyDestroy(&topo);
        CharmmICDestroy(&ic);
      }
      else if (strcmp(RotamerGetType(pRot), "HSD") == 0)
      {
        int ndx1, ndx2;
        RotamerFindAtom(&rot, "CD2", &ndx1);
        RotamerFindAtom(&rot, "ND1", &ndx2);
        XYZ xyz1 = RotamerGetAtom(&rot, ndx1)->xyz;
        XYZ xyz2 = RotamerGetAtom(&rot, ndx2)->xyz;
        RotamerGetAtom(&rot, ndx1)->xyz = xyz2;
        RotamerGetAtom(&rot, ndx2)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz2);
        XYZArraySet(&rot.xyzs, ndx2, &xyz1);
        RotamerFindAtom(&rot, "NE2", &ndx1);
        RotamerFindAtom(&rot, "CE1", &ndx2);
        xyz1 = RotamerGetAtom(&rot, ndx1)->xyz;
        xyz2 = RotamerGetAtom(&rot, ndx2)->xyz;
        RotamerGetAtom(&rot, ndx1)->xyz = xyz2;
        RotamerGetAtom(&rot, ndx2)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz2);
        XYZArraySet(&rot.xyzs, ndx2, &xyz1);
        ResTopology topo;
        CharmmIC ic;
        ResTopologyCreate(&topo);
        CharmmICCreate(&ic);
        ResToposGet(resiTopos, "HSD", &topo);
        ResTopologyFindCharmmIC(&topo, "HD1", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&rot, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &xyz1);
        RotamerFindAtom(&rot, "HD1", &ndx1);
        RotamerGetAtom(&rot, ndx1)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz1);
        ResTopologyDestroy(&topo);
        CharmmICDestroy(&ic);
      }
      else if (strcmp(RotamerGetType(pRot), "HSE") == 0)
      {
        int ndx1, ndx2;
        RotamerFindAtom(&rot, "CD2", &ndx1);
        RotamerFindAtom(&rot, "ND1", &ndx2);
        XYZ xyz1 = RotamerGetAtom(&rot, ndx1)->xyz;
        XYZ xyz2 = RotamerGetAtom(&rot, ndx2)->xyz;
        RotamerGetAtom(&rot, ndx1)->xyz = xyz2;
        RotamerGetAtom(&rot, ndx2)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz2);
        XYZArraySet(&rot.xyzs, ndx2, &xyz1);
        RotamerFindAtom(&rot, "NE2", &ndx1);
        RotamerFindAtom(&rot, "CE1", &ndx2);
        xyz1 = RotamerGetAtom(&rot, ndx1)->xyz;
        xyz2 = RotamerGetAtom(&rot, ndx2)->xyz;
        RotamerGetAtom(&rot, ndx1)->xyz = xyz2;
        RotamerGetAtom(&rot, ndx2)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz2);
        XYZArraySet(&rot.xyzs, ndx2, &xyz1);
        ResTopology topo;
        CharmmIC ic;
        ResTopologyCreate(&topo);
        CharmmICCreate(&ic);
        ResToposGet(resiTopos, "HSE", &topo);
        ResTopologyFindCharmmIC(&topo, "HE2", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&rot, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&rot, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &xyz1);
        RotamerFindAtom(&rot, "HE2", &ndx1);
        RotamerGetAtom(&rot, ndx1)->xyz = xyz1;
        XYZArraySet(&rot.xyzs, ndx1, &xyz1);
        ResTopologyDestroy(&topo);
        CharmmICDestroy(&ic);
      }
      RotamersAdd(pSet, &rot);
      RotamerDestroy(&rot);
    }
  }

  return Success;
}


int ProteinSiteExpandHydroxylRotamers(Structure* pStruct, int chnNdx, int resNdx, ResTopos* resTopos)
{
  Chain* pChain = StructureGetChain(pStruct, chnNdx);
  if (pChain->type == Type_Chain_Protein)
  {
    DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    int rotNum = RotamersGetNum(pSet);
    ResTopology topo;
    CharmmIC ic;
    ResTopologyCreate(&topo);
    CharmmICCreate(&ic);
    if ( RotamersFindRep(pSet, "SER") != NULL)
    {
      ResToposGet(resTopos, "SER", &topo);
      ResTopologyFindCharmmIC(&topo, "HG", &ic);
      double icParaX = ic.icParam[2];
      int addNum = 0;
      Rotamer rot;
      RotamerCreate(&rot);
      for (int j = 0; j < rotNum; j++)
      {
        Rotamer* pRot = RotamersGet(pSet, j);
        if (strcmp(RotamerGetType(pRot), "SER") == 0)
        {
          RotamerRestore(pRot, pSet);
          int atNdx;
          RotamerFindAtom(pRot, "HG", &atNdx);
          RotamerCopy(&rot, pRot);
          RotamerExtract(pRot);
          RotamersAdd(pSet, &rot);
          for (int k = 0; k < EXPANDED_ROT_SER; k++)
          {
            ic.icParam[2] = icParaX + 2.0 * PI * (k + 1) / (EXPANDED_ROT_SER + 1);
            if (ic.icParam[2] > PI) ic.icParam[2] -= 2 * PI;
            GetFourthAtom(&RotamerGetAtomByName(&rot, "CA")->xyz,
              &RotamerGetAtomByName(&rot, "CB")->xyz,
              &RotamerGetAtomByName(&rot, "OG")->xyz,
              ic.icParam,
              &RotamerGetAtomByName(&rot, "HG")->xyz);
            XYZArraySet(&rot.xyzs, atNdx, &RotamerGetAtomByName(&rot, "HG")->xyz);
            // check if the rotamer should be added
            BOOL isAccepted = TRUE;
            for (int kk = 0; kk < RotamerGetAtomNum(&rot); kk++)
            {
              Atom* pAtomK = RotamerGetAtom(&rot, kk);
              if (strcmp(AtomGetName(pAtomK), "HG") != 0)
              {
                continue;
              }
              for (int ss = 0; ss < RotamerGetAtomNum(&rot); ss++)
              {
                Atom* pAtomS = RotamerGetAtom(&rot, ss);
                if (strcmp(AtomGetName(pAtomK), AtomGetName(pAtomS)) == 0)
                {
                  continue;
                }
                int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtomK), AtomGetName(pAtomS), RotamerGetBonds(&rot));
                if (bondType == 14 || bondType == 15)
                {
                  double distance = XYZDist(&pAtomK->xyz, &pAtomS->xyz);
                  //if(distance < (pAtomK->vdw_radius+pAtomS->vdw_radius)*0.75){
                  if (distance < 2.25)
                  {
                    isAccepted = FALSE;
                    break;
                  }
                }
              }
              if (isAccepted == FALSE)
              {
                break;
              }
            }
            if (isAccepted == TRUE)
            {
              RotamersAdd(pSet, &rot);
              addNum++;
            }
          }
        }
      }
      RotamerDestroy(&rot);
    }

    if (RotamersFindRep(pSet, "THR") != NULL)
    {
      ResToposGet(resTopos, "THR", &topo);
      ResTopologyFindCharmmIC(&topo, "HG1", &ic);
      double icParaX = ic.icParam[2];
      int addNum = 0;
      Rotamer rot;
      RotamerCreate(&rot);
      for (int j = 0; j < rotNum; j++)
      {
        Rotamer* pRot = RotamersGet(pSet, j);
        if (strcmp(RotamerGetType(pRot), "THR") == 0)
        {
          RotamerRestore(pRot, pSet);
          int atNdx;
          RotamerFindAtom(pRot, "HG1", &atNdx);
          RotamerCopy(&rot, pRot);
          RotamerExtract(pRot);
          for (int k = 0; k < EXPANDED_ROT_THR; k++)
          {
            BOOL isAccepted = TRUE;
            ic.icParam[2] = icParaX + 2.0 * PI * (k + 1) / (EXPANDED_ROT_THR + 1);
            if (ic.icParam[2] > PI) ic.icParam[2] -= 2 * PI;
            GetFourthAtom(&RotamerGetAtomByName(&rot, "CA")->xyz,
              &RotamerGetAtomByName(&rot, "CB")->xyz,
              &RotamerGetAtomByName(&rot, "OG1")->xyz,
              ic.icParam,
              &RotamerGetAtomByName(&rot, "HG1")->xyz);
            XYZArraySet(&rot.xyzs, atNdx, &RotamerGetAtomByName(&rot, "HG1")->xyz);
            for (int kk = 0; kk < RotamerGetAtomNum(&rot); kk++)
            {
              Atom* pAtomK = RotamerGetAtom(&rot, kk);
              if (strcmp(AtomGetName(pAtomK), "HG1") != 0)
              {
                continue;
              }
              for (int ss = 0; ss < RotamerGetAtomNum(&rot); ss++)
              {
                Atom* pAtomS = RotamerGetAtom(&rot, ss);
                if (strcmp(AtomGetName(pAtomK), AtomGetName(pAtomS)) == 0)
                {
                  continue;
                }
                int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtomK), AtomGetName(pAtomS), RotamerGetBonds(&rot));
                if (bondType == 14 || bondType == 15)
                {
                  double distance = XYZDist(&pAtomK->xyz, &pAtomS->xyz);
                  //if(distance < (pAtomK->vdw_radius+pAtomS->vdw_radius)*0.75){
                  if (distance < 2.25)
                  {
                    isAccepted = FALSE;
                    break;
                  }
                }
              }
              if (isAccepted == FALSE)
              {
                break;
              }
            }
            if (isAccepted == TRUE)
            {
              RotamersAdd(pSet, &rot);
              addNum++;
            }
          }

        }
      }
      RotamerDestroy(&rot);
    }

    if (RotamersFindRep(pSet, "TYR") != NULL)
    {
      ResToposGet(resTopos, "TYR", &topo);
      ResTopologyFindCharmmIC(&topo, "HH", &ic);
      double icParaX = ic.icParam[2];
      int addNum = 0;
      Rotamer rot;
      RotamerCreate(&rot);
      for (int j = 0; j < rotNum; j++)
      {
        Rotamer* pRot = RotamersGet(pSet, j);
        if (strcmp(RotamerGetType(pRot), "TYR") == 0)
        {
          RotamerRestore(pRot, pSet);
          int atNdx = -1;
          RotamerFindAtom(pRot, "HH", &atNdx);
          RotamerCopy(&rot, pRot);
          RotamerExtract(pRot);
          for (int k = 0; k < EXPANDED_ROT_TYR; k++)
          {
            ic.icParam[2] = icParaX + 2.0 * PI * (k + 1) / (EXPANDED_ROT_TYR + 1);
            if (ic.icParam[2] > PI) ic.icParam[2] -= 2 * PI;
            GetFourthAtom(&RotamerGetAtomByName(&rot, "CE1")->xyz,
              &RotamerGetAtomByName(&rot, "CZ")->xyz,
              &RotamerGetAtomByName(&rot, "OH")->xyz,
              ic.icParam,
              &RotamerGetAtomByName(&rot, "HH")->xyz);
            XYZArraySet(&rot.xyzs, atNdx, &RotamerGetAtomByName(&rot, "HH")->xyz);
            RotamersAdd(pSet, &rot);
            addNum++;
          }
        }

      }
      RotamerDestroy(&rot);
    }
    ResTopologyDestroy(&topo);
    CharmmICDestroy(&ic);
  }
  return Success;
}

int StructureGenerateWildtypeRotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain* pChn = StructureGetChain(pStruct, i);
    if (ChainGetType(pChn) != Type_Chain_Protein)
    {
      continue;
    }
    if (strstr(DES_CHAINS, ChainGetName(pChn)) == NULL)
    {
      continue;
    }
    for (int j = 0; j < ChainGetResidueNum(pChn); j++)
    {
      Residue* pRes = ChainGetResidue(pChn, j);
      if (pRes->desType != Type_DesType_Fixed)
      {
        continue;
      }
      ResidueSetDesignType(pRes, Type_DesType_Repackable);
      ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resiTopos);
      if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE)
      {
        ProteinSiteBuildNativeRotamer(pStruct, i, j, resiTopos);
      }
      if (FLAG_ROTATE_HYDROXYL == TRUE)
      {
        ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resiTopos);
      }
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain* pChn = StructureGetChain(pStruct, i);
    if (ChainGetType(pChn) != Type_Chain_Protein)
    {
      continue;
    }
    if (strstr(DES_CHAINS, ChainGetName(pChn)) != NULL)
    {
      continue;
    }
    for (int j = 0; j < ChainGetResidueNum(pChn); j++)
    {
      Residue* pRes = ChainGetResidue(pChn, j);
      if (pRes->desType != Type_DesType_Fixed)
      {
        continue;
      }
      Atom* pAtomCA1 = ResidueGetAtomByName(pRes, "CA");
      BOOL interResi = FALSE;
      for (int k = 0; k < StructureGetChainNum(pStruct); k++)
      {
        if (k == i)
        {
          continue;
        }
        Chain* pChnK = StructureGetChain(pStruct, k);
        if (ChainGetType(pChnK) != Type_Chain_Protein)
        {
          continue;
        }
        if (strstr(DES_CHAINS, ChainGetName(pChnK)) == NULL)
        {
          continue;
        }
        for (int s = 0;s < ChainGetResidueNum(pChnK);s++)
        {
          Residue* pResKS = ChainGetResidue(pChnK, s);
          Atom* pAtomCA2 = ResidueGetAtomByName(pResKS, "CA");
          if (XYZDist(&pAtomCA2->xyz, &pAtomCA1->xyz) > 15.0)
          {
            continue;
          }
          if (AtomArrayCalcMinDistance(&pRes->atoms, &pResKS->atoms) < CUT_PPI_DIST_SHELL1)
          {
            interResi = TRUE;
            break;
          }
        }
        if (interResi == TRUE)
        {
          break;
        }
      }
      if (interResi == TRUE)
      {
        ResidueSetDesignType(pRes, Type_DesType_Repackable);
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStruct, i, j, rotlib, atomParams, resiTopos);
        if (FLAG_USE_INPUT_SC == TRUE && pRes->isSCIntact == TRUE)
        {
          ProteinSiteBuildNativeRotamer(pStruct, i, j, resiTopos);
        }
        if (FLAG_ROTATE_HYDROXYL == TRUE)
        {
          ProteinSiteExpandHydroxylRotamers(pStruct, i, j, resiTopos);
        }
      }
    }
  }

  return Success;
}


int SelfEnergyGenerate(Structure* pStruct, char* file)
{
  int result = Success;
  char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
  FILE* pOut = fopen(file, "w");
  if (pOut == NULL)
  {
    sprintf(errMsg, "in file %s line %d, cannot create file %s to write", __FILE__, __LINE__, file);
    result = IOError;
    TraceError(errMsg, result);
    return result;
  }

  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++)
  {
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    Chain* pChn = StructureGetChain(pStruct, pSite->chnNdx);
    for (int j = 0; j < RotamersGetNum(pSet); j++)
    {
      Rotamer* pRot = RotamersGet(pSet, j);
      RotamerRestore(pRot, pSet);
      double energyTerms[MAX_ENERGY_TERM] = { 0 };
      if (ChainGetType(pChn) == Type_Chain_Protein)
      {
        for (int a = 0; a < StructureGetChainNum(pStruct); a++)
        {
          Chain* pChnA = StructureGetChain(pStruct, a);
          if (a == pSite->chnNdx)
          { // the same chain
            for (int b = 0; b < ChainGetResidueNum(pChnA); b++)
            {
              Residue* pResAB = ChainGetResidue(pChnA, b);
              if (b == pSite->resNdx)
              { // the same residue
                AminoAcidReferenceEnergy(RotamerGetType(pRot), energyTerms);
                EnergyIntraRotamer(pRot, energyTerms);
              }
              else
              { // different residues
                if (ChainGetType(pChnA) != Type_Chain_SmallMol)
                {
                  if (pResAB->desType == Type_DesType_Fixed)
                  {
                    EnergyRotamerAndFixedResidueSameChain(pRot, pResAB, energyTerms);
                  }
                  else if ( pResAB->desType == Type_DesType_Repackable ||
                            pResAB->desType == Type_DesType_Mutable ||
                            pResAB->desType == Type_DesType_Catalytic )
                  {
                    EnergyRotamerAndDesignResidueSameChain(pRot, pResAB, energyTerms);
                  }
                }
              }
            }
          }
          else
          { // different chains
            if (ChainGetType(pChnA) != Type_Chain_SmallMol)
            {
              for (int b = 0; b < ChainGetResidueNum(pChnA); b++)
              {
                Residue* pResAB = ChainGetResidue(pChnA, b);
                if (pResAB->desType == Type_DesType_Fixed)
                {
                  EnergyRotamerAndFixedResidueDiffChain(pRot, pResAB, energyTerms);
                }
                else if ( pResAB->desType == Type_DesType_Repackable ||
                          pResAB->desType == Type_DesType_Mutable ||
                          pResAB->desType == Type_DesType_Catalytic )
                {
                  EnergyRotamerAndDesignResidueDiffChain(pRot, pResAB, energyTerms);
                }
              }
            }
            else
            { // fixed small molecule
              for (int b = 0; b < ChainGetResidueNum(pChnA); b++)
              {
                Residue* pResAB = ChainGetResidue(pChnA, b);
                if (pResAB->desType == Type_DesType_Fixed)
                {
                  EnergyRotamerAndFixedLigResidue(pRot, pResAB, energyTerms);
                }
                else if ( pResAB->desType == Type_DesType_Repackable ||
                          pResAB->desType == Type_DesType_Mutable ||
                          pResAB->desType == Type_DesType_Catalytic )
                {
                  EnergyRotamerAndDesignResidueDiffChain(pRot, pResAB, energyTerms);
                }
              }
            }
          }
        }
      }

      else if (ChainGetType(pChn) == Type_Chain_SmallMol)
      {
        for (int a = 0; a < StructureGetChainNum(pStruct); a++)
        {
          Chain* pChnA = StructureGetChain(pStruct, a);
          if (a == pSite->chnNdx)
          { // the same chain
            pRot->selfEnergy += pRot->vdwInternal;
          }
          else
          { // different chains
            for (int b = 0; b < ChainGetResidueNum(pChnA); b++)
            {
              Residue* pResAB = ChainGetResidue(pChnA, b);
              if (pResAB->desType == Type_DesType_Fixed)
              {
                EnergyLigRotamerAndFixedResidue(pRot, pResAB, energyTerms);
              }
              else if ( pResAB->desType == Type_DesType_Repackable ||
                        pResAB->desType == Type_DesType_Mutable ||
                        pResAB->desType == Type_DesType_Catalytic )
              {
                EnergyLigRotamerAndDesignResidue(pRot, pResAB, energyTerms);
              }
            }
          }
        }
      }
      EnergyTermWeighting(energyTerms);
      pRot->selfEnergy += energyTerms[0];
      fprintf(pOut, "%d %d %f\n", i, j, pRot->selfEnergy);
      RotamerExtract(pRot);
    }
  }
  fclose(pOut);

  return result;
}


int SelfEnergyGenerate2(Structure* pStruct, AAppTable* pAAppTable, RamaTable* pRamaTable, char* file)
{
  FILE* pOut = fopen(file, "w");
  if (pOut == NULL)
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, IOError);
    return IOError;
  }

  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++)
  {
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    Chain* pChn = StructureGetChain(pStruct, pSite->chnNdx);
    for (int j = 0; j < RotamersGetNum(pSet); j++)
    {
      Rotamer* pRot = RotamersGet(pSet, j);
      RotamerRestore(pRot, pSet);
      double energyTerms[MAX_ENERGY_TERM] = { 0 };
      double energyTermsBind[MAX_ENERGY_TERM] = { 0 };
      if (ChainGetType(pChn) == Type_Chain_Protein)
      {
        for (int a = 0; a < StructureGetChainNum(pStruct); a++)
        {
          Chain* pChnA = StructureGetChain(pStruct, a);
          if (a == pSite->chnNdx)
          { // the same chain
            for (int b = 0; b < ChainGetResidueNum(pChnA); b++)
            {
              Residue* pResAB = ChainGetResidue(pChnA, b);
              if (b == pSite->resNdx)
              { // the same residue
                AminoAcidReferenceEnergy(RotamerGetType(pRot), energyTerms);
                EnergyIntraRotamer(pRot, energyTerms);
                RotamerPropensityAndRamachandranEnergy(pRot, pResAB, pAAppTable, pRamaTable, energyTerms);
                RotamerDunbrackEnergy(pRot, energyTerms);
              }
              else
              {
                if ( ChainGetType(pChnA) == Type_Chain_Protein ||
                     ChainGetType(pChnA) == Type_Chain_DNA ||
                     ChainGetType(pChnA) == Type_Chain_RNA ||
                     ChainGetType(pChnA) == Type_Chain_Water )
                {
                  if (pResAB->desType == Type_DesType_Fixed)
                  {
                    EnergyRotamerAndFixedResidueSameChain(pRot, pResAB, energyTerms);
                  }
                  else if ( pResAB->desType == Type_DesType_Repackable ||
                            pResAB->desType == Type_DesType_Mutable ||
                            pResAB->desType == Type_DesType_Catalytic ||
                            pResAB->desType == Type_DesType_NatRot )
                  {
                    EnergyRotamerAndDesignResidueSameChain(pRot, pResAB, energyTerms);
                  }
                }
              }
            }
          }
          else
          { // different chains
            if ( ChainGetType(pChnA) == Type_Chain_Protein ||
                 ChainGetType(pChnA) == Type_Chain_DNA ||
                 ChainGetType(pChnA) == Type_Chain_RNA ||
                 ChainGetType(pChnA) == Type_Chain_Water )
            {
              for (int b = 0; b < ChainGetResidueNum(pChnA); b++)
              {
                Residue* pResAB = ChainGetResidue(pChnA, b);
                if (pResAB->desType == Type_DesType_Fixed)
                {
                  EnergyRotamerAndFixedResidueDiffChain(pRot, pResAB, energyTerms);
                  if ( (strstr(DES_CHAINS, ChainGetName(pChnA)) == NULL && strstr(DES_CHAINS, ChainGetName(pChn)) != NULL) ||
                       (strstr(DES_CHAINS, ChainGetName(pChnA)) != NULL && strstr(DES_CHAINS, ChainGetName(pChn)) == NULL))
                  {
                    if (FLAG_DESIGN_PPI == TRUE)
                    {
                      EnergyRotamerAndFixedResidueDiffChain(pRot, pResAB, energyTermsBind);
                    }
                  }
                }
                else if ( pResAB->desType == Type_DesType_Repackable ||
                          pResAB->desType == Type_DesType_Mutable ||
                          pResAB->desType == Type_DesType_Catalytic ||
                          pResAB->desType == Type_DesType_NatRot)
                {
                  EnergyRotamerAndDesignResidueDiffChain(pRot, pResAB, energyTerms);
                  if ( (strstr(DES_CHAINS, ChainGetName(pChnA)) == NULL && strstr(DES_CHAINS, ChainGetName(pChn)) != NULL) ||
                       (strstr(DES_CHAINS, ChainGetName(pChnA)) != NULL && strstr(DES_CHAINS, ChainGetName(pChn)) == NULL))
                  {
                    if (FLAG_DESIGN_PPI == TRUE)
                    {
                      EnergyRotamerAndDesignResidueDiffChain(pRot, pResAB, energyTermsBind);
                    }
                  }
                }
              }
            }
            else if (ChainGetType(pChnA) == Type_Chain_SmallMol)
            { // chainA is small molecule
              for (int b = 0; b < ChainGetResidueNum(pChnA); b++)
              {
                Residue* pResAB = ChainGetResidue(pChnA, b);
                if (pResAB->desType == Type_DesType_Fixed)
                {
                  EnergyRotamerAndFixedLigResidue(pRot, pResAB, energyTerms);
                  if (FLAG_DESIGN_PLI == TRUE || FLAG_DESIGN_ENZ == TRUE)
                  {
                    EnergyRotamerAndFixedLigResidue(pRot, pResAB, energyTermsBind);
                  }
                }
              }
            }
          }
        }
      }
      else if (ChainGetType(pChn) == Type_Chain_SmallMol)
      {
        for (int a = 0; a < StructureGetChainNum(pStruct); a++)
        {
          Chain* pChnA = StructureGetChain(pStruct, a);
          if (a == pSite->chnNdx)
          { // the same chain
            pRot->selfEnergy += pRot->vdwInternal;
            // do not add vdwBackbone because it will be re-calculated
          }
          else
          { // different chains
            for (int b = 0; b < ChainGetResidueNum(pChnA); b++)
            {
              Residue* pResAB = ChainGetResidue(pChnA, b);
              if (pResAB->desType == Type_DesType_Fixed)
              {
                EnergyLigRotamerAndFixedResidue(pRot, pResAB, energyTerms);
                EnergyLigRotamerAndFixedResidue(pRot, pResAB, energyTermsBind);
              }
              else if ( pResAB->desType == Type_DesType_Repackable ||
                        pResAB->desType == Type_DesType_Mutable ||
                        pResAB->desType == Type_DesType_Catalytic ||
                        pResAB->desType == Type_DesType_NatRot )
              {
                EnergyLigRotamerAndDesignResidue(pRot, pResAB, energyTerms);
                EnergyLigRotamerAndDesignResidue(pRot, pResAB, energyTermsBind);
              }
            }
          }
        }
      }
      EnergyTermWeighting(energyTerms);
      pRot->selfEnergy += energyTerms[0];
      EnergyTermWeighting(energyTermsBind);
      pRot->selfEnergyBin += energyTermsBind[0];
      fprintf(pOut, "%d %d %f %f\n", i, j, pRot->selfEnergy, pRot->selfEnergyBin);
      RotamerExtract(pRot);
    }
  }
  fclose(pOut);

  return Success;
}


int SelfEnergyReadAndCheck(Structure* pStruct, RotamerList* pList, char* file)
{
  FILE* pIn = fopen(file, "r");
  if (pIn == NULL)
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, IOError);
    return IOError;
  }

  char buffer[MAX_LEN_ONE_LINE_CONTENT+1];
  while (fgets(buffer, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    int siteNdx = -1;
    int rotNdx = -1;
    double selfTot = 1000.0;
    double selfBin = 1000.0;
    sscanf(buffer, "%d %d %lf %lf\n", &siteNdx, &rotNdx, &selfTot, &selfBin);
    DesignSite* pSite = StructureGetDesignSite(pStruct, siteNdx);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    Rotamer* pRot = RotamersGet(pSet, rotNdx);
    pRot->selfEnergy = selfTot;
    pRot->selfEnergyBin = selfBin;
  }

  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++)
  {
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    for (int j = 0; j < RotamersGetNum(pSet); j++)
    {
      Rotamer* pRot = RotamersGet(pSet, j);
      if (pList->remainFlag[i][j] == FALSE)
      {
        continue;
      }
      for (int k = 0; k < RotamersGetNum(pSet); k++)
      {
        Rotamer* pRotK = RotamersGet(pSet, k);
        if (pList->remainFlag[i][k] == FALSE)
        {
          continue;
        }
        if (ResidueGetDesignType(pSite->pRes) == Type_DesType_SmallMol)
        {
          if (pRot->selfEnergy - pRotK->selfEnergy > SELF_ENERGY_THRESHOLD_SMALLMOL)
          {
            pList->remainFlag[i][j] = FALSE;
            break;
          }
        }
        else
        {
          if ( (pRot->selfEnergy - pRotK->selfEnergy > SELF_ENERGY_THRESHOLD_DIFF_AA_RES_TYPE && RotamerAndRotamerInSameType(pRot, pRotK) == FALSE) ||
               (pRot->selfEnergy - pRotK->selfEnergy > SELF_ENERGY_THRESHOLD_SAME_AA_RES_TYPE && RotamerAndRotamerInSameType(pRot, pRotK) == TRUE) )
          {
            pList->remainFlag[i][j] = FALSE;
            break;
          }
        }
      }
    }
  }
  fclose(pIn);

  return Success;
}


int RotamerListCreateFromStructure(RotamerList* pThis, Structure* pStruct)
{
  pThis->desSiteNum = StructureGetDesignSiteNum(pStruct);
  pThis->rotNum = (int*)malloc(sizeof(int) * pThis->desSiteNum);
  pThis->remainFlag = (BOOL**)malloc(sizeof(BOOL*) * pThis->desSiteNum);
  for (int i = 0; i < pThis->desSiteNum; i++)
  {
    pThis->rotNum[i] = RotamersGetNum(DesignSiteGetRotamers(StructureGetDesignSite(pStruct, i)));
    pThis->remainFlag[i] = (BOOL*)malloc(sizeof(BOOL) * pThis->rotNum[i]);
    for (int j = 0; j < pThis->rotNum[i]; j++)
    {
      pThis->remainFlag[i][j] = TRUE;
    }
  }
  pThis->remainRotNum = (int*)malloc(sizeof(int) * pThis->desSiteNum);
  for (int i = 0; i < pThis->desSiteNum; i++)
  {
    pThis->remainRotNum[i] = 0;
    for (int j = 0; j < pThis->rotNum[i]; j++)
    {
      if (pThis->remainFlag[i][j])
      {
        pThis->remainRotNum[i]++;
      }
    }
  }
  return Success;
}


void RotamerListDestroy(RotamerList* pList)
{
  for (int i = 0; i < pList->desSiteNum; i++) free(pList->remainFlag[i]);
  free(pList->remainFlag);
  free(pList->rotNum);
  free(pList->remainRotNum);
}


int RotamerListCopy(RotamerList* pThis, RotamerList* pOther)
{
  RotamerListDestroy(pThis);
  pThis->desSiteNum = pOther->desSiteNum;
  pThis->rotNum = (int*)malloc(pThis->desSiteNum * sizeof(int));
  memcpy(pThis->rotNum, pOther->rotNum, sizeof(int) * pThis->desSiteNum);
  pThis->remainFlag = (BOOL**)malloc(pThis->desSiteNum * sizeof(BOOL*));
  for (int i = 0; i < pThis->desSiteNum; i++)
  {
    pThis->remainFlag[i] = (BOOL*)malloc(pThis->rotNum[i] * sizeof(BOOL));
    memcpy(pThis->remainFlag[i], pOther->remainFlag[i], sizeof(BOOL) * pThis->rotNum[i]);
  }
  pThis->remainRotNum = (int*)malloc(pThis->desSiteNum * sizeof(int));
  memcpy(pThis->remainRotNum, pOther->remainRotNum, sizeof(int) * pThis->desSiteNum);
  return Success;
}


int RotamerListRead(RotamerList* pList, char* file)
{
  int result = Success;
  char errMsg[MAX_LEN_ERR_MSG+1];
  FILE* pIn = fopen(file, "r");
  if (pIn == NULL)
  {
    sprintf(errMsg, "in file %s line %d, cannot read file %s", __FILE__, __LINE__, file);
    result = IOError;
    TraceError(errMsg, result);
    return result;
  }

  char line[MAX_LEN_ONE_LINE_CONTENT+1];
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    int i, j;
    char flag[MAX_LEN_ONE_LINE_CONTENT+1];
    sscanf(line, "%d %d %s", &i, &j, flag);
    if (i >= pList->desSiteNum)
    {
      sprintf(errMsg, "in file %s line %d, structure does not have site #%d at this line:\n%s", __FILE__, __LINE__, i, line);
      result = ValueError;
      TraceError(errMsg, result);
      return result;
    }
    if (j >= pList->rotNum[i])
    {
      sprintf(errMsg, "in file %s line %d, structure does not have rotamer #%d on site #%d at this line:\n%s", __FILE__, __LINE__, j, i, line);
      result = ValueError;
      TraceError(errMsg, result);
      return result;
    }
    if (strcmp(flag, "TRUE") == 0)
    {
      pList->remainFlag[i][j] = TRUE;
    }
    else if (strcmp(flag, "FALSE") == 0)
    {
      pList->remainFlag[i][j] = FALSE;
    }
    else
    {
      sprintf(errMsg, "in file %s line %d, incorrect file format at this line:\n%s", __FILE__, __LINE__, line);
      result = FormatError;
      TraceError(errMsg, result);
      return result;
    }
  }

  for (int i = 0; i < pList->desSiteNum; i++)
  {
    pList->remainRotNum[i] = 0;
    for (int j = 0; j < pList->rotNum[i]; j++)
    {
      if (pList->remainFlag[i][j])
      {
        pList->remainRotNum[i]++;
      }
    }
  }

  return Success;
}


int RotamerListWrite(RotamerList* pList, char* file)
{
  int result;
  char errMsg[MAX_LEN_ERR_MSG+1];
  FILE* pOut = fopen(file, "w");
  if (pOut == NULL)
  {
    result = IOError;
    sprintf(errMsg, "in file %s line %d, cannot read file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0; i < pList->desSiteNum; i++)
  {
    for (int j = 0; j < pList->rotNum[i]; j++)
    {
      fprintf(pOut, "%d\t%d\t%s\n", i, j,
        pList->remainFlag[i][j] ? "TRUE" : "FALSE");
    }
  }
  fclose(pOut);
  return Success;
}



int RotamerListShow(RotamerList* pList)
{
  int totRotNum = 0;
  for (int i = 0; i < pList->desSiteNum; i++)
  {
    int remainRotNum = 0;
    printf("site %d, %d rotamers:\n", i, pList->rotNum[i]);
    for (int j = 0; j < pList->rotNum[i]; j++)
    {
      if (pList->remainFlag[i][j])
      {
        printf("rotamer %d\n", j);
        remainRotNum++;
      }
    }
    totRotNum += remainRotNum;
    printf("%d rotamers left\n", remainRotNum);
  }
  printf("number of total rotamers:\t%d\n", totRotNum);
  return Success;
}


int StructureShowDesignSites2(Structure* pStruct, RotamerList* pList)
{
  int result = Success;
  char errMsg[MAX_LEN_ERR_MSG+1];
  double space = 0.0;
  int totRotNum = 0;
  if (StructureGetDesignSiteNum(pStruct) != pList->desSiteNum) {
    sprintf(errMsg, "in file %s line %d, number of design site in structure (%d) is not equal to that in rotamer list (%d)",
      __FILE__, __LINE__, StructureGetDesignSiteNum(pStruct), pList->desSiteNum);
    result = ValueError;
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0; i < StructureGetDesignSiteNum(pStruct); i++) {
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    int remainRotNum = 0;
    if (RotamersGetNum(pSet) != pList->rotNum[i]) {
      sprintf(errMsg, "in file %s line %d, number of rotamers (%d) in site %d of structure is not equal to that in rotamer list (%d)",
        __FILE__, __LINE__, RotamersGetNum(pSet), i, pList->rotNum[i]);
      result = ValueError;
      TraceError(errMsg, result);
      return result;
    }
    for (int j = 0; j < RotamersGetNum(pSet); j++) {
      Rotamer* pRot = RotamersGet(pSet, j);
      if (pList->remainFlag[i][j] == TRUE) remainRotNum++;
    }
    printf("site %3d: %s %4d %3s, %6d rotamers, ",
      i, ResidueGetChainName(pSite->pRes), ResidueGetPosInChain(pSite->pRes), ResidueGetName(pSite->pRes), remainRotNum);
    totRotNum += remainRotNum;
    if (remainRotNum > 0) space += log((double)remainRotNum) / log(10.0);
    switch (pSite->pRes->desType) {
      case Type_DesType_Catalytic:
        printf("Catalytic\n"); break;
      case Type_DesType_Fixed:
        printf("Fixed\n"); break;
      case Type_DesType_Mutable:
        printf("Mutable\n"); break;
      case Type_DesType_SmallMol:
        printf("Smallmol\n"); break;
      case Type_DesType_Repackable:
        printf("Repackable\n"); break;
      case Type_DesType_NatRot:
        printf("Fixed_Natrot\n"); break;
      default:
        break;
    }
  }
  fprintf(stdout, "number of total rotamers: %d, conformation space: 1e^%.1f\n", totRotNum, space);

  return Success;
}

int RotamerOriginalIndexGet(RotamerList* pList, int siteNdx, int rotNdx, int* trueRotNdx)
{
  int j, flag = 0;
  for (j = 0; j < pList->rotNum[siteNdx]; j++)
  {
    if (pList->remainFlag[siteNdx][j] == FALSE)
    {
      continue;
    }
    if (flag == rotNdx)
    {
      break;
    }
    flag++;
  }
  if (j == pList->rotNum[siteNdx])
  {
    return DataNotExistError;
  }
  else
  {
    *trueRotNdx = j;
  }

  return Success;
}


int RotamerReducedIndexGet(RotamerList* pList, int siteNdx, int trueRotNdx, int* reducedRotNdx)
{
  int j, flag = 0;
  for (j = 0; j < pList->rotNum[siteNdx]; j++)
  {
    if (pList->remainFlag[siteNdx][j] == FALSE)
    {
      continue;
    }
    if (j == trueRotNdx)
    {
      break;
    }
    flag++;
  }
  if (j == pList->rotNum[siteNdx])
  {
    return DataNotExistError;
  }
  else
  {
    *reducedRotNdx = flag;
  }

  return Success;
}


int RotamerListCreateFromFile(RotamerList* pList, char* file)
{
  int result;
  char errMsg[MAX_LEN_ERR_MSG+1];
  FILE* pIn = fopen(file, "r");
  if (pIn == NULL)
  {
    sprintf(errMsg, "in file %s line %d, cannot open file %s", __FILE__, __LINE__, file);
    result = IOError;
    TraceError(errMsg, result);
    return result;
  }

  int desSiteNum = 0;
  char line[MAX_LEN_ONE_LINE_CONTENT+1];
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    int i, j;
    char flag[MAX_LEN_ONE_LINE_CONTENT+1];
    if (sscanf(line, "%d %d %s", &i, &j, flag) != 3)
    {
      sprintf(errMsg, "in file %s line %d, wrong format in file %s line\n%s", __FILE__, __LINE__, file, line);
      result = IOError;
      TraceError(errMsg, result);
      return result;
    }
    if (i >= desSiteNum) desSiteNum++;
  }

  IntArray remainRotNumArray;
  IntArrayCreate(&remainRotNumArray, desSiteNum);
  for (int i = 0; i < IntArrayGetLength(&remainRotNumArray); i++)
  {
    IntArraySet(&remainRotNumArray, i, 0);
  }

  fseek(pIn, 0, SEEK_SET);
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    int i, j;
    char flag[MAX_LEN_ONE_LINE_CONTENT+1];
    sscanf(line, "%d %d %s", &i, &j, flag);
    if (j >= IntArrayGet(&remainRotNumArray, i))
    {
      IntArraySet(&remainRotNumArray, i, j + 1);
    }
  }
  fseek(pIn, 0, SEEK_SET);

  pList->desSiteNum = desSiteNum;
  pList->rotNum = (int*)malloc(sizeof(int) * pList->desSiteNum);
  pList->remainFlag = (BOOL**)malloc(sizeof(BOOL*) * pList->desSiteNum);
  for (int i = 0; i < pList->desSiteNum; i++)
  {
    pList->rotNum[i] = IntArrayGet(&remainRotNumArray, i);
    pList->remainFlag[i] = (BOOL*)malloc(sizeof(BOOL) * pList->rotNum[i]);
    for (int j = 0; j < pList->rotNum[i]; j++)
    {
      pList->remainFlag[i][j] = TRUE;
    }
  }
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pIn))
  {
    int i, j;
    char flag[MAX_LEN_ONE_LINE_CONTENT+1];
    sscanf(line, "%d %d %s", &i, &j, flag);
    if (strcmp(flag, "TRUE") == 0)
    {
      pList->remainFlag[i][j] = TRUE;
    }
    else if (strcmp(flag, "FALSE") == 0)
    {
      pList->remainFlag[i][j] = FALSE;
    }
  }
  fclose(pIn);
  IntArrayDestroy(&remainRotNumArray);

  return Success;
}

