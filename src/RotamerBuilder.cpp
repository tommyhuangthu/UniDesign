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
#pragma warning(disable:6308)
#pragma warning(disable:28182)
#pragma warning(disable:26812)
#pragma warning(disable:6387)

#include "RotamerBuilder.h"
#include "EnergyFunction.h"
#include <string.h>

extern double CUT_PPI_DIST_SHELL1;
extern double CUT_PPI_DIST_SHELL2;
extern double CUT_PLI_DIST_SHELL1;
extern double CUT_PLI_DIST_SHELL2;


extern BOOL FLAG_MONOMER;
extern BOOL FLAG_PPI;
extern BOOL FLAG_PROT_LIG;
extern BOOL FLAG_ENZYME;

extern BOOL FLAG_USE_INPUT_SC;
extern BOOL FLAG_ROTATE_HYDROXYL;
extern BOOL FLAG_WILDTYPE_ONLY;
extern char DES_CHAINS[MAX_LEN_ONE_LINE_CONTENT + 1];

extern BOOL FLAG_EXCL_CYS_ROTS;

#define DEAL_WITH_PROTEIN_ROTAMERS_BBIND
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// functions used to deal with sidechain rots, sidechain repacking and protein design
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ProteinSiteBuildAllRotamers(Structure* pStruct, int chnNdx, int resNdx, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos)
{
  int result = Success;
  DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pSite != NULL)
  {
    DesignSiteRemoveRotamers(pSite);
  }
  else
  {
    (pStruct->desSiteCount)++;
    pStruct->designSites = (DesignSite*)realloc(pStruct->designSites, sizeof(DesignSite) * pStruct->desSiteCount);
    DesignSiteCreate(&pStruct->designSites[pStruct->desSiteCount - 1]);
    pSite = StructureGetDesignSite(pStruct, pStruct->desSiteCount - 1);
    Chain* pDestChain = StructureGetChain(pStruct, chnNdx);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resNdx);
    pSite->pRes = pDestResidue;
    pSite->chnNdx = chnNdx;
    pSite->resNdx = resNdx;
  }

  // set design types - 20 AA types;
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
  result = RotamerSetOfProteinGenerate(DesignSiteGetRotamers(pSite), pSite->pRes, &designTypes, &patchTypes, rotlib, atomParams, resiTopos);
  ResidueSetDesignType(pSite->pRes, Type_DesType_Mutable);
  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}


int ProteinSiteBuildMutatedRotamers(Structure* pStruct, int chnNdx, int resNdx, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, StringArray* pDesignTypes, StringArray* pPatchTypes)
{
  int result = Success;
  DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pSite != NULL)
  {
    DesignSiteRemoveRotamers(pSite);
  }
  else
  {
    (pStruct->desSiteCount)++;
    pStruct->designSites = (DesignSite*)realloc(pStruct->designSites, sizeof(DesignSite) * pStruct->desSiteCount);
    DesignSiteCreate(&pStruct->designSites[pStruct->desSiteCount - 1]);
    pSite = StructureGetDesignSite(pStruct, pStruct->desSiteCount - 1);
    Chain* pDestChain = StructureGetChain(pStruct, chnNdx);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resNdx);
    pSite->pRes = pDestResidue;
    pSite->chnNdx = chnNdx;
    pSite->resNdx = resNdx;
  }

  result = RotamerSetOfProteinGenerate(DesignSiteGetRotamers(pSite), pSite->pRes, pDesignTypes, pPatchTypes, rotlib, atomParams, resiTopos);
  ResidueSetDesignType(pSite->pRes, Type_DesType_Mutable);
  return result;
}


int ProteinSiteBuildWildtypeRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos)
{
  int result = Success;
  DesignSite* pSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if (pSite != NULL)
  {
    DesignSiteRemoveRotamers(pSite);
  }
  else
  {
    (pThis->desSiteCount)++;
    pThis->designSites = (DesignSite*)realloc(pThis->designSites, sizeof(DesignSite) * pThis->desSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->desSiteCount - 1]);
    pSite = StructureGetDesignSite(pThis, pThis->desSiteCount - 1);
    Chain* pDestChn = StructureGetChain(pThis, chainIndex);
    Residue* pDestRes = ChainGetResidue(pDestChn, resiIndex);
    pSite->pRes = pDestRes;
    pSite->chnNdx = chainIndex;
    pSite->resNdx = resiIndex;
  }

  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  StringArrayAppend(&designTypes, ResidueGetName(pSite->pRes));
  StringArrayAppend(&patchTypes, "");
  if (strcmp(ResidueGetName(pSite->pRes), "HSD") == 0)
  {
    StringArrayAppend(&designTypes, "HSE");
    StringArrayAppend(&patchTypes, "");
  }
  else if (strcmp(ResidueGetName(pSite->pRes), "HSE") == 0)
  {
    StringArrayAppend(&designTypes, "HSD");
    StringArrayAppend(&patchTypes, "");
  }

  result = RotamerSetOfProteinGenerate(DesignSiteGetRotamers(pSite), pSite->pRes, &designTypes, &patchTypes, rotlib, atomParams, resiTopos);
  ResidueSetDesignType(pSite->pRes, Type_DesType_Repackable);

  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}


int ProteinSiteBuildSpecifiedRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, Type_ResidueDesignType type)
{
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if (pCurrentDesignSite != NULL)
  {
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else
  {
    (pThis->desSiteCount)++;
    pThis->designSites = (DesignSite*)realloc(pThis->designSites, sizeof(DesignSite) * pThis->desSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->desSiteCount - 1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->desSiteCount - 1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pRes = pDestResidue;
    pCurrentDesignSite->chnNdx = chainIndex;
    pCurrentDesignSite->resNdx = resiIndex;
  }

  StringArray designTypes;
  StringArray patchTypes;
  StringArrayCreate(&designTypes);
  StringArrayCreate(&patchTypes);
  if (type == Type_DesType_Mutable)
  {
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
  }
  else if (type == Type_DesType_Repackable || type == Type_DesType_Catalytic)
  {
    StringArrayAppend(&designTypes, ResidueGetName(pCurrentDesignSite->pRes)); StringArrayAppend(&patchTypes, "");
  }
  result = RotamerSetOfProteinGenerate(DesignSiteGetRotamers(pCurrentDesignSite), pCurrentDesignSite->pRes, &designTypes, &patchTypes, rotlib, atomParams, resiTopos);
  StringArrayDestroy(&patchTypes);
  StringArrayDestroy(&designTypes);

  return result;
}


int StructureGenerateSpecifiedProteinRotamers(Structure* pStructure, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* resfile)
{
  FileReader fr;
  FileReaderCreate(&fr, resfile);
  char buffer[MAX_LEN_ONE_LINE_CONTENT + 1];
  while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
  {
    char keyword[MAX_LEN_ONE_LINE_CONTENT + 1];
    ExtractFirstStringFromSourceString(keyword, buffer);
    if (strcmp(keyword, "SITES_CATALYTIC_START") == 0)
    {
      StringArray strings;
      StringArrayCreate(&strings);
      while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
      {
        StringArraySplitString(&strings, buffer, ' ');
        if (strcmp(StringArrayGet(&strings, 0), "SITES_CATALYTIC_END") == 0)
        {
          break;
        }
        char* chnname = StringArrayGet(&strings, 0);
        int posInChain = atoi(StringArrayGet(&strings, 1));
        Chain* pChain = StructureFindChainByName(pStructure, chnname);
        int resiIndex = -1;
        ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
        if (resiIndex != -1)
        {
          ChainGetResidue(pChain, resiIndex)->desType = Type_DesType_Catalytic;
        }
      }
      StringArrayDestroy(&strings);
    }
    else if (strcmp(keyword, "SITES_DESIGN_START") == 0)
    {
      StringArray strings;
      StringArrayCreate(&strings);
      while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
      {
        StringArraySplitString(&strings, buffer, ' ');
        if (strcmp(StringArrayGet(&strings, 0), "SITES_DESIGN_END") == 0)
        {
          break;
        }
        char* chnname = StringArrayGet(&strings, 0);
        int posInChain = atoi(StringArrayGet(&strings, 1));
        Chain* pChain = StructureFindChainByName(pStructure, chnname);
        int resiIndex = -1;
        ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
        if (resiIndex != -1)
        {
          ChainGetResidue(pChain, resiIndex)->desType = Type_DesType_Mutable;
        }
      }
      StringArrayDestroy(&strings);
    }
    else if (strcmp(keyword, "SITES_REPACK_START") == 0)
    {
      StringArray strings;
      StringArrayCreate(&strings);
      while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
      {
        StringArraySplitString(&strings, buffer, ' ');
        if (strcmp(StringArrayGet(&strings, 0), "SITES_REPACK_END") == 0)
        {
          break;
        }
        char* chnname = StringArrayGet(&strings, 0);
        int posInChain = atoi(StringArrayGet(&strings, 1));
        Chain* pChain = StructureFindChainByName(pStructure, chnname);
        int resiIndex = -1;
        ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
        if (resiIndex != -1)
        {
          ChainGetResidue(pChain, resiIndex)->desType = Type_DesType_Repackable;
        }
      }
      StringArrayDestroy(&strings);
    }
  }
  FileReaderDestroy(&fr);

  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
    {
      Residue* pResi = ChainGetResidue(pChainI, j);
      if (pResi->desType != Type_DesType_Fixed)
      {
        //printf("generate rots for residue %d\n", j);
        ProteinSiteBuildSpecifiedRotamers(pStructure, i, j, rotlib, atomParams, resiTopos, pResi->desType);
        if (FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
        if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
      }
    }
  }

  return Success;
}


int StructureGenerateAllRotamers(Structure* pStructure, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos)
{
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL) continue;
    for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
    {
      Residue* pResi = ChainGetResidue(pChainI, j);
      if (pResi->desType != Type_DesType_Fixed) continue;
      ResidueSetDesignType(pResi, Type_DesType_Mutable);
      ProteinSiteBuildSpecifiedRotamers(pStructure, i, j, rotlib, atomParams, resiTopos, pResi->desType);
      if (FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
      if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
    }
  }

  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL) continue;
    for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
    {
      Residue* pResi = ChainGetResidue(pChainI, j);
      if (pResi->desType != Type_DesType_Fixed) continue;
      Atom* pAtomCA1 = ResidueGetAtomByName(pResi, "CA");
      BOOL interResi = FALSE;
      for (int k = 0; k < StructureGetChainCount(pStructure); k++)
      {
        if (k == i)continue;
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) != Type_Chain_Protein)continue;
        if (strstr(DES_CHAINS, ChainGetName(pChainK)) == NULL)continue;
        for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
        {
          Residue* pResidueKS = ChainGetResidue(pChainK, s);
          Atom* pAtomCA2 = ResidueGetAtomByName(pResidueKS, "CA");
          if (XYZDistance(&pAtomCA2->xyz, &pAtomCA1->xyz) > 15.0) continue;
          if (AtomArrayCalcMinDistance(&pResi->atoms, &pResidueKS->atoms) < CUT_PPI_DIST_SHELL1)
          {
            interResi = TRUE;
            break;
          }
        }
        if (interResi) break;
      }
      if (interResi)
      {
        ResidueSetDesignType(pResi, Type_DesType_Repackable);
        ProteinSiteBuildSpecifiedRotamers(pStructure, i, j, rotlib, atomParams, resiTopos, pResi->desType);
        if (FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
        if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
      }
    }
  }

  return Success;
}



#define DEAL_WITH_LIGAND_ROTAMERS

int StructureGetTruncatedBackbone(Structure* pThis, Residue* pSmallMol, double activeSiteRange, BOOL withHydrogen, AtomArray* pBackboneAtoms)
{
  for (int chainIndex = 0; chainIndex < StructureGetChainCount(pThis); chainIndex++)
  {
    Chain* pChain = StructureGetChain(pThis, chainIndex);
    for (int resiIndex = 0; resiIndex < ChainGetResidueCount(pChain); resiIndex++)
    {
      Residue* pResi = ChainGetResidue(pChain, resiIndex);
      if (strcmp(ResidueGetChainName(pResi), ResidueGetChainName(pSmallMol)) == 0 && ResidueGetPosInChain(pResi) == ResidueGetPosInChain(pSmallMol)) continue;
      double minDist = AtomArrayCalcMinDistance(ResidueGetAllAtoms(pResi), ResidueGetAllAtoms(pSmallMol));
      if (minDist > activeSiteRange) continue;
      for (int atomIndex = 0;atomIndex < ResidueGetAtomCount(pResi);atomIndex++)
      {
        Atom* pAtom = ResidueGetAtom(pResi, atomIndex);
        if (AtomIsHydrogen(pAtom) && !withHydrogen) continue;
        if (pAtom->isBBAtom == FALSE && strcmp(AtomGetName(pAtom), "CB") != 0) continue;
        AtomArrayAppend(pBackboneAtoms, pAtom);
      }
    }
  }
  return Success;
}


int StructureDeployCataConsSitePair(Structure* pThis, CataConsSitePair* pCataConsSitePair)
{
  char errMsg[MAX_LEN_ERR_MSG + 1];
  Residue* pSmallMol = NULL;
  RotamerSet* pRotSetOfFirstAndSecondSite[2] = { NULL,NULL };

  int result = StructureFindSmallMol(pThis, &pSmallMol);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, cannot find small molecule", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }

  RotamerSet tempSmallMolRotamerSet;
  RotamerSetCreate(&tempSmallMolRotamerSet);
  Rotamer tempSmallMolRotamer;
  RotamerCreate(&tempSmallMolRotamer);
  RotamerSetType(&tempSmallMolRotamer, ResidueGetName(pSmallMol));
  RotamerAddAtoms(&tempSmallMolRotamer, ResidueGetAllAtoms(pSmallMol));
  RotamerSetAdd(&tempSmallMolRotamerSet, &tempSmallMolRotamer);
  RotamerDestroy(&tempSmallMolRotamer);
  for (int i = 0;i < 2;i++)
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
      pRotSetOfFirstAndSecondSite[i] = &tempSmallMolRotamerSet;
    }
    else
    {
      int chainIndex = -1, resiIndex = -1;
      StructureFindChainIndex(pThis, chainName, &chainIndex);
      ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), posInChain, &resiIndex);
      DesignSite* pDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
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

  RotamerSetDestroy(&tempSmallMolRotamerSet);
  return result;
}


int StructurePlaceSmallMol(Structure* pThis, PlacingRule* pRule, CataConsSitePairArray* pConsArray, int siteCount, DesignSite** ppSites, RotamerSet* pSmallSet)
{
  int result;
  char errMsg[MAX_LEN_ERR_MSG + 1];
  int chnNdx = -1, resNdx = -1;
  StructureFindChainIndex(pThis, PlacingRuleGetChainName(pRule), &chnNdx);
  ChainFindResidueByPosInChain(StructureGetChain(pThis, chnNdx), PlacingRuleGetPosInChain(pRule), &resNdx);
  DesignSite* pStartSite = StructureFindDesignSite(pThis, chnNdx, resNdx);
  if (pStartSite == NULL)
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d, starting site %s%d%s is required for ligand placement", __FILE__, __LINE__,
      PlacingRuleGetChainName(pRule), PlacingRuleGetPosInChain(pRule), PlacingRuleGetResiName(pRule));
    TraceError(errMsg, result);
    return result;
  }

  result = PlacingRulePlaceSmallMol(pRule, pConsArray, pStartSite, siteCount, ppSites, pSmallSet);
  if (FAILED(result))
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d, failed to place ligand rotamers", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }

  return Success;
}


int StructureGenerateSmallMolRotamers(Structure* pThis, char* fileCataCons, char* fileLigPlacing)
{
  int result = Success;
  char errMsg[MAX_LEN_ERR_MSG + 1];

  printf("make ligand conformers for ligand pose sampling in enzyme design\n");

  printf("1. read and deploy catalytic constraints ...\n");
  CataConsSitePairArray cataCons;
  result = CataConsSitePairArrayCreate(&cataCons, fileCataCons);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, failed to read catalytic-constraint file %s", __FILE__, __LINE__, fileCataCons);
    result = ValueError;
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0;i < CataConsSitePairArrayGetCount(&cataCons);i++)
  {
    result = StructureDeployCataConsSitePair(pThis, CataConsSitePairArrayGet(&cataCons, i));
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
  result = PlacingRuleCreate(&placingRule, fileLigPlacing);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, failed to read ligand-placing file %s", __FILE__, __LINE__, fileLigPlacing);
    result = ValueError;
    return result;
  }

  Residue* pSmallMol = NULL;
  result = StructureFindSmallMol(pThis, &pSmallMol);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, cannot find small molecule", __FILE__, __LINE__);
    result = DataNotExistError;
    return result;
  }
  AtomArray truncBackbone;
  AtomArrayCreate(&truncBackbone);
  StructureGetTruncatedBackbone(pThis, pSmallMol, PlacingRuleGetTruncatedBackboneRange(&placingRule), TRUE, &truncBackbone);

  int nRelatedSites = 0;
  DesignSite** ppRelatedSites = NULL;
  for (int i = 0; i < CataConsSitePairArrayGetCount(&cataCons); i++)
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
        StructureFindChainIndex(pThis, pSitePair->chnName1, &chainIndex);
        ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), pSitePair->pos1, &resiIndex);
        pSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
      }
      else
      {
        if (ResidueGetPosInChain(pSmallMol) == pSitePair->pos2 && strcmp(ResidueGetChainName(pSmallMol), pSitePair->chnName2) == 0)
        {
          continue;
        }
        int chainIndex = -1, resiIndex = -1;
        StructureFindChainIndex(pThis, pSitePair->chnName2, &chainIndex);
        ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), pSitePair->pos2, &resiIndex);
        pSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
      }
      //smallmol design site is excluded
      BOOL designsiteExist = FALSE;
      for (int j = 0; j < nRelatedSites; j++)
      {
        DesignSite* pSite1 = ppRelatedSites[j];
        //skip the identical sites
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

  printf("3. start placing ligand rotamers ...\n");
  RotamerSet ligRots;
  RotamerSetCreate(&ligRots);
  result = StructurePlaceSmallMol(pThis, &placingRule, &cataCons, nRelatedSites, ppRelatedSites, &ligRots);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, failed to placing ligand rotamers", __FILE__, __LINE__);
    result = ValueError;
    return result;
  }
  if (RotamerSetGetCount(&ligRots) > 0)
  {
    BOOL smallmolRotSetGenerated = FALSE;
    for (int i = 0;i < nRelatedSites;i++)
    {
      if (ppRelatedSites[i]->pRes == pSmallMol)
      {
        for (int j = 0;j < RotamerSetGetCount(&ligRots);j++)
        {
          RotamerSetAdd(DesignSiteGetRotamers(ppRelatedSites[i]), RotamerSetGet(&ligRots, j));
        }
        smallmolRotSetGenerated = TRUE;
        break;
      }
    }
    if (smallmolRotSetGenerated == FALSE)
    {
      int chnNdx = -1, resNdx = -1;
      StructureFindChainIndex(pThis, ResidueGetChainName(pSmallMol), &chnNdx);
      ChainFindResidueByPosInChain(StructureGetChain(pThis, chnNdx), ResidueGetPosInChain(pSmallMol), &resNdx);
      ProteinSiteAddDesignSite(pThis, chnNdx, resNdx);
      RotamerSetCopy(DesignSiteGetRotamers(StructureGetDesignSite(pThis, pThis->desSiteCount - 1)), &ligRots);
    }
  }

  RotamerSetDestroy(&ligRots);
  PlacingRuleDestroy(&placingRule);
  CataConsSitePairArrayDestroy(&cataCons);
  AtomArrayDestroy(&truncBackbone);
  free(ppRelatedSites);
  ppRelatedSites = NULL;

  return Success;
}


int StructureWriteSmallMolRotamers(Structure* pThis, char* fileSmallMol)
{
  int result;
  char errMsg[MAX_LEN_ERR_MSG + 1];
  Residue* pSmallMol = NULL;
  StructureFindSmallMol(pThis, &pSmallMol);
  if (pSmallMol == NULL)
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d, cannot find the small molecule", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }
  int chainIndex = -1, resiIndex = -1;
  StructureFindChainIndex(pThis, ResidueGetChainName(pSmallMol), &chainIndex);
  ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), ResidueGetPosInChain(pSmallMol), &resiIndex);
  DesignSite* pSmallSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  RotamerSet* pSmallSet = DesignSiteGetRotamers(pSmallSite);
  if (pSmallSet == NULL)
  {
    result = DataNotExistError;
    sprintf(errMsg, "in file %s line %d, cannot find the small-molecule design site", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }

  FILE* outputFile = NULL;
  if (fileSmallMol == NULL || strcmp(fileSmallMol, "") == 0) outputFile = NULL;
  else
  {
    outputFile = fopen(fileSmallMol, "w");
    if (outputFile == NULL)
    {
      sprintf(errMsg, "In file %s line %d, cannot write to file %s", __FILE__, __LINE__, fileSmallMol);
      result = IOError;
      TraceError(errMsg, result);
      return result;
    }
  }

  for (int i = 0; i < RotamerSetGetCount(pSmallSet); i++)
  {
    Rotamer newRot;
    RotamerCreate(&newRot);
    Model(i + 1, outputFile);
    RotamerCopy(&newRot, RotamerSetGet(pSmallSet, i));
    RotamerRestore(&newRot, pSmallSet);
    RotamerShowInPDBFormat(&newRot, "ATOM", ResidueGetChainName(pSmallMol), 1, 1, outputFile);
    RotamerExtract(&newRot);
    // write energy into ligand conformers file
    fprintf(outputFile, "ENERGY INTERNAL: %f BACKBONE: %f\n", newRot.vdwInternal, newRot.vdwBackbone);
    EndModel(outputFile);
    RotamerDestroy(&newRot);
  }
  if (outputFile != NULL) fclose(outputFile);

  return Success;
}


int StructureReadSmallMolRotamers(Structure* pThis, ResiTopoSet* resiTopos, char* smallMolFileName)
{
  char errMsg[MAX_LEN_ERR_MSG + 1];
  Residue* pSmallMol = NULL;
  int result = StructureFindSmallMol(pThis, &pSmallMol);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, cannot find small molecule", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }

  int chainIndex = -1, resiIndex = -1;
  StructureFindChainIndex(pThis, ResidueGetChainName(pSmallMol), &chainIndex);
  ChainFindResidueByPosInChain(StructureGetChain(pThis, chainIndex), ResidueGetPosInChain(pSmallMol), &resiIndex);
  DesignSite* pSmallMolSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if (pSmallMolSite == NULL)
  {
    ProteinSiteAddDesignSite(pThis, chainIndex, resiIndex);
    pSmallMolSite = StructureGetDesignSite(pThis, pThis->desSiteCount - 1);
  }
  ResidueSetDesignType(pSmallMol, Type_DesType_SmallMol);

  FileReader fr;
  result = FileReaderCreate(&fr, smallMolFileName);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, cannot read file %s", __FILE__, __LINE__, smallMolFileName);
    TraceError(errMsg, result);
    FileReaderDestroy(&fr);
    return result;
  }
  while (!FileReaderEndOfFile(&fr))
  {
    char line[MAX_LEN_ONE_LINE_CONTENT + 1];
    char keyword[MAX_LEN_ONE_LINE_CONTENT + 1];
    FileReaderGetNextLine(&fr, line);
    ExtractFirstStringFromSourceString(keyword, line);
    if (strcmp(keyword, "MODEL") != 0) continue;

    Residue tmpRes;
    ResidueCreate(&tmpRes);
    ResidueCopy(&tmpRes, pSmallMol);
    for (int i = 0;i < ResidueGetAtomCount(&tmpRes);i++)
    {
      ResidueGetAtom(&tmpRes, i)->isXyzValid = FALSE;
    }
    if (FAILED(ResidueReadXYZFromPDB(&tmpRes, &fr)))
    {
      sprintf(errMsg, "in file %s line %d, error when reading small-molecule conformers file", __FILE__, __LINE__);
      result = FormatError;
      TraceError(errMsg, result);
      return result;
    }

    ResidueCheckAtomCoordinateValidity(&tmpRes);
    result = ResidueCalcAllAtomXYZ(&tmpRes, resiTopos, NULL, NULL);
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
      BondSetCopy(RotamerGetBonds(&newRot), ResidueGetBonds(&tmpRes));
      RotamerSetChainName(&newRot, ResidueGetChainName(&tmpRes));
      //read the energy from smallmol rotamer file
      newRot.vdwInternal = tmpRes.internalEnergy;
      newRot.vdwBackbone = tmpRes.backboneEnergy;
      //set smallmol internal energy as self energy and recalculate backbone energy later
      //newRot.selfEnergy = tmpRes.internalEnergy;
      RotamerSetPosInChain(&newRot, ResidueGetPosInChain(&tmpRes));
      RotamerSetAdd(DesignSiteGetRotamers(pSmallMolSite), &newRot);
      RotamerDestroy(&newRot);
    }
    ResidueDestroy(&tmpRes);
  }
  FileReaderDestroy(&fr);
  return Success;
}


// this function is used to delete or screen small molecules with correct direct orientation;
// oriFileName is small molecule library file by combining the tmpPDB files together;
// newFileName is the reserved library file after screening;
// screenFileName is the screening rule file;
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
// in each BIND_SITE_GROUP, all the distances between the binding atom and the CA atoms of the residues must be are 
// larger than those between the catalytic atom and the CA atoms of the residues;
// and at least one binding site group must satisfy the above relationship.
int StructureSmallmolOrientationScreen(Structure* pStructure, ResiTopoSet* pResiTopo, char* oriFileName, char* newFileName, char* screenRuleFileName)
{
  typedef struct _BindSite
  {
    char chainName[MAX_LEN_CHAIN_NAME + 1];
    int  posInChain;
    char resiName[MAX_LEN_RES_NAME + 1];
  } BindSite;

  typedef struct _BindSiteGroup
  {
    int       siteNum;
    BindSite* pSites;
  } BindSiteGroup;

  int            totalRotamerCounter, acceptedRotamerCounter;
  char           cataAtomName[MAX_LEN_ATOM_NAME + 1];
  char           bindAtomName[MAX_LEN_ATOM_NAME + 1];
  int            groupCounter;
  char           line[MAX_LEN_ONE_LINE_CONTENT + 1];
  BindSiteGroup* pGroups = NULL;
  FILE* pOutputFile = NULL;
  FILE* pInputFile = NULL;
  StringArray    buffer;
  XYZ            xyzCataAtom;
  XYZ            xyzBindAtom;
  BOOL           cataAtomExist = FALSE;
  BOOL           bindAtomExist = FALSE;
  BOOL           curRotamerAccepted = FALSE;

  FileReader     file;
  FileReaderCreate(&file, screenRuleFileName);
  groupCounter = 0;
  while (!FileReaderEndOfFile(&file))
  {
    BOOL doneInThisGroup = FALSE;
    int siteNum;
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
        groupCounter++;
        pGroups = (BindSiteGroup*)realloc(pGroups, sizeof(BindSiteGroup) * groupCounter);
        siteNum = 0;
        pGroups[groupCounter - 1].siteNum = 0;
        pGroups[groupCounter - 1].pSites = NULL;
        continue;
      }
      else if (strcmp(StringArrayGet(&wordsInLine, 0), "BIND_SITE") == 0)
      {
        siteNum++;
        pGroups[groupCounter - 1].pSites = (BindSite*)realloc(pGroups[groupCounter - 1].pSites, sizeof(BindSite) * siteNum);
        strcpy(pGroups[groupCounter - 1].pSites[siteNum - 1].chainName, StringArrayGet(&wordsInLine, 1));
        strcpy(pGroups[groupCounter - 1].pSites[siteNum - 1].resiName, StringArrayGet(&wordsInLine, 2));
        pGroups[groupCounter - 1].pSites[siteNum - 1].posInChain = atoi(StringArrayGet(&wordsInLine, 3));
        pGroups[groupCounter - 1].siteNum = siteNum;
      }
    }
  }

  char errMsg[MAX_LEN_ERR_MSG + 1];
  printf("reading small-molecule rotamers from file %s\n", oriFileName);
  pInputFile = fopen(oriFileName, "r");
  if (pInputFile == NULL)
  {
    sprintf(errMsg, "in file %s line %d, can not open file %s for reading", __FILE__, __LINE__, oriFileName);
    TraceError(errMsg, IOError);
    return IOError;
  }
  pOutputFile = fopen(newFileName, "w");
  if (pOutputFile == NULL)
  {
    sprintf(errMsg, "in file %s line %d, can not open file %s for reading", __FILE__, __LINE__, newFileName);
    TraceError(errMsg, IOError);
    return IOError;
  }

  totalRotamerCounter = 0;
  acceptedRotamerCounter = 0;
  StringArrayCreate(&buffer);
  while (fgets(line, MAX_LEN_ONE_LINE_CONTENT, pInputFile))
  {
    char keyword[10] = "";
    ExtractTargetStringFromSourceString(keyword, line, 0, 4);
    if (strcmp(keyword, "ENDM") == 0)
    {
      if (cataAtomExist == FALSE || bindAtomExist == FALSE)
      {
        sprintf(errMsg, "in file %s line %d, CataAtom or BindAtom does not exist in file %s", __FILE__, __LINE__, oriFileName);
        TraceError(errMsg, FormatError);
        return FormatError;
      }

      for (int i = 0; i < groupCounter; i++)
      {
        curRotamerAccepted = TRUE;
        for (int j = 0; j < pGroups[i].siteNum; j++)
        {
          Residue* pResidue = ChainGetResidue(StructureFindChainByName(pStructure, pGroups[i].pSites[j].chainName), pGroups[i].pSites[j].posInChain);
          XYZ* pXyzCA = &ResidueGetAtomByName(pResidue, "CA")->xyz;
          if (XYZDistance(&xyzBindAtom, pXyzCA) > XYZDistance(&xyzCataAtom, pXyzCA))
          {
            curRotamerAccepted = FALSE;
            break;
          }
        }
        if (curRotamerAccepted == TRUE)
        {
          fprintf(pOutputFile, "MODEL     %d\n", acceptedRotamerCounter);
          for (int j = 0;j < StringArrayGetCount(&buffer);j++)
          {
            fprintf(pOutputFile, "%s", StringArrayGet(&buffer, j));
          }
          fprintf(pOutputFile, "ENDMDL\n");
          acceptedRotamerCounter++;
          break;
        }
      }

      printf("%d / %d rotamers have been processed      \r", acceptedRotamerCounter, totalRotamerCounter);
      //fflush(stdout);
    }
    else if (strcmp(keyword, "MODE") == 0)
    {
      StringArrayDestroy(&buffer);
      StringArrayCreate(&buffer);
      totalRotamerCounter++;
    }
    else if (strcmp(keyword, "ATOM") == 0)
    {
      char atomName[MAX_LEN_ATOM_NAME + 1];
      ExtractTargetStringFromSourceString(atomName, line, 12, 4);
      if (strcmp(atomName, cataAtomName) == 0)
      {
        cataAtomExist = TRUE;
        ExtractTargetStringFromSourceString(keyword, line, 31, 7);
        xyzCataAtom.X = atof(keyword);
        ExtractTargetStringFromSourceString(keyword, line, 39, 7);
        xyzCataAtom.Y = atof(keyword);
        ExtractTargetStringFromSourceString(keyword, line, 47, 7);
        xyzCataAtom.Z = atof(keyword);
      }

      if (strcmp(atomName, bindAtomName) == 0)
      {
        bindAtomExist = TRUE;
        ExtractTargetStringFromSourceString(keyword, line, 31, 7);
        xyzBindAtom.X = atof(keyword);
        ExtractTargetStringFromSourceString(keyword, line, 39, 7);
        xyzBindAtom.Y = atof(keyword);
        ExtractTargetStringFromSourceString(keyword, line, 47, 7);
        xyzBindAtom.Z = atof(keyword);
      }
      StringArrayAppend(&buffer, line);
    }
    else if (strcmp(keyword, "ENER") == 0 && curRotamerAccepted)
    {
      fprintf(pOutputFile, "%s", line);
    }
  }
  printf("\n");
  fclose(pInputFile);
  fclose(pOutputFile);

  FileReaderDestroy(&file);
  StringArrayDestroy(&buffer);
  for (int j = 0; j < groupCounter; j++)
  {
    free(pGroups[j].pSites);
    pGroups[j].siteNum = 0;
  }
  free(pGroups);

  return Success;
}


#define DEAL_WITH_PROTEIN_ROTAMERS_BBDEP
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//build backbone dependent rots
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int ProteinSiteBuildAllRotamersByBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos)
{
  int result = Success;
  DesignSite* pSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if (pSite != NULL)
  {
    DesignSiteRemoveRotamers(pSite);
  }
  else
  {
    (pThis->desSiteCount)++;
    pThis->designSites = (DesignSite*)realloc(pThis->designSites, sizeof(DesignSite) * pThis->desSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->desSiteCount - 1]);
    pSite = StructureGetDesignSite(pThis, pThis->desSiteCount - 1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pSite->pRes = pDestResidue;
    pSite->chnNdx = chainIndex;
    pSite->resNdx = resiIndex;
  }

  // set design types - 20 AA types;
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


int ProteinSiteBuildMutatedRotamersByBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, StringArray* pDesignTypes, StringArray* pPatchTypes)
{
  int result = Success;
  DesignSite* pCurrentDesignSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if (pCurrentDesignSite != NULL)
  {
    DesignSiteRemoveRotamers(pCurrentDesignSite);
  }
  else
  {
    (pThis->desSiteCount)++;
    pThis->designSites = (DesignSite*)realloc(pThis->designSites, sizeof(DesignSite) * pThis->desSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->desSiteCount - 1]);
    pCurrentDesignSite = StructureGetDesignSite(pThis, pThis->desSiteCount - 1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
    pCurrentDesignSite->pRes = pDestResidue;
    pCurrentDesignSite->chnNdx = chainIndex;
    pCurrentDesignSite->resNdx = resiIndex;
  }

  result = RotamerSetOfProteinGenerateByBBdepRotLib(DesignSiteGetRotamers(pCurrentDesignSite), pCurrentDesignSite->pRes, pDesignTypes, pPatchTypes, rotlib, atomParams, resiTopos);
  ResidueSetDesignType(pCurrentDesignSite->pRes, Type_DesType_Mutable);
  return result;
}


int ProteinSiteBuildWildtypeRotamersByBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos)
{
  int result = Success;
  DesignSite* pSite = StructureFindDesignSite(pThis, chainIndex, resiIndex);
  if (pSite != NULL)
  {
    DesignSiteRemoveRotamers(pSite);
  }
  else
  {
    (pThis->desSiteCount)++;
    pThis->designSites = (DesignSite*)realloc(pThis->designSites, sizeof(DesignSite) * pThis->desSiteCount);
    DesignSiteCreate(&pThis->designSites[pThis->desSiteCount - 1]);
    pSite = StructureGetDesignSite(pThis, pThis->desSiteCount - 1);
    Chain* pDestChain = StructureGetChain(pThis, chainIndex);
    Residue* pResidue = ChainGetResidue(pDestChain, resiIndex);
    pSite->pRes = pResidue;
    pSite->chnNdx = chainIndex;
    pSite->resNdx = resiIndex;
  }

  // set native design types;
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


int ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(Structure* pStruct, int chnNdx, int resNdx, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos)
{
  int result = Success;
  DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pSite != NULL)
  {
    DesignSiteRemoveRotamers(pSite);
  }
  else
  {
    (pStruct->desSiteCount)++;
    pStruct->designSites = (DesignSite*)realloc(pStruct->designSites, sizeof(DesignSite) * pStruct->desSiteCount);
    DesignSiteCreate(&pStruct->designSites[pStruct->desSiteCount - 1]);
    pSite = StructureGetDesignSite(pStruct, pStruct->desSiteCount - 1);
    Chain* pDestChn = StructureGetChain(pStruct, chnNdx);
    Residue* pDestRes = ChainGetResidue(pDestChn, resNdx);
    pSite->pRes = pDestRes;
    pSite->chnNdx = chnNdx;
    pSite->resNdx = resNdx;
  }


  if (ResidueGetDesignType(pSite->pRes) == Type_DesType_Fixed ||
    ResidueGetDesignType(pSite->pRes) == Type_DesType_NatRot)
  {
    ;
  }
  else
  {
    StringArray designTypes;
    StringArray patchTypes;
    StringArrayCreate(&designTypes);
    StringArrayCreate(&patchTypes);
    if (ResidueGetDesignType(pSite->pRes) == Type_DesType_Mutable)
    {
      int len = strlen(pSite->pRes->designAATypes);
      char* desStr = "ACDEFGHIKLMNPQRSTVWY";
      if (!len)
      {
        for (int i = 0;i < 20;i++)
        {
          char aa = desStr[i];
          char rotType[MAX_LEN_RES_NAME + 1];
          AA1ToAA3(aa, rotType);
          if (aa == 'C' && FLAG_EXCL_CYS_ROTS) continue;
          else if (aa == 'H')
          {
            StringArrayAppend(&designTypes, "HSD"); StringArrayAppend(&patchTypes, "");
            StringArrayAppend(&designTypes, "HSE"); StringArrayAppend(&patchTypes, "");
          }
          else if (aa != 'X') { StringArrayAppend(&designTypes, rotType); StringArrayAppend(&patchTypes, ""); }
        }
      }
      else
      {
        for (int i = 0;i < len;i++)
        {
          char aa = pSite->pRes->designAATypes[i];
          char rotType[MAX_LEN_RES_NAME + 1];
          AA1ToAA3(aa, rotType);
          if (aa == 'C' && FLAG_EXCL_CYS_ROTS) continue;
          else if (aa == 'H')
          {
            StringArrayAppend(&designTypes, "HSD"); StringArrayAppend(&patchTypes, "");
            StringArrayAppend(&designTypes, "HSE"); StringArrayAppend(&patchTypes, "");
          }
          else if (aa != 'X') { StringArrayAppend(&designTypes, rotType); StringArrayAppend(&patchTypes, ""); }
        }
      }
    }
    else if (ResidueGetDesignType(pSite->pRes) == Type_DesType_Catalytic)
    {
      int len = strlen(pSite->pRes->designAATypes);
      if (!len)
      {
        StringArrayAppend(&designTypes, ResidueGetName(pSite->pRes)); StringArrayAppend(&patchTypes, "");
      }
      else
      {
        for (int i = 0;i < len;i++)
        {
          char aa = pSite->pRes->designAATypes[i];
          char rotType[MAX_LEN_RES_NAME + 1];
          AA1ToAA3(aa, rotType);
          if (aa == 'H')
          {
            StringArrayAppend(&designTypes, "HSD"); StringArrayAppend(&patchTypes, "");
            StringArrayAppend(&designTypes, "HSE"); StringArrayAppend(&patchTypes, "");
          }
          else if (aa != 'X') { StringArrayAppend(&designTypes, rotType); StringArrayAppend(&patchTypes, ""); }
        }
      }
    }
    else if (ResidueGetDesignType(pSite->pRes) == Type_DesType_Repackable)
    {
      if (IsResidueHistidine(ResidueGetName(pSite->pRes)))
      {
        StringArrayAppend(&designTypes, "HSD"); StringArrayAppend(&patchTypes, "");
        StringArrayAppend(&designTypes, "HSE"); StringArrayAppend(&patchTypes, "");
      }
      else StringArrayAppend(&designTypes, ResidueGetName(pSite->pRes)); StringArrayAppend(&patchTypes, "");
    }
    result = RotamerSetOfProteinGenerateByBBdepRotLib(DesignSiteGetRotamers(pSite), pSite->pRes, &designTypes, &patchTypes, rotlib, atomParams, resiTopos);
    StringArrayDestroy(&patchTypes);
    StringArrayDestroy(&designTypes);
  }

  return result;
}


int StructureBuildAllRotamersByBBdepRotLib(Structure* pStructure, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos)
{
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL) continue;
    for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
    {
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if (pResidue->desType != Type_DesType_Fixed) continue;
      if (FLAG_WILDTYPE_ONLY) ResidueSetDesignType(pResidue, Type_DesType_Repackable);
      else ResidueSetDesignType(pResidue, Type_DesType_Mutable);
      ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
      if (FLAG_USE_INPUT_SC && pResidue->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
      if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
    }
  }

  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL) continue;
    for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
    {
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if (pResidue->desType != Type_DesType_Fixed) continue;
      BOOL interResi = FALSE;
      for (int k = 0; k < StructureGetChainCount(pStructure); k++)
      {
        if (k == i)continue;
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) != Type_Chain_Protein)continue;
        if (strstr(DES_CHAINS, ChainGetName(pChainK)) == 0)continue;
        for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
        {
          Residue* pResidueKS = ChainGetResidue(pChainK, s);
          if (AtomArrayCalcMinDistance(&pResidue->atoms, &pResidueKS->atoms) < CUT_PPI_DIST_SHELL1)
          {
            interResi = TRUE;
            break;
          }
        }
        if (interResi) break;
      }
      if (interResi)
      {
        ResidueSetDesignType(pResidue, Type_DesType_Repackable);
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
        if (FLAG_USE_INPUT_SC && pResidue->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
        if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
      }
    }
  }

  return Success;
}


int StructureBuildPPIRotamersByBBdepRotLib(Structure* pStructure, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos)
{
  if (StructureGetChainCount(pStructure) < 2)
  {
    printf("there is only one chain in the whole structure, no protein-protein interface found\n");
    exit(ValueError);
  }

  IntArray* arrayFlagMutated = (IntArray*)malloc(sizeof(IntArray) * StructureGetChainCount(pStructure));
  IntArray* arrayFlagRotameric = (IntArray*)malloc(sizeof(IntArray) * StructureGetChainCount(pStructure));
  for (int i = 0;i < StructureGetChainCount(pStructure);i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    IntArrayCreate(&arrayFlagMutated[i], pChainI->residueNum);
    IntArrayCreate(&arrayFlagRotameric[i], pChainI->residueNum);
    for (int j = 0;j < pChainI->residueNum;j++)
    {
      IntArraySet(&arrayFlagMutated[i], j, 0);
      IntArraySet(&arrayFlagRotameric[i], j, 0);
    }
  }

  //1. specify mutated and rotameric positions
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    //if(ChainGetType(pChainI)!=Type_Chain_Protein) continue;
    for (int k = i + 1; k < StructureGetChainCount(pStructure); k++)
    {
      Chain* pChainK = StructureGetChain(pStructure, k);
      //if(ChainGetType(pChainK)!=Type_Chain_Protein) continue;
      for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
      {
        Residue* pResiIJ = ChainGetResidue(pChainI, j);
        for (int s = 0; s < ChainGetResidueCount(pChainK); s++)
        {
          Residue* pResiKS = ChainGetResidue(pChainK, s);
          if (AtomArrayCalcMinDistance(&pResiIJ->atoms, &pResiKS->atoms) < CUT_PPI_DIST_SHELL1)
          {
            if (strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL) IntArraySet(&arrayFlagMutated[i], j, 1);
            else IntArraySet(&arrayFlagRotameric[i], j, 1);
            if (strstr(DES_CHAINS, ChainGetName(pChainK)) != NULL) IntArraySet(&arrayFlagMutated[k], s, 1);
            else IntArraySet(&arrayFlagRotameric[k], s, 1);
          }
          else if (AtomArrayCalcMinDistance(&pResiIJ->atoms, &pResiKS->atoms) < CUT_PPI_DIST_SHELL2)
          {
            if (strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL) IntArraySet(&arrayFlagRotameric[i], j, 1);
            if (strstr(DES_CHAINS, ChainGetName(pChainK)) != NULL) IntArraySet(&arrayFlagRotameric[k], s, 1);
          }
        }
      }
    }
  }

  //2. create rots on the design chains
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL)
    {
      for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
      {
        Residue* pResidue = ChainGetResidue(pChainI, j);
        if (pResidue->desType != Type_DesType_Fixed) continue;
        if (IntArrayGet(&arrayFlagMutated[i], j) == 1)
        {
          if (FLAG_WILDTYPE_ONLY) ResidueSetDesignType(pResidue, Type_DesType_Repackable);
          else ResidueSetDesignType(pResidue, Type_DesType_Mutable);
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
          if (FLAG_USE_INPUT_SC && pResidue->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
          if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
        }
        else if (IntArrayGet(&arrayFlagRotameric[i], j) == 1)
        {
          ResidueSetDesignType(pResidue, Type_DesType_Repackable);
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
          if (FLAG_USE_INPUT_SC && pResidue->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
          if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
        }
      }
    }
  }

  //3. create rots on the non-design chains
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL)
    {
      for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
      {
        Residue* pResidue = ChainGetResidue(pChainI, j);
        if (pResidue->desType != Type_DesType_Fixed) continue;
        if (IntArrayGet(&arrayFlagRotameric[i], j) == 1)
        {
          ResidueSetDesignType(pResidue, Type_DesType_Repackable);
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
          if (FLAG_USE_INPUT_SC && pResidue->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
          if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
        }
      }
    }
  }

  for (int i = 0;i < StructureGetChainCount(pStructure);i++)
  {
    IntArrayDestroy(&arrayFlagMutated[i]);
    IntArrayDestroy(&arrayFlagRotameric[i]);
  }
  free(arrayFlagMutated);
  arrayFlagMutated = NULL;
  free(arrayFlagRotameric);
  arrayFlagRotameric = NULL;

  return Success;
}


int StructureBuildResfileRotamersByBBdepRotLib(Structure* pStructure, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* resfile)
{
  FileReader fr;
  int result = FileReaderCreate(&fr, resfile);
  StringArray strings;
  StringArrayCreate(&strings);

  // 1. read sites from file
  if (!FAILED(result))
  {
    char buffer[MAX_LEN_ONE_LINE_CONTENT + 1];
    while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
    {
      char keyword[MAX_LEN_ONE_LINE_CONTENT + 1];
      ExtractFirstStringFromSourceString(keyword, buffer);
      if (!strcmp(keyword, "SITES_CATALYTIC_START"))
      {
        while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
        {
          StringArraySplitString(&strings, buffer, ' ');
          if (!strcmp(StringArrayGet(&strings, 0), "SITES_CATALYTIC_END")) break;
          char* chnName = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStructure, chnName);
          int resiIndex = -1;
          ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
          if (resiIndex != -1)
          {
            if (StringArrayGetCount(&strings) == 3)
            {
              // deal with native rotameric conformer for catalytic design sites
              if (strcmp(StringArrayGet(&strings, 2), "NATROT") == 0)
              {
                ChainGetResidue(pChain, resiIndex)->desType = Type_DesType_NatRot;
              }
              else
              {
                ChainGetResidue(pChain, resiIndex)->desType = Type_DesType_Catalytic;
                strcpy(ChainGetResidue(pChain, resiIndex)->designAATypes, StringArrayGet(&strings, 2));
              }
            }
          }
        }
      }
      else if (!strcmp(keyword, "SITES_DESIGN_START"))
      {
        while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
        {
          StringArraySplitString(&strings, buffer, ' ');
          if (!strcmp(StringArrayGet(&strings, 0), "SITES_DESIGN_END")) break;
          char* chnname = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStructure, chnname);
          int resiIndex = -1;
          ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
          if (resiIndex != -1)
          {
            ResidueSetDesignType(ChainGetResidue(pChain, resiIndex), Type_DesType_Mutable);
            if (StringArrayGetCount(&strings) == 3)
            {
              strcpy(ChainGetResidue(pChain, resiIndex)->designAATypes, StringArrayGet(&strings, 2));
            }
          }
        }
      }
      else if (!strcmp(keyword, "SITES_REPACK_START"))
      {
        while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
        {
          StringArraySplitString(&strings, buffer, ' ');
          if (!strcmp(StringArrayGet(&strings, 0), "SITES_REPACK_END")) break;
          char* chnname = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStructure, chnname);
          int resiIndex = -1;
          ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
          if (resiIndex != -1)
          {
            ResidueSetDesignType(ChainGetResidue(pChain, resiIndex), Type_DesType_Repackable);
          }
        }
      }
      else if (!strcmp(keyword, "SITES_FIX_START"))
      {
        while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
        {
          StringArraySplitString(&strings, buffer, ' ');
          if (!strcmp(StringArrayGet(&strings, 0), "SITES_FIX_END")) break;
          char* chnname = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStructure, chnname);
          int resiIndex = -1;
          ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
          if (resiIndex != -1)
          {
            ResidueSetDesignType(ChainGetResidue(pChain, resiIndex), Type_DesType_Fixed);
          }
        }
      }
    }
    FileReaderDestroy(&fr);
  }
  StringArrayDestroy(&strings);

  // 2. deal with site design type
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    for (int j = 0; j < ChainGetResidueCount(pChainI);j++)
    {
      Residue* pResi = ChainGetResidue(pChainI, j);
      if (ResidueGetDesignType(pResi) == Type_DesType_Fixed)
      {
        continue;
      }
      if (ResidueGetDesignType(pResi) == Type_DesType_Mutable && FLAG_WILDTYPE_ONLY)
      {
        ResidueSetDesignType(pResi, Type_DesType_Repackable);
      }
      if (ResidueGetDesignType(pResi) == Type_DesType_Repackable)
      {
        if (!strcmp(ResidueGetName(pResi), "ALA") || !strcmp(ResidueGetName(pResi), "GLY"))
        {
          ResidueSetDesignType(pResi, Type_DesType_Fixed);
        }
      }
    }
  }

  // 3. build rotamers for each design site
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    for (int j = 0; j < ChainGetResidueCount(pChainI);j++)
    {
      Residue* pResi = ChainGetResidue(pChainI, j);
      if (ResidueGetDesignType(pResi) != Type_DesType_Fixed)
      {
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
        if (FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
        if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
      }
    }
  }

  return Success;
}


int StructureBuildCatalyticRotamersByBBdepRotLib(Structure* pThis, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* resfile)
{
  FileReader fr;
  if (!FAILED(FileReaderCreate(&fr, resfile)))
  {
    BOOL flagCataSite = FALSE;
    char buffer[MAX_LEN_ONE_LINE_CONTENT + 1];
    while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
    {
      char keyword[MAX_LEN_ONE_LINE_CONTENT + 1];
      ExtractFirstStringFromSourceString(keyword, buffer);
      if (!strcmp(keyword, "SITES_CATALYTIC_START"))
      {
        flagCataSite = TRUE;
        StringArray strings;
        StringArrayCreate(&strings);
        while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
        {
          StringArraySplitString(&strings, buffer, ' ');
          if (strcmp(StringArrayGet(&strings, 0), "SITES_CATALYTIC_END") == 0) break;
          char* chnName = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pThis, chnName);
          int resiIndex = -1;
          ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
          if (resiIndex != -1)
          {
            if (StringArrayGetCount(&strings) == 3)
            {
              if (strcmp(StringArrayGet(&strings, 2), "NATROT") == 0)
              {
                ChainGetResidue(pChain, resiIndex)->desType = Type_DesType_NatRot;
              }
              else
              {
                ChainGetResidue(pChain, resiIndex)->desType = Type_DesType_Catalytic;
                strcpy(ChainGetResidue(pChain, resiIndex)->designAATypes, StringArrayGet(&strings, 2));
              }
            }
          }
        }
        StringArrayDestroy(&strings);
      }
    }
    FileReaderDestroy(&fr);

    if (flagCataSite == FALSE)
    {
      char errMsg[MAX_LEN_ERR_MSG + 1];
      sprintf(errMsg, "in file %s line %d, catalytic sites were not designated. Note: user should specify catalytic sites using keyword 'SITES_CATALYTIC_START' "
        "and 'SITES_CATALYTIC_END' in the resfile; in most cases, catalytic site are not designable", __FILE__, __LINE__);
      TraceError(errMsg, FormatError);
      return FormatError;
    }

    for (int i = 0; i < StructureGetChainCount(pThis); i++)
    {
      Chain* pChainI = StructureGetChain(pThis, i);
      if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
      for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
      {
        Residue* pResi = ChainGetResidue(pChainI, j);
        if (ResidueGetDesignType(pResi) == Type_DesType_Fixed) continue;
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pThis, i, j, rotlib, atomParams, resiTopos);
        if (FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pThis, i, j, resiTopos);
        if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pThis, i, j, resiTopos);
      }
    }
  }
  else
  {
    char errMsg[MAX_LEN_ERR_MSG + 1];
    sprintf(errMsg, "in file %s line %d, catalytic sites were not designated. Note: user should specify catalytic sites using keyword 'SITES_CATALYTIC_START' "
      "and 'SITES_CATALYTIC_END' in the resfile; in most cases, a catalytic site is not designable", __FILE__, __LINE__);
    TraceError(errMsg, FormatError);
    return FormatError;
  }

  return Success;
}


int StructureBuildPLIShell1RotamersByBBdepRotLib(Structure* pStructure, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* resfile)
{
  FileReader fr;
  int result = FileReaderCreate(&fr, resfile);
  int mutaSiteCount = 0;
  if (!FAILED(result))
  {
    char buffer[MAX_LEN_ONE_LINE_CONTENT + 1];
    while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
    {
      char keyword[MAX_LEN_ONE_LINE_CONTENT + 1];
      ExtractFirstStringFromSourceString(keyword, buffer);
      if (!strcmp(keyword, "SITES_DESIGN_START"))
      {
        StringArray strings;
        StringArrayCreate(&strings);
        while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
        {
          StringArraySplitString(&strings, buffer, ' ');
          if (!strcmp(StringArrayGet(&strings, 0), "SITES_DESIGN_END")) break;
          char* chain = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStructure, chain);
          int resiIndex = -1;
          ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
          if (resiIndex != -1)
          {
            ChainGetResidue(pChain, resiIndex)->desType = Type_DesType_Mutable;
            if (StringArrayGetCount(&strings) == 3)
            {
              strcpy(ChainGetResidue(pChain, resiIndex)->designAATypes, StringArrayGet(&strings, 2));
            }
          }
          mutaSiteCount++;
        }
        StringArrayDestroy(&strings);
      }
    }
    FileReaderDestroy(&fr);
  }

  if (mutaSiteCount)
  {
    for (int i = 0; i < StructureGetChainCount(pStructure); i++)
    {
      Chain* pChainI = StructureGetChain(pStructure, i);
      if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
      for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
      {
        Residue* pResi = ChainGetResidue(pChainI, j);
        if (ResidueGetDesignType(pResi) == Type_DesType_Mutable)
        {
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
          if (FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
          if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
        }
      }
    }
  }
  else
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];
    Residue* pSmallMol = NULL;
    StructureFindSmallMol(pStructure, &pSmallMol);
    if (pSmallMol == NULL)
    {
      int result = DataNotExistError;
      sprintf(errMsg, "in file %s line %d, cannot find small molecule", __FILE__, __LINE__);
      TraceError(errMsg, result);
      return result;
    }
    for (int i = 0; i < StructureGetChainCount(pStructure); i++)
    {
      Chain* pChainI = StructureGetChain(pStructure, i);
      if (!strcmp(ResidueGetChainName(pSmallMol), ChainGetName(pChainI))) continue;
      if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
      for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
      {
        Residue* pResi = ChainGetResidue(pChainI, j);
        double minDist = AtomArrayCalcMinDistance(&pSmallMol->atoms, &pResi->atoms);
        if (minDist > CUT_PLI_DIST_SHELL1) continue;
        if (ResidueGetDesignType(pResi) != Type_DesType_Fixed) continue;
        ResidueSetDesignType(pResi, Type_DesType_Mutable);
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
        if (FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
        if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
      }
    }
  }

  return Success;
}


int StructureBuildPLIShell2RotamersByBBdepRotLib(Structure* pStructure, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* resfile)
{
  FileReader fr;
  int result = FileReaderCreate(&fr, resfile);
  int rotaSiteCount = 0;
  if (!FAILED(result))
  {
    char buffer[MAX_LEN_ONE_LINE_CONTENT + 1];
    while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
    {
      char keyword[MAX_LEN_ONE_LINE_CONTENT + 1];
      ExtractFirstStringFromSourceString(keyword, buffer);
      if (!strcmp(keyword, "SITES_REPACK_START"))
      {
        StringArray strings;
        StringArrayCreate(&strings);
        while (!FAILED(FileReaderGetNextLine(&fr, buffer)))
        {
          StringArraySplitString(&strings, buffer, ' ');
          if (!strcmp(StringArrayGet(&strings, 0), "SITES_REPACK_END")) break;
          char* chnname = StringArrayGet(&strings, 0);
          int posInChain = atoi(StringArrayGet(&strings, 1));
          Chain* pChain = StructureFindChainByName(pStructure, chnname);
          int resiIndex = -1;
          ChainFindResidueByPosInChain(pChain, posInChain, &resiIndex);
          if (resiIndex != -1) ChainGetResidue(pChain, resiIndex)->desType = Type_DesType_Repackable;
          rotaSiteCount++;
        }
        StringArrayDestroy(&strings);
      }
    }
    FileReaderDestroy(&fr);
  }

  if (rotaSiteCount)
  {
    for (int i = 0; i < StructureGetChainCount(pStructure); i++)
    {
      Chain* pChainI = StructureGetChain(pStructure, i);
      if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
      for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
      {
        Residue* pResi = ChainGetResidue(pChainI, j);
        if (ResidueGetDesignType(pResi) == Type_DesType_Repackable)
        {
          if (!strcmp(ResidueGetName(pResi), "ALA") || !strcmp(ResidueGetName(pResi), "GLY"))
          {
            ResidueSetDesignType(pResi, Type_DesType_Fixed);
            continue;
          }
          ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
          if (FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
          if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
        }
      }
    }
  }
  else
  {
    char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];
    Residue* pSmallMol = NULL;
    StructureFindSmallMol(pStructure, &pSmallMol);
    if (pSmallMol == NULL)
    {
      int result = DataNotExistError;
      sprintf(errMsg, "in file %s line %d, cannot find small molecule", __FILE__, __LINE__);
      TraceError(errMsg, result);
      return result;
    }

    for (int i = 0; i < StructureGetChainCount(pStructure); i++)
    {
      Chain* pChainI = StructureGetChain(pStructure, i);
      if (!strcmp(ResidueGetChainName(pSmallMol), ChainGetName(pChainI))) continue;
      if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
      for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
      {
        Residue* pResi = ChainGetResidue(pChainI, j);
        double minDist = AtomArrayCalcMinDistance(&pSmallMol->atoms, &pResi->atoms);
        if (minDist > CUT_PLI_DIST_SHELL2) continue;
        if (ResidueGetDesignType(pResi) != Type_DesType_Fixed) continue;
        if (!strcmp(ResidueGetName(pResi), "ALA") || !strcmp(ResidueGetName(pResi), "GLY")) continue;
        ResidueSetDesignType(pResi, Type_DesType_Repackable);
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
        if (FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
        if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
      }
    }
  }

  return Success;
}



#define CHECK_ROTAMER_CONF_IN_ROTLIB

BOOL IsNativeRotamerInBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex, ResiTopoSet* pResiTopos, BBdepRotamerLib* pBBdepRotLib, double torsionStd)
{
  char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];
  Chain* pDestChain = StructureGetChain(pThis, chainIndex);
  Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
  if (pDestChain->type == Type_Chain_Protein)
  {
    //do not add rotamer for residue ala and gly
    if (strcmp(ResidueGetName(pDestResidue), "ALA") == 0 || strcmp(ResidueGetName(pDestResidue), "GLY") == 0)
    {
      return TRUE;
    }

    //step1: calculate the X angles for the current residue
    double DUNBRACK_ENERGY = 0.0;
    int torsionCount = 0;
    char resname[MAX_LEN_RES_NAME + 1];
    strcpy(resname, ResidueGetName(pDestResidue));
    if (!strcmp(resname, "ALA")) { torsionCount = 0; }
    else if (!strcmp(resname, "ARG")) { torsionCount = 4; }
    else if (!strcmp(resname, "ASN")) { torsionCount = 2; }
    else if (!strcmp(resname, "ASP")) { torsionCount = 2; }
    else if (!strcmp(resname, "CYS")) { torsionCount = 1; }
    else if (!strcmp(resname, "GLN")) { torsionCount = 3; }
    else if (!strcmp(resname, "GLU")) { torsionCount = 3; }
    else if (!strcmp(resname, "GLY")) { torsionCount = 0; }
    else if (!strcmp(resname, "HSD")) { torsionCount = 2; }
    else if (!strcmp(resname, "HSE")) { torsionCount = 2; }
    else if (!strcmp(resname, "ILE")) { torsionCount = 2; }
    else if (!strcmp(resname, "LEU")) { torsionCount = 2; }
    else if (!strcmp(resname, "LYS")) { torsionCount = 4; }
    else if (!strcmp(resname, "MET")) { torsionCount = 3; }
    else if (!strcmp(resname, "PHE")) { torsionCount = 2; }
    else if (!strcmp(resname, "PRO")) { torsionCount = 2; }
    else if (!strcmp(resname, "SER")) { torsionCount = 1; }
    else if (!strcmp(resname, "THR")) { torsionCount = 1; }
    else if (!strcmp(resname, "TRP")) { torsionCount = 2; }
    else if (!strcmp(resname, "TYR")) { torsionCount = 2; }
    else if (!strcmp(resname, "VAL")) { torsionCount = 1; }

    DoubleArray xangles;
    DoubleArrayCreate(&xangles, 0);
    ResidueTopology resiTop;
    ResidueTopologyCreate(&resiTop);
    ResiTopoSetGet(pResiTopos, resname, &resiTop);
    for (int torsionIndex = 0;torsionIndex < torsionCount;torsionIndex++)
    {
      Type_ProteinAtomOrder desiredAtomBOrder = Type_ProteinAtomOrder_FromInt(torsionIndex);
      Type_ProteinAtomOrder desiredAtomCOrder = Type_ProteinAtomOrder_FromInt(torsionIndex + 1);
      CharmmIC icOfCurrentTorsion;
      CharmmICCreate(&icOfCurrentTorsion);
      BOOL icFound = FALSE;
      for (int icIndex = 0; icIndex < ResidueTopologyGetCharmmICCount(&resiTop); icIndex++)
      {
        Type_ProteinAtomOrder atomBOrder;
        Type_ProteinAtomOrder atomCOrder;
        ResidueTopologyGetCharmmIC(&resiTop, icIndex, &icOfCurrentTorsion);
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
        sprintf(errMsg, "in file %s line %d, cannot find IC of the %dth torsion for residue %s", __FILE__, __LINE__, torsionIndex + 1, resname);
        TraceError(errMsg, DataNotExistError);
        CharmmICDestroy(&icOfCurrentTorsion);
        return DataNotExistError;
      }

      Atom* pAtomA = ResidueGetAtomByName(pDestResidue, CharmmICGetAtomA(&icOfCurrentTorsion));
      Atom* pAtomB = ResidueGetAtomByName(pDestResidue, CharmmICGetAtomB(&icOfCurrentTorsion));
      Atom* pAtomC = ResidueGetAtomByName(pDestResidue, CharmmICGetAtomC(&icOfCurrentTorsion));
      Atom* pAtomD = ResidueGetAtomByName(pDestResidue, CharmmICGetAtomD(&icOfCurrentTorsion));
      double torsion = GetTorsionAngle(&pAtomA->xyz, &pAtomB->xyz, &pAtomC->xyz, &pAtomD->xyz);
      DoubleArrayAppend(&xangles, torsion);

      CharmmICDestroy(&icOfCurrentTorsion);
    }
    ResidueTopologyDestroy(&resiTop);

    //step2: get the phi&psi angles for the current residue
    int binindex = ((int)(pDestResidue->phipsi[0] + 180) / 10) * 36 + (int)(pDestResidue->phipsi[1] + 180) / 10;
    RotLibPhiPsi* pRotLibPhiPsi = &pBBdepRotLib->rotlibphipsis[binindex];
    int rotTypeIndex = -1;
    StringArrayFind(&pRotLibPhiPsi->rotTypes, ResidueGetName(pDestResidue), &rotTypeIndex);
    DoubleArray* pTorsionsArrayForTypeI = pRotLibPhiPsi->torsions[rotTypeIndex];
    DoubleArray* pDeviationsArrayForTypeI = pRotLibPhiPsi->deviations[rotTypeIndex];
    int matchIndex = -1;
    for (int i = 0;i < IntArrayGet(&pRotLibPhiPsi->rotamerCounts, rotTypeIndex);i++)
    {
      DoubleArray* pTorsions = &pTorsionsArrayForTypeI[i];
      DoubleArray* pDeviations = &pDeviationsArrayForTypeI[i];
      BOOL match = TRUE;
      for (int j = 0;j < DoubleArrayGetLength(&xangles);j++)
      {
        double min = DoubleArrayGet(pTorsions, j) - DegToRad(torsionStd);
        double max = DoubleArrayGet(pTorsions, j) + DegToRad(torsionStd);
        double torsion = DoubleArrayGet(&xangles, j);
        double torsionm2pi = torsion - 2 * PI;
        double torsionp2pi = torsion + 2 * PI;
        double torsion2 = torsion;
        if ((strcmp(resname, "PHE") == 0 && j == 1) ||
          (strcmp(resname, "TYR") == 0 && j == 1) ||
          (strcmp(resname, "ASP") == 0 && j == 1) ||
          strcmp(resname, "GLU") == 0 && j == 2)
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
      if (match)
      {
        matchIndex = i;
        break;
      }
    }
    DoubleArrayDestroy(&xangles);

    if (matchIndex != -1)
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }

  }

  return Success;
}

BOOL IsNativeRotamerInBBindRotLib(Structure* pThis, int chainIndex, int resiIndex, ResiTopoSet* pResiTopos, BBindRotamerLib* pBBindRotLib, double torsionStd)
{
  char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];
  Chain* pDestChain = StructureGetChain(pThis, chainIndex);
  Residue* pDestResidue = ChainGetResidue(pDestChain, resiIndex);
  if (pDestChain->type == Type_Chain_Protein)
  {
    //do not add rotamer for residue ala and gly
    if (strcmp(ResidueGetName(pDestResidue), "ALA") == 0 || strcmp(ResidueGetName(pDestResidue), "GLY") == 0)
    {
      return TRUE;
    }

    //step1: calculate the X angles for the current residue
    double DUNBRACK_ENERGY = 0.0;
    int torsionCount = 0;
    char resname[MAX_LEN_RES_NAME + 1];
    strcpy(resname, ResidueGetName(pDestResidue));
    if (!strcmp(resname, "ALA")) { torsionCount = 0; }
    else if (!strcmp(resname, "ARG")) { torsionCount = 4; }
    else if (!strcmp(resname, "ASN")) { torsionCount = 2; }
    else if (!strcmp(resname, "ASP")) { torsionCount = 2; }
    else if (!strcmp(resname, "CYS")) { torsionCount = 1; }
    else if (!strcmp(resname, "GLN")) { torsionCount = 3; }
    else if (!strcmp(resname, "GLU")) { torsionCount = 3; }
    else if (!strcmp(resname, "GLY")) { torsionCount = 0; }
    else if (!strcmp(resname, "HSD")) { torsionCount = 2; }
    else if (!strcmp(resname, "HSE")) { torsionCount = 2; }
    else if (!strcmp(resname, "ILE")) { torsionCount = 2; }
    else if (!strcmp(resname, "LEU")) { torsionCount = 2; }
    else if (!strcmp(resname, "LYS")) { torsionCount = 4; }
    else if (!strcmp(resname, "MET")) { torsionCount = 3; }
    else if (!strcmp(resname, "PHE")) { torsionCount = 2; }
    else if (!strcmp(resname, "PRO")) { torsionCount = 2; }
    else if (!strcmp(resname, "SER")) { torsionCount = 1; }
    else if (!strcmp(resname, "THR")) { torsionCount = 1; }
    else if (!strcmp(resname, "TRP")) { torsionCount = 2; }
    else if (!strcmp(resname, "TYR")) { torsionCount = 2; }
    else if (!strcmp(resname, "VAL")) { torsionCount = 1; }

    DoubleArray xangles;
    DoubleArrayCreate(&xangles, 0);
    ResidueTopology resiTop;
    ResidueTopologyCreate(&resiTop);
    ResiTopoSetGet(pResiTopos, resname, &resiTop);
    for (int torsionIndex = 0;torsionIndex < torsionCount;torsionIndex++)
    {
      Type_ProteinAtomOrder desiredAtomBOrder = Type_ProteinAtomOrder_FromInt(torsionIndex);
      Type_ProteinAtomOrder desiredAtomCOrder = Type_ProteinAtomOrder_FromInt(torsionIndex + 1);
      CharmmIC icOfCurrentTorsion;
      CharmmICCreate(&icOfCurrentTorsion);
      BOOL icFound = FALSE;
      for (int icIndex = 0; icIndex < ResidueTopologyGetCharmmICCount(&resiTop); icIndex++)
      {
        Type_ProteinAtomOrder atomBOrder;
        Type_ProteinAtomOrder atomCOrder;
        ResidueTopologyGetCharmmIC(&resiTop, icIndex, &icOfCurrentTorsion);
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
        sprintf(errMsg, "in file %s line %d, cannot find IC of the %dth torsion for residue %s", __FILE__, __LINE__, torsionIndex + 1, resname);
        TraceError(errMsg, DataNotExistError);
        CharmmICDestroy(&icOfCurrentTorsion);
        return DataNotExistError;
      }

      Atom* pAtomA = ResidueGetAtomByName(pDestResidue, CharmmICGetAtomA(&icOfCurrentTorsion));
      Atom* pAtomB = ResidueGetAtomByName(pDestResidue, CharmmICGetAtomB(&icOfCurrentTorsion));
      Atom* pAtomC = ResidueGetAtomByName(pDestResidue, CharmmICGetAtomC(&icOfCurrentTorsion));
      Atom* pAtomD = ResidueGetAtomByName(pDestResidue, CharmmICGetAtomD(&icOfCurrentTorsion));
      double torsion = GetTorsionAngle(&pAtomA->xyz, &pAtomB->xyz, &pAtomC->xyz, &pAtomD->xyz);
      DoubleArrayAppend(&xangles, torsion);

      CharmmICDestroy(&icOfCurrentTorsion);
    }
    ResidueTopologyDestroy(&resiTop);


    //step2: check the torsions in rotamer library
    int rotTypeIndex = -1, matchIndex = -1;
    StringArrayFind(&pBBindRotLib->resTypeNames, ResidueGetName(pDestResidue), &rotTypeIndex);
    int rotCount = IntArrayGet(&pBBindRotLib->rotCounts, rotTypeIndex);
    DoubleArray torsions;
    DoubleArrayCreate(&torsions, 0);
    for (int i = 0;i < rotCount;i++)
    {
      BBindRotamerLibGet(pBBindRotLib, ResidueGetName(pDestResidue), i, &torsions);
      BOOL match = TRUE;
      for (int j = 0;j < DoubleArrayGetLength(&xangles);j++)
      {
        double min = DoubleArrayGet(&torsions, j) - DegToRad(torsionStd);
        double max = DoubleArrayGet(&torsions, j) + DegToRad(torsionStd);
        double torsion = DoubleArrayGet(&xangles, j);
        double torsionm2pi = torsion - 2 * PI;
        double torsionp2pi = torsion + 2 * PI;
        double torsion2 = torsion;
        if ((strcmp(resname, "PHE") == 0 && j == 1) ||
          (strcmp(resname, "TYR") == 0 && j == 1) ||
          (strcmp(resname, "ASP") == 0 && j == 1) ||
          strcmp(resname, "GLU") == 0 && j == 2)
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
      if (match)
      {
        matchIndex = i;
        break;
      }
    }
    DoubleArrayDestroy(&torsions);
    DoubleArrayDestroy(&xangles);

    if (matchIndex != -1)
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }

  }

  return FALSE;
}


#define GENERAL_FUNCTIONS_FOR_BBDEP_AND_BBIND

int ProteinSiteWriteRotamers(Structure* pStructure, int chainIndex, int resiIndex, char* filepath)
{
  DesignSite* pSite = StructureFindDesignSite(pStructure, chainIndex, resiIndex);
  RotamerSet* pSet = DesignSiteGetRotamers(pSite);
  FILE* pOut = fopen(filepath, "w");
  if (pOut != NULL)
  {
    for (int i = 0; i < RotamerSetGetCount(pSet); i++)
    {
      Rotamer* pRotamer = RotamerSetGet(pSet, i);
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


int ProteinSiteBuildNativeRotamer(Structure* pThis, int chnNdx, int resNdx, ResiTopoSet* pResiTopos)
{
  Chain* pChain = StructureGetChain(pThis, chnNdx);
  Residue* pResi = ChainGetResidue(pChain, resNdx);
  if (!strcmp(ResidueGetName(pResi), "ALA") || !strcmp(ResidueGetName(pResi), "GLY")) return Success;

  if (ChainGetType(pChain) == Type_Chain_Protein)
  {
    DesignSite* pSite = StructureFindDesignSite(pThis, chnNdx, resNdx);
    if (!pSite)
    {
      (pThis->desSiteCount)++;
      pThis->designSites = (DesignSite*)realloc(pThis->designSites, sizeof(DesignSite) * pThis->desSiteCount);
      DesignSiteCreate(&pThis->designSites[pThis->desSiteCount - 1]);
      pSite = StructureGetDesignSite(pThis, pThis->desSiteCount - 1);
      pSite->pRes = pResi;
      pSite->chnNdx = chnNdx;
      pSite->resNdx = resNdx;
    }

    RotamerSet* pSetI = DesignSiteGetRotamers(pSite);
    Rotamer tempRot;
    RotamerCreate(&tempRot);
    RotamerSetType(&tempRot, ResidueGetName(pSite->pRes));
    RotamerSetChainName(&tempRot, ResidueGetChainName(pSite->pRes));
    RotamerSetPosInChain(&tempRot, ResidueGetPosInChain(pSite->pRes));
    Rotamer* pRepresentative = RotamerSetGetRepresentative(pSetI, RotamerGetType(&tempRot));
    if (pRepresentative != NULL)
    {
      AtomArrayCopy(&tempRot.atoms, &pRepresentative->atoms);
      for (int i = 0; i < RotamerGetAtomCount(&tempRot); i++)
      {
        Atom* pAtom = RotamerGetAtom(&tempRot, i);
        pAtom->xyz = ResidueGetAtomByName(pSite->pRes, AtomGetName(pAtom))->xyz;
      }
      BondSetCopy(&tempRot.bonds, &pRepresentative->bonds);
      XYZArrayResize(&tempRot.xyzs, RotamerGetAtomCount(pRepresentative));
      for (int i = 0; i < XYZArrayGetLength(&tempRot.xyzs); i++)
      {
        XYZArraySet(&tempRot.xyzs, i, &AtomArrayGet(&tempRot.atoms, i)->xyz);
      }
      RotamerSetDunbrack(&tempRot, ResidueGetDunbrack(pResi));
      DoubleArrayCopy(&tempRot.Xs, &pResi->Xs);
      RotamerSetAdd(&pSite->rots, &tempRot);
    }
    else
    { // cannot find a representative rotamer
      if (ResidueGetDesignType(pResi) == Type_DesType_NatRot)
      {
        AtomArrayCopy(&tempRot.atoms, &pResi->atoms);
        BondSetCopy(&tempRot.bonds, &pResi->bonds);
        XYZArrayResize(&tempRot.xyzs, AtomArrayGetCount(&tempRot.atoms));
        for (int i = 0; i < XYZArrayGetLength(&tempRot.xyzs); i++)
        {
          XYZArraySet(&tempRot.xyzs, i, &AtomArrayGet(&tempRot.atoms, i)->xyz);
        }
        RotamerSetDunbrack(&tempRot, ResidueGetDunbrack(pResi));
        DoubleArrayCopy(&tempRot.Xs, &pResi->Xs);
        RotamerSetAdd(&pSite->rots, &tempRot);
      }
    }
    RotamerDestroy(&tempRot);

    if (!strcmp(ResidueGetName(pSite->pRes), "HSD"))
    {
      if ((pRepresentative = RotamerSetGetRepresentative(pSetI, "HSE")))
      {
        Residue newResi;
        ResidueCreate(&newResi);
        ResidueSetName(&newResi, "HSE");
        AtomArrayCopy(&newResi.atoms, &pRepresentative->atoms);
        for (int i = 0; i < ResidueGetAtomCount(&newResi); i++)
        {
          Atom* pAtom = ResidueGetAtom(&newResi, i);
          if (!pAtom->isBBAtom && AtomIsHydrogen(pAtom))
          {
            pAtom->isXyzValid = FALSE;
            continue;
          }
          pAtom->xyz = ResidueGetAtomByName(pSite->pRes, AtomGetName(pAtom))->xyz;
        }
        ResidueCalcAllAtomXYZ(&newResi, pResiTopos, NULL, NULL);

        RotamerCreate(&tempRot);
        RotamerSetType(&tempRot, "HSE");
        RotamerSetChainName(&tempRot, ResidueGetChainName(pSite->pRes));
        RotamerSetPosInChain(&tempRot, ResidueGetPosInChain(pSite->pRes));
        AtomArrayCopy(&tempRot.atoms, &pRepresentative->atoms);
        for (int i = 0; i < RotamerGetAtomCount(&tempRot); i++)
        {
          Atom* pAtom = RotamerGetAtom(&tempRot, i);
          pAtom->xyz = ResidueGetAtomByName(&newResi, AtomGetName(pAtom))->xyz;
        }
        BondSetCopy(&tempRot.bonds, &pRepresentative->bonds);

        XYZArrayResize(&tempRot.xyzs, RotamerGetAtomCount(pRepresentative));
        for (int i = 0; i < XYZArrayGetLength(&tempRot.xyzs); i++)
        {
          XYZArraySet(&tempRot.xyzs, i, &AtomArrayGet(&tempRot.atoms, i)->xyz);
        }
        RotamerSetDunbrack(&tempRot, ResidueGetDunbrack(pResi));
        DoubleArrayCopy(&tempRot.Xs, &pResi->Xs);
        RotamerSetAdd(&pSite->rots, &tempRot);

        RotamerDestroy(&tempRot);
        ResidueDestroy(&newResi);
      }
    }
    else if (!strcmp(ResidueGetName(pSite->pRes), "HSE"))
    {
      if ((pRepresentative = RotamerSetGetRepresentative(pSetI, "HSD")))
      {
        Residue newResi;
        ResidueCreate(&newResi);
        ResidueSetName(&newResi, "HSD");
        AtomArrayCopy(&newResi.atoms, &pRepresentative->atoms);
        for (int i = 0; i < ResidueGetAtomCount(&newResi); i++)
        {
          Atom* pAtom = ResidueGetAtom(&newResi, i);
          if (!pAtom->isBBAtom && AtomIsHydrogen(pAtom))
          {
            pAtom->isXyzValid = FALSE;
            continue;
          }
          pAtom->xyz = ResidueGetAtomByName(pSite->pRes, AtomGetName(pAtom))->xyz;
        }
        ResidueCalcAllAtomXYZ(&newResi, pResiTopos, NULL, NULL);

        RotamerCreate(&tempRot);
        RotamerSetType(&tempRot, "HSD");
        RotamerSetChainName(&tempRot, ResidueGetChainName(pSite->pRes));
        RotamerSetPosInChain(&tempRot, ResidueGetPosInChain(pSite->pRes));

        AtomArrayCopy(&tempRot.atoms, &pRepresentative->atoms);
        for (int i = 0; i < RotamerGetAtomCount(&tempRot); i++)
        {
          Atom* pAtom = RotamerGetAtom(&tempRot, i);
          pAtom->xyz = ResidueGetAtomByName(&newResi, AtomGetName(pAtom))->xyz;
        }
        BondSetCopy(&tempRot.bonds, &pRepresentative->bonds);

        XYZArrayResize(&tempRot.xyzs, RotamerGetAtomCount(pRepresentative));
        for (int i = 0; i < XYZArrayGetLength(&tempRot.xyzs); i++)
        {
          XYZArraySet(&tempRot.xyzs, i, &AtomArrayGet(&tempRot.atoms, i)->xyz);
        }
        RotamerSetDunbrack(&tempRot, ResidueGetDunbrack(pResi));
        DoubleArrayCopy(&tempRot.Xs, &pResi->Xs);
        RotamerSetAdd(&pSite->rots, &tempRot);
        RotamerDestroy(&tempRot);
        ResidueDestroy(&newResi);
      }
    }
  }

  return Success;
}


int ProteinSiteBuildFlippedNativeRotamer(Structure* pStructure, int chainIndex, int resiIndex, ResiTopoSet* pResiTopos)
{
  Chain* pChain = StructureGetChain(pStructure, chainIndex);
  Residue* pResidue = ChainGetResidue(pChain, resiIndex);
  if (pChain->type == Type_Chain_Protein)
  {
    if (strcmp(ResidueGetName(pResidue), "ASN") != 0 && strcmp(ResidueGetName(pResidue), "GLN") != 0 &&
      strcmp(ResidueGetName(pResidue), "HSD") != 0 && strcmp(ResidueGetName(pResidue), "HSE") != 0)
    {
      return Success;
    }

    DesignSite* pSite = StructureFindDesignSite(pStructure, chainIndex, resiIndex);
    RotamerSet* pSet = DesignSiteGetRotamers(pSite);
    int rotCount = RotamerSetGetCount(pSet);
    for (int i = 0; i < rotCount; i++)
    {
      Rotamer* pRotamer = RotamerSetGet(pSet, i);
      Rotamer* pRepresentative = RotamerSetGetRepresentative(pSet, RotamerGetType(pRotamer));
      Rotamer tempRotamer;
      RotamerCreate(&tempRotamer);
      RotamerCopy(&tempRotamer, pRepresentative);
      if (strcmp(RotamerGetType(pRotamer), "ASN") == 0)
      {
        int index1, index2;
        RotamerFindAtom(&tempRotamer, "ND2", &index1);
        RotamerFindAtom(&tempRotamer, "OD1", &index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer, index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer, index2)->xyz;
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ2;
        RotamerGetAtom(&tempRotamer, index2)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs, index2, &tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos, "ASN", &resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo, "HD21", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer, "HD21", &index1);
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ1);
        ResidueTopologyFindCharmmIC(&resiTopo, "HD22", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer, "HD22", &index1);
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      else if (strcmp(RotamerGetType(pRotamer), "GLN") == 0)
      {
        int index1, index2;
        RotamerFindAtom(&tempRotamer, "NE2", &index1);
        RotamerFindAtom(&tempRotamer, "OE1", &index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer, index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer, index2)->xyz;
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ2;
        RotamerGetAtom(&tempRotamer, index2)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs, index2, &tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos, "GLN", &resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo, "HE21", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer, "HE21", &index1);
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ1);
        ResidueTopologyFindCharmmIC(&resiTopo, "HE22", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer, "HE22", &index1);
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      else if (strcmp(RotamerGetType(pRotamer), "HSD") == 0)
      {
        int index1, index2;
        RotamerFindAtom(&tempRotamer, "CD2", &index1);
        RotamerFindAtom(&tempRotamer, "ND1", &index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer, index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer, index2)->xyz;
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ2;
        RotamerGetAtom(&tempRotamer, index2)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs, index2, &tempXYZ1);
        RotamerFindAtom(&tempRotamer, "NE2", &index1);
        RotamerFindAtom(&tempRotamer, "CE1", &index2);
        tempXYZ1 = RotamerGetAtom(&tempRotamer, index1)->xyz;
        tempXYZ2 = RotamerGetAtom(&tempRotamer, index2)->xyz;
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ2;
        RotamerGetAtom(&tempRotamer, index2)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs, index2, &tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos, "HSD", &resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo, "HD1", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer, "HD1", &index1);
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      else if (strcmp(RotamerGetType(pRotamer), "HSE") == 0)
      {
        int index1, index2;
        RotamerFindAtom(&tempRotamer, "CD2", &index1);
        RotamerFindAtom(&tempRotamer, "ND1", &index2);
        XYZ tempXYZ1 = RotamerGetAtom(&tempRotamer, index1)->xyz;
        XYZ tempXYZ2 = RotamerGetAtom(&tempRotamer, index2)->xyz;
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ2;
        RotamerGetAtom(&tempRotamer, index2)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs, index2, &tempXYZ1);
        RotamerFindAtom(&tempRotamer, "NE2", &index1);
        RotamerFindAtom(&tempRotamer, "CE1", &index2);
        tempXYZ1 = RotamerGetAtom(&tempRotamer, index1)->xyz;
        tempXYZ2 = RotamerGetAtom(&tempRotamer, index2)->xyz;
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ2;
        RotamerGetAtom(&tempRotamer, index2)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ2);
        XYZArraySet(&tempRotamer.xyzs, index2, &tempXYZ1);
        ResidueTopology resiTopo;
        CharmmIC ic;
        ResidueTopologyCreate(&resiTopo);
        CharmmICCreate(&ic);
        ResiTopoSetGet(pResiTopos, "HSE", &resiTopo);
        ResidueTopologyFindCharmmIC(&resiTopo, "HE2", &ic);
        GetFourthAtom(&RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomA(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomB(&ic))->xyz,
          &RotamerGetAtomByName(&tempRotamer, CharmmICGetAtomC(&ic))->xyz,
          ic.icParam,
          &tempXYZ1);
        RotamerFindAtom(&tempRotamer, "HE2", &index1);
        RotamerGetAtom(&tempRotamer, index1)->xyz = tempXYZ1;
        XYZArraySet(&tempRotamer.xyzs, index1, &tempXYZ1);
        ResidueTopologyDestroy(&resiTopo);
        CharmmICDestroy(&ic);
      }
      RotamerSetAdd(pSet, &tempRotamer);
      RotamerDestroy(&tempRotamer);
    }
  }

  return Success;
}


int ProteinSiteExpandHydroxylRotamers(Structure* pStructure, int chainIndex, int resiIndex, ResiTopoSet* pTopos)
{
  Chain* pChain = StructureGetChain(pStructure, chainIndex);
  if (pChain->type == Type_Chain_Protein)
  {
    DesignSite* pSite = StructureFindDesignSite(pStructure, chainIndex, resiIndex);
    RotamerSet* pSet = DesignSiteGetRotamers(pSite);
    int rotamerCount = RotamerSetGetCount(pSet); // the rotamer set will be expanded, we need to record the current count
    ResidueTopology tops;
    CharmmIC ics;
    ResidueTopologyCreate(&tops);
    CharmmICCreate(&ics);
    if (RotamerSetGetRepresentative(pSet, "SER"))
    {
      ResiTopoSetGet(pTopos, "SER", &tops);
      ResidueTopologyFindCharmmIC(&tops, "HG", &ics);
      double icParaX = ics.icParam[2];
      int addedCount = 0;
      Rotamer tempRot;
      RotamerCreate(&tempRot);
      for (int j = 0; j < rotamerCount; j++)
      {
        Rotamer* pRotamer = RotamerSetGet(pSet, j);
        if (strcmp(RotamerGetType(pRotamer), "SER") == 0)
        {
          RotamerRestore(pRotamer, pSet);
          int atomIndex;
          RotamerFindAtom(pRotamer, "HG", &atomIndex);
          RotamerCopy(&tempRot, pRotamer);
          RotamerExtract(pRotamer);
          RotamerSetAdd(pSet, &tempRot);
          for (int k = 0; k < EXPANDED_ROT_SER; k++)
          {
            ics.icParam[2] = icParaX + 2.0 * PI * (k + 1) / (EXPANDED_ROT_SER + 1);
            if (ics.icParam[2] > PI) ics.icParam[2] -= 2 * PI;
            GetFourthAtom(&RotamerGetAtomByName(&tempRot, "CA")->xyz,
              &RotamerGetAtomByName(&tempRot, "CB")->xyz,
              &RotamerGetAtomByName(&tempRot, "OG")->xyz,
              ics.icParam,
              &RotamerGetAtomByName(&tempRot, "HG")->xyz);
            XYZArraySet(&tempRot.xyzs, atomIndex, &RotamerGetAtomByName(&tempRot, "HG")->xyz);
            // check if the rotamer should be added
            BOOL expandedRotAccepted = TRUE;
            for (int kk = 0; kk < RotamerGetAtomCount(&tempRot); kk++)
            {
              Atom* pAtomK = RotamerGetAtom(&tempRot, kk);
              if (strcmp(AtomGetName(pAtomK), "HG") != 0) continue;
              for (int ss = 0; ss < RotamerGetAtomCount(&tempRot); ss++)
              {
                Atom* pAtomS = RotamerGetAtom(&tempRot, ss);
                if (strcmp(AtomGetName(pAtomK), AtomGetName(pAtomS)) == 0) continue;
                int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtomK), AtomGetName(pAtomS), RotamerGetBonds(&tempRot));
                if (bondType == 14 || bondType == 15)
                {
                  double distance = XYZDistance(&pAtomK->xyz, &pAtomS->xyz);
                  //if(distance < (pAtomK->vdw_radius+pAtomS->vdw_radius)*0.75){
                  if (distance < 2.25)
                  {
                    expandedRotAccepted = FALSE;
                    break;
                  }
                }
              }
              if (expandedRotAccepted == FALSE) break;
            }
            if (expandedRotAccepted == TRUE)
            {
              RotamerSetAdd(pSet, &tempRot);
              addedCount++;
            }
          }
        }
      }
      //printf("Design site (%2d, %4d): %d SER rots expanded\n", chnNdx, ResidueGetPosInChain(pDesignSite->pRes), addedCount);
      RotamerDestroy(&tempRot);
    }
    // for thr rots
    if (RotamerSetGetRepresentative(pSet, "THR") != NULL)
    {
      ResiTopoSetGet(pTopos, "THR", &tops);
      ResidueTopologyFindCharmmIC(&tops, "HG1", &ics);
      double icParaX = ics.icParam[2];
      int addedCount = 0;
      Rotamer tempRot;
      RotamerCreate(&tempRot);
      for (int j = 0; j < rotamerCount; j++)
      {
        Rotamer* pRotamer = RotamerSetGet(pSet, j);
        if (strcmp(RotamerGetType(pRotamer), "THR") == 0)
        {
          RotamerRestore(pRotamer, pSet);
          int atomIndex;
          RotamerFindAtom(pRotamer, "HG1", &atomIndex);
          RotamerCopy(&tempRot, pRotamer);
          RotamerExtract(pRotamer);
          for (int k = 0; k < EXPANDED_ROT_THR; k++)
          {
            BOOL expandedRotAccepted = TRUE;
            ics.icParam[2] = icParaX + 2.0 * PI * (k + 1) / (EXPANDED_ROT_THR + 1);
            if (ics.icParam[2] > PI) ics.icParam[2] -= 2 * PI;
            GetFourthAtom(&RotamerGetAtomByName(&tempRot, "CA")->xyz,
              &RotamerGetAtomByName(&tempRot, "CB")->xyz,
              &RotamerGetAtomByName(&tempRot, "OG1")->xyz,
              ics.icParam,
              &RotamerGetAtomByName(&tempRot, "HG1")->xyz);
            XYZArraySet(&tempRot.xyzs, atomIndex, &RotamerGetAtomByName(&tempRot, "HG1")->xyz);
            //RotamerSetAdd(&tempSet, &tempRot);
            for (int kk = 0; kk < RotamerGetAtomCount(&tempRot); kk++)
            {
              Atom* pAtomK = RotamerGetAtom(&tempRot, kk);
              if (strcmp(AtomGetName(pAtomK), "HG1") != 0) continue;
              for (int ss = 0; ss < RotamerGetAtomCount(&tempRot); ss++)
              {
                Atom* pAtomS = RotamerGetAtom(&tempRot, ss);
                if (strcmp(AtomGetName(pAtomK), AtomGetName(pAtomS)) == 0) continue;
                int bondType = ResidueIntraBondConnectionCheck(AtomGetName(pAtomK), AtomGetName(pAtomS), RotamerGetBonds(&tempRot));
                if (bondType == 14 || bondType == 15)
                {
                  double distance = XYZDistance(&pAtomK->xyz, &pAtomS->xyz);
                  //if(distance < (pAtomK->vdw_radius+pAtomS->vdw_radius)*0.75){
                  if (distance < 2.25)
                  {
                    expandedRotAccepted = FALSE;
                    break;
                  }
                }
              }
              if (expandedRotAccepted == FALSE) break;
            }
            if (expandedRotAccepted == TRUE)
            {
              RotamerSetAdd(pSet, &tempRot);
              addedCount++;
            }
          }

        }
      }
      //printf("Design site (%2d, %4d): %d THR rots expanded\n", chnNdx, ResidueGetPosInChain(pDesignSite->pRes), addedCount);
      RotamerDestroy(&tempRot);
    }
    // for tyr rots
    if (RotamerSetGetRepresentative(pSet, "TYR") != NULL)
    {
      ResiTopoSetGet(pTopos, "TYR", &tops);
      ResidueTopologyFindCharmmIC(&tops, "HH", &ics);
      double icPara_Tyr = ics.icParam[2];
      int addedCount = 0;
      Rotamer tempRot;
      RotamerCreate(&tempRot);
      for (int j = 0; j < rotamerCount; j++)
      {
        Rotamer* pRotamer = RotamerSetGet(pSet, j);
        if (strcmp(RotamerGetType(pRotamer), "TYR") == 0)
        {
          RotamerRestore(pRotamer, pSet);
          int atomIndex = -1;
          RotamerFindAtom(pRotamer, "HH", &atomIndex);
          RotamerCopy(&tempRot, pRotamer);
          RotamerExtract(pRotamer);
          for (int k = 0; k < EXPANDED_ROT_TYR; k++)
          {
            ics.icParam[2] = icPara_Tyr + 2.0 * PI * (k + 1) / (EXPANDED_ROT_TYR + 1);
            if (ics.icParam[2] > PI) ics.icParam[2] -= 2 * PI;
            GetFourthAtom(&RotamerGetAtomByName(&tempRot, "CE1")->xyz,
              &RotamerGetAtomByName(&tempRot, "CZ")->xyz,
              &RotamerGetAtomByName(&tempRot, "OH")->xyz,
              ics.icParam,
              &RotamerGetAtomByName(&tempRot, "HH")->xyz);
            XYZArraySet(&tempRot.xyzs, atomIndex, &RotamerGetAtomByName(&tempRot, "HH")->xyz);
            RotamerSetAdd(pSet, &tempRot);
            addedCount++;
          }
        }

      }
      //printf("Design site (%2d, %4d): %d TYR rots expanded\n", chnNdx, ResidueGetPosInChain(pDesignSite->pRes), addedCount);
      RotamerDestroy(&tempRot);
    }
    ResidueTopologyDestroy(&tops);
    CharmmICDestroy(&ics);
  }
  return Success;
}

#define NOT_FOR_DESIGN

int StructureGenerateWildtypeRotamers(Structure* pStructure, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos)
{
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL) continue;
    for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
    {
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if (pResidue->desType != Type_DesType_Fixed) continue;
      ResidueSetDesignType(pResidue, Type_DesType_Repackable);
      ProteinSiteBuildSpecifiedRotamers(pStructure, i, j, rotlib, atomParams, resiTopos, pResidue->desType);
      if (FLAG_USE_INPUT_SC && pResidue->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
      if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
    }
  }

  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL) continue;
    for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
    {
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if (pResidue->desType != Type_DesType_Fixed) continue;
      Atom* pAtomCA1 = ResidueGetAtomByName(pResidue, "CA");
      BOOL interResi = FALSE;
      for (int k = 0; k < StructureGetChainCount(pStructure); k++)
      {
        if (k == i)continue;
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) != Type_Chain_Protein)continue;
        if (strstr(DES_CHAINS, ChainGetName(pChainK)) == NULL)continue;
        for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
        {
          Residue* pResidueKS = ChainGetResidue(pChainK, s);
          Atom* pAtomCA2 = ResidueGetAtomByName(pResidueKS, "CA");
          if (XYZDistance(&pAtomCA2->xyz, &pAtomCA1->xyz) > 15.0) continue;
          if (AtomArrayCalcMinDistance(&pResidue->atoms, &pResidueKS->atoms) < CUT_PPI_DIST_SHELL1)
          {
            interResi = TRUE;
            break;
          }
        }
        if (interResi) break;
      }
      if (interResi)
      {
        ResidueSetDesignType(pResidue, Type_DesType_Repackable);
        ProteinSiteBuildSpecifiedRotamers(pStructure, i, j, rotlib, atomParams, resiTopos, pResidue->desType);
        if (FLAG_USE_INPUT_SC && pResidue->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
        if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
      }
    }
  }

  return Success;
}


int StructureGenerateWildtypeRotamersByBBdepRotLib(Structure* pStructure, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos)
{
  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChainI)) == NULL) continue;
    for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
    {
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if (pResidue->desType != Type_DesType_Fixed) continue;
      ResidueSetDesignType(pResidue, Type_DesType_Repackable);
      ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
      if (FLAG_USE_INPUT_SC && pResidue->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
      if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
    }
  }

  for (int i = 0; i < StructureGetChainCount(pStructure); i++)
  {
    Chain* pChainI = StructureGetChain(pStructure, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein) continue;
    if (strstr(DES_CHAINS, ChainGetName(pChainI)) != NULL) continue;
    for (int j = 0; j < ChainGetResidueCount(pChainI); j++)
    {
      Residue* pResidue = ChainGetResidue(pChainI, j);
      if (pResidue->desType != Type_DesType_Fixed) continue;
      Atom* pAtomCA1 = ResidueGetAtomByName(pResidue, "CA");
      BOOL interResi = FALSE;
      for (int k = 0; k < StructureGetChainCount(pStructure); k++)
      {
        if (k == i)continue;
        Chain* pChainK = StructureGetChain(pStructure, k);
        if (ChainGetType(pChainK) != Type_Chain_Protein)continue;
        if (strstr(DES_CHAINS, ChainGetName(pChainK)) == NULL)continue;
        for (int s = 0;s < ChainGetResidueCount(pChainK);s++)
        {
          Residue* pResidueKS = ChainGetResidue(pChainK, s);
          Atom* pAtomCA2 = ResidueGetAtomByName(pResidueKS, "CA");
          if (XYZDistance(&pAtomCA2->xyz, &pAtomCA1->xyz) > 15.0) continue;
          if (AtomArrayCalcMinDistance(&pResidue->atoms, &pResidueKS->atoms) < CUT_PPI_DIST_SHELL1)
          {
            interResi = TRUE;
            break;
          }
        }
        if (interResi) break;
      }
      if (interResi)
      {
        ResidueSetDesignType(pResidue, Type_DesType_Repackable);
        ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(pStructure, i, j, rotlib, atomParams, resiTopos);
        if (FLAG_USE_INPUT_SC && pResidue->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
        if (FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
      }
    }
  }

  return Success;
}

