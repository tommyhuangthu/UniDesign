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

#ifndef ROTAMER_BUILDER_H
#define ROTAMER_BUILDER_H

#include "Structure.h"

#define NATROT_DUNBRACK_ENERGY  0.0//kcal
#define EXPANDED_ROT_SER        2
#define EXPANDED_ROT_THR        2
#define EXPANDED_ROT_TYR        1

//deal with BBdep protein rots
int ProteinSiteBuildAllRotamersByBBdepRotLib(Structure* pStruct, int chnNdx, int resNdx, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos);
int ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(Structure* pStruct, int chnNdx, int resNdx, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos);
int ProteinSiteBuildMutatedRotamersByBBdepRotLib(Structure* pStruct, int chnNdx, int resNdx, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos, StringArray* pDesTypes, StringArray* pPatchTypes);
int ProteinSiteBuildWildtypeRotamersByBBdepRotLib(Structure* pStruct, int chnNdx, int resNdx, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos);
int StructureBuildAllRotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos);
int StructureBuildResfileRotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos, char* resfile);
int StructureGenerateWildtypeRotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos);
int StructureBuildCatalyticRotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos, char* resfile);
int StructureBuildPLIShell1RotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos, char* resfile);
int StructureBuildPLIShell2RotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos, char* resfile);
int StructureBuildPPIRotamersByBBdepRotLib(Structure* pStruct, BBdepRotLib* rotlib, AtomParams* atomParams, ResTopos* resiTopos);
BOOL IsNativeRotamerInBBdepRotLib(Structure* pStruct, int chnNdx, int resNdx, ResTopos* pResiTopos, BBdepRotLib* pBBdepRotLib, double cut);

//deal with small molecule
int StructureGetTruncatedBackbone(Structure* pStruct, Residue* pSmallMol, double range, BOOL withH, AtomArray* pBackbone);
int StructureGetTruncatedBackbone2(Structure* pStruct, Residue* pSmallMol, BOOL withH, AtomArray* pBackbone);
int StructureDeployCataConsSitePair(Structure* pStruct, CataConsSitePair* pCataConsSitePair);
int StructurePlaceSmallMol(Structure* pStruct, PlacingRule* pRule, CataConsSitePairArray* pCataArray, int siteNum, DesignSite** ppSites, Rotamers* pSmallSet);
int StructureGenerateSmallMolRotamers(Structure* pStruct, char* cataConsFile, char* ligplaceFile);
int StructureReadSmallMolRotamers(Structure* pStruct, ResTopos* resiTopos, char* smallmolFile);
int StructureWriteSmallMolRotamers(Structure* pStruct, char* smallmolFile);
int StructureSmallmolOrientationScreen(Structure* pStruct, ResTopos* resTopos, char* oldposeFile, char* newposeFile, char* orientRuleFile);

int ProteinSiteWriteRotamers(Structure* pStruct, int chnNdx, int resNdx, char* file);
int ProteinSiteBuildNativeRotamer(Structure* pStruct, int chnNdx, int resNdx, ResTopos* resiTopos);
int ProteinSiteBuildFlippedNativeRotamer(Structure* pStruct, int chnNdx, int resNdx, ResTopos* resiTopos);
int ProteinSiteExpandHydroxylRotamers(Structure* pStruct, int chnNdx, int resNdx, ResTopos* resiTopos);


#define  SELF_ENERGY_THRESHOLD_SMALLMOL           30
#define  SELF_ENERGY_THRESHOLD_SAME_AA_RES_TYPE   15
#define  SELF_ENERGY_THRESHOLD_DIFF_AA_RES_TYPE   30

typedef struct _RotamerList
{
  int desSiteNum;
  int* rotNum;
  int* remainRotNum;
  BOOL** remainFlag;
} RotamerList;

int RotamerListCreateFromStructure(RotamerList* pList, Structure* pStruct);
int RotamerListCreateFromFile(RotamerList* pList, char* file);
void RotamerListDestroy(RotamerList* pList);
int RotamerListCopy(RotamerList* pThis, RotamerList* pOther);
int RotamerListRead(RotamerList* pList, char* file);
int RotamerListWrite(RotamerList* pList, char* file);
int RotamerListShow(RotamerList* pList);

int RotamerOriginalIndexGet(RotamerList* pList, int siteNdx, int rotNdx, int* trueRotNdx);
int RotamerReducedIndexGet(RotamerList* pList, int siteNdx, int trueRotNdx, int* reducedRotNdx);

int SelfEnergyGenerate(Structure* pStruct, char* file);
int SelfEnergyGenerate2(Structure* pStruct, AAppTable* pAAppTable, RamaTable* pRamaTable, char* file);
int SelfEnergyReadAndCheck(Structure* pStruct, RotamerList* pList, char* file);

int StructureShowDesignSites2(Structure* pStruct, RotamerList* pList);



#endif