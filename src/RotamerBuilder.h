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

#define NATIVE_ROTAMER_DUNBRACK 0.0//kcal

#define EXPANDED_ROT_SER  5
#define EXPANDED_ROT_THR  5
#define EXPANDED_ROT_TYR  1

//deal with BBind protein rots
int ProteinSiteBuildAllRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int ProteinSiteBuildSpecifiedRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, Type_ResidueDesignType type);
int ProteinSiteBuildMutatedRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, StringArray* pDesignTypes, StringArray* pPatchTypes);
int ProteinSiteBuildWildtypeRotamers(Structure* pThis, int chainIndex, int resiIndex, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int StructureGenerateSpecifiedProteinRotamers(Structure* pThis, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* resfile);
int StructureGenerateAllRotamers(Structure* pThis, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int StructureGenerateWildtypeRotamers(Structure* pThis, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
BOOL IsNativeRotamerInBBindRotLib(Structure* pThis, int chainIndex, int resiIndex, ResiTopoSet* pResiTopos, BBindRotamerLib* pBBindRotLib, double torsionStd);

//deal with BBdep protein rots
int ProteinSiteBuildAllRotamersByBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int ProteinSiteBuildSpecifiedRotamersByBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int ProteinSiteBuildMutatedRotamersByBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, StringArray* pDesignTypes, StringArray* pPatchTypes);
int ProteinSiteBuildWildtypeRotamersByBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int StructureBuildAllRotamersByBBdepRotLib(Structure* pThis, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int StructureBuildResfileRotamersByBBdepRotLib(Structure* pThis, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* resfile);
int StructureGenerateWildtypeRotamersByBBdepRotLib(Structure* pThis, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);
int StructureBuildCatalyticRotamersByBBdepRotLib(Structure* pThis, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* resfile);
int StructureBuildPLIShell1RotamersByBBdepRotLib(Structure* pThis, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* resfile);
int StructureBuildPLIShell2RotamersByBBdepRotLib(Structure* pThis, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* resfile);
int StructureBuildPPIRotamersByBBdepRotLib(Structure* pStructure, BBdepRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos);

BOOL IsNativeRotamerInBBdepRotLib(Structure* pThis, int chainIndex, int resiIndex, ResiTopoSet* pResiTopos, BBdepRotamerLib* pBBdepRotLib, double torsionStd);

//deal with small molecule
int StructureGetTruncatedBackbone(Structure* pThis, Residue* pSmallMol, double activeSiteRange, BOOL withHydrogen, AtomArray* pBackboneAtoms);
int StructureGetTruncatedBackboneNew(Structure* pThis, Residue* pSmallMol, BOOL withHydrogen, AtomArray* pBackboneAtoms);
int StructureDeployCataConsSitePair(Structure* pThis, CataConsSitePair* pCataConsSitePair);
int StructurePlaceSmallMol(Structure* pThis, PlacingRule* pPlacingRule, CataConsSitePairArray* pCataConsCollection, int relatedProteinSiteCount, DesignSite** relatedProteinSites, RotamerSet* pSmallMolRotSet);
int StructureGenerateSmallMolRotamers(Structure* pThis, char* cataConsFile, char* placingRuleFile);
int StructureReadSmallMolRotamers(Structure* pThis, ResiTopoSet* resiTopos, char* smallMolFile);
int StructureWriteSmallMolRotamers(Structure* pThis, char* smallMolFile);
int StructureSmallmolOrientationScreen(Structure* pStructure, ResiTopoSet* pResiTopo, char* initialPoseFile, char* newPoseFile, char* screenRuleFile);

// general function for BBdep and BBind protein rots, and ligand rots
int ProteinSiteWriteRotamers(Structure* pStructure, int chainIndex, int resiIndex, char* rotamerFilePath);
int ProteinSiteBuildNativeRotamer(Structure* pThis, int chainIndex, int resiIndex, ResiTopoSet* pResiTopos);
int ProteinSiteBuildFlippedNativeRotamer(Structure* pStructure, int chainIndex, int resiIndex, ResiTopoSet* pResiTopos);
int ProteinSiteExpandHydroxylRotamers(Structure* pStructure, int chainIndex, int resiIndex, ResiTopoSet* pTopos);


#endif