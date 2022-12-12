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

#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "Chain.h"
#include "Rotamer.h"
#include "DesignSite.h"
#include "EnergyFunction.h"
#include "SmallMol.h"


typedef struct _Structure
{
  char name[MAX_LEN_STRUCTURE_NAME + 1];
  Chain* chains;
  DesignSite* designSites;
  int chainNum;
  int desSiteCount;
} Structure;

int StructureCreate(Structure* pThis);
int StructureDestroy(Structure* pThis);
char* StructureGetName(Structure* pThis);
int StructureSetName(Structure* pThis, char* newName);
int StructureGetChainCount(Structure* pThis);
Chain* StructureGetChain(Structure* pThis, int index);
Chain* StructureFindChainByName(Structure* pThis, char* chainName);
int StructureFindChainIndex(Structure* pThis, char* chainName, int* index);
int StructureFindSmallMol(Structure* pThis, Residue** ppSmallMol);
int StructureAddChain(Structure* pThis, Chain* newChain);
int StructureDeleteChain(Structure* pThis, char* chainName);
int StructureShowInPDBFormat(Structure* pThis, FILE* pFile);

// deal with design sites
int StructureGetDesignSiteCount(Structure* pStructure);
DesignSite* StructureGetDesignSite(Structure* pStructure, int index);
int StructureFindDesignSiteIndexByChainNameAndPosInChain(Structure* pStructure, char* chainName, int posInChain);
int StructureGetDesignSiteIndex(Structure* pStructure, int chainIndex, int resiIndex);
DesignSite* StructureFindDesignSite(Structure* pStructure, int chainIndex, int resiIndex);
DesignSite* StructureFindDesignSiteByChainName(Structure* pStructure, char* chainName, int posInChain);

int StructureShowDesignSites(Structure* pThis);

int ProteinSiteAddDesignSite(Structure* pThis, int chainIndex, int resiIndex);
int ProteinSiteRemoveDesignSite(Structure* pThis, int chainIndex, int resiIndex);
int StructureRemoveAllDesignSites(Structure* pThis);

// other functions
int ChainComputeResiduePosition(Structure* pStructure, int chainIndex);
int StructureComputeResiduePosition(Structure* pStructure);

// debuggers
int StructureShowAtomParameter(Structure* pStructure);
int StructureShowBondInformation(Structure* pStructure);

int StructureReadPDB(Structure* pStructure, char* pdbFile, AtomParamsSet* pAtomParams, ResiTopoSet* pResiTopos);
int StructureReadMol2(Structure* pStructure, char* mol2file, AtomParamsSet* pAtomParams, ResiTopoSet* pTopos);

int StructureCalcProteinResidueSidechainTorsion(Structure* pThis, ResiTopoSet* pResiTopos);

int StructureCalcPhiPsi(Structure* pStructure);
int StructureCopy(Structure* pThis, Structure* pOther);
Residue* StructureFindNucleotidePair(Structure* pStructure, Residue* pQuery);
int StructureMutateNucleotide(Residue* pQuery, char oneLetterName, Structure* pStructure, AtomParamsSet* pAtomParaSet, ResiTopoSet* pResiTopos);
int StructureCalcAminoAcidPropensityAndRamaEnergy(Structure* pStructure, AAppTable* pAAPP, RamaTable* pRama);
int StructureCalcAminoAcidDunbrackEnergy(Structure* pStructure, BBdepRotamerLib* pBBdep);

//deal with nucleic acid

#endif // STRUCTURE_H
