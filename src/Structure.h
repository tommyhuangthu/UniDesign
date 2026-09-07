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

#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "Chain.h"
#include "Rotamer.h"
#include "DesignSite.h"
#include "EnergyFunction.h"
#include "SmallMol.h"


typedef struct _Structure
{
  char name[MAX_LEN_STRUCTURE_NAME+1];
  Chain* chains;
  DesignSite* desSites;
  int chnNum;
  int desSiteNum;
} Structure;

int StructureCreate(Structure* pStruct);
int StructureDestroy(Structure* pStruct);
char* StructureGetName(Structure* pStruct);
int StructureSetName(Structure* pStruct, char* newName);
int StructureGetChainNum(Structure* pStruct);
Chain* StructureGetChain(Structure* pStruct, int index);
Chain* StructureFindChainByName(Structure* pStruct, char* chainName);
int StructureFindChainIndex(Structure* pStruct, char* chainName, int* index);
int StructureFindSmallMol(Structure* pStruct, Residue** ppSmallMol);
int StructureAddChain(Structure* pStruct, Chain* newChain);
int StructureDeleteChain(Structure* pStruct, char* chainName);
int StructureShowInPDBFormat(Structure* pStruct, FILE* pFile);

int StructureGetDesignSiteNum(Structure* pStruct);
DesignSite* StructureGetDesignSite(Structure* pStruct, int index);
int StructureFindDesignSiteIndexByChainNameAndPosInChain(Structure* pStruct, char* chainName, int posInChain);
int StructureGetDesignSiteIndex(Structure* pStruct, int chainIndex, int resiIndex);
DesignSite* StructureFindDesignSite(Structure* pStruct, int chainIndex, int resiIndex);
DesignSite* StructureFindDesignSiteByChainName(Structure* pStruct, char* chainName, int posInChain);

int StructureShowDesignSites(Structure* pStruct);

int ProteinSiteAddDesignSite(Structure* pStruct, int chainIndex, int resiIndex);
int ProteinSiteRemoveDesignSite(Structure* pStruct, int chainIndex, int resiIndex);
int StructureRemoveAllDesignSites(Structure* pStruct);

int ChainComputeResiduePosition(Structure* pStruct, int chainIndex);
int StructureComputeResiduePosition(Structure* pStruct);

int StructureShowAtomParameter(Structure* pStruct);
int StructureShowBondInformation(Structure* pStruct);

int StructureReadPDB(Structure* pStruct, char* pdbfile, AtomParams* pParams, ResTopos* pTops);
int StructureReadMol2(Structure* pStruct, char* mol2file, AtomParams* pParams, ResTopos* pTops);

int StructureCalcProteinResidueSidechainTorsion(Structure* pStruct, ResTopos* pTops);

int StructureCalcPhiPsi(Structure* pStruct);
int StructureCopy(Structure* pStruct1, Structure* pStruct2);
Residue* StructureFindNucleotidePair(Structure* pStruct, Residue* pRes);
int StructureMutateNucleotide(Residue* pRes, char oneLetterName, Structure* pStruct, AtomParams* pParams, ResTopos* pTops);
int StructureCalcAminoAcidPropensityAndRamaEnergy(Structure* pStruct, AAppTable* pAAPP, RamaTable* pRama);
int StructureCalcAminoAcidDunbrackEnergy(Structure* pStruct, BBdepRotLib* pLib);


#endif // STRUCTURE_H
