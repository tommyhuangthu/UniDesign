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

#ifndef PROGRAM_FUNCTION_H
#define PROGRAM_FUNCTION_H

#include "Structure.h"
#include "EnergyFunction.h"


int PrintHelp();
int PrintVersion();
int PrintAdvertisement();

int ComputeStructureStability(Structure* pStructure, AAppTable* pAAppTable, RamaTable* pRama, double energyTerms[MAX_ENERGY_TERM]);
int ComputeStructureStabilityByBBdepRotLib(Structure* pStructure, AAppTable* pAAppTable, RamaTable* pRama, BBdepRotamerLib* pRotLib, double energyTerms[MAX_ENERGY_TERM]);
int ComputeStructureStabilityByBBdepRotLib2(Structure* pStructure, AAppTable* pAAppTable, RamaTable* pRama, char* dunlibfile, double energyTerms[MAX_ENERGY_TERM]);

int ComputeBinding(Structure* pStructure);
int ComputeBindingWithChainSplitting(Structure* pStructure, char split1[], char split2[]);

int BuildMutant(Structure* pStructure, char* mutantfile, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* pdbid);
int BuildMutantByBBdepRotLib(Structure* pStructure, char* mutantfile, BBdepRotamerLib* pBBdepRotLib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* pdbid);

int RepairStructure(Structure* pStructure, BBindRotamerLib* rotlib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* pdbid);
int RepairStructureByBBdepRotLib(Structure* pStructure, BBdepRotamerLib* pBBdepRotLib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* pdbid);
int EnergyMinimizationByBBdepRotLib(Structure* pStructure, BBdepRotamerLib* pBBdepRotLib, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* pdbid);

int AddPolarHydrogen(Structure* pStructure, char* pdbid);
int OptimizeHydrogen(Structure* pStructure, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* pdbid);

int ComputeResidueInteractionWithFixedEnvironment(Structure* pStructure, int chainIndex, int residueIndex);
int ComputeAllRotamersEnergy(Structure* pStructure, BBindRotamerLib* rotlib, AAppTable* pAAppTable, RamaTable* pRama, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* pdbid);
int ComputeWildtypeRotamersEnergy(Structure* pStructure, BBindRotamerLib* rotlib, AAppTable* pAAppTable, RamaTable* pRama, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* pdbid);
int ComputeRotamersEnergyByBBdepRotLib(Structure* pStructure, BBdepRotamerLib* rotlib, AAppTable* pAAppTable, RamaTable* pRama, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* pdbid);
int ComputeWildtypeRotamersEnergyByBBdepRotLib(Structure* pStructure, BBdepRotamerLib* rotlib, AAppTable* pAAppTable, RamaTable* pRama, AtomParamsSet* atomParams, ResiTopoSet* resiTopos, char* pdbid);

int ShowPhiPsi(Structure* pStructure, char* phipsifile);

int CheckRotamerInBBindRotLib(Structure* pStructure, BBindRotamerLib* pRotLib, ResiTopoSet* pTopos, double cutoff, char* pdbid);
int CheckRotamerInBBdepRotLib(Structure* pStructure, BBdepRotamerLib* pRotLib, ResiTopoSet* pTopos, double cutoff, char* pdbid);
int FindMinRmsdRotFromRotLib(Structure* pStructure, char* pdbid);
int CompareSidechainsOf2Structures(Structure* pStructure, Structure* pStructure2, FILE* pTorsion, FILE* pRmsd);

int CheckClash0(Structure* pStructure, double clashRatio);
int CheckClash1(Structure* pStructure, double clashRatio);
int CheckClash2(Structure* pStructure, double clashRatio);

int FindInterfaceResidues(Structure* pStructure);
int FindInterfaceResiduesWithChainSplitting(Structure* pStructure, char split1[], char split2[]);
int FindCoreResidues(Structure* pStructure);
int FindSurfaceResidues(Structure* pStructure);
int FindIntermediateResidues(Structure* pStructure);

int StructureGetAminoAcidComposition(Structure* pStructure, int* aas);

int SelectResiduesWithin(Structure* pStructure);


#endif //PROGRAM_FUNCTION_H