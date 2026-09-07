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

#ifndef PROGRAM_FUNCTION_H
#define PROGRAM_FUNCTION_H

#include "Structure.h"
#include "EnergyFunction.h"

int PrintProgramHelp();
int PrintProgramVersion();
int PrintProgramInfo();
int ComputeStability(Structure * pStruct, AAppTable * pAAppTable, RamaTable * pRama, BBdepRotLib * pLib, double energyTerms[MAX_ENERGY_TERM]);
int ComputeBinding(Structure* pStruct);
int ComputeBinding2(Structure* pStruct, char split1[], char split2[]);
int BuildMutant(Structure * pStruct, char * mutfile, BBdepRotLib * pLib, AtomParams * pParams, ResTopos * pTops);
int RepairStructure(Structure* pStruct, BBdepRotLib* pLib, AtomParams* pParams, ResTopos* pTops);
int EnergyMinimization(Structure* pStruct, BBdepRotLib* pLib, AtomParams* pParams, ResTopos* pTops);
int AddPolarHydrogen(Structure* pStruct);
int OptimizeHydrogen(Structure* pStruct, AtomParams* pParams, ResTopos* pTops);
int ComputeResidueInteractionWithFixedEnvironment(Structure * pStruct, AAppTable * pAAPP, RamaTable * pRama, BBdepRotLib * pLib, int chnNdx, int resNdx);

int ComputeRotamersEnergy(Structure* pStruct, BBdepRotLib* pLib, AAppTable* pAAPP, RamaTable* pRama, AtomParams* pParams, ResTopos* pTops);
int ComputeWildtypeRotamersEnergy(Structure* pStruct, BBdepRotLib* pLib, AAppTable* pAAPP, RamaTable* pRama, AtomParams* pParams, ResTopos* pTops);

int WritePhiPsi(Structure* pStruct, char* file);

int CheckRotamerInBBdepRotLib(Structure* pStruct, BBdepRotLib* pLib, ResTopos* pTops, double cutoff);
int FindMinRmsdRotFromRotLib(Structure* pStruct);
int CompareTwoStructureProteinSidechains(Structure* pStruct, Structure* pStruct2, FILE* pTorsion, FILE* pRmsd);

int CheckClash0(Structure* pStruct, double clashRatio);
int CheckClash1(Structure* pStruct, double clashRatio);
int CheckClash2(Structure* pStruct, double clashRatio);

int FindInterfaceResidues(Structure* pStruct);
int FindInterfaceResidues2(Structure* pStruct, char split1[], char split2[]);
int FindCoreResidues(Structure* pStruct);
int FindSurfaceResidues(Structure* pStruct);
int FindIntermediateResidues(Structure* pStruct);

int ComputeAminoAcidComposition(Structure* pStruct, int* aas);

int SelectResiduesInRange(Structure* pStruct);

int StructureFindMutations(Structure* pStruct, Structure* pStruct2, char* mutStr);


#endif //PROGRAM_FUNCTION_H

