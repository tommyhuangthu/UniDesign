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

#ifndef ROTAMER_H
#define ROTAMER_H

#include "Atom.h"
#include "Residue.h"

typedef struct _Rotamer
{
  AtomArray atoms;
  Bonds bonds;
  XYZArray xyzs;
  char type[MAX_LEN_RES_NAME+1];
  char chainName[MAX_LEN_CHAIN_NAME+1];
  int  posInChain;
  double vdwBackbone;
  double vdwInternal;
  double selfEnergy;
  double selfEnergyBin;
  double dunbrack;
  DoubleArray Xs; //side-chain torsions
} Rotamer;

int RotamerCreate(Rotamer* pThis);
int RotamerDestroy(Rotamer* pThis);
int RotamerCopy(Rotamer* pThis, Rotamer* pOther);
char* RotamerGetType(Rotamer* pThis);
int RotamerSetType(Rotamer* pThis, char* newType);
int RotamerCopyAtomXYZ(Rotamer* pThis, XYZArray* pNewXYZ);
int RotamerAtomsCopyXYZFromArray(Rotamer* pThis, XYZArray* pNewXYZ);
char* RotamerGetChainName(Rotamer* pThis);
int RotamerSetChainName(Rotamer* pThis, char* newChainname);
int RotamerGetPosInChain(Rotamer* pThis);
int RotamerSetPosInChain(Rotamer* pThis, int newPosInChain);
double RotamerGetDunbrack(Rotamer* pThis);
int RotamerSetDunbrack(Rotamer* pThis, double dun);
int RotamerShow(Rotamer* pThis);
BOOL RotamerAndRotamerInSameType(Rotamer* pThis, Rotamer* pOther);
BOOL RotamerAndResidueInSameType(Rotamer* pThis, Residue* pOther);

// These functions below are valid only when Rotamer is restored
int RotamerGetAtomNum(Rotamer* pThis);
Atom* RotamerGetAtom(Rotamer* pThis, int index);
Atom* RotamerGetAtomByName(Rotamer* pThis, char* atomName);
int RotamerFindAtom(Rotamer* pThis, char* atomName, int* index);
int RotamerAddAtoms(Rotamer* pThis, AtomArray* pNewAtoms);
Bonds* RotamerGetBonds(Rotamer* pThis);

typedef enum _Type_ProteinAtomOrder
{
  Type_ProteinAtomOrder_Alpha = 0,
  Type_ProteinAtomOrder_Beta,
  Type_ProteinAtomOrder_Gamme,
  Type_ProteinAtomOrder_Delta,
  Type_ProteinAtomOrder_Epsilon,
  Type_ProteinAtomOrder_Zeta,
  Type_ProteinAtomOrder_Other
} Type_ProteinAtomOrder;

int Type_ProteinAtomOrder_ToInt(Type_ProteinAtomOrder order);
Type_ProteinAtomOrder Type_ProteinAtomOrder_FromInt(int order);
Type_ProteinAtomOrder Type_ProteinAtomOrder_JudgedByAtomName(char* atomName);

int RotamerOfProteinInitAtomsAndBonds_Charmm19(Rotamer* pThis, Residue* pResi, AtomParams* atomParams, ResTopos* resiTopos);
int RotamerOfProteinInitAtomsAndBonds_Charmm22(Rotamer* pThis, Residue* pResi, AtomParams* atomParams, ResTopos* resiTopos);
int RotamerOfProteinPatch(Rotamer* pThis, char* patchType, AtomParams* atomParams, ResTopos* resiTopo);
int RotamerOfProteinCalcXYZ(Rotamer* pThis, Residue* pResi, char* patchName, DoubleArray* torsions, ResTopos* resiTopos);
int RotamerOfProteinGenerate(Rotamer* pThis, Residue* pResi, char* rotamerType, char* patchType, DoubleArray* torsions, AtomParams* atomParams, ResTopos* resiTopos);
int RotamerShowInPDBFormat(Rotamer* pThis, char* header, char* chainName, int atomIndex, int resiIndex, FILE* pFile);
int RotamerShowAtomParameter(Rotamer* pThis);


double RotamerAndRotamerSidechainRMSD(Rotamer* pThis, Rotamer* pOther);
double RotamerAndResidueSidechainRMSD(Rotamer* pThis, Residue* pOther);



typedef struct _Rotamers
{
  Rotamer* rots;
  Rotamer* reps; // representative rotamers
  int num;
  int capacity;
  int repNum;
} Rotamers;

int RotamersCreate(Rotamers* pThis);
int RotamersDestroy(Rotamers* pThis);
int RotamersCopy(Rotamers* pThis, Rotamers* pOther);
int RotamersGetNum(Rotamers* pThis);
Rotamer* RotamersGet(Rotamers* pThis, int index);
Rotamer* RotamersGetRepByIndex(Rotamers* pThis, int index);
Rotamer* RotamersFindRep(Rotamers* pThis, char* type);
int RotamersAdd(Rotamers* pThis, Rotamer* pNewRotamer);
int RotamersShow(Rotamers* pThis, FILE* pFile);

int RotamerExtract(Rotamer* pThis);
int RotamerRestore(Rotamer* pThis, Rotamers* pRots);
int RotamersGetRepNum(Rotamers* pThis);
int RotamerShowBondInformation(Rotamer* pThis);


typedef struct _RotLibPhiPsi
{
  StringArray rotTypes;
  IntArray    rotNums;
  DoubleArray** tors;
  DoubleArray** devs;
  DoubleArray* prob;
  double phi;
  double psi;
}RotLibPhiPsi;

int RotLibPhiPsiCreate(RotLibPhiPsi* pThis);
int RotLibPhiPsiDestroy(RotLibPhiPsi* pThis);
int RotLibPhiPsiGetNum(RotLibPhiPsi* pThis, char* typeName);
int RotLibPhiPsiGet(RotLibPhiPsi* pThis, char* typeName, int index, DoubleArray* pDestTorsion, double* probability);


typedef struct _BBdepRotLib
{
  int num; // number of Phi-
  RotLibPhiPsi* pps; // Phi-Psi's
}BBdepRotLib;


int BBdepRotamerLibCreate(BBdepRotLib* pRotLib, char* rotlibfile);
int BBdepRotamerLibCreate2(BBdepRotLib* pRotLib, char* binlibfile);
int BBdepRotamerLibDestroy(BBdepRotLib* pThis);
int RotamerOfProteinGenerateByBBdepRot(Rotamer* pThis, Residue* pResi, char* rotamerType, char* patchType, DoubleArray* torsions, AtomParams* atomParams, ResTopos* resiTopos);
int RotamerSetOfProteinGenerateByBBdepRotLib(Rotamers* pThis, Residue* pResi, StringArray* designTypes, StringArray* patchTypes, BBdepRotLib* bbrotlib, AtomParams* atomParams, ResTopos* resiTopo);
int RotamerCalcDunbrackEnergy(Rotamer* pThis, double probability);


BOOL RotamerAndResidueWithSimilarTorsions(Rotamer* pThis, Residue* pOther, double cutoff);
int ResidueCalcSidechainTorsion(Residue* pThis, ResTopos* pTopos);
Rotamer* RotamerSetFindFirstGivenTypeRotamer(Rotamers* pThis, char* type);

BOOL RotamerIsSymmetricalCheck(Rotamer* pRotamer);
int SymmetricalRotamerGenerate(Rotamer* pSymmetrical, Rotamer* pRotamer);

int BBdepRotamerLibFromText2Binary(char* textlibfile, char* binlibfile);
int BBdepRotamerLibFromBinary2Text(char* binlibfile, char* textlibfile);

#endif //ROTAMER_H
