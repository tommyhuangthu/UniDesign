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

#ifndef RESIDUE_TOPOLOGY_H
#define RESIDUE_TOPOLOGY_H

#include "Utility.h"
#include "Atom.h"

typedef enum _Type_Bond
{
  Type_Bond_Single,
  Type_Bond_double,
  Type_Bond_Triple,
  Type_Bond_None
}Type_Bond;

typedef struct _Bond
{
  char atomFromName[MAX_LEN_ATOM_NAME+1];
  char atomToName[MAX_LEN_ATOM_NAME+1];
  Type_Bond type;
} Bond;

int BondCreate(Bond* pThis);
int BondDestroy(Bond* pThis);
int BondCopy(Bond* pThis, Bond* pOther);
char* BondGetFromName(Bond* pThis);
int BondSetFromName(Bond* pThis, char* from);
char* BondGetToName(Bond* pThis);
int BondSetToName(Bond* pThis, char* to);
Type_Bond BondGetType(Bond* pThis);
int BondSetType(Bond* pThis, Type_Bond newType);
int BondShow(Bond* pThis);

typedef struct _Bonds
{
  int num;
  Bond* bonds;
} Bonds;

int BondsCreate(Bonds* pThis);
int BondsDestroy(Bonds* pThis);
int BondsCopy(Bonds* pThis, Bonds* pOther);
int BondsGetNum(Bonds* pThis);
int BondsAdd(Bonds* pThis, char* atom1, char* atom2, Type_Bond bondType);
int BondsRemove(Bonds* pThis, char* atom1, char* atom2);
Type_Bond BondsFind(Bonds* pThis, char* atom1, char* atom2);
Bond* BondsGet(Bonds* pThis, int index);
int BondsShow(Bonds* pThis);

typedef struct _CharmmIC
{
  double icParam[5];
  char atomNames[4][MAX_LEN_ATOM_NAME+1];
  BOOL torFlag;
  // Proper dihedral angle (A, B, C, D): five parameters are R(AB), Theta(ABC), Phi(ABCD), Theta(BCD), R(CD)
  // Improper dihedral angle (A, B, *C, D); five parameters are R(AC), Theta(ACB), Phi(ABCD), Theta(BCD), R(CD)
}CharmmIC;

int CharmmICCreate(CharmmIC* pThis);
int CharmmICCreateFromStringArray(CharmmIC* pThis, StringArray* params);
int CharmmICDestroy(CharmmIC* pThis);
int CharmmICCopy(CharmmIC* pThis, CharmmIC* pOther);
char* CharmmICGetAtomA(CharmmIC* pThis);
char* CharmmICGetAtomB(CharmmIC* pThis);
char* CharmmICGetAtomC(CharmmIC* pThis);
char* CharmmICGetAtomD(CharmmIC* pThis);
double* CharmmICGetICParams(CharmmIC* pThis);
int CharmmICSetTorsionAngle(CharmmIC* pThis, double newTorsionAngle);
BOOL CharmmICGetTorsionProperFlag(CharmmIC* pThis);
int CharmmICCalcXYZ(CharmmIC* pThis, AtomArray* atomArray, XYZ* pDestXYZ);
int CharmmICShow(CharmmIC* pThis);

typedef struct _ResTopology
{
  StringArray atoms;
  StringArray dels;
  Bonds bonds;
  CharmmIC* ics;
  char resName[MAX_LEN_RES_NAME+1];
  int icNum;
} ResTopology;

int ResTopologyCreate(ResTopology* pThis);
int ResTopologyCreateFromFileReader(ResTopology* pThis, FileReader* fr);
int ResTopologyDestroy(ResTopology* pThis);
int ResTopologyCopy(ResTopology* pThis, ResTopology* pOther);

char* ResTopologyGetName(ResTopology* pThis);
int ResTopologySetName(ResTopology* pThis, char* newName);

StringArray* ResTopologyGetAtoms(ResTopology* pThis);
StringArray* ResTopologyGetDeletes(ResTopology* pThis);
Bonds* ResTopologyGetBonds(ResTopology* pThis);

int ResTopologyGetCharmmICNum(ResTopology* pThis);
int ResTopologyGetCharmmIC(ResTopology* pThis, int index, CharmmIC* ic);
int ResTopologyFindCharmmIC(ResTopology* pThis, char* atomDName, CharmmIC* ic);
int ResTopologyFindCharmmICIndex(ResTopology* pThis, char* atomDName, int* ndx);
int ResTopologyAddCharmmIC(ResTopology* pThis, CharmmIC* ic);
int ResTopologyShow(ResTopology* pThis);



typedef struct _ResTopos
{
  int num;
  ResTopology* topos;
} ResTopos;

int ResToposCreate(ResTopos* pThis);
int ResToposDestroy(ResTopos* pThis);
int ResToposCopy(ResTopos* pThis, ResTopos* pOther);
int ResToposGet(ResTopos* pThis, char* resName, ResTopology* pTopo);
int ResToposFindIndexByName(ResTopos* pThis, char* resName, int* ndx);
int ResToposAdd(ResTopos* pThis, ResTopology* pTopo);
int ResToposAddFromFile(ResTopos* pThis, char* file);
int ResToposRead(ResTopos* pThis, char* file);
int ResToposShow(ResTopos* pThis);
#endif