/*******************************************************************************************************************************
Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

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

#ifndef SMALLMOL_PAR_AND_TOPO_H
#define SMALLMOL_PAR_AND_TOPO_H

#include "Utility.h"
#include "Atom.h"
#include "Residue.h"
#include "SmallMolEEF1.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


int FindAtomC(Atom *pAtomD,AtomArray *pAtomArray,IntArray *icFlags,IntArray *bondsFromToAndType,AtomArray *pAtomCArray);
int FindAtomB(Atom *pAtomC,Atom* pAtomD,AtomArray *pAtoms,IntArray *icFlags,IntArray *bondsFromToAndType,AtomArray *pAtomBArray);
int GenerateSmallMolTopologyFromMol2(char* mol2file,char* topfile,char* firstA,char *firstB,char* firstC);
int GetAttachedAtomTypeNum(AtomArray* pAtoms,int bonds[][2],int bondNum,int atomIndex,int attachAtomTypeNum[10],int *attachCount,int attachIndex[10]);
int GenerateSmallMolParameterFromMol2(char* mol2file,char* parfile);
int GenerateParameterAndTopologyFromMol2(char* mol2file,char* parfile,char* topfile);


#endif // SMALLMOL_PAR_AND_TOPO_H
