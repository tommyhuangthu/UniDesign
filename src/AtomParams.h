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

#ifndef ATOMPARAMS_H
#define ATOMPARAMS_H

#include "Atom.h"

typedef struct _AtomParams
{
  StringArray resNames;
  Atom** atoms;
  int* atomNum;
} AtomParams;

int AtomParamsCreate(AtomParams* pThis);
int AtomParamsDestroy(AtomParams* pThis);
int AtomParamsAdd(AtomParams* pThis, char* resName, Atom* pAtom);
int AtomParamsAddFromFile(AtomParams* pThis, char* file);
int AtomParamsGetResidueNum(AtomParams* pThis);
int AtomParamsGetResidueName(AtomParams* pThis, int index, char* resName);
int AtomParamsGetAtomNumOfResidue(AtomParams* pThis, char* resName, int* num);
int AtomParamsGetAtomParamOfResidue(AtomParams* pThis, char* resName, int ndx, Atom* pAtom);
int AtomParamsGetAtomParamByResidueName(AtomParams* pThis, char* resName, char* atomName, Atom* pAtom);
int AtomParamsRead(AtomParams* pThis, char* file);

#endif