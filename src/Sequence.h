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

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "Structure.h"
#include "EnergyMatrix.h"
#include "EnergyFunction.h"

#define  MAX_SEQ_LEN  2000

typedef struct _Sequence
{
  double        etot;
  double        ephy;
  double        ebin;
  double        eevo;
  double        ebpf;
  int           desSiteCount;
  IntArray      rotNdxs;
  int           numOfUnsatisfiedCons;
  char          fileToSaveThisSeq[MAX_LEN_FILE_NAME + 1];
} Sequence;

int SequenceCreate(Sequence* pThis);
int SequenceDestroy(Sequence* pThis);
int SequenceCopy(Sequence* pThis, Sequence* pOther);
int SequenceWriteDesignRotamer(Sequence* pThis, Structure* pStructure, int index, FILE* pFile);
int SequenceWriteDesignFasta(Sequence* pThis, Structure* pStructure, int index, FILE* pFile);

int DesignShowMinEnergyDesignStructure(Structure* pStructure, Sequence* pSeq, char* pdbfile);
int DesignShowMinEnergyDesignSites(Structure* pStructure, Sequence* pSeq, char* pdbfile);
int DesignShowMinEnergyDesignLigand(Structure* pStructure, Sequence* sequence, char* mol2file);
int StructureGetWholeSequence(Structure* pStructure, Sequence* sequence, char* seq);
int StructureWriteFirstLigandRotamerIntoMol2(Structure* pStructure, char* mol2file);

#endif // SEQUENCE_H
