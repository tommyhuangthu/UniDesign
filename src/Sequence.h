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

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "Structure.h"
#include "EnergyFunction.h"

#define  MAX_LEN_SEQ  3000

typedef struct _Sequence
{
  double        etot;
  double        ephy;
  double        ebin;
  double        eevo;
  double        ebpf;
  int           desSiteNum;
  IntArray      rotNdxs;
  int           unmetConsNum;
  int           mutNum;
  BOOL          seqExist;
  char          pdbToSaveSeq[MAX_LEN_FILE_NAME+1];
} Sequence;

int SequenceCreate(Sequence* pThis);
int SequenceDestroy(Sequence* pThis);
int SequenceCopy(Sequence* pThis, Sequence* pOther);
int SequenceWriteRotamerIndex(Sequence* pThis, Structure* pStruct, int index, FILE* pFile);
int SequenceWriteDesignSummaryToFile(Sequence* pThis, Structure* pStruct, int ndx, FILE* pFile);
int SequenceWriteDesignSequenceToFile(Sequence* pSeq, Structure* pStruct, int ndx, FILE* pFile);
int SequenceWriteStructureInPDB(Sequence* pSeq, Structure* pStruct, char* file);
int SequenceWriteAllSitesInPDB(Sequence* pSeq, Structure* pStruct, char* file);
int SequenceWriteMutableSitesInPDB(Sequence* pSeq, Structure* pStruct, char* file);
int StructureWriteFirstLigandRotamerIntoMol2(Structure* pStruct, char* file);
int SequenceWriteLigandInMol2(Sequence* pSeq, Structure* pStruct, char* file);
int StructureGetWholeSequence(Sequence* pSeq, Structure* pStruct, char* seqstr);

#endif // SEQUENCE_H
