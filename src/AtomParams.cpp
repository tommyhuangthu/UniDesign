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

#include "AtomParams.h"
#include <string.h>

int AtomParamsCreate(AtomParams* pThis)
{
  StringArrayCreate(&pThis->resNames);
  pThis->atomNum = NULL;
  pThis->atoms = NULL;
  return Success;
}


int AtomParamsDestroy(AtomParams* pThis)
{
  for (int i = 0;i < AtomParamsGetResidueNum(pThis);i++)
  {
    for (int j = 0;j < pThis->atomNum[i];j++)
    {
      AtomDestroy(&pThis->atoms[i][j]);
    }
    free(pThis->atoms[i]);
  }
  free(pThis->atoms);
  free(pThis->atomNum);
  StringArrayDestroy(&pThis->resNames);
  return Success;
}


int AtomParamsAdd(AtomParams* pThis, char* resName, Atom* pAtom)
{
  if (strlen(resName) > MAX_LEN_RES_NAME)
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    int result = ValueError;
    sprintf(errMsg, "in file %s line %d, residue name is too long", __FILE__, __LINE__);
    TraceError(errMsg, result);
    return result;
  }

  int resNdx;
  int atomNdx;
  int result = StringArrayFind(&pThis->resNames, resName, &resNdx);

  if (FAILED(result))
  {
    int newResNum = AtomParamsGetResidueNum(pThis) + 1;
    pThis->atoms = (Atom**)realloc(pThis->atoms, sizeof(Atom*) * newResNum);
    pThis->atoms[newResNum - 1] = NULL;
    StringArrayAppend(&pThis->resNames, resName);
    pThis->atomNum = (int*)realloc(pThis->atomNum, sizeof(int) * newResNum);
    pThis->atomNum[newResNum - 1] = 0;
    resNdx = AtomParamsGetResidueNum(pThis) - 1;
  }

  for (atomNdx = 0;atomNdx < pThis->atomNum[resNdx];atomNdx++)
  {
    if (strcmp(pThis->atoms[resNdx][atomNdx].name, pAtom->name) == 0)
    {
      break;
    }
  }
  if (atomNdx == pThis->atomNum[resNdx])
  {
    int newAtomNum = pThis->atomNum[resNdx] + 1;
    pThis->atoms[resNdx] = (Atom*)realloc(pThis->atoms[resNdx], sizeof(Atom) * newAtomNum);
    pThis->atomNum[resNdx]++;
  }

  AtomCopy(&pThis->atoms[resNdx][atomNdx], pAtom);

  return Success;
}


int AtomParamsAddFromFile(AtomParams* pThis, char* filePath)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  FileReader file;
  int result = FileReaderCreate(&file, filePath);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, failed to read parameter file %s", __FILE__, __LINE__, filePath);
    TraceError(errMsg, result);
    return result;
  }

  char line[MAX_LEN_ONE_LINE_CONTENT+1];
  while (!FAILED(FileReaderGetNextLine(&file, line)))
  {
    Atom atom;
    StringArray params;
    AtomCreate(&atom);
    StringArrayCreate(&params);
    StringArraySplitString(&params, line, ' ');
    char* resiName = StringArrayGet(&params, 0);
    if (resiName == NULL || FAILED(AtomSetParamsByStringArray(&atom, &params)) || FAILED(AtomParamsAdd(pThis, resiName, &atom)))
    {
      sprintf(errMsg, "in file %s line %d, wrong format in file %s line \n%s", __FILE__, __LINE__, filePath, line);
      result = FormatError;
      TraceError(errMsg, result);
      return result;
    }
    StringArrayDestroy(&params);
    AtomDestroy(&atom);
  }

  FileReaderDestroy(&file);
  return Success;
}


int AtomParamsGetResidueNum(AtomParams* pThis)
{
  return StringArrayGetNum(&pThis->resNames);
}


int AtomParamsGetResidueName(AtomParams* pThis, int index, char* residueName)
{
  char* name = StringArrayGet(&pThis->resNames, index);
  if (!name) return DataNotExistError;
  else
  {
    strcpy(residueName, name);
    return Success;
  }
}


int AtomParamsGetAtomNumOfResidue(AtomParams* pThis, char* residueName, int* num)
{
  int pos;
  int result = StringArrayFind(&pThis->resNames, residueName, &pos);
  if (FAILED(result)) return result;
  *num = pThis->atomNum[pos];
  return Success;
}


int AtomParamsGetAtomParamOfResidue(AtomParams* pThis, char* residueName, int index, Atom* pDestAtom)
{
  int pos;
  int result = StringArrayFind(&pThis->resNames, residueName, &pos);
  if (FAILED(result)) return result;
  if (index < 0 || index >= pThis->atomNum[pos])
  {
    return IndexError;
  }
  AtomCopy(pDestAtom, &pThis->atoms[pos][index]);
  return Success;
}


int AtomParamsGetAtomParamByResidueName(AtomParams* pThis, char* residueName, char* atomName, Atom* pDestAtom)
{
  int pos, atomIndex;
  int result = StringArrayFind(&pThis->resNames, residueName, &pos);
  if (FAILED(result))  return result;
  for (atomIndex = 0;atomIndex < pThis->atomNum[pos];atomIndex++)
  {
    if (strcmp(pThis->atoms[pos][atomIndex].name, atomName) == 0) break;
  }
  if (atomIndex == pThis->atomNum[pos])
  { // atom not found
    return DataNotExistError;
  }
  AtomCopy(pDestAtom, &pThis->atoms[pos][atomIndex]);
  return Success;
}


int AtomParamsRead(AtomParams* pAtomParam, char* file)
{
  char errMsg[MAX_LEN_ERR_MSG+1];
  int result = AtomParamsAddFromFile(pAtomParam, file);
  if (FAILED(result))
  {
    sprintf(errMsg, "in file %s line %d, failed to read parameter file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, result);
    return result;
  }
  return Success;
}

