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

#include "Chain.h"
#include <string.h>

int ChainTypeConvertFromString(char* typeName, Type_Chain* type)
{
  if (strcmp(typeName, "PROTEIN") == 0)
  {
    *type = Type_Chain_Protein;
  }
  else if (strcmp(typeName, "SMALLMOL") == 0)
  {
    *type = Type_Chain_SmallMol;
  }
  else if (strcmp(typeName, "METALION") == 0)
  {
    *type = Type_Chain_MetalIon;
  }
  else if (strcmp(typeName, "DNA") == 0)
  {
    *type = Type_Chain_DNA;
  }
  else if (strcmp(typeName, "RNA") == 0)
  {
    *type = Type_Chain_RNA;
  }
  else if (strcmp(typeName, "WATER") == 0)
  {
    *type = Type_Chain_Water;
  }
  else
  {
    *type = Type_Chain_Unknown;
  }

  return Success;
}


Type_Chain ChainTypeIdentifiedFromResidueName(char* resiName)
{
  if (!strcmp(resiName, "ALA")
    || !strcmp(resiName, "ARG")
    || !strcmp(resiName, "ASN")
    || !strcmp(resiName, "ASP")
    || !strcmp(resiName, "CYS")
    || !strcmp(resiName, "GLN")
    || !strcmp(resiName, "GLU")
    || !strcmp(resiName, "GLY")
    || !strcmp(resiName, "HIS")
    || !strcmp(resiName, "HSD")
    || !strcmp(resiName, "HSE")
    || !strcmp(resiName, "HSP")
    || !strcmp(resiName, "ILE")
    || !strcmp(resiName, "LEU")
    || !strcmp(resiName, "LYS")
    || !strcmp(resiName, "MET")
    || !strcmp(resiName, "PHE")
    || !strcmp(resiName, "PRO")
    || !strcmp(resiName, "SER")
    || !strcmp(resiName, "THR")
    || !strcmp(resiName, "TRP")
    || !strcmp(resiName, "TYR")
    || !strcmp(resiName, "VAL"))
  {
    return Type_Chain_Protein;
  }
  else if (!strcmp(resiName, "DA")
    || !strcmp(resiName, "DG")
    || !strcmp(resiName, "DC")
    || !strcmp(resiName, "DT")
    || !strcmp(resiName, "DN"))
  {
    return Type_Chain_DNA;
  }
  else if (!strcmp(resiName, "A")
    || !strcmp(resiName, "G")
    || !strcmp(resiName, "C")
    || !strcmp(resiName, "U")
    || !strcmp(resiName, "N"))
  {
    return Type_Chain_RNA;
  }
  else if (!strcmp(resiName, "HOH")
    || !strcmp(resiName, "H2O")
    || !strcmp(resiName, "WAT"))
  {
    return Type_Chain_Water;
  }
  else
  {
    return Type_Chain_Unknown;
  }
}


int ChainCreate(Chain* pThis)
{
  ChainSetName(pThis, "");
  pThis->resNum = 0;
  pThis->residues = NULL;
  return Success;
}


int ChainDestroy(Chain* pThis)
{
  for (int i = 0;i < pThis->resNum;i++)
  {
    ResidueDestroy(&pThis->residues[i]);
  }
  free(pThis->residues);
  pThis->residues = NULL;
  pThis->resNum = 0;
  return Success;
}


int ChainCopy(Chain* pThis, Chain* pOther)
{
  ChainDestroy(pThis);
  ChainSetName(pThis, pOther->name);
  ChainSetType(pThis, pOther->type);

  pThis->resNum = pOther->resNum;
  pThis->residues = (Residue*)malloc(sizeof(Residue) * pThis->resNum);
  for (int i = 0;i < pThis->resNum;i++)
  {
    int result;
    ResidueCreate(&pThis->residues[i]);
    result = ResidueCopy(&pThis->residues[i], &pOther->residues[i]);
    if (FAILED(result))
      return result;
  }
  return Success;
}


char* ChainGetName(Chain* pThis)
{
  return pThis->name;
}


int ChainSetName(Chain* pThis, char* newName)
{
  if (strlen(newName) > MAX_LEN_CHAIN_NAME)
  {
    char errMsg[MAX_LEN_ERR_MSG+1];
    int result = ValueError;
    sprintf(errMsg, "in file %s line %d, chain name %s is too long", __FILE__, __LINE__, newName);
    TraceError(errMsg, result);
    return result;
  }
  strcpy(pThis->name, newName);
  return Success;
}


Type_Chain ChainGetType(Chain* pThis)
{
  return pThis->type;
}


int ChainSetType(Chain* pThis, Type_Chain newType)
{
  pThis->type = newType;
  return Success;
}


int ChainGetResidueNum(Chain* pThis)
{
  return pThis->resNum;
}


Residue* ChainGetResidue(Chain* pThis, int index)
{
  if (index < 0 || index >= ChainGetResidueNum(pThis)) return NULL;
  else return &pThis->residues[index];
}


int ChainInsertResidue(Chain* pThis, int index, Residue* pNewRes)
{
  if (index<0 || index>pThis->resNum) return IndexError;
  int newNum = pThis->resNum + 1;
  pThis->residues = (Residue*)realloc(pThis->residues, sizeof(Residue) * newNum);
  ResidueCreate(&pThis->residues[newNum - 1]);
  for (int i = newNum - 1;i > index;i--)
  {
    ResidueCopy(&pThis->residues[i], &pThis->residues[i - 1]);
    ResidueSetPosInChain(&pThis->residues[i], pThis->residues[i - 1].posInChain);
  }
  ResidueCopy(&pThis->residues[index], pNewRes);
  ResidueSetChainName(&pThis->residues[index], pThis->name);
  ResidueSetPosInChain(&pThis->residues[index], pNewRes->posInChain);
  pThis->resNum = newNum;

  return Success;
}


int ChainRemoveResidue(Chain* pThis, int index)
{
  if (index < 0 || index >= pThis->resNum) return IndexError;
  for (int i = index;i < pThis->resNum - 1;i++)
  {
    ResidueCopy(&pThis->residues[i], &pThis->residues[i + 1]);
    ResidueSetPosInChain(&pThis->residues[i], ResidueGetPosInChain(&pThis->residues[i + 1]));
  }
  ResidueDestroy(&pThis->residues[pThis->resNum - 1]);
  (pThis->resNum)--;
  return Success;
}


int ChainAppendResidue(Chain* pThis, Residue* pNewResi)
{
  return ChainInsertResidue(pThis, pThis->resNum, pNewResi);
}


int ChainCalcAllAtomXYZ(Chain* pThis, ResTopos* topos)
{
  for (int i = 0;i < ChainGetResidueNum(pThis);i++)
  {
    Residue* prevResi = ChainGetResidue(pThis, i - 1);
    Residue* nextResi = ChainGetResidue(pThis, i + 1);
    Residue* curResi = ChainGetResidue(pThis, i);
    int result = ResidueCalcAllAtomXYZ(curResi, topos, prevResi, nextResi);
    if (FAILED(result))
    {
      char errMsg[MAX_LEN_ERR_MSG+1];
      sprintf(errMsg, "in file %s line %d, failed to calculate some atoms' coordinates", __FILE__, __LINE__);
      TraceError(errMsg, result);
      return result;
    }
  }
  return Success;
}


int ChainShowInPDBFormat(Chain* pThis, int resiIndex, int atomNdx, FILE* pOut)
{
  char header[7];
  char chainName[2];
  if (pThis->type == Type_Chain_Protein)
  {
    strcpy(header, "ATOM  ");
    strcpy(chainName, " ");
    chainName[0] = pThis->name[strlen(pThis->name) - 1];
  }
  else
  {
    strcpy(header, "HETATM");
    strcpy(chainName, " ");
  }
  int num = atomNdx;
  for (int i = 0;i < pThis->resNum;i++)
  {
    ResidueShowInPDBFormat(&pThis->residues[i], header, chainName, num, ResidueGetPosInChain(&pThis->residues[i]) + resiIndex, pOut);
    num += ResidueGetAtomNum(&pThis->residues[i]);
  }
  return Success;
}


int ChainFindResidueIndexByPosInChain(Chain* pThis, int posInchain, int* index)
{
  for (int i = 0; i < ChainGetResidueNum(pThis); ++i)
  {
    Residue* pResi = ChainGetResidue(pThis, i);
    if (ResidueGetPosInChain(pResi) == posInchain)
    {
      *index = i;
      return Success;
    }
  }
  return DataNotExistError;
}


int ChainShowAtomParameter(Chain* pThis)
{
  for (int i = 0; i < ChainGetResidueNum(pThis); i++)
  {
    ResidueShowAtomParameter(ChainGetResidue(pThis, i));
  }
  return Success;
}


int ChainShowBondInformation(Chain* pThis)
{
  printf("------------------Bonds information for Chain %s--------------\n", ChainGetName(pThis));
  for (int i = 0; i < ChainGetResidueNum(pThis); i++)
  {
    Residue* pResidue = ChainGetResidue(pThis, i);
    printf("Residue %s PosInChain %d:\n", ResidueGetName(pResidue), ResidueGetPosInChain(pResidue));
    ResidueShowBondInformation(pResidue);
  }
  return Success;
}

