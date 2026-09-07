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

#include "Sequence.h"
#include "Structure.h"
#include "EnergyFunction.h"
#include <time.h>
#include <string.h>
#include <ctype.h>

extern double TORSION_DEVIATION_THRESHOLD;
extern char MOL2[MAX_LEN_FILE_NAME+1];
extern char DES_CHAINS[MAX_LEN_ONE_LINE_CONTENT+1];

extern BOOL FLAG_DESIGN_MON;
extern BOOL FLAG_DESIGN_PPI;
extern BOOL FLAG_DESIGN_PLI;
extern BOOL FLAG_DESIGN_ENZ;
extern BOOL FLAG_WRITE_MUTATIONS;

int SequenceCreate(Sequence* pThis)
{
  pThis->desSiteNum = 0;
  pThis->etot = 0;
  pThis->eevo = 0;
  pThis->ephy = 0;
  pThis->ebin = 0;
  pThis->ebpf = 0;
  pThis->mutNum = -1;
  pThis->seqExist = FALSE;
  pThis->unmetConsNum = 0;
  IntArrayCreate(&pThis->rotNdxs, 0);
  strcpy(pThis->pdbToSaveSeq, "None");
  return Success;
}


int SequenceDestroy(Sequence* pThis)
{
  pThis->desSiteNum = 0;
  pThis->etot = 0;
  pThis->eevo = 0;
  pThis->ephy = 0;
  pThis->ebin = 0;
  pThis->ebpf = 0;
  pThis->mutNum = -1;
  pThis->seqExist = FALSE;
  pThis->unmetConsNum = 0;
  IntArrayDestroy(&pThis->rotNdxs);
  return Success;
}


int SequenceCopy(Sequence* pThis, Sequence* pOther)
{
  pThis->desSiteNum = pOther->desSiteNum;
  pThis->etot = pOther->etot;
  pThis->ephy = pOther->ephy;
  pThis->ebin = pOther->ebin;
  pThis->eevo = pOther->eevo;
  pThis->ebpf = pOther->ebpf;
  pThis->mutNum = pOther->mutNum;
  pThis->seqExist = pOther->seqExist;
  pThis->unmetConsNum = pOther->unmetConsNum;
  IntArrayCopy(&pThis->rotNdxs, &pOther->rotNdxs);
  strcpy(pThis->pdbToSaveSeq, pOther->pdbToSaveSeq);
  return Success;
}


int SequenceWriteRotamerIndex(Sequence* pSeq, Structure* pStruct, int ndx, FILE* pFile)
{
  BOOL firstchain = TRUE;
  int iniChainNdx = -1;
  for (int i = 0;i < StructureGetChainNum(pStruct);i++) 
  {
    Chain* pChain = StructureGetChain(pStruct, i);
    for (int j = 0;j < ChainGetResidueNum(pChain);j++) 
    {
      int siteNdx = StructureGetDesignSiteIndex(pStruct, i, j);
      if (siteNdx == -1)
      {
        continue;
      }
      DesignSite* pSite = StructureGetDesignSite(pStruct, siteNdx);
      Rotamers* pSet = DesignSiteGetRotamers(pSite);
      int rotNdx = IntArrayGet(&pSeq->rotNdxs, siteNdx);
      if (i != iniChainNdx) 
      {
        if (firstchain) 
        {
          fprintf(pFile, "%d", rotNdx);
          firstchain = FALSE;
        }
        else
        {
          fprintf(pFile, ";%d", rotNdx);
        }
        iniChainNdx = i;
      }
      else
      {
        fprintf(pFile, ",%d", rotNdx);
      }
    }
  }
  fprintf(pFile, " %d\n", ndx);
  return Success;
}


int SequenceWriteDesignSummaryToFile(Sequence* pSeq, Structure* pStruct, int ndx, FILE* pFile)
{
  BOOL firstchain = TRUE;
  int iniChainNdx = -1;
  int mutResNum = 0;
  int sameResNum = 0;
  char mutStr[2*MAX_LEN_FILE_NAME+1] = "";
  BOOL firstMut = TRUE;
  for (int i = 0; i < StructureGetChainNum(pStruct);i++)
  {
    Chain* pChain = StructureGetChain(pStruct, i);
    if (!strstr(DES_CHAINS, ChainGetName(pChain)))
    {
      continue;
    }
    for (int j = 0;j < ChainGetResidueNum(pChain);j++)
    {
      Residue* pRes = ChainGetResidue(pChain, j);
      char res = AA3ToAA1(ResidueGetName(pRes));
      int siteNdx = StructureGetDesignSiteIndex(pStruct, i, j);
      if (siteNdx == -1) 
      { // current site is fixed (not a design site)
        continue;
      }
      DesignSite* pSite = StructureGetDesignSite(pStruct, siteNdx);
      Rotamers* pSet = DesignSiteGetRotamers(pSite);
      int rotIdx = IntArrayGet(&pSeq->rotNdxs, siteNdx);
      Rotamer* pRot = RotamersGet(pSet, rotIdx);
      char rot = AA3ToAA1(RotamerGetType(pRot));
      if (ChainGetType(pChain) == Type_Chain_Protein && ResidueGetDesignType(pRes) == Type_DesType_Mutable) 
      {
        if (rot == res)
        {
          sameResNum++;
        }
        else 
        {
          if (FLAG_DESIGN_PLI || FLAG_DESIGN_ENZ || FLAG_WRITE_MUTATIONS) 
          {
            if (firstMut == TRUE) 
            {
              sprintf(mutStr, "%s%c%c%d%c", mutStr, res, ResidueGetChainName(pRes)[0], ResidueGetPosInChain(pRes), rot);
              firstMut = FALSE;
            }
            else 
            {
              sprintf(mutStr, "%s,%c%c%d%c", mutStr, res, ResidueGetChainName(pRes)[0], ResidueGetPosInChain(pRes), rot);
            }
          }
        }
        mutResNum++;
      }
    }
  }

  if(strcmp(mutStr, "") == 0)
  {
    strcpy(mutStr, "Wildtype");
  }

  double seq_identity = (mutResNum > 0 ? (double)sameResNum / mutResNum : 1.0);
  if (FLAG_DESIGN_PLI || FLAG_DESIGN_ENZ || FLAG_WRITE_MUTATIONS) 
  {
    fprintf(pFile, "%s\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%s\n",
      pSeq->pdbToSaveSeq, seq_identity, pSeq->etot, pSeq->eevo, pSeq->ephy, pSeq->ebin, pSeq->unmetConsNum, pSeq->mutNum, pSeq->seqExist, mutStr);
  }
  else 
  {
    fprintf(pFile, "%s\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%s\n",
      pSeq->pdbToSaveSeq, seq_identity, pSeq->etot, pSeq->eevo, pSeq->ephy, pSeq->ebin, pSeq->unmetConsNum, pSeq->mutNum, pSeq->seqExist, "mutations_not_shown");
  }

  return Success;
}


int SequenceWriteDesignSequenceToFile(Sequence* pSeq, Structure* pStruct, int ndx, FILE* pFile)
{
  BOOL firstchain = TRUE;
  int iniChainNdx = -1;
  char fullSeq[MAX_LEN_ONE_LINE_CONTENT+1] = "";
  char mutSeq[MAX_LEN_ONE_LINE_CONTENT+1] = "";
  for (int i = 0; i < StructureGetChainNum(pStruct); i++) 
  {
    Chain* pChain = StructureGetChain(pStruct, i);
    if (!strstr(DES_CHAINS, ChainGetName(pChain)))
    {
      continue;
    }
    if (ChainGetType(pChain) != Type_Chain_Protein)
    {
      continue;
    }
    for (int j = 0; j < ChainGetResidueNum(pChain); j++) 
    {
      Residue* pRes = ChainGetResidue(pChain, j);
      char res = AA3ToAA1(ResidueGetName(pRes));
      int siteNdx = StructureGetDesignSiteIndex(pStruct, i, j);
      if (siteNdx == -1) 
      { // current site is fixed (not design site)
        if (i != iniChainNdx) 
        {
          if (firstchain) 
          {
            sprintf(fullSeq, "%s%c", fullSeq, res);
            firstchain = FALSE;
          }
          else
          {
            sprintf(fullSeq, "%s;%c", fullSeq, res);
          }
          iniChainNdx = i;
        }
        else 
        {
          sprintf(fullSeq, "%s%c", fullSeq, res);
        }
        continue;
      }
      DesignSite* pSite = StructureGetDesignSite(pStruct, siteNdx);
      Rotamers* pSet = DesignSiteGetRotamers(pSite);
      int rotIdx = IntArrayGet(&pSeq->rotNdxs, siteNdx);
      Rotamer* pRot = RotamersGet(pSet, rotIdx);
      char rot = AA3ToAA1(RotamerGetType(pRot));
      if (i != iniChainNdx) 
      {
        if (firstchain) 
        {
          sprintf(fullSeq, "%s%c", fullSeq, rot);
          if (ResidueGetDesignType(pRes) == Type_DesType_Mutable)
          { 
            sprintf(mutSeq, "%s%c", mutSeq, rot);
          }
          firstchain = FALSE;
        }
        else 
        {
          sprintf(fullSeq, "%s;%c", fullSeq, rot);
          if (ResidueGetDesignType(pRes) == Type_DesType_Mutable)
          {
            sprintf(mutSeq, "%s;%c", mutSeq, rot);
          }
        }
        iniChainNdx = i;
      }
      else {
        sprintf(fullSeq, "%s%c", fullSeq, rot);
        if (ResidueGetDesignType(pRes) == Type_DesType_Mutable) 
        {
          sprintf(mutSeq, "%s%c", mutSeq, rot);
        }
      }
    }
  }
  if (FLAG_DESIGN_PLI || FLAG_DESIGN_ENZ || FLAG_WRITE_MUTATIONS)
  {
    fprintf(pFile, "%s\t%s\t%s\n", pSeq->pdbToSaveSeq, fullSeq, mutSeq);
  }
  else
  {
    fprintf(pFile, "%s\t%s\t%s\n", pSeq->pdbToSaveSeq, fullSeq, "N.A.");
  }

  return Success;
}


int SequenceWriteStructureInPDB(Sequence* pSeq, Structure* pStruct, char* file)
{
  int result = Success;
  char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
  FILE* pFile = fopen(file, "w");
  if (pFile == NULL) 
  {
    result = IOError;
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0; i < StructureGetChainNum(pStruct);i++) 
  {
    Chain* pChain = StructureGetChain(pStruct, i);
    for (int j = 0;j < ChainGetResidueNum(pChain);j++) 
    {
      Residue* pRes = ChainGetResidue(pChain, j);
      if (pRes->desType == Type_DesType_Fixed) 
      {
        if (ChainGetType(pChain) == Type_Chain_Protein || ChainGetType(pChain) == Type_Chain_DNA || ChainGetType(pChain) == Type_Chain_RNA) 
        {
          ResidueShowInPDBFormat(pRes, "ATOM", ResidueGetChainName(pRes), 1, ResidueGetPosInChain(pRes), pFile);
        }
        else if (ChainGetType(pChain) != Type_Chain_SmallMol) 
        {
          ResidueShowInPDBFormat(pRes, "HETATM", ResidueGetChainName(pRes), 1, ResidueGetPosInChain(pRes), pFile);
        }
      }
      else 
      {
        int siteIdx = StructureGetDesignSiteIndex(pStruct, i, j);
        DesignSite* pSite = StructureGetDesignSite(pStruct, siteIdx);
        Rotamers* pSet = DesignSiteGetRotamers(pSite);
        int rotIdx = IntArrayGet(&pSeq->rotNdxs, siteIdx);
        Rotamer* pRot = RotamersGet(pSet, rotIdx);
        RotamerRestore(pRot, pSet);
        if (ChainGetType(pChain) == Type_Chain_Protein || ChainGetType(pChain) == Type_Chain_DNA || ChainGetType(pChain) == Type_Chain_RNA) 
        {
          RotamerShowInPDBFormat(pRot, "ATOM", RotamerGetChainName(pRot), 1, RotamerGetPosInChain(pRot), pFile);
        }
        else if (ChainGetType(pChain) != Type_Chain_SmallMol) 
        {
          RotamerShowInPDBFormat(pRot, "HETATM", RotamerGetChainName(pRot), 1, RotamerGetPosInChain(pRot), pFile);
        }
        RotamerExtract(pRot);
      }
    }
    AddChainTER(pFile);
  }
  fclose(pFile);
  return result;
}

int SequenceWriteAllSitesInPDB(Sequence* pSeq, Structure* pStruct, char* file)
{
  int result = Success;
  char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
  FILE* pFile = fopen(file, "w");
  if (pFile == NULL) 
  {
    result = IOError;
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0;i < StructureGetDesignSiteNum(pStruct);i++) 
  {
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    Rotamers* pSet = DesignSiteGetRotamers(pSite);
    int rotIdx = IntArrayGet(&pSeq->rotNdxs, i);
    Rotamer* pRotamer = RotamersGet(pSet, rotIdx);
    RotamerRestore(pRotamer, pSet);
    RotamerShowInPDBFormat(pRotamer, "ATOM", RotamerGetChainName(pRotamer), 1, RotamerGetPosInChain(pRotamer), pFile);
    RotamerExtract(pRotamer);
  }
  fclose(pFile);
  return result;
}


int SequenceWriteMutableSitesInPDB(Sequence* pSeq, Structure* pStruct, char* file)
{
  int result = Success;
  char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
  FILE* pFile = fopen(file, "w");
  if (pFile == NULL) 
  {
    result = IOError;
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0;i < StructureGetDesignSiteNum(pStruct);i++) 
  {
    DesignSite* pSite = StructureGetDesignSite(pStruct, i);
    if (ResidueGetDesignType(pSite->pRes) == Type_DesType_Mutable || ResidueGetDesignType(pSite->pRes) == Type_DesType_SmallMol) 
    {
      Rotamers* pSet = DesignSiteGetRotamers(pSite);
      int rotIdx = IntArrayGet(&pSeq->rotNdxs, i);
      Rotamer* pRotamer = RotamersGet(pSet, rotIdx);
      RotamerRestore(pRotamer, pSet);
      RotamerShowInPDBFormat(pRotamer, "ATOM", RotamerGetChainName(pRotamer), 1, RotamerGetPosInChain(pRotamer), pFile);
      RotamerExtract(pRotamer);
    }
  }
  fclose(pFile);
  return result;
}



int SequenceWriteLigandInMol2(Sequence* pSeq, Structure* pStruct, char* file)
{
  int result = Success;
  char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];

  StringArray lines;
  StringArrayCreate(&lines);
  FILE* pfile = fopen(MOL2, "r");
  char buffer[MAX_LEN_ONE_LINE_CONTENT+1];
  int start = -1, end = -1;
  int num = 0;
  while (fgets(buffer, MAX_LEN_ONE_LINE_CONTENT, pfile)) 
  {
    StringArrayAppend(&lines, buffer);
    if (strstr(buffer, "@<TRIPOS>ATOM") != NULL) 
    {
      start = num + 1;
    }
    else if (strstr(buffer, "@<TRIPOS>BOND") != NULL) 
    {
      end = num - 1;
    }
    num++;
  }
  fclose(pfile);

  FILE* pMol2 = fopen(file, "w");
  if (pMol2 == NULL) 
  {
    result = IOError;
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, result);
    return result;
  }

  for (int i = 0; i < StructureGetChainNum(pStruct);i++) 
  {
    Chain* pChain = StructureGetChain(pStruct, i);
    for (int j = 0;j < ChainGetResidueNum(pChain);j++) 
    {
      Residue* pRes = ChainGetResidue(pChain, j);
      if (ResidueGetDesignType(pRes) == Type_DesType_Fixed) 
      {
        if (ChainGetType(pChain) == Type_Chain_SmallMol) 
        {
          for (int l = 0;l < StringArrayGetNum(&lines);l++) 
          {
            if (l<start || l>end) 
            {
              fprintf(pMol2, "%s", StringArrayGet(&lines, l));
            }
            else 
            {
              char line[MAX_LEN_ONE_LINE_CONTENT+1];
              strcpy(line, StringArrayGet(&lines, l));
              StringArray words;
              StringArrayCreate(&words);
              StringArraySplitString(&words, line, ' ');
              XYZ xyz = ResidueGetAtom(pRes, l - start)->xyz;
              fprintf(pMol2, "%7s %-10s %10.3f %10.3f %10.3f %-10s %s %s %10s\n",
                StringArrayGet(&words, 0),
                StringArrayGet(&words, 1),
                xyz.X,
                xyz.Y,
                xyz.Z,
                StringArrayGet(&words, 5),
                StringArrayGet(&words, 6),
                StringArrayGet(&words, 7),
                StringArrayGet(&words, 8));
              StringArrayDestroy(&words);
            }
          }
        }
      }
      else 
      {
        for (int k = 0;k < pSeq->desSiteNum;k++) 
        {
          DesignSite* pSite = StructureGetDesignSite(pStruct, k);
          Rotamers* pSet = DesignSiteGetRotamers(pSite);
          int rotNdx = IntArrayGet(&pSeq->rotNdxs, k);
          Rotamer* pOptRot = RotamersGet(pSet, rotNdx);
          if (pSite->chnNdx == i && pSite->resNdx == j) 
          {
            RotamerRestore(pOptRot, pSet);
            if (ChainGetType(pChain) == Type_Chain_SmallMol) 
            {
              for (int l = 0;l < StringArrayGetNum(&lines);l++) 
              {
                if (l<start || l>end) 
                {
                  fprintf(pMol2, "%s", StringArrayGet(&lines, l));
                }
                else 
                {
                  char line[MAX_LEN_ONE_LINE_CONTENT+1];
                  strcpy(line, StringArrayGet(&lines, l));
                  StringArray words;
                  StringArrayCreate(&words);
                  StringArraySplitString(&words, line, ' ');
                  XYZ xyz = RotamerGetAtom(pOptRot, l - start)->xyz;
                  fprintf(pMol2, "%7s %-10s %10.3f %10.3f %10.3f %-10s %s %s %10s\n",
                    StringArrayGet(&words, 0),
                    StringArrayGet(&words, 1),
                    xyz.X,
                    xyz.Y,
                    xyz.Z,
                    StringArrayGet(&words, 5),
                    StringArrayGet(&words, 6),
                    StringArrayGet(&words, 7),
                    StringArrayGet(&words, 8));
                  StringArrayDestroy(&words);
                }
              }
            }
            RotamerExtract(pOptRot);
            break;
          }
        }
      }
    }
  }
  fclose(pMol2);
  StringArrayDestroy(&lines);

  return result;
}


int StructureGetWholeSequence(Sequence* pSeq, Structure* pStruct, char* seqstr)
{
  int index = 0;
  for (int i = 0;i < StructureGetChainNum(pStruct);i++) 
  {
    Chain* pChain = StructureGetChain(pStruct, i);
    if (strcmp(DES_CHAINS, ChainGetName(pChain))) 
    {
      continue;
    }
    for (int j = 0;j < ChainGetResidueNum(pChain);j++) 
    {
      Residue* pResi = ChainGetResidue(pChain, j);
      if (ResidueGetDesignType(pResi) == Type_DesType_Fixed) 
      {
        seqstr[index++] = AA3ToAA1(ResidueGetName(pResi));
      }
      else 
      {
        int siteIdx = StructureGetDesignSiteIndex(pStruct, i, j);
        DesignSite* pSite = StructureGetDesignSite(pStruct, siteIdx);
        Rotamers* pSet = DesignSiteGetRotamers(pSite);
        int rotIdx = IntArrayGet(&pSeq->rotNdxs, siteIdx);
        Rotamer* pRot = RotamersGet(pSet, rotIdx);
        seqstr[index++] = AA3ToAA1(RotamerGetType(pRot));
      }
    }
  }
  seqstr[index] = '\0';
  return Success;
}


int StructureWriteFirstLigandRotamerIntoMol2(Structure* pStruct, char* file)
{
  int result = Success;
  char errMsg[MAX_LEN_ONE_LINE_CONTENT+1];
  StringArray lines;
  StringArrayCreate(&lines);
  FILE* pfile = fopen(MOL2, "r");
  char buffer[MAX_LEN_ONE_LINE_CONTENT+1];
  int start = 100000, end = -1;
  int num = 0;
  int atomend = 0;
  while (fgets(buffer, MAX_LEN_ONE_LINE_CONTENT, pfile)) 
  {
    StringArrayAppend(&lines, buffer);
    if (strstr(buffer, "@<TRIPOS>ATOM") != NULL) 
    {
      start = num + 1;
    }
    else if (strstr(buffer, "@<TRIPOS>") != NULL && num > start && atomend == 0) 
    {
      end = num - 1;
      atomend = 1;
    }
    num++;
  }
  fclose(pfile);

  FILE* pMol2 = fopen(file, "w");
  if (pMol2 == NULL) 
  {
    result = IOError;
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, file);
    TraceError(errMsg, result);
    return result;
  }

  for (int i = 0; i < StructureGetChainNum(pStruct);i++) 
  {
    Chain* pChain = StructureGetChain(pStruct, i);
    for (int j = 0;j < ChainGetResidueNum(pChain);j++) 
    {
      Residue* pRes = ChainGetResidue(pChain, j);
      if (pRes->desType == Type_DesType_Fixed) 
      {
        if (pChain->type == Type_Chain_SmallMol) 
        {
          for (int l = 0;l < StringArrayGetNum(&lines);l++) 
          {
            if (l<start || l>end) 
            {
              fprintf(pMol2, "%s", StringArrayGet(&lines, l));
            }
            else 
            {
              char line[MAX_LEN_ONE_LINE_CONTENT+1];
              strcpy(line, StringArrayGet(&lines, l));
              StringArray words;
              StringArrayCreate(&words);
              StringArraySplitString(&words, line, ' ');
              XYZ xyz = ResidueGetAtom(pRes, l - start)->xyz;
              fprintf(pMol2, "%7s %-10s %10.3f %10.3f %10.3f %-10s %s %s %10s\n",
                StringArrayGet(&words, 0),
                StringArrayGet(&words, 1),
                xyz.X,
                xyz.Y,
                xyz.Z,
                StringArrayGet(&words, 5),
                StringArrayGet(&words, 6),
                StringArrayGet(&words, 7),
                StringArrayGet(&words, 8));
              StringArrayDestroy(&words);
            }
          }
        }
      }
      else 
      {
        for (int k = 0;k < pStruct->desSiteNum;k++) 
        {
          DesignSite* pSite = StructureGetDesignSite(pStruct, k);
          Rotamers* pSet = DesignSiteGetRotamers(pSite);
          Rotamer* pOptRot = RotamersGet(pSet, 0);
          if (pSite->chnNdx == i && pSite->resNdx == j) 
          {
            RotamerRestore(pOptRot, pSet);
            if (ChainGetType(pChain) == Type_Chain_SmallMol) 
            {
              for (int l = 0;l < StringArrayGetNum(&lines);l++) 
              {
                if (l<start || l>end) 
                {
                  fprintf(pMol2, "%s", StringArrayGet(&lines, l));
                }
                else 
                {
                  char line[MAX_LEN_ONE_LINE_CONTENT+1];
                  strcpy(line, StringArrayGet(&lines, l));
                  StringArray words;
                  StringArrayCreate(&words);
                  StringArraySplitString(&words, line, ' ');
                  XYZ xyz = RotamerGetAtom(pOptRot, l - start)->xyz;
                  fprintf(pMol2, "%7s %-10s %10.3f %10.3f %10.3f %-10s %s %s %10s\n",
                    StringArrayGet(&words, 0),
                    StringArrayGet(&words, 1),
                    xyz.X,
                    xyz.Y,
                    xyz.Z,
                    StringArrayGet(&words, 5),
                    StringArrayGet(&words, 6),
                    StringArrayGet(&words, 7),
                    StringArrayGet(&words, 8));
                  StringArrayDestroy(&words);
                }
              }
            }
            RotamerExtract(pOptRot);
            break;
          }
        }
      }
    }
  }
  fclose(pMol2);
  StringArrayDestroy(&lines);

  return result;
}