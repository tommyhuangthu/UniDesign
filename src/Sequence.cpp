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

extern double CUT_TORSION_DEVIATION;
extern char MOL2[MAX_LEN_FILE_NAME + 1];
extern char DES_CHAINS[MAX_LEN_ONE_LINE_CONTENT + 1];

int SequenceCreate(Sequence* pThis)
{
  pThis->desSiteCount = 0;
  pThis->etot = 0;
  pThis->eevo = 0;
  pThis->ephy = 0;
  pThis->ebin = 0;
  pThis->ebpf = 0;
  IntArrayCreate(&pThis->rotNdxs, 0);
  pThis->numOfUnsatisfiedCons = 0;
  strcpy(pThis->fileToSaveThisSeq, "None");
  return Success;
}


int SequenceDestroy(Sequence* pThis)
{
  pThis->desSiteCount = 0;
  pThis->etot = 0;
  pThis->eevo = 0;
  pThis->ephy = 0;
  pThis->ebin = 0;
  pThis->ebpf = 0;
  IntArrayDestroy(&pThis->rotNdxs);
  pThis->numOfUnsatisfiedCons = 0;
  return Success;
}


int SequenceCopy(Sequence* pThis, Sequence* pOther)
{
  pThis->desSiteCount = pOther->desSiteCount;
  pThis->etot = pOther->etot;
  pThis->ephy = pOther->ephy;
  pThis->ebin = pOther->ebin;
  pThis->eevo = pOther->eevo;
  pThis->ebpf = pOther->ebpf;
  IntArrayCopy(&pThis->rotNdxs, &pOther->rotNdxs);
  pThis->numOfUnsatisfiedCons = pOther->numOfUnsatisfiedCons;
  strcpy(pThis->fileToSaveThisSeq, pOther->fileToSaveThisSeq);
  return Success;
}


int SequenceWriteDesignRotamer(Sequence* pThis, Structure* pStructure, int seqIdx, FILE* pFile)
{
  BOOL firstchain = TRUE;
  int flagIndex = -1;
  for (int i = 0;i < StructureGetChainCount(pStructure);i++)
  {
    Chain* pChain = StructureGetChain(pStructure, i);
    for (int j = 0;j < ChainGetResidueCount(pChain);j++)
    {
      int siteIdx = StructureGetDesignSiteIndex(pStructure, i, j);
      if (siteIdx == -1) continue;
      DesignSite* pSite = StructureGetDesignSite(pStructure, siteIdx);
      RotamerSet* pSet = DesignSiteGetRotamers(pSite);
      int rotIdx = IntArrayGet(&pThis->rotNdxs, siteIdx);
      if (i != flagIndex)
      {
        if (firstchain)
        {
          fprintf(pFile, "%d", rotIdx);
          firstchain = FALSE;
        }
        else fprintf(pFile, ";%d", rotIdx);
        flagIndex = i;
      }
      else fprintf(pFile, ",%d", rotIdx);
    }
  }
  fprintf(pFile, " %d\n", seqIdx);
  return Success;
}


int SequenceWriteDesignFasta(Sequence* pThis, Structure* pStructure, int seqIdx, FILE* pFile)
{
  BOOL firstchain = TRUE;
  int chainFlag = -1;
  int totResCount = 0;
  int samResCount = 0;
  for (int i = 0; i < StructureGetChainCount(pStructure);i++)
  {
    Chain* pChain = StructureGetChain(pStructure, i);
    if (!strstr(DES_CHAINS, ChainGetName(pChain))) continue;
    for (int j = 0;j < ChainGetResidueCount(pChain);j++)
    {
      Residue* pResi = ChainGetResidue(pChain, j);
      char res = AA3ToAA1(ResidueGetName(pResi));
      int siteIdx = StructureGetDesignSiteIndex(pStructure, i, j);
      if (siteIdx == -1)
      {
        if (i != chainFlag)
        {
          if (firstchain)
          {
            if (ChainGetType(pChain) == Type_Chain_Protein)
            {
              fprintf(pFile, "%c", res);
            }
            firstchain = FALSE;
          }
          else
          {
            if (ChainGetType(pChain) == Type_Chain_Protein)
            {
              fprintf(pFile, ";%c", res);
            }
          }
          chainFlag = i;
        }
        else
        {
          if (ChainGetType(pChain) == Type_Chain_Protein)
          {
            fprintf(pFile, "%c", res);
          }
        }
        continue;
      }

      DesignSite* pSite = StructureGetDesignSite(pStructure, siteIdx);
      RotamerSet* pSet = DesignSiteGetRotamers(pSite);
      int rotIdx = IntArrayGet(&pThis->rotNdxs, siteIdx);
      Rotamer* pRot = RotamerSetGet(pSet, rotIdx);
      char rot = AA3ToAA1(RotamerGetType(pRot));

      if (i != chainFlag)
      {
        if (firstchain)
        {
          if (ChainGetType(pChain) == Type_Chain_Protein)
          {
            fprintf(pFile, "%c", rot);
          }
          firstchain = FALSE;
        }
        else
        {
          if (ChainGetType(pChain) == Type_Chain_Protein)
          {
            fprintf(pFile, ";%c", rot);
          }
        }
        chainFlag = i;
      }
      else
      {
        if (ChainGetType(pChain) == Type_Chain_Protein)
        {
          fprintf(pFile, "%c", rot);
        }
      }

      if (ChainGetType(pChain) == Type_Chain_Protein && ResidueGetDesignType(pResi) == Type_DesType_Mutable)
      {
        if (rot == res) samResCount++;
        totResCount++;
      }
    }
  }
  double seqid = (totResCount > 0 ? (double)samResCount / totResCount : 1.0);
  fprintf(pFile, " %50s %13.6f %13.6f %13.6f %13.6f %13.6f %13d\n",
    pThis->fileToSaveThisSeq, seqid, pThis->etot, pThis->eevo, pThis->ephy, pThis->ebin, pThis->numOfUnsatisfiedCons);

  return Success;
}


int DesignShowMinEnergyDesignStructure(Structure* pStructure, Sequence* pSeq, char* pdbfile)
{
  int result = Success;
  char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];
  FILE* pFile = fopen(pdbfile, "w");
  if (pFile == NULL)
  {
    result = IOError;
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, pdbfile);
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0; i < StructureGetChainCount(pStructure);i++)
  {
    Chain* pChain = StructureGetChain(pStructure, i);
    for (int j = 0;j < ChainGetResidueCount(pChain);j++)
    {
      Residue* pResi = ChainGetResidue(pChain, j);
      if (pResi->desType == Type_DesType_Fixed)
      {
        if (ChainGetType(pChain) == Type_Chain_Protein || ChainGetType(pChain) == Type_Chain_DNA || ChainGetType(pChain) == Type_Chain_RNA)
        {
          ResidueShowInPDBFormat(pResi, "ATOM", ResidueGetChainName(pResi), 1, ResidueGetPosInChain(pResi), pFile);
        }
        else if (ChainGetType(pChain) != Type_Chain_SmallMol)
        {
          ResidueShowInPDBFormat(pResi, "HETATM", ResidueGetChainName(pResi), 1, ResidueGetPosInChain(pResi), pFile);
        }
      }
      else
      {
        int siteIdx = StructureGetDesignSiteIndex(pStructure, i, j);
        DesignSite* pSite = StructureGetDesignSite(pStructure, siteIdx);
        RotamerSet* pSet = DesignSiteGetRotamers(pSite);
        int rotIdx = IntArrayGet(&pSeq->rotNdxs, siteIdx);
        Rotamer* pRot = RotamerSetGet(pSet, rotIdx);
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

int DesignShowMinEnergyDesignSites(Structure* pStructure, Sequence* sequence, char* pdbfile)
{
  int result = Success;
  char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];
  FILE* pFile = fopen(pdbfile, "w");
  if (pFile == NULL)
  {
    result = IOError;
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, pdbfile);
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0;i < StructureGetDesignSiteCount(pStructure);i++)
  {
    DesignSite* pSite = StructureGetDesignSite(pStructure, i);
    RotamerSet* pSet = DesignSiteGetRotamers(pSite);
    int rotIdx = IntArrayGet(&sequence->rotNdxs, i);
    Rotamer* pRotamer = RotamerSetGet(pSet, rotIdx);
    RotamerRestore(pRotamer, pSet);
    RotamerShowInPDBFormat(pRotamer, "ATOM", RotamerGetChainName(pRotamer), 1, RotamerGetPosInChain(pRotamer), pFile);
    RotamerExtract(pRotamer);
  }
  fclose(pFile);
  return result;
}


int DesignShowMinEnergyDesignMutableSites(Structure* pStructure, Sequence* sequence, char* pdbfile)
{
  int result = Success;
  char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];
  FILE* pFile = fopen(pdbfile, "w");
  if (pFile == NULL)
  {
    result = IOError;
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, pdbfile);
    TraceError(errMsg, result);
    return result;
  }
  for (int i = 0;i < StructureGetDesignSiteCount(pStructure);i++)
  {
    DesignSite* pSite = StructureGetDesignSite(pStructure, i);
    if (pSite->pRes->desType == Type_DesType_Mutable || pSite->pRes->desType == Type_DesType_SmallMol)
    {
      RotamerSet* pSet = DesignSiteGetRotamers(pSite);
      int rotIdx = IntArrayGet(&sequence->rotNdxs, i);
      Rotamer* pRotamer = RotamerSetGet(pSet, rotIdx);
      RotamerRestore(pRotamer, pSet);
      RotamerShowInPDBFormat(pRotamer, "ATOM", RotamerGetChainName(pRotamer), 1, RotamerGetPosInChain(pRotamer), pFile);
      RotamerExtract(pRotamer);
    }
  }
  fclose(pFile);
  return result;
}



int DesignShowMinEnergyDesignLigand(Structure* pStructure, Sequence* sequence, char* mol2file)
{
  int result = Success;
  char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];

  StringArray lines;
  StringArrayCreate(&lines);
  FILE* pfile = fopen(MOL2, "r");
  char buffer[MAX_LEN_ONE_LINE_CONTENT + 1];
  int start = -1, end = -1;
  int counter = 0;
  while (fgets(buffer, MAX_LEN_ONE_LINE_CONTENT, pfile))
  {
    StringArrayAppend(&lines, buffer);
    if (strstr(buffer, "@<TRIPOS>ATOM") != NULL)
    {
      start = counter + 1;
    }
    else if (strstr(buffer, "@<TRIPOS>BOND") != NULL)
    {
      end = counter - 1;
    }
    counter++;
  }
  fclose(pfile);

  FILE* pMol2 = fopen(mol2file, "w");
  if (pMol2 == NULL)
  {
    result = IOError;
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, mol2file);
    TraceError(errMsg, result);
    return result;
  }

  for (int i = 0; i < StructureGetChainCount(pStructure);i++)
  {
    Chain* pChain = StructureGetChain(pStructure, i);
    for (int j = 0;j < ChainGetResidueCount(pChain);j++)
    {
      Residue* pResi = ChainGetResidue(pChain, j);
      if (pResi->desType == Type_DesType_Fixed)
      {
        if (pChain->type == Type_Chain_SmallMol)
        {
          for (int l = 0;l < StringArrayGetCount(&lines);l++)
          {
            if (l<start || l>end)
            {
              fprintf(pMol2, "%s", StringArrayGet(&lines, l));
            }
            else
            {
              char line[MAX_LEN_ONE_LINE_CONTENT + 1];
              strcpy(line, StringArrayGet(&lines, l));
              StringArray words;
              StringArrayCreate(&words);
              StringArraySplitString(&words, line, ' ');
              XYZ xyz = ResidueGetAtom(pResi, l - start)->xyz;
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
        for (int k = 0;k < sequence->desSiteCount;k++)
        {
          DesignSite* pDesignSite = StructureGetDesignSite(pStructure, k);
          RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
          int rotamerIndex = IntArrayGet(&sequence->rotNdxs, k);
          Rotamer* pOptRotamer = RotamerSetGet(pRotamerSet, rotamerIndex);
          if (pDesignSite->chnNdx == i && pDesignSite->resNdx == j)
          {
            RotamerRestore(pOptRotamer, pRotamerSet);
            if (pChain->type == Type_Chain_SmallMol)
            {
              for (int l = 0;l < StringArrayGetCount(&lines);l++)
              {
                if (l<start || l>end)
                {
                  fprintf(pMol2, "%s", StringArrayGet(&lines, l));
                }
                else
                {
                  char line[MAX_LEN_ONE_LINE_CONTENT + 1];
                  strcpy(line, StringArrayGet(&lines, l));
                  StringArray words;
                  StringArrayCreate(&words);
                  StringArraySplitString(&words, line, ' ');
                  XYZ xyz = RotamerGetAtom(pOptRotamer, l - start)->xyz;
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
            RotamerExtract(pOptRotamer);
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


int StructureGetWholeSequence(Structure* pStructure, Sequence* sequence, char* seq)
{
  int index = 0;
  for (int i = 0;i < StructureGetChainCount(pStructure);i++)
  {
    Chain* pChain = StructureGetChain(pStructure, i);
    if (strcmp(DES_CHAINS, ChainGetName(pChain)))continue;
    for (int j = 0;j < ChainGetResidueCount(pChain);j++)
    {
      Residue* pResi = ChainGetResidue(pChain, j);
      if (ResidueGetDesignType(pResi) == Type_DesType_Fixed)
      {
        seq[index++] = AA3ToAA1(ResidueGetName(pResi));
      }
      else
      {
        int siteIdx = StructureGetDesignSiteIndex(pStructure, i, j);
        DesignSite* pSite = StructureGetDesignSite(pStructure, siteIdx);
        RotamerSet* pSet = DesignSiteGetRotamers(pSite);
        int rotIdx = IntArrayGet(&sequence->rotNdxs, siteIdx);
        Rotamer* pRot = RotamerSetGet(pSet, rotIdx);
        seq[index++] = AA3ToAA1(RotamerGetType(pRot));
      }
    }
  }
  seq[index] = '\0';
  return Success;
}


int StructureWriteFirstLigandRotamerIntoMol2(Structure* pStructure, char* mol2file)
{
  int result = Success;
  char errMsg[MAX_LEN_ONE_LINE_CONTENT + 1];

  StringArray lines;
  StringArrayCreate(&lines);
  FILE* pfile = fopen(MOL2, "r");
  char buffer[MAX_LEN_ONE_LINE_CONTENT + 1];
  int start = 100000, end = -1;
  int counter = 0;
  int atomend = 0;
  while (fgets(buffer, MAX_LEN_ONE_LINE_CONTENT, pfile))
  {
    StringArrayAppend(&lines, buffer);
    if (strstr(buffer, "@<TRIPOS>ATOM") != NULL)
    {
      start = counter + 1;
    }
    else if (strstr(buffer, "@<TRIPOS>") != NULL && counter > start && atomend == 0)
    {
      end = counter - 1;
      atomend = 1;
    }
    counter++;
  }
  fclose(pfile);

  FILE* pMol2 = fopen(mol2file, "w");
  if (pMol2 == NULL)
  {
    result = IOError;
    sprintf(errMsg, "in file %s line %d, cannot write to file %s", __FILE__, __LINE__, mol2file);
    TraceError(errMsg, result);
    return result;
  }

  for (int i = 0; i < StructureGetChainCount(pStructure);i++)
  {
    Chain* pChain = StructureGetChain(pStructure, i);
    for (int j = 0;j < ChainGetResidueCount(pChain);j++)
    {
      Residue* pResi = ChainGetResidue(pChain, j);
      if (pResi->desType == Type_DesType_Fixed)
      {
        if (pChain->type == Type_Chain_SmallMol)
        {
          for (int l = 0;l < StringArrayGetCount(&lines);l++)
          {
            if (l<start || l>end)
            {
              fprintf(pMol2, "%s", StringArrayGet(&lines, l));
            }
            else
            {
              char line[MAX_LEN_ONE_LINE_CONTENT + 1];
              strcpy(line, StringArrayGet(&lines, l));
              StringArray words;
              StringArrayCreate(&words);
              StringArraySplitString(&words, line, ' ');
              XYZ xyz = ResidueGetAtom(pResi, l - start)->xyz;
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
        for (int k = 0;k < pStructure->desSiteCount;k++)
        {
          DesignSite* pDesignSite = StructureGetDesignSite(pStructure, k);
          RotamerSet* pRotamerSet = DesignSiteGetRotamers(pDesignSite);
          Rotamer* pOptRotamer = RotamerSetGet(pRotamerSet, 0);
          if (pDesignSite->chnNdx == i && pDesignSite->resNdx == j)
          {
            RotamerRestore(pOptRotamer, pRotamerSet);
            if (pChain->type == Type_Chain_SmallMol)
            {
              for (int l = 0;l < StringArrayGetCount(&lines);l++)
              {
                if (l<start || l>end)
                {
                  fprintf(pMol2, "%s", StringArrayGet(&lines, l));
                }
                else
                {
                  char line[MAX_LEN_ONE_LINE_CONTENT + 1];
                  strcpy(line, StringArrayGet(&lines, l));
                  StringArray words;
                  StringArrayCreate(&words);
                  StringArraySplitString(&words, line, ' ');
                  XYZ xyz = RotamerGetAtom(pOptRotamer, l - start)->xyz;
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
            RotamerExtract(pOptRotamer);
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