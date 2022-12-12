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

#include "RotamerOptimizer.h"
#include "EnergyFunction.h"
#include <string.h>

int ProteinSiteOptimizeRotamer(Structure* pStructure, int chainIndex, int resiIndex)
{
  DesignSite* pDesignSite = StructureFindDesignSite(pStructure, chainIndex, resiIndex);
  if (pDesignSite == NULL) return Success;
  RotamerSet* pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue* pDesign = pDesignSite->pRes;

  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue** ppSurroundingResidues = NULL;
  for (int i = 0; i < pStructure->chainNum; i++)
  {
    Chain* pChainI = pStructure->chains + i;
    for (int j = 0; j < pChainI->residueNum; j++)
    {
      Residue* pResi2 = pChainI->residues + j;
      if (strcmp(pDesign->name, pResi2->name) == 0 && pDesign->posInChain == pResi2->posInChain)
      {
        continue;
      }
      else if (AtomArrayCalcMinDistance(&pDesign->atoms, &pResi2->atoms) < ENERGY_DISTANCE_CUTOFF)
      {
        surroundingResiNum++;
        ppSurroundingResidues = (Residue**)realloc(ppSurroundingResidues, sizeof(Residue*) * surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum - 1] = pResi2;
      }
    }
  }

  //step 2: calculate the energy between the rots of the design site
  double minEnergy = 1000.0;
  for (int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++)
  {
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    Rotamer* pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pDesign);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);
    AminoAcidReferenceEnergy(tempResidue.name, energyTerms);
    EnergyIntraResidue(&tempResidue, energyTerms);
    for (int is = 0; is < surroundingResiNum; is++)
    {
      Residue* pResIS = ppSurroundingResidues[is];
      if (strcmp(ResidueGetChainName(&tempResidue), ResidueGetChainName(pResIS)) == 0)
      {
        if (tempResidue.posInChain == pResIS->posInChain - 1)
        {
          EnergyResidueAndNextResidue(&tempResidue, pResIS, energyTerms);
        }
        else if (tempResidue.posInChain == pResIS->posInChain + 1)
        {
          EnergyResidueAndNextResidue(pResIS, &tempResidue, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueSameChain(&tempResidue, pResIS, energyTerms);
        }
      }
      else
      {
        if (ChainGetType(StructureFindChainByName(pStructure, ResidueGetChainName(pResIS))) == Type_Chain_SmallMol)
        {
          EnergyResidueAndLigandResidue(&tempResidue, pResIS, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueDiffChain(&tempResidue, pResIS, energyTerms);
        }
      }
    }

    EnergyTermWeighting(energyTerms);
    if (energyTerms[0] < minEnergy)
    {
      minEnergy = energyTerms[0];
      ResidueCopy(pDesign, &tempResidue);
    }
    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  if (ppSurroundingResidues != NULL)
  {
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }
  return Success;
}


int ProteinSiteOptimizeRotamerHBondEnergy(Structure* pStructure, int chainIndex, int resiIndex)
{
  DesignSite* pDesignSite = StructureFindDesignSite(pStructure, chainIndex, resiIndex);
  if (pDesignSite == NULL) return Success;
  RotamerSet* pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue* pDesign = pDesignSite->pRes;

  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue** ppSurroundingResidues = NULL;
  for (int i = 0; i < pStructure->chainNum; i++)
  {
    Chain* pChainI = pStructure->chains + i;
    for (int j = 0; j < pChainI->residueNum; j++)
    {
      Residue* pResi2 = pChainI->residues + j;
      if (strcmp(pDesign->name, pResi2->name) == 0 && pDesign->posInChain == pResi2->posInChain)
      {
        continue;
      }
      else if (AtomArrayCalcMinDistance(&pDesign->atoms, &pResi2->atoms) < ENERGY_DISTANCE_CUTOFF)
      {
        surroundingResiNum++;
        ppSurroundingResidues = (Residue**)realloc(ppSurroundingResidues, sizeof(Residue*) * surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum - 1] = pResi2;
      }
    }
  }

  //step 2: calculate the energy between the rots of the design site
  double minEnergy = 1000.0;
  for (int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++)
  {
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    Rotamer* pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pDesign);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);
    AminoAcidReferenceEnergy(tempResidue.name, energyTerms);
    EnergyIntraResidue(&tempResidue, energyTerms);
    for (int is = 0; is < surroundingResiNum; is++)
    {
      Residue* pResIS = ppSurroundingResidues[is];
      if (strcmp(ResidueGetChainName(&tempResidue), ResidueGetChainName(pResIS)) == 0)
      {
        if (tempResidue.posInChain == pResIS->posInChain - 1)
        {
          EnergyResidueAndNextResidue(&tempResidue, pResIS, energyTerms);
        }
        else if (tempResidue.posInChain == pResIS->posInChain + 1)
        {
          EnergyResidueAndNextResidue(pResIS, &tempResidue, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueSameChain(&tempResidue, pResIS, energyTerms);
        }
      }
      else
      {
        if (ChainGetType(StructureFindChainByName(pStructure, ResidueGetChainName(pResIS))) == Type_Chain_SmallMol)
        {
          EnergyResidueAndLigandResidue(&tempResidue, pResIS, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueDiffChain(&tempResidue, pResIS, energyTerms);
        }
      }
    }

    EnergyTermWeighting(energyTerms);
    //only consider the non-intra residue hbond energy
    double hbenergy;
    for (int i = 41;i <= 49;i++) hbenergy += energyTerms[i];
    for (int i = 61;i <= 69;i++) hbenergy += energyTerms[i];
    for (int i = 81;i <= 89;i++) hbenergy += energyTerms[i];

    if (hbenergy < minEnergy)
    {
      minEnergy = hbenergy;
      ResidueCopy(pDesign, &tempResidue);
    }
    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  if (ppSurroundingResidues != NULL)
  {
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }

  return Success;
}

int ProteinSiteCalcRotamersEnergy(Structure* pStructure, AAppTable* pAAppTable, RamaTable* pRama, int chainIndex, int resiIndex, FILE* fp)
{
  DesignSite* pDesignSite = StructureFindDesignSite(pStructure, chainIndex, resiIndex);
  if (pDesignSite == NULL) return Success;
  RotamerSet* pRotSet = DesignSiteGetRotamers(pDesignSite);
  Residue* pResidue = pDesignSite->pRes;

  //step 1: find out residues within 5 angstroms to the design site of interest;
  int surroundingResiNum = 0;
  Residue** ppSurroundingResidues = NULL;
  for (int i = 0; i < pStructure->chainNum; i++)
  {
    Chain* pChainI = pStructure->chains + i;
    for (int j = 0; j < pChainI->residueNum; j++)
    {
      Residue* pResi2 = pChainI->residues + j;
      if (strcmp(pResidue->name, pResi2->name) == 0 && pResidue->posInChain == pResi2->posInChain)
      {
        continue;
      }
      else if (AtomArrayCalcMinDistance(&pResidue->atoms, &pResi2->atoms) < ENERGY_DISTANCE_CUTOFF)
      {
        surroundingResiNum++;
        ppSurroundingResidues = (Residue**)realloc(ppSurroundingResidues, sizeof(Residue*) * surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum - 1] = pResi2;
      }
    }
  }

  // step 2: calculate the energy between the rots of the design site
  for (int ir = 0; ir < RotamerSetGetCount(pRotSet); ir++)
  {
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    Rotamer* pRotIR = RotamerSetGet(pRotSet, ir);
    RotamerRestore(pRotIR, pRotSet);

    Residue tempResidue;
    ResidueCreate(&tempResidue);
    ResidueCopy(&tempResidue, pResidue);
    ResidueSetName(&tempResidue, pRotIR->type);
    AtomArrayCopy(&tempResidue.atoms, &pRotIR->atoms);
    BondSetCopy(&tempResidue.bonds, &pRotIR->bonds);

    //intra energy, add the statistical energy
    AminoAcidReferenceEnergy(tempResidue.name, energyTerms);
    EnergyIntraResidue(&tempResidue, energyTerms);
    RotamerPropensityAndRamachandranEnergy(pRotIR, pResidue, pAAppTable, pRama, energyTerms);
    RotamerDunbrackEnergy(pRotIR, energyTerms);
    for (int is = 0; is < surroundingResiNum; is++)
    {
      Residue* pResIS = ppSurroundingResidues[is];
      if (strcmp(ResidueGetChainName(&tempResidue), ResidueGetChainName(pResIS)) == 0)
      {
        if (tempResidue.posInChain == pResIS->posInChain - 1)
        {
          EnergyResidueAndNextResidue(&tempResidue, pResIS, energyTerms);
        }
        else if (tempResidue.posInChain == pResIS->posInChain + 1)
        {
          EnergyResidueAndNextResidue(pResIS, &tempResidue, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueSameChain(&tempResidue, pResIS, energyTerms);
        }
      }
      else
      {
        if (ChainGetType(StructureFindChainByName(pStructure, ResidueGetName(pResIS))) == Type_Chain_SmallMol)
        {
          EnergyResidueAndLigandResidue(&tempResidue, pResIS, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueDiffChain(&tempResidue, pResIS, energyTerms);
        }
      }
    }

    //disable the weighting
    EnergyTermWeighting(energyTerms);
    if (RotamerAndResidueInSameType(pRotIR, pResidue))
    {
      double rmsd = RotamerAndResidueSidechainRMSD(pRotIR, pResidue);
      if (rmsd < 0.005)
      {//this should be a crystal rotamer
        fprintf(fp, "POS %3d XAL", RotamerGetPosInChain(pRotIR));
        //if the energy for monomer & PPI is too large, skip the crystal rotamer
        BOOL energyTooBig = FALSE;
        for (int index = 1;index <= 70;index++)
        {
          if (energyTerms[index] > 30.0 || energyTerms[index] < -30.0)
          {
            energyTooBig = TRUE;
            break;
          }
        }
        if (energyTooBig)
        {
          char errMsg[MAX_LEN_ERR_MSG + 1];
          sprintf(errMsg, "in file %s line %d, the energy value of native rotamer on residue %s%d%s is too high", __FILE__, __LINE__, ResidueGetChainName(&tempResidue), ResidueGetPosInChain(&tempResidue), ResidueGetName(&tempResidue));
          TraceError(errMsg, Warning);
        }

      }
      else
      {
        fprintf(fp, "POS %3d NAT", RotamerGetPosInChain(pRotIR));
      }
    }
    else
    {
      fprintf(fp, "POS %3d %s", RotamerGetPosInChain(pRotIR), RotamerGetType(pRotIR));
    }
    //20 reference energy
    for (int index = 1; index <= 20; index++)
    {
      fprintf(fp, " %8.2f", energyTerms[index]);
    }
    //intraR vdw->desol
    for (int index = 21; index <= 25; index++)
    {
      fprintf(fp, " %8.2f", energyTerms[index]);
    }
    //intraR hbonds
    for (int index = 26; index <= 28; index++)
    {
      fprintf(fp, " %8.2f", energyTerms[index]);
    }
    //intraR aa propensity, rama and dunbrack
    for (int index = 91; index <= 93; index++)
    {
      fprintf(fp, " %8.2f", energyTerms[index]);
    }
    //interS vdw->desol
    for (int index = 31; index <= 36; index++)
    {
      fprintf(fp, " %8.2f", energyTerms[index]);
    }
    //interS hbonds
    for (int index = 41; index <= 49; index++)
    {
      fprintf(fp, " %8.2f", energyTerms[index]);
    }
    //interD vdw->desol
    for (int index = 51; index <= 56; index++)
    {
      fprintf(fp, " %8.2f", energyTerms[index]);
    }
    //interD hbonds
    for (int index = 61; index <= 69; index++)
    {
      fprintf(fp, " %8.2f", energyTerms[index]);
    }
    //prot-lig vdw->desol
    for (int index = 71; index <= 75; index++)
    {
      fprintf(fp, " %8.2f", energyTerms[index]);
    }
    //prot-lig hbonds
    for (int index = 81; index <= 86; index++)
    {
      fprintf(fp, " %8.2f", energyTerms[index]);
    }
    fprintf(fp, "\n");

    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResidue);
  }

  if (ppSurroundingResidues != NULL)
  {
    free(ppSurroundingResidues);
    ppSurroundingResidues = NULL;
  }

  return Success;
}


int ProteinSiteOptimizeRotamerWithBBdepRotLib(Structure* pStructure, int chainIndex, int resiIndex, BBdepRotamerLib* pBBdepRotLib)
{
  DesignSite* pSite = StructureFindDesignSite(pStructure, chainIndex, resiIndex);
  if (pSite == NULL) return Success;
  RotamerSet* pSet = DesignSiteGetRotamers(pSite);
  Residue* pResi = pSite->pRes;

  int surroundingResiNum = 0;
  Residue** ppSurroundingResidues = NULL;
  for (int i = 0; i < pStructure->chainNum; i++)
  {
    Chain* pChainI = pStructure->chains + i;
    for (int j = 0; j < pChainI->residueNum; j++)
    {
      Residue* pResi2 = pChainI->residues + j;
      if (strcmp(pResi->name, pResi2->name) == 0 && pResi->posInChain == pResi2->posInChain) continue;
      else if (AtomArrayCalcMinDistance(&pResi->atoms, &pResi2->atoms) < ENERGY_DISTANCE_CUTOFF)
      {
        surroundingResiNum++;
        ppSurroundingResidues = (Residue**)realloc(ppSurroundingResidues, sizeof(Residue*) * surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum - 1] = pResi2;
      }
    }
  }

  double minEnergy = 1000.0;
  for (int ir = 0; ir < RotamerSetGetCount(pSet); ir++)
  {
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    Rotamer* pRotIR = RotamerSetGet(pSet, ir);
    RotamerRestore(pRotIR, pSet);

    Residue tempResi;
    ResidueCreate(&tempResi);
    ResidueCopy(&tempResi, pResi);
    ResidueSetName(&tempResi, RotamerGetType(pRotIR));
    AtomArrayCopy(ResidueGetAllAtoms(&tempResi), &pRotIR->atoms);
    BondSetCopy(ResidueGetBonds(&tempResi), RotamerGetBonds(pRotIR));
    AminoAcidReferenceEnergy(ResidueGetName(&tempResi), energyTerms);
    EnergyIntraResidue(&tempResi, energyTerms);
    RotamerDunbrackEnergy(pRotIR, energyTerms);
    ResidueSetDunbrack(&tempResi, RotamerGetDunbrack(pRotIR));
    DoubleArrayCopy(&tempResi.Xs, &pRotIR->Xs);
    for (int is = 0; is < surroundingResiNum; is++)
    {
      Residue* pResIS = ppSurroundingResidues[is];
      if (!strcmp(ResidueGetChainName(&tempResi), ResidueGetChainName(pResIS)))
      {
        if (tempResi.posInChain == pResIS->posInChain - 1)
        {
          EnergyResidueAndNextResidue(&tempResi, pResIS, energyTerms);
        }
        else if (tempResi.posInChain == pResIS->posInChain + 1)
        {
          EnergyResidueAndNextResidue(pResIS, &tempResi, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueSameChain(&tempResi, pResIS, energyTerms);
        }
      }
      else
      {
        if (ChainGetType(StructureFindChainByName(pStructure, ResidueGetChainName(pResIS))) == Type_Chain_SmallMol)
        {
          EnergyResidueAndLigandResidue(&tempResi, pResIS, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueDiffChain(&tempResi, pResIS, energyTerms);
        }
      }
    }

    EnergyTermWeighting(energyTerms);
    if (energyTerms[0] < minEnergy)
    {
      minEnergy = energyTerms[0];
      ResidueCopy(pResi, &tempResi);
    }

    RotamerExtract(pRotIR);
    ResidueDestroy(&tempResi);
  }
  free(ppSurroundingResidues);
  ppSurroundingResidues = NULL;

  return Success;
}
