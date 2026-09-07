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

#include "RotamerOptimizer.h"
#include <string.h>

int ProteinSiteOptimizeRotamer(Structure* pStruct, int chnNdx, int resNdx)
{
  DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pSite == NULL)
  {
    return Success;
  }
  Rotamers* pSet = DesignSiteGetRotamers(pSite);
  Residue* pRes = pSite->pRes;

  int numCloseRes = 0;
  Residue** ppCloseRes = NULL;
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain* pChn = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChn); j++)
    {
      Residue* pRes2 = ChainGetResidue(pChn, j);
      if ( strcmp(ResidueGetName(pRes), ResidueGetName(pRes2)) == 0 && 
           ResidueGetPosInChain(pRes) == ResidueGetPosInChain(pRes2) )
      {
        continue;
      }
      else if (AtomArrayCalcMinDistance(&pRes->atoms, &pRes2->atoms) < ENERGY_DISTANCE_CUTOFF)
      {
        numCloseRes++;
        ppCloseRes = (Residue**)realloc(ppCloseRes, sizeof(Residue*) * numCloseRes);
        ppCloseRes[numCloseRes - 1] = pRes2;
      }
    }
  }

  double minEnergy = 1000.0;
  for (int ir = 0; ir < RotamersGetNum(pSet); ir++)
  {
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    Rotamer* pRotIR = RotamersGet(pSet, ir);
    RotamerRestore(pRotIR, pSet);

    Residue res;
    ResidueCreate(&res);
    ResidueCopy(&res, pRes);
    ResidueSetName(&res, RotamerGetType(pRotIR));
    AtomArrayCopy(ResidueGetAllAtoms(&res), &pRotIR->atoms);
    BondsCopy(ResidueGetBonds(&res), RotamerGetBonds(pRotIR));
    AminoAcidReferenceEnergy(ResidueGetName(&res), energyTerms);
    EnergyIntraResidue(&res, energyTerms);
    for (int is = 0; is < numCloseRes; is++)
    {
      Residue* pResIS = ppCloseRes[is];
      if (strcmp(ResidueGetChainName(&res), ResidueGetChainName(pResIS)) == 0)
      {
        if (ResidueGetPosInChain(&res) == ResidueGetPosInChain(pResIS) - 1)
        {
          EnergyResidueAndNextResidue(&res, pResIS, energyTerms);
        }
        else if (ResidueGetPosInChain(&res) == ResidueGetPosInChain(pResIS) + 1)
        {
          EnergyResidueAndNextResidue(pResIS, &res, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueSameChain(&res, pResIS, energyTerms);
        }
      }
      else
      {
        if (ChainGetType(StructureFindChainByName(pStruct, ResidueGetChainName(pResIS))) == Type_Chain_SmallMol)
        {
          EnergyResidueAndLigandResidue(&res, pResIS, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueDiffChain(&res, pResIS, energyTerms);
        }
      }
    }

    EnergyTermWeighting(energyTerms);
    if (energyTerms[0] < minEnergy)
    {
      minEnergy = energyTerms[0];
      ResidueCopy(pRes, &res);
    }
    RotamerExtract(pRotIR);
    ResidueDestroy(&res);
  }

  if (ppCloseRes != NULL)
  {
    free(ppCloseRes);
    ppCloseRes = NULL;
  }
  return Success;
}


int ProteinSiteOptimizeRotamerHBondEnergy(Structure* pStruct, int chnNdx, int resNdx)
{
  DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pSite == NULL)
  {
    return Success;
  }
  Rotamers* pSet = DesignSiteGetRotamers(pSite);
  Residue* pRes = pSite->pRes;

  int numCloseRes = 0;
  Residue** ppCloseRes = NULL;
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain* pChn = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChn); j++)
    {
      Residue* pRes2 = ChainGetResidue(pChn, j);
      if (strcmp(ResidueGetName(pRes), ResidueGetName(pRes2)) == 0 &&
        ResidueGetPosInChain(pRes) == ResidueGetPosInChain(pRes2))
      {
        continue;
      }
      else if (AtomArrayCalcMinDistance(&pRes->atoms, &pRes2->atoms) < ENERGY_DISTANCE_CUTOFF)
      {
        numCloseRes++;
        ppCloseRes = (Residue**)realloc(ppCloseRes, sizeof(Residue*) * numCloseRes);
        ppCloseRes[numCloseRes - 1] = pRes2;
      }
    }
  }

  double minEnergy = 1000.0;
  for (int ir = 0; ir < RotamersGetNum(pSet); ir++)
  {
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    Rotamer* pRotIR = RotamersGet(pSet, ir);
    RotamerRestore(pRotIR, pSet);

    Residue res;
    ResidueCreate(&res);
    ResidueCopy(&res, pRes);
    ResidueSetName(&res, RotamerGetType(pRotIR));
    AtomArrayCopy(ResidueGetAllAtoms(&res), &pRotIR->atoms);
    BondsCopy(ResidueGetBonds(&res), RotamerGetBonds(pRotIR));
    AminoAcidReferenceEnergy(ResidueGetName(&res), energyTerms);
    EnergyIntraResidue(&res, energyTerms);
    for (int is = 0; is < numCloseRes; is++)
    {
      Residue* pResIS = ppCloseRes[is];
      if (strcmp(ResidueGetChainName(&res), ResidueGetChainName(pResIS)) == 0)
      {
        if (ResidueGetPosInChain(&res) == ResidueGetPosInChain(pResIS) - 1)
        {
          EnergyResidueAndNextResidue(&res, pResIS, energyTerms);
        }
        else if (ResidueGetPosInChain(&res) == ResidueGetPosInChain(pResIS) + 1)
        {
          EnergyResidueAndNextResidue(pResIS, &res, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueSameChain(&res, pResIS, energyTerms);
        }
      }
      else
      {
        if (ChainGetType(StructureFindChainByName(pStruct, ResidueGetChainName(pResIS))) == Type_Chain_SmallMol)
        {
          EnergyResidueAndLigandResidue(&res, pResIS, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueDiffChain(&res, pResIS, energyTerms);
        }
      }
    }

    EnergyTermWeighting(energyTerms);

    // only consider the non-intra residue hbond energy
    double hbenergy;
    for (int i = 41; i <= 49; i++) hbenergy += energyTerms[i];
    for (int i = 61; i <= 69; i++) hbenergy += energyTerms[i];
    for (int i = 81; i <= 89; i++) hbenergy += energyTerms[i];

    if (hbenergy < minEnergy)
    {
      minEnergy = hbenergy;
      ResidueCopy(pRes, &res);
    }
    RotamerExtract(pRotIR);
    ResidueDestroy(&res);
  }

  if (ppCloseRes != NULL)
  {
    free(ppCloseRes);
    ppCloseRes = NULL;
  }

  return Success;
}

int ProteinSiteCalcRotamersEnergy(Structure* pStruct, AAppTable* pAAppTable, RamaTable* pRama, int chnNdx, int resNdx, FILE* pOut)
{
  DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pSite == NULL)
  {
    return Success;
  }
  Rotamers* pSet = DesignSiteGetRotamers(pSite);
  Residue* pRes = pSite->pRes;

  int numCloseRes = 0;
  Residue** ppCloseRes = NULL;
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain* pChn = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChn); j++)
    {
      Residue* pRes2 = ChainGetResidue(pChn, j);
      if (strcmp(ResidueGetName(pRes), ResidueGetName(pRes2)) == 0 &&
        ResidueGetPosInChain(pRes) == ResidueGetPosInChain(pRes2))
      {
        continue;
      }
      else if (AtomArrayCalcMinDistance(&pRes->atoms, &pRes2->atoms) < ENERGY_DISTANCE_CUTOFF)
      {
        numCloseRes++;
        ppCloseRes = (Residue**)realloc(ppCloseRes, sizeof(Residue*) * numCloseRes);
        ppCloseRes[numCloseRes - 1] = pRes2;
      }
    }
  }

  for (int ir = 0; ir < RotamersGetNum(pSet); ir++)
  {
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    Rotamer* pRotIR = RotamersGet(pSet, ir);
    RotamerRestore(pRotIR, pSet);

    Residue res;
    ResidueCreate(&res);
    ResidueCopy(&res, pRes);
    ResidueSetName(&res, RotamerGetType(pRotIR));
    AtomArrayCopy(ResidueGetAllAtoms(&res), &pRotIR->atoms);
    BondsCopy(ResidueGetBonds(&res), RotamerGetBonds(pRotIR));
    AminoAcidReferenceEnergy(ResidueGetName(&res), energyTerms);
    EnergyIntraResidue(&res, energyTerms);
    RotamerPropensityAndRamachandranEnergy(pRotIR, pRes, pAAppTable, pRama, energyTerms);
    RotamerDunbrackEnergy(pRotIR, energyTerms);
    for (int is = 0; is < numCloseRes; is++)
    {
      Residue* pResIS = ppCloseRes[is];
      if (strcmp(ResidueGetChainName(&res), ResidueGetChainName(pResIS)) == 0)
      {
        if (ResidueGetPosInChain(&res) == ResidueGetPosInChain(pResIS) - 1)
        {
          EnergyResidueAndNextResidue(&res, pResIS, energyTerms);
        }
        else if (ResidueGetPosInChain(&res) == ResidueGetPosInChain(pResIS) + 1)
        {
          EnergyResidueAndNextResidue(pResIS, &res, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueSameChain(&res, pResIS, energyTerms);
        }
      }
      else
      {
        if (ChainGetType(StructureFindChainByName(pStruct, ResidueGetName(pResIS))) == Type_Chain_SmallMol)
        {
          EnergyResidueAndLigandResidue(&res, pResIS, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueDiffChain(&res, pResIS, energyTerms);
        }
      }
    }

    EnergyTermWeighting(energyTerms);
    if (RotamerAndResidueInSameType(pRotIR, pRes))
    {
      double rmsd = RotamerAndResidueSidechainRMSD(pRotIR, pRes);
      if (rmsd < 0.005)
      { // must be the crystal rotamer
        fprintf(pOut, "POS %3d XAL", RotamerGetPosInChain(pRotIR));
        // if the energy for monomer & PPI is too large, skip the crystal rotamer
        BOOL energyTooHigh = FALSE;
        for (int index = 1; index <= 70; index++)
        {
          if (energyTerms[index] > 30.0 || energyTerms[index] < -30.0)
          {
            energyTooHigh = TRUE;
            break;
          }
        }
        if (energyTooHigh == TRUE)
        {
          char errMsg[MAX_LEN_ERR_MSG+1];
          sprintf(errMsg, "in file %s line %d, the energy value of native rotamer on residue %s%d%s is too high", 
            __FILE__, __LINE__, ResidueGetChainName(&res), ResidueGetPosInChain(&res), ResidueGetName(&res));
          TraceError(errMsg, Warning);
        }

      }
      else
      {
        fprintf(pOut, "POS %3d NAT", RotamerGetPosInChain(pRotIR));
      }
    }
    else
    {
      fprintf(pOut, "POS %3d %s", RotamerGetPosInChain(pRotIR), RotamerGetType(pRotIR));
    }
    //20 reference energy
    for (int index = 1; index <= 20; index++)
    {
      fprintf(pOut, " %8.2f", energyTerms[index]);
    }
    //intraR vdw->desol
    for (int index = 21; index <= 25; index++)
    {
      fprintf(pOut, " %8.2f", energyTerms[index]);
    }
    //intraR hbonds
    for (int index = 26; index <= 28; index++)
    {
      fprintf(pOut, " %8.2f", energyTerms[index]);
    }
    //intraR aa propensity, rama and dunbrack
    for (int index = 91; index <= 93; index++)
    {
      fprintf(pOut, " %8.2f", energyTerms[index]);
    }
    //interS vdw->desol
    for (int index = 31; index <= 36; index++)
    {
      fprintf(pOut, " %8.2f", energyTerms[index]);
    }
    //interS hbonds
    for (int index = 41; index <= 49; index++)
    {
      fprintf(pOut, " %8.2f", energyTerms[index]);
    }
    //interD vdw->desol
    for (int index = 51; index <= 56; index++)
    {
      fprintf(pOut, " %8.2f", energyTerms[index]);
    }
    //interD hbonds
    for (int index = 61; index <= 69; index++)
    {
      fprintf(pOut, " %8.2f", energyTerms[index]);
    }
    //prot-lig vdw->desol
    for (int index = 71; index <= 75; index++)
    {
      fprintf(pOut, " %8.2f", energyTerms[index]);
    }
    //prot-lig hbonds
    for (int index = 81; index <= 86; index++)
    {
      fprintf(pOut, " %8.2f", energyTerms[index]);
    }
    fprintf(pOut, "\n");

    RotamerExtract(pRotIR);
    ResidueDestroy(&res);
  }

  if (ppCloseRes != NULL)
  {
    free(ppCloseRes);
    ppCloseRes = NULL;
  }

  return Success;
}


int ProteinSiteOptimizeRotamer(Structure* pStruct, int chnNdx, int resNdx, BBdepRotLib* rotlib)
{
  DesignSite* pSite = StructureFindDesignSite(pStruct, chnNdx, resNdx);
  if (pSite == NULL)
  {
    return Success;
  }
  Rotamers* pSet = DesignSiteGetRotamers(pSite);
  Residue* pRes = pSite->pRes;

  int numCloseRes = 0;
  Residue** ppCloseRes = NULL;
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain* pChn = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChn); j++)
    {
      Residue* pRes2 = ChainGetResidue(pChn, j);
      if (strcmp(ResidueGetName(pRes), ResidueGetName(pRes2)) == 0 &&
        ResidueGetPosInChain(pRes) == ResidueGetPosInChain(pRes2))
      {
        continue;
      }
      else if (AtomArrayCalcMinDistance(&pRes->atoms, &pRes2->atoms) < ENERGY_DISTANCE_CUTOFF)
      {
        numCloseRes++;
        ppCloseRes = (Residue**)realloc(ppCloseRes, sizeof(Residue*) * numCloseRes);
        ppCloseRes[numCloseRes - 1] = pRes2;
      }
    }
  }

  double minEnergy = 1000.0;
  for (int ir = 0; ir < RotamersGetNum(pSet); ir++)
  {
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    Rotamer* pRotIR = RotamersGet(pSet, ir);
    RotamerRestore(pRotIR, pSet);

    Residue res;
    ResidueCreate(&res);
    ResidueCopy(&res, pRes);
    ResidueSetName(&res, RotamerGetType(pRotIR));
    AtomArrayCopy(ResidueGetAllAtoms(&res), &pRotIR->atoms);
    BondsCopy(ResidueGetBonds(&res), RotamerGetBonds(pRotIR));
    AminoAcidReferenceEnergy(ResidueGetName(&res), energyTerms);
    EnergyIntraResidue(&res, energyTerms);
    RotamerDunbrackEnergy(pRotIR, energyTerms);
    ResidueSetDunbrack(&res, RotamerGetDunbrack(pRotIR));
    DoubleArrayCopy(&res.Xs, &pRotIR->Xs);
    for (int is = 0; is < numCloseRes; is++)
    {
      Residue* pResIS = ppCloseRes[is];
      if (!strcmp(ResidueGetChainName(&res), ResidueGetChainName(pResIS)))
      {
        if (ResidueGetPosInChain(&res) == ResidueGetPosInChain(pResIS) - 1)
        {
          EnergyResidueAndNextResidue(&res, pResIS, energyTerms);
        }
        else if (ResidueGetPosInChain(&res) == ResidueGetPosInChain(pResIS) + 1)
        {
          EnergyResidueAndNextResidue(pResIS, &res, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueSameChain(&res, pResIS, energyTerms);
        }
      }
      else
      {
        if (ChainGetType(StructureFindChainByName(pStruct, ResidueGetChainName(pResIS))) == Type_Chain_SmallMol)
        {
          EnergyResidueAndLigandResidue(&res, pResIS, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueDiffChain(&res, pResIS, energyTerms);
        }
      }
    }

    EnergyTermWeighting(energyTerms);
    if (energyTerms[0] < minEnergy)
    {
      minEnergy = energyTerms[0];
      ResidueCopy(pRes, &res);
    }
    RotamerExtract(pRotIR);
    ResidueDestroy(&res);
  }

  if (ppCloseRes != NULL)
  {
    free(ppCloseRes);
    ppCloseRes = NULL;
  }

  return Success;
}

