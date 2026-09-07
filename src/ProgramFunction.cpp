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
#pragma warning(disable : 6031)
#pragma warning(disable : 4476)

#include "ProgramFunction.h"
#include "RotamerBuilder.h"
#include "RotamerOptimizer.h"
//#include "Evolution.h"
#include <string.h>
#include <ctype.h>

extern BOOL FLAG_DESIGN_MON;
extern BOOL FLAG_DESIGN_PPI;
extern BOOL FLAG_DESIGN_ENZ;
extern BOOL FLAG_DESIGN_PLI;
extern BOOL FLAG_USE_INPUT_SC;
extern BOOL FLAG_ROTATE_HYDROXYL;
extern double TORSION_DEVIATION_THRESHOLD;
extern double ROTLIB_LOW_PROB_CUTOFF;
extern char PROGRAM_NAME[MAX_LEN_FILE_NAME + 1];
extern char PROGRAM_VERSION[MAX_LEN_FILE_NAME + 1];
extern int CUT_NUM_CB_CORE;
extern int CUT_NUM_CB_SURF;
extern double CUT_PPI_DIST_SHELL1;
extern int MAX_NUM_OF_RUNS;
extern char REFERENCE_RESIDUES[MAX_LEN_ONE_LINE_CONTENT + 1];
extern double DISTANCE_TO_REFERENCE_RESIDUES;
extern char EXCLUDE_RESIDUES[MAX_LEN_ONE_LINE_CONTENT + 1];

int PrintProgramHelp()
{
  printf(
      "NAME \n\n"
      "  %s - a light weight program for computational protein design and modeling \n\n"
      "SYNOPSIS\n\n"
      "  %s [OPTIONS] ...\n\n"
      "DESCRIPTION\n\n"
      "  -v, --version \n"
      "    print out program version \n\n"
      "  -h, -?, --help \n"
      "    print out help information \n\n"
      "  --command <arg> \n"
      "    supported commands are: \n"
      "    1. Energy calculation commands, including: \n"
      "      'ComputeStability', 'ComputeBinding', 'ComputeEvoScore', 'ComputeResEnergy', 'ComputeResPairEnergy',  'ComputeRotamersEnergy'. \n"
      "    2. Structure modeling commands, including: \n"
      "      'RepairStructure', 'Minimize', 'BuildMutant', 'AddPolarHydrogen', 'OptimizeHydrogen', and 'ProteinDesign' used with flag \n"
      "      '--wildtype_only' for protein sidechain repacking. \n"
      "    3. Protein Design commands: \n"
      "      'ProteinDesign' used with one of the four flags (--design_mon, --design_ppi, --design_pli, or --design_enz) to design new \n"
      "      protein sequences or variants (for protein engineering) \n"
      "    4. Ligand pose generatiion and analysis commands: "
      "      'MakeLigParamAndTopo', 'MakeLigPoses', 'AnalyzeLigPoses', 'ScreenLigPoses' \n"
      "    5. Protein structure checking and assessment commands, including: \n"
      "      'CheckClash0', 'CheckClash1', 'CheckClash2', 'CheckResiInRotLib', 'CompareSideChain', 'FindCoreResidue', 'FindInterfaceResidue', \n"
      "      'FindIntermediateResidue', 'FindSurfaceResidue', 'GetPhiPsi', 'GetResiMinRMSD', 'GetResiBfactor', 'ShowResiComposition' \n"
      "    6. Sequence-based feature prediction commands: \n"
      "      'PhiPsiPred', 'SSPred', 'SAPred' \n"
      "    7. Other commands: \n"
      "      'OptimizeWeight' \n"
      "                               \n"
      "  commands can be used in conjunction with one or more options below: \n"
      "    --pdb  <arg> \n"
      "                            <arg> specifies an input pdb file. \n"
      "    --split_chains  <arg>\n"
      "                            <arg> specifies the chain splitting manner for a system with no less than three chains (e.g. 'AB,C', 'AB,CDE'). \n"
      "    --mutant_file  <arg> \n"
      "                            <arg> specifies an input file that records mutants. \n"
      "    --use_input_sc \n"
      "                            Flag for using the input protein's amino-acid sidechains as rotamers. \n"
      "    --ppi_shell1  <arg> \n"
      "                            <arg> specifies the distance cutoff for the first shell of PPI (default: 5 Ang). \n"
      "    --ppi_shell2  <arg> \n"
      "                            <arg> specifies the distance cutoff for the second shell of PPI (default: 8 Ang). \n"
      "    --pli_shell1  <arg> \n"
      "                            <arg> specifies the distance cutoff for the first shell of protein-ligand interaction (default: 5 Ang). \n"
      "    --pli_shell2  <arg> \n"
      "                            <arg> specifies the distance cutoff for the second shell of protein-ligand interaction (default: 8 Ang). \n"
      "    --evolution \n"
      "                            This flag turns on the evolutionary score for protein design. \n"
      "    --design_mon \n"
      "                            This flag turns on monomer protein design. \n"
      "    --design_ppi \n"
      "                            This flag turns on protein-protein interaction design. \n"
      "    --design_pli \n"
      "                            This flag turns on protein-ligand interaction design. \n"
      "    --design_enz \n"
      "                            This flag turns on small molecule ligand-involved enzyme design. \n"
      "    --design_chains  <arg> \n"
      "                            <arg> specifies protein chains to be redesigned (e.g. 'A', 'AB').\n"
      "    --wildtype_only \n"
      "                            This flag specifies that rotamers should be generated for native amino acid types. \n"
      "    --interface_only\n"
      "                            This flag specifies that only interface residues will be designed. \n"
      "    --resfile  <arg> \n"
      "                            <arg> specifies an input file that set residues for protein design. \n"
      "    --seed_from_wildtype\n"
      "                            turn on the flag for using native sequence as a starting point for protein design\n"
      "    --exclude_cys_rots\n"
      "                            exclude the Cystenine rotamers for protein design\n"
      "    --mol2  <arg> \n"
      "                            <arg> specifies an input mol2 file (Note: atomic charge should be already included). \n"
      "    --rotate_hydroxyl  <arg> \n"
      "                            <arg> = 'true'/'false'. This option turns on/off the flag for rotating the hydroxyl groups of Ser, \n"
      "                            Thr, and Tyr to generate extra rotamers to improve sampling (default: true). \n"
      "    --xdeviation  <arg> \n"
      "                            <arg> specifies a value threshold (in degrees) to determine the similarity of two torsional angles. \n"
      "    --wread  <arg> \n"
      "                            <arg> specifies the file for reading energy weights \n"
      "    --pdblist  <arg> \n"
      "                            <arg> specifies a file that records a list of protein IDs, one in each line. \n"
      "    --compare_pdb  <arg> \n"
      "                            <arg> specifies another input pdb file, whose model will be compared to the structure specified with '--pdb'. \n"
      "                            This option is used with command 'CompareSideChain' to compare two proteins' sidechain similarity. \n"
      "    --clash_ratio  <arg> \n"
      "                            <arg> specifies a threshold value for command 'CheckClash0', 'CheckClash1', and 'CheckClash2' (default: 0.75). \n"
      "    --ntraj  <arg> \n"
      "                            <arg> specifies the number of protein design trajectories (default: 1). \n"
      "    --low_prob_rot_cutoff  <arg> \n"
      "                            <arg> specifies a threshold value for excluding low-probability rotamers. \n"
      "                            Suggested values: 0.01~0.03 (default: 0.03). \n"
      "    --seq_file  <arg> \n"
      "                            <arg> specifies an input protein fasta file. Protein sequence should be written in a single line. \n"
      "    --all_evo_terms \n"
      "                            turn on the flag for using all-term evolution (including PSSM, SA, SS, and psi-phi) for protein design. \n"
      "    --wprof  <arg> \n"
      "                            <arg> is a float value (default: 1.0)\n"

      "    --ncut_cb_core  <arg>\n"
      "                            <arg> is an integer cutoff value for the number of CB atoms for a CORE (buried) residue (default: 20)\n"
      "    --ncut_cb_surf  <arg>\n"
      "                            <arg> is an integer cutoff value for the number of CB atoms for a SURFACE (exposed) residue (default: 15)\n"

      "    --init_3atoms  <arg>\n"
      "                            <arg> is the name of three atoms divided by commas (e.g. C1,C2,C3)\n"
      "                            this option is only used with the command GenLigParamAndTopo\n"
      "    --read_lig_poses  <arg>\n"
      "                            <arg> is the PDB file from which ligand rotameric poses will be read\n"
      "    --write_lig_poses  <arg>\n"
      "                            <arg> is a PDB file into which ligand rotameric poses will be written\n"
      "    --by_orientation  <arg>\n"
      "                            This option is used with commend 'ScreenLigPoses'. <arg> is a file that defines the rules for screening ligand poses based on the ligand's orientations \n"
      "    --by_bbvdw_perc  <arg>\n"
      "                            <arg> is a value between 0 and 1 that specifies the top percentile threshold (e.g., 0.25) \n"
      "                            for selecting ligand poses whose internalVDW and backboneVDW scores both fall within this range\n"
      "    --by_invdw_perc  <arg>\n"
      "                            <arg> is a value between 0 and 1 that specifies the top percentile threshold (e.g., 0.25) \n"
      "                            for selecting ligand poses whose internalVDW and backboneVDW scores both fall within this range\n"
      "    --by_rmsd  <arg>\n"
      "                            <arg> is a value between 0 and 1 that specifies the top percentile threshold (e.g., 0.25) \n"
      "                            for selecting ligand poses whose internalVDW and backboneVDW scores both fall within this range\n"
      "    --lig_param  <arg> \n"
      "    --lig_topo <arg> \n"
      "    --lig_catacons <arg>  \n"
      "    --lig_placing <arg> \n"
      "    --num_of_runs <arg>  \n"
      "    --within_residues <arg>  \n"
      "    --within_range <arg> \n"
      "    --ntraj_start_ndx <arg> \n"
      "    --wbind <arg> \n"
      "    --residue <arg> \n"
      "    --residue_pair <arg> \n"
      "    --exclude_residues <arg> \n"
      "    --min_muts <arg> \n"
      "    --max_muts <arg>  \n"
      "    --penalize_identical_seqs \n"
      "    --min_temperature  \n"
      "    --max_temperature  \n"
      "    --num_of_samc_cycles \n"
      "    --num_of_metropolis_steps \n"
      "    --samc_temperature_decrease_factor \n"
      "    --write_mutations \n"
      "\n"
      ,PROGRAM_NAME, PROGRAM_NAME);
  return Success;
}

int PrintProgramVersion()
{
  printf("%s %s\n", PROGRAM_NAME, PROGRAM_VERSION);
  return Success;
}

int PrintProgramInfo()
{
  printf(
      "####################################################################################################################\n"
      "#                                       %s %s\n"
      "#  A command-line-based tool featuring protein structure modeling and protein design, including:\n"
      "#    *) de novo protein sequence design and protein sequence redesign (ProteinDesign)\n"
      "#       can be used to design monomer proteins, protein-protein interactions, protein-ligand interactions, \n"
      "#       protein-nucleic acid interactions, and enzymes\n"
      "#    *) protein side-chain packing (ProteinDesign in conjuncion with the option --wildtype_only)\n"
      "#    *) calculate protein fold stability score (ComputeStability)\n"
      "#    *) calculate protein-protein binding energy score (ComputeBinding)\n"
      "#    *) repair missing protein side chains (RepairStructure)\n"
      "#    *) protein structure energy minimization (Minimize)\n"
      "#    *) build mutant strucure models (BuildMutant)\n"
      "#    *) add polar hydrogen atoms (AddPolarHydrogen)\n"
      "#    *) optimize hydrogen's position to maximize hydrogen bonding networks (OptimizeHydrogen)\n"
      "#  Copyright (c) Xiaoqiang Huang\n"
      "#  Contact: xiaoqiah@outlook.com\n"
      "####################################################################################################################\n",
      PROGRAM_NAME, PROGRAM_VERSION);
  return Success;
}

int ComputeStability(Structure *pStruct, AAppTable *pAAPP, RamaTable *pRama, BBdepRotLib* pLib, double energyTerms[MAX_ENERGY_TERM])
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChainI = StructureGetChain(pStruct, i);
    for (int ir = 0; ir < ChainGetResidueNum(pChainI); ir++)
    {
      Residue *pResIR = ChainGetResidue(pChainI, ir);
      if (ChainGetType(pChainI) == Type_Chain_Protein)
      {
        AminoAcidReferenceEnergy(pResIR->name, energyTerms);
        EnergyIntraResidue(pResIR, energyTerms);
        AminoAcidPropensityAndRamachandranEnergy(pResIR, pAAPP, pRama);
        AminoAcidDunbrackEnergy(pResIR, pLib);
        energyTerms[91] += pResIR->aapp;
        energyTerms[92] += pResIR->rama;
        energyTerms[93] += pResIR->dunbrack;
      }
      for (int is = ir + 1; is < ChainGetResidueNum(pChainI); is++)
      {
        Residue *pResIS = ChainGetResidue(pChainI, is);
        if (ResidueGetPosInChain(pResIR) + 1 == ResidueGetPosInChain(pResIS))
        {
          EnergyResidueAndNextResidue(pResIR, pResIS, energyTerms);
        }
        else
        {
          EnergyResidueAndOtherResidueSameChain(pResIR, pResIS, energyTerms);
        }
      }
      for (int k = i + 1; k < StructureGetChainNum(pStruct); k++)
      {
        Chain *pChainK = StructureGetChain(pStruct, k);
        for (int ks = 0; ks < ChainGetResidueNum(pChainK); ks++)
        {
          Residue *pResKS = ChainGetResidue(pChainK, ks);
          if (ChainGetType(pChainI) == Type_Chain_SmallMol)
          {
            EnergyResidueAndLigandResidue(pResKS, pResIR, energyTerms);
          }
          else if (ChainGetType(pChainK) == Type_Chain_SmallMol)
          {
            EnergyResidueAndLigandResidue(pResIR, pResKS, energyTerms);
          }
          else
          {
            EnergyResidueAndOtherResidueDiffChain(pResIR, pResKS, energyTerms);
          }
        }
      }
    }
  }

  EnergyTermWeighting(energyTerms);
  printf("\nStability energy details:\n");
  EnergyTermShowComplex(energyTerms);
  return Success;
}

int ComputeBinding(Structure *pStructure)
{
  if (StructureGetChainNum(pStructure) > 2)
  {
    printf("Your structure has more than two protein chains, and you should specify how to split chains "
           "before computing the binding energy. Otherwise, the interactions between any chain pair is calculated\n");
  }
  else if (StructureGetChainNum(pStructure) <= 1)
  {
    printf("Your structure has less than or equal to one chain, binding energy cannot be calculated\n");
    return Warning;
  }

  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChainI = StructureGetChain(pStructure, i);
    for (int k = i + 1; k < StructureGetChainNum(pStructure); k++)
    {
      Chain *pChainK = StructureGetChain(pStructure, k);
      double energyTerms[MAX_ENERGY_TERM] = {0};
      for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
      {
        Residue *pResIJ = ChainGetResidue(pChainI, j);
        for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
        {
          Residue *pResKS = ChainGetResidue(pChainK, s);
          if (ChainGetType(pChainI) == Type_Chain_SmallMol)
          {
            EnergyResidueAndLigandResidue(pResKS, pResIJ, energyTerms);
          }
          else if (ChainGetType(pChainK) == Type_Chain_SmallMol)
          {
            EnergyResidueAndLigandResidue(pResIJ, pResKS, energyTerms);
          }
          else
          {
            EnergyResidueAndOtherResidueDiffChain(pResIJ, pResKS, energyTerms);
          }
        }
      }
      EnergyTermWeighting(energyTerms);
      printf("Structure binding energy details between chain %s and chain %s:\n", ChainGetName(pChainI), ChainGetName(pChainK));
      EnergyTermShowComplex(energyTerms);
    }
  }

  return Success;
}

int ComputeBinding2(Structure *pStructure, char split1[], char split2[])
{
  double energyTerms[MAX_ENERGY_TERM] = {0};
  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChainI = StructureGetChain(pStructure, i);
    for (int k = i + 1; k < StructureGetChainNum(pStructure); k++)
    {
      Chain *pChainK = StructureGetChain(pStructure, k);
      if ( (strstr(split1, ChainGetName(pChainI)) != NULL && strstr(split1, ChainGetName(pChainK)) == NULL) || 
      ((strstr(split1, ChainGetName(pChainI)) == NULL && strstr(split1, ChainGetName(pChainK)) != NULL)) )
      {
        for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
        {
          Residue *pResIJ = ChainGetResidue(pChainI, j);
          for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
          {
            Residue *pResKS = ChainGetResidue(pChainK, s);
            if (ChainGetType(pChainI) == Type_Chain_SmallMol)
            {
              EnergyResidueAndLigandResidue(pResKS, pResIJ, energyTerms);
            }
            else if (ChainGetType(pChainK) == Type_Chain_SmallMol)
            {
              EnergyResidueAndLigandResidue(pResIJ, pResKS, energyTerms);
            }
            else
            {
              EnergyResidueAndOtherResidueDiffChain(pResIJ, pResKS, energyTerms);
            }
          }
        }
      }
    }
  }
  EnergyTermWeighting(energyTerms);
  printf("Binding energy details between chain(s) %s and chain(s) %s "
  "(DG_bind = DG(stability,complex) - DG(stability,%s) - DG(stability,%s):\n", split1, split2, split1, split2);
  EnergyTermShowComplex(energyTerms);
  return Success;
}

int BuildMutant(Structure *pStruct, char *mutfile, BBdepRotLib *pLib, AtomParams *pParams, ResTopos *pTops)
{
  FileReader fr;
  FileReaderCreate(&fr, mutfile);
  int mutNum = FileReaderGetLineNum(&fr);
  if (mutNum <= 0)
  {
    printf("in file %s line %d, no mutation was found in file %s\n", __FILE__, __LINE__, mutfile);
    exit(DataNotExistError);
  }
  StringArray *mutants = (StringArray *)malloc(sizeof(StringArray) * mutNum);
  char line[MAX_LEN_ONE_LINE_CONTENT + 1];
  int mutNdx = 0;
  while (!FAILED(FileReaderGetNextLine(&fr, line)))
  {
    StringArrayCreate(&mutants[mutNdx]);
    StringArraySplitString(&mutants[mutNdx], line, ',');
    char lastMutant[MAX_LEN_ONE_LINE_CONTENT + 1];
    int lastMutNdx = StringArrayGetNum(&mutants[mutNdx]) - 1;
    strcpy(lastMutant, StringArrayGet(&mutants[mutNdx], lastMutNdx));
    // deal with the last char of the last single mutant
    if ((!isdigit(lastMutant[strlen(lastMutant) - 1])) && !isalpha(lastMutant[strlen(lastMutant) - 1]))
    {
      lastMutant[strlen(lastMutant) - 1] = '\0';
    }
    StringArraySet(&mutants[mutNdx], lastMutNdx, lastMutant);
    mutNdx++;
  }
  FileReaderDestroy(&fr);

  for (int mutNdx = 0; mutNdx < mutNum; mutNdx++)
  {
    Structure tempStruct;
    StructureCreate(&tempStruct);
    StructureCopy(&tempStruct, pStruct);
    Structure tempStruct2;
    StructureCreate(&tempStruct2);
    StructureCopy(&tempStruct2, pStruct);

    // 1.deal with mutations on the nucleotide
    IntArray flagNucleotides;
    IntArrayCreate(&flagNucleotides, StringArrayGetNum(&mutants[mutNdx]));
    for (int mutation = 0; mutation < StringArrayGetNum(&mutants[mutNdx]); mutation++)
    {
      char mutstr[10];
      char aa1, chain, aa2;
      int posInChain;
      strcpy(mutstr, StringArrayGet(&mutants[mutNdx], mutation));
      sscanf(mutstr, "%c%c%d%c", &aa1, &chain, &posInChain, &aa2);
      int chainIndex = -1, resIndex = -1;
      char chainName[MAX_LEN_CHAIN_NAME];
      chainName[0] = chain;
      chainName[1] = '\0';
      StructureFindChainIndex(&tempStruct, chainName, &chainIndex);
      if (chainIndex == -1)
      {
        printf("in file %s line %d, cannot find mutation %s due to an incorrect chain id\n", __FILE__, __LINE__, mutstr);
        exit(ValueError);
      }
      ChainFindResidueIndexByPosInChain(StructureGetChain(&tempStruct, chainIndex), posInChain, &resIndex);
      if (resIndex == -1)
      {
        printf("in file %s line %d, cannot find mutation %s due to an incorrect position\n", __FILE__, __LINE__, mutstr);
        exit(ValueError);
      }
      Chain *pChain = StructureFindChainByName(&tempStruct, chainName);
      if (ChainGetType(pChain) == Type_Chain_DNA || ChainGetType(pChain) == Type_Chain_RNA)
      {
        Residue *pTargetNucleotide = ChainGetResidue(StructureGetChain(&tempStruct, chainIndex), resIndex);
        StructureMutateNucleotide(pTargetNucleotide, aa2, &tempStruct, pParams, pTops);
        IntArraySet(&flagNucleotides, mutation, 1);
      }
    }

    // 2. for each mutant, build the rotamer-tree
    IntArray mutatedArray, rotamericArray;
    IntArrayCreate(&mutatedArray, 0);
    IntArrayCreate(&rotamericArray, 0);
    for (int mutation = 0; mutation < StringArrayGetNum(&mutants[mutNdx]); mutation++)
    {
      // skip the nucleotide mutation
      if (IntArrayGet(&flagNucleotides, mutation) == 1)
        continue;
      char mutstr[10];
      char aa1, chn, aa2;
      int posInChain;
      strcpy(mutstr, StringArrayGet(&mutants[mutNdx], mutation));
      sscanf(mutstr, "%c%c%d%c", &aa1, &chn, &posInChain, &aa2);
      int chnNdx = -1, resNdx = -1;
      char chnName[MAX_LEN_CHAIN_NAME];
      chnName[0] = chn;
      chnName[1] = '\0';
      StructureFindChainIndex(&tempStruct, chnName, &chnNdx);
      if (chnNdx == -1)
      {
        printf("in file %s line %d, cannot find mutation %s\n", __FILE__, __LINE__, mutstr);
        exit(ValueError);
      }
      ChainFindResidueIndexByPosInChain(StructureGetChain(&tempStruct, chnNdx), posInChain, &resNdx);
      if (resNdx == -1)
      {
        printf("in file %s line %d, cannot find mutation %s\n", __FILE__, __LINE__, mutstr);
        exit(ValueError);
      }
      char mutAAType[MAX_LEN_RES_NAME];
      AA1ToAA3(aa2, mutAAType);
      StringArray mutRotType, mutPatch;
      StringArrayCreate(&mutRotType);
      StringArrayCreate(&mutPatch);
      StringArrayAppend(&mutRotType, mutAAType);
      StringArrayAppend(&mutPatch, "");
      // for histidine, the default mutAAType is HSD, we need to add HSE
      if (aa2 == 'H')
      {
        StringArrayAppend(&mutRotType, "HSE");
        StringArrayAppend(&mutPatch, "");
      }
      ProteinSiteBuildMutatedRotamersByBBdepRotLib(&tempStruct, chnNdx, resNdx, pLib, pParams, pTops, &mutRotType, &mutPatch);
      StringArrayDestroy(&mutRotType);
      StringArrayDestroy(&mutPatch);

      IntArrayAppend(&mutatedArray, chnNdx);
      IntArrayAppend(&mutatedArray, resNdx);
      IntArrayAppend(&rotamericArray, chnNdx);
      IntArrayAppend(&rotamericArray, resNdx);

      // create a reference model for the wild-type
      StringArray wtRotType, wtPatch;
      StringArrayCreate(&wtRotType);
      StringArrayCreate(&wtPatch);
      char wtAAType[MAX_LEN_RES_NAME];
      AA1ToAA3(aa1, wtAAType);
      StringArrayAppend(&wtRotType, wtAAType);
      StringArrayAppend(&wtPatch, "");
      ProteinSiteBuildMutatedRotamersByBBdepRotLib(&tempStruct2, chnNdx, resNdx, pLib, pParams, pTops, &wtRotType, &wtPatch);
      StringArrayDestroy(&mutRotType);
      StringArrayDestroy(&mutPatch);
    }
    IntArrayDestroy(&flagNucleotides);

    // 3. build rots for surrounding residues
    for (int i = 0; i < IntArrayGetLength(&mutatedArray); i += 2)
    {
      int chnNdx = IntArrayGet(&mutatedArray, i);
      int resNdx = IntArrayGet(&mutatedArray, i + 1);
      Residue *pResi1 = ChainGetResidue(StructureGetChain(&tempStruct, chnNdx), resNdx);
      for (int j = 0; j < StructureGetChainNum(&tempStruct); j++)
      {
        Chain *pChain = StructureGetChain(&tempStruct, j);
        if (ChainGetType(pChain) != Type_Chain_Protein)
          continue;
        for (int k = 0; k < ChainGetResidueNum(pChain); k++)
        {
          Residue *pResi2 = ChainGetResidue(pChain, k);
          if (AtomArrayCalcMinDistance(&pResi1->atoms, &pResi2->atoms) < ENERGY_DISTANCE_CUTOFF)
          {
            if (pResi2->desType == Type_DesType_Fixed)
            {
              ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&tempStruct, j, k, pLib, pParams, pTops);
              ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&tempStruct2, j, k, pLib, pParams, pTops);
              if (pResi2->isSCIntact)
              {
                ProteinSiteBuildNativeRotamer(&tempStruct, j, k, pTops);
                ProteinSiteBuildNativeRotamer(&tempStruct2, j, k, pTops);
              }
              IntArrayAppend(&rotamericArray, j);
              IntArrayAppend(&rotamericArray, k);
            }
          }
        }
      }
    }

    // 4.optimization rots sequentially
    printf("Building mutation model %d, the following sites will be optimized:\n", mutNdx + 1);
    printf("chnIndex resIndex (both of them starts from zero on the chain)\n");
    for (int i = 0; i < IntArrayGetLength(&rotamericArray); i += 2)
    {
      printf("%8d %8d\n", IntArrayGet(&rotamericArray, i), IntArrayGet(&rotamericArray, i + 1));
    }
    for (int run = 0; run < MAX_NUM_OF_RUNS; run++)
    {
      printf("optimization run %d ... \n", run + 1);
      for (int i = 0; i < IntArrayGetLength(&rotamericArray); i += 2)
      {
        int chnNdx = IntArrayGet(&rotamericArray, i);
        int resNdx = IntArrayGet(&rotamericArray, i + 1);
        ProteinSiteOptimizeRotamer(&tempStruct, chnNdx, resNdx, pLib);
        ProteinSiteOptimizeRotamer(&tempStruct2, chnNdx, resNdx, pLib);
      }
    }
    IntArrayDestroy(&mutatedArray);
    IntArrayDestroy(&rotamericArray);

    // delete rots for previous mutant
    StructureRemoveAllDesignSites(&tempStruct);
    StructureRemoveAllDesignSites(&tempStruct2);

    // finally, write the mutant and reference WT model into files
    char mut_model_file[MAX_LEN_ONE_LINE_CONTENT + 1];
    char wt_model_file[MAX_LEN_ONE_LINE_CONTENT + 1];
    sprintf(mut_model_file, "UD_model_%04d.pdb", mutNdx + 1);
    sprintf(wt_model_file, "UD_model_%04d_WT.pdb", mutNdx + 1);
    FILE *pFile = fopen(mut_model_file, "w");
    fprintf(pFile, "REMARK file generated by module <BuildMutant> with backbone-dependent rotamer library\n");
    StructureShowInPDBFormat(&tempStruct, pFile);
    fclose(pFile);

    pFile = fopen(wt_model_file, "w");
    fprintf(pFile, "REMARK file generated by module <BuildMutant> with backbone-dependent rotamer library\n");
    StructureShowInPDBFormat(&tempStruct2, pFile);
    fclose(pFile);

    StructureDestroy(&tempStruct);
    StructureDestroy(&tempStruct2);
  }

  return Success;
}

int RepairStructure(Structure *pStruct, BBdepRotLib *pLib, AtomParams *pParams, ResTopos *pTopos)
{
  for (int run = 0; run < MAX_NUM_OF_RUNS; run++)
  {
    printf("Repairing structure run %d ... \n", run + 1);
    for (int i = 0; i < StructureGetChainNum(pStruct); ++i)
    {
      Chain *pChain = StructureGetChain(pStruct, i);
      if (ChainGetType(pChain) != Type_Chain_Protein)
      {
        continue;
      }
      for (int j = 0; j < ChainGetResidueNum(pChain); j++)
      {
        Residue *pResi = ChainGetResidue(pChain, j);
        if (pResi->isSCIntact == TRUE)
        {
          if (strcmp(ResidueGetName(pResi), "ALA") == 0 || strcmp(ResidueGetName(pResi), "GLY") == 0)
            continue;
          // skip CYS which may form disulfide bonds
          if (strcmp(ResidueGetName(pResi), "CYS") == 0)
            continue;
          if (strcmp(ResidueGetName(pResi), "ASN") == 0 || strcmp(ResidueGetName(pResi), "GLN") == 0 || strcmp(ResidueGetName(pResi), "HSD") == 0 || strcmp(ResidueGetName(pResi), "HSE") == 0)
          {
            printf("flip %s%d%s sidechain to optimize H-bonds\n", ResidueGetChainName(pResi), ResidueGetPosInChain(pResi), ResidueGetName(pResi));
            ProteinSiteBuildNativeRotamer(pStruct, i, j, pTopos);
            ProteinSiteBuildFlippedNativeRotamer(pStruct, i, j, pTopos);
            ProteinSiteOptimizeRotamerHBondEnergy(pStruct, i, j);
            ProteinSiteRemoveDesignSite(pStruct, i, j);
          }
          else if (strcmp(ResidueGetName(pResi), "SER") == 0 || strcmp(ResidueGetName(pResi), "THR") == 0 || strcmp(ResidueGetName(pResi), "TYR") == 0)
          {
            printf("rotate %s%d%s hydroxyl group to optimize H-bonds\n", ResidueGetChainName(pResi), ResidueGetPosInChain(pResi), ResidueGetName(pResi));
            ProteinSiteBuildNativeRotamer(pStruct, i, j, pTopos);
            ProteinSiteExpandHydroxylRotamers(pStruct, i, j, pTopos);
            ProteinSiteOptimizeRotamerHBondEnergy(pStruct, i, j);
            ProteinSiteRemoveDesignSite(pStruct, i, j);
          }
        }
        else
        {
          printf("repair %s%d%s sidechain to complete atoms\n", ResidueGetChainName(pResi), ResidueGetPosInChain(pResi), ResidueGetName(pResi));
          ResidueCalcAllAtomXYZ(pResi, pTopos, ChainGetResidue(pChain, j - 1), ChainGetResidue(pChain, j + 1));
          //ResidueCalcSidechainTorsion(pResi,pTopos);
          ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStruct, i, j, pLib, pParams, pTopos);
          ProteinSiteExpandHydroxylRotamers(pStruct, i, j, pTopos);
          ProteinSiteOptimizeRotamer(pStruct, i, j, pLib);
          ProteinSiteRemoveDesignSite(pStruct, i, j);
        }
      }
    }
  }

  FILE *pFile = fopen("UD_RepairedStructure.pdb", "w");
  fprintf(pFile, "REMARK file generated by module <RepairStructure>\n");
  StructureShowInPDBFormat(pStruct, pFile);
  fclose(pFile);

  return Success;
}

int EnergyMinimization(Structure *pStructure, BBdepRotLib *pBBdepRotLib, AtomParams *atomParams, ResTopos *resiTopos)
{
  for (int run = 0; run < MAX_NUM_OF_RUNS; run++)
  {
    printf("Minimzation run %d ... \n", run + 1);
    for (int i = 0; i < StructureGetChainNum(pStructure); ++i)
    {
      Chain *pChain = StructureGetChain(pStructure, i);
      if (ChainGetType(pChain) != Type_Chain_Protein)
        continue;
      for (int j = 0; j < ChainGetResidueNum(pChain); j++)
      {
        Residue *pResi = ChainGetResidue(pChain, j);
        if (pResi->isSCIntact == TRUE)
        {
          if (strcmp(ResidueGetName(pResi), "ALA") == 0 || strcmp(ResidueGetName(pResi), "GLY") == 0)
            continue;
          // skip CYS which may form disulfide bonds
          if (strcmp(ResidueGetName(pResi), "CYS") == 0)
            continue;
          if (strcmp(ResidueGetName(pResi), "ASN") == 0 || strcmp(ResidueGetName(pResi), "GLN") == 0 || strcmp(ResidueGetName(pResi), "HSD") == 0 || strcmp(ResidueGetName(pResi), "HSE") == 0)
          {
            printf("flip %s%d%s to optimize H-bonds\n", ResidueGetChainName(pResi), ResidueGetPosInChain(pResi), ResidueGetName(pResi));
            ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
            ProteinSiteBuildFlippedNativeRotamer(pStructure, i, j, resiTopos);
            ProteinSiteOptimizeRotamerHBondEnergy(pStructure, i, j);
            ProteinSiteRemoveDesignSite(pStructure, i, j);
          }
          else if (strcmp(ResidueGetName(pResi), "SER") == 0 || strcmp(ResidueGetName(pResi), "THR") == 0 || strcmp(ResidueGetName(pResi), "TYR") == 0)
          {
            printf("rotate %s%d%s hydroxyl group to optimize H-bonds\n", ResidueGetChainName(pResi), ResidueGetPosInChain(pResi), ResidueGetName(pResi));
            ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
            ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
            ProteinSiteOptimizeRotamerHBondEnergy(pStructure, i, j);
            ProteinSiteRemoveDesignSite(pStructure, i, j);
          }
          printf("optimize %s%d%s sidechain conformation\n", ResidueGetChainName(pResi), ResidueGetPosInChain(pResi), ResidueGetName(pResi));
          ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure, i, j, pBBdepRotLib, atomParams, resiTopos);
          ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
          ProteinSiteOptimizeRotamer(pStructure, i, j, pBBdepRotLib);
          ProteinSiteRemoveDesignSite(pStructure, i, j);
        }
        else
        {
          printf("repair %s%d%s sidechain and optimize conformation\n", ResidueGetChainName(pResi), ResidueGetPosInChain(pResi), ResidueGetName(pResi));
          ResidueCalcAllAtomXYZ(pResi, resiTopos, ChainGetResidue(pChain, j - 1), ChainGetResidue(pChain, j + 1));
          // ResidueCalcSidechainTorsion(pResi,resiTopos);
          ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure, i, j, pBBdepRotLib, atomParams, resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
          ProteinSiteOptimizeRotamer(pStructure, i, j, pBBdepRotLib);
          ProteinSiteRemoveDesignSite(pStructure, i, j);
        }
      }
    }
  }

  FILE *pOut = fopen("UD_sidechainEnergyMinimized.pdb", "w");
  fprintf(pOut, "REMARK file generated by module <Minimize>\n");
  StructureShowInPDBFormat(pStructure, pOut);
  fclose(pOut);

  return Success;
}

int AddPolarHydrogen(Structure *pStructure)
{
  FILE *pOut = fopen("UD_polarHydrogensAdded.pdb", "w");
  fprintf(pOut, "REMARK file generated by module <AddPolarHydrogen>\n");
  StructureShowInPDBFormat(pStructure, pOut);
  fclose(pOut);
  return Success;
}

int OptimizeHydrogen(Structure *pStructure, AtomParams *atomParams, ResTopos *resiTopos)
{
  for (int run = 0; run < MAX_NUM_OF_RUNS; run++)
  {
    printf("Optimizing hydrogen run %d \n", run + 1);
    for (int i = 0; i < StructureGetChainNum(pStructure); ++i)
    {
      Chain *pChain = StructureGetChain(pStructure, i);
      for (int j = 0; j < ChainGetResidueNum(pChain); j++)
      {
        Residue *pResi = ChainGetResidue(pChain, j);
        if (strcmp(ResidueGetName(pResi), "SER") == 0 ||
            strcmp(ResidueGetName(pResi), "THR") == 0 ||
            strcmp(ResidueGetName(pResi), "TYR") == 0)
        {
          printf("rotating hydroxyl group of residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi), ResidueGetPosInChain(pResi), AA3ToAA1(ResidueGetName(pResi)));
          if (pResi->isSCIntact)
            ProteinSiteBuildNativeRotamer(pStructure, i, j, resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure, i, j, resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure, i, j);
          ProteinSiteRemoveDesignSite(pStructure, i, j);
        }
      }
    }
  }

  // output the struture with optimized polar hydrogens
  FILE *pFile = fopen("UD_polarHydrogensOptimized.pdb", "w");
  fprintf(pFile, "REMARK file generated by module <OptimizeHydrogen>\n");
  StructureShowInPDBFormat(pStructure, pFile);
  fclose(pFile);

  return Success;
}

int ComputeResidueInteractionWithFixedEnvironment(Structure *pStruct, AAppTable *pAAppTable, RamaTable *pRama, BBdepRotLib* pBBdepRotLib, int chnNdx, int resNdx)
{
  Chain *pChainI = StructureGetChain(pStruct, chnNdx);
  Residue *pResIR = ChainGetResidue(pChainI, resNdx);
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;

  StringArray resStrs;
  StringArrayCreate(&resStrs);
  if (strcmp(EXCLUDE_RESIDUES, "") != 0)
  {
    StringArraySplitString(&resStrs, EXCLUDE_RESIDUES, ',');
  }
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChainI = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
    {
      Residue *pResi2 = ChainGetResidue(pChainI, j);
      if (strcmp(ResidueGetName(pResIR), ResidueGetName(pResi2)) == 0 && ResidueGetPosInChain(pResIR) == ResidueGetPosInChain(pResi2))
      {
        continue;
      }
      else if (AtomArrayCalcMinDistance(&pResIR->atoms, &pResi2->atoms) < ENERGY_DISTANCE_CUTOFF)
      {
        char str[10] = "";
        sprintf(str, "%c%d", ResidueGetChainName(pResi2)[0], ResidueGetPosInChain(pResi2));
        int exist = -1;
        StringArrayFind(&resStrs, str, &exist);
        if (exist >= 0)
        {
          continue;
        }
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue *) * surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum - 1] = pResi2;
      }
    }
  }
  StringArrayDestroy(&resStrs);

  double energyTerms[MAX_ENERGY_TERM] = {0};
  if (ChainGetType(pChainI) == Type_Chain_Protein)
  {
    AminoAcidReferenceEnergy(ResidueGetName(pResIR), energyTerms);
    EnergyIntraResidue(pResIR, energyTerms);
    AminoAcidPropensityAndRamachandranEnergy(pResIR, pAAppTable, pRama);
    AminoAcidDunbrackEnergy(pResIR, pBBdepRotLib);
    energyTerms[91] += pResIR->aapp;
    energyTerms[92] += pResIR->rama;
    energyTerms[93] += pResIR->dunbrack;
  }
  for (int is = 0; is < surroundingResiNum; is++)
  {
    Residue *pResIS = ppSurroundingResidues[is];
    if (!strcmp(ResidueGetChainName(pResIR), ResidueGetChainName(pResIS)))
    {
      if (ResidueGetPosInChain(pResIR) == ResidueGetPosInChain(pResIS) - 1)
      {
        EnergyResidueAndNextResidue(pResIR, pResIS, energyTerms);
      }
      else if (ResidueGetPosInChain(pResIR) == ResidueGetPosInChain(pResIS) + 1)
      {
        EnergyResidueAndNextResidue(pResIS, pResIR, energyTerms);
      }
      else
      {
        EnergyResidueAndOtherResidueSameChain(pResIR, pResIS, energyTerms);
      }
    }
    else
    {
      if (ChainGetType(StructureFindChainByName(pStruct, ResidueGetChainName(pResIS))) == Type_Chain_SmallMol)
      {
        EnergyResidueAndLigandResidue(pResIR, pResIS, energyTerms);
      }
      else
      {
        EnergyResidueAndOtherResidueDiffChain(pResIR, pResIS, energyTerms);
      }
    }
  }

  EnergyTermWeighting(energyTerms);
  EnergyTermShowComplex(energyTerms);

  return Success;
}

int FindInterfaceResidues(Structure *pStructure)
{
  IntArray *flagArrays = (IntArray *)malloc(sizeof(IntArray) * StructureGetChainNum(pStructure));
  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChainI = StructureGetChain(pStructure, i);
    IntArrayCreate(&flagArrays[i], ChainGetResidueNum(pChainI));
    for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
    {
      IntArraySet(&flagArrays[i], j, 0);
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChainI = StructureGetChain(pStructure, i);
    for (int k = i + 1; k < StructureGetChainNum(pStructure); k++)
    {
      Chain *pChainK = StructureGetChain(pStructure, k);
      for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
      {
        Residue *pResiIJ = ChainGetResidue(pChainI, j);
        for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
        {
          Residue *pResiKS = ChainGetResidue(pChainK, s);
          if (AtomArrayCalcMinDistance(&pResiIJ->atoms, &pResiKS->atoms) < CUT_PPI_DIST_SHELL1)
          {
            IntArraySet(&flagArrays[i], j, 1);
            IntArraySet(&flagArrays[k], s, 1);
          }
        }
      }
    }
  }

  printf("Interface residues: \n");
  for (int j = 0; j < StructureGetChainNum(pStructure); j++)
  {
    Chain *pChain = StructureGetChain(pStructure, j);
    if (ChainGetType(pChain) == Type_Chain_Protein || ChainGetType(pChain) == Type_Chain_DNA || ChainGetType(pChain) == Type_Chain_RNA)
    {
      for (int i = 0; i < ChainGetResidueNum(pChain); i++)
      {
        Residue *pResidue = ChainGetResidue(pChain, i);
        if (IntArrayGet(&flagArrays[j], i) == 1)
        {
          printf("%s %4d %c\n", ResidueGetChainName(pResidue), ResidueGetPosInChain(pResidue), AA3ToAA1(ResidueGetName(pResidue)));
        }
      }
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    IntArrayDestroy(&flagArrays[i]);
  }
  free(flagArrays);
  flagArrays = NULL;

  return Success;
}

int FindInterfaceResidues2(Structure *pStructure, char split1[], char split2[])
{
  IntArray *flagArrays = (IntArray *)malloc(sizeof(IntArray) * StructureGetChainNum(pStructure));
  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChainI = StructureGetChain(pStructure, i);
    IntArrayCreate(&flagArrays[i], ChainGetResidueNum(pChainI));
    for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
    {
      IntArraySet(&flagArrays[i], j, 0);
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChainI = StructureGetChain(pStructure, i);
    for (int k = i + 1; k < StructureGetChainNum(pStructure); k++)
    {
      Chain *pChainK = StructureGetChain(pStructure, k);
      if ((strstr(split1, ChainGetName(pChainI)) != NULL && strstr(split1, ChainGetName(pChainK)) == NULL) || ((strstr(split1, ChainGetName(pChainI)) == NULL && strstr(split1, ChainGetName(pChainK)) != NULL)))
      {
        for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
        {
          Residue *pResiIJ = ChainGetResidue(pChainI, j);
          for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
          {
            Residue *pResiKS = ChainGetResidue(pChainK, s);
            if (AtomArrayCalcMinDistance(&pResiIJ->atoms, &pResiKS->atoms) < CUT_PPI_DIST_SHELL1)
            {
              IntArraySet(&flagArrays[i], j, 1);
              IntArraySet(&flagArrays[k], s, 1);
            }
          }
        }
      }
    }
  }

  printf("Interface residues: \n");
  for (int j = 0; j < StructureGetChainNum(pStructure); j++)
  {
    Chain *pChain = StructureGetChain(pStructure, j);
    if (ChainGetType(pChain) == Type_Chain_Protein || ChainGetType(pChain) == Type_Chain_DNA || ChainGetType(pChain) == Type_Chain_RNA)
    {
      for (int i = 0; i < ChainGetResidueNum(pChain); i++)
      {
        Residue *pResidue = ChainGetResidue(pChain, i);
        if (IntArrayGet(&flagArrays[j], i) == 1)
        {
          printf("%s %4d %c\n", ResidueGetChainName(pResidue), ResidueGetPosInChain(pResidue), AA3ToAA1(ResidueGetName(pResidue)));
        }
      }
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    IntArrayDestroy(&flagArrays[i]);
  }
  free(flagArrays);
  flagArrays = NULL;

  return Success;
}

int FindCoreResidues(Structure *pStructure)
{
  StructureComputeResiduePosition(pStructure);
  printf("Core residues:\n");
  for (int j = 0; j < StructureGetChainNum(pStructure); j++)
  {
    Chain *pChain = StructureGetChain(pStructure, j);
    if (ChainGetType(pChain) == Type_Chain_Protein)
    {
      for (int i = 0; i < ChainGetResidueNum(pChain); i++)
      {
        Residue *pResidue = ChainGetResidue(pChain, i);
        if (pResidue->nCbIn10A > CUT_NUM_CB_CORE)
        {
          printf("%s %4d %c\n", ResidueGetChainName(pResidue), ResidueGetPosInChain(pResidue), AA3ToAA1(ResidueGetName(pResidue)));
        }
      }
    }
  }
  return Success;
}

int FindSurfaceResidues(Structure *pStructure)
{
  StructureComputeResiduePosition(pStructure);
  printf("Surface residues:\n");
  for (int j = 0; j < StructureGetChainNum(pStructure); j++)
  {
    Chain *pChain = StructureGetChain(pStructure, j);
    if (ChainGetType(pChain) == Type_Chain_Protein)
    {
      for (int i = 0; i < ChainGetResidueNum(pChain); i++)
      {
        Residue *pResidue = ChainGetResidue(pChain, i);
        if (pResidue->nCbIn10A < CUT_NUM_CB_SURF)
        {
          printf("%s %4d %c\n", ResidueGetChainName(pResidue), ResidueGetPosInChain(pResidue), AA3ToAA1(ResidueGetName(pResidue)));
        }
      }
    }
  }
  return Success;
}

int FindIntermediateResidues(Structure *pStructure)
{
  StructureComputeResiduePosition(pStructure);
  printf("Intermediate residues:\n");
  for (int j = 0; j < StructureGetChainNum(pStructure); j++)
  {
    Chain *pChain = StructureGetChain(pStructure, j);
    if (ChainGetType(pChain) == Type_Chain_Protein)
    {
      for (int i = 0; i < ChainGetResidueNum(pChain); i++)
      {
        Residue *pResidue = ChainGetResidue(pChain, i);
        if (pResidue->nCbIn10A <= CUT_NUM_CB_CORE && pResidue->nCbIn10A >= CUT_NUM_CB_SURF)
        {
          printf("%s %4d %c\n", ResidueGetChainName(pResidue), ResidueGetPosInChain(pResidue), AA3ToAA1(ResidueGetName(pResidue)));
        }
      }
    }
  }
  return Success;
}

int WritePhiPsi(Structure *pStruct, char *filename)
{
  FILE *fout = NULL;
  if (filename == NULL)
    fout = stdout;
  else
    fout = fopen(filename, "w");
  fprintf(fout, "#phi-psi angles of protein residues: \n");
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChain = StructureGetChain(pStruct, i);
    if (ChainGetType(pChain) != Type_Chain_Protein)
      continue;
    for (int j = 0; j < ChainGetResidueNum(pChain); j++)
    {
      Residue *pResi1 = ChainGetResidue(pChain, j);
      fprintf(fout, "%s %s %4d %8.3f %8.3f\n", ResidueGetChainName(pResi1), ResidueGetName(pResi1), ResidueGetPosInChain(pResi1), pResi1->phipsi[0], pResi1->phipsi[1]);
    }
  }
  if (fout != stdout)
    fclose(fout);

  return Success;
}

int ComputeRotamersEnergy(Structure *pStruct, BBdepRotLib *pLib, AAppTable *pAAPP, RamaTable *pRama, AtomParams *pParams, ResTopos *pTops)
{
  FILE *pFile = fopen("UD_rotamers_energy.txt", "w");
  for (int i = 0; i < StructureGetChainNum(pStruct); ++i)
  {
    Chain *pChain = StructureGetChain(pStruct, i);
    if (ChainGetType(pChain) != Type_Chain_Protein)
      continue;
    for (int j = 0; j < ChainGetResidueNum(pChain); j++)
    {
      Residue *pRes = ChainGetResidue(pChain, j);
      if (FLAG_DESIGN_PPI)
      {
        for (int k = 0; k < StructureGetChainNum(pStruct); k++)
        {
          if (k == i)
            continue;
          Chain *pChainK = StructureGetChain(pStruct, k);
          for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
          {
            Residue *pResKS = ChainGetResidue(pChainK, s);
            if (AtomArrayCalcMinDistance(&pRes->atoms, &pResKS->atoms) < ENERGY_DISTANCE_CUTOFF)
            {
              ProteinSiteBuildAllRotamersByBBdepRotLib(pStruct, i, j, pLib, pParams, pTops);
              if (FLAG_USE_INPUT_SC && pRes->isSCIntact)
                ProteinSiteBuildNativeRotamer(pStruct, i, j, pTops);
              if (FLAG_ROTATE_HYDROXYL)
                ProteinSiteExpandHydroxylRotamers(pStruct, i, j, pTops);
              ProteinSiteCalcRotamersEnergy(pStruct, pAAPP, pRama, i, j, pFile);
              ProteinSiteRemoveDesignSite(pStruct, i, j);
            }
          }
        }
      }
      else if (FLAG_DESIGN_ENZ || FLAG_DESIGN_PLI)
      {
        Residue *pSmallMol = NULL;
        StructureFindSmallMol(pStruct, &pSmallMol);
        if (AtomArrayCalcMinDistance(&pRes->atoms, &pSmallMol->atoms) < ENERGY_DISTANCE_CUTOFF)
        {
          ProteinSiteBuildAllRotamersByBBdepRotLib(pStruct, i, j, pLib, pParams, pTops);
          if (FLAG_USE_INPUT_SC && pRes->isSCIntact)
            ProteinSiteBuildNativeRotamer(pStruct, i, j, pTops);
          if (FLAG_ROTATE_HYDROXYL)
            ProteinSiteExpandHydroxylRotamers(pStruct, i, j, pTops);
          ProteinSiteCalcRotamersEnergy(pStruct, pAAPP, pRama, i, j, pFile);
          ProteinSiteRemoveDesignSite(pStruct, i, j);
        }
      }
      else
      {
        ProteinSiteBuildAllRotamersByBBdepRotLib(pStruct, i, j, pLib, pParams, pTops);
        if (FLAG_USE_INPUT_SC && pRes->isSCIntact)
          ProteinSiteBuildNativeRotamer(pStruct, i, j, pTops);
        if (FLAG_ROTATE_HYDROXYL)
          ProteinSiteExpandHydroxylRotamers(pStruct, i, j, pTops);
        ProteinSiteCalcRotamersEnergy(pStruct, pAAPP, pRama, i, j, pFile);
        ProteinSiteRemoveDesignSite(pStruct, i, j);
      }
    }
  }
  fclose(pFile);

  return Success;
}

int ComputeWildtypeRotamersEnergy(Structure *pStruct, BBdepRotLib *pLib, AAppTable *pAAPP, RamaTable *pRama, AtomParams *pParams, ResTopos *pTops)
{
  FILE *pFile = fopen("UD_rotamers_energy.txt", "w");
  for (int i = 0; i < StructureGetChainNum(pStruct); ++i)
  {
    Chain *pChainI = StructureGetChain(pStruct, i);
    if (ChainGetType(pChainI) != Type_Chain_Protein)
      continue;
    for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
    {
      Residue *pResIJ = ChainGetResidue(pChainI, j);
      if (FLAG_DESIGN_PPI)
      {
        // check if residue is an interface residue
        for (int k = 0; k < StructureGetChainNum(pStruct); k++)
        {
          if (k == i)
            continue;
          Chain *pChainK = StructureGetChain(pStruct, k);
          for (int s = 0; s < ChainGetResidueNum(pChainK); s++)
          {
            Residue *pResKS = ChainGetResidue(pChainK, s);
            if (AtomArrayCalcMinDistance(&pResIJ->atoms, &pResKS->atoms) < ENERGY_DISTANCE_CUTOFF)
            {
              ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStruct, i, j, pLib, pParams, pTops);
              if (FLAG_USE_INPUT_SC && pResIJ->isSCIntact)
                ProteinSiteBuildNativeRotamer(pStruct, i, j, pTops);
              if (FLAG_ROTATE_HYDROXYL)
                ProteinSiteExpandHydroxylRotamers(pStruct, i, j, pTops);
              ProteinSiteCalcRotamersEnergy(pStruct, pAAPP, pRama, i, j, pFile);
              ProteinSiteRemoveDesignSite(pStruct, i, j);
            }
          }
        }
      }
      else if (FLAG_DESIGN_ENZ || FLAG_DESIGN_PLI)
      {
        Residue *pSmallMol = NULL;
        StructureFindSmallMol(pStruct, &pSmallMol);
        if (AtomArrayCalcMinDistance(&pResIJ->atoms, &pSmallMol->atoms) < ENERGY_DISTANCE_CUTOFF)
        {
          ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStruct, i, j, pLib, pParams, pTops);
          if (FLAG_USE_INPUT_SC && pResIJ->isSCIntact)
            ProteinSiteBuildNativeRotamer(pStruct, i, j, pTops);
          if (FLAG_ROTATE_HYDROXYL)
            ProteinSiteExpandHydroxylRotamers(pStruct, i, j, pTops);
          ProteinSiteCalcRotamersEnergy(pStruct, pAAPP, pRama, i, j, pFile);
          ProteinSiteRemoveDesignSite(pStruct, i, j);
        }
      }
      else
      {
        ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStruct, i, j, pLib, pParams, pTops);
        if (FLAG_USE_INPUT_SC && pResIJ->isSCIntact)
          ProteinSiteBuildNativeRotamer(pStruct, i, j, pTops);
        if (FLAG_ROTATE_HYDROXYL)
          ProteinSiteExpandHydroxylRotamers(pStruct, i, j, pTops);
        ProteinSiteCalcRotamersEnergy(pStruct, pAAPP, pRama, i, j, pFile);
        ProteinSiteRemoveDesignSite(pStruct, i, j);
      }
    }
  }
  fclose(pFile);

  return Success;
}

int CheckRotamerInBBdepRotLib(Structure *pStruct, BBdepRotLib *pLib, ResTopos *pTops, double cutoff)
{
  FILE *pFile = fopen("UD_torsion_recovery.txt", "w");
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChain = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChain); j++)
    {
      Residue *pRes = ChainGetResidue(pChain, j);
      if (strcmp(ResidueGetName(pRes), "ALA") != 0 && strcmp(ResidueGetName(pRes), "GLY") != 0)
      {
        BOOL result = IsNativeRotamerInBBdepRotLib(pStruct, i, j, pTops, pLib, cutoff);
        if (result)
          fprintf(pFile, "%s %s 1\n", ChainGetName(pChain), ResidueGetName(pRes));
        else
          fprintf(pFile, "%s %s 0\n", ChainGetName(pChain), ResidueGetName(pRes));
      }
      else
        fprintf(pFile, "%s %s 1\n", ChainGetName(pChain), ResidueGetName(pRes));
    }
  }
  fclose(pFile);
  return Success;
}

int FindMinRmsdRotFromRotLib(Structure *pStruct)
{
  FILE *pFile = fopen("UD_MinRMSDRotamers.txt", "w");
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChain = StructureGetChain(pStruct, i);
    for (int j = 0; j < ChainGetResidueNum(pChain); j++)
    {
      Residue *pRes = ChainGetResidue(pChain, j);
      DesignSite *pSite = StructureFindDesignSite(pStruct, i, j);
      if (pSite != NULL)
      {
        double minRMSD = 1e8;
        Rotamers *pSet = DesignSiteGetRotamers(pSite);
        for (int k = 0; k < RotamersGetNum(pSet); k++)
        {
          Rotamer *pRot = RotamersGet(pSet, k);
          RotamerRestore(pRot, pSet);
          double rmsd = RotamerAndResidueSidechainRMSD(pRot, pRes);
          if (RotamerIsSymmetricalCheck(pRot))
          {
            Rotamer tempRot;
            RotamerCreate(&tempRot);
            SymmetricalRotamerGenerate(&tempRot, pRot);
            double rmsd2 = RotamerAndResidueSidechainRMSD(&tempRot, pRes);
            if (rmsd2 < rmsd)
              rmsd = rmsd2;
            RotamerDestroy(&tempRot);
          }
          if (rmsd < minRMSD)
            minRMSD = rmsd;
          RotamerExtract(pRot);
        }
        fprintf(pFile, "%s %s %f\n", ChainGetName(pChain), ResidueGetName(pRes), minRMSD);
      }
    }
  }
  fclose(pFile);
  return Success;
}

int CompareTwoStructureProteinSidechains(Structure *pStruct, Structure *pStruct2, FILE *pTorsion, FILE *pRmsd)
{
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChain = StructureGetChain(pStruct, i);
    Chain *pChain2 = StructureGetChain(pStruct2, i);
    for (int j = 0; j < ChainGetResidueNum(pChain); j++)
    {
      Residue *pRes = ChainGetResidue(pChain, j);
      Residue *pRes2 = ChainGetResidue(pChain2, j);
      char res = AA3ToAA1(ResidueGetName(pRes));
      char res2 = AA3ToAA1(ResidueGetName(pRes2));
      if (ChainGetType(pChain) == Type_Chain_Protein)
      {
        if (res != 'A' && res != 'G')
        {
          double rmsd = ResidueAndResidueSidechainRMSD(pRes, pRes2);
          if (ResidueIsSymmetricalCheck(pRes))
          {
            Residue tempRes;
            ResidueCreate(&tempRes);
            SymmetricalResidueGenerate(&tempRes, pRes);
            double rmsd2 = ResidueAndResidueSidechainRMSD(&tempRes, pRes2);
            if (rmsd2 < rmsd)
              rmsd = rmsd2;
            ResidueDestroy(&tempRes);
          }
          fprintf(pRmsd, "%s %s %f\n", ChainGetName(pChain), ResidueGetName(pRes), rmsd);
          IntArray similarityArray;
          IntArrayCreate(&similarityArray, 0);
          fprintf(pTorsion, "%s %s", ChainGetName(pChain), ResidueGetName(pRes));
          if (ResidueAndResidueCheckTorsionSimilarity(pRes, pRes2, TORSION_DEVIATION_THRESHOLD, &similarityArray))
            fprintf(pTorsion, " 1");
          else
            fprintf(pTorsion, " 0");
          // calculate torsion angle-wise similarity
          for (int k = 0; k < IntArrayGetLength(&similarityArray); k++)
          {
            if (IntArrayGet(&similarityArray, k) == 1)
              fprintf(pTorsion, " 1");
            else if (IntArrayGet(&similarityArray, k) == 2)
              fprintf(pTorsion, " 0");
          }
          fprintf(pTorsion, "\n");
          IntArrayDestroy(&similarityArray);
        }
        else
        {
          double rmsd = 0.0;
          fprintf(pRmsd, "%s %s %f\n", ChainGetName(pChain), ResidueGetName(pRes), rmsd);
          fprintf(pTorsion, "%s %s 1\n", ChainGetName(pChain), ResidueGetName(pRes));
        }
      }
    }
  }
  return Success;
}

int CheckClash0(Structure *pStructure, double clashRatio)
{
  printf("checking structure clashes with atom pairwise distance < %.2f*(ri+rj)\n", clashRatio);
  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChainI = StructureGetChain(pStructure, i);

    for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
    {
      Residue *pResiIJ = ChainGetResidue(pChainI, j);
      for (int j2 = j + 1; j2 < ChainGetResidueNum(pChainI); j2++)
      {
        Residue *pResiIJ2 = ChainGetResidue(pChainI, j2);
        for (int atom1 = 0; atom1 < ResidueGetAtomNum(pResiIJ); atom1++)
        {
          Atom *pAtom1 = ResidueGetAtom(pResiIJ, atom1);
          if (pAtom1->isBBAtom || AtomIsHydrogen(pAtom1))
            continue;
          for (int atom2 = 0; atom2 < ResidueGetAtomNum(pResiIJ2); atom2++)
          {
            Atom *pAtom2 = ResidueGetAtom(pResiIJ2, atom2);
            if (pAtom2->isBBAtom || AtomIsHydrogen(pAtom2))
              continue;
            double dist = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
            if (dist < clashRatio * (pAtom1->vdw_rad + pAtom2->vdw_rad))
            {
              printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                     AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                     AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiIJ2), AtomGetName(pAtom2),
                     dist);
            }
          }
        }
      }
    }

    for (int k = i + 1; k < StructureGetChainNum(pStructure); k++)
    {
      Chain *pChainK = StructureGetChain(pStructure, k);
      for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
      {
        Residue *pResiIJ = ChainGetResidue(pChainI, j);
        for (int l = 0; l < ChainGetResidueNum(pChainK); l++)
        {
          Residue *pResiKL = ChainGetResidue(pChainK, l);
          for (int atom1 = 0; atom1 < ResidueGetAtomNum(pResiIJ); atom1++)
          {
            Atom *pAtom1 = ResidueGetAtom(pResiIJ, atom1);
            if (pAtom1->isBBAtom || AtomIsHydrogen(pAtom1))
              continue;
            for (int atom2 = 0; atom2 < ResidueGetAtomNum(pResiKL); atom2++)
            {
              Atom *pAtom2 = ResidueGetAtom(pResiKL, atom2);
              if (pAtom2->isBBAtom || AtomIsHydrogen(pAtom2))
                continue;
              double dist = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
              if (dist < clashRatio * (pAtom1->vdw_rad + pAtom2->vdw_rad))
              {
                if (ChainGetType(pChainI) == Type_Chain_SmallMol)
                {
                  printf("ligand residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                         AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                         AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiKL), AtomGetName(pAtom2),
                         dist);
                }
                else if (ChainGetType(pChainK) == Type_Chain_SmallMol)
                {
                  printf("protein residue %s%d%s atom %s and ligand residue %s%d%s atom %s has clash, dist=%f\n",
                         AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                         AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiKL), AtomGetName(pAtom2),
                         dist);
                }
                else
                {
                  printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                         AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                         AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiKL), AtomGetName(pAtom2),
                         dist);
                }
              }
            }
          }
        }
      }
    }
  }

  return Success;
}

int CheckClash1(Structure *pStructure, double clashRatio)
{
  printf("checking structure clashes with atom pairwise distance < %.2f*(ri+rj)\n", clashRatio);
  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChainI = StructureGetChain(pStructure, i);

    for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
    {
      Residue *pResiIJ = ChainGetResidue(pChainI, j);
      for (int j2 = j + 1; j2 < ChainGetResidueNum(pChainI); j2++)
      {
        Residue *pResiIJ2 = ChainGetResidue(pChainI, j2);
        for (int atom1 = 0; atom1 < ResidueGetAtomNum(pResiIJ); atom1++)
        {
          Atom *pAtom1 = ResidueGetAtom(pResiIJ, atom1);
          if (pAtom1->isBBAtom || AtomIsHydrogen(pAtom1))
            continue;
          for (int atom2 = 0; atom2 < ResidueGetAtomNum(pResiIJ2); atom2++)
          {
            Atom *pAtom2 = ResidueGetAtom(pResiIJ2, atom2);
            if (pAtom2->isBBAtom || AtomIsHydrogen(pAtom2))
              continue;
            double dist = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
            if (dist < clashRatio * (pAtom1->vdw_rad + pAtom2->vdw_rad))
            {
              printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                     AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                     AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiIJ2), AtomGetName(pAtom2),
                     dist);
            }
          }
        }
      }
    }

    for (int k = i + 1; k < StructureGetChainNum(pStructure); k++)
    {
      Chain *pChainK = StructureGetChain(pStructure, k);
      for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
      {
        Residue *pResiIJ = ChainGetResidue(pChainI, j);
        for (int l = 0; l < ChainGetResidueNum(pChainK); l++)
        {
          Residue *pResiKL = ChainGetResidue(pChainK, l);
          for (int atom1 = 0; atom1 < ResidueGetAtomNum(pResiIJ); atom1++)
          {
            Atom *pAtom1 = ResidueGetAtom(pResiIJ, atom1);
            if (pAtom1->isBBAtom || AtomIsHydrogen(pAtom1))
              continue;
            for (int atom2 = 0; atom2 < ResidueGetAtomNum(pResiKL); atom2++)
            {
              Atom *pAtom2 = ResidueGetAtom(pResiKL, atom2);
              if (pAtom2->isBBAtom || AtomIsHydrogen(pAtom2))
                continue;
              double dist = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
              if (dist < clashRatio * (pAtom1->vdw_rad + pAtom2->vdw_rad))
              {
                if (ChainGetType(pChainI) == Type_Chain_SmallMol)
                {
                  printf("ligand residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                         AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                         AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiKL), AtomGetName(pAtom2),
                         dist);
                }
                else if (ChainGetType(pChainK) == Type_Chain_SmallMol)
                {
                  printf("protein residue %s%d%s atom %s and ligand residue %s%d%s atom %s has clash, dist=%f\n",
                         AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                         AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiKL), AtomGetName(pAtom2),
                         dist);
                }
                else
                {
                  printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                         AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                         AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiKL), AtomGetName(pAtom2),
                         dist);
                }
              }
            }
          }
        }
      }
    }
  }

  return Success;
}

int CheckClash2(Structure *pStructure, double clashRatio)
{
  printf("checking structure clashes with atom pairwise distance < %.2f*(ri+rj)\n", clashRatio);
  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChainI = StructureGetChain(pStructure, i);

    for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
    {
      Residue *pResiIJ = ChainGetResidue(pChainI, j);
      for (int j2 = j + 1; j2 < ChainGetResidueNum(pChainI); j2++)
      {
        Residue *pResiIJ2 = ChainGetResidue(pChainI, j2);
        for (int atom1 = 0; atom1 < ResidueGetAtomNum(pResiIJ); atom1++)
        {
          Atom *pAtom1 = ResidueGetAtom(pResiIJ, atom1);
          if (pAtom1->isBBAtom || AtomIsHydrogen(pAtom1))
            continue;
          for (int atom2 = 0; atom2 < ResidueGetAtomNum(pResiIJ2); atom2++)
          {
            Atom *pAtom2 = ResidueGetAtom(pResiIJ2, atom2);
            if (pAtom2->isBBAtom || AtomIsHydrogen(pAtom2))
              continue;
            double dist = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
            if (dist < clashRatio * (pAtom1->vdw_rad + pAtom2->vdw_rad))
            {
              printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                     AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                     AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiIJ2), AtomGetName(pAtom2),
                     dist);
            }
          }
        }
      }
    }

    for (int k = i + 1; k < StructureGetChainNum(pStructure); k++)
    {
      Chain *pChainK = StructureGetChain(pStructure, k);
      for (int j = 0; j < ChainGetResidueNum(pChainI); j++)
      {
        Residue *pResiIJ = ChainGetResidue(pChainI, j);
        for (int l = 0; l < ChainGetResidueNum(pChainK); l++)
        {
          Residue *pResiKL = ChainGetResidue(pChainK, l);
          for (int atom1 = 0; atom1 < ResidueGetAtomNum(pResiIJ); atom1++)
          {
            Atom *pAtom1 = ResidueGetAtom(pResiIJ, atom1);
            if (pAtom1->isBBAtom || AtomIsHydrogen(pAtom1))
              continue;
            for (int atom2 = 0; atom2 < ResidueGetAtomNum(pResiKL); atom2++)
            {
              Atom *pAtom2 = ResidueGetAtom(pResiKL, atom2);
              if (pAtom2->isBBAtom || AtomIsHydrogen(pAtom2))
                continue;
              double dist = XYZDist(&pAtom1->xyz, &pAtom2->xyz);
              if (dist < clashRatio * (pAtom1->vdw_rad + pAtom2->vdw_rad))
              {
                if (ChainGetType(pChainI) == Type_Chain_SmallMol)
                {
                  printf("ligand residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                         AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                         AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiKL), AtomGetName(pAtom2),
                         dist);
                }
                else if (ChainGetType(pChainK) == Type_Chain_SmallMol)
                {
                  printf("protein residue %s%d%s atom %s and ligand residue %s%d%s atom %s has clash, dist=%f\n",
                         AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                         AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiKL), AtomGetName(pAtom2),
                         dist);
                }
                else
                {
                  printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                         AtomGetChainName(pAtom1), AtomGetPosInChain(pAtom1), ResidueGetName(pResiIJ), AtomGetName(pAtom1),
                         AtomGetChainName(pAtom2), AtomGetPosInChain(pAtom2), ResidueGetName(pResiKL), AtomGetName(pAtom2),
                         dist);
                }
              }
            }
          }
        }
      }
    }
  }

  return Success;
}

int ComputeAminoAcidComposition(Structure *pStructure, int *aas)
{
  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChain = StructureGetChain(pStructure, i);
    for (int j = 0; j < ChainGetResidueNum(pChain); j++)
    {
      Residue *pResi = ChainGetResidue(pChain, j);
      if (strcmp(ResidueGetName(pResi), "ALA") == 0)
        aas[0]++;
      else if (strcmp(ResidueGetName(pResi), "CYS") == 0)
        aas[1]++;
      else if (strcmp(ResidueGetName(pResi), "ASP") == 0)
        aas[2]++;
      else if (strcmp(ResidueGetName(pResi), "GLU") == 0)
        aas[3]++;
      else if (strcmp(ResidueGetName(pResi), "PHE") == 0)
        aas[4]++;
      else if (strcmp(ResidueGetName(pResi), "GLY") == 0)
        aas[5]++;
      else if (strcmp(ResidueGetName(pResi), "HSD") == 0)
        aas[6]++;
      else if (strcmp(ResidueGetName(pResi), "HSE") == 0)
        aas[6]++;
      else if (strcmp(ResidueGetName(pResi), "HIS") == 0)
        aas[6]++;
      else if (strcmp(ResidueGetName(pResi), "ILE") == 0)
        aas[7]++;
      else if (strcmp(ResidueGetName(pResi), "LYS") == 0)
        aas[8]++;
      else if (strcmp(ResidueGetName(pResi), "LEU") == 0)
        aas[9]++;
      else if (strcmp(ResidueGetName(pResi), "MET") == 0)
        aas[10]++;
      else if (strcmp(ResidueGetName(pResi), "ASN") == 0)
        aas[11]++;
      else if (strcmp(ResidueGetName(pResi), "PRO") == 0)
        aas[12]++;
      else if (strcmp(ResidueGetName(pResi), "GLN") == 0)
        aas[13]++;
      else if (strcmp(ResidueGetName(pResi), "ARG") == 0)
        aas[14]++;
      else if (strcmp(ResidueGetName(pResi), "SER") == 0)
        aas[15]++;
      else if (strcmp(ResidueGetName(pResi), "THR") == 0)
        aas[16]++;
      else if (strcmp(ResidueGetName(pResi), "VAL") == 0)
        aas[17]++;
      else if (strcmp(ResidueGetName(pResi), "TRP") == 0)
        aas[18]++;
      else if (strcmp(ResidueGetName(pResi), "TYR") == 0)
        aas[19]++;
      else if (strcmp(ResidueGetName(pResi), "A") == 0)
        aas[20]++;
      else if (strcmp(ResidueGetName(pResi), "C") == 0)
        aas[21]++;
      else if (strcmp(ResidueGetName(pResi), "G") == 0)
        aas[22]++;
      else if (strcmp(ResidueGetName(pResi), "U") == 0)
        aas[23]++;
      else if (strcmp(ResidueGetName(pResi), "DA") == 0)
        aas[24]++;
      else if (strcmp(ResidueGetName(pResi), "DC") == 0)
        aas[25]++;
      else if (strcmp(ResidueGetName(pResi), "DG") == 0)
        aas[26]++;
      else if (strcmp(ResidueGetName(pResi), "DT") == 0)
        aas[27]++;
    }
  }
  printf("\nResidue composition of the structure:\n");
  printf("ALA     = %4d\n", aas[0]);
  printf("CYS     = %4d\n", aas[1]);
  printf("ASP     = %4d\n", aas[2]);
  printf("GLU     = %4d\n", aas[3]);
  printf("PHE     = %4d\n", aas[4]);
  printf("GLY     = %4d\n", aas[5]);
  printf("HIS     = %4d\n", aas[6]);
  printf("ILE     = %4d\n", aas[7]);
  printf("LYS     = %4d\n", aas[8]);
  printf("LEU     = %4d\n", aas[9]);
  printf("MET     = %4d\n", aas[10]);
  printf("ASN     = %4d\n", aas[11]);
  printf("PRO     = %4d\n", aas[12]);
  printf("GLN     = %4d\n", aas[13]);
  printf("ARG     = %4d\n", aas[14]);
  printf("SER     = %4d\n", aas[15]);
  printf("THR     = %4d\n", aas[16]);
  printf("VAL     = %4d\n", aas[17]);
  printf("TRP     = %4d\n", aas[18]);
  printf("TYR     = %4d\n", aas[19]);
  printf("A(RNA)  = %4d\n", aas[20]);
  printf("C(RNA)  = %4d\n", aas[21]);
  printf("G(RNA)  = %4d\n", aas[22]);
  printf("U(RNA)  = %4d\n", aas[23]);
  printf("DA(DNA) = %4d\n", aas[24]);
  printf("DC(DNA) = %4d\n", aas[25]);
  printf("DG(DNA) = %4d\n", aas[26]);
  printf("DT(DNA) = %4d\n", aas[27]);

  return Success;
}

int SelectResiduesInRange(Structure *pStructure)
{
  StringArray chain_and_positions;
  StringArrayCreate(&chain_and_positions);
  StringArraySplitString(&chain_and_positions, REFERENCE_RESIDUES, ',');

  Residue **pResiduesWithin = (Residue **)malloc(sizeof(Residue *) * StringArrayGetNum(&chain_and_positions));
  for (int res = 0; res < StringArrayGetNum(&chain_and_positions); res++)
  {
    char chName[MAX_LEN_CHAIN_NAME + 1];
    int pos;
    sscanf(StringArrayGet(&chain_and_positions, res), "%c%d", &chName[0], &pos);
    chName[1] = '\0';
    Chain *pChain = StructureFindChainByName(pStructure, chName);
    int resNdx = -1;
    ChainFindResidueIndexByPosInChain(pChain, pos, &resNdx);
    if (resNdx != -1)
    {
      pResiduesWithin[res] = ChainGetResidue(pChain, resNdx);
    }
  }

  for (int i = 0; i < StructureGetChainNum(pStructure); i++)
  {
    Chain *pChain = StructureGetChain(pStructure, i);
    for (int j = 0; j < ChainGetResidueNum(pChain); j++)
    {
      Residue *pResi = ChainGetResidue(pChain, j);
      BOOL resIsAwithinRes = FALSE;
      BOOL resIsWithinRange = FALSE;
      for (int k = 0; k < StringArrayGetNum(&chain_and_positions); k++)
      {
        Residue *pResiWith = pResiduesWithin[k];
        if (pResi == pResiWith)
        {
          resIsAwithinRes = TRUE;
          break;
        }
        else
        {
          if (AtomArrayCalcMinDistance(&pResi->atoms, &pResiWith->atoms) < DISTANCE_TO_REFERENCE_RESIDUES)
          {
            resIsWithinRange = TRUE;
            break;
          }
        }
      }
      if (resIsAwithinRes == FALSE && resIsWithinRange == TRUE)
      {
        printf("%s %4d %c\n", ResidueGetChainName(pResi), ResidueGetPosInChain(pResi), AA3ToAA1(ResidueGetName(pResi)));
      }
    }
  }

  return Success;
}

int StructureFindMutations(Structure *pStruct, Structure *pStruct2, char *mutStr)
{
  char strs[2 * MAX_LEN_FILE_NAME + 1] = "";
  BOOL firstMut = TRUE;
  for (int i = 0; i < StructureGetChainNum(pStruct); i++)
  {
    Chain *pChain = StructureGetChain(pStruct, i);
    Chain *pChain2 = StructureGetChain(pStruct2, i);
    if (ChainGetType(pChain) == Type_Chain_Protein)
    {
      for (int j = 0; j < ChainGetResidueNum(pChain); j++)
      {
        Residue *pRes = ChainGetResidue(pChain, j);
        Residue *pRes2 = ChainGetResidue(pChain2, j);
        char res = AA3ToAA1(ResidueGetName(pRes));
        char res2 = AA3ToAA1(ResidueGetName(pRes2));
        if (res != res2)
        {
          if (firstMut == TRUE)
          {
            sprintf(strs, "%s%c%c%d%c", strs, res, ResidueGetChainName(pRes)[0], ResidueGetPosInChain(pRes), res2);
            firstMut = FALSE;
          }
          else
            sprintf(strs, "%s,%c%c%d%c", strs, res, ResidueGetChainName(pRes)[0], ResidueGetPosInChain(pRes), res2);
        }
      }
    }
  }
  strcpy(mutStr, strs);
  return Success;
}