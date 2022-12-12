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

#pragma warning(disable:4996)
#pragma warning(disable:6031)
#pragma warning(disable:26812)
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "Getopt.h"
#include "ProgramFunction.h"
#include "RotamerBuilder.h"
#include "Structure.h"
//#include "EnergyMatrix.h"
#include "Evolution.h"
//#include "DEE.h"
#include "ProteinDesign.h"
#include "WeightOpt.h"
#include "SmallMolEEF1.h"
#include "SmallMolParAndTopo.h"
#include "EvoAminoName.h"
#include "ProgramPreprocess.h"

#define  PROGRAM_PARAMETERS

char PROGRAM_PATH[MAX_LEN_ONE_LINE_CONTENT + 1] = ".";
char PROGRAM_NAME[MAX_LEN_FILE_NAME + 1] = "UniDesign";
char PROGRAM_VERSION[MAX_LEN_FILE_NAME + 1] = "1.0";

//////////////////////////////////////////////////////////////
// 1. PARAMETER, TOPOLOGY and ROTAMER LIBRARY
//////////////////////////////////////////////////////////////
// Atom Parameter abd Residue Topology
char FILE_ATOMPAR[MAX_LEN_FILE_NAME + 1] = "library/toppar/param_charmm19_lk.prm";
char FILE_TOPO[MAX_LEN_FILE_NAME + 1] = "library/toppar/top_polh19.inp";
char LIG_ATOMPAR[MAX_LEN_FILE_NAME + 1] = "LIG_PARAM.prm";
char LIG_TOPFILE[MAX_LEN_FILE_NAME + 1] = "LIG_TOPO.inp";
// Statistical Energy Terms and Weights
char FILE_AAPROPENSITY[MAX_LEN_FILE_NAME + 1] = "library/eterms/aapropensity.nrg";
char FILE_RAMACHANDRAN[MAX_LEN_FILE_NAME + 1] = "library/eterms/ramachandran.nrg";
char FILE_WEIGHT_READ[MAX_LEN_FILE_NAME + 1] = "wread/weight_all1.wgt";
// Rotamer Library, includes:
// Backbone-dependent library: dun2010bb.lib, dun2010bb1per.lib, and dun2010bb3per.lib in TEXT format or ALLbbdep.bin in BINARY format
// Backbone-independent library: honig984lib, honig3222.lib, honig7421.lib, and honig11810.lib in TEXT format
char FILE_ROTLIB[MAX_LEN_FILE_NAME + 1] = "library/rotlib/dun2010bb3per.lib";
// By default, the ALLbbdep.bin is used with a rotamer probability cutoff of <=0.03
char FILE_ROTLIB_BIN[MAX_LEN_FILE_NAME + 1] = "library/rotlib/ALLbbdep.bin";
double CUT_EXCL_LOW_ROT_PROB = 0.03; // 0.01~0.03 suggested
// User can select a rotamer library
BOOL FLAG_USER_ROTLIB = FALSE;
char USER_ROTLIB_NAME[MAX_LEN_FILE_NAME + 1] = "";

// Evolution and PSSM
char TGT_PRF[MAX_LEN_FILE_NAME + 1] = "prf.txt"; // Positiion-Specific Scoring Matrix file
char TGT_MSA[MAX_LEN_FILE_NAME + 1] = "msa.txt"; // Multiple Sequence Alignment file
char TGT_SA[MAX_LEN_FILE_NAME + 1] = "sa.txt";  // Solvent Accessibility file
char TGT_SS[MAX_LEN_FILE_NAME + 1] = "ss.txt";  // Secondary Structure file
char TGT_SEQ[MAX_LEN_FILE_NAME + 1] = "seq.txt"; // Sequence file in single-line plain text
char TGT_PHIPSI[MAX_LEN_FILE_NAME + 1] = "phi-psi.txt"; // Main-chain Phi and Psi angles for input backbone
//store profile into memory
float** PROT_PROFILE = NULL;      // Store PSSM into a 2-D float** array
double WGT_PROFILE = 1.0;       // Energy Weight of Evolutionary Profile
double WGT_BIND = 1.0;       // Energy Weight of Binding interaction (default: 1.0)

// for ENZYME design
char FILE_CATACONS[MAX_LEN_FILE_NAME + 1] = "LIG_CATACONS.txt";
char FILE_LIG_PLACEMENT[MAX_LEN_FILE_NAME + 1] = "LIG_PLACEMENT.txt";

// protein design OUTPUT
char FILE_SELF_ENERGY[MAX_LEN_FILE_NAME + 1] = "selfenergy.txt";
char ROT_LIST_ALL[MAX_LEN_FILE_NAME + 1] = "rotlist.txt";
char ROT_LIST_SEC[MAX_LEN_FILE_NAME + 1] = "rotlistSEC.txt";
char ROT_INDEX_FILE[MAX_LEN_FILE_NAME + 1] = "desrots";
char SEQ_FILE[MAX_LEN_FILE_NAME + 1] = "desseqs";
char BEST_SEQ[MAX_LEN_FILE_NAME + 1] = "bestseq";
char BEST_STRUCT[MAX_LEN_FILE_NAME + 1] = "beststruct";
char BEST_DESSITE[MAX_LEN_FILE_NAME + 1] = "bestsites";
char BEST_MOL2[MAX_LEN_FILE_NAME + 1] = "bestlig";


///////////////////////////////////////////////////////
// 2. INPUT
///////////////////////////////////////////////////////
// Pdb
char PDB[MAX_LEN_FILE_NAME + 1] = "2pjyAB_aligned.pdb";
char PDBPATH[MAX_LEN_FILE_NAME + 1] = ".";
char PDBNAME[MAX_LEN_FILE_NAME + 1] = "pdbname.pdb";
char PDBID[MAX_LEN_FILE_NAME + 1] = "pdb";
// Mol2 for ligand
char MOL2[MAX_LEN_FILE_NAME + 1] = "ligand.mol2";

///////////////////////////////////////////////////////
// Parameters for ProteinDesign
char DES_CHAINS[MAX_LEN_ONE_LINE_CONTENT + 1] = "AE";
int  NTRAJ = 1;
int  NTRAJ_START_NDX = 1;
// parameters for PPI design
double CUT_PPI_DIST_SHELL1 = 5.0;
double CUT_PPI_DIST_SHELL2 = 8.0;
// parameters for ENZYME and PROTEIN-LIGAND-INTERACTION design
double CUT_PLI_DIST_SHELL1 = 5.0;
double CUT_PLI_DIST_SHELL2 = 8.0;
// Set default desType
Type_ResidueDesignType DEFAULT_RESIDUE_DESIGNTYPE = Type_DesType_Fixed;

// Parameters for EnzymeDesign
// user should specify three initial atoms for creating small-molecule topology
char INI_ATOM1[MAX_LEN_ATOM_NAME + 1] = "O39";
char INI_ATOM2[MAX_LEN_ATOM_NAME + 1] = "C20";
char INI_ATOM3[MAX_LEN_ATOM_NAME + 1] = "S17";

// parameters for ComputeBinding
BOOL FLAG_CHAIN_SPLIT = TRUE;
char SPLIT_CHAINS[MAX_LEN_FILE_NAME + 1] = "AB,C";
char SPLIT_PART1[MAX_LEN_FILE_NAME + 1] = "AB";
char SPLIT_PART2[MAX_LEN_FILE_NAME + 1] = "C";

// parameters for ComputeResiEnergy
char RESI[MAX_LEN_CHAIN_NAME + 1] = "C8";
char EXCLUDE_RESI[MAX_LEN_ONE_LINE_CONTENT + 1] = "";

// parameters for ComputeResiPairEnergy
char RESI_PAIR[2 * MAX_LEN_CHAIN_NAME + 1] = "";

// Parameters for BuildMutant
char MUTANT_FILE[MAX_LEN_FILE_NAME + 1] = "mutant_file.txt";
int MAX_NUM_OF_RUNS = 4; // default: 4 runs

// parameters for CheckResiInRotLib
double CUT_TORSION_DEVIATION = 20.0;

// parameters for CheckClash
double CUT_CLASH_RATIO = 0.6;

// parameters for FindCoreResidue (n.CB>=20)
int CUT_NUM_CB_CORE = 20;
// parameters for FindSurfaceResidue (n.CB<=15)
int CUT_NUM_CB_SURF = 15;

// parameters for MakeLigConformers, ScreenLigConformers, and ENZYME design
BOOL FLAG_LIG_CONFORMERS = TRUE;
char FILE_LIG_CONFORMERS_IN[MAX_LEN_FILE_NAME + 1] = "LIG_CONFORMERS.pdb";
char FILE_LIG_CONFORMERS_OUT[MAX_LEN_FILE_NAME + 1] = "LIG_CONFORMERS2.pdb";
BOOL FLAG_LIG_SCREEN_BY_ORIENTATION = FALSE;
char LIGSCREEN_ORITENTATION_FILE[MAX_LEN_ATOM_NAME + 1] = "";
BOOL FLAG_LIG_SCREEN_BY_TOPVDW = FALSE;
double LIGSCREEN_TOPVDW_PERCENT = 0.5;
BOOL FLAG_LIG_SCREEN_BY_RMSD = TRUE;
double LIGSCREEN_RMSD_CUTOFF = 0.5;

// Parameters for CompareSidechain
char PDB2[MAX_LEN_FILE_NAME + 1] = "1aho2.pdb";

// Parameters for OptimizeWeight
char PDBLIST[MAX_LEN_FILE_NAME + 1] = "pdblist.txt";

// Parameters for energy score normalization
int PROT_LEN_NORM = 100;

#define PROGRAM_FLAGS

BOOL FLAG_PDB = TRUE;
BOOL FLAG_MOL2 = FALSE;

BOOL FLAG_MONOMER = TRUE;
BOOL FLAG_PPI = FALSE;
BOOL FLAG_PROT_LIG = FALSE;
BOOL FLAG_ENZYME = FALSE;

// score function (default: physics)
BOOL FLAG_PHYSICS = TRUE;
BOOL FLAG_EVOLUTION = FALSE;
BOOL FLAG_EVOPHIPSI = FALSE;

// restrictions for ProteinDesign
BOOL FLAG_BBDEP_ROTLIB = TRUE;  // specify a rotamer library for side-chain sampling (default: Dunbrack 2010 bb-dep rotlib)
BOOL FLAG_USE_INPUT_SC = TRUE;  // input structure's side chains used as rots (default: Yes)
BOOL FLAG_ROTATE_HYDROXYL = TRUE;  // rotate hydrogens for Ser, Thr, and Tyr (default: Yes)
BOOL FLAG_WILDTYPE_ONLY = FALSE; // set optional amino-acid types as wild-type (equivalent to side-chain repacking)
BOOL FLAG_INTERFACE_ONLY = FALSE; // design interface residues only
BOOL FLAG_EXCL_CYS_ROTS = FALSE; // Cysteine rots excluded (default: No)
// designate design sites and optional amino-acid types using a resfile
char FILE_RESFILE[MAX_LEN_FILE_NAME + 1] = "RESFILE.txt";
BOOL FLAG_RESFILE = TRUE;

// start a protein design job from the native sequence (input from the given PDB)
BOOL FLAG_DESIGN_FROM_NATAA = FALSE; // default: start from random; thus this flag is set to FALSE

// flag for reading & writing hydrogen atoms
BOOL FLAG_READ_HYDROGEN = TRUE;
BOOL FLAG_WRITE_HYDROGEN = TRUE;

// select residues within a range to a set of target residues
char WITHIN_RESIDUES[MAX_LEN_ONE_LINE_CONTENT + 1];
double WITHIN_RANGE = 9.0;

#define COMMANDS_AND_OPTIONS
const int N_CMD = 200;//max count of supported commands
char SUPPORTED_COMMANDS[][30] = {
  "ComputeBinding",
  "ComputeEvolutionScore",
  "ComputeResiEnergy",
  "ComputeRotamersEnergy",
  "ComputeStability",
  "ComputeResiPairEnergy",

  "ProteinDesign",

  "AddPolarHydrogen",
  "BuildMutant",
  "CalcPhiPsi",
  "CalcResiMinRMSD",
  "CheckClash0",
  "CheckClash1",
  "CheckClash2",
  "CheckResiInRotLib",
  "CompareSideChain",
  "FindCoreResidue",
  "FindIntermediateResidue",
  "FindInterfaceResidue",
  "FindSurfaceResidue",
  "Minimize",
  "OptimizeHydrogen",
  "RepairStructure",
  "ShowResiComposition",

  "PredPhiPsi",
  "PredSecStr",
  "PredSolAcc",

  "CalcResiBfactor",
  "OptimizeWeight",

  "MakeLigParamAndTopo",
  "MakeLigConformers",
  "AnalyzeLigConformers",
  "ScreenLigConformers",

  "WriteFirstLigRotIntoMol2",

  "SelectResiduesWithin",

  NULL,
};


//supported options
const char* SHORT_OPTS = "?hv";
struct option LONG_OPTS[] = {
  {"help",                no_argument,       NULL,  1},
  {"version",             no_argument,       NULL,  2},
  {"command",             required_argument, NULL,  3},
  {"pdb",                 required_argument, NULL,  4},
  {"split_chains",        required_argument, NULL,  5},
  {"mutant_file",         required_argument, NULL,  6},
  {"bbdep",               required_argument, NULL,  7},
  {"use_input_sc",        required_argument, NULL,  8},
  {"ppi_shell1",          required_argument, NULL,  9},
  {"ppi_shell2",          required_argument, NULL, 10},
  {"evolution",           no_argument,       NULL, 11},
  {"physics",             no_argument,       NULL, 12},
  {"monomer",             no_argument,       NULL, 13},
  {"ppint",               no_argument,       NULL, 14},
  {"protlig",             no_argument,       NULL, 15},
  {"enzyme",              no_argument,       NULL, 16},
  {"design_chains",       required_argument, NULL, 17},
  {"mol2",                required_argument, NULL, 18},
  {"rotate_hydroxyl",     required_argument, NULL, 19},
  {"xdeviation",          required_argument, NULL, 20},
  {"wread",               required_argument, NULL, 21},
  {"pdblist",             required_argument, NULL, 22},
  {"wildtype_only",       no_argument      , NULL, 23},
  {"rotlib",              required_argument, NULL, 24},
  {"pdb2",                required_argument, NULL, 25},
  {"pli_shell1",          required_argument, NULL, 26},
  {"pli_shell2",          required_argument, NULL, 27},
  {"clash_ratio",         required_argument, NULL, 28},
  {"ntraj",               required_argument, NULL, 29},
  {"excl_low_prob",       required_argument, NULL, 30},
  {"interface_only",      no_argument,       NULL, 32},
  {"seq",                 required_argument, NULL, 33},
  {"evo_allterms",        no_argument,       NULL, 34},
  {"wprof",               required_argument, NULL, 35},
  {"resfile",             required_argument, NULL, 36},
  {"seed_from_nat",       no_argument,       NULL, 38},
  {"excl_rot_cys",        no_argument,       NULL, 40},
  {"show_hydrogen",       required_argument, NULL, 41},
  {"ncut_cb_core",        required_argument, NULL, 42},
  {"ncut_cb_surf",        required_argument, NULL, 43},
  {"init_3atoms",         required_argument, NULL, 44}, // designate three initial atoms for generating ligand topology
  {"read_lig_conformers",   required_argument, NULL, 45}, // read ligand conformers from file
  {"write_lig_conformers",  required_argument, NULL, 46}, // write ligand conformers to file
  {"scrn_by_ornt",        required_argument, NULL, 47}, // screen ligand rots by orientation
  {"scrn_by_vdwpct",      required_argument, NULL, 48}, // screen ligand rots with both internalVDW and backboneVDW ranked in a low percentile, e.g. 25%
  {"scrn_by_rmsd",        required_argument, NULL, 49}, // screen ligand rots with an RMSD cutoff
  {"init_rotype",         required_argument, NULL, 50}, // initialize rotamer type for design (allaa, allaaxc, nataa, natrot)
  {"lig_atomparam",       required_argument, NULL, 51},
  {"lig_topology",        required_argument, NULL, 52},
  {"cata_cons",           required_argument, NULL, 53},
  {"num_of_runs",         required_argument, NULL, 54},
  {"within_residues",     required_argument, NULL, 55},
  {"within_range",        required_argument, NULL, 56},
  {"ntraj_start_ndx",     required_argument, NULL, 57},
  {"wbind",               required_argument, NULL, 58},
  {"resi",                required_argument, NULL, 59},
  {"resi_pair",           required_argument, NULL, 60},
  {"exclude",             required_argument, NULL, 61},
  {NULL,                  no_argument,       NULL,  0},
};


BOOL CheckCommandName(char* queryname)
{
  for (int i = 0; i < N_CMD; i++)
  {
    if (SUPPORTED_COMMANDS[i] == NULL) break;
    else if (!strcmp(queryname, SUPPORTED_COMMANDS[i])) return TRUE;
  }
  return FALSE;
}


int ParseOptions(int argc, char** argv)
{
  return getopt_long(argc, argv, SHORT_OPTS, LONG_OPTS, NULL);;
}


int main(int argc, char* argv[])
{
  char errMsg[MAX_LEN_ERR_MSG + 1];
  //char *cmdname = "ProteinDesign";
  //char* cmdname = "MakeLigParamAndTopo";
  //char* cmdname = "MakeLigConformers";
  //char* cmdname = "AnalyzeLigConformers";
  //char* cmdname = "ScreenLigConformers";
  //char* cmdname = "ComputeResiEnergy";
  //char* cmdname = "ComputeStability";
  //char* cmdname = "ComputeResiPairEnergy";
  char* cmdname = "FindInterfaceResidue";

  clock_t timestart = clock();
  setvbuf(stdout, NULL, _IONBF, 0);
  PrintAdvertisement();

  // set filepaths 
  ExtractPathAndName(argv[0], PROGRAM_PATH, PROGRAM_NAME);
  sprintf(FILE_ATOMPAR, "%s/library/toppar/param_charmm19_lk.prm", PROGRAM_PATH);
  sprintf(FILE_TOPO, "%s/library/toppar/top_polh19.inp", PROGRAM_PATH);
  sprintf(FILE_ROTLIB, "%s/library/rotlib/dun2010bb3per.lib", PROGRAM_PATH);
  sprintf(FILE_AAPROPENSITY, "%s/library/eterms/aapropensity.nrg", PROGRAM_PATH);
  sprintf(FILE_RAMACHANDRAN, "%s/library/eterms/ramachandran.nrg", PROGRAM_PATH);
  sprintf(FILE_WEIGHT_READ, "%s/wread/weight_all1.wgt", PROGRAM_PATH);

  while (TRUE)
  {
    int opt = ParseOptions(argc, argv);
    if (opt == -1) break;
    switch (opt)
    {
    case '?':
      PrintHelp();
      exit(Success);
    case 'h':
      PrintHelp();
      exit(Success);
    case 'v':
      PrintVersion();
      exit(Success);
    case 1:
      PrintHelp();
      exit(Success);
    case 2:
      PrintVersion();
      exit(Success);
    case 3:
      cmdname = optarg;
      if (!CheckCommandName(cmdname))
      {
        printf("command %s is not supported, program exits\n", cmdname);
        exit(ValueError);
      }
      else
      {
        printf("command %s is supported\n", cmdname);
      }
      break;
    case 4:
      FLAG_PDB = TRUE;
      strcpy(PDB, optarg);
      break;
    case 5:
      strcpy(SPLIT_CHAINS, optarg);
      sscanf(SPLIT_CHAINS, "%[^,],%s", SPLIT_PART1, SPLIT_PART2);
      FLAG_CHAIN_SPLIT = TRUE;
      // check if the two parts contain identical chains
      for (int i = 0; i < (int)strlen(SPLIT_PART1); i++)
      {
        char tmp[2] = { SPLIT_PART1[i],'\0' };
        if (strstr(SPLIT_PART2, tmp) != NULL)
        {
          printf("the two splitting parts should NOT have identical chains, program exits\n");
          exit(FormatError);
        }
      }
      break;
    case 6:
      strcpy(MUTANT_FILE, optarg);
      break;
    case 7:
      if (strcmp(optarg, "yes") == 0)
      {
        printf("backbone-dependent rotamer library enabled by user\n");
        FLAG_BBDEP_ROTLIB = TRUE;
      }
      else if (strcmp(optarg, "no") == 0)
      {
        printf("backbone-independent rotamer library enabled by user\n");
        FLAG_BBDEP_ROTLIB = FALSE;
      }
      else
      {
        printf("backbone-dependent rotamer library enabled by default\n");
        FLAG_BBDEP_ROTLIB = TRUE;
      }
      break;
    case 8:
      if (strcmp(optarg, "yes") == 0)
      {
        FLAG_USE_INPUT_SC = TRUE;
      }
      else if (strcmp(optarg, "no") == 0)
      {
        FLAG_USE_INPUT_SC = FALSE;
      }
      else
      {
        FLAG_USE_INPUT_SC = TRUE;
      }
      break;
    case 9:
      CUT_PPI_DIST_SHELL1 = atof(optarg);
      break;
    case 10:
      CUT_PPI_DIST_SHELL2 = atof(optarg);
      break;
    case 11:
      FLAG_EVOLUTION = TRUE;
      break;
    case 12:
      FLAG_PHYSICS = TRUE;
      break;
    case 13:
      FLAG_MONOMER = TRUE;
      FLAG_PPI = FALSE;
      FLAG_PROT_LIG = FALSE;
      FLAG_ENZYME = FALSE;
      break;
    case 14:
      FLAG_PPI = TRUE;
      FLAG_MONOMER = FALSE;
      FLAG_PROT_LIG = FALSE;
      FLAG_ENZYME = FALSE;
      break;
    case 15:
      FLAG_PROT_LIG = TRUE;
      FLAG_PPI = FALSE;
      FLAG_MONOMER = FALSE;
      FLAG_ENZYME = FALSE;
      break;
    case 16:
      FLAG_ENZYME = TRUE;
      FLAG_PROT_LIG = FALSE;
      FLAG_PPI = FALSE;
      FLAG_MONOMER = FALSE;
      break;
    case 17:
      strcpy(DES_CHAINS, optarg);
      break;
    case 18:
      strcpy(MOL2, optarg);
      FLAG_MOL2 = TRUE;
      break;
    case 19:
      if (strcmp(optarg, "yes") == 0) FLAG_ROTATE_HYDROXYL = TRUE;
      else if (strcmp(optarg, "no") == 0) FLAG_ROTATE_HYDROXYL = FALSE;
      else FLAG_ROTATE_HYDROXYL = TRUE;
      break;
    case 20:
      CUT_TORSION_DEVIATION = atof(optarg);
      break;
    case 21:
      strcpy(FILE_WEIGHT_READ, optarg);
      break;
    case 22:
      strcpy(PDBLIST, optarg);
      break;
    case 23:
      FLAG_WILDTYPE_ONLY = TRUE;
      break;
    case 24:
      FLAG_USER_ROTLIB = TRUE;
      strcpy(USER_ROTLIB_NAME, optarg);
      break;
    case 25:
      strcpy(PDB2, optarg);
      break;
    case 26:
      CUT_PLI_DIST_SHELL1 = atof(optarg);
      break;
    case 27:
      CUT_PLI_DIST_SHELL2 = atof(optarg);
      break;
    case 28:
      CUT_CLASH_RATIO = atof(optarg);
      break;
    case 29:
      NTRAJ = atoi(optarg);
      break;
    case 30:
      CUT_EXCL_LOW_ROT_PROB = atof(optarg);
      break;
    case 32:
      FLAG_INTERFACE_ONLY = TRUE;
      break;
    case 33:
      strcpy(TGT_SEQ, optarg);
      break;
    case 34:
      FLAG_EVOPHIPSI = TRUE;
      break;
    case 35:
      WGT_PROFILE = atof(optarg);
      break;
    case 36:
      strcpy(FILE_RESFILE, optarg);
      FLAG_RESFILE = TRUE;
      break;
    case 38:
      FLAG_DESIGN_FROM_NATAA = TRUE;
      break;
    case 40:
      FLAG_EXCL_CYS_ROTS = TRUE;
      break;
    case 41:
      if (strcmp(optarg, "yes") == 0)
      {
        FLAG_WRITE_HYDROGEN = TRUE;
      }
      else if (strcmp(optarg, "no") == 0)
      {
        FLAG_WRITE_HYDROGEN = FALSE;
      }
      else
      {
        FLAG_WRITE_HYDROGEN = TRUE;
      }
      break;
    case 42:
      CUT_NUM_CB_CORE = atoi(optarg);
      break;
    case 43:
      CUT_NUM_CB_SURF = atoi(optarg);
      break;
    case 44:
      StringArray atoms;
      StringArrayCreate(&atoms);
      StringArraySplitString(&atoms, optarg, ',');
      if (StringArrayGetCount(&atoms) < 3)
      {
        sprintf(errMsg, "in file %s line %d, please input three initial atoms for creating ligand topology, program exits", __FILE__, __LINE__);
        TraceError(errMsg, ValueError);
        exit(ValueError);
      }
      strcpy(INI_ATOM1, StringArrayGet(&atoms, 0));
      strcpy(INI_ATOM2, StringArrayGet(&atoms, 1));
      strcpy(INI_ATOM3, StringArrayGet(&atoms, 2));
      StringArrayDestroy(&atoms);
      break;
    case 45:
      FLAG_LIG_CONFORMERS = TRUE;
      strcpy(FILE_LIG_CONFORMERS_IN, optarg);
      break;
    case 46:
      strcpy(FILE_LIG_CONFORMERS_OUT, optarg);
      break;
    case 47:
      FLAG_LIG_SCREEN_BY_ORIENTATION = TRUE;
      strcpy(LIGSCREEN_ORITENTATION_FILE, optarg);
      break;
    case 48:
      FLAG_LIG_SCREEN_BY_TOPVDW = TRUE;
      LIGSCREEN_TOPVDW_PERCENT = atof(optarg);
      break;
    case 49:
      FLAG_LIG_SCREEN_BY_RMSD = TRUE;
      LIGSCREEN_RMSD_CUTOFF = atof(optarg);
      break;
    case 50:
      if (!strcmp(optarg, "natro")) DEFAULT_RESIDUE_DESIGNTYPE = Type_DesType_Fixed;
      else if (!strcmp(optarg, "nataa")) DEFAULT_RESIDUE_DESIGNTYPE = Type_DesType_Repackable;
      else if (!strcmp(optarg, "allaa")) DEFAULT_RESIDUE_DESIGNTYPE = Type_DesType_Mutable;
      else DEFAULT_RESIDUE_DESIGNTYPE = Type_DesType_Fixed;
      break;
    case 51:
      strcpy(LIG_ATOMPAR, optarg);
      break;
    case 52:
      strcpy(LIG_TOPFILE, optarg);
      break;
    case 53:
      strcpy(FILE_CATACONS, optarg);
      break;
    case 54:
      MAX_NUM_OF_RUNS = atoi(optarg);
      break;
    case 55:
      strcpy(WITHIN_RESIDUES, optarg);
      break;
    case 56:
      WITHIN_RANGE = atof(optarg);
      break;
    case 57:
      NTRAJ_START_NDX = atoi(optarg);
      break;
    case 58:
      WGT_BIND = atof(optarg);
      break;
    case 59:
      strcpy(RESI, optarg);
      break;
    case 60:
      strcpy(RESI_PAIR, optarg);
      break;
    case 61:
      strcpy(EXCLUDE_RESI, optarg);
      break;
    default:
      sprintf(errMsg, "in file %s line %d, unknown option, program exits", __FILE__, __LINE__);
      TraceError(errMsg, ValueError);
      exit(ValueError);
      break;
    }
  }

  // Early termination if no command detected
  if (!strcmp(cmdname, ""))
  {
    printf("no command name detected, please select a command for execution, program exits\n");
    printf("run './%s -h' to show help document\n", PROGRAM_NAME);
    exit(Success);
  }


  // quick check and file name settings
  if (FLAG_PDB)
  {
    ExtractPathAndName(PDB, PDBPATH, PDBNAME);
    if (!strcmp(PDBPATH, "") || !strcmp(PDBNAME, ""))
    {
      sprintf(errMsg, "in file %s line %d, no pdb was identified, program exits\n", __FILE__, __LINE__);
      TraceError(errMsg, ValueError);
      exit(ValueError);
    }
    GetPDBID(PDBNAME, PDBID);
    sprintf(FILE_SELF_ENERGY, "%s_selfenergy.txt", PDBID);
    sprintf(ROT_LIST_ALL, "%s_rotlist.txt", PDBID);
    sprintf(ROT_LIST_SEC, "%s_rotlistSEC.txt", PDBID);
    sprintf(ROT_INDEX_FILE, "%s_desrots", PDBID);
    sprintf(SEQ_FILE, "%s_desseqs", PDBID);
    sprintf(BEST_SEQ, "%s_bestseqs", PDBID);
    sprintf(BEST_STRUCT, "%s_beststruct", PDBID);
    sprintf(BEST_DESSITE, "%s_bestsites", PDBID);
    sprintf(BEST_MOL2, "%s_bestlig", PDBID);
  }

  if (strcmp(cmdname, "CheckClash1") == 0) sprintf(FILE_ATOMPAR, "%s/library/param_checkclash1.prm", PROGRAM_PATH);
  else if (strcmp(cmdname, "CheckClash2") == 0) sprintf(FILE_ATOMPAR, "%s/library/param_checkclash2.prm", PROGRAM_PATH);

  // read atom parameters
  AtomParamsSet atomParam;
  AtomParamsSetCreate(&atomParam);
  if (FAILED(AtomParameterRead(&atomParam, FILE_ATOMPAR)))
  {
    sprintf(errMsg, "in file %s line %d, failed to parse atom parameter file %s", __FILE__, __LINE__, FILE_ATOMPAR);
    TraceError(errMsg, IOError);
    exit(IOError);
  }

  // read residue topology
  ResiTopoSet resiTopo;
  ResiTopoSetCreate(&resiTopo);
  if (FAILED(ResiTopoSetRead(&resiTopo, FILE_TOPO)))
  {
    sprintf(errMsg, "in file %s line %d, failed to parse topology file %s", __FILE__, __LINE__, FILE_TOPO);
    TraceError(errMsg, IOError);
    exit(IOError);
  }

  Structure structure;
  StructureCreate(&structure);
  if (FLAG_PDB)
  {
    if (FAILED(StructureReadPDB(&structure, PDB, &atomParam, &resiTopo)))
    {
      sprintf(errMsg, "in file %s line %d, failed to parse pdb file %s", __FILE__, __LINE__, PDB);
      TraceError(errMsg, IOError);
      exit(IOError);
    }
    else
    {
      printf("read pdb file %s\n", PDB);
      StructureCalcPhiPsi(&structure);
    }
  }

  // read mol2 and generate parameter and topology files
  if (FLAG_MOL2)
  {
    if (!strcmp(cmdname, "MakeLigParamAndTopo"))
    {
      if (FAILED(GenerateSmallMolParameterAndTopologyFromMol2(MOL2, LIG_ATOMPAR, LIG_TOPFILE, INI_ATOM1, INI_ATOM2, INI_ATOM3)))
      {
        sprintf(errMsg, "in file %s line %d, failed to generate ligand atom parameter or ligand topology from mol2 %s, program exits", __FILE__, __LINE__, MOL2);
        TraceError(errMsg, ValueError);
        exit(ValueError);
      }
      else
      {
        printf("ligand atom parameter file '%s' and topology file '%s' has been generated\n", LIG_ATOMPAR, LIG_TOPFILE);
        exit(Success);
      }
    }

    if (FAILED(AtomParameterRead(&atomParam, LIG_ATOMPAR)))
    {
      sprintf(errMsg, "in file %s line %d, failed to read ligand atom parameters from %s, program exits. please check!", __FILE__, __LINE__, LIG_ATOMPAR);
      TraceError(errMsg, IOError);
      exit(IOError);
    }
    if (FAILED(ResiTopoSetRead(&resiTopo, LIG_TOPFILE)))
    {
      sprintf(errMsg, "in file %s line %d, failed to read ligand topology from %s, program exits. please check!", __FILE__, __LINE__, LIG_TOPFILE);
      TraceError(errMsg, IOError);
      exit(IOError);
    }

    if (FAILED(StructureReadMol2(&structure, MOL2, &atomParam, &resiTopo)))
    {
      sprintf(errMsg, "in file %s line %d, failed to parse mol2 file %s", __FILE__, __LINE__, MOL2);
      TraceError(errMsg, IOError);
      exit(IOError);
    }
    else
    {
      printf("read mol2 file %s\n", MOL2);
    }
  }

  int resiCount = 0;
  for (int i = 0;i < StructureGetChainCount(&structure);i++)
  {
    if (ChainGetType(StructureGetChain(&structure, i)) == Type_Chain_Protein)
    {
      resiCount += ChainGetResidueCount(StructureGetChain(&structure, i));
    }
  }
  PROT_LEN_NORM = resiCount > 0 ? resiCount : 1;

  // set the default rotamer library
  if (FLAG_USER_ROTLIB)
  {
    if (!(strcmp(USER_ROTLIB_NAME, "honig984") == 0 || strcmp(USER_ROTLIB_NAME, "honig3222") == 0 || strcmp(USER_ROTLIB_NAME, "honig7421") == 0 || strcmp(USER_ROTLIB_NAME, "honig11810") == 0 ||
      strcmp(USER_ROTLIB_NAME, "dun2010bb") == 0 || strcmp(USER_ROTLIB_NAME, "dun2010bb1per") == 0 || strcmp(USER_ROTLIB_NAME, "dun2010bb3per") == 0))
    {
      printf("please choose one of the following rotamer libraries (dun2010bb3per recommended):\n");
      printf("  dun2010bb\n  dun2010bb1per\n  dun2010bb3per\n");
      printf("  honig984\n  honig3222\n  honig7421\n  honig11810\n");
      exit(ValueError);
    }
    FLAG_BBDEP_ROTLIB = strstr(USER_ROTLIB_NAME, "honig") != NULL ? FALSE : TRUE;
    sprintf(FILE_ROTLIB, "%s/library/rotlib/%s.lib", PROGRAM_PATH, USER_ROTLIB_NAME);
  }
  else
  {
    if (FLAG_BBDEP_ROTLIB)
    {
      sprintf(FILE_ROTLIB, "%s/library/rotlib/dun2010bb3per.lib", PROGRAM_PATH);
      sprintf(FILE_ROTLIB_BIN, "%s/library/rotlib/ALLbbdep.bin", PROGRAM_PATH);
    }
    else sprintf(FILE_ROTLIB, "%s/library/rotlib/honig984.lib", PROGRAM_PATH);
  }

  // read energy weights from file
  if (FAILED(EnergyWeightRead(FILE_WEIGHT_READ)))
  {
    printf("warning! failed to parse energy weight file. all energy weights are set to 1\n");
  }
  else
  {
    printf("read energy weight file %s\n", FILE_WEIGHT_READ);
  }

  // deploy evolution information
  if (FLAG_EVOLUTION && strlen(DES_CHAINS) <= 1)
  {
    StructureDeployEvolutionInfo(&structure);
  }

  ////////////////////////////////////////////
  //            MAIN FUNCTION
  ////////////////////////////////////////////
#define MAIN_FUNC_PROTEIN_DESIGN
  if (!strcmp(cmdname, "ProteinDesign"))
  {
    // set protein chains for design
    if (!strcmp(DES_CHAINS, ""))
    {
      int proteinChainCount = 0;
      for (int i = 0;i < StructureGetChainCount(&structure);i++)
      {
        Chain* pChain = StructureGetChain(&structure, i);
        if (ChainGetType(pChain) == Type_Chain_Protein) proteinChainCount++;
      }
      if (proteinChainCount >= 1)
      {
        strcpy(DES_CHAINS, ChainGetName(StructureGetChain(&structure, 0)));
      }
      else
      {
        sprintf(errMsg, "in file %s line %d, no protein chain is identified and protein design fails", __FILE__, __LINE__);
        TraceError(errMsg, ValueError);
        exit(ValueError);
      }
    }

    AAppTable aapptable;
    RamaTable ramatable;
    BBdepRotamerLib bbrotlib;
    AApropensityTableReadFromFile(&aapptable, FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable, FILE_RAMACHANDRAN);
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    StructureCalcAminoAcidPropensityAndRamaEnergy(&structure, &aapptable, &ramatable);
    StructureCalcAminoAcidDunbrackEnergy(&structure, &bbrotlib);
    if (FLAG_PPI && FLAG_INTERFACE_ONLY)
    {
      if (FLAG_RESFILE) StructureBuildResfileRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
      else StructureBuildPPIRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo);
    }
    else if (FLAG_PROT_LIG)
    {
      StructureBuildPLIShell1RotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
      StructureBuildPLIShell2RotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
    }
    else if (FLAG_ENZYME)
    {
      if (FLAG_RESFILE)
      {
        StructureBuildCatalyticRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
        StructureBuildPLIShell1RotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
        StructureBuildPLIShell2RotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
      }
      else
      {
        sprintf(errMsg, "in file %s line %d, catalytic sites must be specified in RESFILE, otherwise catalytic sites may be mutated by the program\n"
          "this situation should be avoided in a meaningful enzyme design process", __FILE__, __LINE__);
        TraceError(errMsg, DataNotExistError);
        return DataNotExistError;
      }
    }
    else
    {
      if (FLAG_RESFILE)
      {
        StructureBuildResfileRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
      }
      else StructureBuildAllRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo);
    }
    BBdepRotamerLibDestroy(&bbrotlib);

    //deal with ligand rots if applicable
    if (FLAG_PROT_LIG || FLAG_ENZYME)
    {
      if (FLAG_LIG_CONFORMERS && !FAILED(StructureReadSmallMolRotamers(&structure, &resiTopo, FILE_LIG_CONFORMERS_IN)))
      {
        printf("read ligand conformers from %s\n", FILE_LIG_CONFORMERS_IN);
      }
      else
      {
        printf("use the ligand conformation in mol2 file for design\n");
        FILE* pOut = fopen(FILE_LIG_CONFORMERS_OUT, "w");
        Model(1, pOut);
        Residue* pSmallMol = NULL;
        StructureFindSmallMol(&structure, &pSmallMol);
        AtomArrayShowInPDBFormat(ResidueGetAllAtoms(pSmallMol), "ATOM", ResidueGetName(pSmallMol), ResidueGetChainName(pSmallMol), 1, ResidueGetPosInChain(pSmallMol), pOut);
        EndModel(pOut);
        fclose(pOut);
        StructureReadSmallMolRotamers(&structure, &resiTopo, FILE_LIG_CONFORMERS_OUT);
      }
    }

    StructureShowDesignSites(&structure);
    SelfEnergyGenerate2(&structure, &aapptable, &ramatable, FILE_SELF_ENERGY);
    RotamerList rotList;
    RotamerListCreateFromStructure(&rotList, &structure);
    RotamerListWrite(&rotList, ROT_LIST_ALL);
    SelfEnergyReadAndCheck(&structure, &rotList, FILE_SELF_ENERGY);
    RotamerListWrite(&rotList, ROT_LIST_SEC);
    RotamerListRead(&rotList, ROT_LIST_SEC);
    StructureShowDesignSitesAfterRotamerDelete(&structure, &rotList);
    SimulatedAnnealing(&structure, &rotList);
    RotamerListDestroy(&rotList);
  }
  else if (!strcmp(cmdname, "ComputeStability"))
  {
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable, FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable, FILE_RAMACHANDRAN);
    if (FLAG_BBDEP_ROTLIB)
    {
      ComputeStructureStabilityByBBdepRotLib2(&structure, &aapptable, &ramatable, FILE_ROTLIB_BIN, energyTerms);
    }
    else
    {
      ComputeStructureStability(&structure, &aapptable, &ramatable, energyTerms);
    }
  }
  else if (!strcmp(cmdname, "ComputeBinding"))
  {
    if (StructureGetChainCount(&structure) <= 2 || !FLAG_CHAIN_SPLIT)
    {
      ComputeBinding(&structure);
    }
    else ComputeBindingWithChainSplitting(&structure, SPLIT_PART1, SPLIT_PART2);
  }
  else if (!strcmp(cmdname, "RepairStructure"))
  {
    if (FLAG_BBDEP_ROTLIB)
    {
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
      StructureCalcAminoAcidDunbrackEnergy(&structure, &bbrotlib);
      RepairStructureByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, PDBID);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else
    {
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib, FILE_ROTLIB);
      RepairStructure(&structure, &rotlib, &atomParam, &resiTopo, PDBID);
      BBindRotamerLibDestroy(&rotlib);
    }
  }
  else if (!strcmp(cmdname, "Minimize"))
  {
    if (FLAG_BBDEP_ROTLIB)
    {
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
      StructureCalcAminoAcidDunbrackEnergy(&structure, &bbrotlib);
      EnergyMinimizationByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, PDBID);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else
    {
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib, FILE_ROTLIB);
      RepairStructure(&structure, &rotlib, &atomParam, &resiTopo, PDBID);
      BBindRotamerLibDestroy(&rotlib);
    }
  }
  else if (!strcmp(cmdname, "BuildMutant"))
  {
    if (FLAG_BBDEP_ROTLIB)
    {
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
      StructureCalcAminoAcidDunbrackEnergy(&structure, &bbrotlib);
      BuildMutantByBBdepRotLib(&structure, MUTANT_FILE, &bbrotlib, &atomParam, &resiTopo, PDBID);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else
    {
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib, FILE_ROTLIB);
      BuildMutant(&structure, MUTANT_FILE, &rotlib, &atomParam, &resiTopo, PDBID);
      BBindRotamerLibDestroy(&rotlib);
    }
  }
  else if (!strcmp(cmdname, "ComputeResiEnergy"))
  {
    if (strcmp(RESI, "") == 0)
    {
      printf("residue was not specified. Please use '--resi=ABB' to specify the residue, where X is the chain ID, YY is the residue position in the chain\n");
      exit(Warning);
    }
    char chnName[2] = "";
    chnName[0] = RESI[0];
    int pos = atoi(RESI + 1);
    int chnNdx = -100;
    StructureFindChainIndex(&structure, chnName, &chnNdx);
    Chain* pChain = NULL;
    if (chnNdx != -100)
    {
      pChain = StructureGetChain(&structure, chnNdx);
    }
    int resNdx = -1;
    ChainFindResidueByPosInChain(pChain, pos, &resNdx);
    Residue* pResi = NULL;
    if (resNdx != -1)
    {
      pResi = ChainGetResidue(pChain, resNdx);
    }
    BBdepRotamerLib bbrotlib;
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    StructureCalcAminoAcidDunbrackEnergy(&structure, &bbrotlib);
    printf("residue %s%d%s energy details:\n", ChainGetName(pChain), ResidueGetPosInChain(pResi), ResidueGetName(pResi));
    ComputeResidueInteractionWithFixedEnvironment(&structure, chnNdx, resNdx);
    BBdepRotamerLibDestroy(&bbrotlib);
  }

  else if (!strcmp(cmdname, "ComputeResiPairEnergy"))
  {
    if (strcmp(RESI_PAIR, "") == 0)
    {
      printf("residue was not specified. Please use '--resi_pair=ABB,CDD' to specify the residue, where A/C is the chain ID, BB/DD is the residue position in the chain\n");
      exit(Warning);
    }
    char part1[MAX_LEN_CHAIN_NAME + 1] = "";
    char part2[MAX_LEN_CHAIN_NAME + 1] = "";
    sscanf(RESI_PAIR, "%[^,],%s", part1, part2);
    char chnName1[2] = "";
    char chnName2[2] = "";
    chnName1[0] = part1[0];
    chnName2[0] = part2[0];
    int pos1 = atoi(part1 + 1);
    int pos2 = atoi(part2 + 1);
    int chnNdx1 = -100, chnNdx2 = -100;
    StructureFindChainIndex(&structure, chnName1, &chnNdx1);
    StructureFindChainIndex(&structure, chnName2, &chnNdx2);
    Chain* pChain1 = NULL;
    Chain* pChain2 = NULL;
    if (chnNdx1 != -100 && chnNdx2 != -100)
    {
      pChain1 = StructureGetChain(&structure, chnNdx1);
      pChain2 = StructureGetChain(&structure, chnNdx2);
    }
    int resNdx1 = -1, resNdx2 = -1;
    ChainFindResidueByPosInChain(pChain1, pos1, &resNdx1);
    ChainFindResidueByPosInChain(pChain2, pos2, &resNdx2);
    Residue* pResi1 = NULL;
    Residue* pResi2 = NULL;
    if (resNdx1 != -1 && resNdx2 != -1)
    {
      pResi1 = ChainGetResidue(pChain1, resNdx1);
      pResi2 = ChainGetResidue(pChain2, resNdx2);
    }
    double energyTerms[MAX_ENERGY_TERM] = { 0 };
    if (strcmp(chnName1, chnName2) == 0)
    { // same chain
      if (pos1 == pos2 + 1 || pos2 == pos1 + 1)
      {
        EnergyResidueAndNextResidue(pResi1, pResi2, energyTerms);
      }
      else
      {
        EnergyResidueAndOtherResidueSameChain(pResi1, pResi2, energyTerms);
      }
    }
    else
    { // different chains
      if (ChainGetType(pChain1) == Type_Chain_SmallMol)
      {
        EnergyResidueAndLigandResidue(pResi2, pResi1, energyTerms);
      }
      else if (ChainGetType(pChain2) == Type_Chain_SmallMol)
      {
        EnergyResidueAndLigandResidue(pResi1, pResi2, energyTerms);
      }
      else
      {
        EnergyResidueAndOtherResidueDiffChain(pResi1, pResi2, energyTerms);
      }
    }
    EnergyTermWeighting(energyTerms);
    printf("interaction energy between residue %s%d%s and residue %s%d%s:\n",
      ChainGetName(pChain1), ResidueGetPosInChain(pResi1), ResidueGetName(pResi1),
      ChainGetName(pChain2), ResidueGetPosInChain(pResi2), ResidueGetName(pResi2));
    EnergyTermShowComplex(energyTerms);
  }

  else if (!strcmp(cmdname, "ComputeRotamersEnergy"))
  {
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable, FILE_AAPROPENSITY);
    RamaTableReadFromFile(&ramatable, FILE_RAMACHANDRAN);
    if (FLAG_BBDEP_ROTLIB)
    {
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
      if (FLAG_WILDTYPE_ONLY) ComputeWildtypeRotamersEnergyByBBdepRotLib(&structure, &bbrotlib, &aapptable, &ramatable, &atomParam, &resiTopo, PDBID);
      else ComputeRotamersEnergyByBBdepRotLib(&structure, &bbrotlib, &aapptable, &ramatable, &atomParam, &resiTopo, PDBID);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else
    {
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib, FILE_ROTLIB);
      if (FLAG_WILDTYPE_ONLY) ComputeWildtypeRotamersEnergy(&structure, &rotlib, &aapptable, &ramatable, &atomParam, &resiTopo, PDBID);
      else ComputeAllRotamersEnergy(&structure, &rotlib, &aapptable, &ramatable, &atomParam, &resiTopo, PDBID);
      BBindRotamerLibDestroy(&rotlib);
    }
  }
  
  else if (!strcmp(cmdname, "ShowResiComposition"))
  {
    int aas[30] = { 0 };
    StructureGetAminoAcidComposition(&structure, aas);
  }
  
  else if (!strcmp(cmdname, "OptimizeHydrogen"))
  {
    OptimizeHydrogen(&structure, &atomParam, &resiTopo, PDBID);
  }
  
  else if (!strcmp(cmdname, "AddPolarHydrogen"))
  {
    AddPolarHydrogen(&structure, PDBID);
  }
  
  else if (!strcmp(cmdname, "FindInterfaceResidue"))
  {
    if (StructureGetChainCount(&structure) < 2)
    {
      printf("there is only one chain in the whole structure, no protein-protein interface found\n");
      exit(ValueError);
    }
    else if (StructureGetChainCount(&structure) == 2)
    {
      FindInterfaceResidues(&structure);
    }
    else
    {
      printf("there are >= 3 chains in the whole structure, please split chains before calculating interface residues\n");
      FindInterfaceResiduesWithChainSplitting(&structure, SPLIT_PART1, SPLIT_PART2);
    }
  }
  
  else if (!strcmp(cmdname, "SelectResiduesWithin"))
  {
    SelectResiduesWithin(&structure);
  }
  
  else if (!strcmp(cmdname, "MakeLigConformers"))
  {
    if (!strcmp(FILE_LIG_CONFORMERS_OUT, ""))
    {
      sprintf(errMsg, "in file %s line %d, no file designated to output ligand conformers, e.g. please use:\n"
        "--write_lig_conformers=LIG_CONFORMERS.pdb to write the ligand rotamers into a file named 'LIG_CONFORMERS.pdb'", __FILE__, __LINE__);
      TraceError(errMsg, InvalidInputError);
      exit(InvalidInputError);
    }

    BBdepRotamerLib bbrotlib;
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    StructureBuildResfileRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
    BBdepRotamerLibDestroy(&bbrotlib);
    StructureShowDesignSites(&structure);
    StructureGenerateSmallMolRotamers(&structure, FILE_CATACONS, FILE_LIG_PLACEMENT);
    StructureWriteSmallMolRotamers(&structure, FILE_LIG_CONFORMERS_OUT);
    StructureShowDesignSites(&structure);
  }
  
  else if (!strcmp(cmdname, "AnalyzeLigConformers"))
  {
    if (!strcmp(FILE_LIG_CONFORMERS_IN, ""))
    {
      sprintf(errMsg, "in file %s line %d, ligand conformers to be analyzed was not specified, e.g. please use:\n"
        "--read_lig_conformers=LIG_CONFORMERS.pdb to specify the file for saving ligand conformers", __FILE__, __LINE__);
      TraceError(errMsg, InvalidInputError);
      exit(InvalidInputError);
    }
    Residue* pSmallmol = NULL;
    StructureFindSmallMol(&structure, &pSmallmol);
    AnalyzeSmallMolRotamers(FILE_LIG_CONFORMERS_IN, pSmallmol);
  }
  
  else if (!strcmp(cmdname, "ScreenLigConformers"))
  {
    if (FLAG_LIG_SCREEN_BY_ORIENTATION)
    {
      StructureSmallmolOrientationScreen(&structure, &resiTopo, FILE_LIG_CONFORMERS_IN, FILE_LIG_CONFORMERS_OUT, LIGSCREEN_ORITENTATION_FILE);
    }
    else if (FLAG_LIG_SCREEN_BY_TOPVDW)
    {
      SmallMolRotamersGetBothHighRankOfBackboneVdwAndInternalVdw(FILE_LIG_CONFORMERS_IN, FILE_LIG_CONFORMERS_OUT, LIGSCREEN_TOPVDW_PERCENT);
    }
    else if (FLAG_LIG_SCREEN_BY_RMSD)
    {
      ScreenSmallmolRotamersByRMSD(FILE_LIG_CONFORMERS_IN, FILE_LIG_CONFORMERS_OUT, LIGSCREEN_RMSD_CUTOFF);
    }
  }
  
  else if (!strcmp(cmdname, "ComputeEvolutionScore"))
  {
    double evoscore = EvolutionScorePrfFromFile(TGT_SEQ);
    printf("evoscore: %f\n", evoscore);
  }
  
  else if (!strcmp(cmdname, "OptimizeWeight"))
  {
    PNATAA_WeightOptByGradientDescent(PDBLIST);
  }
  
  else if (!strcmp(cmdname, "CalcPhiPsi"))
  {
    char PHI_PSI_FILE[MAX_LEN_ONE_LINE_CONTENT + 1] = "";
    sprintf(PHI_PSI_FILE, "%s_phipsi.txt", PDBID);
    ShowPhiPsi(&structure, PHI_PSI_FILE);
  }
  
  else if (!strcmp(cmdname, "CheckResiInRotLib"))
  {
    if (FLAG_BBDEP_ROTLIB)
    {
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
      CheckRotamerInBBdepRotLib(&structure, &bbrotlib, &resiTopo, CUT_TORSION_DEVIATION, PDBID);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else
    {
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib, FILE_ROTLIB);
      CheckRotamerInBBindRotLib(&structure, &rotlib, &resiTopo, CUT_TORSION_DEVIATION, PDBID);
      BBindRotamerLibDestroy(&rotlib);
    }
  }
  
  else if (!strcmp(cmdname, "CalcResiMinRMSD"))
  {
    if (FLAG_BBDEP_ROTLIB)
    {
      BBdepRotamerLib bbrotlib;
      BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
      StructureGenerateWildtypeRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo);
      BBdepRotamerLibDestroy(&bbrotlib);
    }
    else
    {
      BBindRotamerLib rotlib;
      BBindRotamerLibCreate(&rotlib, FILE_ROTLIB);
      StructureGenerateWildtypeRotamers(&structure, &rotlib, &atomParam, &resiTopo);
      BBindRotamerLibDestroy(&rotlib);
    }
    FindMinRmsdRotFromRotLib(&structure, PDBID);
  }
  
  else if (!strcmp(cmdname, "CompareSideChain"))
  {
    Structure newStruct;
    StructureCreate(&newStruct);
    StructureReadPDB(&newStruct, PDB2, &atomParam, &resiTopo);
    StructureCalcPhiPsi(&newStruct);
    StructureCalcProteinResidueSidechainTorsion(&newStruct, &resiTopo);
    char RMSDFILE[MAX_LEN_FILE_NAME + 1];
    char TORSIONFILE[MAX_LEN_FILE_NAME + 1];
    sprintf(RMSDFILE, "%s_cmprmsd.txt", PDBID);
    sprintf(TORSIONFILE, "%s_cmptorsion.txt", PDBID);
    FILE* pTorsion = fopen(TORSIONFILE, "w");
    FILE* pRMSD = fopen(RMSDFILE, "w");
    CompareSidechainsOf2Structures(&structure, &newStruct, pTorsion, pRMSD);
    fclose(pTorsion);
    fclose(pRMSD);
    StructureDestroy(&newStruct);
  }
  
  else if (!strcmp(cmdname, "CheckClash0"))
  {
    CheckClash0(&structure, CUT_CLASH_RATIO);
  }
  
  else if (!strcmp(cmdname, "CheckClash1"))
  {
    CheckClash1(&structure, CUT_CLASH_RATIO);
  }
  
  else if (!strcmp(cmdname, "CheckClash2"))
  {
    CheckClash2(&structure, CUT_CLASH_RATIO);
  }
  
  else if (!strcmp(cmdname, "PredSecStr"))
  {
    SSPred(TGT_SEQ);
  }
  
  else if (!strcmp(cmdname, "PredSolAcc"))
  {
    SAPred(TGT_SEQ);
  }
  
  else if (!strcmp(cmdname, "PredPhiPsi"))
  {
    PhiPsiPred(TGT_SEQ);
  }
  
  else if (!strcmp(cmdname, "CalcResiBfactor"))
  {
    char FILE_BFACTOR[MAX_LEN_FILE_NAME + 1];
    sprintf(FILE_BFACTOR, "%s_bfactor.txt", PDBID);
    FILE* pFile = fopen(FILE_BFACTOR, "w");
    for (int i = 0;i < StructureGetChainCount(&structure);i++)
    {
      Chain* pChain = StructureGetChain(&structure, i);
      if (ChainGetType(pChain) != Type_Chain_Protein) continue;
      for (int j = 0;j < ChainGetResidueCount(pChain);j++)
      {
        Residue* pResi = ChainGetResidue(pChain, j);
        double bfac = ResidueGetAverageBfactor(pResi);
        fprintf(pFile, "%3s %8.3f\n", ResidueGetName(pResi), bfac);
      }
    }
    fclose(pFile);
  }
  
  else if (!strcmp(cmdname, "FindCoreResidue"))
  {
    FindCoreResidues(&structure);
  }
  
  else if (!strcmp(cmdname, "FindSurfaceResidue"))
  {
    FindSurfaceResidues(&structure);
  }
  
  else if (!strcmp(cmdname, "FindIntermediateResidue"))
  {
    FindIntermediateResidues(&structure);
  }
  
  else if (!strcmp(cmdname, "WriteFirstLigRotIntoMol2"))
  {
    StructureReadSmallMolRotamers(&structure, &resiTopo, FILE_LIG_CONFORMERS_IN);
    StructureWriteFirstLigandRotamerIntoMol2(&structure, "TheFirstLigRotamer.mol2");
  }
  
  else
  {
    printf("unknown command name: %s, program exits\n", cmdname);
    exit(ValueError);
  }

  if (FLAG_EVOLUTION)
  {
    Chain* pChain = StructureFindChainByName(&structure, DES_CHAINS);
    int len = ChainGetResidueCount(pChain);
    if (PROT_PROFILE)
    {
      for (int i = 0;i < len;i++) delete[] PROT_PROFILE[i];
      delete[] PROT_PROFILE;
    }
  }
  ResiTopoSetDestroy(&resiTopo);
  AtomParamsSetDestroy(&atomParam);
  StructureDestroy(&structure);

  clock_t timeend = clock();
  SpentTimeShow(timestart, timeend);

  return Success;
}



