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
#pragma warning(disable : 4996)
#pragma warning(disable : 6031)
#pragma warning(disable : 26812)
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include "Getopt.h"
#include "ProgramFunction.h"
#include "RotamerBuilder.h"
#include "Structure.h"
#include "Evolution.h"
#include "ProteinDesign.h"
#include "WeightOpt.h"
#include "SmallMolParAndTopo.h"
#include "ProgramPreprocess.h"

// Paths, Parameters, and Paths that should NOT be changed by low-level users
char PROGRAM_PATH[MAX_LEN_FILE_NAME + 1] = ".";
char PROGRAM_NAME[MAX_LEN_FILE_NAME + 1] = "UniDesign";
char PROGRAM_VERSION[MAX_LEN_FILE_NAME + 1] = "v1.3";
char FILE_ATOMPARAM[MAX_LEN_FILE_NAME + 1] = "library/toppar/param_charmm19_lk.prm";
char FILE_TOPO[MAX_LEN_FILE_NAME + 1] = "library/toppar/top_polh19.inp";
char FILE_ENERGY_AAPP[MAX_LEN_FILE_NAME + 1] = "library/eterms/aapropensity.nrg";
char FILE_ENERGY_RAMA[MAX_LEN_FILE_NAME + 1] = "library/eterms/ramachandran.nrg";
char FILE_WEIGHT_READ[MAX_LEN_FILE_NAME + 1] = "wread/weight_all1.wgt";
char FILE_ROTLIB[MAX_LEN_FILE_NAME + 1] = "library/rotlib/dun2010bb3per.lib";
char FILE_ROTLIB_BIN[MAX_LEN_FILE_NAME + 1] = "library/rotlib/ALLbbdep.bin";
double ROTLIB_LOW_PROB_CUTOFF = 0.03; // 0.01-0.03 suggested
char USER_ROTLIB_NAME[MAX_LEN_FILE_NAME + 1] = "";
char TGT_PRF[MAX_LEN_FILE_NAME + 1] = "prf.txt";
char TGT_MSA[MAX_LEN_FILE_NAME + 1] = "msa.txt";
char TGT_SA[MAX_LEN_FILE_NAME + 1] = "sa.txt";
char TGT_SS[MAX_LEN_FILE_NAME + 1] = "ss.txt";
char TGT_SEQ[MAX_LEN_FILE_NAME + 1] = "seq.txt";
char TGT_PHIPSI[MAX_LEN_FILE_NAME + 1] = "phi-psi.txt";
float **PROT_PROFILE = NULL;
char FILE_SELF_ENERGY[MAX_LEN_FILE_NAME + 1] = "selfenergy.txt";
char FILE_ROTLIST[MAX_LEN_FILE_NAME + 1] = "rotlist.txt";
char FILE_ROTLIST_SEC[MAX_LEN_FILE_NAME + 1] = "rotlistSEC.txt";
char TAG_BEST_DESIGNS[MAX_LEN_FILE_NAME + 1] = "designs.tsv";
char TAG_BEST_SEQS[MAX_LEN_FILE_NAME + 1] = "seqs.tsv";
char TAG_BEST_MODEL[MAX_LEN_FILE_NAME + 1] = "mod";
char TAG_BEST_LIG_MOL2[MAX_LEN_FILE_NAME + 1] = "lig";
double TORSION_DEVIATION_THRESHOLD = 20.0;
double CLASH_RATIO_THRESHOLD = 0.75;    // clash is counted if dist(i, j) < clash_ratio * (vdw_r_i + vdw_r_j)
double CLASH_MIN_DISTANCE = 2.0;       // (Angstrom): clash if dist(i, j) < 2.0
double CLASH_DISTANCE_THRESHOLD = 0.4; // (Angstrom): clash if dist(i, j) < vdw_r_i + vdw_r_j - 0.4
int CUT_NUM_CB_CORE = 20;
int CUT_NUM_CB_SURF = 15;
int PROT_LEN_NORM = 100;

double SA_T_INITIAL = 50;
double SA_T_FINAL = 0.001;
int SA_CYCLE = 2;
double SA_DECREASE_FACTOR = 0.8;
int NUM_METROPOLIS_STEPS = 20000;

// Paths, Parameters, and Paths that can be modified by users
char PDB[MAX_LEN_FILE_NAME + 1] = "PROTEIN.pdb";
char PDB2[MAX_LEN_FILE_NAME + 1] = "PROTEIN2.pdb";
char PDBPATH[MAX_LEN_FILE_NAME + 1] = ".";
char PDBNAME[MAX_LEN_FILE_NAME + 1] = "PROTEIN.pdb";
char PDBID[MAX_LEN_FILE_NAME + 1] = "PROTEIN";
char MOL2[MAX_LEN_FILE_NAME + 1] = "LIG.mol2";

char DES_CHAINS[10] = "A";
int NTRAJ = 1;
int NTRAJ_START_NDX = 1;

double WGT_PROFILE = 1.0;
double WGT_BIND = 1.0;
double WGT_CATA_CONS = 5.0;
int MIN_MUT_NUM = 0;
int MAX_MUT_NUM = 1000;
double WGT_MUT_NUM = 1000.0;  // strong penalty to disfavor out-of-range number of mutations
double WGT_SAME_SEQ = 1000.0; // strong penalty to disfavor identical designs

double CUT_PPI_DIST_SHELL1 = 5.0;
double CUT_PPI_DIST_SHELL2 = 8.0;
double CUT_PLI_DIST_SHELL1 = 4.5;
double CUT_PLI_DIST_SHELL2 = 7.5;

char INI_ATOM1[MAX_LEN_ATOM_NAME + 1] = "C1";
char INI_ATOM2[MAX_LEN_ATOM_NAME + 1] = "C2";
char INI_ATOM3[MAX_LEN_ATOM_NAME + 1] = "C3";

char SPLIT_CHAINS[10] = "AB,C";
char SPLIT_PART1[10] = "AB";
char SPLIT_PART2[10] = "C";

char RESIDUE[10] = "A1"; // the first character indicates the chain ID, and the digits indicate the residue position on the chain
char EXCLUDE_RESIDUES[MAX_LEN_FILE_NAME + 1] = "";
char RESIDUE_PAIR[20] = "";

char MUTANT_FILE[MAX_LEN_FILE_NAME + 1] = "MUTANT_FILE.txt";
int MAX_NUM_OF_RUNS = 4;

// parameters for MakeLigPoses, ScreenLigPoses, and enzyme design
char FILE_LIG_POSES_READ[MAX_LEN_FILE_NAME + 1] = "LIG_POSES.pdb";
char FILE_LIG_POSES_WRITE[MAX_LEN_FILE_NAME + 1] = "LIG_POSES.pdb";
char FILE_PICK_POSE_BY_ORITENTATION[MAX_LEN_FILE_NAME + 1] = "PICK_LIG_POSES_BY_ORIENTATION_RULE.txt";
double PICK_POSE_BBVDW_RANK_CUTOFF = 0.5; // default: 0.5 (50%)
double PICK_POSE_INVDW_RANK_CUTOFF = 0.5; // default: 0.5 (50%)
double PICK_POSE_RMSD_CUTOFF = 0.5;       // default: 0.5 Angstrom

char FILE_PDBLIST[MAX_LEN_FILE_NAME + 1] = "PDB_LIST.txt";

char FILE_LIG_PARAM[MAX_LEN_FILE_NAME + 1] = "LIG_PARAM.prm";
char FILE_LIG_TOPO[MAX_LEN_FILE_NAME + 1] = "LIG_TOPO.inp";
char FILE_LIG_CATACONS[MAX_LEN_FILE_NAME + 1] = "LIG_CATACONS.txt";
char FILE_LIG_PLACING[MAX_LEN_FILE_NAME + 1] = "LIG_PLACING.txt";
char FILE_RESFILE[MAX_LEN_FILE_NAME + 1] = "RESFILE.txt";
char OUTPUT_PREFIX[MAX_LEN_FILE_NAME + 1] = "UD"; // UD stands for UniDesign

char REFERENCE_RESIDUES[MAX_LEN_FILE_NAME + 1];
double DISTANCE_TO_REFERENCE_RESIDUES = 9.0;

// Flags to control UniDesign workflow
BOOL FLAG_PDB = FALSE;                          // default: FALSE
BOOL FLAG_MOL2 = FALSE;                         // default: FALSE
BOOL FLAG_CHAIN_SPLIT = FALSE;                  // default: FALSE
BOOL FLAG_USE_LIG_POSES = FALSE;                // default: FALSE
BOOL FLAG_BY_ORIENTATION = FALSE;               // default: FALSE
BOOL FLAG_BY_VDW_PERCENTILE = FALSE;            // default: FALSE
BOOL FLAG_BY_RMSD = FALSE;                      // default: FALSE
BOOL FLAG_DESIGN_MON = FALSE;                   // default: FALSE
BOOL FLAG_DESIGN_PPI = FALSE;                   // default: FALSE
BOOL FLAG_DESIGN_PLI = FALSE;                   // default: FALSE
BOOL FLAG_DESIGN_ENZ = FALSE;                   // default: FALSE
BOOL FLAG_EVOLUTION = FALSE;                    // default: FALSE
BOOL FLAG_EVOPHIPSI = FALSE;                    // default: FALSE
BOOL FLAG_USE_INPUT_SC = TRUE;                  // default: TRUE
BOOL FLAG_USE_RESFILE = FALSE;                  // default: FALSE
BOOL FLAG_ROTATE_HYDROXYL = TRUE;               // default: TRUE
BOOL FLAG_WILDTYPE_ONLY = FALSE;                // default: FALSE
BOOL FLAG_INTERFACE_ONLY = FALSE;               // default: FALSE
BOOL FLAG_EXCL_CYS_ROTS = FALSE;                // default: FALSE
BOOL FLAG_DESIGN_FROM_NATAA = FALSE;            // default: FALSE
BOOL FLAG_READ_HYDROGEN = TRUE;                 // default: TRUE
BOOL FLAG_PENALIZE_IDENTICAL_SEQUENCES = FALSE; // default: FALSE
BOOL FLAG_RESTRICT_NUM_OF_MUTATIONS = FALSE;    // default: FALSE
BOOL FLAG_WRITE_MUTATIONS = FALSE;              // default: FALSE

const int N_CMD = 200;
char SUPPORTED_COMMANDS[][30] = {
    "ComputeBinding",
    "ComputeEvolutionScore",
    "ComputeResEnergy",
    "ComputeResPairEnergy",
    "ComputeRotamersEnergy",
    "ComputeStability",
    "ProteinDesign",
    "AddPolarHydrogen",
    "BuildMutant",
    "CalcPhiPsi",
    "CalcResMinRMSD",
    "CheckClash0",
    "CheckClash1",
    "CheckClash2",
    "CheckResInRotLib",
    "CompareProtSideChain",
    "FindCoreRes",
    "FindIntermediateRes",
    "FindInterfaceRes",
    "FindSurfaceRes",
    "Minimize",
    "OptimizeHydrogen",
    "RepairStructure",
    "ShowResComposition",
    "PredPhiPsi",
    "PredSS",
    "PredSA",
    "CalcResBfactor",
    "OptimizeWeight",
    "MakeLigParamAndTopo",
    "MakeLigPoses",
    "AnalyzeLigPoses",
    "ScreenLigPoses",
    "WriteFirstLigPoseToMol2",
    "SelectResWithin",
    "WriteFasta",
    "FindMutations",
    NULL,
};

const char *SHORT_OPTS = "?hv";
struct option LONG_OPTS[] = {
    {"help", no_argument, NULL, 1},
    {"version", no_argument, NULL, 2},
    {"command", required_argument, NULL, 3},
    {"pdb", required_argument, NULL, 4},
    {"split_chains", required_argument, NULL, 5},
    {"mutant_file", required_argument, NULL, 6},
    {"use_input_sc", required_argument, NULL, 8},
    {"ppi_shell1", required_argument, NULL, 9},
    {"ppi_shell2", required_argument, NULL, 10},
    {"pli_shell1", required_argument, NULL, 26},
    {"pli_shell2", required_argument, NULL, 27},
    {"evolution", no_argument, NULL, 11},
    {"design_mon", no_argument, NULL, 13},
    {"design_ppi", no_argument, NULL, 14},
    {"design_pli", no_argument, NULL, 15},
    {"design_enz", no_argument, NULL, 16},
    {"design_chains", required_argument, NULL, 17},
    {"wildtype_only", no_argument, NULL, 23},
    {"interface_only", no_argument, NULL, 32},
    {"resfile", required_argument, NULL, 36},
    {"seed_from_wildtype", no_argument, NULL, 38},
    {"exclude_cys_rots", no_argument, NULL, 40},
    {"mol2", required_argument, NULL, 18},
    {"rotate_hydroxyl", required_argument, NULL, 19},
    {"xdeviation", required_argument, NULL, 20},
    {"wread", required_argument, NULL, 21},
    {"pdblist", required_argument, NULL, 22},
    {"compare_pdb", required_argument, NULL, 25},
    {"clash_ratio", required_argument, NULL, 28},
    {"ntraj", required_argument, NULL, 29},
    {"low_prob_rot_cutoff", required_argument, NULL, 30},
    {"seq_file", required_argument, NULL, 33},
    {"all_evo_terms", no_argument, NULL, 34},
    {"wprof", required_argument, NULL, 35},

    {"prefix", required_argument, NULL, 37},

    {"ncut_cb_core", required_argument, NULL, 42},
    {"ncut_cb_surf", required_argument, NULL, 43},

    {"init_3atoms", required_argument, NULL, 44}, 
    {"read_lig_poses", required_argument, NULL, 45}, 
    {"write_lig_poses", required_argument, NULL, 46},  
    {"by_orientation", required_argument, NULL, 47}, 
    {"by_bbvdw_perc", required_argument, NULL, 48}, 
    {"by_invdw_perc", required_argument, NULL, 480}, 
    {"by_rmsd", required_argument, NULL, 49}, 
    {"lig_param", required_argument, NULL, 51},
    {"lig_topo", required_argument, NULL, 52},
    {"lig_catacons", required_argument, NULL, 53},
    {"num_of_runs", required_argument, NULL, 54},
    {"within_residues", required_argument, NULL, 55},
    {"within_range", required_argument, NULL, 56},
    {"ntraj_start_ndx", required_argument, NULL, 57},
    {"wbind", required_argument, NULL, 58},
    {"residue", required_argument, NULL, 59},
    {"residue_pair", required_argument, NULL, 60},
    {"exclude_residues", required_argument, NULL, 61},
    {"lig_placing", required_argument, NULL, 62},
    // 'val == 63' cannot be used, as the compilation always reports error C2196: case value '63' already used
    {"min_muts", required_argument, NULL, 64}, 
    {"max_muts", required_argument, NULL, 65}, 
    {"penalize_exist_seqs", no_argument, NULL, 66}, 
    {"min_temperature", required_argument, NULL, 67}, 
    {"max_temperature", required_argument, NULL, 68}, 
    {"num_of_samc_cycles", required_argument, NULL, 69}, 
    {"metropolis_steps", required_argument, NULL, 70}, 
    {"samc_temperature_decay_rate", required_argument, NULL, 71}, 
    {"write_mutations", no_argument, NULL, 72}, 
    {NULL, no_argument, NULL, 0},
};

static BOOL CheckCommandName(const char *cmdname)
{
  for (int i = 0; i < N_CMD; i++)
  {
    if (SUPPORTED_COMMANDS[i] == NULL)
      break;
    else if (strcmp(cmdname, SUPPORTED_COMMANDS[i]) == 0)
      return TRUE;
  }
  return FALSE;
}

static int ParseOptions(int argc, char **argv)
{
  return getopt_long(argc, argv, SHORT_OPTS, LONG_OPTS, NULL);
}

static int SetProteinDesignType(int code)
{
  FLAG_DESIGN_PPI = FALSE;
  FLAG_DESIGN_MON = FALSE;
  FLAG_DESIGN_PLI = FALSE;
  FLAG_DESIGN_ENZ = FALSE;
  switch (code)
  {
  case 14:
    FLAG_DESIGN_PPI = TRUE;
    break;
  case 15:
    FLAG_DESIGN_PLI = TRUE;
    break;
  case 16:
    FLAG_DESIGN_ENZ = TRUE;
    break;
  default:
    FLAG_DESIGN_MON = TRUE;
    break;
  }
  return Success;
}

int PrintCommandLine(int argc, char *argv[])
{
  printf("your input command line: ");
  for (int i = 0; i < argc; i++)
  {
    printf("%s ", argv[i]);
  }
  printf("\n");
  return Success;
}

int main(int argc, char *argv[])
{
  const char *cmdname = "";
  char errMsg[MAX_LEN_ERR_MSG + 1];
  clock_t timestart = clock();
  setvbuf(stdout, NULL, _IONBF, 0);
  PrintProgramInfo();
  PrintCommandLine(argc, argv);

  // Set default filepaths for parameters and weights
  ExtractPathAndName(argv[0], PROGRAM_PATH, PROGRAM_NAME);
  sprintf(FILE_ATOMPARAM, "%s/library/toppar/param_charmm19_lk.prm", PROGRAM_PATH);
  sprintf(FILE_TOPO, "%s/library/toppar/top_polh19.inp", PROGRAM_PATH);
  sprintf(FILE_ROTLIB, "%s/library/rotlib/dun2010bb3per.lib", PROGRAM_PATH);
  sprintf(FILE_ENERGY_AAPP, "%s/library/eterms/aapropensity.nrg", PROGRAM_PATH);
  sprintf(FILE_ENERGY_RAMA, "%s/library/eterms/ramachandran.nrg", PROGRAM_PATH);
  sprintf(FILE_WEIGHT_READ, "%s/wread/weight_all1.wgt", PROGRAM_PATH);

  while (TRUE)
  {
    int code = ParseOptions(argc, argv);
    if (code == -1)
    {
      break;
    }
    switch (code)
    {
    case '?':
      PrintProgramHelp();
      exit(Success);
    case 'h':
      PrintProgramHelp();
      exit(Success);
    case 'v':
      PrintProgramVersion();
      exit(Success);
    case 1:
      PrintProgramHelp();
      exit(Success);
    case 2:
      PrintProgramVersion();
      exit(Success);
    case 3:
      cmdname = optarg;
      if (CheckCommandName(cmdname) == FALSE)
      {
        printf("command '%s' is not supported, program exits\n", cmdname);
        exit(ValueError);
      }
      else
      {
        printf("command '%s' is supported\n", cmdname);
      }
      break;
    case 4: // pdb
      FLAG_PDB = TRUE;
      strcpy(PDB, optarg);
      break;
    case 5: // split_chains
      strcpy(SPLIT_CHAINS, optarg);
      sscanf(SPLIT_CHAINS, "%[^,],%s", SPLIT_PART1, SPLIT_PART2);
      FLAG_CHAIN_SPLIT = TRUE;
      for (int i = 0; i < (int)strlen(SPLIT_PART1); i++)
      {
        char tmp[2] = {SPLIT_PART1[i], '\0'};
        if (strstr(SPLIT_PART2, tmp) != NULL)
        {
          printf("for option 'split_chains', the two splitted parts should not have identical chains, program exits\n");
          exit(FormatError);
        }
      }
      break;
    case 6: // mutant_file
      strcpy(MUTANT_FILE, optarg);
      break;
    case 8: // use_input_sc
      if (strcmp(optarg, "true") == 0)
      {
        FLAG_USE_INPUT_SC = FALSE;
      }
      else
      {
        FLAG_USE_INPUT_SC = TRUE;
      }
      break;
    case 9: // ppi_shell1
      CUT_PPI_DIST_SHELL1 = atof(optarg);
      break;
    case 10: // ppi_shell2
      CUT_PPI_DIST_SHELL2 = atof(optarg);
      break;
    case 11:
      FLAG_EVOLUTION = TRUE;
      break;
    case 13:
    case 14:
    case 15:
    case 16:
      SetProteinDesignType(code);
      break;
    case 17:
      strcpy(DES_CHAINS, optarg);
      break;
    case 18:
      strcpy(MOL2, optarg);
      FLAG_MOL2 = TRUE;
      break;
    case 19:
      if (strcmp(optarg, "no") == 0)
      {
        FLAG_ROTATE_HYDROXYL = FALSE;
      }
      else
      {
        FLAG_ROTATE_HYDROXYL = TRUE;
      }
      break;
    case 20:
      TORSION_DEVIATION_THRESHOLD = atof(optarg);
      break;
    case 21:
      strcpy(FILE_WEIGHT_READ, optarg);
      break;
    case 22:
      strcpy(FILE_PDBLIST, optarg);
      break;
    case 23:
      FLAG_WILDTYPE_ONLY = TRUE;
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
      CLASH_RATIO_THRESHOLD = atof(optarg);
      break;
    case 29:
      NTRAJ = atoi(optarg);
      break;
    case 30:
      ROTLIB_LOW_PROB_CUTOFF = atof(optarg);
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
      FLAG_USE_RESFILE = TRUE;
      break;
    case 38:
      FLAG_DESIGN_FROM_NATAA = TRUE;
      break;
    case 40:
      FLAG_EXCL_CYS_ROTS = TRUE;
      break;
    case 42:
      CUT_NUM_CB_CORE = atoi(optarg);
      break;
    case 43:
      CUT_NUM_CB_SURF = atoi(optarg);
      break;
    case 44: // init_3atoms
      StringArray atoms;
      StringArrayCreate(&atoms);
      StringArraySplitString(&atoms, optarg, ',');
      if (StringArrayGetNum(&atoms) < 3)
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
      FLAG_USE_LIG_POSES = TRUE;
      strcpy(FILE_LIG_POSES_READ, optarg);
      break;
    case 46:
      strcpy(FILE_LIG_POSES_WRITE, optarg);
      break;
    case 47:
      FLAG_BY_ORIENTATION = TRUE;
      strcpy(FILE_PICK_POSE_BY_ORITENTATION, optarg);
      break;
    case 48:
      FLAG_BY_VDW_PERCENTILE = TRUE;
      PICK_POSE_BBVDW_RANK_CUTOFF = atof(optarg);
      break;
    case 480:
      FLAG_BY_VDW_PERCENTILE = TRUE;
      PICK_POSE_INVDW_RANK_CUTOFF = atof(optarg);
      break;
    case 49:
      FLAG_BY_RMSD = TRUE;
      PICK_POSE_RMSD_CUTOFF = atof(optarg);
      break;
    case 51:
      strcpy(FILE_LIG_PARAM, optarg);
      break;
    case 52:
      strcpy(FILE_LIG_TOPO, optarg);
      break;
    case 53:
      strcpy(FILE_LIG_CATACONS, optarg);
      break;
    case 54:
      MAX_NUM_OF_RUNS = atoi(optarg);
      break;
    case 55:
      strcpy(REFERENCE_RESIDUES, optarg);
      break;
    case 56:
      DISTANCE_TO_REFERENCE_RESIDUES = atof(optarg);
      break;
    case 57:
      NTRAJ_START_NDX = atoi(optarg);
      break;
    case 58:
      WGT_BIND = atof(optarg);
      break;
    case 59:
      strcpy(RESIDUE, optarg);
      break;
    case 60:
      strcpy(RESIDUE_PAIR, optarg);
      break;
    case 61:
      strcpy(EXCLUDE_RESIDUES, optarg);
      break;
    case 62:
      strcpy(FILE_LIG_PLACING, optarg);
      break;
    case 37:
      strcpy(OUTPUT_PREFIX, optarg);
      break;
    case 64:
      MIN_MUT_NUM = atoi(optarg);
      FLAG_RESTRICT_NUM_OF_MUTATIONS = TRUE;
      break;
    case 65:
      MAX_MUT_NUM = atoi(optarg);
      FLAG_RESTRICT_NUM_OF_MUTATIONS = TRUE;
      break;
    case 66:
      FLAG_PENALIZE_IDENTICAL_SEQUENCES = TRUE;
      break;
    case 67:
      SA_T_FINAL = atof(optarg);
      break;
    case 68:
      SA_T_INITIAL = atof(optarg);
      break;
    case 69:
      SA_CYCLE = atoi(optarg);
      break;
    case 70:
      NUM_METROPOLIS_STEPS = atoi(optarg);
      break;
    case 71:
      SA_DECREASE_FACTOR = atof(optarg);
      break;
    case 72:
      FLAG_WRITE_MUTATIONS = TRUE;
      break;
    default:
      sprintf(errMsg, "in file %s line %d, unknown option, program exits", __FILE__, __LINE__);
      TraceError(errMsg, ValueError);
      exit(ValueError);
      break;
    }
  }

  // terminate if no command detected
  if (strcmp(cmdname, "") == 0)
  {
    printf("no command name detected, please choose a command for execution, program exits\n");
    printf("run './%s -h' for help\n", PROGRAM_NAME);
    exit(Success);
  }
  else if (strcmp(cmdname, "CheckClash1") == 0)
  {
    sprintf(FILE_ATOMPARAM, "%s/library/param_checkclash1.prm", PROGRAM_PATH);
  }
  else if (strcmp(cmdname, "CheckClash2") == 0)
  {
    sprintf(FILE_ATOMPARAM, "%s/library/param_checkclash2.prm", PROGRAM_PATH);
  }

  // set file names
  sprintf(FILE_SELF_ENERGY, "%s_selfenergy.txt", OUTPUT_PREFIX);
  sprintf(FILE_ROTLIST, "%s_rotlist.txt", OUTPUT_PREFIX);
  sprintf(FILE_ROTLIST_SEC, "%s_rotlistSEC.txt", OUTPUT_PREFIX);
  sprintf(TAG_BEST_DESIGNS, "%s_designs.tsv", OUTPUT_PREFIX);
  sprintf(TAG_BEST_SEQS, "%s_seqs.tsv", OUTPUT_PREFIX);
  sprintf(TAG_BEST_MODEL, "%s_mod", OUTPUT_PREFIX);
  sprintf(TAG_BEST_LIG_MOL2, "%s_lig", OUTPUT_PREFIX);

  // read atom parameters
  AtomParams atomParam;
  AtomParamsCreate(&atomParam);
  if (FAILED(AtomParamsRead(&atomParam, FILE_ATOMPARAM)))
  {
    AtomParamsDestroy(&atomParam);
    sprintf(errMsg, "in file %s line %d, failed to read parameter file %s", __FILE__, __LINE__, FILE_ATOMPARAM);
    TraceError(errMsg, IOError);
    exit(IOError);
  }

  // read residue topology
  ResTopos resiTopo;
  ResToposCreate(&resiTopo);
  if (FAILED(ResToposRead(&resiTopo, FILE_TOPO)))
  {
    ResToposDestroy(&resiTopo);
    sprintf(errMsg, "in file %s line %d, failed to read topology file %s", __FILE__, __LINE__, FILE_TOPO);
    TraceError(errMsg, IOError);
    exit(IOError);
  }

  Structure structure;
  StructureCreate(&structure);
  if (FLAG_PDB == TRUE)
  {
    if (FAILED(StructureReadPDB(&structure, PDB, &atomParam, &resiTopo)))
    {
      sprintf(errMsg, "in file %s line %d, failed to read pdb file %s", __FILE__, __LINE__, PDB);
      TraceError(errMsg, IOError);
      exit(IOError);
    }
    else
    {
      printf("PDB file %s has been read\n", PDB);
      StructureCalcPhiPsi(&structure);
    }
  }

  // read mol2 and generate parameter and topology files
  if (FLAG_MOL2 == TRUE)
  {
    if (strcmp(cmdname, "MakeLigParamAndTopo") == 0)
    {
      if (FAILED(GenerateSmallMolParameterAndTopologyFromMol2(MOL2, FILE_LIG_PARAM, FILE_LIG_TOPO, INI_ATOM1, INI_ATOM2, INI_ATOM3)))
      {
        sprintf(errMsg, "in file %s line %d, failed to generate ligand atom parameter or ligand topology from mol2 %s, program exits", __FILE__, __LINE__, MOL2);
        TraceError(errMsg, ValueError);
        exit(ValueError);
      }
      else
      {
        printf("ligand atom parameter file '%s' and topology file '%s' were generated\n", FILE_LIG_PARAM, FILE_LIG_TOPO);
        exit(Success);
      }
    }

    if (FAILED(AtomParamsRead(&atomParam, FILE_LIG_PARAM)))
    {
      sprintf(errMsg, "in file %s line %d, failed to read ligand atom parameters from %s, program exits. please check!", 
      __FILE__, __LINE__, FILE_LIG_PARAM);
      TraceError(errMsg, IOError);
      exit(IOError);
    }

    if (FAILED(ResToposRead(&resiTopo, FILE_LIG_TOPO)))
    {
      sprintf(errMsg, "in file %s line %d, failed to read ligand topology from %s, program exits. please check!", 
      __FILE__, __LINE__, FILE_LIG_TOPO);
      TraceError(errMsg, IOError);
      exit(IOError);
    }

    if (FAILED(StructureReadMol2(&structure, MOL2, &atomParam, &resiTopo)))
    {
      sprintf(errMsg, "in file %s line %d, failed to read mol2 file %s", __FILE__, __LINE__, MOL2);
      TraceError(errMsg, IOError);
      exit(IOError);
    }
    else
    {
      printf("MOL2 file %s has been read\n", MOL2);
    }
  }

  int totResNum = 0;
  for (int i = 0; i < StructureGetChainNum(&structure); i++)
  {
    if (ChainGetType(StructureGetChain(&structure, i)) == Type_Chain_Protein)
    {
      totResNum += ChainGetResidueNum(StructureGetChain(&structure, i));
    }
  }
  PROT_LEN_NORM = totResNum > 0 ? totResNum : 1;

  // set the default rotamer library
  sprintf(FILE_ROTLIB, "%s/library/rotlib/dun2010bb3per.lib", PROGRAM_PATH);
  sprintf(FILE_ROTLIB_BIN, "%s/library/rotlib/ALLbbdep.bin", PROGRAM_PATH);

  // read energy weights from file
  if (FAILED(EnergyWeightRead(FILE_WEIGHT_READ)))
  {
    printf("failed to read weights from file. all energy weights are set to 1\n");
  }

  // apply evolution information
  if (FLAG_EVOLUTION == TRUE && strlen(DES_CHAINS) <= 1)
  {
    StructureDeployEvolutionInfo(&structure);
  }

  if (strcmp(cmdname, "ProteinDesign") == 0)
  {
    // set protein chains for design
    if (strcmp(DES_CHAINS, "") == 0)
    {
      int protChnNum = 0;
      for (int i = 0; i < StructureGetChainNum(&structure); i++)
      {
        Chain *pChain = StructureGetChain(&structure, i);
        if (ChainGetType(pChain) == Type_Chain_Protein)
        {
          protChnNum++;
        }
      }
      if (protChnNum >= 1)
      {
        strcpy(DES_CHAINS, ChainGetName(StructureGetChain(&structure, 0)));
      }
      else
      {
        sprintf(errMsg, "in file %s line %d, protein chain(s) %s were not identified and job failed; program exit", __FILE__, __LINE__, DES_CHAINS);
        TraceError(errMsg, ValueError);
        exit(ValueError);
      }
    }

    AAppTable aapptable;
    RamaTable ramatable;
    BBdepRotLib bbrotlib;
    AApropensityTableReadFromFile(&aapptable, FILE_ENERGY_AAPP);
    RamaTableReadFromFile(&ramatable, FILE_ENERGY_RAMA);
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    StructureCalcAminoAcidPropensityAndRamaEnergy(&structure, &aapptable, &ramatable);
    StructureCalcAminoAcidDunbrackEnergy(&structure, &bbrotlib);
    if (FLAG_DESIGN_PPI == TRUE && FLAG_INTERFACE_ONLY == TRUE)
    {
      if (FLAG_USE_RESFILE == TRUE)
        StructureBuildResfileRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
      else
        StructureBuildPPIRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo);
    }
    else if (FLAG_DESIGN_PLI == TRUE)
    {
      StructureBuildPLIShell1RotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
      StructureBuildPLIShell2RotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
    }
    else if (FLAG_DESIGN_ENZ == TRUE)
    {
      if (FLAG_USE_RESFILE == TRUE)
      {
        StructureBuildCatalyticRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
        StructureBuildPLIShell1RotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
        StructureBuildPLIShell2RotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
      }
      else
      {
        sprintf(errMsg, "in file %s line %d, catalytic sites must be specified in RESFILE, otherwise catalytic sites may be mutated by the program\n"
                        "this situation should be avoided in a meaningful enzyme design process",
                __FILE__, __LINE__);
        TraceError(errMsg, DataNotExistError);
        return DataNotExistError;
      }
    }
    else
    {
      if (FLAG_USE_RESFILE == TRUE)
      {
        StructureBuildResfileRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
      }
      else
      {
        StructureBuildAllRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo);
      }
    }
    BBdepRotamerLibDestroy(&bbrotlib);

    // deal with ligand rots if applicable
    if (FLAG_DESIGN_PLI == TRUE || FLAG_DESIGN_ENZ == TRUE)
    {
      if (FLAG_USE_LIG_POSES && !FAILED(StructureReadSmallMolRotamers(&structure, &resiTopo, FILE_LIG_POSES_READ)))
        printf("read ligand poses from %s\n", FILE_LIG_POSES_READ);
      else
      {
        printf("use the ligand pose in mol2 file for design\n");
        FILE *pOut = fopen(FILE_LIG_POSES_WRITE, "w");
        Model(1, pOut);
        Residue *pSmallMol = NULL;
        StructureFindSmallMol(&structure, &pSmallMol);
        AtomArrayShowInPDBFormat(ResidueGetAllAtoms(pSmallMol), "ATOM", ResidueGetName(pSmallMol), ResidueGetChainName(pSmallMol), 1, ResidueGetPosInChain(pSmallMol), pOut);
        EndModel(pOut);
        fclose(pOut);
        StructureReadSmallMolRotamers(&structure, &resiTopo, FILE_LIG_POSES_WRITE);
      }
    }

    RotamerList rotList;
    RotamerListCreateFromStructure(&rotList, &structure);
    RotamerListWrite(&rotList, FILE_ROTLIST);
    StructureShowDesignSites2(&structure, &rotList);
    SelfEnergyGenerate2(&structure, &aapptable, &ramatable, FILE_SELF_ENERGY);
    SelfEnergyReadAndCheck(&structure, &rotList, FILE_SELF_ENERGY);
    RotamerListWrite(&rotList, FILE_ROTLIST_SEC);
    RotamerListRead(&rotList, FILE_ROTLIST_SEC);
    StructureShowDesignSites2(&structure, &rotList);
    SimulatedAnnealing(&structure, &rotList);
    RotamerListDestroy(&rotList);
  }

  else if (strcmp(cmdname, "ComputeStability") == 0)
  {
    double energyTerms[MAX_ENERGY_TERM] = {0};
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable, FILE_ENERGY_AAPP);
    RamaTableReadFromFile(&ramatable, FILE_ENERGY_RAMA);
    BBdepRotLib bbrotlib;
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    ComputeStability(&structure, &aapptable, &ramatable, &bbrotlib, energyTerms);
    BBdepRotamerLibDestroy(&bbrotlib);
  }

  else if (strcmp(cmdname, "ComputeBinding") == 0)
  {
    if (StructureGetChainNum(&structure) <= 2 || FLAG_CHAIN_SPLIT == FALSE)
      ComputeBinding(&structure);
    else
      ComputeBinding2(&structure, SPLIT_PART1, SPLIT_PART2);
  }

  else if (strcmp(cmdname, "RepairStructure") == 0)
  {
    BBdepRotLib bbrotlib;
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    StructureCalcAminoAcidDunbrackEnergy(&structure, &bbrotlib);
    RepairStructure(&structure, &bbrotlib, &atomParam, &resiTopo);
    BBdepRotamerLibDestroy(&bbrotlib);
  }

  else if (strcmp(cmdname, "Minimize") == 0)
  {
    BBdepRotLib bbrotlib;
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    StructureCalcAminoAcidDunbrackEnergy(&structure, &bbrotlib);
    EnergyMinimization(&structure, &bbrotlib, &atomParam, &resiTopo);
    BBdepRotamerLibDestroy(&bbrotlib);
  }

  else if (strcmp(cmdname, "BuildMutant") == 0)
  {
    BBdepRotLib bbrotlib;
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    StructureCalcAminoAcidDunbrackEnergy(&structure, &bbrotlib);
    BuildMutant(&structure, MUTANT_FILE, &bbrotlib, &atomParam, &resiTopo);
    BBdepRotamerLibDestroy(&bbrotlib);
  }

  else if (strcmp(cmdname, "ComputeResEnergy") == 0)
  {
    if (strcmp(RESIDUE, "") == 0)
    {
      printf("residue was not specified. Please use '--residue=ABB' to specify the residue, "
        "where A is the chain ID, BB is the residue position in the chain\n");
      exit(Warning);
    }

    char chainName[2] = "";
    chainName[0] = RESIDUE[0];
    int pos = atoi(RESIDUE + 1);
    int chainIndex = -100;
    StructureFindChainIndex(&structure, chainName, &chainIndex);
    Chain *pChain = NULL;
    if (chainIndex != -100)
      pChain = StructureGetChain(&structure, chainIndex);
    int resNdx = -1;
    ChainFindResidueIndexByPosInChain(pChain, pos, &resNdx);
    Residue *pResi = NULL;
    if (resNdx != -1)
      pResi = ChainGetResidue(pChain, resNdx);
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable, FILE_ENERGY_AAPP);
    RamaTableReadFromFile(&ramatable, FILE_ENERGY_RAMA);
    BBdepRotLib bbrotlib;
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    printf("residue %s%d%s energy details:\n", ChainGetName(pChain), ResidueGetPosInChain(pResi), ResidueGetName(pResi));
    ComputeResidueInteractionWithFixedEnvironment(&structure, &aapptable, &ramatable, &bbrotlib, chainIndex, resNdx);
    BBdepRotamerLibDestroy(&bbrotlib);
  }

  else if (!strcmp(cmdname, "ComputeResPairEnergy"))
  {
    if (strcmp(RESIDUE_PAIR, "") == 0)
    {
      printf("residue was not specified. Please use '--residue_pair=ABB,CDD' to specify the residue, "
        "where A/C is the chain ID, BB/DD is the residue position in the chain\n");
      exit(Warning);
    }
    char part1[MAX_LEN_CHAIN_NAME + 1] = "";
    char part2[MAX_LEN_CHAIN_NAME + 1] = "";
    sscanf(RESIDUE_PAIR, "%[^,],%s", part1, part2);
    char chainName1[2] = "";
    char chainName2[2] = "";
    chainName1[0] = part1[0];
    chainName2[0] = part2[0];
    int pos1 = atoi(part1 + 1);
    int pos2 = atoi(part2 + 1);
    int chainIndex1 = -100, chainIndex2 = -100;
    StructureFindChainIndex(&structure, chainName1, &chainIndex1);
    StructureFindChainIndex(&structure, chainName2, &chainIndex2);
    Chain *pChain1 = NULL;
    Chain *pChain2 = NULL;
    if (chainIndex1 != -100 && chainIndex2 != -100)
    {
      pChain1 = StructureGetChain(&structure, chainIndex1);
      pChain2 = StructureGetChain(&structure, chainIndex2);
    }
    int resNdx1 = -1, resNdx2 = -1;
    ChainFindResidueIndexByPosInChain(pChain1, pos1, &resNdx1);
    ChainFindResidueIndexByPosInChain(pChain2, pos2, &resNdx2);
    Residue *pRes1 = NULL;
    Residue *pRes2 = NULL;
    if (resNdx1 != -1 && resNdx2 != -1)
    {
      pRes1 = ChainGetResidue(pChain1, resNdx1);
      pRes2 = ChainGetResidue(pChain2, resNdx2);
    }
    double energyTerms[MAX_ENERGY_TERM] = {0};
    if (strcmp(chainName1, chainName2) == 0)
    { // same chain
      if (pos1 == pos2 + 1 || pos2 == pos1 + 1)
        EnergyResidueAndNextResidue(pRes1, pRes2, energyTerms);
      else
        EnergyResidueAndOtherResidueSameChain(pRes1, pRes2, energyTerms);
    }
    else
    { // different chains
      if (ChainGetType(pChain1) == Type_Chain_SmallMol)
        EnergyResidueAndLigandResidue(pRes2, pRes1, energyTerms);
      else if (ChainGetType(pChain2) == Type_Chain_SmallMol)
        EnergyResidueAndLigandResidue(pRes1, pRes2, energyTerms);
      else
        EnergyResidueAndOtherResidueDiffChain(pRes1, pRes2, energyTerms);
    }
    printf("interaction energy between residue %s%d%s and residue %s%d%s:\n",
           ChainGetName(pChain1), ResidueGetPosInChain(pRes1), ResidueGetName(pRes1),
           ChainGetName(pChain2), ResidueGetPosInChain(pRes2), ResidueGetName(pRes2));
    EnergyTermWeighting(energyTerms);
    EnergyTermShowComplex(energyTerms);
  }

  else if (strcmp(cmdname, "ComputeRotamersEnergy") == 0)
  {
    AAppTable aapptable;
    RamaTable ramatable;
    AApropensityTableReadFromFile(&aapptable, FILE_ENERGY_AAPP);
    RamaTableReadFromFile(&ramatable, FILE_ENERGY_RAMA);
    BBdepRotLib bbrotlib;
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    if (FLAG_WILDTYPE_ONLY == TRUE)
      ComputeWildtypeRotamersEnergy(&structure, &bbrotlib, &aapptable, &ramatable, &atomParam, &resiTopo);
    else
      ComputeRotamersEnergy(&structure, &bbrotlib, &aapptable, &ramatable, &atomParam, &resiTopo);
    BBdepRotamerLibDestroy(&bbrotlib);
  }

  else if (strcmp(cmdname, "ShowResComposition") == 0)
  {
    int aas[30] = {0};
    ComputeAminoAcidComposition(&structure, aas);
  }

  else if (strcmp(cmdname, "OptimizeHydrogen") == 0)
  {
    OptimizeHydrogen(&structure, &atomParam, &resiTopo);
  }

  else if (strcmp(cmdname, "AddPolarHydrogen") == 0)
  {
    AddPolarHydrogen(&structure);
  }

  else if (strcmp(cmdname, "FindInterfaceRes") == 0)
  {
    if (StructureGetChainNum(&structure) < 2)
    {
      printf("there is only one chain in the whole structure, no protein-protein interface found\n");
      exit(ValueError);
    }
    else if (StructureGetChainNum(&structure) == 2)
    {
      FindInterfaceResidues(&structure);
    }
    else
    {
      printf("there are >= 3 chains in the whole structure, please split chains before calculating interface residues\n");
      FindInterfaceResidues2(&structure, SPLIT_PART1, SPLIT_PART2);
    }
  }

  else if (strcmp(cmdname, "SelectResWithin") == 0)
  {
    SelectResiduesInRange(&structure);
  }

  else if (strcmp(cmdname, "MakeLigPoses") == 0)
  {
    if (strcmp(FILE_LIG_POSES_WRITE, "") == 0)
    {
      sprintf(errMsg, "in file %s line %d, no file is specified to save ligand poses, please use e.g.:\n"
                      "--write_lig_poses=LIG_POSES.pdb to save the ligand poses into a file named 'LIG_POSES.pdb'",
              __FILE__, __LINE__);
      TraceError(errMsg, InvalidInputError);
      exit(InvalidInputError);
    }

    BBdepRotLib bbrotlib;
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    StructureBuildResfileRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo, FILE_RESFILE);
    BBdepRotamerLibDestroy(&bbrotlib);
    StructureShowDesignSites(&structure);
    StructureGenerateSmallMolRotamers(&structure, FILE_LIG_CATACONS, FILE_LIG_PLACING);
    StructureWriteSmallMolRotamers(&structure, FILE_LIG_POSES_WRITE);
    StructureShowDesignSites(&structure);
  }

  else if (strcmp(cmdname, "AnalyzeLigPoses") == 0)
  {
    if (strcmp(FILE_LIG_POSES_READ, "") == 0)
    {
      sprintf(errMsg, "in file %s line %d, the file recording ligand poses for analysis is not specified, please use e.g.:\n"
                      "--read_lig_poses=LIG_POSES.pdb to read ligand poses from the file named 'LIG_POSES.pdb'",
              __FILE__, __LINE__);
      TraceError(errMsg, InvalidInputError);
      exit(InvalidInputError);
    }
    Residue *pSmallmol = NULL;
    StructureFindSmallMol(&structure, &pSmallmol);
    AnalyzeSmallMolRotamers(FILE_LIG_POSES_READ, pSmallmol);
  }

  else if (strcmp(cmdname, "ScreenLigPoses") == 0)
  {
    if (FLAG_BY_ORIENTATION)
      StructureSmallmolOrientationScreen(&structure, &resiTopo, FILE_LIG_POSES_READ, FILE_LIG_POSES_WRITE, FILE_PICK_POSE_BY_ORITENTATION);
    else if (FLAG_BY_VDW_PERCENTILE)
      SmallMolRotamersGetBothHighRankOfBackboneVdwAndInternalVdw(FILE_LIG_POSES_READ, FILE_LIG_POSES_WRITE, PICK_POSE_BBVDW_RANK_CUTOFF, PICK_POSE_INVDW_RANK_CUTOFF);
    else if (FLAG_BY_RMSD)
      ScreenSmallmolRotamersByRMSD(FILE_LIG_POSES_READ, FILE_LIG_POSES_WRITE, PICK_POSE_RMSD_CUTOFF);
  }

  else if (strcmp(cmdname, "ComputeEvolutionScore") == 0)
  {
    double evoscore = EvolutionScorePrfFromFile(TGT_SEQ);
    printf("evoscore: %f\n", evoscore);
  }

  else if (strcmp(cmdname, "OptimizeWeight") == 0)
  {
    PNATAA_WeightOptByGradientDescent(FILE_PDBLIST);
  }

  else if (strcmp(cmdname, "CalcPhiPsi") == 0)
  {
    char PHI_PSI_FILE[MAX_LEN_FILE_NAME + 1] = "";
    sprintf(PHI_PSI_FILE, "%s_phipsi.txt", PDBID);
    WritePhiPsi(&structure, PHI_PSI_FILE);
  }

  else if (strcmp(cmdname, "CheckResInRotLib") == 0)
  {
    BBdepRotLib bbrotlib;
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    CheckRotamerInBBdepRotLib(&structure, &bbrotlib, &resiTopo, TORSION_DEVIATION_THRESHOLD);
    BBdepRotamerLibDestroy(&bbrotlib);
  }

  else if (strcmp(cmdname, "CalcResMinRMSD") == 0)
  {
    BBdepRotLib bbrotlib;
    BBdepRotamerLibCreate2(&bbrotlib, FILE_ROTLIB_BIN);
    StructureGenerateWildtypeRotamersByBBdepRotLib(&structure, &bbrotlib, &atomParam, &resiTopo);
    BBdepRotamerLibDestroy(&bbrotlib);
    FindMinRmsdRotFromRotLib(&structure);
  }

  else if (strcmp(cmdname, "CompareSideChain") == 0)
  {
    Structure newStruct;
    StructureCreate(&newStruct);
    StructureReadPDB(&newStruct, PDB2, &atomParam, &resiTopo);
    StructureCalcPhiPsi(&newStruct);
    StructureCalcProteinResidueSidechainTorsion(&newStruct, &resiTopo);
    FILE *pTorsion = fopen("UD_rmsdComparison.txt", "w");
    FILE *pRMSD = fopen("UD_torsionComparison.txt", "w");
    CompareTwoStructureProteinSidechains(&structure, &newStruct, pTorsion, pRMSD);
    fclose(pTorsion);
    fclose(pRMSD);
    StructureDestroy(&newStruct);
  }

  else if (strcmp(cmdname, "CheckClash0") == 0)
  {
    CheckClash0(&structure, CLASH_RATIO_THRESHOLD);
  }

  else if (strcmp(cmdname, "CheckClash1") == 0)
  {
    CheckClash1(&structure, CLASH_RATIO_THRESHOLD);
  }

  else if (strcmp(cmdname, "CheckClash2") == 0)
  {
    CheckClash2(&structure, CLASH_RATIO_THRESHOLD);
  }

  else if (strcmp(cmdname, "PredSS") == 0)
  {
    SSPred(TGT_SEQ);
  }

  else if (strcmp(cmdname, "PredSA") == 0)
  {
    SAPred(TGT_SEQ);
  }

  else if (strcmp(cmdname, "PredPhiPsi") == 0)
  {
    PhiPsiPred(TGT_SEQ);
  }

  else if (strcmp(cmdname, "CalcResBfactor") == 0)
  {
    char FILE_BFACTOR[MAX_LEN_FILE_NAME + 1];
    sprintf(FILE_BFACTOR, "%s_bfactor.txt", PDBID);
    FILE *pFile = fopen(FILE_BFACTOR, "w");
    for (int i = 0; i < StructureGetChainNum(&structure); i++)
    {
      Chain *pChain = StructureGetChain(&structure, i);
      if (ChainGetType(pChain) != Type_Chain_Protein)
        continue;
      for (int j = 0; j < ChainGetResidueNum(pChain); j++)
      {
        Residue *pResi = ChainGetResidue(pChain, j);
        double bfac = ResidueGetAverageBfactor(pResi);
        fprintf(pFile, "%3s %8.3f\n", ResidueGetName(pResi), bfac);
      }
    }
    fclose(pFile);
  }

  else if (strcmp(cmdname, "FindCoreRes") == 0)
  {
    FindCoreResidues(&structure);
  }

  else if (strcmp(cmdname, "FindSurfaceRes") == 0)
  {
    FindSurfaceResidues(&structure);
  }

  else if (strcmp(cmdname, "FindIntermediateRes") == 0)
  {
    FindIntermediateResidues(&structure);
  }

  else if (strcmp(cmdname, "WriteFirstLigPoseToMol2") == 0)
  {
    StructureReadSmallMolRotamers(&structure, &resiTopo, FILE_LIG_POSES_READ);
    StructureWriteFirstLigandRotamerIntoMol2(&structure, "FirstLigPose.mol2");
  }

  else if (strcmp(cmdname, "WriteFasta") == 0)
  {
    char seqFileName[MAX_LEN_FILE_NAME + 1];
    sprintf(seqFileName, "%s_seq.fa", OUTPUT_PREFIX);
    FILE *pFile = fopen(seqFileName, "w");
    for (int i = 0; i < StructureGetChainNum(&structure); i++)
    {
      Chain *pChn = StructureGetChain(&structure, i);
      fprintf(pFile, ">%s%s\n", PDBID, ChainGetName(pChn));
      for (int j = 0; j < ChainGetResidueNum(pChn); j++)
      {
        Residue *pRes = ChainGetResidue(pChn, j);
        if (pRes->isSCIntact == TRUE)
          fprintf(pFile, "%c", tolower(AA3ToAA1(ResidueGetName(pRes))));
        else
          fprintf(pFile, "%c", AA3ToAA1(ResidueGetName(pRes)));
      }
      fprintf(pFile, "\n");
    }
    fclose(pFile);
  }

  else if (strcmp(cmdname, "FindMutations") == 0)
  {
    Structure newStruct;
    StructureCreate(&newStruct);
    StructureReadPDB(&newStruct, PDB2, &atomParam, &resiTopo);
    StructureCalcPhiPsi(&newStruct);
    StructureCalcProteinResidueSidechainTorsion(&newStruct, &resiTopo);
    char mutStr[2 * MAX_LEN_FILE_NAME + 1] = "";
    StructureFindMutations(&structure, &newStruct, mutStr);
    printf("mutations = %s\n", mutStr);
    StructureDestroy(&newStruct);
  }

  else
  {
    printf("unknown command name: %s, program exits\n", cmdname);
    exit(ValueError);
  }

  if (FLAG_EVOLUTION == TRUE)
  {
    Chain *pChain = StructureFindChainByName(&structure, DES_CHAINS);
    int len = ChainGetResidueNum(pChain);
    if (PROT_PROFILE != NULL)
    {
      for (int i = 0; i < len; i++)
      {
        delete[] PROT_PROFILE[i];
      }
      delete[] PROT_PROFILE;
    }
  }
  ResToposDestroy(&resiTopo);
  AtomParamsDestroy(&atomParam);
  StructureDestroy(&structure);

  clock_t timeend = clock();
  SpentTimeShow(timestart, timeend);

  return Success;
}
