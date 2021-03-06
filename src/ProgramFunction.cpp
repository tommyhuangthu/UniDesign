/*******************************************************************************************************************************
Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

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

#include "ProgramFunction.h"
#include "RotamerBuilder.h"
#include "RotamerOptimizer.h"
#include "Evolution.h"
#include <string.h>
#include <ctype.h>

extern BOOL FLAG_MONOMER;
extern BOOL FLAG_PPI;
extern BOOL FLAG_ENZYME;
extern BOOL FLAG_PROT_LIG;

extern BOOL FLAG_SHOW_HYDROGEN;

extern BOOL FLAG_USE_INPUT_SC;
extern BOOL FLAG_ROTATE_HYDROXYL;
extern double CUT_TORSION_DEVIATION;

extern double CUT_EXCL_LOW_ROT_PROB;
extern char PROGRAM_NAME[MAX_LENGTH_FILE_NAME+1];
extern char PROGRAM_VERSION[MAX_LENGTH_FILE_NAME+1];

extern int CUT_NUM_CB_CORE;
extern int CUT_NUM_CB_SURF;

extern double CUT_PPI_DIST_SHELL1;



int PrintHelp(){
  printf(  
    "Usage: %s [OPTIONS]\n\n"  
    "OPTIONS:\n"
    "   -v              print version\n"
    "   --version       print version\n"
    "   -h              print help message\n"
    "   -?              print help message\n"
    "   --help          print help message\n"  
    "   --command=arg   choose your command:\n"
    "                   [Energy calculation module]\n"
    "                     ComputeStability\n"
    "                     ComputeBinding\n"
    "                     ComputeEvolutionScore\n"
    "                     ComputeResiEnergy\n"
    "                     ComputeRotamersEnergy\n"
    "                   [Structure modeling module]\n"
    "                     RepairStructure\n"
    "                     Minimize\n"
    "                     BuildMutant\n"
    "                     AddPolarHydrogen\n"
    "                     OptimizeHydrogen\n"
    "                     ProteinDesign\n"
    "                   [ProteinDesign and side-chain packing module]\n"
    "                     monomer design: ProteinDesign --monomer\n"
    "                     PPI design    : ProteinDesign --ppint\n"
    "                     PLI design    : ProteinDesign --protlig\n"
    "                     enzyme design : ProteinDesign --enzyme\n"
    "                     packing       : ProteinDesign --wildtype_only\n"
    "                   [Enzyme and PLI design module]\n"
    "                     GenLigParamAndTopo\n"
    "                     MakeLigEnsemble\n"
    "                     AnalyzeLigEnsemble\n"
    "                     ScreenLigEnsemble\n"
    "                   [Structure checking and assessment module]\n"
    "                     CheckClash0\n"
    "                     CheckClash1\n"
    "                     CheckClash2\n"
    "                     CheckResiInRotLib\n"
    "                     CompareSideChain\n"
    "                     FindCoreResidue\n"
    "                     FindInterfaceResidue\n"
    "                     FindIntermediateResidue\n"
    "                     FindSurfaceResidue\n"
    "                     GetPhiPsi\n"
    "                     GetResiMinRMSD\n"
    "                     GetResiBfactor\n"
    "                     ShowResiComposition\n"
    "                   [Sequence-based feature prediction module]\n"
    "                     PhiPsiPred\n"
    "                     SSPred\n"
    "                     SAPred\n"
    "                   [Other modules]\n"
    "                     OptimizeWeight\n"
    "\n"
    "command can be used in conjunction with one or more following options:\n\n"
    "   --pdb=arg                 arg is a standard readable PDB file\n"
    "   --split_chains=arg        arg specifies chain splitting for ComputeBinding with >=3 protein chains (e.g. AB,C or AB,CDE)\n"
    "   --mutant_file=arg         arg is a file recording one or more mutants\n"
    "   --bbdep=arg               arg = yes/no, turning on/off the flag for using bb-dep rotamer library (default: yes)\n"
    "   --use_input_sc            use input protein structure's side chains as rotamers\n"
    "   --ppi_shell1=arg          arg is the distance cutoff for the 1st shell of PPI (default: 5.0 Angstroms)\n"
    "   --ppi_shell2=arg          arg is the distance cutoff for the 2nd shell of PPI (default: 8.0 Angstroms)\n"
    "   --evolution               protein design using both physics and evolution scores\n"
    "   --physics                 protein design using physical scoring function (EvoEF) only (default)\n"
    "   --monomer                 the designing protein is treated as a monomer\n"
    "   --ppint                   protein-protein interaction design\n"
    "   --protlig                 protein-ligand interaction design\n"
    "   --enzyme                  enzyme design (extended version of PLI design) with catalytic constraints\n"
    "   --design_chains=arg       arg is the designing chain(s)'s identifier (e.g. A, or AB)\n"
    "   --mol2=arg                arg is a standard readable mol2 file, atomic charge should be included\n"
    "   --rotate_hydroxyl=arg     arg = yes/no, turning on/off the flag for rotating Ser, Thr, and Tyr hydroxyl hydrogen (default: on)\n"
    "   --xdeviation=arg          arg is a float value cutoff for side-chain Chi angle deviation\n"
    "   --wread=arg               arg is a energy-weight file\n"
    "   --pdblist=arg             arg is a file recording a list of PDB IDs, each in one line\n"
    "   --wildtype_only\n"        
    "   --rotlib=arg              arg is the name of a rotamer lib\n"
    "   --pdb2=arg                arg is the 2nd PDB file for comparing side-chains by command CompareSideChain\n"
    "                             this flag makes design very slow (default: turned off)\n"
    "   --pli_shell1=arg          arg is the distance cutoff for the 1st shell of protein-ligand interaction (default: 5.0 Angstroms)\n"
    "   --pli_shell2=arg          arg is the distance cutoff for the 2nd shell of protein-ligand interaction (default: 8.0 Angstroms)\n"
    "   --clash_ratio=arg         arg is a float value cutoff for the command CheckClash[0-2] (default: 0.6)\n"
    "   --ntraj=arg               arg is an integer for the number of independent protein design trajectories (default: 1)\n"
    "   --excl_low_prob=arg       arg is a flat value cutoff for excluding low-probability rotamers (default: 0.03), 0~0.05 suggested\n"
    "   --interface_only\n"       
    "   --seq=arg                 arg is a single-line plain-text FASTA protein sequence file\n"
    "   --evo_allterms            turn on the flag for using all-term evolution (including PSSM, SA, SS, and psi-phi) for protein design\n"
    "   --wprof=arg               arg is a float value (default: 1.0)\n"
    "   --resfile=arg             arg is the name of residue-restriction file\n"
    "   --seed_from_nat           turn on the flag for using native sequence as a starting point for protein design\n"
    "   --excl_rot_cys            exclude the Cystenine rotamers for protein design\n"
    "   --show_hydrogen=arg       arg = yes, write polar hydrogen's coordinates into PDB\n"
    "                             arg = no, do not write hydrogens into PDB\n"
    "   --ncut_cb_core=arg        arg is an integer cutoff value for the number of CB atoms for a CORE (buried) residue (default: 20)\n"
    "   --ncut_cb_surf=arg        arg is an integer cutoff value for the number of CB atoms for a SURFACE (exposed) residue (default: 15)\n"
    "   --init_3atoms=arg         arg is the name of three atoms divided by commas (e.g. C1,C2,C3)\n"
    "                             this option is only used with the command GenLigParamAndTopo\n"
    "   --read_lig_ensemble=arg   read ligand poses from the file arg, a PDB file recording multiple ligand poses\n"
    "   --write_lig_ensemble=arg  write ligand poses to the file arg, a PDB file for recording multiple ligand poses\n"
    "   --scrn_by_ornt=arg        screen ligand conformers using the rule defined in the arg file\n"
    "   --scrn_by_vdwpct=arg      screen ligand conformers with both internalVDW and backboneVDW ranked in a low percentile, e.g. 25%\n"
    "   --scrn_by_rmsd=arg        screen ligand conformers using an RMSD cutoff value (default: 1.0)\n"
    "\n\n",
    PROGRAM_NAME);
  return Success;
}

int PrintVersion(){
  printf("%s v%s\n",PROGRAM_NAME,PROGRAM_VERSION);
  return Success;
}


int PrintAdvertisement(){
  printf(
    "#############################################################################################################\n"
    "                                       %s v%s\n"
    "  A command-line-based tool featuring protein structure modeling and protein design, including:\n\n"
    "    *) de novo protein sequence design and protein sequence redesign (ProteinDesign)\n"
    "       can be used to design monomer proteins, protein-protein interactions, protein-ligand interactions, \n"
    "       protein-nucleic acid interactions, and enzymes\n"
    "    *) protein side-chain packing (ProteinDesign in conjuncion with the option --wildtype_only)\n"
    "    *) calculate protein fold stability score (ComputeStability)\n"
    "    *) calculate protein-protein binding energy score (ComputeBinding)\n"
    "    *) repair missing protein side chains (RepairStructure)\n"
    "    *) protein structure energy minimization (Minimize)\n"
    "    *) build mutant strucure models (BuildMutant)\n"
    "    *) add polar hydrogen atoms (AddPolarHydrogen)\n"
    "    *) optimize hydrogen's position to maximize hydrogen bonding networks (OptimizeHydrogen)\n"
    "\n"
    "  Copyright (c) 2020 Xiaoqiang Huang\n"
    "  tommyhuangthu@foxmail.com; xiaoqiah@umich.edu\n"
    "  Dept. of Computational Medicine & Bioinformatics\n"
    "  Medical School\n"
    "  University of Michigan\n"
    "############################################################################################################\n\n",
    PROGRAM_NAME,PROGRAM_VERSION);
  return Success;
}


///////////////////////////////////////////////////////////////////////////////////////////
//MAIN functions
//////////////////////////////////////////////////////////////////////////////////////////
int ComputeChainStability(Structure *pStructure,int chainIndex,double *energyTerms){
  Chain *pChainI = StructureGetChain(pStructure, chainIndex);
  for(int ir=0; ir<ChainGetResidueCount(pChainI); ir++){
    Residue *pResIR = ChainGetResidue(pChainI,ir);
    AminoAcidReferenceEnergy(pResIR->name, energyTerms);
    EnergyIntraResidue(pResIR,energyTerms);
    for(int is=ir+1; is<ChainGetResidueCount(pChainI); is++){
      Residue *pResIS = ChainGetResidue(pChainI,is);
      if(ResidueGetPosInChain(pResIR)+1==ResidueGetPosInChain(pResIS)) EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
      else EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
    }
  }

  EnergyTermWeighting(energyTerms);
  printf("Chain %s energy details:\n", ChainGetName(pChainI));
  EnergyTermShowMonomer(energyTerms);

  return Success;
}


int ComputeStructureStability(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRama,double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    Chain *pChainI=StructureGetChain(pStructure,i);
    for(int ir=0; ir<ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      if(ChainGetType(pChainI)==Type_Chain_Protein){
        AminoAcidReferenceEnergy(ResidueGetName(pResIR),energyTerms);
        EnergyIntraResidue(pResIR,energyTerms);
        AminoAcidPropensityAndRamachandranEnergy(pResIR,pAAppTable,pRama);
        energyTerms[91]+=pResIR->aapp;
        energyTerms[92]+=pResIR->rama;
      }
      for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
        Residue *pResIS = ChainGetResidue(pChainI,is);
        if(ResidueGetPosInChain(pResIR)+1==ResidueGetPosInChain(pResIS)) EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
        else EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
      }
      for(int k = i+1; k < StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        for(int ks = 0; ks < ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          if(ChainGetType(pChainI)==Type_Chain_SmallMol){
            EnergyResidueAndLigResidue(pResKS,pResIR,energyTerms);
          }
          else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
            EnergyResidueAndLigResidue(pResIR,pResKS,energyTerms);
          }
          else{
            EnergyResidueAndOtherResidueDiffChain(pResIR,pResKS,energyTerms);
          }
        }
      }
    }
  }

  EnergyTermWeighting(energyTerms);
  printf("\nStructure energy details:\n");
  EnergyTermShowComplex(energyTerms);
  
  return Success;
}


int ComputeStructureStabilityByBBdepRotLib(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRama,BBdepRotamerLib* pRotLib,double energyTerms[MAX_ENERGY_TERM]){
  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int ir=0; ir<ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      if(ChainGetType(pChainI)==Type_Chain_Protein){
        AminoAcidReferenceEnergy(ResidueGetName(pResIR),energyTerms);
        EnergyIntraResidue(pResIR,energyTerms);
        AminoAcidPropensityAndRamachandranEnergy(pResIR,pAAppTable,pRama);
        AminoAcidDunbrackEnergy(pResIR,pRotLib);
        energyTerms[91]+=pResIR->aapp;
        energyTerms[92]+=pResIR->rama;
        energyTerms[93]+=ResidueGetDunbrack(pResIR);
      }
      for(int is=ir+1; is<ChainGetResidueCount(pChainI); is++){
        Residue *pResIS = ChainGetResidue(pChainI,is);
        if(ResidueGetPosInChain(pResIR)+1==ResidueGetPosInChain(pResIS)) EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
        else EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
      }
      for(int k=i+1; k<StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        for(int ks=0; ks<ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          if(ChainGetType(pChainI)==Type_Chain_SmallMol) EnergyResidueAndLigResidue(pResKS,pResIR,energyTerms);
          else if(ChainGetType(pChainK)==Type_Chain_SmallMol) EnergyResidueAndLigResidue(pResIR,pResKS,energyTerms);
          else EnergyResidueAndOtherResidueDiffChain(pResIR,pResKS,energyTerms);
        }
      }
    }
  }

  EnergyTermWeighting(energyTerms);
  printf("\nStructure energy details:\n");
  EnergyTermShowComplex(energyTerms);
  return Success;
}


int ComputeStructureStabilityByBBdepRotLib2(Structure *pStructure,AAppTable* pAAppTable,RamaTable* pRama,char* binLib,double energyTerms[MAX_ENERGY_TERM]){
  // ACDEFGHIKLMNPQRSTVWY, only for regular amino acid
  int xcount[20]={0,1,2,3,2,0,2,2,4,2,3,2,2,3,4,1,1,1,2,2};
  int nrot[20]={0,3,18,54,18,0,36,9,73,9,27,36,2,108,75,3,3,3,36,18};
  int lrot[20]={0,129,111,240,448,0,294,330,348,339,421,75,466,132,0,468,471,528,474,510};
  FILE* pFileLib=fopen(binLib,"rb");
  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    if(ChainGetType(pChainI)==Type_Chain_Protein){
      for(int ir=0;ir<ChainGetResidueCount(pChainI);ir++){
        Residue* pResiIR=ChainGetResidue(pChainI,ir);
        if(pResiIR->isSCIntact){// calculate dunbrack energy for the amino acids with complete side chains
          int aaIdx=ThreeLetterAAGetIndex(ResidueGetName(pResiIR));
          if(aaIdx==0 || aaIdx==5) continue;//skip Ala and Gly
          int nchi=xcount[aaIdx];
          int binIdx=((int)(pResiIR->phipsi[0]+180)/10)*36+(int)(pResiIR->phipsi[1]+180)/10;
          fseek(pFileLib,(1296*lrot[aaIdx]+binIdx*nrot[aaIdx])*36,SEEK_SET);
          int matchIdx=-1;
          double pMin=10;
          double pMatch;
          for(int rotIdx=0;rotIdx<nrot[aaIdx];rotIdx++){
            float p=0.;
            fread((char*)&p,sizeof(char),4,pFileLib);
            if(p<CUT_EXCL_LOW_ROT_PROB) break;
            if(p<pMin){pMin=p;}
            float val;
            double libtorsions[4]={0};
            double libdeviations[4]={0};
            for(int k=0;k<nchi;k++){
              fread((char*)&val,sizeof(char),4,pFileLib);
              libtorsions[k]=(float)val;
            }
            fseek(pFileLib,(4-nchi)*4,SEEK_CUR);
            for(int k=0;k<nchi;k++){
              fread((char*)&val,sizeof(char),4,pFileLib);
              libdeviations[k]=(float)val;
            }
            fseek(pFileLib,(4-nchi)*4,SEEK_CUR);

            BOOL match=TRUE;
            for(int torIdx=0;torIdx<DoubleArrayGetLength(&pResiIR->Xs);torIdx++){
              double min=DegToRad(libtorsions[torIdx])-DegToRad(libdeviations[torIdx]);
              double max=DegToRad(libtorsions[torIdx])+DegToRad(libdeviations[torIdx]);
              double torsion=DoubleArrayGet(&pResiIR->Xs,torIdx);
              double torsionm2pi=torsion-2*PI;
              double torsionp2pi=torsion+2*PI;
              double torsion2=torsion;
              if((strcmp(ResidueGetName(pResiIR),"PHE")==0 && torIdx==1)|| (strcmp(ResidueGetName(pResiIR),"TYR")==0 && torIdx==1)||
                (strcmp(ResidueGetName(pResiIR),"ASP")==0 && torIdx==1)|| strcmp(ResidueGetName(pResiIR),"GLU")==0 && torIdx==2){
                  torsion2=torsion+PI;
                  torsion2=torsion>0?torsion-PI:torsion2;
              }
              double torsion2m2pi=torsion2-2*PI;
              double torsion2p2pi=torsion2+2*PI;
              if(!((torsion<=max && torsion>=min) || (torsionm2pi<=max && torsionm2pi>=min) || (torsionp2pi<=max && torsionp2pi>=min) ||
                (torsion2<=max && torsion2>=min) || (torsion2m2pi<=max && torsion2m2pi>=min) || (torsion2p2pi<=max && torsion2p2pi>=min))){
                  match=TRUE;
              }
            }
            if(match){
              matchIdx=rotIdx;
              pMatch=p;
              break;
            }
          }
          double pDelta=1e-7;
          if(matchIdx!=-1) energyTerms[93]+=-log(pMatch+pDelta);
          else energyTerms[93]+=-log(pMin+pDelta);
        }
      }
    }
  }
  fclose(pFileLib);

  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int ir=0; ir<ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      if(ChainGetType(pChainI)==Type_Chain_Protein){
        AminoAcidReferenceEnergy(pResIR->name, energyTerms);
        EnergyIntraResidue(pResIR,energyTerms);
        AminoAcidPropensityAndRamachandranEnergy(pResIR,pAAppTable,pRama);
        energyTerms[91]+=pResIR->aapp;
        energyTerms[92]+=pResIR->rama;
      }
      for(int is=ir+1;is<ChainGetResidueCount(pChainI);is++){
        Residue *pResIS = ChainGetResidue(pChainI,is);
        if(ResidueGetPosInChain(pResIR)+1==ResidueGetPosInChain(pResIS)) EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
        else EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
      }
      for(int k=i+1; k<StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        for(int ks=0; ks<ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          if(ChainGetType(pChainI)==Type_Chain_SmallMol) EnergyResidueAndLigResidue(pResKS,pResIR,energyTerms);
          else if(ChainGetType(pChainK)==Type_Chain_SmallMol) EnergyResidueAndLigResidue(pResIR,pResKS,energyTerms);
          else EnergyResidueAndOtherResidueDiffChain(pResIR,pResKS,energyTerms);
        }
      }
    }
  }

  EnergyTermWeighting(energyTerms);
  printf("\nStructure energy details:\n");
  EnergyTermShowComplex(energyTerms);
  return Success;
}


int ComputeBinding(Structure *pStructure){
  if(StructureGetChainCount(pStructure)>2){
    printf("Your structure has more than two protein chains, and you should specify how to split chains "
      "before computing the binding energy. Otherwise, the interactions between any chain pair is calculated\n");
  }
  else if(StructureGetChainCount(pStructure)<=1){
    printf("Your structure has less than or equal to one chain, binding energy cannot be calculated\n");
    return Warning;
  }

  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChainI=StructureGetChain(pStructure,i);
    for(int k=i+1; k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      double energyTerms[MAX_ENERGY_TERM]={0};
      for(int j=0;j<ChainGetResidueCount(pChainI);j++){
        Residue* pResIJ=ChainGetResidue(pChainI,j);
        for(int s=0;s<ChainGetResidueCount(pChainK);s++){
          Residue* pResKS=ChainGetResidue(pChainK,s);
          if(ChainGetType(pChainI)==Type_Chain_SmallMol){
            EnergyResidueAndLigResidue(pResKS,pResIJ,energyTerms);
          }
          else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
            EnergyResidueAndLigResidue(pResIJ,pResKS,energyTerms);
          }
          else{
            EnergyResidueAndOtherResidueDiffChain(pResIJ,pResKS,energyTerms);
          }
        }
      }
      EnergyTermWeighting(energyTerms);
      // energy terms are weighted during the calculation, don't weight them for the difference
      printf("Binding energy details between chain %s and chain %s:\n",ChainGetName(pChainI),ChainGetName(pChainK));
      EnergyTermShowComplex(energyTerms);
    }
  }

  return Success;
}


int ComputeBindingByChainSplitting(Structure *pStructure,char split1[], char split2[]){
  double energyTerms[MAX_ENERGY_TERM]={0};
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChainI=StructureGetChain(pStructure,i);
    for(int k=i+1;k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      if((strstr(split1,ChainGetName(pChainI))!=NULL && strstr(split2,ChainGetName(pChainK))!=NULL)||
        ((strstr(split2,ChainGetName(pChainI))!=NULL && strstr(split1,ChainGetName(pChainK))!=NULL))){
          for(int j=0;j<ChainGetResidueCount(pChainI);j++){
            Residue* pResIJ=ChainGetResidue(pChainI,j);
            for(int s=0;s<ChainGetResidueCount(pChainK);s++){
              Residue* pResKS=ChainGetResidue(pChainK,s);
              if(ChainGetType(pChainI)==Type_Chain_SmallMol){
                EnergyResidueAndLigResidue(pResKS,pResIJ,energyTerms);
              }
              else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
                EnergyResidueAndLigResidue(pResIJ,pResKS,energyTerms);
              }
              else{
                EnergyResidueAndOtherResidueDiffChain(pResIJ,pResKS,energyTerms);
              }
            }
          }
      }
    }
  }
  EnergyTermWeighting(energyTerms);
  printf("Binding energy details between chain(s) %s and chain(s) %s (DG_bind = DG(stability,complex) - DG(stability,%s) - DG(stability,%s):\n",split1,split2,split1,split2);
  EnergyTermShowComplex(energyTerms);
  return Success;
}


//this function is used to build the structure model of mutations
int BuildMutant(Structure* pStructure, char* mutantfile, BBindRotamerLib* rotlib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  FileReader fr;
  FileReaderCreate(&fr, mutantfile);
  int mutantcount = FileReaderGetLineCount(&fr);
  if(mutantcount<=0){
    printf("There is no mutant found in the mutant file\n");
    return DataNotExistError;
  }

  StringArray* mutants = (StringArray*)malloc(sizeof(StringArray)*mutantcount);
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  int mutantIndex=0;
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArrayCreate(&mutants[mutantIndex]);
    StringArraySplitString(&mutants[mutantIndex], line, ',');
    char lastMutant[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    int lastmutindex = StringArrayGetCount(&mutants[mutantIndex])-1;
    strcpy(lastMutant, StringArrayGet(&mutants[mutantIndex], lastmutindex));
    //deal with the last char of the last single mutant
    if((!isdigit(lastMutant[strlen(lastMutant)-1])) && !isalpha(lastMutant[strlen(lastMutant)-1])){
      lastMutant[strlen(lastMutant)-1] = '\0';
    }
    StringArraySet(&mutants[mutantIndex], lastmutindex, lastMutant);
    mutantIndex++;
  }
  FileReaderDestroy(&fr);

  for(int mutantIndex = 0; mutantIndex < mutantcount; mutantIndex++){
    Structure tempStruct;
    StructureCreate(&tempStruct);
    StructureCopy(&tempStruct,pStructure);
    //for each mutant, build the rotamer-tree
    IntArray mutantArray,rotamersArray;
    IntArrayCreate(&mutantArray,0);
    IntArrayCreate(&rotamersArray,0);
    for(int cycle=0; cycle<StringArrayGetCount(&mutants[mutantIndex]); cycle++){
      char mutstr[10];
      char aa1, chn, aa2;
      int posInChain;
      strcpy(mutstr, StringArrayGet(&mutants[mutantIndex], cycle));
      sscanf(mutstr, "%c%c%d%c", &aa1, &chn, &posInChain, &aa2);
      int chainIndex = -1, residueIndex = -1;
      char chainname[MAX_LENGTH_CHAIN_NAME]; chainname[0] = chn; chainname[1] = '\0';
      StructureFindChainIndex(&tempStruct,chainname,&chainIndex);
      if(chainIndex==-1){
        printf("in file %s line %d, cannot find mutation %s\n", __FILE__,__LINE__,mutstr);
        exit(ValueError);
      }
      ChainFindResidueByPosInChain(StructureGetChain(&tempStruct, chainIndex), posInChain, &residueIndex);
      if(residueIndex==-1){
        printf("in file %s line %d, cannot find mutation %s\n", __FILE__,__LINE__,mutstr);
        exit(ValueError);
      }
      char mutaatype[MAX_LENGTH_RESIDUE_NAME];
      OneLetterAAToThreeLetterAA(aa2, mutaatype);
      StringArray designType, patchType;
      StringArrayCreate(&designType);
      StringArrayCreate(&patchType);
      // for histidine, the default mutaatype is HSD, we need to add HSE
      StringArrayAppend(&designType, mutaatype); StringArrayAppend(&patchType, "");
      if(aa2=='H'){StringArrayAppend(&designType, "HSE"); StringArrayAppend(&patchType, "");}
      ProteinSiteBuildMutatedRotamers(&tempStruct,chainIndex,residueIndex,rotlib,atomParams,resiTopos,&designType,&patchType);
      IntArrayAppend(&mutantArray, chainIndex);
      IntArrayAppend(&mutantArray, residueIndex);
      IntArrayAppend(&rotamersArray,chainIndex);
      IntArrayAppend(&rotamersArray,residueIndex);
      StringArrayDestroy(&designType);
      StringArrayDestroy(&patchType);
    }

    // for each mutant, find the surrounding residues and build the wild-type rotamer-tree
    for(int ii=0; ii<IntArrayGetLength(&mutantArray); ii+=2){
      int chainIndex = IntArrayGet(&mutantArray,ii);
      int resiIndex = IntArrayGet(&mutantArray,ii+1);
      Residue *pResi1 = ChainGetResidue(StructureGetChain(&tempStruct, chainIndex), resiIndex);
      for(int j = 0; j < StructureGetChainCount(&tempStruct); ++j){
        Chain* pChain = StructureGetChain(&tempStruct,j);
        for(int k=0; k<ChainGetResidueCount(pChain); k++){
          Residue* pResi2 = ChainGetResidue(pChain,k);
          if(AtomArrayCalcMinDistance(&pResi1->atoms,&pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
            if(pResi2->designType==Type_ResidueDesignType_Fixed){
              ProteinSiteBuildWildtypeRotamers(&tempStruct,j,k,rotlib,atomParams,resiTopos);
              if(pResi2->isSCIntact) ProteinSiteBuildNativeRotamer(&tempStruct,j,k,resiTopos);
              IntArrayAppend(&rotamersArray,j);
              IntArrayAppend(&rotamersArray,k);
            }
          }
        }
      }
    }

    // optimization rotamers sequentially
    printf("Building mutation model %d, the following sites will be optimized:\n",mutantIndex+1);
    printf("chnIndex resIndex (both of them starts from zero on the chain)\n");
    for(int ii=0;ii<IntArrayGetLength(&rotamersArray);ii+=2){
      printf("%8d %8d\n",IntArrayGet(&rotamersArray,ii),IntArrayGet(&rotamersArray,ii+1));
    }
    for(int cycle=0; cycle<10; cycle++){
      printf("optimization cycle %d ... \n",cycle+1);
      for(int ii=0; ii<IntArrayGetLength(&rotamersArray); ii+=2){
        int chainIndex = IntArrayGet(&rotamersArray, ii);
        int resiIndex = IntArrayGet(&rotamersArray, ii+1);
        ProteinSiteOptimizeRotamer(&tempStruct, chainIndex, resiIndex);
      }
    }
    IntArrayDestroy(&mutantArray);
    IntArrayDestroy(&rotamersArray);
    //remember to delete rotamers for previous mutant
    StructureRemoveAllDesignSites(&tempStruct);

    char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    if(pdbid!=NULL)
      sprintf(modelfile,"%s_Model_%04d.pdb",pdbid,mutantIndex+1);
    else
      sprintf(modelfile,"Model_%04d.pdb",mutantIndex+1);
    FILE* pf=fopen(modelfile,"w");
    fprintf(pf,"REMARK file generated by module <BuildMutant>\n");
    StructureShowInPDBFormat(&tempStruct,pf);
    fclose(pf);
    StructureDestroy(&tempStruct);
  }

  return Success;
}


int BuildMutantByBBdepRotLib(Structure* pStructure,char* mutantfile,BBdepRotamerLib* pBBdepRotLib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  FileReader fr;
  FileReaderCreate(&fr, mutantfile);
  int mutantcount = FileReaderGetLineCount(&fr);
  if(mutantcount<=0){
    printf("in file %s line %d, no mutation was found in file %s\n",__FILE__,__LINE__,mutantfile);
    exit(DataNotExistError);
  }
  StringArray* mutants = (StringArray*)malloc(sizeof(StringArray)*mutantcount);
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  int mutantIndex=0;
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArrayCreate(&mutants[mutantIndex]);
    StringArraySplitString(&mutants[mutantIndex], line, ',');
    char lastMutant[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    int lastmutindex = StringArrayGetCount(&mutants[mutantIndex])-1;
    strcpy(lastMutant, StringArrayGet(&mutants[mutantIndex], lastmutindex));
    //deal with the last char of the last single mutant
    if((!isdigit(lastMutant[strlen(lastMutant)-1])) && !isalpha(lastMutant[strlen(lastMutant)-1])){
      lastMutant[strlen(lastMutant)-1] = '\0';
    }
    StringArraySet(&mutants[mutantIndex], lastmutindex, lastMutant);
    mutantIndex++;
  }
  FileReaderDestroy(&fr);

  for(int mutantIndex = 0;mutantIndex<mutantcount;mutantIndex++){
    Structure tempStruct;
    StructureCreate(&tempStruct);
    StructureCopy(&tempStruct,pStructure);
    // 1.deal with mutations on the nucleotide
    IntArray flagNucleotides;
    IntArrayCreate(&flagNucleotides,StringArrayGetCount(&mutants[mutantIndex]));
    for(int posIndex=0; posIndex<StringArrayGetCount(&mutants[mutantIndex]); posIndex++){
      char mutstr[10];
      char aa1,chn,aa2;
      int posInChain;
      strcpy(mutstr,StringArrayGet(&mutants[mutantIndex],posIndex));
      sscanf(mutstr,"%c%c%d%c",&aa1,&chn,&posInChain,&aa2);
      int chIdx = -1, resIdx = -1;
      char chName[MAX_LENGTH_CHAIN_NAME]; chName[0]=chn; chName[1]='\0';
      StructureFindChainIndex(&tempStruct,chName,&chIdx);
      if(chIdx==-1){
        printf("in file %s line %d, cannot find mutation %s\n", __FILE__,__LINE__,mutstr);
        exit(ValueError);
      }
      ChainFindResidueByPosInChain(StructureGetChain(&tempStruct, chIdx), posInChain, &resIdx);
      if(resIdx==-1){
        printf("in file %s line %d, cannot find mutation %s\n", __FILE__,__LINE__,mutstr);
        exit(ValueError);
      }
      Chain* pChain=StructureFindChainByName(&tempStruct,chName);
      if(ChainGetType(pChain)==Type_Chain_DNA || ChainGetType(pChain)==Type_Chain_RNA){
        Residue* pTargetNucleotide=ChainGetResidue(StructureGetChain(&tempStruct, chIdx),resIdx);
        Residue* pPairedNucleotide=StructureFindNucleotidePair(&tempStruct,pTargetNucleotide);
        StructureMutateNucleotide(pTargetNucleotide,aa2,&tempStruct,atomParams,resiTopos);
        if(pPairedNucleotide!=NULL){
          char aa22;
          if(ChainGetType(StructureFindChainByName(&tempStruct,ResidueGetChainName(pPairedNucleotide)))==Type_Chain_DNA){
            if(aa2=='a') aa22='t';
            else if(aa2=='g') aa22='c';
            else if(aa2=='c') aa22='g';
            else if(aa2=='t') aa22='a';
            else if(aa2=='u') aa22='a';
            else if(aa2=='n') aa22='n';
          }
          else if(ChainGetType(StructureFindChainByName(&tempStruct,ResidueGetChainName(pPairedNucleotide)))==Type_Chain_RNA){
            if(aa2=='a') aa22='u';
            else if(aa2=='g') aa22='c';
            else if(aa2=='c') aa22='g';
            else if(aa2=='t') aa22='a';
            else if(aa2=='u') aa22='a';
            else if(aa2=='n') aa22='n';
          }
          StructureMutateNucleotide(pPairedNucleotide,aa22,&tempStruct,atomParams,resiTopos);
        }
        IntArraySet(&flagNucleotides,posIndex,1);
      }
    }

    // 2.for each mutant, build the rotamer-tree
    IntArray mutatedArray,rotamericArray;
    IntArrayCreate(&mutatedArray,0);
    IntArrayCreate(&rotamericArray,0);
    for(int posIndex=0; posIndex<StringArrayGetCount(&mutants[mutantIndex]); posIndex++){
      if(IntArrayGet(&flagNucleotides,posIndex)==1) continue;
      char mutstr[10];
      char aa1,chn,aa2;
      int posInChain;
      strcpy(mutstr,StringArrayGet(&mutants[mutantIndex],posIndex));
      sscanf(mutstr,"%c%c%d%c",&aa1,&chn,&posInChain,&aa2);
      int chainIndex = -1, residueIndex = -1;
      char chainname[MAX_LENGTH_CHAIN_NAME]; chainname[0]=chn; chainname[1]='\0';
      StructureFindChainIndex(&tempStruct,chainname,&chainIndex);
      if(chainIndex==-1){
        printf("in file %s line %d, cannot find mutation %s\n", __FILE__,__LINE__,mutstr);
        exit(ValueError);
      }
      ChainFindResidueByPosInChain(StructureGetChain(&tempStruct, chainIndex), posInChain, &residueIndex);
      if(residueIndex==-1){
        printf("in file %s line %d, cannot find mutation %s\n", __FILE__,__LINE__,mutstr);
        exit(ValueError);
      }
      char mutaatype[MAX_LENGTH_RESIDUE_NAME];
      OneLetterAAToThreeLetterAA(aa2, mutaatype);
      StringArray designType, patchType;
      StringArrayCreate(&designType);
      StringArrayCreate(&patchType);
      //for histidine, the default mutaatype is HSD, we need to add HSE
      StringArrayAppend(&designType, mutaatype); StringArrayAppend(&patchType, "");
      if(aa2=='H'){StringArrayAppend(&designType, "HSE"); StringArrayAppend(&patchType, "");}
      ProteinSiteBuildMutatedRotamersByBBdepRotLib(&tempStruct,chainIndex,residueIndex,pBBdepRotLib,atomParams,resiTopos,&designType,&patchType);
      IntArrayAppend(&mutatedArray, chainIndex);
      IntArrayAppend(&mutatedArray, residueIndex);
      IntArrayAppend(&rotamericArray,chainIndex);
      IntArrayAppend(&rotamericArray,residueIndex);
      StringArrayDestroy(&designType);
      StringArrayDestroy(&patchType);
    }
    IntArrayDestroy(&flagNucleotides);

    // 3.build rotamers for surrounding residues
    for(int ii=0; ii<IntArrayGetLength(&mutatedArray); ii+=2){
      int chainIndex = IntArrayGet(&mutatedArray,ii);
      int resiIndex = IntArrayGet(&mutatedArray,ii+1);
      Residue *pResi1 = ChainGetResidue(StructureGetChain(&tempStruct, chainIndex), resiIndex);
      for(int j = 0; j < StructureGetChainCount(&tempStruct); ++j){
        Chain* pChain = StructureGetChain(&tempStruct,j);
        if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
        for(int k=0; k<ChainGetResidueCount(pChain); k++){
          Residue* pResi2 = ChainGetResidue(pChain,k);
          if(AtomArrayCalcMinDistance(&pResi1->atoms,&pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
            if(pResi2->designType==Type_ResidueDesignType_Fixed){
              ProteinSiteBuildWildtypeRotamersByBBdepRotLib(&tempStruct,j,k,pBBdepRotLib,atomParams,resiTopos);
              if(pResi2->isSCIntact) ProteinSiteBuildNativeRotamer(&tempStruct,j,k,resiTopos);
              IntArrayAppend(&rotamericArray,j);
              IntArrayAppend(&rotamericArray,k);
            }
          }
        }
      }
    }

    // 4.optimization rotamers sequentially
    printf("Building mutation model %d, the following sites will be optimized:\n",mutantIndex+1);
    printf("chnIndex resIndex (both of them starts from zero on the chain)\n");
    for(int ii=0;ii<IntArrayGetLength(&rotamericArray);ii+=2){
      printf("%8d %8d\n",IntArrayGet(&rotamericArray,ii),IntArrayGet(&rotamericArray,ii+1));
    }
    for(int posIndex=0; posIndex<10; posIndex++){
      printf("optimization cycle %d ... \n",posIndex+1);
      for(int ii=0; ii<IntArrayGetLength(&rotamericArray); ii+=2){
        int chainIndex = IntArrayGet(&rotamericArray, ii);
        int resiIndex = IntArrayGet(&rotamericArray, ii+1);
        ProteinSiteOptimizeRotamerWithBBdepRotLib(&tempStruct,chainIndex,resiIndex,pBBdepRotLib);
      }
    }
    IntArrayDestroy(&mutatedArray);
    IntArrayDestroy(&rotamericArray);
    //remember to delete rotamers for previous mutant
    StructureRemoveAllDesignSites(&tempStruct);

    char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    if(pdbid!=NULL) sprintf(modelfile,"%s_Model_%04d.pdb",pdbid,mutantIndex+1);
    else sprintf(modelfile,"Model_%04d.pdb",mutantIndex+1);
    FILE* pf=fopen(modelfile,"w");
    fprintf(pf,"REMARK file generated by module <BuildMutant>\n");
    StructureShowInPDBFormat(&tempStruct,pf);
    fclose(pf);
    StructureDestroy(&tempStruct);
  }

  return Success;
}


int RepairStructure(Structure* pStructure, BBindRotamerLib* pBBindRotLib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  for(int cycle=0; cycle<3; cycle++){
    printf("Repairing PDB cycle %d ... \n",cycle+1);
    for(int i=0; i<StructureGetChainCount(pStructure); ++i){
      Chain* pChain = StructureGetChain(pStructure, i);
      if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
      for(int j=0; j<ChainGetResidueCount(pChain); j++){
        Residue* pResi = ChainGetResidue(pChain, j);
        if(pResi->isSCIntact){
          if(strcmp(ResidueGetName(pResi),"ALA")==0 || strcmp(ResidueGetName(pResi),"GLY")==0) continue;
          //skip CYS which may form disulfide bonds
          if(strcmp(ResidueGetName(pResi),"CYS")==0) continue;
          if(strcmp(ResidueGetName(pResi),"ASN")==0||strcmp(ResidueGetName(pResi),"GLN")==0||strcmp(ResidueGetName(pResi),"HSD")==0||strcmp(ResidueGetName(pResi),"HSE")==0){
            printf("We now flip residue %s%d%s to optimize hydrogen bonds\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ResidueGetName(pResi));
            ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
            ProteinSiteBuildFlippedNativeRotamer(pStructure,i,j,resiTopos);
            ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
            ProteinSiteRemoveDesignSite(pStructure,i,j);
          }
          else if(strcmp(ResidueGetName(pResi),"SER")==0 || strcmp(ResidueGetName(pResi),"THR")==0 || strcmp(ResidueGetName(pResi),"TYR")==0){
            printf("We now rotate hydroxyl group of residue %s%d%s to optimize hydrogen bonds\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ResidueGetName(pResi));
            ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
            ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
            ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
            ProteinSiteRemoveDesignSite(pStructure,i,j);
          }
        }
        else{
          printf("We now optimize side-chain conformation of residue %s%d%s\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ResidueGetName(pResi));
          ResidueCalcAllAtomXYZ(pResi,resiTopos,ChainGetResidue(pChain,j-1),ChainGetResidue(pChain,j+1));
          ResidueCalcSidechainTorsion(pResi,resiTopos);
          ProteinSiteBuildWildtypeRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamer(pStructure,i,j);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
        
      }
    }
  }

  //output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_Repaired.pdb",pdbid);}
  else{strcpy(modelfile,"Repaired.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK file generated by module <RepairStructure>\n");
  StructureShowInPDBFormat(pStructure,pf);
  fclose(pf);

  return Success;
}


int RepairStructureByBBdepRotLib(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  for(int cycle=0; cycle<3; cycle++){
    printf("Repairing PDB cycle %d ... \n",cycle+1);
    for(int i=0; i<StructureGetChainCount(pStructure); ++i){
      Chain* pChain = StructureGetChain(pStructure, i);
      if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
      for(int j=0; j<ChainGetResidueCount(pChain); j++){
        Residue* pResi = ChainGetResidue(pChain, j);
        if(pResi->isSCIntact){
          if(strcmp(ResidueGetName(pResi),"ALA")==0 || strcmp(ResidueGetName(pResi),"GLY")==0) continue;
          //skip CYS which may form disulfide bonds
          if(strcmp(ResidueGetName(pResi),"CYS")==0) continue;
          if(strcmp(ResidueGetName(pResi),"ASN")==0||strcmp(ResidueGetName(pResi),"GLN")==0||strcmp(ResidueGetName(pResi),"HSD")==0||strcmp(ResidueGetName(pResi),"HSE")==0){
            printf("We now flip residue %s%d%s to optimize hydrogen bonds\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ResidueGetName(pResi));
            ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
            ProteinSiteBuildFlippedNativeRotamer(pStructure,i,j,resiTopos);
            ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
            ProteinSiteRemoveDesignSite(pStructure,i,j);
          }
          else if(strcmp(ResidueGetName(pResi),"SER")==0 || strcmp(ResidueGetName(pResi),"THR")==0 || strcmp(ResidueGetName(pResi),"TYR")==0){
            printf("We now rotate hydroxyl group of residue %s%d%s to optimize hydrogen bonds\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ResidueGetName(pResi));
            ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
            ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
            ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
            ProteinSiteRemoveDesignSite(pStructure,i,j);
          }
       }
        else{
          printf("We now optimize side-chain conformation of residue %s%d%s\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ResidueGetName(pResi));
          ResidueCalcAllAtomXYZ(pResi,resiTopos,ChainGetResidue(pChain,j-1),ChainGetResidue(pChain,j+1));
          //ResidueCalcSidechainTorsion(pResi,resiTopos);
          ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerWithBBdepRotLib(pStructure,i,j,pBBdepRotLib);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
    }
  }

  //output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_Repaired.pdb",pdbid);}
  else{strcpy(modelfile,"Repaired.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK file generated by module <RepairStructure>\n");
  StructureShowInPDBFormat(pStructure,pf);
  fclose(pf);

  return Success;
}


int EnergyMinimizationByBBdepRotLib(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  for(int cycle=0; cycle<3; cycle++){
    printf("Minimzation cycle %d ... \n",cycle+1);
    for(int i=0; i<StructureGetChainCount(pStructure); ++i){
      Chain* pChain = StructureGetChain(pStructure,i);
      if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
      for(int j=0; j<ChainGetResidueCount(pChain); j++){
        Residue* pResi = ChainGetResidue(pChain,j);
        if(pResi->isSCIntact){
          if(strcmp(ResidueGetName(pResi),"ALA")==0 || strcmp(ResidueGetName(pResi),"GLY")==0) continue;
          //skip CYS which may form disulfide bonds
          if(strcmp(ResidueGetName(pResi),"CYS")==0) continue;
          if(strcmp(ResidueGetName(pResi),"ASN")==0||strcmp(ResidueGetName(pResi),"GLN")==0||strcmp(ResidueGetName(pResi),"HSD")==0||strcmp(ResidueGetName(pResi),"HSE")==0){
            printf("We now flip residue %s%d%s to optimize hydrogen bonds\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ResidueGetName(pResi));
            ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
            ProteinSiteBuildFlippedNativeRotamer(pStructure,i,j,resiTopos);
            ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
            ProteinSiteRemoveDesignSite(pStructure,i,j);
          }
          else if(strcmp(ResidueGetName(pResi),"SER")==0 || strcmp(ResidueGetName(pResi),"THR")==0 || strcmp(ResidueGetName(pResi),"TYR")==0){
            printf("We now rotate hydroxyl group of residue %s%d%s to optimize hydrogen bonds\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ResidueGetName(pResi));
            ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
            ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
            ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
            ProteinSiteRemoveDesignSite(pStructure,i,j);
          }
          printf("We now optimize side-chain conformation of residue %s%d%s\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ResidueGetName(pResi));
          ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
          ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerWithBBdepRotLib(pStructure,i,j,pBBdepRotLib);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
        else{
          printf("We now optimize side-chain conformation of residue %s%d%s\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ResidueGetName(pResi));
          ResidueCalcAllAtomXYZ(pResi,resiTopos,ChainGetResidue(pChain,j-1),ChainGetResidue(pChain,j+1));
          //ResidueCalcSidechainTorsion(pResi,resiTopos);
          ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerWithBBdepRotLib(pStructure,i,j,pBBdepRotLib);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
    }
  }

  //output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_Minimized.pdb",pdbid);}
  else{strcpy(modelfile,"Minimized.pdb");}
  FILE* pOut=fopen(modelfile,"w");
  fprintf(pOut,"REMARK file generated by module <Minimize>\n");
  StructureShowInPDBFormat(pStructure,pOut);
  fclose(pOut);

  return Success;
}


int AddPolarHydrogen(Structure* pStructure, char* pdbid){
  FLAG_SHOW_HYDROGEN=TRUE;
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_PolarH.pdb",pdbid);}
  else{strcpy(modelfile,"PolarH.pdb");}
  FILE* pOut=fopen(modelfile,"w");
  fprintf(pOut,"REMARK file generated by module <AddHydrogen>\n");
  StructureShowInPDBFormat(pStructure,pOut);
  fclose(pOut);
  return Success;
}


int OptimizeHydrogen(Structure* pStructure, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  for(int cycle=0; cycle<1; cycle++){
    printf("Optimizing hydrogen cycle %d ... \n",cycle+1);
    for(int i=0; i<StructureGetChainCount(pStructure); ++i){
      Chain* pChain = StructureGetChain(pStructure, i);
      for(int j=0; j<ChainGetResidueCount(pChain); j++){
        Residue* pResi = ChainGetResidue(pChain, j);
        if(strcmp(ResidueGetName(pResi),"SER")==0 || strcmp(ResidueGetName(pResi),"THR")==0 || strcmp(ResidueGetName(pResi),"TYR")==0){
          printf("We will rotate hydroxyl group of residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          if(pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
    }
  }

  //output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_OptH.pdb",pdbid);}
  else{strcpy(modelfile,"OptH.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK file generated by module <OptimizeHydrogen>\n");
  StructureShowInPDBFormat(pStructure,pf);
  fclose(pf);

  return Success;
}


int ComputeResidueInteractionWithFixedEnvironment(Structure *pStructure, int chainIndex, int residueIndex){
  Chain *pChainI=StructureGetChain(pStructure, chainIndex);
  Residue *pResIR= ChainGetResidue(pChainI, residueIndex);
  int surroundingResiNum = 0;
  Residue **ppSurroundingResidues = NULL;
  for(int i=0; i<StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int j=0; j<ChainGetResidueCount(pChainI); j++){
      Residue *pResi2 = ChainGetResidue(pChainI,j);
      if(!strcmp(ResidueGetName(pResIR),ResidueGetName(pResi2)) && ResidueGetPosInChain(pResIR)==ResidueGetPosInChain(pResi2)) continue;
      else if(AtomArrayCalcMinDistance(&pResIR->atoms,&pResi2->atoms)<ENERGY_DISTANCE_CUTOFF){
        surroundingResiNum++;
        ppSurroundingResidues = (Residue **)realloc(ppSurroundingResidues, sizeof(Residue*)*surroundingResiNum);
        ppSurroundingResidues[surroundingResiNum-1] = pResi2;
      }
    }
  }

  double energyTerms[MAX_ENERGY_TERM]={0};
  AminoAcidReferenceEnergy(ResidueGetName(pResIR),energyTerms);
  EnergyIntraResidue(pResIR,energyTerms);
  for(int is=0; is<surroundingResiNum; is++){
    Residue *pResIS = ppSurroundingResidues[is];
    if(!strcmp(ResidueGetChainName(pResIR),ResidueGetChainName(pResIS))){
      if(ResidueGetPosInChain(pResIR)==ResidueGetPosInChain(pResIS)-1) EnergyResidueAndNextResidue(pResIR,pResIS,energyTerms);
      else if(ResidueGetPosInChain(pResIR)==ResidueGetPosInChain(pResIS)+1) EnergyResidueAndNextResidue(pResIS,pResIR,energyTerms);
      else EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,energyTerms);
    }
    else{
      if(ChainGetType(StructureFindChainByName(pStructure,ResidueGetName(pResIS)))==Type_Chain_SmallMol) EnergyResidueAndLigResidue(pResIR,pResIS,energyTerms);
      else EnergyResidueAndOtherResidueDiffChain(pResIR,pResIS,energyTerms);
    }
  }

  return Success;
}


//////////////////////////////////////////////////////////////////////////////
//The following are new program functions
//////////////////////////////////////////////////////////////////////////////
int ComputeAllRotamersEnergy(Structure* pStructure,BBindRotamerLib* pBBindRotLib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  char energyfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  sprintf(energyfile,"%s_rotenergy.txt",pdbid);
  FILE* fp=fopen(energyfile,"w");
  for(int i=0; i<StructureGetChainCount(pStructure); ++i){
    Chain* pChain = StructureGetChain(pStructure, i);
    if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
    for(int j=0; j<ChainGetResidueCount(pChain); j++){
      Residue* pResi = ChainGetResidue(pChain, j);
      if(FLAG_PPI){
        //check if residue is an interface residue
        for(int k=0;k<StructureGetChainCount(pStructure);k++){
          if(k==i) continue;
          Chain* pChainK=StructureGetChain(pStructure,k);
          for(int s=0;s<ChainGetResidueCount(pChainK);s++){
            Residue* pResiKS=ChainGetResidue(pChainK,s);
            if(AtomArrayCalcMinDistance(&pResi->atoms,&pResiKS->atoms)<ENERGY_DISTANCE_CUTOFF){
              ProteinSiteBuildAllRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
              if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
              if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
              ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
              ProteinSiteRemoveDesignSite(pStructure,i,j);
            }
          }
        }
      }
      else if(FLAG_ENZYME || FLAG_PROT_LIG){
        Residue* pSmallMol=NULL;
        StructureFindSmallMol(pStructure,&pSmallMol);
        if(AtomArrayCalcMinDistance(&pResi->atoms,&pSmallMol->atoms)<ENERGY_DISTANCE_CUTOFF){
          ProteinSiteBuildAllRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
          if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
          if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
      else{
        ProteinSiteBuildAllRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
        if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
        if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
        ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
        ProteinSiteRemoveDesignSite(pStructure,i,j);
      }
    }
  }
  fclose(fp);

  return Success;
}


int ComputeWildtypeRotamersEnergy(Structure* pStructure,BBindRotamerLib* pBBindRotLib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  char energyfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  sprintf(energyfile,"%s_rotenergy.txt",pdbid);
  FILE* fp=fopen(energyfile,"w");
  for(int i=0; i<StructureGetChainCount(pStructure); ++i){
    Chain* pChain = StructureGetChain(pStructure, i);
    if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
    for(int j=0; j<ChainGetResidueCount(pChain); j++){
      Residue* pResi = ChainGetResidue(pChain, j);
      if(FLAG_PPI){
        //check if residue is an interface residue
        for(int k=0;k<StructureGetChainCount(pStructure);k++){
          if(k==i) continue;
          Chain* pChainK=StructureGetChain(pStructure,k);
          for(int s=0;s<ChainGetResidueCount(pChainK);s++){
            Residue* pResiKS=ChainGetResidue(pChainK,s);
            if(AtomArrayCalcMinDistance(&pResi->atoms,&pResiKS->atoms)<ENERGY_DISTANCE_CUTOFF){
              ProteinSiteBuildWildtypeRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
              if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
              if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
              ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
              ProteinSiteRemoveDesignSite(pStructure,i,j);
            }
          }
        }
      }
      else if(FLAG_ENZYME || FLAG_PROT_LIG){
        Residue* pSmallMol=NULL;
        StructureFindSmallMol(pStructure,&pSmallMol);
        if(AtomArrayCalcMinDistance(&pResi->atoms,&pSmallMol->atoms)<ENERGY_DISTANCE_CUTOFF){
          ProteinSiteBuildWildtypeRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
          if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
          if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
      else{
        ProteinSiteBuildWildtypeRotamers(pStructure,i,j,pBBindRotLib,atomParams,resiTopos);
        if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
        if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
        ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
        ProteinSiteRemoveDesignSite(pStructure,i,j);
      }
    }
  }
  fclose(fp);

  return Success;
}


int FindInterfaceResidues(Structure *pStructure){
  if(StructureGetChainCount(pStructure)<2){
    printf("there is only one chain in the whole structure, no protein-protein interface found\n");
    exit(ValueError);
  }

  IntArray* chainArrays=(IntArray*)malloc(sizeof(IntArray)*StructureGetChainCount(pStructure));
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChainI=StructureGetChain(pStructure,i);
    IntArrayCreate(&chainArrays[i],pChainI->residueNum);
    for(int j=0;j<pChainI->residueNum;j++){
      IntArraySet(&chainArrays[i],j,0);
    }
  }

  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = StructureGetChain(pStructure, i);
    for(int k = i+1; k < pStructure->chainNum; k++){
      Chain *pChainK = StructureGetChain(pStructure, k);
      for(int j = 0; j < pChainI->residueNum; j++){
        Residue* pResiIJ = ChainGetResidue(pChainI, j);
        //if(IntArrayGet(&chainArrays[i],j)==1)continue;
        for(int s = 0; s < pChainK->residueNum; s++){
          Residue* pResiKS = ChainGetResidue(pChainK, s);
          //if(IntArrayGet(&chainArrays[k],s)==1)continue;
          if(AtomArrayCalcMinDistance(&pResiIJ->atoms, &pResiKS->atoms) < CUT_PPI_DIST_SHELL1){
            IntArraySet(&chainArrays[i], j, 1);
            IntArraySet(&chainArrays[k], s, 1);
          }
        }
      }
    }
  }

  printf("Interface residues: \n");
  for(int j=0;j<StructureGetChainCount(pStructure);j++){
    Chain* pChain=StructureGetChain(pStructure,j);
    if(ChainGetType(pChain)==Type_Chain_Protein || ChainGetType(pChain)==Type_Chain_DNA || ChainGetType(pChain)==Type_Chain_RNA){
      for(int i=0; i<ChainGetResidueCount(pChain); i++){
        Residue* pResidue=ChainGetResidue(pChain,i);
        if(IntArrayGet(&chainArrays[j],i)==1){
          printf("%s %4d\n",ResidueGetChainName(pResidue),ResidueGetPosInChain(pResidue));
        }
      }
    }
  }

  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    IntArrayDestroy(&chainArrays[i]);
  }
  free(chainArrays);
  chainArrays=NULL;

  return Success;
}



int FindCoreResidues(Structure* pStructure){
  StructureComputeResiduePosition(pStructure);
  printf("Core residues:\n");
  for(int j=0;j<StructureGetChainCount(pStructure);j++){
    Chain* pChain=StructureGetChain(pStructure,j);
    if(ChainGetType(pChain)==Type_Chain_Protein){
      for(int i=0; i<ChainGetResidueCount(pChain); i++){
        Residue* pResidue=ChainGetResidue(pChain,i);
        if(pResidue->nCbIn10A>CUT_NUM_CB_CORE){
          printf("%s %4d\n",ResidueGetChainName(pResidue),ResidueGetPosInChain(pResidue));
        }
      }
    }
  }
  return Success;
}


int FindSurfaceResidues(Structure* pStructure){
  StructureComputeResiduePosition(pStructure);
  printf("Surface residues:\n");
  for(int j=0;j<StructureGetChainCount(pStructure);j++){
    Chain* pChain=StructureGetChain(pStructure,j);
    if(ChainGetType(pChain)==Type_Chain_Protein){
      for(int i=0; i<ChainGetResidueCount(pChain); i++){
        Residue* pResidue=ChainGetResidue(pChain,i);
        if(pResidue->nCbIn10A<CUT_NUM_CB_SURF){
          printf("%s %4d\n",ResidueGetChainName(pResidue),ResidueGetPosInChain(pResidue));
        }
      }
    }
  }
  return Success;
}

int FindIntermediateResidues(Structure* pStructure){
  StructureComputeResiduePosition(pStructure);
  printf("Intermediate residues:\n");
  for(int j=0;j<StructureGetChainCount(pStructure);j++){
    Chain* pChain=StructureGetChain(pStructure,j);
    if(ChainGetType(pChain)==Type_Chain_Protein){
      for(int i=0; i<ChainGetResidueCount(pChain); i++){
        Residue* pResidue=ChainGetResidue(pChain,i);
        if(pResidue->nCbIn10A<=CUT_NUM_CB_CORE && pResidue->nCbIn10A>=CUT_NUM_CB_SURF){
          printf("%s %4d\n",ResidueGetChainName(pResidue),ResidueGetPosInChain(pResidue));
        }
      }
    }
  }
  return Success;
}


int ShowPhiPsi(Structure* pStructure,char* phipsifile){
  FILE* fout=NULL;
  if(phipsifile==NULL) fout=stdout;
  else fout=fopen(phipsifile,"w");
  fprintf(fout,"#phi-psi angles of protein residues: \n");
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi1=ChainGetResidue(pChain,j);
      fprintf(fout,"%s %s %4d %8.3f %8.3f\n",ResidueGetChainName(pResi1),ResidueGetName(pResi1),ResidueGetPosInChain(pResi1),pResi1->phipsi[0],pResi1->phipsi[1]);
    }
  }
  if(fout != stdout) fclose(fout);

  return Success;
}


int ComputeRotamersEnergyByBBdepRotLib(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  char energyfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  sprintf(energyfile,"%s_rotenergy.txt",pdbid);
  FILE* pFile=fopen(energyfile,"w");
  for(int i=0; i<StructureGetChainCount(pStructure); ++i){
    Chain* pChain = StructureGetChain(pStructure, i);
    if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
    for(int j=0; j<ChainGetResidueCount(pChain); j++){
      Residue* pResi = ChainGetResidue(pChain, j);
      if(FLAG_PPI){
        for(int k=0;k<StructureGetChainCount(pStructure);k++){
          if(k==i) continue;
          Chain* pChainK=StructureGetChain(pStructure,k);
          for(int s=0;s<ChainGetResidueCount(pChainK);s++){
            Residue* pResiKS=ChainGetResidue(pChainK,s);
            if(AtomArrayCalcMinDistance(&pResi->atoms,&pResiKS->atoms)<ENERGY_DISTANCE_CUTOFF){
              ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
              if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
              if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
              ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,pFile);
              ProteinSiteRemoveDesignSite(pStructure,i,j);
            }
          }
        }
      }
      else if(FLAG_ENZYME || FLAG_PROT_LIG){
        Residue* pSmallMol=NULL;
        StructureFindSmallMol(pStructure,&pSmallMol);
        if(AtomArrayCalcMinDistance(&pResi->atoms,&pSmallMol->atoms)<ENERGY_DISTANCE_CUTOFF){
          ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
          if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
          if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,pFile);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
      else{
        ProteinSiteBuildAllRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
        if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
        if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
        ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,pFile);
        ProteinSiteRemoveDesignSite(pStructure,i,j);
      }
    }
  }
  fclose(pFile);

  return Success;
}


int ComputeWildtypeRotamersEnergyByBBdepRotLib(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,AAppTable* pAAppTable,RamaTable* pRama,AtomParamsSet* atomParams,ResiTopoSet* resiTopos,char* pdbid){
  char energyfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  sprintf(energyfile,"%s_rotenergy.txt",pdbid);
  FILE* fp=fopen(energyfile,"w");
  for(int i=0; i<StructureGetChainCount(pStructure); ++i){
    Chain* pChain = StructureGetChain(pStructure, i);
    if(ChainGetType(pChain)!=Type_Chain_Protein) continue;
    for(int j=0; j<ChainGetResidueCount(pChain); j++){
      Residue* pResi = ChainGetResidue(pChain, j);
      if(FLAG_PPI){
        //check if residue is an interface residue
        for(int k=0;k<StructureGetChainCount(pStructure);k++){
          if(k==i) continue;
          Chain* pChainK=StructureGetChain(pStructure,k);
          for(int s=0;s<ChainGetResidueCount(pChainK);s++){
            Residue* pResiKS=ChainGetResidue(pChainK,s);
            if(AtomArrayCalcMinDistance(&pResi->atoms,&pResiKS->atoms)<ENERGY_DISTANCE_CUTOFF){
              ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
              if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
              if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
              ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
              ProteinSiteRemoveDesignSite(pStructure,i,j);
            }
          }
        }
      }
      else if(FLAG_ENZYME || FLAG_PROT_LIG){
        Residue* pSmallMol=NULL;
        StructureFindSmallMol(pStructure,&pSmallMol);
        if(AtomArrayCalcMinDistance(&pResi->atoms,&pSmallMol->atoms)<ENERGY_DISTANCE_CUTOFF){
          ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
          if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
          if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
          ProteinSiteRemoveDesignSite(pStructure,i,j);
        }
      }
      else{
        ProteinSiteBuildWildtypeRotamersByBBdepRotLib(pStructure,i,j,pBBdepRotLib,atomParams,resiTopos);
        if(FLAG_USE_INPUT_SC && pResi->isSCIntact) ProteinSiteBuildNativeRotamer(pStructure,i,j,resiTopos);
        if(FLAG_ROTATE_HYDROXYL) ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
        ProteinSiteCalcRotamersEnergy(pStructure,pAAppTable,pRama,i,j,fp);
        ProteinSiteRemoveDesignSite(pStructure,i,j);
      }
    }
  }
  fclose(fp);

  return Success;
}



int CheckRotamerInBBindRotLib(Structure* pStructure,BBindRotamerLib* pBBindRotLib,ResiTopoSet* pTopos,double cutoff,char* pdbid){
  char FILE_TORSION[MAX_LENGTH_FILE_NAME+1];
  sprintf(FILE_TORSION,"%s_torsionrecover.txt",pdbid);
  FILE* pf=fopen(FILE_TORSION,"w");  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResidue=ChainGetResidue(pChain,j);
      if(strcmp(ResidueGetName(pResidue),"ALA")!=0 && strcmp(ResidueGetName(pResidue),"GLY")!=0){
        BOOL result=IsNativeRotamerInBBindRotLib(pStructure,i,j,pTopos,pBBindRotLib,cutoff);
        if(result){
          fprintf(pf,"%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResidue));
        }
        else{
          fprintf(pf,"%s %s 0\n",ChainGetName(pChain),ResidueGetName(pResidue));
        }
      }
      else{
        fprintf(pf,"%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResidue));
      }
    }
  }
  fclose(pf);
  return Success;
}


int CheckRotamerInBBdepRotLib(Structure* pStructure,BBdepRotamerLib* pBBdepRotLib,ResiTopoSet* pTopos,double cutoff,char* pdbid){
  char FILE_TORSION[MAX_LENGTH_FILE_NAME+1];
  sprintf(FILE_TORSION,"%s_torsionrecover.txt",pdbid);
  FILE* pf=fopen(FILE_TORSION,"w");
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResidue=ChainGetResidue(pChain,j);
      if(strcmp(ResidueGetName(pResidue),"ALA")!=0 && strcmp(ResidueGetName(pResidue),"GLY")!=0){
        BOOL result=IsNativeRotamerInBBdepRotLib(pStructure,i,j,pTopos,pBBdepRotLib,cutoff);
        if(result){
          fprintf(pf,"%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResidue));
        }
        else{
          fprintf(pf,"%s %s 0\n",ChainGetName(pChain),ResidueGetName(pResidue));
        }
      }
      else{
        fprintf(pf,"%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResidue));
      }
    }
  }
  fclose(pf);
  return Success;
}


int FindMinRmsdRotFromRotLib(Structure* pStructure,char* pdbid){
  char FILE_RMSD[MAX_LENGTH_FILE_NAME+1];
  sprintf(FILE_RMSD,"%s_minrmsd.txt",pdbid);
  FILE* pf=fopen(FILE_RMSD,"w");
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResidue=ChainGetResidue(pChain,j);
      DesignSite* pSite=StructureFindDesignSite(pStructure,i,j);
      if(pSite!=NULL){
        double minRMSD=1e8;
        RotamerSet* pSet=DesignSiteGetRotamers(pSite);
        for(int k=0;k<RotamerSetGetCount(pSet);k++){
          Rotamer* pRotamer=RotamerSetGet(pSet,k);
          RotamerRestore(pRotamer,pSet);
          double rmsd = RotamerAndResidueSidechainRMSD(pRotamer,pResidue);
          if(RotamerIsSymmetricalCheck(pRotamer)){
            Rotamer tempRot;
            RotamerCreate(&tempRot);
            SymmetricalRotamerGenerate(&tempRot,pRotamer);
            double rmsd2=RotamerAndResidueSidechainRMSD(&tempRot,pResidue);
            if(rmsd2<rmsd) rmsd=rmsd2;
            RotamerDestroy(&tempRot);
          }
          if(rmsd<minRMSD) minRMSD=rmsd;
          RotamerExtract(pRotamer);
        }
        fprintf(pf,"%s %s %f\n",ChainGetName(pChain),ResidueGetName(pResidue),minRMSD);
      }
    }
  }
  fclose(pf);
  return Success;
}

int CompareSidechainsOf2Structures(Structure*pStructure,Structure* pStructure2,FILE* pTorsion,FILE* pRmsd){
  for(int i=0; i<StructureGetChainCount(pStructure);i++){
    Chain* pChain=StructureGetChain(pStructure,i);
    Chain* pChain2=StructureGetChain(pStructure2,i);
    for(int j=0;j<ChainGetResidueCount(pChain);j++){
      Residue* pResi=ChainGetResidue(pChain,j);
      Residue* pResi2=ChainGetResidue(pChain2,j);
      char res=ThreeLetterAAToOneLetterAA(ResidueGetName(pResi));
      char res2=ThreeLetterAAToOneLetterAA(ResidueGetName(pResi2));
      if(ChainGetType(pChain)==Type_Chain_Protein){
        if(res != 'A' && res != 'G'){
          double rmsd = ResidueAndResidueSidechainRMSD(pResi,pResi2);
          if(ResidueIsSymmetricalCheck(pResi)){
            Residue tempRot;
            ResidueCreate(&tempRot);
            SymmetricalResidueGenerate(&tempRot,pResi);
            double rmsd2=ResidueAndResidueSidechainRMSD(&tempRot,pResi2);
            if(rmsd2<rmsd) rmsd=rmsd2;
            ResidueDestroy(&tempRot);
          }
          fprintf(pRmsd, "%s %s %f\n",ChainGetName(pChain),ResidueGetName(pResi),rmsd);
#if 0
          if(ResidueAndResidueAllTorsionsAreSimilar(pResi,pResi2,TORSION_DEVIATION_CUTOFF)){
            fprintf(pTorsion, "%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResi));
          }
          else{
            fprintf(pTorsion, "%s %s 0\n",ChainGetName(pChain),ResidueGetName(pResi));
          }
#endif
#if 1
          IntArray simArray;
          IntArrayCreate(&simArray,0);
          //BOOL sim=TRUE;
          fprintf(pTorsion,"%s %s",ChainGetName(pChain),ResidueGetName(pResi));
          if(ResidueAndResidueCheckTorsionSimilarity(pResi,pResi2,CUT_TORSION_DEVIATION,&simArray)){
            fprintf(pTorsion," 1");
          }
          else{
            fprintf(pTorsion," 0");
          }
          for(int k=0;k<IntArrayGetLength(&simArray);k++){
            if(IntArrayGet(&simArray,k)==1 /*&& sim*/){
              fprintf(pTorsion," 1");
            }
            else if(IntArrayGet(&simArray,k)==2){
              fprintf(pTorsion," 0");
              //sim=FALSE;
            }
            /*else{
              fprintf(pTorsion," -");
            }*/
          }
          fprintf(pTorsion,"\n");
          IntArrayDestroy(&simArray);
#endif
        }
        else{
          double rmsd=0.0;
          fprintf(pRmsd, "%s %s %f\n",ChainGetName(pChain),ResidueGetName(pResi),rmsd);
          fprintf(pTorsion, "%s %s 1\n",ChainGetName(pChain),ResidueGetName(pResi));
        }
      }
    }
  }
  return Success;
}

int CheckClash0(Structure *pStructure, double clashRatio){
  printf("checking structure clashes with atom pairwise distance < %.2f*(ri+rj)\n",clashRatio);
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain *pChainI = StructureGetChain(pStructure, i);

    for(int j=0;j<ChainGetResidueCount(pChainI);j++){
      Residue* pResiIJ=ChainGetResidue(pChainI,j);
      for(int j2=j+1;j2<ChainGetResidueCount(pChainI);j2++){
        Residue* pResiIJ2=ChainGetResidue(pChainI,j2);
        for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
          Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
          if(pAtom1->isBBAtom || AtomIsHydrogen(pAtom1)) continue;
          for(int atom2=0;atom2<ResidueGetAtomCount(pResiIJ2);atom2++){
            Atom* pAtom2=ResidueGetAtom(pResiIJ2,atom2);
            if(pAtom2->isBBAtom || AtomIsHydrogen(pAtom2)) continue;
            double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
            if(dist<clashRatio*(pAtom1->vdw_radius+pAtom2->vdw_radius)){
              printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiIJ2),AtomGetName(pAtom2),
                dist);
            }
          }
        }
      }
    }

    for(int k=i+1;k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      for(int j=0;j<ChainGetResidueCount(pChainI);j++){
        Residue* pResiIJ=ChainGetResidue(pChainI,j);
        for(int l=0;l<ChainGetResidueCount(pChainK);l++){
          Residue* pResiKL=ChainGetResidue(pChainK,l);
          for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
            Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
            if(pAtom1->isBBAtom || AtomIsHydrogen(pAtom1)) continue;
            for(int atom2=0;atom2<ResidueGetAtomCount(pResiKL);atom2++){
              Atom* pAtom2=ResidueGetAtom(pResiKL,atom2);
              if(pAtom2->isBBAtom || AtomIsHydrogen(pAtom2)) continue;
              double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
              if(dist<clashRatio*(pAtom1->vdw_radius+pAtom2->vdw_radius)){
                if(ChainGetType(pChainI)==Type_Chain_SmallMol){
                  printf("ligand residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
                  printf("protein residue %s%d%s atom %s and ligand residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else{
                  printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
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

int CheckClash1(Structure *pStructure, double clashRatio){
  printf("checking structure clashes with atom pairwise distance < %.2f*(ri+rj)\n",clashRatio);
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain *pChainI = StructureGetChain(pStructure, i);

    for(int j=0;j<ChainGetResidueCount(pChainI);j++){
      Residue* pResiIJ=ChainGetResidue(pChainI,j);
      for(int j2=j+1;j2<ChainGetResidueCount(pChainI);j2++){
        Residue* pResiIJ2=ChainGetResidue(pChainI,j2);
        for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
          Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
          if(pAtom1->isBBAtom || AtomIsHydrogen(pAtom1)) continue;
          for(int atom2=0;atom2<ResidueGetAtomCount(pResiIJ2);atom2++){
            Atom* pAtom2=ResidueGetAtom(pResiIJ2,atom2);
            if(pAtom2->isBBAtom || AtomIsHydrogen(pAtom2)) continue;
            double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
            if(dist<clashRatio*(pAtom1->vdw_radius+pAtom2->vdw_radius)){
              printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiIJ2),AtomGetName(pAtom2),
                dist);
            }
          }
        }
      }
    }

    for(int k=i+1;k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      for(int j=0;j<ChainGetResidueCount(pChainI);j++){
        Residue* pResiIJ=ChainGetResidue(pChainI,j);
        for(int l=0;l<ChainGetResidueCount(pChainK);l++){
          Residue* pResiKL=ChainGetResidue(pChainK,l);
          for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
            Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
            if(pAtom1->isBBAtom || AtomIsHydrogen(pAtom1)) continue;
            for(int atom2=0;atom2<ResidueGetAtomCount(pResiKL);atom2++){
              Atom* pAtom2=ResidueGetAtom(pResiKL,atom2);
              if(pAtom2->isBBAtom || AtomIsHydrogen(pAtom2)) continue;
              double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
              if(dist<clashRatio*(pAtom1->vdw_radius+pAtom2->vdw_radius)){
                if(ChainGetType(pChainI)==Type_Chain_SmallMol){
                  printf("ligand residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
                  printf("protein residue %s%d%s atom %s and ligand residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else{
                  printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
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


int CheckClash2(Structure *pStructure, double clashRatio){
  printf("checking structure clashes with atom pairwise distance < %.2f*(ri+rj)\n",clashRatio);
  for(int i=0;i<StructureGetChainCount(pStructure);i++){
    Chain *pChainI = StructureGetChain(pStructure, i);

    for(int j=0;j<ChainGetResidueCount(pChainI);j++){
      Residue* pResiIJ=ChainGetResidue(pChainI,j);
      for(int j2=j+1;j2<ChainGetResidueCount(pChainI);j2++){
        Residue* pResiIJ2=ChainGetResidue(pChainI,j2);
        for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
          Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
          if(pAtom1->isBBAtom || AtomIsHydrogen(pAtom1)) continue;
          for(int atom2=0;atom2<ResidueGetAtomCount(pResiIJ2);atom2++){
            Atom* pAtom2=ResidueGetAtom(pResiIJ2,atom2);
            if(pAtom2->isBBAtom || AtomIsHydrogen(pAtom2)) continue;
            double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
            if(dist<clashRatio*(pAtom1->vdw_radius+pAtom2->vdw_radius)){
              printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiIJ2),AtomGetName(pAtom2),
                dist);
            }
          }
        }
      }
    }

    for(int k=i+1;k<StructureGetChainCount(pStructure);k++){
      Chain* pChainK=StructureGetChain(pStructure,k);
      for(int j=0;j<ChainGetResidueCount(pChainI);j++){
        Residue* pResiIJ=ChainGetResidue(pChainI,j);
        for(int l=0;l<ChainGetResidueCount(pChainK);l++){
          Residue* pResiKL=ChainGetResidue(pChainK,l);
          for(int atom1=0;atom1<ResidueGetAtomCount(pResiIJ);atom1++){
            Atom* pAtom1=ResidueGetAtom(pResiIJ,atom1);
            if(pAtom1->isBBAtom || AtomIsHydrogen(pAtom1)) continue;
            for(int atom2=0;atom2<ResidueGetAtomCount(pResiKL);atom2++){
              Atom* pAtom2=ResidueGetAtom(pResiKL,atom2);
              if(pAtom2->isBBAtom || AtomIsHydrogen(pAtom2)) continue;
              double dist=XYZDistance(&pAtom1->xyz,&pAtom2->xyz);
              if(dist<clashRatio*(pAtom1->vdw_radius+pAtom2->vdw_radius)){
                if(ChainGetType(pChainI)==Type_Chain_SmallMol){
                  printf("ligand residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else if(ChainGetType(pChainK)==Type_Chain_SmallMol){
                  printf("protein residue %s%d%s atom %s and ligand residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
                    dist);
                }
                else{
                  printf("protein residue %s%d%s atom %s and protein residue %s%d%s atom %s has clash, dist=%f\n",
                    AtomGetChainName(pAtom1),AtomGetPosInChain(pAtom1),ResidueGetName(pResiIJ),AtomGetName(pAtom1),
                    AtomGetChainName(pAtom2),AtomGetPosInChain(pAtom2),ResidueGetName(pResiKL),AtomGetName(pAtom2),
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




