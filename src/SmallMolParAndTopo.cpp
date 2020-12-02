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

#include "Utility.h"
#include "Atom.h"
#include "Residue.h"
#include "SmallMolEEF1.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


int FindAtomC(Atom *pAtomD,AtomArray *pAtomArray,IntArray *icFlags,IntArray *bondsFromToAndType,AtomArray *pAtomCArray){
	AtomArrayDestroy(pAtomCArray);
	AtomArrayCreate(pAtomCArray);
	for(int i=0; i<IntArrayGetLength(bondsFromToAndType); i+=3){
		int atomIndex = IntArrayGet(bondsFromToAndType, i);
		Atom *pAtom = AtomArrayGet(pAtomArray,atomIndex);
		if(strcmp(AtomGetName(pAtomD),AtomGetName(pAtom)) == 0){
			int atomCIndex = IntArrayGet(bondsFromToAndType,i+1);
			if(IntArrayGet(icFlags, atomCIndex) == 0) continue;
			AtomArrayAppend(pAtomCArray,AtomArrayGet(pAtomArray, atomCIndex));
		}

		atomIndex = IntArrayGet(bondsFromToAndType, i+1);
		pAtom = AtomArrayGet(pAtomArray,atomIndex);
		if(strcmp(AtomGetName(pAtomD),AtomGetName(pAtom)) == 0){
			int atomCIndex = IntArrayGet(bondsFromToAndType,i);
			if(IntArrayGet(icFlags, atomCIndex) == 0) continue;
			AtomArrayAppend(pAtomCArray, AtomArrayGet(pAtomArray, atomCIndex));
		}
	}

	return Success;
}


int FindAtomB(Atom *pAtomC,Atom* pAtomD,AtomArray *pAtoms,IntArray *icFlags,IntArray *bondsFromToAndType,AtomArray *pAtomBArray){
	AtomArrayDestroy(pAtomBArray);
	AtomArrayCreate(pAtomBArray);
	for(int i=0; i<IntArrayGetLength(bondsFromToAndType); i+=3){
		int atomIndex = IntArrayGet(bondsFromToAndType, i);
		Atom *pAtom = AtomArrayGet(pAtoms, atomIndex);
		if(strcmp(AtomGetName(pAtomC),AtomGetName(pAtom)) == 0){
			int atomCIndex = IntArrayGet(bondsFromToAndType,i+1);
			if(IntArrayGet(icFlags, atomCIndex) == 0) continue;
			if(strcmp(AtomGetName(pAtomD), AtomGetName(AtomArrayGet(pAtoms,atomCIndex))) == 0) continue;
			AtomArrayAppend(pAtomBArray, AtomArrayGet(pAtoms,atomCIndex));
		}

		atomIndex = IntArrayGet(bondsFromToAndType, i+1);
		pAtom = AtomArrayGet(pAtoms, atomIndex);
		if(strcmp(AtomGetName(pAtomC), AtomGetName(pAtom)) == 0){
			int atomCIndex = IntArrayGet(bondsFromToAndType,i);
			if(IntArrayGet(icFlags, atomCIndex) == 0) continue;
			if(strcmp(AtomGetName(pAtomD), AtomGetName(AtomArrayGet(pAtoms, atomCIndex))) == 0) continue;
			AtomArrayAppend(pAtomBArray, AtomArrayGet(pAtoms, atomCIndex));
		}
	}

	return Success;
}


int GenerateSmallMolTopologyFromMol2(char* mol2file,char* topfile,char* iniAtom1,char *iniAtom2,char* iniAtom3){
  int result = Success;
  char errMsg[MAX_LENGTH_ERR_MSG+1];
	char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
	char resiName[MAX_LENGTH_RESIDUE_NAME+1]="";
  BOOL readingAtom = FALSE;
  BOOL readingBond = FALSE;

  FILE* pMol2=fopen(mol2file,"r");
  AtomArray atoms;
  AtomArrayCreate(&atoms);
  IntArray bondsFromToAndType;
  IntArrayCreate(&bondsFromToAndType, 0);
  Atom atom;
  AtomCreate(&atom);
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pMol2)!=NULL){
    char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1]="";
    sscanf(line,"%s",keyword);
    if(strcmp(keyword,"@<TRIPOS>ATOM")==0){
      readingAtom = TRUE;
      readingBond = FALSE;
      continue;
    }
    else if(strcmp(keyword,"@<TRIPOS>BOND")==0){
      readingAtom = FALSE;
      readingBond = TRUE;
      continue;
    }
    else if(keyword[0] == '@'){
      readingAtom = readingBond = FALSE;
      continue;
    }
    if(readingAtom){
      int atomId;
      char mol2Resi[MAX_LENGTH_RESIDUE_NAME+1];
      sscanf(line,"%d %s %lf %lf %lf %s %d %s %lf",&atomId,atom.name,&atom.xyz.X,&atom.xyz.Y,&atom.xyz.Z,atom.type,&atom.posInChain,mol2Resi,&atom.charge);
      if(strcmp(resiName,"")==0){
        strcpy(resiName,mol2Resi);
        if(strlen(resiName)>3) resiName[3]='\0';
        if(LigandResidueNameConflictWithAminoAcid(resiName)){
          strcpy(resiName,"LIG");
        }
      }
      AtomArrayAppend(&atoms,&atom);
    }
    if(readingBond){
      int id, from, to, type;
      sscanf(line,"%d %d %d %d",&id,&from,&to,&type);
      IntArrayAppend(&bondsFromToAndType, from-1);
      IntArrayAppend(&bondsFromToAndType, to-1);
      IntArrayAppend(&bondsFromToAndType, type);
    }
  }
  AtomDestroy(&atom);
  fclose(pMol2);

	// write atoms and bonds
  FILE *pFileTopo = fopen(topfile,"w");
  if(pFileTopo==NULL){
    sprintf(errMsg,"in file %s line %d, cannot write to file %s",__FILE__,__LINE__,topfile);
    result = IOError;
    TraceError(errMsg,result);
    return result;
  }
	fprintf(pFileTopo, "RESI %-3.3s\n", resiName);
	for(int i=0; i<AtomArrayGetCount(&atoms); i++){
		fprintf(pFileTopo, "ATOM  %s\n", AtomGetName(AtomArrayGet(&atoms, i)));
	}
	for(int i=0; i<IntArrayGetLength(&bondsFromToAndType); i+=3){
		if(IntArrayGet(&bondsFromToAndType, i+2) == 1) fprintf(pFileTopo, "BOND  ");
		else if(IntArrayGet(&bondsFromToAndType, i+2) == 3) fprintf(pFileTopo, "TRIPLE  ");
		else fprintf(pFileTopo, "DOUBLE  ");
    fprintf(pFileTopo,"%s %s\n",AtomGetName(AtomArrayGet(&atoms, IntArrayGet(&bondsFromToAndType,i))), AtomGetName(AtomArrayGet(&atoms, IntArrayGet(&bondsFromToAndType, i+1))));
	}
	// write ICs
  fprintf(pFileTopo, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	fprintf(pFileTopo, "! Internal Coordinate correlates with small-molecule placement for enzyme design\n");
	fprintf(pFileTopo, "! It depends on three initial placing atoms on small molecule. If the three atoms\n");
  fprintf(pFileTopo, "! were not provided by user, the first three atoms from mol2 file will be used.\n");
	fprintf(pFileTopo, "! Anyway, it is very important to check the ICs manually!\n");
	fprintf(pFileTopo, "! The format of an IC entry is: \n");
  fprintf(pFileTopo, "! IC atomA atomB  atomC atomD distAB angleABC angleABCD angleBCD distCD    or:\n");
	fprintf(pFileTopo, "! IC atomA atomB *atomC atomD distAC angleACB angleABCD angleBCD distCD\n");
  fprintf(pFileTopo, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

  Atom* pIniAtomA = AtomArrayGetByName(&atoms,iniAtom1);
  Atom* pIniAtomB = AtomArrayGetByName(&atoms,iniAtom2);
  Atom* pIniAtomC = AtomArrayGetByName(&atoms,iniAtom3);
  if(pIniAtomA == NULL || pIniAtomB == NULL || pIniAtomC == NULL){
    pIniAtomA = AtomArrayGet(&atoms, 0);
    pIniAtomB = AtomArrayGet(&atoms, 1);
    pIniAtomC = AtomArrayGet(&atoms, 2);
  }


  IntArray icFlags;
  IntArrayCreate(&icFlags, AtomArrayGetCount(&atoms));
  for(int i=0; i<IntArrayGetLength(&icFlags); i++){
    Atom* pAtom = AtomArrayGet(&atoms,i);
    if(pAtom==pIniAtomA || pAtom==pIniAtomB || pAtom==pIniAtomC) IntArraySet(&icFlags,i,1);
    else IntArraySet(&icFlags,i,0);
  }

  BOOL allICsCalculated = FALSE;
  while(!allICsCalculated){
    allICsCalculated = TRUE;
    for(int i=0; i<AtomArrayGetCount(&atoms); i++){
      if(IntArrayGet(&icFlags, i) == 1) continue;
      BOOL icFoundFlag = FALSE;
      Atom *pAtomD = AtomArrayGet(&atoms, i);
      AtomArray atomCArray;
      AtomArrayCreate(&atomCArray);
      // find atomCs bonded with atomD, where atomC's IC should have been calculated
      FindAtomC(pAtomD,&atoms,&icFlags,&bondsFromToAndType,&atomCArray);
      for(int j=0; j<AtomArrayGetCount(&atomCArray); j++){
        Atom *pAtomC = AtomArrayGet(&atomCArray, j);
        AtomArray atomBArray;
        AtomArrayCreate(&atomBArray);
        FindAtomB(pAtomC,pAtomD,&atoms,&icFlags,&bondsFromToAndType,&atomBArray);
        if(AtomArrayGetCount(&atomBArray) >= 2){ // improper IC
          Atom* pAtomA = AtomArrayGet(&atomBArray,0);
          Atom* pAtomB = AtomArrayGet(&atomBArray,1);
          XYZ AB = XYZDifference(&pAtomA->xyz,&pAtomB->xyz);
          XYZ BC = XYZDifference(&pAtomB->xyz,&pAtomC->xyz);
          XYZ CD = XYZDifference(&pAtomC->xyz,&pAtomD->xyz);
          XYZ AC = XYZDifference(&pAtomA->xyz,&pAtomC->xyz);
          double dAC=XYZDistance(&pAtomA->xyz,&pAtomC->xyz);
          double aACB=XYZAngle(&AC,&BC);
          double xABCD=GetTorsionAngle(&pAtomA->xyz,&pAtomB->xyz,&pAtomC->xyz,&pAtomD->xyz);
          double aBCD=PI-XYZAngle(&BC,&CD);
          double dCD=XYZDistance(&pAtomC->xyz,&pAtomD->xyz);
          fprintf(pFileTopo, "IC   %-7.7s %-7.7s *%-7.7s %-7.7s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
            AtomGetName(pAtomA),AtomGetName(pAtomB),AtomGetName(pAtomC),AtomGetName(pAtomD),dAC,RadToDeg(aACB),RadToDeg(xABCD),RadToDeg(aBCD),dCD);
          IntArraySet(&icFlags, i, 1);
          icFoundFlag = TRUE;
        }
        else{
          for(int k = 0; k < AtomArrayGetCount(&atomBArray); k++){
            Atom* pAtomB = AtomArrayGet(&atomBArray,k);
            AtomArray atomAArray;
            AtomArrayCreate(&atomAArray);
            // the method for finding aomA is identical to that of finding atomB
            FindAtomB(pAtomB,pAtomC,&atoms,&icFlags,&bondsFromToAndType,&atomAArray);
            if(AtomArrayGetCount(&atomAArray) == 0){
              sprintf(errMsg,"in file %s line %d, cannot find proper IC for atom %s",__FILE__,__LINE__,AtomGetName(pAtomD));
              result=DataNotExistError;
              TraceError(errMsg,result);
              return result;
            }
            else{
              Atom* pAtomA = AtomArrayGet(&atomAArray,0);
              XYZ AB=XYZDifference(&pAtomA->xyz,&pAtomB->xyz);
              XYZ BC=XYZDifference(&pAtomB->xyz,&pAtomC->xyz);
              XYZ CD = XYZDifference(&pAtomC->xyz,&pAtomD->xyz);
              XYZ AC = XYZDifference(&pAtomA->xyz,&pAtomC->xyz);
              double dAB=XYZDistance(&pAtomA->xyz,&pAtomB->xyz);
              double aABC=PI-XYZAngle(&AB,&BC);
              double xABCD=GetTorsionAngle(&pAtomA->xyz,&pAtomB->xyz,&pAtomC->xyz,&pAtomD->xyz);
              double aBCD=PI-XYZAngle(&BC,&CD);
              double dCD=XYZDistance(&pAtomC->xyz,&pAtomD->xyz);
              fprintf(pFileTopo, "IC   %-7.7s %-7.7s  %-7.7s %-7.7s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
                AtomGetName(pAtomA),AtomGetName(pAtomB),AtomGetName(pAtomC),AtomGetName(pAtomD),dAB,RadToDeg(aABC),RadToDeg(xABCD),RadToDeg(aBCD),dCD);
              IntArraySet(&icFlags, i, 1);
              icFoundFlag = TRUE;
            }
            AtomArrayDestroy(&atomAArray);
            if(icFoundFlag == TRUE) break;
          }
        }
        AtomArrayDestroy(&atomBArray);
        if(icFoundFlag == TRUE) break;
      }
      AtomArrayDestroy(&atomCArray);
      if(icFoundFlag == TRUE) break;
    }

    for(int i=0; i<IntArrayGetLength(&icFlags); i++){
      if(IntArrayGet(&icFlags, i) == 0){
        allICsCalculated = FALSE;
        break;
      }
    }
  }
  fclose(pFileTopo);
  IntArrayDestroy(&icFlags);
	IntArrayDestroy(&bondsFromToAndType);
  AtomArrayDestroy(&atoms);

	return Success;
}


int GetAttachedAtomTypeNum(AtomArray* pAtoms,int bonds[][2],int bondNum,int atomIndex,int attachAtomTypeNum[10],int *attachCount,int attachIndex[10]){
  //attach[ 0]:C  attach[ 1]:H   attach[ 2]:O   attach[ 3]:N  attach[ 4]:P  attach[ 5]:S
  //attach[ 6];F  attach[ 7]:Cl  attach[ 8]:Br  attach[ 9]:I  attach[10]:B  attach[11]:Fe/Zn/X,metal
  for(int i=0;i<bondNum;i++){
    if(atomIndex==bonds[i][0]-1){
      attachIndex[*attachCount]=bonds[i][1]-1;
      (*attachCount)++;
      Atom* pAtom = AtomArrayGet(pAtoms,bonds[i][1]-1);
      if(pAtom->name[0]=='C') attachAtomTypeNum[0]++;
      else if(pAtom->name[0]=='H')attachAtomTypeNum[1]++;
      else if(pAtom->name[0]=='O')attachAtomTypeNum[2]++;
      else if(pAtom->name[0]=='N')attachAtomTypeNum[3]++;
      else attachAtomTypeNum[4]++;
    }
    else if(atomIndex==bonds[i][1]-1){
      attachIndex[*attachCount]=bonds[i][0]-1;
      (*attachCount)++;
      Atom* pAtom = AtomArrayGet(pAtoms,bonds[i][0]-1);
      if(pAtom->name[0]=='C') attachAtomTypeNum[0]++;
      else if(pAtom->name[0]=='H')attachAtomTypeNum[1]++;
      else if(pAtom->name[0]=='O')attachAtomTypeNum[2]++;
      else if(pAtom->name[0]=='N')attachAtomTypeNum[3]++;
      else attachAtomTypeNum[4]++;
    }
  }

  return Success;
}


int GenerateSmallMolParameterFromMol2(char* mol2file,char* parfile){
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  BOOL readingAtom = FALSE;
  BOOL readingBond = FALSE;
  char resiName[MAX_LENGTH_RESIDUE_NAME+1]="";
  int bondNum=0;
  int bonds[500][2];

  FILE* pMol2=fopen(mol2file,"r");

  AtomArray atoms;
  AtomArrayCreate(&atoms);
  Atom atom;
  AtomCreate(&atom);
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pMol2)!=NULL){
    char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1]="";
    sscanf(line,"%s",keyword);
    if(strcmp(keyword,"@<TRIPOS>ATOM")==0){
      readingAtom = TRUE;
      readingBond = FALSE;
      continue;
    }
    else if(strcmp(keyword,"@<TRIPOS>BOND")==0){
      readingAtom = FALSE;
      readingBond = TRUE;
      continue;
    }
    else if(keyword[0] == '@'){
      readingAtom = readingBond = FALSE;
      continue;
    }
    if(readingAtom){
      int atomId;
      char mol2Resi[MAX_LENGTH_RESIDUE_NAME+1];
      sscanf(line,"%d %s %lf %lf %lf %s %d %s %lf",&atomId,atom.name,&atom.xyz.X,&atom.xyz.Y,&atom.xyz.Z,atom.type,&atom.posInChain,mol2Resi,&atom.charge);
      if(strcmp(resiName,"")==0){
        strcpy(resiName,mol2Resi);
        if(strlen(resiName)>3) resiName[3]='\0';
        if(LigandResidueNameConflictWithAminoAcid(resiName)){
          strcpy(resiName,"LIG");
        }
      }
      AtomArrayAppend(&atoms,&atom);
    }
    if(readingBond){
      int id;
      char type[5];
      sscanf(line,"%d %d %d %s",&id,&bonds[bondNum][0],&bonds[bondNum][1],type);
      bondNum++;
    }
  }
  AtomDestroy(&atom);

  //create atom parameter for ligand
  for(int i=0; i<AtomArrayGetCount(&atoms); i++){
    int attachAtomTypeNum[10];
    int attachCount=0;
    int attachIndex[10];
    for(int j=0;j<10;j++) attachAtomTypeNum[j]=0;
    for(int j=0;j<10;j++) attachIndex[j]=-1;
    Atom* pAtom=AtomArrayGet(&atoms,i);
    pAtom->isBBAtom=FALSE;
    strcpy(pAtom->hbB2,"-");
    strcpy(pAtom->hbDorB,"-");
    strcpy(pAtom->hbHorA,"-");
    pAtom->hybridType=Type_AtomHybridType_None;
    pAtom->polarity=Type_AtomPolarity_NPAliphatic;

    if(pAtom->name[0]=='C' && strstr(pAtom->type,"C.")!=NULL){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      if(attachAtomTypeNum[1]>=3) pAtom->EEF1_atType=EEF1_AtomType_CH3E;
      else if(attachAtomTypeNum[1]==2) pAtom->EEF1_atType=EEF1_AtomType_CH2E;
      else if(attachAtomTypeNum[1]==1){
        pAtom->EEF1_atType=EEF1_AtomType_CH1E;
        if(strcmp(pAtom->type,"C.ar")==0){
          pAtom->EEF1_atType=EEF1_AtomType_CR1E;
        }
      }
      else{
        if(strcmp(pAtom->type,"C.ar")==0) pAtom->EEF1_atType=EEF1_AtomType_CR;
        else if(strcmp(pAtom->type,"C.cat")==0) pAtom->EEF1_atType=EEF1_AtomType_Carg;
        else if(strcmp(pAtom->type,"C.2")==0){
          pAtom->EEF1_atType=EEF1_AtomType_CR;
          if(attachAtomTypeNum[2]>0){
            pAtom->EEF1_atType=EEF1_AtomType_CO;
          }
          for(int k=0;k<attachCount;k++){
            Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
            if(strcmp(pAttachAtom->type,"O.co2")==0 || pAttachAtom->EEF1_atType==EEF1_AtomType_OOC){
              pAtom->EEF1_atType=EEF1_AtomType_COO;
              break;
            }
          }
        }
        else if(strcmp(pAtom->type,"C.3")==0){
          if(attachCount==4) pAtom->EEF1_atType=EEF1_AtomType_CH1E;
          else if(attachCount==3) pAtom->EEF1_atType=EEF1_AtomType_CH1E;
          else if(attachCount==2) pAtom->EEF1_atType=EEF1_AtomType_CH2E;
          else if(attachCount==1) pAtom->EEF1_atType=EEF1_AtomType_CH3E;
        }
        else pAtom->EEF1_atType=EEF1_AtomType_CR;
      }
    }
    else if(pAtom->name[0]=='H'){
      pAtom->EEF1_atType=EEF1_AtomType_Hnpl;
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      if(attachAtomTypeNum[2]>0 || attachAtomTypeNum[3]>0){
        pAtom->EEF1_atType=EEF1_AtomType_Hpol;
        pAtom->polarity=Type_AtomPolarity_P;
        strcpy(pAtom->hbHorA,"H");
        strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
      }
    }
    else if(pAtom->name[0]=='O'){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      pAtom->polarity=Type_AtomPolarity_P;
      if(strcmp(pAtom->type,"O.co2")==0){
        pAtom->EEF1_atType=EEF1_AtomType_OOC;
        pAtom->polarity=Type_AtomPolarity_C;
        strcpy(pAtom->hbHorA,"A");
        pAtom->hybridType=Type_AtomHybridType_SP2;
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(pAttachAtom->name[0]!='H'){
            strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
            break;
          }
        }
      }//COO-
      else if(strcmp(pAtom->type,"O.2")==0){
        pAtom->EEF1_atType=EEF1_AtomType_OC;
        strcpy(pAtom->hbHorA,"A");
        pAtom->hybridType=Type_AtomHybridType_SP2;
        if(attachAtomTypeNum[1]==1){
          pAtom->EEF1_atType=EEF1_AtomType_OH1;
        }
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(pAttachAtom->name[0]!='H'){
            strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
            break;
          }
        }
      }
      else if(strcmp(pAtom->type,"O.3")==0){
        pAtom->EEF1_atType=EEF1_AtomType_Oest;
        if(attachAtomTypeNum[1]==1){
          pAtom->EEF1_atType=EEF1_AtomType_OH1;
          strcpy(pAtom->hbHorA,"A");
          pAtom->hybridType=Type_AtomHybridType_SP3;
          for(int k=0;k<attachCount;k++){
            Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
            if(pAttachAtom->name[0]!='H'){
              strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
              break;
            }
          }
        }
        else if(attachAtomTypeNum[1]==2){
          pAtom->EEF1_atType=EEF1_AtomType_OH2;
          strcpy(pAtom->hbHorA,"A");
          pAtom->hybridType=Type_AtomHybridType_SP3;
          for(int k=0;k<attachCount;k++){
            Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
            if(pAttachAtom->name[0]!='H'){
              strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
              break;
            }
          }
        }
      }
      else{
        pAtom->EEF1_atType=EEF1_AtomType_Oest;
      }
    }
    else if(pAtom->name[0]=='N' && strstr(pAtom->type,"N.")!=NULL){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      pAtom->polarity=Type_AtomPolarity_P;
      if(strcmp(pAtom->type,"N.4")==0) pAtom->polarity=Type_AtomPolarity_C;
      if(attachAtomTypeNum[1]==3){
        pAtom->EEF1_atType=EEF1_AtomType_NH3;
        pAtom->hybridType=Type_AtomHybridType_SP3;
      }
      else if(attachAtomTypeNum[1]==2){
        pAtom->EEF1_atType=EEF1_AtomType_NH2;
        pAtom->hybridType=Type_AtomHybridType_SP3;
        if(strcmp(pAtom->type,"N.2")==0) pAtom->hybridType=Type_AtomHybridType_SP2;
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(strcmp(pAttachAtom->type,"C.cat")==0 || pAttachAtom->EEF1_atType==EEF1_AtomType_Carg){
            pAtom->EEF1_atType=EEF1_AtomType_Narg;
            break;
          }
        }
      }
      else if(attachAtomTypeNum[1]==1){
        pAtom->EEF1_atType=EEF1_AtomType_NH1;
        pAtom->hybridType=Type_AtomHybridType_SP3;
        if(strcmp(pAtom->type,"N.2")==0) pAtom->hybridType=Type_AtomHybridType_SP2;
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(strcmp(pAttachAtom->type,"C.cat")==0 || pAttachAtom->EEF1_atType==EEF1_AtomType_Carg){
            pAtom->EEF1_atType=EEF1_AtomType_Narg;
            break;
          }
        }
      }
      else{
        if(strcmp(pAtom->type,"N.ar")==0 || strcmp(pAtom->type,"N.2")==0){
          pAtom->EEF1_atType=EEF1_AtomType_NR;
          if(attachCount<3){
            strcpy(pAtom->hbHorA,"A");
            pAtom->hybridType=Type_AtomHybridType_SP2;
            for(int k=0;k<attachCount;k++){
              Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
              if(pAttachAtom->name[0]!='H'){
                strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
                break;
              }
            }
          }
          else pAtom->EEF1_atType=EEF1_AtomType_Npro;
        }
        else if(strcmp(pAtom->type,"N.am")==0) pAtom->EEF1_atType=EEF1_AtomType_Npro;
        else if(strcmp(pAtom->type,"N.3")==0 || strcmp(pAtom->type,"N.pl3")==0) pAtom->EEF1_atType=EEF1_AtomType_Npro;
        else if(attachCount==3) pAtom->EEF1_atType=EEF1_AtomType_Npro;
        else if(strcmp(pAtom->type,"N.1")==0 && attachCount==1) pAtom->EEF1_atType=EEF1_AtomType_NR;
        else pAtom->EEF1_atType=EEF1_AtomType_NR;
      }
    }
    else if(pAtom->name[0]=='P'){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      pAtom->EEF1_atType=EEF1_AtomType_P;
    }
    else if(pAtom->name[0]=='S'){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      if(attachAtomTypeNum[1]==1) pAtom->EEF1_atType=EEF1_AtomType_SH1E;
      else pAtom->EEF1_atType=EEF1_AtomType_S;
    }
    else if(pAtom->name[0]=='B' && strcmp(pAtom->type,"B")==0) pAtom->EEF1_atType=EEF1_AtomType_B;
    else if(pAtom->name[0]=='F' && strcmp(pAtom->type,"F")==0) pAtom->EEF1_atType=EEF1_AtomType_F;
    else if(pAtom->name[0]=='C' && strcmp(pAtom->type,"Cl")==0) pAtom->EEF1_atType=EEF1_AtomType_Cl;
    else if(pAtom->name[0]=='B' && strcmp(pAtom->type,"Br")==0) pAtom->EEF1_atType=EEF1_AtomType_Br;
    else if(pAtom->name[0]=='I' && strcmp(pAtom->type,"I")==0) pAtom->EEF1_atType=EEF1_AtomType_I;
    else if(pAtom->name[0]=='F' && strcmp(pAtom->type,"Fe")==0) pAtom->EEF1_atType=EEF1_AtomType_Fe;
    else if(pAtom->name[0]=='Z' && strcmp(pAtom->type,"F")==0) pAtom->EEF1_atType=EEF1_AtomType_Zn;
    else if(pAtom->name[0]=='A' && strcmp(pAtom->type,"Al")==0) pAtom->EEF1_atType=EEF1_AtomType_Al;
    else if(pAtom->name[0]=='N' && strcmp(pAtom->type,"Na")==0) pAtom->EEF1_atType=EEF1_AtomType_Na;
    else if(pAtom->name[0]=='K' && strcmp(pAtom->type,"K")==0) pAtom->EEF1_atType=EEF1_AtomType_K;
    else if(pAtom->name[0]=='C' && strcmp(pAtom->type,"Ca")==0) pAtom->EEF1_atType=EEF1_AtomType_Ca;
    else if(pAtom->name[0]=='M' && strcmp(pAtom->type,"Mg")==0) pAtom->EEF1_atType=EEF1_AtomType_Mg;
    else pAtom->EEF1_atType=EEF1_AtomType_Other;

    AssignAtomParameterByEEF1Type(pAtom,pAtom->EEF1_atType);
  }

  //output the parameter file
  FILE* pFilePara=fopen(parfile,"w");
  fprintf(pFilePara,"!%s\n",resiName);
  fprintf(pFilePara,"!Resi   Name   Type   Backbone   Polar   epsilon   Rmin     charge  HbH/A  HbD/B  HbB2  Hybrid  DG_free  Volume  Lambda\n");
  for(int i=0; i<AtomArrayGetCount(&atoms); i++){
    Atom* pAtom=AtomArrayGet(&atoms,i);
    char polarity[10];
    char hybrid[5];
    if(pAtom->polarity==Type_AtomPolarity_P) strcpy(polarity,"P");
    else if(pAtom->polarity==Type_AtomPolarity_C) strcpy(polarity,"C");
    else if(pAtom->polarity==Type_AtomPolarity_NPAliphatic) strcpy(polarity,"NP1");
    else strcpy(polarity,"NP2");
    if(pAtom->hybridType==Type_AtomHybridType_SP) strcpy(hybrid,"SP1");
    else if(pAtom->hybridType==Type_AtomHybridType_SP2) strcpy(hybrid,"SP2");
    else if(pAtom->hybridType==Type_AtomHybridType_SP3) strcpy(hybrid,"SP3");
    else strcpy(hybrid,"-");
    fprintf(pFilePara,"%-6s  %-6s %-6s %-10s %-7s %-8.4f %7.4f  %7.4f    %-5s  %-5s  %-5s %-4s   %7.2f     %4.1f    %3.1f\n",
      resiName,pAtom->name,pAtom->type,"N",polarity,
      pAtom->vdw_epsilon,pAtom->vdw_radius,pAtom->charge,
      pAtom->hbHorA,pAtom->hbDorB,pAtom->hbB2,hybrid,
      pAtom->EEF1_freeDG,pAtom->EEF1_volume,pAtom->EEF1_lamda_);
  }
  fclose(pFilePara);
  AtomArrayDestroy(&atoms);

  return Success;
}



int GenerateParameterAndTopologyFromMol2(char* mol2file,char* parfile,char* topfile){
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  BOOL readingAtom = FALSE;
  BOOL readingBond = FALSE;
  char resiName[MAX_LENGTH_RESIDUE_NAME+1]="";
  int bondNum=0;
  int bonds[500][2];

  AtomArray atoms;
  Atom atom;
  AtomArrayCreate(&atoms);
  AtomCreate(&atom);
  FILE* pMol2=fopen(mol2file,"r");
  while(fgets(line,MAX_LENGTH_ONE_LINE_IN_FILE,pMol2)!=NULL){
    char keyword[MAX_LENGTH_ONE_LINE_IN_FILE+1]="";
    sscanf(line,"%s",keyword);
    if(strcmp(keyword,"@<TRIPOS>ATOM")==0){
      readingAtom = TRUE;
      readingBond = FALSE;
      continue;
    }
    else if(strcmp(keyword,"@<TRIPOS>BOND")==0){
      readingAtom = FALSE;
      readingBond = TRUE;
      continue;
    }
    else if(keyword[0] == '@'){
      readingAtom = readingBond = FALSE;
      continue;
    }
    if(readingAtom){
      int atomId;
      char mol2Resi[MAX_LENGTH_RESIDUE_NAME+1];
      sscanf(line,"%d %s %lf %lf %lf %s %d %s %lf",&atomId,atom.name,&atom.xyz.X,&atom.xyz.Y,&atom.xyz.Z,atom.type,&atom.posInChain,mol2Resi,&atom.charge);
      if(strcmp(resiName,"")==0){
        strcpy(resiName,mol2Resi);
        if(strlen(resiName)>3) resiName[3]='\0';
        if(LigandResidueNameConflictWithAminoAcid(resiName)){
          strcpy(resiName,"LIG");
        }
      }
      AtomArrayAppend(&atoms,&atom);
    }
    if(readingBond){
      int id;
      char type[5];
      sscanf(line,"%d %d %d %s",&id,&bonds[bondNum][0],&bonds[bondNum][1],type);
      bondNum++;
    }
  }
  AtomDestroy(&atom);

  //create atom parameter for ligand
  for(int i=0; i<AtomArrayGetCount(&atoms); i++){
    int attachAtomTypeNum[10];
    int attachCount=0;
    int attachIndex[10];
    for(int j=0;j<10;j++) attachAtomTypeNum[j]=0;
    for(int j=0;j<10;j++) attachIndex[j]=-1;
    Atom* pAtom=AtomArrayGet(&atoms,i);
    pAtom->isBBAtom=FALSE;
    strcpy(pAtom->hbB2,"-");
    strcpy(pAtom->hbDorB,"-");
    strcpy(pAtom->hbHorA,"-");
    pAtom->hybridType=Type_AtomHybridType_None;
    pAtom->polarity=Type_AtomPolarity_NPAliphatic;

    if(pAtom->name[0]=='C' && strstr(pAtom->type,"C.")!=NULL){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      if(attachAtomTypeNum[1]>=3) pAtom->EEF1_atType=EEF1_AtomType_CH3E;
      else if(attachAtomTypeNum[1]==2) pAtom->EEF1_atType=EEF1_AtomType_CH2E;
      else if(attachAtomTypeNum[1]==1){
        pAtom->EEF1_atType=EEF1_AtomType_CH1E;
        if(strcmp(pAtom->type,"C.ar")==0){
          pAtom->EEF1_atType=EEF1_AtomType_CR1E;
        }
      }
      else{
        if(strcmp(pAtom->type,"C.ar")==0) pAtom->EEF1_atType=EEF1_AtomType_CR;
        else if(strcmp(pAtom->type,"C.cat")==0) pAtom->EEF1_atType=EEF1_AtomType_Carg;
        else if(strcmp(pAtom->type,"C.2")==0){
          pAtom->EEF1_atType=EEF1_AtomType_CR;
          if(attachAtomTypeNum[2]>0){
            pAtom->EEF1_atType=EEF1_AtomType_CO;
          }
          for(int k=0;k<attachCount;k++){
            Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
            if(strcmp(pAttachAtom->type,"O.co2")==0 || pAttachAtom->EEF1_atType==EEF1_AtomType_OOC){
              pAtom->EEF1_atType=EEF1_AtomType_COO;
              break;
            }
          }
        }
        else if(strcmp(pAtom->type,"C.3")==0){
          if(attachCount==4) pAtom->EEF1_atType=EEF1_AtomType_CH1E;
          else if(attachCount==3) pAtom->EEF1_atType=EEF1_AtomType_CH1E;
          else if(attachCount==2) pAtom->EEF1_atType=EEF1_AtomType_CH2E;
          else if(attachCount==1) pAtom->EEF1_atType=EEF1_AtomType_CH3E;
        }
        else pAtom->EEF1_atType=EEF1_AtomType_CR;
      }
    }
    else if(pAtom->name[0]=='H'){
      pAtom->EEF1_atType=EEF1_AtomType_Hnpl;
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      if(attachAtomTypeNum[2]>0 || attachAtomTypeNum[3]>0){
        pAtom->EEF1_atType=EEF1_AtomType_Hpol;
        pAtom->polarity=Type_AtomPolarity_P;
        strcpy(pAtom->hbHorA,"H");
        strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
      }
    }
    else if(pAtom->name[0]=='O'){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      pAtom->polarity=Type_AtomPolarity_P;
      if(strcmp(pAtom->type,"O.co2")==0){
        pAtom->EEF1_atType=EEF1_AtomType_OOC;
        pAtom->polarity=Type_AtomPolarity_C;
        strcpy(pAtom->hbHorA,"A");
        pAtom->hybridType=Type_AtomHybridType_SP2;
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(pAttachAtom->name[0]!='H'){
            strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
            break;
          }
        }
      }//COO-
      else if(strcmp(pAtom->type,"O.2")==0){
        pAtom->EEF1_atType=EEF1_AtomType_OC;
        strcpy(pAtom->hbHorA,"A");
        pAtom->hybridType=Type_AtomHybridType_SP2;
        if(attachAtomTypeNum[1]==1){
          pAtom->EEF1_atType=EEF1_AtomType_OH1;
        }
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(pAttachAtom->name[0]!='H'){
            strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
            break;
          }
        }
      }
      else if(strcmp(pAtom->type,"O.3")==0){
        pAtom->EEF1_atType=EEF1_AtomType_Oest;
        if(attachAtomTypeNum[1]==1){
          pAtom->EEF1_atType=EEF1_AtomType_OH1;
          strcpy(pAtom->hbHorA,"A");
          pAtom->hybridType=Type_AtomHybridType_SP3;
          for(int k=0;k<attachCount;k++){
            Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
            if(pAttachAtom->name[0]!='H'){
              strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
              break;
            }
          }
        }
        else if(attachAtomTypeNum[1]==2){
          pAtom->EEF1_atType=EEF1_AtomType_OH2;
          strcpy(pAtom->hbHorA,"A");
          pAtom->hybridType=Type_AtomHybridType_SP3;
          for(int k=0;k<attachCount;k++){
            Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
            if(pAttachAtom->name[0]!='H'){
              strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
              break;
            }
          }
        }
      }
      else{
        pAtom->EEF1_atType=EEF1_AtomType_Oest;
      }
    }
    else if(pAtom->name[0]=='N' && strstr(pAtom->type,"N.")!=NULL){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      pAtom->polarity=Type_AtomPolarity_P;
      if(strcmp(pAtom->type,"N.4")==0) pAtom->polarity=Type_AtomPolarity_C;
      if(attachAtomTypeNum[1]==3){
        pAtom->EEF1_atType=EEF1_AtomType_NH3;
        pAtom->hybridType=Type_AtomHybridType_SP3;
      }
      else if(attachAtomTypeNum[1]==2){
        pAtom->EEF1_atType=EEF1_AtomType_NH2;
        pAtom->hybridType=Type_AtomHybridType_SP3;
        if(strcmp(pAtom->type,"N.2")==0) pAtom->hybridType=Type_AtomHybridType_SP2;
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(strcmp(pAttachAtom->type,"C.cat")==0 || pAttachAtom->EEF1_atType==EEF1_AtomType_Carg){
            pAtom->EEF1_atType=EEF1_AtomType_Narg;
            break;
          }
        }
      }
      else if(attachAtomTypeNum[1]==1){
        pAtom->EEF1_atType=EEF1_AtomType_NH1;
        pAtom->hybridType=Type_AtomHybridType_SP3;
        if(strcmp(pAtom->type,"N.2")==0) pAtom->hybridType=Type_AtomHybridType_SP2;
        for(int k=0;k<attachCount;k++){
          Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
          if(strcmp(pAttachAtom->type,"C.cat")==0 || pAttachAtom->EEF1_atType==EEF1_AtomType_Carg){
            pAtom->EEF1_atType=EEF1_AtomType_Narg;
            break;
          }
        }
      }
      else{
        if(strcmp(pAtom->type,"N.ar")==0 || strcmp(pAtom->type,"N.2")==0){
          pAtom->EEF1_atType=EEF1_AtomType_NR;
          if(attachCount<3){
            strcpy(pAtom->hbHorA,"A");
            pAtom->hybridType=Type_AtomHybridType_SP2;
            for(int k=0;k<attachCount;k++){
              Atom* pAttachAtom=AtomArrayGet(&atoms,attachIndex[k]);
              if(pAttachAtom->name[0]!='H'){
                strcpy(pAtom->hbDorB,AtomGetName(AtomArrayGet(&atoms,attachIndex[0])));
                break;
              }
            }
          }
          else pAtom->EEF1_atType=EEF1_AtomType_Npro;
        }
        else if(strcmp(pAtom->type,"N.am")==0) pAtom->EEF1_atType=EEF1_AtomType_Npro;
        else if(strcmp(pAtom->type,"N.3")==0 || strcmp(pAtom->type,"N.pl3")==0) pAtom->EEF1_atType=EEF1_AtomType_Npro;
        else if(attachCount==3) pAtom->EEF1_atType=EEF1_AtomType_Npro;
        else if(strcmp(pAtom->type,"N.1")==0 && attachCount==1) pAtom->EEF1_atType=EEF1_AtomType_NR;
        else pAtom->EEF1_atType=EEF1_AtomType_NR;
      }
    }
    else if(pAtom->name[0]=='P'){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      pAtom->EEF1_atType=EEF1_AtomType_P;
    }
    else if(pAtom->name[0]=='S'){
      GetAttachedAtomTypeNum(&atoms,bonds,bondNum,i,attachAtomTypeNum,&attachCount,attachIndex);
      if(attachAtomTypeNum[1]==1) pAtom->EEF1_atType=EEF1_AtomType_SH1E;
      else pAtom->EEF1_atType=EEF1_AtomType_S;
    }
    else if(pAtom->name[0]=='B' && strcmp(pAtom->type,"B")==0) pAtom->EEF1_atType=EEF1_AtomType_B;
    else if(pAtom->name[0]=='F' && strcmp(pAtom->type,"F")==0) pAtom->EEF1_atType=EEF1_AtomType_F;
    else if(pAtom->name[0]=='C' && strcmp(pAtom->type,"Cl")==0) pAtom->EEF1_atType=EEF1_AtomType_Cl;
    else if(pAtom->name[0]=='B' && strcmp(pAtom->type,"Br")==0) pAtom->EEF1_atType=EEF1_AtomType_Br;
    else if(pAtom->name[0]=='I' && strcmp(pAtom->type,"I")==0) pAtom->EEF1_atType=EEF1_AtomType_I;
    else if(pAtom->name[0]=='F' && strcmp(pAtom->type,"Fe")==0) pAtom->EEF1_atType=EEF1_AtomType_Fe;
    else if(pAtom->name[0]=='Z' && strcmp(pAtom->type,"F")==0) pAtom->EEF1_atType=EEF1_AtomType_Zn;
    else if(pAtom->name[0]=='A' && strcmp(pAtom->type,"Al")==0) pAtom->EEF1_atType=EEF1_AtomType_Al;
    else if(pAtom->name[0]=='N' && strcmp(pAtom->type,"Na")==0) pAtom->EEF1_atType=EEF1_AtomType_Na;
    else if(pAtom->name[0]=='K' && strcmp(pAtom->type,"K")==0) pAtom->EEF1_atType=EEF1_AtomType_K;
    else if(pAtom->name[0]=='C' && strcmp(pAtom->type,"Ca")==0) pAtom->EEF1_atType=EEF1_AtomType_Ca;
    else if(pAtom->name[0]=='M' && strcmp(pAtom->type,"Mg")==0) pAtom->EEF1_atType=EEF1_AtomType_Mg;
    else pAtom->EEF1_atType=EEF1_AtomType_Other;

    AssignAtomParameterByEEF1Type(pAtom,pAtom->EEF1_atType);
  }

  //output the parameter file
  FILE* pFilePara=fopen(parfile,"w");
  fprintf(pFilePara,"!%s\n",resiName);
  fprintf(pFilePara,"!Resi   Name   Type   Backbone   Polar   epsilon   Rmin     charge  HbH/A  HbD/B  HbB2  Hybrid  DG_free  Volume  Lambda\n");
  for(int i=0; i<AtomArrayGetCount(&atoms); i++){
    Atom* pAtom=AtomArrayGet(&atoms,i);
    char polarity[10];
    char hybrid[5];
    if(pAtom->polarity==Type_AtomPolarity_P) strcpy(polarity,"P");
    else if(pAtom->polarity==Type_AtomPolarity_C) strcpy(polarity,"C");
    else if(pAtom->polarity==Type_AtomPolarity_NPAliphatic) strcpy(polarity,"NP1");
    else strcpy(polarity,"NP2");
    if(pAtom->hybridType==Type_AtomHybridType_SP) strcpy(hybrid,"SP1");
    else if(pAtom->hybridType==Type_AtomHybridType_SP2) strcpy(hybrid,"SP2");
    else if(pAtom->hybridType==Type_AtomHybridType_SP3) strcpy(hybrid,"SP3");
    else strcpy(hybrid,"-");
    fprintf(pFilePara,"%-6s  %-6s %-6s %-10s %-7s %-8.4f %7.4f  %5.2f    %-5s  %-5s  %-5s %-4s   %7.2f     %4.1f    %3.1f\n",
      resiName,pAtom->name,pAtom->type,"N",polarity,
      pAtom->vdw_epsilon,pAtom->vdw_radius,pAtom->charge,
      pAtom->hbHorA,pAtom->hbDorB,pAtom->hbB2,hybrid,
      pAtom->EEF1_freeDG,pAtom->EEF1_volume,pAtom->EEF1_lamda_);
  }
  fclose(pFilePara);

  //output the topology file
  FILE* pFileTopo=fopen(topfile,"w");
  fprintf(pFileTopo,"RESI %s\n",resiName);
  for(int i=0;i<AtomArrayGetCount(&atoms);i++){
    fprintf(pFileTopo,"ATOM %s\n",AtomGetName(AtomArrayGet(&atoms,i)));
  }
  fprintf(pFileTopo,"\n");
  for(int i=0;i<bondNum;i++){
    int from=bonds[i][0];
    int to=bonds[i][1];
    fprintf(pFileTopo,"BOND %s %s\n",AtomGetName(AtomArrayGet(&atoms,from-1)),AtomGetName(AtomArrayGet(&atoms,to-1)));
  }
  fprintf(pFileTopo,"\n");
  //calculate charmmIC for hydrogen atoms
  for(int i=0;i<AtomArrayGetCount(&atoms);i++){
    Atom* pAtom=AtomArrayGet(&atoms,i);
    if(AtomIsHydrogen(pAtom)==FALSE) continue;
    //get atomC
    int atomcIndex=-1;
    for(int j=0;j<bondNum;j++){
      if(bonds[j][0]-1==i){atomcIndex=bonds[j][1]-1; break;}
      else if(bonds[j][1]-1==i){atomcIndex=bonds[j][0]-1; break;}
    }
    //get atomB's
    int atomBnum=0;
    int atomBindex[3]={-1,-1,-1};
    for(int j=0;j<bondNum;j++){
      if(bonds[j][0]-1==atomcIndex && bonds[j][1]-1!=i && AtomIsHydrogen(AtomArrayGet(&atoms,bonds[j][1]-1))==FALSE){atomBindex[atomBnum++]=bonds[j][1]-1;}
      else if(bonds[j][1]-1==atomcIndex && bonds[j][0]-1!=i && AtomIsHydrogen(AtomArrayGet(&atoms,bonds[j][0]-1))==FALSE){atomBindex[atomBnum++]=bonds[j][0]-1;}
    }
    //check if improper dihedral angle can be used
    if(atomBnum>=2){
      Atom* pAtomA=AtomArrayGet(&atoms,atomBindex[0]);
      Atom* pAtomB=AtomArrayGet(&atoms,atomBindex[1]);
      Atom* pAtomC=AtomArrayGet(&atoms,atomcIndex);
      Atom* pAtomD=AtomArrayGet(&atoms,i);
      XYZ AB=XYZDifference(&pAtomA->xyz,&pAtomB->xyz);
      XYZ BC=XYZDifference(&pAtomB->xyz,&pAtomC->xyz);
      XYZ CD = XYZDifference(&pAtomC->xyz,&pAtomD->xyz);
      XYZ AC = XYZDifference(&pAtomA->xyz,&pAtomC->xyz);

      double dAC=XYZDistance(&pAtomA->xyz,&pAtomC->xyz);
      double aACB=XYZAngle(&AC,&BC);
      double xABCD=GetTorsionAngle(&pAtomA->xyz,&pAtomB->xyz,&pAtomC->xyz,&pAtomD->xyz);
      double aBCD=PI-XYZAngle(&BC,&CD);
      double dCD=XYZDistance(&pAtomC->xyz,&pAtomD->xyz);
      fprintf(pFileTopo,"IC   %-7.7s %-7.7s *%-7.7s %-7.7s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
        AtomGetName(pAtomA),AtomGetName(pAtomB),AtomGetName(pAtomC),AtomGetName(pAtomD),dAC,RadToDeg(aACB),RadToDeg(xABCD),RadToDeg(aBCD),dCD);
    }
    else{
      //get atomA's
      int atomaIndex=-1;
      for(int j=0;j<bondNum;j++){
        if(bonds[j][0]-1==atomBindex[0] && bonds[j][1]-1!=atomcIndex){
          atomaIndex=bonds[j][1]-1;
          break;
        }
        else if(bonds[j][1]-1==atomBindex[0] && bonds[j][0]-1!=atomcIndex){
          atomaIndex=bonds[j][0]-1;
          break;
        }
      }
      Atom* pAtomA=AtomArrayGet(&atoms,atomaIndex);
      Atom* pAtomB=AtomArrayGet(&atoms,atomBindex[0]);
      Atom* pAtomC=AtomArrayGet(&atoms,atomcIndex);
      Atom* pAtomD=AtomArrayGet(&atoms,i);
      XYZ AB=XYZDifference(&pAtomA->xyz,&pAtomB->xyz);
      XYZ BC=XYZDifference(&pAtomB->xyz,&pAtomC->xyz);
      XYZ CD = XYZDifference(&pAtomC->xyz,&pAtomD->xyz);
      XYZ AC = XYZDifference(&pAtomA->xyz,&pAtomC->xyz);

      double dAB=XYZDistance(&pAtomA->xyz,&pAtomB->xyz);
      double aABC=PI-XYZAngle(&AB,&BC);
      double xABCD=GetTorsionAngle(&pAtomA->xyz,&pAtomB->xyz,&pAtomC->xyz,&pAtomD->xyz);
      double aBCD=PI-XYZAngle(&BC,&CD);
      double dCD=XYZDistance(&pAtomC->xyz,&pAtomD->xyz);
      fprintf(pFileTopo,"IC   %-7.7s %-7.7s  %-7.7s %-7.7s  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
        AtomGetName(pAtomA),AtomGetName(pAtomB),AtomGetName(pAtomC),AtomGetName(pAtomD),dAB,RadToDeg(aABC),RadToDeg(xABCD),RadToDeg(aBCD),dCD);
    }
  }
  fclose(pFileTopo);

  AtomArrayDestroy(&atoms);

  return Success;
}
