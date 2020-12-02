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

#ifndef SMALL_MOL_H
#define SMALL_MOL_H

#include "DesignSite.h"
#include "Rotamer.h"

#define DISTURBANCE_IN_RANGE_CHECK 1e-3 
#define MAX_COUNT_CHECK_MULTI_CONS 50
#define SMALLMOL_ATOM_SYMBOL       '$'
#define SECOND_CATACON_SITE_SYMBOL '$'
#define PSEUDO_ATOM_SYMBOL         '+'

typedef enum _Type_CataCons{
    Type_CataCons_NewPseudoAtom,
    Type_CataCons_Distance,
    Type_CataCons_Angle,
    Type_CataCons_Torsion
} Type_CataCons;


typedef enum _Type_CataConsAtomLocation{
    Type_CataConsAtomLocation_FirstSite,
    Type_CataConsAtomLocation_SecondSite,
    Type_CataConsAtomLocation_Pseudo,
    Type_CataConsAtomLocation_Undefined
} Type_CataConsAtomLocation;


typedef struct _CataConsItem{
	Type_CataCons type;
	char atomName[4][MAX_LENGTH_ATOM_NAME+1];
	Type_CataConsAtomLocation atomLoc[4];
	int  atomIndex[4];
	double min;
	double max;
} CataConsItem;


int CataConsItemShow(CataConsItem* pThis);

typedef struct _CataConsGroup{
	char site1RotType[MAX_LENGTH_RESIDUE_NAME+1];
	char site2RotType[MAX_LENGTH_RESIDUE_NAME+1];
	int nconsItems;
	CataConsItem* consItems;
	int npseudoAtoms;
	XYZ* pseudoAtoms;
} CataConsGroup;

int CataConsGroupDeploy(CataConsGroup* pThis,Rotamer* pFirstSiteRotamer,Rotamer* pSecondSiteRotamer);

typedef enum _Type_SitePairCons{
  Type_SitePairCons_Covalent,
  Type_SitePairCons_Hbond,
  Type_SitePairCons_SaltBridge,
  Type_SitePairCons_VDW,
  Type_SitePairCons_PiPiStack,
  Type_SitePairCons_TStack,
  Type_SitePairCons_Other
}Type_SitePairCons;

typedef struct _CataConsSitePair{
  char chName1[MAX_LENGTH_CHAIN_NAME+1];
  int  pos1;
  char resiName1[MAX_LENGTH_RESIDUE_NAME+1]; 
	char chName2[MAX_LENGTH_CHAIN_NAME+1];
	int  pos2;
	char resiName2[MAX_LENGTH_RESIDUE_NAME+1];
	int ngroups;
  CataConsGroup* groups;
  BOOL deployedFlag;
  Type_SitePairCons pairConsType;
} CataConsSitePair;


int CataConsSitePairCreate(CataConsSitePair* pThis,FileReader* file);
void CataConsSitePairDestroy(CataConsSitePair* pThis);
BOOL CataConsSitePairGetDeployedFlag(CataConsSitePair* pThis);
int CataConsSitePairDeploy(CataConsSitePair* pThis,RotamerSet* pFirstSiteRotSet,RotamerSet* pSecondSiteRotSet);
int CataConsSitePairShow(CataConsSitePair* pThis);


typedef struct _CataConsSitePairArray{
 int count;
 CataConsSitePair* sites;
} CataConsSitePairArray;


int CataConsSitePairArrayCreate(CataConsSitePairArray* pThis,char* cataConsFile);
void CataConsSitePairArrayDestroy(CataConsSitePairArray* pThis);
int CataConsSitePairArrayGetCount(CataConsSitePairArray* pThis);
CataConsSitePair* CataConsSitePairArrayGet(CataConsSitePairArray* pThis,int index);
CataConsSitePair* CataConsSitePairArrayFind(CataConsSitePairArray* pThis,char* chnName1,int pos1,char* resiName1,char* chnName2,int pos2,char* resiName2);

int CataConsSitePairArrayShow(CataConsSitePairArray* pThis);
int CataConsSitePairArrayTester(char* filename);

BOOL CataConsItemCheck(CataConsItem* pThis,XYZArray* pOnFirstSite, XYZArray* pOnSecondSite,XYZ* pseudoAtoms);
BOOL CataConsGroupCheck(CataConsGroup* pThis,Rotamer* pOnFirstSite,Rotamer* pOnSecondSite);
BOOL CataConsSitePairCheck(CataConsSitePair* pThis,Rotamer* pOnFirstSite,Rotamer* pOnSecondSite);

//Various Placing Actions
typedef enum _Type_PlacingAction{
  Type_PlacingAction_Load,
  Type_PlacingAction_Evaluate,
  Type_PlacingAction_Variate,
  Type_PlacingAction_Calc,
  Type_PlacingAction_CheckCataCons,
	Type_PlacingAction_CheckMultiCons,
  Type_PlacingAction_CheckVDW_Backbone,
	Type_PlacingAction_CheckVDW_Internal,
  Type_PlacingAction_CheckRMSD
} Type_PlacingAction;


typedef enum _Type_PlacingActionEvalParam{
  Type_PlacingActionEvaluate_Distance,
  Type_PlacingActionEvaluate_HalfDistance,
  Type_PlacingActionEvaluate_Angle,
  Type_PlacingActionEvaluate_Torsion
} Type_PlacingActionEvalParam;


typedef struct _CheckMultiConsStep{
	int designSite;
	int countOfConsToCheckAtThisStep;
	int consToCheckAtThisStep[MAX_COUNT_CHECK_MULTI_CONS];
} CheckMultiConsStep;


typedef struct _PlacingAction{
  Type_PlacingAction actionType;

  //Fields if actionType is "LOAD"
  IntArray load_atoms;

  //Fields if actionType is "VARIATE" 
  int variate_param;
  double variate_from;
  double variate_to;
  double variate_increment;

  //Fields if actionType is "CALC" 
  IntArray calc_atoms;
  IntArray calc_params;

  //Fields if actionType is "CHECK_CATA_CONS" 
  char checkCataCons_firstSiteChainName[MAX_LENGTH_CHAIN_NAME+1];
  char checkCataCons_firstSiteResidueName[MAX_LENGTH_RESIDUE_NAME+1];
  int  checkCataCons_firstSitePosInChain;
	char checkCataCons_secondSiteChainName[MAX_LENGTH_CHAIN_NAME+1];
	char checkCataCons_secondSiteResidueName[MAX_LENGTH_RESIDUE_NAME+1];
	int  checkCataCons_secondSitePosInChain;
	CataConsSitePair* checkCataCons_pCataCon;
	int checkCataCons_pSite[2];

	//Fields if actionType is "CHECK_MULTI_CONS"
	int checkMultiCons_cataConsCount;
	StringArray checkMultiCons_firstSiteChainNames;
	StringArray checkMultiCons_firstSiteResidueNames;
	IntArray    checkMultiCons_firstSitePosInChains;
	StringArray checkMultiCons_secondSiteChainNames;
	StringArray checkMultiCons_secondSiteResidueNames;
	IntArray    checkMultiCons_secondSitePosInChains;
	CataConsSitePair* checkMultiCons_pCataCons[MAX_COUNT_CHECK_MULTI_CONS];
	int checkMultiCons_siteIndexes[MAX_COUNT_CHECK_MULTI_CONS][2];
	int checkMultiCons_rotamerCountOnEachSite[MAX_COUNT_CHECK_MULTI_CONS][2];
	BOOL** checkMultiCons_predeterminedConsRelations[MAX_COUNT_CHECK_MULTI_CONS];
	int checkMultiCons_stepCount;
	CheckMultiConsStep checkMultiCons_steps[MAX_COUNT_CHECK_MULTI_CONS*2];
	
  //Fields if actionType is "CHECK_VDW_BACKBONE"
  double checkVDW_backbone_maxAllowed;
  BOOL checkVDW_backbone_withHydrogen;
  double checkVDW_backbone_activeRange;
  double checkVDW_backbone_totalPotential;
	IntArray checkVDW_backbone_smallmolAtomHasXyz;

	//Fields if actionType is "CHECK_VDW_INTERNAL"
	double checkVDW_internal_maxAllowed;
	BOOL checkVDW_internal_withHydrogen;
	double checkVDW_internal_totalPotential;
	int checkVDW_internal_smallMolAtomCount;
	BOOL** checkVDW_internal_smallMolAtom13bondedMatrix;
	IntArray checkVDW_internal_smallmolAtomHasXyz;

  //Fields if actionType is "CHECK_RMSD" 
	IntArray checkRMSD_smallmolAtomHasXyz;
  double checkRMSD_minDifference;
  BOOL checkRMSD_withHydrogen;

  //Fields if actionType is "EVALUATE" 
  int evaluate_param;
  Type_PlacingActionEvalParam evaluate_type;
  IntArray evaluate_atoms; 
} PlacingAction;


int PlacingActionCreate(PlacingAction* pThis,Type_PlacingAction type);
int PlacingActionDestroy(PlacingAction* pThis);
BOOL PlacingActionValidParamName(char* paramName);
int PlacingActionRead_CALC(PlacingAction* pThis,StringArray* pAtomNames,XYZArray* pAtomXYZs,StringArray* pParamNames,DoubleArray* pParams,StringArray* pContent);
//The following functions begin with "PlacingActionRead" are all sub-routines of "PlacingRuleReadFile()" 
int PlacingActionRead_CHECK_CATA_CONS(PlacingAction* pThis,StringArray* pContent);
int PlacingActionRead_CHECK_MULTI_CONS(PlacingAction* pThis,StringArray* pContent);
int PlacingActionRead_CHECK_RMSD(PlacingAction* pThis,StringArray* pContent);
int PlacingActionRead_CHECK_VDW_BACKBONE(PlacingAction* pThis,StringArray* pContent);
int PlacingActionRead_CHECK_VDW_INTERNAL(PlacingAction* pThis,StringArray* pContent);
int PlacingActionRead_EVALUATE(PlacingAction* pThis,StringArray* pAtomNames,StringArray* pParamNames,DoubleArray* pParams,StringArray* pContent);
int PlacingActionRead_LOAD(PlacingAction* pThis,StringArray* pAtomNames,XYZArray* pAtomXYZs,StringArray* pContent);
int PlacingActionRead_VARIATE(PlacingAction* pThis,StringArray* pParamNames,DoubleArray* pParams,StringArray* pContent);


typedef struct _PlacingRule{
  char chainName[MAX_LENGTH_CHAIN_NAME+1];
  char residueName[MAX_LENGTH_RESIDUE_NAME+1];
  int  posInChain;
  char rotamerType[MAX_LENGTH_RESIDUE_NAME+1];
  StringArray atomNames;
  XYZArray    atomXYZs;
  StringArray paramNames;
  DoubleArray params;
  PlacingAction* actions;
  int actionCount;

	// these two fields have no meaning before the rule is deployed;
  IntArray atomPosOnSmallMol;
  XYZArray smallMolAtomXYZs;
	Residue* pSmallMol;
  AtomArray* pTruncatedBackbone;

	// a flag for recording if the rule has been deployed.
  BOOL deployedFlag;

	double vdwBackbone;
	double vdwInternal;
} PlacingRule;

int PlacingRuleCreate(PlacingRule* pThis,char* fileName);
void PlacingRuleDestroy(PlacingRule* pThis);
BOOL PlacingRuleGetDeployedFlag(PlacingRule* pThis);
char* PlacingRuleGetResiName(PlacingRule* pThis);
int PlacingRuleGetPosInChain(PlacingRule* pThis);
char* PlacingRuleGetChainName(PlacingRule* pThis);
char* PlacingRuleGetRotamerType(PlacingRule* pThis);
double PlacingRuleGetTruncatedBackboneRange(PlacingRule* pThis);
//The following three functions are all huge
//PlacingRuleReadFile() will call PlacingActionRead_XXX to read every PlacingAction 
int	PlacingRuleReadFile(PlacingRule* pThis,FileReader* pFileReader);
//PlacingRuleDeploy() will call PlacingActionDeploy_XXX to deploy every PlacingAction
int PlacingRuleDeploy(PlacingRule* pThis,
							  Residue* pSmallMol,
							  CataConsSitePairArray* pCataConsArray,
							  int relatedProteinSiteCount,
							  DesignSite** relatedProteinSites,
							  AtomArray* pTruncatedBackbone);
//The following functions begin with "PlacingActionDeploy_XXX" are all sub-routines of "PlacingRuleDepoly()" 
int PlacingActionDeploy_CheckCataCons(PlacingAction* pThis,
									  Residue* pSmallMol,
									  CataConsSitePairArray* pCataConsArray,
									  char* startingSiteChainName,
									  int startingSitePosInChain,
									  int relatedProteinSiteCount,
									  DesignSite* relatedProteinSites);
int PlacingActionDeploy_CheckMultiCons(PlacingAction* pThis,
									   Residue* pSmallMol,
									   CataConsSitePairArray* pCataConsArray,
									   char* startingSiteChainName,
									   char* startingSiteResidueName,
									   int startingSitePosInChain,
									   int relatedProteinSiteCount,
									   DesignSite* relatedProteinSites);
int PlacingActionDeploy_CheckRMSD(PlacingAction* pThis,IntArray* pSmallmolAtomGetXyzByThisStep);
int PlacingActionDeploy_CheckVDWBackbone(PlacingAction* pThis,IntArray* pSmallmolAtomGetXyzByThisStep);
int PlacingActionDeploy_CheckVDWInternal(PlacingAction* pThis,Residue* pSmallMol,
										 IntArray* pSmallmolAtomGetXyzByThisStep);
int PlacingActionDeploy_Calc(PlacingAction* pThis,
							 IntArray* pAtomPosOnSmallMol,
							 IntArray* pSmallmolAtomGetXyzByThisStep);

int PlacingActionDeploy_Load(PlacingAction* pThis,
							 IntArray* pAtomPosOnSmallMol,
							 IntArray* pSmallmolAtomGetXyzByThisStep);
//PlacingRuleProcess() will call the following sub-routines to execute every PlacingAction 
int PlacingRuleProcess(	PlacingRule* pThis,
								Rotamer* pProteinRotamerOnStartingSite,
								CataConsSitePairArray* pCataConsArray,
								int relatedProteinSiteCount,
								DesignSite* relatedProteinSites,
								RotamerSet* pSmallMolRotSetForOutput,
								int step);
//The following functions begin with "PlacingRuleProcess_" are all sub-routines of "PlacingRuleProcess()" 
int PlacingRuleProcess_CALC(PlacingRule* pThis,PlacingAction* pAction);
BOOL PlacingRuleProcess_CHECK_CATA_CONS(PlacingRule* pThis,
										  Rotamer* pProteinRotamerOnStartingSite,										  
										  int relatedProteinSiteCount,
										  DesignSite* relatedProteinSites,
										  PlacingAction* pAction);
BOOL PlacingRuleProcess_CHECK_MULTI_CONS_Recurrence(PlacingAction* pAction,
													  int relatedProteinSiteCount,
													  DesignSite** relatedProteinSites,
													  int* rotamersSelectedOnPreviousSteps,
													  Rotamer* pSmallMolRotamer,
													  Rotamer* pStartingProteinRotamer,
													  int currentStepIndex);
BOOL PlacingRuleProcess_CHECK_MULTI_CONS( PlacingRule* pThis,
											Rotamer* pProteinRotamerOnStartingSite,										  
											int relatedProteinSiteCount,
											DesignSite** relatedProteinSites,
											PlacingAction* pAction);
BOOL PlacingRuleProcess_CHECK_RMSD(PlacingRule* pThis,RotamerSet* pSmallMolRotSet,PlacingAction* pAction);
BOOL PlacingRuleProcess_CHECK_VDW_BACKBONE(PlacingRule* pThis,PlacingAction* pAction);
BOOL PlacingRuleProcess_CHECK_VDW_INTERNAL(PlacingRule* pThis,PlacingAction* pAction);
int PlacingRuleProcess_EVALUATE(PlacingRule* pThis,PlacingAction* pAction);
int PlacingRuleProcess_LOAD(PlacingRule* pThis,
							   Rotamer* pProteinRotamerOnStartingSite,
							   PlacingAction* pAction);

//PlacingRulePlaceSmallMol() is the ultimate user interface in order to place small molecules 
int PlacingRulePlaceSmallMol(PlacingRule* pThis,
								CataConsSitePairArray* pCataConsArray,
								DesignSite* pStartingSite,
								int relatedProteinDesignSiteCount,
								DesignSite** relatedProteinDesignSites,
								RotamerSet* pSmallMolRotSetForOutput);

//Methods used for debugging 
int PlacingRuleShowAtom(PlacingRule* pThis,int index);
int PlacingRuleShowAtomDistances(PlacingRule* pThis);
int PlacingRuleShowParam(PlacingRule* pThis,int index);
int PlacingRuleShowAction(PlacingRule* pThis,int index);
int PlacingRuleShow(PlacingRule* pThis);
int PlacingRuleTester(char* filename);


// functions for screening and analyzing small-molecule rotamers
int ScreenSmallmolRotamersByRMSD(char* oriFileName,char* outputFileName,double rmsdThresold);
int AnalyzeSmallMolRotamers(char* oriFileName,Residue* pOriginalSmallMol);
int AnalyzeSmallMolRotamersForSpecifiedAtoms(char* oriFileName,Residue* pNativeSmallMol,char* specificAtomFile);
int CompareByInternalEnergy(const void *a, const void *b);
int CompareByBackboneEnergy(const void *a, const void *b);
int SmallMolRotamersGetBothHighRankOfBackboneVdwAndInternalVdw(char *oldRotamersFile,char* newRotamersFile,double percent);

#endif  //SMALL_MOL_H 
