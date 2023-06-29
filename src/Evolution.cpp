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


/******************************************************************************/
/*                                                                            */
/*                                xtd-seqali                                  */
/*                                                                            */
/* Computes the optimum alignment between two amino acid sequences (S1 & S2). */
/* The score of each pair (i,j) of positions takes into account not only      */
/* the similarity between the amino acids but also other features that        */
/* characterize those positions.                                              */
/*                                                                            */
/* In the current implementation, the score of each pair (i,j) of positions   */
/* is given by the secondary structures at positions i and j in S1 and S2,    */
/* respectively, and by the value of S1[i] at position j in the PSSM of S2.   */
/*                                                                            */
/* ARGUMENTS:                                                                 */
/* ---> S1:  path to the file containing amino acid sequence S1.              */
/*           The sequence is all on one line, starts at the start of the      */
/*           file, and is immediately followed by a line terminator.          */
/* ---> PF2: path to the file containing the PSSM of S2.                      */
/*           The PSSM can be output either by program makemat from the        */
/*           Impala suite, or by program mkprf (see notes).                   */
/* ---> SS1: path to the file containing the secondary structure of S1.       */
/*           The secondary structure is a horizontal sequence, with one       */
/*           character per residue position. It starts at the start of the    */
/*           file and is immediately followed by a line terminator.           */
/*           The i-th character represents the secondary structure at the     */
/*           i-th residue position (i=0,...,N-1, where N is the number of     */
/*           residues).            .                                          */
/* ---> SS2: path to the file containing the secondary structure of S2.       */
/*           The file format is the same as for SS1.                          */
/* ---> WSS: weight of the secondary structure in the alignment.              */
/*                                                                            */
/* NOTES:                                                                     */
/* - the PSSM output by program mkprf consists of a header line, N data lines,*/
/*   and a footer line.                                                       */
/*     The header line contains 21 columns, with one amino acid type per      */
/*   column, expressed in one-letter code. The last column, labeled '-'       */
/*   refers to the gap. The columns are separated by a certain number of      */
/*   blanks and the last column is immediately followd by a line terminator.  */
/*     The i-th data line contains the profile data for the i-th residue      */
/*   position, for i=0,...,N-1, where N is the number of residue positions.   */
/*   The j-th column is the score of the j-th amino acid type at that         */
/*   position, for j=0,...,19. The 20-th column is the score of the gap at    */
/*   that position. The 20-th column is immediately followed by a line        */
/*   terminator.                                                              */
/*     The footer line is identical to the header line. After its line        */
/*   terminator the file terminates.                                          */
/*                                                                            */
/******************************************************************************/
#pragma warning(disable:4996)
#pragma warning(disable:4305)
#pragma warning(disable:4244)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "EvoSeqAlign.h"
#include "EvoGetSeq2SA.h"
#include "EvoGetSeq2SS.h"
#include "EvoGetPhiPsi.h"
#include "EvoUtility.h"
#include "ErrorTracker.h"
#include "EvoAminoName.h"
#include "Utility.h"

float wss = 1.58;
float wsa = 2.45;
float wang = 1.00;

using namespace CharSeq;

extern float** PROT_PROFILE;
extern char PROGRAM_PATH[MAX_LEN_ONE_LINE_CONTENT + 1];
extern char TGT_PRF[MAX_LEN_FILE_NAME + 1];
extern char TGT_SA[MAX_LEN_FILE_NAME + 1];
extern char TGT_SS[MAX_LEN_FILE_NAME + 1];
extern char TGT_SEQ[MAX_LEN_FILE_NAME + 1];
extern char TGT_PHIPSI[MAX_LEN_FILE_NAME + 1];

namespace Text
{
  int ncharLine(FILE* fp);
}


float EvolutionScoreAllFromFile(char* seqfile)
{
  char evolutiondir[500] = "";
  sprintf(evolutiondir, "%s/evolution", PROGRAM_PATH);
  printf("evolution parameter file path is: %s\n", evolutiondir);

  FILE* fp = fopen(seqfile, "r");
  if (fp == NULL)
  {
    printf("in file %s line %d, cannot read file %s\n", __FILE__, __LINE__, seqfile);
    return IOError;
  }
  SequenceData dsInfo;
  dsInfo.len = Text::ncharLine(fp);
  dsInfo.seq = new char[dsInfo.len + 1];
  fseek(fp, 0, SEEK_SET);
  fscanf(fp, "%s", dsInfo.seq);
  fclose(fp);


  dsInfo.ss1 = new char[dsInfo.len + 1];
  SSPrediction::getSeq2SS SS(dsInfo.ss1, dsInfo.seq, evolutiondir);
  dsInfo.sa1 = new char[dsInfo.len + 1];
  SAPrediction::getSeq2SA SA(dsInfo.sa1, dsInfo.seq, dsInfo.ss1, evolutiondir);

  //phi-psi prediction
  PhiPsiPrediction::getPhiPsi(dsInfo, evolutiondir);

  dsInfo.ss2 = new char[dsInfo.len + 1];
  fp = fopen(TGT_SS, "r");
  fscanf(fp, "%s", dsInfo.ss2);
  fclose(fp);
  dsInfo.sa2 = new char[dsInfo.len + 1];
  fp = fopen(TGT_SA, "r");
  fscanf(fp, "%s", dsInfo.sa2);
  fclose(fp);


  wsa /= (dsInfo.len * 1.0);    // weight for solvent accessibility
  wss /= (dsInfo.len * 1.0);    // weight for secondary structure
  wang /= (dsInfo.len * 180.0);  // weight for phi-psi prediction

  CharSeq::SeqAlign sa(dsInfo, TGT_PRF, TGT_PHIPSI, wss, wsa, wang);

  float max = sa.printMtxInfo();
  printf("%f\n", max);

  delete[] dsInfo.ss1;  delete[] dsInfo.ss2;
  delete[] dsInfo.sa1;  delete[] dsInfo.sa2;
  delete[] dsInfo.seq;

  return max;
}


float EvolutionScorePrfFromFile(char* seqfile)
{
  char evolutiondir[500] = "";
  sprintf(evolutiondir, "%s/evolution", PROGRAM_PATH);
  //printf("evolution parameter file path is: %s\n", evolutiondir);

  FILE* fp = fopen(seqfile, "r");
  if (fp == NULL)
  {
    printf("in file %s line %d, cannot read file %s\n", __FILE__, __LINE__, seqfile);
    return IOError;
  }
  SequenceData dsInfo;
  dsInfo.len = Text::ncharLine(fp);
  dsInfo.seq = new char[dsInfo.len + 1];
  fseek(fp, 0, SEEK_SET);
  fscanf(fp, "%s", dsInfo.seq);
  fclose(fp);

  //read profile
  float** prf;
  fp = fopen(TGT_PRF, "r");
  int len = dsInfo.len;
  prf = new float* [len];
  for (int i = 0; i < len; ++i)
    prf[i] = new float[AMINOS];
  AminoName map[AMINOS];
  char ami[2];
  for (int i = 0; i < AMINOS; ++i)
  {
    fscanf(fp, "%s", ami);
    map[i] = charToAmino(ami[0]);
  }
  fscanf(fp, "%*s");
  for (int i = 0; i < len; ++i)
  {
    fscanf(fp, "%*s");
    for (int j = 0; j < AMINOS; ++j)
    {
      fscanf(fp, "%f", &prf[i][map[j]]);
    }
    fscanf(fp, "%*s");
  }
  fclose(fp);

  //calculate profile score without alignment
  float score = 0;
  for (int i = 0; i < len; i++)
  {
    AminoName map = charToAmino(dsInfo.seq[i]);
    score += prf[i][map];
  }

  for (int i = 0; i < len; i++) delete[] prf[i];
  delete[] prf;

  return -1.0 * score;
}


float EvolutionScoreAllFromSeq(char* seq)
{
  char evolutiondir[500] = "";
  sprintf(evolutiondir, "%s/evolution", PROGRAM_PATH);
  //printf("evolution parameter file path is: %s\n", evolutiondir);

  SequenceData dsInfo;
  dsInfo.len = strlen(seq);
  dsInfo.seq = new char[dsInfo.len + 1];
  strcpy(dsInfo.seq, seq);

  dsInfo.ss1 = new char[dsInfo.len + 1];
  SSPrediction::getSeq2SS SS(dsInfo.ss1, dsInfo.seq, evolutiondir);
  dsInfo.sa1 = new char[dsInfo.len + 1];
  SAPrediction::getSeq2SA SA(dsInfo.sa1, dsInfo.seq, dsInfo.ss1, evolutiondir);

  //phi-psi prediction
  PhiPsiPrediction::getPhiPsi(dsInfo, evolutiondir);

  dsInfo.ss2 = new char[dsInfo.len + 1];
  FILE* fp = fopen(TGT_SS, "r");
  fscanf(fp, "%s", dsInfo.ss2);
  fclose(fp);
  dsInfo.sa2 = new char[dsInfo.len + 1];
  fp = fopen(TGT_SA, "r");
  fscanf(fp, "%s", dsInfo.sa2);
  fclose(fp);


  wsa /= (dsInfo.len * 1.0);    // weight for solvent accessibility
  wss /= (dsInfo.len * 1.0);    // weight for secondary structure
  wang /= (dsInfo.len * 180.0);  // weight for phi-psi prediction

  CharSeq::SeqAlign sa(dsInfo, TGT_PRF, TGT_PHIPSI, wss, wsa, wang);

  float max = sa.printMtxInfo();
  //printf("%f\n",max);

  delete[] dsInfo.ss1;  delete[] dsInfo.ss2;
  delete[] dsInfo.sa1;  delete[] dsInfo.sa2;
  delete[] dsInfo.seq;

  return -1.0 * max;
}


float EvolutionScorePrfSSSAFromSeq(char* seq)
{
  char evolutiondir[500] = "";
  sprintf(evolutiondir, "%s/evolution", PROGRAM_PATH);
  //printf("evolution parameter file path is: %s\n", evolutiondir);

  SequenceData dsInfo;
  dsInfo.len = strlen(seq);
  dsInfo.seq = new char[dsInfo.len + 1];
  strcpy(dsInfo.seq, seq);

  dsInfo.ss1 = new char[dsInfo.len + 1];
  SSPrediction::getSeq2SS SS(dsInfo.ss1, dsInfo.seq, evolutiondir);
  dsInfo.sa1 = new char[dsInfo.len + 1];
  SAPrediction::getSeq2SA SA(dsInfo.sa1, dsInfo.seq, dsInfo.ss1, evolutiondir);

  //phi-psi prediction
  //PhiPsiPrediction::getPhiPsi(dsInfo, evolutiondir);

  dsInfo.ss2 = new char[dsInfo.len + 1];
  FILE* fp = fopen(TGT_SS, "r");
  fscanf(fp, "%s", dsInfo.ss2);
  fclose(fp);
  dsInfo.sa2 = new char[dsInfo.len + 1];
  fp = fopen(TGT_SA, "r");
  fscanf(fp, "%s", dsInfo.sa2);
  fclose(fp);


  wsa /= (dsInfo.len * 1.0);    // weight for solvent accessibility
  wss /= (dsInfo.len * 1.0);    // weight for secondary structure
  wang /= (dsInfo.len * 180.0);  // weight for phi-psi prediction

  //read profile
  float** prf;
  fp = fopen(TGT_PRF, "r");
  prf = new float* [dsInfo.len];
  for (int i = 0; i < dsInfo.len; ++i)
    prf[i] = new float[AMINOS];
  AminoName map[AMINOS];
  char ami[2];
  for (int i = 0; i < AMINOS; ++i)
  {
    fscanf(fp, "%s", ami);
    map[i] = charToAmino(ami[0]);
  }
  fscanf(fp, "%*s");
  for (int i = 0; i < dsInfo.len; ++i)
  {
    fscanf(fp, "%*s");
    for (int j = 0; j < AMINOS; ++j)
    {
      fscanf(fp, "%f", &prf[i][map[j]]);
    }
    fscanf(fp, "%*s");
  }
  fclose(fp);

  //calculate profile score without alignment
  float score = 0;
  for (int i = 0; i < dsInfo.len; i++)
  {
    AminoName map = charToAmino(seq[i]);
    score += prf[i][map];
    if (dsInfo.ss1[i] == dsInfo.ss2[i]) score += wss * 1.0;
    else if (dsInfo.ss1[i] == '3' && dsInfo.ss2[i] == '3') score += 0;
    else score += wss * (-1.0);
    if (dsInfo.sa1[i] == dsInfo.sa2[i]) score += wsa * 1.0;
    else if (dsInfo.sa1[i] == '2' && dsInfo.sa2[i] == '2') score += 0;
    else score += wsa * (-1.0);
  }

  delete[] dsInfo.ss1;  delete[] dsInfo.ss2;
  delete[] dsInfo.sa1;  delete[] dsInfo.sa2;
  delete[] dsInfo.seq;

  return -1.0 * score;
}



float EvolutionScorePrfFromSeq(char* seq)
{
  char evolutiondir[500] = "";
  sprintf(evolutiondir, "%s/evolution", PROGRAM_PATH);
  //printf("evolution parameter file path is: %s\n", evolutiondir);

  int len = strlen(seq);
  //read profile
  float** prf;
  FILE* fp = NULL;
  fp = fopen(TGT_PRF, "r");
  prf = new float* [len];
  for (int i = 0; i < len; ++i)
    prf[i] = new float[AMINOS];
  AminoName map[AMINOS];
  char ami[2];
  for (int i = 0; i < AMINOS; ++i)
  {
    fscanf(fp, "%s", ami);
    map[i] = charToAmino(ami[0]);
  }
  fscanf(fp, "%*s");
  for (int i = 0; i < len; ++i)
  {
    fscanf(fp, "%*s");
    for (int j = 0; j < AMINOS; ++j)
    {
      fscanf(fp, "%f", &prf[i][map[j]]);
    }
    fscanf(fp, "%*s");
  }
  fclose(fp);

  //calculate profile score without alignment
  float score = 0;
  for (int i = 0; i < len; i++)
  {
    AminoName map = charToAmino(seq[i]);
    score += prf[i][map];
  }

  for (int i = 0; i < len; i++) delete[] prf[i];
  delete[] prf;

  return -1.0 * score;
}


float EvolutionScoreFromPSSMWithoutAlignment(char* seq)
{
  //calculate profile score without alignment
  int len = strlen(seq);
  float score = 0;
  for (int i = 0; i < len; i++)
  {
    AminoName map = charToAmino(seq[i]);
    score += PROT_PROFILE[i][map];
  }

  return -1.0 * score;
}


int SSPred(char* seqfile)
{
  char evolutiondir[500] = "";
  sprintf(evolutiondir, "%s/evolution", PROGRAM_PATH);
  printf("evolution parameter file path is: %s\n", evolutiondir);

  FILE* fp = fopen(seqfile, "r");
  if (fp == NULL)
  {
    printf("In file %s line %d, cannot read file %s\n", __FILE__, __LINE__, seqfile);
    return IOError;
  }
  SequenceData dsInfo;
  dsInfo.len = Text::ncharLine(fp);
  dsInfo.seq = new char[dsInfo.len + 1];
  fseek(fp, 0, SEEK_SET);
  fscanf(fp, "%s\n", dsInfo.seq);
  fclose(fp);


  dsInfo.ss1 = new char[dsInfo.len + 1];
  SSPrediction::getSeq2SS SS(dsInfo.ss1, dsInfo.seq, evolutiondir);
  dsInfo.ss1[dsInfo.len + 1] = '\0';
  FILE* pFile = fopen("SSpred_result.txt", "w");
  for (int i = 0;i < dsInfo.len;i++)
  {
    if (dsInfo.ss1[i] == '1') fprintf(pFile, "H");
    else if (dsInfo.ss1[i] == '2') fprintf(pFile, "E");
    else fprintf(pFile, "C");
  }
  fprintf(pFile, "\n");
  fclose(pFile);

  delete[] dsInfo.ss1;
  delete[] dsInfo.seq;

  return Success;
}


int SAPred(char* seqfile)
{
  char evolutiondir[500] = "";
  sprintf(evolutiondir, "%s/evolution", PROGRAM_PATH);

  FILE* fp = fopen(seqfile, "r");
  if (fp == NULL)
  {
    printf("in file %s line %d, cannot read file %s\n", __FILE__, __LINE__, seqfile);
    return IOError;
  }
  SequenceData dsInfo;
  dsInfo.len = Text::ncharLine(fp);
  dsInfo.seq = new char[dsInfo.len + 1];
  fseek(fp, 0, SEEK_SET);
  fscanf(fp, "%s\n", dsInfo.seq);
  fclose(fp);

  dsInfo.ss1 = new char[dsInfo.len + 1];
  SSPrediction::getSeq2SS SS(dsInfo.ss1, dsInfo.seq, evolutiondir);
  dsInfo.sa1 = new char[dsInfo.len + 1];
  SAPrediction::getSeq2SA SA(dsInfo.sa1, dsInfo.seq, dsInfo.ss1, evolutiondir);
  dsInfo.sa1[dsInfo.len + 1] = '\0';
  FILE* pFile = fopen("SApred_result.txt", "w");
  for (int i = 0;i < dsInfo.len;i++)
  {
    if (dsInfo.sa1[i] == '1') fprintf(pFile, "B");
    else if (dsInfo.sa1[i] == '2') fprintf(pFile, "I");
    else fprintf(pFile, "E");
  }
  fprintf(pFile, "\n");
  fclose(pFile);

  delete[] dsInfo.ss1;
  delete[] dsInfo.sa1;
  delete[] dsInfo.seq;

  return Success;
}

int PhiPsiPred(char* seqfile)
{
  char evolutiondir[500] = "";
  sprintf(evolutiondir, "%s/evolution", PROGRAM_PATH);

  FILE* fp = fopen(seqfile, "r");
  if (fp == NULL)
  {
    printf("In file %s line %d, cannot read file %s\n", __FILE__, __LINE__, seqfile);
    return IOError;
  }
  SequenceData dsInfo;
  dsInfo.len = Text::ncharLine(fp);
  dsInfo.seq = new char[dsInfo.len + 1];
  fseek(fp, 0, SEEK_SET);
  fscanf(fp, "%s", dsInfo.seq);
  fclose(fp);

  dsInfo.ss1 = new char[dsInfo.len + 1];
  SSPrediction::getSeq2SS SS(dsInfo.ss1, dsInfo.seq, evolutiondir);
  dsInfo.sa1 = new char[dsInfo.len + 1];
  SAPrediction::getSeq2SA SA(dsInfo.sa1, dsInfo.seq, dsInfo.ss1, evolutiondir);
  //phi-psi prediction
  PhiPsiPrediction::getPhiPsi(dsInfo, evolutiondir);

  delete[] dsInfo.ss1;
  delete[] dsInfo.sa1;
  delete[] dsInfo.seq;

  return Success;
}