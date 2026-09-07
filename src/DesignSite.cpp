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

#include "DesignSite.h"

int DesignSiteCreate(DesignSite* pThis)
{
  pThis->pRes = NULL;
  RotamersCreate(&pThis->rots);
  pThis->chnNdx = -1;
  pThis->resNdx = -1;
  return Success;
}
int DesignSiteDestroy(DesignSite* pThis)
{
  pThis->pRes->desType = Type_DesType_Fixed;
  pThis->pRes = NULL;
  RotamersDestroy(&pThis->rots);
  pThis->chnNdx = -1;
  pThis->resNdx = -1;
  return Success;
}
Rotamers* DesignSiteGetRotamers(DesignSite* pThis)
{
  return &pThis->rots;
}
int DesignSiteShow(DesignSite* pThis, FILE* pFile)
{
  printf("Chain %s residue %s %d has %d rotamers\n", ResidueGetChainName(pThis->pRes), ResidueGetName(pThis->pRes), ResidueGetPosInChain(pThis->pRes), RotamersGetNum(&pThis->rots));
  return Success;
}


char* DesignSiteGetChainName(DesignSite* pThis)
{
  return pThis->pRes->chainName;
}
char* DesignSiteGetResiName(DesignSite* pThis)
{
  return pThis->pRes->name;
}
int DesignSiteGetPosInChain(DesignSite* pThis)
{
  return pThis->pRes->posInChain;
}

Residue* DesignSiteGetResidue(DesignSite* pThis)
{
  return pThis->pRes;
}

int DesignSiteShowRepresentativeRotamerAtomParam(DesignSite* pThis)
{
  Rotamers* pSet = DesignSiteGetRotamers(pThis);
  for (int i = 0; i < RotamersGetRepNum(pSet); ++i)
  {
    Rotamer* pRotamer = RotamersGetRepByIndex(pSet, i);
    printf("Rotamer %s Atom Parameter:\n", RotamerGetType(pRotamer));
    RotamerShowAtomParameter(pRotamer);
  }
  return Success;
}


int DesignSiteShowRepresentativeRotamerBondInfo(DesignSite* pThis)
{
  Rotamers* pSet = DesignSiteGetRotamers(pThis);
  for (int i = 0; i < RotamersGetRepNum(pSet); ++i)
  {
    Rotamer* pRotamer = RotamersGetRepByIndex(pSet, i);
    printf("Rotamer %s Bond information:\n", RotamerGetType(pRotamer));
    RotamerShowBondInformation(pRotamer);
  }
  return Success;
}

int DesignSiteRemoveRotamers(DesignSite* pThis)
{
  RotamersDestroy(&pThis->rots);
  return RotamersCreate(&pThis->rots);
}


int DesignSiteCopy(DesignSite* pThis, DesignSite* pOther)
{
  pThis->chnNdx = pOther->chnNdx;
  pThis->resNdx = pOther->resNdx;
  pThis->pRes = pOther->pRes;
  RotamersCopy(&pThis->rots, &pOther->rots);
  return Success;
}