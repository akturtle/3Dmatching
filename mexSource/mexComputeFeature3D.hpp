#define D 3

void computeFeatureSimple( double* pP1, int i , int j, int k , double* pF);
void computeFeature( double* pP1 , int nP1 , double* pP2 , int nP2 ,
                      int* pT1 , int nT1 , double* pF1 , double* pF2);

void computeFeature( double* pP1 , int nP1 , double* pP2 , int nP2 ,
                      int* pT1 , int nT1 , double* pF1 , double* pF2)
{ 
  const int nFeature=3;
  for(int t=0;t<nT1;t++)
  {
    computeFeatureSimple(pP1,pT1[t*3],pT1[t*3+1],pT1[t*3+2],pF1+t*nFeature);
  }
  
  for(int i=0;i<nP2;i++)
    for(int j=0;j<nP2;j++)
      for(int k=0;k<nP2;k++)
        computeFeatureSimple(pP2,i,j,k,pF2+((i*nP2+j)*nP2+k)*nFeature);
  
}

void computeFeatureSimple( double* pP1, int i , int j, int k , double* pF)
{ 
  const int nFeature=3;
  double vecX[nFeature];
  double vecY[nFeature];
  double vecZ[nFeature];
  int ind[nFeature];
  ind[0]=i;ind[1]=j;ind[2]=k;
  double n;
  if((ind[0]==ind[1])||(ind[0]==ind[2])||(ind[1]==ind[2]))
  {
    pF[0]=pF[1]=pF[2]=-10;
    return;
  }
  for(int f=0;f<nFeature;f++)
  {
    vecX[f]=pP1[ind[((f+1)%3)]*2]-pP1[ind[f]*2];
    vecY[f]=pP1[ind[((f+1)%3)]*2+1]-pP1[ind[f]*2+1];
    vecZ[f]=pP1[ind[((f+1)%3)]*2+2]-pP1[ind[f]*2+2];
    double norm=sqrt(vecX[f]*vecX[f]+vecY[f]*vecY[f]+vecZ[f]*vecZ[f]);
    if(norm!=0)
    {
      vecX[f]/=norm;
      vecY[f]/=norm;
      vecZ[f]/=norm;
    }else{
      vecX[f]=0;
      vecY[f]=0;
    }
  }
  for(int f=0;f<nFeature;f++)
  {
    double x,y,z;
      //cross product get the area of unit vector
    x=vecX[f]*vecZ[((f+1)%3)]-vecX[((f+1)%3)]*vecZ[f];
    y=vecY[f]*vecX[((f+1)%3)]-vecY[((f+1)%3)]*vecX[f];
    z=vecY[f]*vecZ[((f+1)%3)]-vecY[((f+1)%3)]*vecZ[f];
    pF[f] = x*x+y*y+z*z;
    //features: unit vector and Z axis[0 0 1] triangle area 
    pF[f+D]=vecX[f]*vecX[f]+vecY[f]+vecY[f]*vecY[f] ;
  }
}