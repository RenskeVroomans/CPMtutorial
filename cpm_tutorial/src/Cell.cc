#include "Cell.hh"


Cell::Cell()
{
  
  neighbours.clear();
  vol=0;
  tarvol=0;
  surf=0;
  tarsurf=0;
  age=0;
  
  EmptyMoments();
  iyy=0.;
  ixx=0.;
  ixy=0.;
  rhs1=0.;
  rhs2=0.;
  tveci=0.;
  tvecj=0.;
  previ=0.;
  prevj=0.;
}

Cell::~Cell()
{
  //printf("destructing cell\n");
}

void Cell::CreateCell(int cid, int type)
{
  id=cid;
  age=0;
  celltype=type;
  neighbours.clear();
  EmptyMoments();
  iyy=0.;
  ixx=0.;
  ixy=0.;
  rhs1=0.;
  rhs2=0.;
  tveci=0.;
  tvecj=0.;
}

Cell::Cell(const Cell &obj)
{
  //copy vars
  id=0;
  age=0;
  celltype=obj.celltype;
  vol=obj.vol;
  tarvol=obj.tarvol;
  surf=obj.surf;
  tarsurf=obj.tarsurf;
  
  
  Sxx=obj.Sxx;
  Syy=obj.Syy;
  Sx=obj.Sx;
  Sxy=obj.Sxy;
  Sy=obj.Sy;
  iyy=obj.iyy;
  ixx=obj.ixx;
  ixy=obj.ixy;
  rhs1=obj.rhs1;
  rhs2=obj.rhs2;
  
  neighbours.clear();
}

Cell& Cell::operator=(const Cell& celltocopy)
{
  id=0; //note: id needs to be set separately!
  age=0;
  celltype=celltocopy.celltype;
  celltype=celltocopy.celltype;
  vol=celltocopy.vol;
  tarvol=celltocopy.tarvol;
  surf=celltocopy.surf;
  tarsurf=celltocopy.tarsurf;
  
  Sxx=celltocopy.Sxx;
  Syy=celltocopy.Syy;
  Sx=celltocopy.Sx;
  Sxy=celltocopy.Sxy;
  Sy=celltocopy.Sy;
  iyy=celltocopy.iyy;
  ixx=celltocopy.ixx;
  ixy=celltocopy.ixy;
  rhs1=celltocopy.rhs1;
  rhs2=celltocopy.rhs2;
  
  neighbours.clear();
  return *this;
}

void Cell::StartTarVec()
{
   
    prevj=meanj;
    previ=meani;
}

void Cell::UpdateTarVec()
{
    tvecj=meanj-prevj;
    tveci=meani-previ;
    tvecj/=hypot(tveci,tvecj);
    tveci/=hypot(tveci,tvecj);
    prevj=meanj;
    previ=meani;
}

void Cell::UpdateMoments(int i,int j, int add)
{
 if(add>0){
    Sx+=j;
    Sy+=i;
    Sxx+=(j*j);
    Syy+=(i*i);
    Sxy+=(i*j);
  }

  else{
    Sx-=j;
    Sy-=i;
    Sxx-=(j*j);
    Syy-=(i*i);
    Sxy-=(i*j);
  }
}
void Cell::EmptyMoments(void)
{
  Sx=0;
  Sy=0;
  Sxx=0;
  Syy=0;
  Sxy=0;
  
}
void Cell::UpdateCellPixel(int i, int j, int add)
{
  //volume should be updated already!
  UpdateMoments(i,j,add);
  if(vol>0){
    meani=Sy/vol;
    meanj=Sx/vol;
  }
  else{
    printf("Cell %d has died\n",id);
    return;
  }

  UpdateInertiaTensor();
  
}

void Cell::UpdateCell(void)//assumes all Moments are already up to date: useful for initialization
{
  if(vol>1){
    meani=Sy/vol;
    meanj=Sx/vol;
    
    UpdateInertiaTensor();
  }
 
}
void Cell::UpdateInertiaTensor(void)
{ 
  long double temp1;
  // inertia tensor (constructed from the raw momenta)
  if(vol>1){
    temp1=Sx*Sx;
    iyy=Sxx-temp1/(long double)vol;
     temp1=Sy*Sy;
    ixx=Syy-temp1/(long double)vol;
   
    temp1=Sx*Sy;
    ixy=-Sxy+temp1/(long double)vol;
  
    rhs1=(ixx+iyy)*0.5;
    rhs2=sqrt( (ixx-iyy)*(ixx-iyy)+4*ixy*ixy )*0.5;
  }

}

double Cell::ReturnMajorAxisLength(void)//taken from Roeland Merks' code
{    
  UpdateInertiaTensor();
  
  double lambda_b=rhs1+rhs2;
  
  return 4*sqrt(lambda_b/vol);
}

double Cell::ReturnMinorAxisLength(void)//taken from Roeland Merks' code
{  
  UpdateInertiaTensor();  
  double lambda_a=rhs1-rhs2;
   
  if(lambda_a<-0.0001){
    printf("ReturnMinorAxisLength, cell %d: error: small eigenvalue negative\n",id);
    printf("rhs1: %.2Lf, rhs2: %.2Lf, vol: %d\n",rhs1,rhs2,vol);
    printf("Sxx %Lf,Sx %Lf, Sxy %Lf, Syy %Lf, Sy %Lf\n",Sxx, Sx, Sxy, Syy, Sy);

    return -1;
  }
  return 4*sqrt(lambda_a/vol);
}

double Cell::ReturnLongval(void)//taken from Roeland Merks' code
{  
  UpdateInertiaTensor();  

  return rhs1+rhs2;
 
}
double Cell::ReturnShortval(void)//taken from Roeland Merks' code
{  
  UpdateInertiaTensor();  

  return rhs1-rhs2;
 
}

//Careful! vector belonging to large eigenvalue points along MINOR axis and vice versa!!
FPosition  Cell::ReturnLongVec(void)
{
  FPosition longvec;
  double hyp;
  UpdateInertiaTensor();  
  longvec.xx=ixy;
  longvec.yy=rhs1-rhs2-ixx;

  hyp=sqrt(longvec.xx*longvec.xx+longvec.yy*longvec.yy);
  if(hyp<0.001){
    printf("longaxis of cell %d does not exist. Something is wrong!!\n",id);//now we calculate components of shortaxis and do perpendicular instead.
    longvec.xx=-rhs1+rhs2-ixx;
    longvec.yy=ixy;
    hyp=sqrt(longvec.xx*longvec.xx+longvec.yy*longvec.yy);
    if((hyp-1.00)*(hyp-1.00)>0.1){
      longvec.xx=longvec.xx/hyp;
      longvec.yy=longvec.yy/hyp;
    }
  }
  else if((hyp-1.00)*(hyp-1.00)>0.1){
    longvec.xx=longvec.xx/hyp;
    longvec.yy=longvec.yy/hyp;
  }
  return longvec;
}


FPosition  Cell::ReturnShortVec(void)
{
  FPosition shortvec;
  double hyp;
  UpdateInertiaTensor();  
  shortvec.xx=ixy;
  shortvec.yy=rhs1+rhs2-ixx;
  hyp=sqrt(shortvec.xx*shortvec.xx+shortvec.yy*shortvec.yy);
  if(hyp>1.00){//normalize the vector to length 1
   shortvec.xx=shortvec.xx/hyp;
   shortvec.yy=shortvec.yy/hyp;
  }
  return shortvec;
}

FPosition Cell::DetermineMomentofInertia(int i, int j, int add)
{
  
  UpdateInertiaTensor(); 

  double inertold=ixx+iyy;

  long double tx,ty,txx,tyy,txy;
  
  int newvol;

  FPosition oldnew;

  newvol=vol+add;
  tx=Sx+(add*j);
  ty=Sy+(add*i);
  txx=Sxx+(add*j*j);
  txy=Sxy+(add*i*j);
  tyy=Syy+(add*i*i); 
  
  tx=Sx+(add*j);
  ty=Sy+(add*i);
  txx=Sxx+(add*j*j);
  txy=Sxy+(add*i*j);
  tyy=Syy+(add*i*i);
 
  iyy=txx-tx*tx/(long double)newvol;
  ixx=tyy-ty*ty/(long double)newvol;
  ixy=-txy+tx*ty/(long double)newvol;

  double inertnew=ixx+iyy;

  oldnew.xx=inertold;
  oldnew.yy=inertnew;

  return oldnew; 

}

FPosition Cell::ReturnChangeInAspect(int i, int j, int add)
{
  double initial, newval;
  long double tx,ty,txx,tyy,txy;
  double ths1, ths2;
  int newvol;

  FPosition oldnew;

  //assume the info is updated
  UpdateInertiaTensor();  
  initial=(rhs1-rhs2)/(rhs1+rhs2);

  newvol=vol+add;
  tx=Sx+(add*j);
  ty=Sy+(add*i);
  txx=Sxx+(add*j*j);
  txy=Sxy+(add*i*j);
  tyy=Syy+(add*i*i);
 
  iyy=txx-tx*tx/(long double)newvol;
  ixx=tyy-ty*ty/(long double)newvol;
  ixy=-txy+tx*ty/(long double)newvol;
  
  ths1=(ixx+iyy)*0.5;
  ths2=sqrt( (ixx-iyy)*(ixx-iyy)+4*ixy*ixy )*0.5;
  
  newval=(ths1-ths2)/(ths1+ths2);
  
  oldnew.xx=initial;
  oldnew.yy=newval;

  return oldnew;
}


