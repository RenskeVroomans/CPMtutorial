#include "Agent.hh"

double Agent::DeltaSurfHamiltonian(int I1,int J1,int I2,int J2)
{
  // double Jtemp1,Jtemp2;
  double DeltaHam=0;
  double DeltaSurf=0;
  double DeltaVol=0;
  double DeltaJ=0;
  double temp1=0,temp2;
  int idA=0;
  int idB=0;
  int idN=0;
  int NB;
  int typeA, typeB, typeN;

  //int posNB[8][2];

  idA=CellIdGrid[I1][J1];
  idB=CellIdGrid[I2][J2];

  typeA=CellArray[idA].celltype;
  typeB=CellArray[idB].celltype;

  if(!neighbourhoodsize){
    printf("DeltaSurfHamiltonian:error: no neighbourhood defined. Use FindNeighbourhood first.\n");
    exit(1);
  }
  
  //compute DeltaJ
  for(NB=0;NB<neighbourhoodsize;NB++)
    {    
      idN=CellIdGrid[I1+wideneighbourhood[NB].yy][J1+wideneighbourhood[NB].xx];
      typeN=CellArray[idN].celltype;
      if(idN!=idA)//now inside same cell, after change in different cell
	{
	  temp1-=JTable[typeN][typeA]; //typeN, typeA
	}
      if(idN!=idB)//now inside different cell, after change in same cell
	{
	  temp1+=JTable[typeN][typeB];//typeN, typeB
	}
     
      DeltaJ+=temp1;
      temp1=0;
    }
 
 
  //compute DeltaVol
  if(idA!=0)
    {
      temp1=(CellArray[idA].vol-CellArray[idA].tarvol);
      temp2=temp1-1;
      temp1=temp1*temp1;
      temp2=temp2*temp2;
      DeltaVol+=labdavol*(temp2-temp1);
    }
  if(idB!=0)
    {
      temp1=(CellArray[idB].vol-CellArray[idB].tarvol);
      temp2=temp1+1;
      temp1=temp1*temp1;
      temp2=temp2*temp2;
      DeltaVol+=labdavol*(temp2-temp1);
    }
  
  
  //Compute DeltaSurf
  int delsurfpos=0,delsurfneigh=0;
  if(surfaceconstraint){
    //printf("here\n");
    for(NB=0;NB<neighbourhoodsize;NB++)
      {
	idN=CellIdGrid[I1+wideneighbourhood[NB].yy][J1+wideneighbourhood[NB].xx];
	if(idN==idA)//now inside same cell, after change in different cell
	  delsurfpos++;
	else
	  delsurfpos--;
	
	if(idN==idB)//now inside different cell, after change in same cell
	  delsurfneigh--;
	else
	  delsurfneigh++;
	
     }
    if(idA)//if not medium
      DeltaSurf+=labdasurf*(delsurfpos*delsurfpos+delsurfpos*2*(CellArray[idA].surf-CellArray[idA].tarsurf));
    if(idB)//if not medium
      DeltaSurf+=labdasurf*(delsurfneigh*delsurfneigh+delsurfneigh*2*(CellArray[idB].surf-CellArray[idB].tarsurf));
  }

  DeltaHam=DeltaJ+DeltaVol+DeltaSurf;

  return DeltaHam;
}