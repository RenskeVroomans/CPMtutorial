#include "Agent.hh"
#include "Misc.hh"

void Agent::ShuffleField(vector<int> &shuffle)
{
  int rand, temp;
  int LL=L-10;
  int WW=W-10;
  
  for(int i=0; i<LL*WW;i++)
    shuffle.push_back(i);
  
  for(int i=LL*WW-1; i>0;i--)
  {
    rand=(int)(uniform()*(i+1));
    temp=shuffle[i];
    shuffle[i]=shuffle[rand];
    shuffle[rand]=temp;
  }

//   for(int i=0; i<LL*WW; i++)
//   {
//     int found=0;
//     for(int j=0; j<LL*WW; j++)
//     {
//       if(shuffle[j]==i)
//       {
// 	found=1;
// 	break;
//       }
//     }
//     if(!found) 
//       cout<<"oh no! "<<i <<" not found\n";
//   }
}

void Agent::PosAndNeigh(int num, int *i, int *j, int *nbi, int *nbj)
{
  (*i)=num/(W-10)+5;
  (*j)=num%(W-10)+5;
  int array[8]={0,1,2,3,5,6,7,8};
  
  int rand=(int)(uniform()*8);
  (*nbi)=array[rand]/3+(*i)-1;
  (*nbj)=array[rand]%3+(*j)-1;
  
}

void Agent::randomCouple(int *i, int *j, int *nbi, int *nbj)
{
  static int first_time=0;
  static int NBpos[8][2];
  static int LL=L-10;
  static int WW=W-10;

  if(first_time==0)
    {
      NBpos[0][0]=-1;
      NBpos[0][1]=-1;
      NBpos[1][0]=0;
      NBpos[1][1]=-1;
      NBpos[2][0]=1;
      NBpos[2][1]=-1;
      NBpos[3][0]=1;
      NBpos[3][1]=0;
      NBpos[4][0]=1;
      NBpos[4][1]=1;
      NBpos[5][0]=0;
      NBpos[5][1]=1;
      NBpos[6][0]=-1;
      NBpos[6][1]=1;
      NBpos[7][0]=-1;
      NBpos[7][1]=0;
      first_time=1;
    }

  int z=(int)(uniform()*(LL*WW*8));
  (*i)=z/(WW*8);
  (*j)=(z-(*i)*WW*8)/8;
  int offset=z-(*i)*WW*8-(*j)*8;
  (*i)+=5;
  (*j)+=5;
  (*nbi)=(*i)+NBpos[offset][0];
  (*nbj)=(*j)+NBpos[offset][1];

  //printf("%i %i and %i %i\n",*i,*j,*nbi,*nbj);
}
void Agent::randomMaskCouple(int *i, int *j, int *nbi, int *nbj, int imin, int imax, int jmin, int jmax)
{
  static int first_time=0;
  static int NBpos[8][2];
  int LL=(imax-imin+1);
  int WW=(jmax-jmin+1);

  if(first_time==0)
    {
      NBpos[0][0]=-1;
      NBpos[0][1]=-1;
      NBpos[1][0]=0;
      NBpos[1][1]=-1;
      NBpos[2][0]=1;
      NBpos[2][1]=-1;
      NBpos[3][0]=1;
      NBpos[3][1]=0;
      NBpos[4][0]=1;
      NBpos[4][1]=1;
      NBpos[5][0]=0;
      NBpos[5][1]=1;
      NBpos[6][0]=-1;
      NBpos[6][1]=1;
      NBpos[7][0]=-1;
      NBpos[7][1]=0;
      first_time=1;
    }

  int z=(int)(uniform()*(LL*WW*8));
  (*i)=z/(WW*8);
  (*j)=(z-(*i)*WW*8)/8;
  int offset=z-(*i)*WW*8-(*j)*8;
  (*i)+=imin;
  (*j)+=jmin;
  (*nbi)=(*i)+NBpos[offset][0];
  (*nbj)=(*j)+NBpos[offset][1];

  //printf("%i %i and %i %i\n",*i,*j,*nbi,*nbj);
}

void Agent::randomLocalCouple(int *i, int *j, int *ii,int *jj, int *nbii,int *nbjj)
{
  static int first_time=0;
  static int NBpos[8][2];

  if(first_time==0)
    {
      NBpos[0][0]=-1;
      NBpos[0][1]=-1;
      NBpos[1][0]=0;
      NBpos[1][1]=-1;
      NBpos[2][0]=1;
      NBpos[2][1]=-1;
      NBpos[3][0]=1;
      NBpos[3][1]=0;
      NBpos[4][0]=1;
      NBpos[4][1]=1;
      NBpos[5][0]=0;
      NBpos[5][1]=1;
      NBpos[6][0]=-1;
      NBpos[6][1]=1;
      NBpos[7][0]=-1;
      NBpos[7][1]=0;
      first_time=1;
    }
 
  int z=(int)(uniform()*8*8);
  int offset1=z/8;
  int offset2=z-offset1*8;
  *ii=(*i)+NBpos[offset1][0];
  *jj=(*j)+NBpos[offset1][1];
  *nbii=(*ii)+NBpos[offset2][0];
  *nbjj=(*jj)+NBpos[offset2][1];
}

int Agent::CPMAttempt(int i,int j,int ii,int jj,double extravalue)
{
  double DeltaHam;
  double c;

  DeltaHam=DeltaSurfHamiltonian(i,j,ii,jj);

  DeltaHam+=extravalue;

  if(DeltaHam>0)
    c=exp(-DeltaHam*invT);
  else
    c=1;

  
  if(uniform()<c)//copy action is performed
    return TRUE;
  else
    return FALSE;
    
}

void Agent::DoCPMUpdate(int i,int j,int ii,int jj)
{
  int posNB[8][2];
  int NB;
  int idA;
  int idB;
  int idN;

  posNB[0][0]=i-1;
  posNB[0][1]=j-1;
  posNB[1][0]=i+0;
  posNB[1][1]=j-1;
  posNB[2][0]=i+1;
  posNB[2][1]=j-1;
  posNB[3][0]=i+1;
  posNB[3][1]=j+0;
  posNB[4][0]=i+1;
  posNB[4][1]=j+1;
  posNB[5][0]=i+0;
  posNB[5][1]=j+1;
  posNB[6][0]=i-1;
  posNB[6][1]=j+1;
  posNB[7][0]=i-1;
  posNB[7][1]=j+0;


  //Store cell ids of original configuration
  idA=CellIdGrid[i][j];
  idB=CellIdGrid[ii][jj];
  
  //update volumes
  if(idA!=0)
    CellArray[idA].vol--;
  if(idB!=0)
    CellArray[idB].vol++;

  
  if(neighbourhoodsize && surfaceconstraint){
    for(NB=0;NB<neighbourhoodsize;NB++)
      {
	idN=CellIdGrid[i+wideneighbourhood[NB].yy][j+wideneighbourhood[NB].xx]; 
	if(idN==idA)//now inside same cell, after change in different cell
	  {
	    if(idA!=0 && idB!=0)
	      {
		CellArray[idA].surf++;
		CellArray[idB].surf++;
	      }
	    else if(idA!=0)
	      CellArray[idA].surf++;
	    else if(idB!=0)
	      CellArray[idB].surf++;
	  }
	else if(idN==idB)//now inside different cell, after change in same cell
	  {
	    if(idA!=0 && idB!=0)
	      {
		CellArray[idA].surf--;
		CellArray[idB].surf--;
	      }
	    else if(idA!=0)
	      CellArray[idA].surf--;
	    else if(idB!=0)
	      CellArray[idB].surf--;
	    
	  }
	else//now in different cells, after change in different different cells
	  {
	    if(idA!=0 && idB!=0)
	      {
	      CellArray[idA].surf--;
	      CellArray[idB].surf++;
	      }
	    else if(idA!=0)
	      CellArray[idA].surf--;
	    else if(idB!=0)
	      CellArray[idB].surf++;
	  }
      }
    
  }
  


  //update mean positions and other cell info:
  // DetermineSingleCellinfo(idA);
  //DetermineSingleCellinfo(idB);
}

void Agent::InitContactLength()
{
  int i, j,k;
  int celltype, celltypeneigh;
  int boundary=0;
  int totalperimeter=0;
  int sigma, sigmaneigh;
  map<int, Cell>::iterator it;
  
  for(it=CellArray.begin(); it!=CellArray.end(); ++it)
  {
    it->second.clearNeighbours();
  }
  
  for(i=0;i<L;i++)
  {
    for(j=0;j<W;j++)
    {
      if((sigma=CellIdGrid[i][j]) && (celltype=CellArray[CellIdGrid[i][j]].celltype)) //focus is on a cell
      {
	for(k=0; k<neighbourhoodsize; k++)//go through neighbourhood of the pixel
	{
	  if((sigmaneigh=CellIdGrid[i+wideneighbourhood[k].yy][j+wideneighbourhood[k].xx])!=sigma)//medium can also be a neighbour!
	  {
	    totalperimeter++;
	    CellArray[sigma].updateNeighbour(sigmaneigh,1);
	  }
	}
      }
    }
  }
  
}

//after divisions, neighbours need to be re-established. !!! need to have the old position of the cell
// so call this function before UpdateCell in the division functions.
void Agent::RecountContactLength(int id, int newid)
{
  int i, j,k;
  int sigma, sigmaneigh;

  if(id==0 || newid==0)
  {
    printf("RecountContactLength: error, updating medium for contacts\n");
    exit(1);
  }

  set<int> store;
  set<int> store2;
  double scale=3.;
  ///find the section of the field over which to iterate (taken quite big)
  int minL=(int)max(CellArray[id].meani-scale*CellArray[id].ReturnMajorAxisLength(), 4.);
  int maxL=(int)min(CellArray[id].meani+scale*CellArray[id].ReturnMajorAxisLength(), (double)L-5);
  int minW=(int)max(CellArray[id].meanj-scale*CellArray[id].ReturnMajorAxisLength(), 4.);
  int maxW=(int)min(CellArray[id].meanj+scale*CellArray[id].ReturnMajorAxisLength(), (double)W-5);
 
  ///go through the neighbours of the old cell, and remove the contact with ID
  map<int, int>::iterator it;
  int checkbounds=0;
  int count;
  for(it=CellArray[id].neighbours.begin(),count=0; it!=CellArray[id].neighbours.end(); ++it, count++)
  {
    checkbounds+=it->second;
    store.insert(it->first);
    if(it->first)
      CellArray[it->first].setNeighbour(id, 0); //this neighbour removes cell ID from its contacts
      
  }
  CellArray[id].clearNeighbours();
  
  int check=0;
  int checkbounds2=0;
  ///reset the contact by iterating over the field
  for(i=minL;i<=maxL;i++)
  {
    for(j=minW;j<=maxW;j++)
    {
      if((sigma=CellIdGrid[i][j])==id || sigma==newid) //focus is on the cell or its daughter
      {
	for(k=0; k<neighbourhoodsize; k++)//go through neighbourhood of the pixel
	{
	  if((sigmaneigh=CellIdGrid[i+wideneighbourhood[k].yy][j+wideneighbourhood[k].xx])!=sigma){
	    check+=CellArray[sigma].updateNeighbour(sigmaneigh,1);
	    if(sigmaneigh!=id && sigmaneigh!=newid)
	    {
	      checkbounds2++;
	      store2.insert(sigmaneigh);
	    }
	    if(sigmaneigh!=id && sigmaneigh!=newid && sigmaneigh){//also update the contacts of the neighbour, unless it is one of the two focus cells:otherwise you update those twice
	      check+=CellArray[sigmaneigh].updateNeighbour(sigma,1);
	    }
	  }
	  if(check)
	  {
	    printf("error in RecountContactlength: updating pixels sigma %d and sigmaneigh %d\n", sigma, sigmaneigh);
          }
	}
	if(i==minL || i==maxL || j==minW || j==maxW )
	  printf("RecountContactLength: caution: cell in question at border of field. cell: %d, pos: %d %d\n", id, i, j);
      }
    }
  }
  
  if(checkbounds2!=checkbounds)
  {
    printf("error in RecountContactLength: wrong nr of boundary pixels! divided %d to %d,  checkbounds %d, checkbounds2 %d\n",id, newid, checkbounds, checkbounds2);
    for(i=minL;i<maxL;i++)
    {
      for(j=minW;j<maxW;j++)
      {
	printf("%d\t", CellIdGrid[i][j]);
      }
      printf("\n");
    }
   // exit(1);
  }
  if(store!=store2)
  {
    printf("error in RecountContactLength: different registered neighbours!\n");
    ///for debugging:
//     set<int> diff1;
//     printf("in store 1 but not store 2: ");
//     set_difference(store.begin(), store.end(), store2.begin(), store2.end(),inserter(diff1, diff1.begin()));
//   
//     for (set<int>::iterator dit=diff1.begin(); dit!=diff1.end(); ++dit) printf("%d\t", (*dit));
//     printf("\n");
//     
//     set<int> diff2;
//     printf("in store 2 but not store 1: ");
//     set_difference(store2.begin(), store2.end(), store.begin(), store.end(),inserter(diff2, diff2.begin()));
//   
//     for (set<int>::iterator dit=diff2.begin(); dit!=diff2.end(); ++dit) printf("%d\t", (*dit));
//     printf("\n");
    InitContactLength();
  }
  

 
}