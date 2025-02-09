#include "Agent.hh"

double   Agent::labdavol;
double   Agent::labdasurf;
double   Agent::invT;
int      Agent::targetvolume;
int      Agent::targetsurface;
int      Agent::surfaceconstraint; //flag to see if a surfaceconstraint is used.
int      Agent::L;               //grid length
int      Agent::W;               //grid width
int      Agent::neighbourhoodsize; //if neighbourhoodsize is given as 0, the Moore neighbourhood is used. Necessary when using DeltaSurfHamiltonian
int      Agent::MaxNrCells;      //256
Position* Agent::wideneighbourhood;  //array holding the positions of the wider neighbourhood used for hamiltonian stuff
int      Agent::InitNrCells;
int      Agent::NrDevSteps;
int**    Agent::JTable;
int      Agent::picinterval;
int      Agent::zoom;
int      Agent::NrCellTypes;
int      Agent::CellPlace;
int      Agent::perstime;
double   Agent::mu;

Agent::Agent()
{
  CellArray.clear();
  piccounter=0;
  piccounter2=0;
  anrcells_=0;
  
  CreateField();
  
 
  
}

Agent::~Agent()
{
  CellArray.clear();
 
  DestroyField();
  
 
   free(JTable[0]);
   free(JTable);
}



//place initial cells in the field; if you need info on contacts between cells, use InitContactLength afterwards
void Agent::PlaceCellsInGrid()
{
  int type;
  int countfail=0;
  int succ;
  ///create medium as special cell type
  CellArray.insert(make_pair(0, Cell()));
  CellArray[0].CreateCell(0,0);
  
  //place cells at random spot in the field
  if (!CellPlace)
  {
    int zz, ii,jj;
    static int LL=L-14;
    static int WW=W-14;
    int ll=0;
    while(ll<InitNrCells)
    {
        zz=(int)(uniform()*(LL*WW*8));
        ii=zz/(WW*8);
        jj=((zz-ii*WW*8)/8);
        ii+=7;
        jj+=7;
        type=int(uniform()*double(NrCellTypes))+1;
        succ=PlaceOneCell(type,ll+1,ii, jj);
        if(succ)
        { 
            countfail++;
        }
        else
        {
            ll++;
        }
        if (countfail>1000)
        {
            printf("PlaceCellsInGrid:error: could not place all cells due to overlap\n Consider increasing grid size\n");
            exit(1);
        }
    }
  }
  else  //place cells in clump
  {
    int counter=1;
    //determine block size
    while(InitNrCells-counter*counter>0){
      counter++;
    }
    int iborder=(L-7*counter)/2;
    int jborder=(W-7*counter)/2;

    if(iborder<7 ||jborder<7){
      printf("PlaceCellsInGrid:error: too many cells, they cannot be placed\n");
      return;
    }

    int placed=1;
    for(int p=0;p<counter;p++){
      
      if(placed>InitNrCells)
	break;
      
      for(int q=0;q<counter;q++){
	
	if(placed>InitNrCells)
	  break;

	int ny=iborder+2+7*p;
	int nx=jborder+2+7*q;
        
        type=int(uniform()*double(NrCellTypes))+1;
        if(PlaceOneCell(type,placed,ny, nx))
            exit(1);
        placed++;
      }
    }
  }
  
  //initialise cell contacts
  //InitContactLength();
  
  
}

void Agent::DevelopAgent(char *show)//run the CPM dynamics + whatever else needs to happen meanwhile (divisions, network updates, pictures) 
{
  map<int, Cell>::iterator it;
  
  PlaceCellsInGrid();

  for(int step=0; step<NrDevSteps; step++)
  {
   
   
       
   //do migration after cells have grown to target size
   if (step==300)
   {
     for(it=CellArray.begin(); it!=CellArray.end(); ++it)
     {
       it->second.StartTarVec();  
       
     }
   }
   
   if(step>300 && mu>0.0001)
   {   
     
    for(it=CellArray.begin(); it!=CellArray.end(); ++it)
    {
       if((it->first%perstime)==(step%perstime))
       {
         it->second.UpdateTarVec();  
       }
     }
   }
   
   
   //make pictures of the field during development, if the target directory is defined
   if(strlen(show)>1 &&  !(step%picinterval))
   {
     Snapshot(zoom,show);
     
   }
    
     /// update the CPM
   UpdateAgent_Random();
    UpdateCellAge();
  }
 
}

void Agent::UpdateCellAge()
{
  map<int, Cell>::iterator it;
  for(it=CellArray.begin(); it!=CellArray.end(); ++it)
  {
    it->second.age++;
  }
  
}


double Agent::Extra(int i,int j, int ii, int jj)//here you can expand the Hamiltonian with more functions
{

  int idA, idB;
  double ai, aj;
  
  double DeltaTar=0;
  
  //cell id's of the cell into which you are copying (A), and the cell which is copying (B)
  idA=CellIdGrid[i][j];
  idB=CellIdGrid[ii][jj];

  
  //Joost's method
    if(mu>0.0001){
      if(idA){
	ai=i-CellArray[idA].meani;
	aj=j-CellArray[idA].meanj;
	DeltaTar+=mu*(ai*CellArray[idA].tveci+aj*CellArray[idA].tvecj)/hypot(ai,aj);
      }
      if(idB){
	aj=j-CellArray[idB].meanj;
	ai=i-CellArray[idB].meani;
	DeltaTar-=mu*(ai*CellArray[idB].tveci+aj*CellArray[idB].tvecj)/hypot(ai,aj);
      }
    }
  
 
 
  return DeltaTar;
}

void Agent::After(int i,int j, int ii, int jj)//extra book keeping if you need it.
{
  int idA, idB;
  idA=CellIdGrid[i][j];
  idB=CellIdGrid[ii][jj];

  ///use this if positional and/or shape info is required (e.g. if you need divisions); take care of dead cells.
  if(idA){
    if(CellArray[idA].vol<1)
    {
      //printf("cell %d has no vol!\n",idA);
      CellArray.erase(idA);
    }
    else
      CellArray[idA].UpdateCellPixel(i,j,-1);
  }
  if(idB){
    if(CellArray[idB].vol<1)
    {
      //printf("cell %d has no vol!\n",idB);
      CellArray.erase(idB);
    }
    else
      CellArray[idB].UpdateCellPixel(i,j,1);
    
  }
  

  
  ///use this if you need info about contacts with neighbours
   // int check=0;
   // int point;
//   for(int k=0; k<neighbourhoodsize; k++)
//   {
//     point=CellIdGrid[i+wideneighbourhood[k].yy][j+wideneighbourhood[k].xx];
//     
//     if(point!=idA)
//     {
//       if(idA && CellArray.count(idA))
// 	check+=CellArray[idA].updateNeighbour(point, -1);
//       if(point)
// 	check+=CellArray[point].updateNeighbour(idA, -1);
//     }
//     if(point!=idB)
//     {
//       if(idB && CellArray.count(idB))
// 	check+=CellArray[idB].updateNeighbour(point, 1);
//       if(point)
// 	check+=CellArray[point].updateNeighbour(idB, 1);
//     }
//     if(check)
//     {
//       printf("error in After: wrongly updating neighbours of copy event idA %d and idB %d\n", idA, idB);
//       printf("agent nr %d\n", agentid);
//     }
//   }
 
  
}

int Agent::PlaceOneCell(int type, int id, int ipos, int jpos)//cell type, cell id, i position, j position
{
  // int stripelength=20; //nr of cells over the length of a stripe
  int i, j,m;
   
  if(ipos<7 || ipos>L-7 || jpos<7 || jpos>W-7)
  {
       printf("too close to rim\n");
       return 1;
  }
  
  //check if spot is taken.
  for(i=-2;i<=2;i++)
  {
    for(j=-2;j<=2;j++)
    {

      
      if(CellIdGrid[ipos+i][jpos+j])
      {
// 	printf("Agent:PlaceOneCell is asked to place cell %d on top of another. exiting\n",id);
// 	exit(1);
        return 1;  
      }
    }
  }
  
  if(!CellArray.count(id))
  {
    CellArray.insert(make_pair(id, Cell()));
    CellArray[id].CreateCell(id, type);
//     printf("created cell\n");
    
  }
  else
  {
    printf("Agent:PlaceOneCell: cell of this id already exists\n");
    //exit(1);
    return 1;
  }
  
  for(i=-2;i<=2;i++)
  {
    for(j=-2;j<=2;j++)
    {
      CellIdGrid[ipos+i][jpos+j]=id;
      CellArray[id].UpdateMoments(ipos+i,jpos+j,1);
      CellArray[id].vol++;
    }
  }
  
  CellArray[id].tarvol=targetvolume;
  
  if(surfaceconstraint){
    for(i=-2;i<=2;i++)
      for(j=-2;j<=2;j++)
	for(m=0;m<neighbourhoodsize;m++){
	  if(CellIdGrid[ipos+i+wideneighbourhood[m].yy][jpos+j+wideneighbourhood[m].xx]!=CellIdGrid[ipos+i][jpos+j])
	    CellArray[id].surf++;
	}
    CellArray[id].tarsurf=targetsurface;
  }
  
  CellArray[id].UpdateCell();
  anrcells_++;
  
  return 0;
}
 

void Agent::CreateAgentManyCells(int number, int type, int clumped) //creates a number of cells of a certain type
{
  int i, j, k, m,p,q, l=0;
  int ny, nx;
  int occupied=0, attempts; 
  
  if (clumped){
    
    int counter=1;
    //determine block size
    while(number-counter*counter>0){
      counter++;
    }
    int iborder=(L-7*counter)/2;
    int jborder=(W-7*counter)/2;

    if(iborder<5 ||jborder<5){
      printf("CreateAgentManyCells:error: too many cells, they cannot be placed\n");
      return;
    }

    int placed=1;
    for(p=0;p<counter;p++){
      
      if(placed>number)
	break;
      
      for(q=0;q<counter;q++){
	
	if(placed>number)
	  break;

	ny=iborder+2+7*p;
	nx=jborder+2+7*q;
	
	if(!CellArray.count(placed))
	{
	  CellArray.insert(make_pair(placed, Cell()));
	  CellArray[placed].CreateCell(placed, type);
	  
	}
	else
	{
	  printf("Agent:CreateAgentManyCells: cell of this id already exists. exiting...\n");
	  exit(1);
	}
	
	for(i=-2;i<=2;i++)
	  for(j=-2;j<=2;j++){
	    CellIdGrid[ny+i][nx+j]=placed;
	    CellArray[placed].UpdateMoments(ny+i,nx+j,1);
	    CellArray[placed].vol++;
	  }


	CellArray[placed].tarvol=targetvolume;

	if(surfaceconstraint){
	  for(i=-2;i<=2;i++)
	    for(j=-2;j<=2;j++)
	      for(m=0;m<neighbourhoodsize;m++){
		if(CellIdGrid[ny+i+wideneighbourhood[m].yy][nx+j+wideneighbourhood[m].xx]!=CellIdGrid[ny+i][nx+j])
		  CellArray[placed].surf++;
	      }
	  CellArray[placed].tarsurf=targetsurface;
	}
	CellArray[placed].UpdateCell();
	placed++;
	anrcells_++;
      }
    }
  }


  else{
    //first test which cell ids have not been used yet
    int start=1;
    while(CellArray[start].tarvol){
      start++;
    }

    for(k=start;k<=number+start; k++){
      
      CellArray[k].CreateCell(k,type); //give the id, type and delay

      do{
	occupied=0;
	
	ny=(int)(uniform()*((double)(L-14))); //pick a random position, taking the unupdated border into account
	nx=(int)(uniform()*((double)(W-14)));
	ny+=5+2;
	nx+=5+2;
	
	//printf("y=%d, x=%d\n",ny,nx);
	
	for(i=-2;i<=2;i++){
	  if(ny+i>=L || ny+i<0){
	    occupied=1;
	    continue;
	  }
	  for(j=-2;j<=2;j++){
	    if(nx+j>=W || nx+j<0){
	      occupied=1;
	      continue;
	    }
	    occupied+=CellIdGrid[ny+i][nx+j];//check whether this position is empty
	    // printf("occupied:%d, counter=%d\n",occupied,l);
	    l++;
	  }
	}
	l=0;
	
	if(!occupied){ //if not, fill with new cell. Update volume, targetvolume and if necessary, (target) surface
	
	for(i=-2;i<=2;i++)
	  for(j=-2;j<=2;j++){
	    CellIdGrid[ny+i][nx+j]=k;
	    CellArray[k].UpdateMoments(ny+i,nx+j,1);
	    CellArray[k].vol++;
	  }
	  
	  CellArray[k].tarvol=targetvolume;
	
	if(surfaceconstraint){
	  for(i=-2;i<=2;i++)
	    for(j=-2;j<=2;j++)
	      for(m=0;m<neighbourhoodsize;m++){
		if(CellIdGrid[ny+i+wideneighbourhood[m].yy][nx+j+wideneighbourhood[m].xx]!=CellIdGrid[ny+i][nx+j])
		  CellArray[k].surf++;
	      }
	      CellArray[k].tarsurf=targetsurface;
	}
	CellArray[j].UpdateCell();
	
	anrcells_++;
	}
      
      else{
	attempts++;
	if (attempts > 100) {
	  printf("cannot place cell with number %d! exiting...!\n",k);
	  exit(1);
	}
      }
    }while(occupied);
    attempts=0;
    occupied=0;

    }
  }
}


void Agent::UpdateAgent_Random()//Do one monte carlo update step. Random order. Note that it is possible that not all pixels are updated, while others are updated 2x. This should not matter for dynamics.
{
  int I;
  int i,j;
  int ii,jj;
  int idA;
  int idB;
  double extravalue=0;

  for(I=0;I<L*W;I++)
    {
      //draw random site i,j, and random neighbor site
      randomCouple(&i,&j,&ii,&jj);     
      idA=CellIdGrid[i][j];
      idB=CellIdGrid[ii][jj];
      if(idA!=idB)//points belong to different cells, or cell and medium
	{
	  extravalue=Extra(i,j,ii,jj);
 
	  if(CPMAttempt(i,j,ii,jj,extravalue)==TRUE){//copy action is performed
	    DoCPMUpdate(i,j,ii,jj);
	    After(i,j,ii,jj);
	    //update ids in the plane: was first in DoCPMUpdate, but that is not very convenient: After should also have the original situation
	    CellIdGrid[i][j]=CellIdGrid[ii][jj];
	  }

	}
    }
}

int Agent::UpdateAgent_Mask()//only update small part of the field, where the tissue is.
{
  int I;
  int i,j;
  int ii,jj;
  int idA;
  int idB;
  double extravalue=0;
  int imin, imax, jmin, jmax;
  //find center of the tissue and its length: only update within those boundaries.
  if(DetermineCollectioninfo(0)==2)
  {
    return 1;
  }
  else
  {    
    imin=max(CellArray[0].meani-CellArray[0].ReturnMajorAxisLength()/2-10, 5.);
    imax=min(CellArray[0].meani+CellArray[0].ReturnMajorAxisLength()/2+10, (double)L-6.);
    jmin=max(CellArray[0].meanj-CellArray[0].ReturnMajorAxisLength()/2-10, 5.);
    jmax=min(CellArray[0].meanj+CellArray[0].ReturnMajorAxisLength()/2+10, (double)W-6.);
    
    int pixelnr=(imax-imin+1)*(jmax-jmin+1);
    
    for(I=0;I<pixelnr;I++)
    {
      //draw random site i,j, and random neighbor site
      randomMaskCouple(&i,&j,&ii,&jj, imin, imax, jmin, jmax);     
      idA=CellIdGrid[i][j];
      idB=CellIdGrid[ii][jj];
      if(idA!=idB)//points belong to different cells, or cell and medium
	{
	  extravalue=Extra(i,j,ii,jj);
 
	  if(CPMAttempt(i,j,ii,jj,extravalue)==TRUE){//copy action is performed
	    DoCPMUpdate(i,j,ii,jj);
	    After(i,j,ii,jj);
	    //update ids in the plane: was first in DoCPMUpdate, but that is not very convenient: After should also have the original situation
	    CellIdGrid[i][j]=CellIdGrid[ii][jj];
	  }

	}
    }
    return 0;
  }
}

void Agent::UpdateAgent_Walker()//deprecated; tries to do clever updates by following boundaries. may have scary side-effects.
{ 
  int I;
  int i,j;
  int NRSTEPS;
  int NRLOST;
  int lost;
  int boundary_;
  int Xi,Xj;
  int NBXi,NBXj;
  int Yi,Yj;
  int NBYi,NBYj;
  int here,there;//randomly chosen point and neighbour
  int Lx=5*(int)sqrt(L);
  int Wx=12*(int)sqrt(W);

  here=0;
  there=0;
  boundary_=FALSE;
  for(I=0;I<Lx;I++)//for(I=0;I<10;I++) 
    {
      NRSTEPS=Wx;//25;//NRSTEPS*NRLOST is max nr steps without starting with new random walker
      NRLOST=8;//64;//max nr of attempts to find new boundary point in locality of current one
      // let's find a new boundary point
      while(boundary_==FALSE) 
	{
	  randomCouple(&Xi,&Xj,&NBXi,&NBXj);// pick a global random point and neighbour
	  here=CellIdGrid[Xi][Xj];
	  there=CellIdGrid[NBXi][NBXj];
	  if(here!=there) 
	    boundary_ = TRUE;
	}
      // boundary found, now we know here and there
      // so let's travel along boundary
      lost = NRLOST;
      Yi=Xi;
      Yj=Xj;
      NBYi=NBXi;
      NBYj=NBXj;
      i = 0;
      j = NRSTEPS * NRLOST; 
      while( i != j ) 
	{
	  if( here != there ) // boundary point
	    {
	      lost = NRLOST;//start over with counting available local attempts
	      Xi=Yi;
	      Xj=Yj;
	      NBXi=NBYi;
	      NBXj=NBYj;
	      if(CPMAttempt(Xi,Xj,NBXi,NBXj,0)==TRUE) 
		{
		  DoCPMUpdate(Xi,Xj,NBXi,NBXj);
		}
	    } 
	  else // not a boundary point, countdown max number of local attempts
	    --lost;
	  if( lost == 0 ) // stop if boundary lost (maximum nr of local attempts reached)
	    j = i;
	  else // next attempt to find point on local boundary
	    {
	      randomLocalCouple(&Xi,&Xj,&Yi,&Yj,&NBYi,&NBYj);
	      here = CellIdGrid[Yi][Yj];
	      there =CellIdGrid[NBYi][NBYj];
	      i++;
	    }
	}
      //if( j != NRSTEPS * NRLOST )//find out if you left loop becaus of boundary lost or not 
      boundary_ = false;
    }
}