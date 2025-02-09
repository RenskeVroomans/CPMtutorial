#include <cstdarg>
#include "Agent.hh"
#include "Cell.hh"



void Agent::CreateField()
{
  int i,j;
//   if (CellIdGrid!=NULL){
//     printf("CreateField:warning: memory for CellIdGrid has already been allocated.\n");
//     return;
//   }
//   CellIdGrid=(int **)calloc((size_t)L,sizeof(int *));
//   if (CellIdGrid==NULL){
//     fprintf(stderr,"error in memory allocation CellIdGrid\n");
//     exit(1);
//   }
// 
//   CellIdGrid[0]=(int *)calloc((size_t)(L*W),sizeof(int));
// 
//   if (CellIdGrid[0]==NULL){
//     fprintf(stderr,"error in memory allocation CellIdGrid\n");
//     exit(1);
//   }
//   for(i=1,j=L; i<j; i++){
//     CellIdGrid[i]=CellIdGrid[i-1]+W;
//   }
  
  
  for(i=0; i<L; i++)
  {
    CellIdGrid.push_back(vector<int>());
    
    for(j=0; j<W; j++)
    {
      CellIdGrid[i].push_back(0);
    }
  }

}

void Agent::DestroyField()
{
  //free(*CellIdGrid);
  //free(CellIdGrid);
  
}



void Agent::FindNeighbourhood()
{
  int i,j,k;

  int dist=1;
  int radius=1;
  int number=0;
  int found=0;
  int level=0;
  int border=0;
  
  for(k=1; level<neighbourhoodsize; k++){
    for(i=-dist; i<=dist;i++){
      for(j=-dist; j<=dist;j++){
	if(i*i+j*j==radius){
	  found=1;
	  number++;
	}
	if(found && i==dist)
	  border=1;
      }
    }
    if(found){
      if(border)
	dist++;
      level++;
      found=0;
      border=0;
    }
    radius++;
  }
  
  wideneighbourhood=(Position *)calloc((size_t) number,sizeof(Position));
  if(wideneighbourhood==NULL){
    printf("FindNeighbourhood: could not allocate memory for wideneighbourhood array. Using Moore neighbourhood\n");
    return;
  }
  dist=1;
  radius=1;
  number=0;
  found=0;
  level=0;
  border=0;
  
  for(k=1; level<neighbourhoodsize; k++){
    for(i=-dist; i<=dist;i++){
      for(j=-dist; j<=dist;j++){
	if(i*i+j*j==radius){
	  found=1;
	  wideneighbourhood[number].xx=j;
	  wideneighbourhood[number].yy=i;
	  number++;
	}
	if(found && i==dist)
	  border=1;
      }
    }
    if(found){
      if(border)
	dist++;
      level++;
      found=0;
      border=0;
    }
    radius++;
  }
  
  neighbourhoodsize=number;

}

/*function which allocates and fills the JTable. Assumes the variable arguments come in the format int type1, int type2, int value.*/
int Agent::FillJTable(int nrtypes, int nrentries, ...)
{
  int i,j;
  /* allocate memory for the JTable */
  if (JTable!=NULL){
    printf("AllocateJTable:warning: memory for JTable has already been allocated.\n");
  }
  else{
    JTable=(int **)calloc((size_t)nrtypes,sizeof(int *));
    if (JTable==NULL){
      fprintf(stderr,"error in memory allocation JTable\n");
      exit(1);
    }
    
    JTable[0]=(int *)calloc((size_t)(nrtypes*nrtypes),sizeof(int));
    
    if (JTable[0]==NULL){
      fprintf(stderr,"error in memory allocation JTable\n");
      exit(1);
    }
    for(i=1,j=nrtypes; i<j; i++){
      JTable[i]=JTable[i-1]+nrtypes;
    }
  }

  /* in this section, the JTable is filled*/
  va_list J_list; //list of given arguments
  int type1, type2, value;

  va_start(J_list,nrentries*3);

  /* if(((nrentries+1)%3)!=0){
    printf("FillJtable:incorrect argument list. Please check and try again\n");
    return 1;
    }*/
  /*the number of entries required to fill the table depends on the nr of celltypes and follows the triangular numbers (-1 for the medium/medium). The equation to obtain the right number of entries is 0.5*n(n+1), where n is the nr of celltypes.*/
 
  int reqnr=(int)(0.5*nrtypes*(nrtypes+1));
  
  if(nrentries!=reqnr && nrentries!=reqnr-1){
    printf("FillJtable:incorrect size of argument list. Please check and try again\n");
    return 1;
  }


  for(i=0;i<(int)(nrentries);i++){
    type1=va_arg(J_list,int);
    type2=va_arg(J_list,int);
    value=va_arg(J_list,int);
   
    //   printf("FillJtable:incorrect argument list. Please check and try again\n");
    //return 1;
   
    JTable[type1][type2]=JTable[type2][type1]=value;
  }

  if(JTable[0][0]!=0){
    printf("FillJtable:nonzero J value for Medium/Medium given. Set to zero.\n");
    JTable[0][0]=0;
  }

  va_end(J_list);

  return 0;
   

}

