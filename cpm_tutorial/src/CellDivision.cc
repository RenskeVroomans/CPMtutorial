#include "Agent.hh"
#include "Misc.hh"



  //PCA 
  //http://en.wikipedia.org/wiki/Principal_component_analysis
  //5 Algorithm #1: the covariance method

  // 5.1 Organize the data set
  // data are set of N vectors representing N observations of the M variables.
  // write as M * N matrix M=2(i coord, j coord) N=nr of points in cell
  // ->goes in matrix AA

  // 5.2 Calculate the empirical mean row vector
  // u[m] = 1/N \sum_{n=1}^N X[m,n] : mean for each of M variables: mean i, mean j
  // ->goes in matrix BB


  // 5.3 Calculate the deviations from the mean
  // Subtract row vector u[m] from each column of  matrix X.
  // Store mean-subtracted data in the M × N matrix B.
  // B = X - u[m] h (h vector with all 1-s) 
  // -> difference of AA and BB, result goes in matrix AA again: gsl_matrix_sub

  // 5.4 Find the covariance matrix
  // Find the M × M empirical covariance matrix C from the outer product of matrix B with itself:
  // C = Exp [ B outerproduct B]
  //   = Exp [B product Bconjugatetranspose]
  //   = 1/N B product Bconjugatetranspose
  //   = 1/N B product B transpose
  // Transpose is just turning rows to columns and columns to rows
  // Matrix C will just be 2 by 2 matrix
  // -> matrix CC is transpose of matrix AA: gsl_matrix_transpose_memcpy
  // -> product of AA times CC gives matrix DD: gsl_blas_dgemm
  // -> DD is divided by N, result again in DD: gsl_matrix_scale


  // 5.5 Find the eigenvectors and eigenvalues of the covariance matrix C
  // Remember labda*labda-(a+d)labda+(ad-bc)=0 (characteristic polynomial, a,b,c,d entries of C)
  // So labda=((a+d)+/-sqrt((a+d)*(a+d)-4*(ad-bc)))/2*(a+d)
  // (a-labda)*x+b*y=0 en c*x+(d-labda)*y=0, fill in found eigenvalue gives eigenvector
  // -> :gsl_eigen_symmv
  // The eigenvectors are guaranteed to be mutually orthogonal and normalised to unit magnitude. 

  // 5.6 Rearrange the eigenvectors and eigenvalues
  // -> :gsl_eigen_symmv_sort: sort them on magnitude and in descending order




//this function is a derivative of the functions in Cell.cc: it determines the principle components of the collection of cells. It sort of assumes they form a coherent blob. 
//It updates the Moments in the memory space for the medium (cell 0). In this particular case, these moments need to be emptied first; 
//they will not be updated in the Hamiltonian itself.
int Agent::DetermineCollectioninfo(int connected)
{
  int i,j;
  int totalvol=0;
  //we use the elements of Medium here to do the calculations
  for(i=1;i<=anrcells_;i++){
    totalvol+=CellArray[i].vol;
  }

  CellArray[0].EmptyMoments();
  CellArray[0].vol=0;
 //in matrix A column 0 i coordinates, column 1 j coordinates of the real points of the cell
  for(i=0;i<L;i++)
  {
    for(j=0;j<W;j++)
    {
      // if connected is 1, a cell only counts when it's connected to tissue. Assume that when a cell is connected it has more than 1 neigh (at least medium + 1 cell)
      if(CellIdGrid[i][j] && (!connected || CellArray[CellIdGrid[i][j]].neighbours.size()>1) ) 
      {
	CellArray[0].UpdateMoments(i,j,1);
	CellArray[0].vol++;
      }
    }
  }  
  
  CellArray[0].UpdateCell();
  
  if(totalvol!=CellArray[0].vol && !connected)
  {
    printf("DetermineCollectionInfo: error: added volume of cell is not equal to counted volume of non-zero spins\n");
    return 1;
  }
  if(CellArray[0].vol==0)
  {
    //printf("DetermineCollectionInfo: error:no volume counted!\n"); //don't sqeak, just pass info on to kill off the agent
    return 2;
  }
  
  return 0;
}

void Agent::DivideSingleVolume(int id, int cleave)
{
  int i;
  int j;
  int k;
  double aa;
  double bb;
  int newid;
  double ran;
  FPosition shortvec;

  if(CellArray[id].vol>1 && CellArray[id].tarvol) //cell exists and has a target volume -> not in apoptosis
  {
    newid=anrcells_+1;
    CellArray[newid]=CellArray[id];
    CellArray[newid].id=newid;
    anrcells_++;
        
    shortvec=CellArray[id].ReturnShortVec();
    if((int)(1000*shortvec.xx)!=0) //used newid previously, which is silly
	aa=(shortvec.yy/shortvec.xx);
    else
      aa=1000;
    bb=CellArray[id].meani-CellArray[id].meanj*aa;//was: CellArray[id].meanj-CellArray[id].meani*aa +0.5
    
    ran=uniform();
    for(i=0;i<L;i++)
      for(j=0;j<W;j++)
      {
	
	if(CellIdGrid[i][j]==id)
	{
	  // if(ran<0.5)//randomise which half of the cell gets the new id. may not be necessary
	  // {
	    if (i>((int)(bb+aa*(double)j))) //this is the new cell; was j>(bb+aa*(double)i)
	     { 
	       CellArray[id].vol--;
	       CellArray[id].UpdateMoments(i,j,-1);
	       CellIdGrid[i][j]=newid;
	       if(cleave){ //when we cleave the cell should not grow again. so we update the tarvol to  its actual vol.
		 CellArray[id].tarvol--;
	       }
	     }
	     else //this is the old cell: remove data from new cell
	     {
	       CellArray[newid].vol--;
	       CellArray[newid].UpdateMoments(i,j,-1);
	       if(cleave)
		 CellArray[newid].tarvol--;
	     }
	}
      }
      
      //InitContactLength();
      RecountContactLength(id, newid);
      CellArray[id].UpdateCell();
      CellArray[newid].UpdateCell();
      
      
      if(neighbourhoodsize && surfaceconstraint){
	//need to recompute the surface and targetsurface
	CellArray[id].surf=0;
	CellArray[newid].surf=0;
	if(cleave){
	  CellArray[newid].tarsurf=0;
	  CellArray[id].tarsurf=0;
	}
	else
	  CellArray[newid].tarsurf=CellArray[id].tarsurf;
	
	for(i=1;i<L-1;i++)
	  for(j=1;j<W-1;j++)
	  {
	    if(CellIdGrid[i][j]==id || CellIdGrid[i][j]==newid)//just split cells there?
		{
		  for(k=0; k<neighbourhoodsize; k++){
		    {
		      if(CellIdGrid[i+wideneighbourhood[k].yy][j+wideneighbourhood[k].xx]!=CellIdGrid[i][j])//nr of neighbour points not part of same cell equals nr of boundary surfaces
			{
			  CellArray[CellIdGrid[i][j]].surf++;
			  if(cleave)
			    CellArray[CellIdGrid[i][j]].tarsurf++;
			}
		    }
		  }
		}
	  }
      }
  }
}

///may do it completely wrong!!
void Agent::DivideCellOverAxis(int id, int cleave,double xvec,double yvec)
{
  int i;
  int j;
  int k;
  double aa;
  double bb;
  int newid;
  double ran;
  FPosition shortvec;

  if(CellArray[id].vol>1  && CellArray[id].tarvol)
  {
    newid=anrcells_+1;
    CellArray[newid]=CellArray[id];
    CellArray[newid].id=newid;
    anrcells_++;
        
    if((int)(1000*xvec)!=0) //used newid previously, which is silly
	aa=(yvec/xvec);
    else
      aa=1000;
    bb=CellArray[id].meanj-CellArray[id].meani*aa+0.5;
    
    ran=uniform();
    for(i=0;i<L;i++)
      for(j=0;j<W;j++)
      {
	
	if(CellIdGrid[i][j]==id)
	{
	  // if(ran<0.5)//randomise which half of the cell gets the new id. may not be necessary
	  // {
	    if (j>((int)(bb+aa*(double)i))) //this is the new cell; remove data from old cell
	     { 
	       CellArray[id].vol--;
	       CellArray[id].UpdateMoments(i,j,-1);
	       CellIdGrid[i][j]=newid;
	       if(cleave){ //when we cleave the cell should not grow again. so we update the tarvol to  its actual vol.
		 CellArray[id].tarvol--;
	       }
	     }
	     else //this is the old cell: remove data from new cell
	     {
	       CellArray[newid].vol--;
	       CellArray[newid].UpdateMoments(i,j,-1);
	       if(cleave)
		 CellArray[newid].tarvol--;
	     }
	}
      }
      
      
      RecountContactLength(id, newid);
      
      CellArray[id].UpdateCell();
      CellArray[newid].UpdateCell();
      
      
      if(neighbourhoodsize && surfaceconstraint){
	//need to recompute the surface and targetsurface
	CellArray[id].surf=0;
	CellArray[newid].surf=0;
	if(cleave){
	  CellArray[newid].tarsurf=0;
	  CellArray[id].tarsurf=0;
	}
	else
	  CellArray[newid].tarsurf=CellArray[id].tarsurf;
	
	for(i=1;i<L-1;i++)
	  for(j=1;j<W-1;j++)
	  {
	    if(CellIdGrid[i][j]==id || CellIdGrid[i][j]==newid)//just split cells there?
		{
		  for(k=0; k<neighbourhoodsize; k++){
		    {
		      if(CellIdGrid[i+wideneighbourhood[k].yy][j+wideneighbourhood[k].xx]!=CellIdGrid[i][j])//nr of neighbour points not part of same cell equals nr of boundary surfaces
			{
			  CellArray[CellIdGrid[i][j]].surf++;
			  if(cleave)
			    CellArray[CellIdGrid[i][j]].tarsurf++;
			}
		    }
		  }
		}
	  }
      }
    }
}

void Agent::DivideAllVolumes(int cleave)
{
  int i,j,k;
  int I;
  map<int, Cell>::iterator it;
  int maxid=anrcells_;
  int newid;
  double aa;
  double bb;
  double ran;
  FPosition shortvec;

  for(it=CellArray.begin(); it!=CellArray.end(); ++it)
  {
    if(it->second.vol>1 && it->second.tarvol && it->second.id>0 && it->second.id<=maxid)//don't divide the medium, don't divide newer cells.
    {
      shortvec=it->second.ReturnShortVec();
      cout << "Short axis: xx="<<shortvec.xx<<", yy="<<shortvec.yy<<endl;
      newid=anrcells_+1;
      CellArray[newid]=it->second;//copy the old cell with all its data into the new
      CellArray[newid].id=newid;
      anrcells_++;
                  
      if((int)(1000*shortvec.xx)!=0)
	aa=shortvec.yy/shortvec.xx;
      else
	aa=1000;
      bb=it->second.meani-it->second.meanj*aa;
      
      
      //printf("cell %i meani %f meanj %f aa %f bb %f\n",I,CellArray[I].meani,CellArray[I].meanj,aa,bb);
      //ran=uniform();
      for(i=0;i<L;i++)
	for(j=0;j<W;j++)
	{
	  
	  if(CellIdGrid[i][j]==it->second.id)
	  {
	   // if(ran<0.5)//randomise which half of the cell gets the new id. may not be necessary
	   // {
	     if (i>((int)(bb+aa*(double)j))) //this is the new cell; remove data from old cell
	     { 
	       it->second.vol--;
	       it->second.UpdateMoments(i,j,-1);
	       CellIdGrid[i][j]=newid;
	       if(cleave){ //when we cleave the cell should not grow again. so we update the tarvol to  its actual vol.
		 it->second.tarvol--;
	       }
	     }
	     else //this is the old cell: remove data from new cell
	     {
	       CellArray[newid].vol--;
	       CellArray[newid].UpdateMoments(i,j,-1);
	       if(cleave)
		 CellArray[newid].tarvol--;
	     }
	    }
	}
	
	
	RecountContactLength(it->second.id, newid);//if not working with neighbours, disable this! costly.
	it->second.UpdateCell();
	CellArray[newid].UpdateCell();
	
	//need to recompute the surface and targetsurface
	//can be done smarter: surface for one is usually also surface for other
	if(neighbourhoodsize &&surfaceconstraint){
	
	it->second.surf=0;
	CellArray[newid].surf=0;
	if(cleave){
	  it->second.tarsurf=0;
	  CellArray[newid].tarsurf=0;
	}
	else
	  CellArray[newid].tarsurf=it->second.tarsurf;
      }
    }
  }    
	
  for(i=1;i<L-1;i++)
    for(j=1;j<W-1;j++)
      {
	if(CellIdGrid[i][j]!=0)//cell there?
	  {
	    for(k=0; k<neighbourhoodsize; k++){
	      
	      if(CellIdGrid[i+wideneighbourhood[k].yy][j+wideneighbourhood[k].xx]!=CellIdGrid[i][j])//nr of neighbour points not part of same cell equals nr of boundary surfaces
		{
		  
		  CellArray[CellIdGrid[i][j]].surf++;
		  if(cleave)
		    CellArray[CellIdGrid[i][j]].tarsurf++;
		}
	    }
	  }
      }
	

//   for(I=0;I<=anrcells_;I++)
//     {
//       printf("cellid %i cellvol %i celltarvol %i cellsurf %i celltarsurf %i\n",I,CellArray[I].vol,CellArray[I].tarvol,CellArray[I].surf,CellArray[I].tarsurf);
//     }
  

}

void Agent::DivideSingleCell(int id)
{
   DivideSingleVolume(id,0);
}

void Agent::CleaveSingleCell(int id)
{
  DivideSingleVolume(id,1);
}

void Agent::DivideAllCells()
{
  DivideAllVolumes(0);
}

void Agent::CleaveAllCells()
{
  DivideAllVolumes(1);
}
