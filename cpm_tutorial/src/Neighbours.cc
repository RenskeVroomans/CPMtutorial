#include "Cell.hh"

void Cell::setNeighbour(int neighbour, int boundarylength)
{
  
  if(boundarylength==0)//remove this neighbour 
    neighbours.erase(neighbour);
  else
    neighbours[neighbour]=boundarylength; //if the element is already present, the boundarylength will be modified, otherwise a new element will be created.
    
}

int Cell::returnBoundaryLength(int cell)
{
 if(neighbours.count(cell))
   return neighbours[cell];
  
 return 0;

}

void Cell::clearNeighbours()
{
  neighbours.clear(); 
}

int Cell::updateNeighbour(int cell, int modification)
{
  if(!neighbours.count(cell) && modification<0)
  {  
    printf("Cell.updateNeighbour: error: negatively updating contact of cell %d with nonexisting neighbour %d\n",id,cell);
    return 1;
  }
  else if(!neighbours.count(cell))
    neighbours[cell]=modification;
  
  else if(neighbours.count(cell))
  {
    neighbours[cell]+=modification;
    if(neighbours[cell]==0)//remove this neighbour 
    {
      neighbours.erase(cell);
    }
    else if(neighbours[cell]<0)
    {
      neighbours.erase(cell);
      printf("Cell.updateNeighbour: error: updating contact of cell %d with neighbour %d to negative value\n",id,cell);
      return 2;
    }
  }
  return 0;
}