#ifndef AgentHeader
#define AgentHeader

#include "Header.hh"
#include "Cell.hh"
#include <sys/stat.h>

class Agent
{
public:
   
  /*the load of variables which were first in Basic.cc*/
  vector< vector<int> > CellIdGrid;        //was previously CellIdGrid[L][W], but remember that L and W are no longer constants
  map <int, Cell> CellArray;
  
  int anrcells_;
  
  static int **JTable;
  //map<int, map <int, int> > JTable;
  
  /* static class parameters: the same for every agent. Definition found in Agent.cc */
  static double   labdavol;
  static double   labdasurf;
  static double   invT;
  static int      targetvolume;
  static int      targetsurface;
  static int      surfaceconstraint; //flag to see if a surfaceconstraint is used.
  static int      L;               //grid length
  static int      W;               //grid width
  static int      neighbourhoodsize; //if neighbourhoodsize is given as 0, the Moore neighbourhood is used. Necessary when using DeltaSurfHamiltonian
  static int      MaxNrCells;      //256
  static Position *wideneighbourhood;  //array holding the positions of the wider neighbourhood used for hamiltonian stuff
  static int      InitNrCells;
  static int      NrDevSteps;
  static int      picinterval;
  static int      NrCellTypes;
  static int      CellPlace;
  static int      perstime;
  static int      zoom;
  static double   mu;
  
   /************************************/
  //agent.cc
  Agent();
  ~Agent();

  void PlaceCellsInGrid(); //start up field, place cells.
  void DevelopAgent(char *desdir); //run dynamics

  void UpdateCellAge();
   
  double Extra(int i,int j, int ii, int jj); //additions to Hamiltonian
  void After(int i,int j, int ii, int jj); //bookkeeping after updating a pixel
  
  int PlaceOneCell(int type, int id, int ipos, int jpos);
  void CreateAgentManyCells(int number,int type, int clumped);
  void UpdateAgent_Random();//one MCS //double (*Extra)(int,int,int,int),void (*After)(int,int,int,int)
  void UpdateAgent_Walker();//deprecated
  int UpdateAgent_Mask();
  
  //basic.cc : to create a field and cell struct and table with Jvalues.
  void CreateField();
  void DestroyField();
  static void FindNeighbourhood();
  static int FillJTable(int nrtypes,int nrentries,...);
    
  //CellAdhesion
  double DeltaSurfHamiltonian(int I1,int J1,int I2,int J2);
  
  //CellLattice
  void ShuffleField(vector<int> &shuffle);
  void PosAndNeigh(int num, int *i, int *j, int *nbi, int *nbj);
  void randomCouple(int *i, int *j, int *nbi, int *nbj);
  void randomMaskCouple(int *i, int *j, int *nbi, int *nbj, int imin, int imax, int jmin, int jmax);
  void randomLocalCouple(int *i, int *j, int *ii,int *jj, int *nbii,int *nbjj);
  int CPMAttempt(int i,int j,int ii,int jj,double extravalue);
  void DoCPMUpdate(int i,int j,int ii,int jj);
  void InitContactLength(); //initialises contact info for cells
  void RecountContactLength(int id, int newid);//re-counts contacts for cells that divided.


  //CellDivision
  //difference between cleavage and division: targetvolume is halved in the first case.
  void DivideSingleCell(int id);
  void DivideSingleVolume(int id, int cleave);
  void CleaveSingleCell(int id);
  int DetermineCollectioninfo(int connected);//to determine overall elongation of the blob. 
  void DivideAllCells();
  void CleaveAllCells();
  void DivideAllVolumes(int cleave);
  void DivideCellOverAxis(int id, int cleave,double xvec,double yvec);
 
  //Graphics
  void Snapshot(int zoom, char *dirname); //colours CPM cells according to type
  int ColourGradient(double **field,int zoom, double maxval, char *dirname); //in case you have a float plane you want to visualise
  void printField(int time);
  
  //BlobInfo
  //functions to print info about cells in a blob
  void PrintCellPositions(FILE *file, int tijd);
  void PrintCellEccentricities(FILE *file, int tijd);
  void PrintAxes(FILE *file, int tijd, double vecx, double vecy);
    
      
private:
  int piccounter;
  int piccounter2;
};

#endif
