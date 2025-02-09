#ifndef CellHeader
#define CellHeader

#include "Header.hh"

class Cell
{
public:
  int id;
  int age;
  int celltype;
  int vol;
  int tarvol;
  int surf;
  int tarsurf;
  
   
  map<int, int>neighbours; //stores neighbouring cells(ID) and the amount of membrane contact
  
 
  double meani;
  double meanj;
 
  double tveci;
  double tvecj;
  
  //stores old pos of cells for target vector calculations
  double previ;
  double prevj;
  
  Cell();
  ~Cell();
  void CreateCell(int cid, int celltype);
  Cell& operator=(const Cell& celltocopy); //assignment operator
  Cell(const Cell &obj);//copy constructor
   
  //for cell migration
  void UpdateTarVec(void);
  void StartTarVec();
  //Cell Shape functions
  void UpdateMoments(int i, int j, int add);
  void EmptyMoments(void);
  void UpdateCellPixel(int i, int j, int add);
  void UpdateCell(void);
  double ReturnMajorAxisLength(void);
  double ReturnMinorAxisLength(void);
  FPosition ReturnChangeInAspect(int i, int j, int add); 
  FPosition DetermineMomentofInertia(int i, int j, int add);
  double ReturnLongval(void);
  double ReturnShortval(void);
  FPosition ReturnLongVec(void); 
  FPosition ReturnShortVec(void);

   //cell neighbours
  void setNeighbour(int neighbour, int boundarylength);
  void clearNeighbours();
  int returnBoundaryLength(int cell);
  int updateNeighbour(int cell, int modification);
  
protected:
  long double Sxx;
  long double Syy; 
  long double Sx; 
  long double Sxy;
  long double Sy;

  long double iyy, ixx, ixy;
  long double rhs1, rhs2;
  void UpdateInertiaTensor(void);

};



#endif
