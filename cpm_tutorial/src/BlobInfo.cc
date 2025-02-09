#include "Header.hh"
#include "Misc.hh"
#include ""Agent.hh"


void PrintCellPositions(FILE *file, int tijd) //does what it says.
{
  int i;

  for(i=1;i<=anrcells_;i++){
    fprintf(file, "%d\t%d\t%3.2lf\t%3.2lf\t0.00\n",i, tijd,CellArray[i].meanj,CellArray[i].meani);
  }
}

void PrintCellEccentricities(FILE *file, int tijd) //eccentricity is a frequently used measure for cell elongation between 0 and 1. 1 is a 1d line, 0 is perfectly spherical. Could also just print ratio of axes.
{

  int i;

  fprintf(file,"%d\t",tijd);
  for(i=0;i<=anrcells_;i++){
    fprintf(file,"%2.3lf\t",sqrt(1-CellArray[i].ReturnShortval()/CellArray[i].ReturnLongval()));
  }
  fprintf(file,"\n");


}


//if you have a connected tissue, this function will print the direction and length of its main axes. Make sure to use DetermineCollecioninfo before!!
void PrintAxes(FILE *file, int tijd, double vecx, double vecy)
{
  double dot, veclength,angle;
  FPosition longvec;
  longvec=CellArray[0].ReturnLongVec();
  dot=vecx*longvec.xx+vecy*longvec.yy;
  veclength=sqrt(longvec.xx*longvec.xx+longvec.yy*longvec.yy);

  angle=dot/veclength;

  if(angle<-0.0001)
    angle=-angle;
  angle=acos(angle);
  angle*=180/M_PI;

  fprintf(file,"%d\t%3.3lf\t%3.3lf\t%3.3lf\t%3.3lf\t%3.3lf\n",tijd,CellArray[0].ReturnMajorAxisLength(),CellArray[0].ReturnMinorAxisLength()/CellArray[0].ReturnMajorAxisLength()/*sqrt(1-CellArray[0].ReturnShortval()/CellArray[0].ReturnLongval())*/,CellArray[0].meanj,CellArray[0].meani,angle);
 
}
