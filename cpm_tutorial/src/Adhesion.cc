/* this program incorporates (and can combine) two mechanisms of convergent extension: blokjes and targetvector. If so desired, a growth zone is implemented from a third celltype, which divides and differentiates.*/

#include "Header.hh"
#include "Agent.hh"
#include "Misc.hh"
#include <sys/stat.h>

#define MEDIUM     0
#define CELL       1
#define CELLTYPES  2

#define GENE       0
#define ADHESION   1

/**********************************************************/
/*********** externally defined variables *****************/
/**********************************************************/

extern int **CellIdGrid;
extern Cell *CellArray;
extern int anrcells_;
extern double labdavol;
extern double labdasurf;
extern  double invT;
extern  int L;//grid length=maximum length of animal
extern  int W;//grid width=maximum width of animal
extern  int MaxNrCells;
extern int firstvol;
extern int seed;
extern int neighbourhoodsize;
extern double yieldzero;
extern double alpha;
extern int surfaceconstraint; 
extern int targetvolume, targetsurface;
extern int **JTable;
extern Position *wideneighbourhood;
extern int genenr;
extern double timestep;


/**********************************************************/
/**************** program variables ***********************/
/**********************************************************/

Agent *A;

struct stat sb;        //allows to check for the existence of directories
int order;             //neighbourhoodorder
dsfmt_t dsfmt;
char dirname[2000];
char subdirname[2000];
char picturedirname[2000];


/* files to write to */
FILE *celldistdata;     //contains motilitydata of cells and their distance to the center/longaxis
FILE *cellposdata;      //position of cells at each time point
FILE *fielddata;
FILE *genedata;

int    **boundarylength;     //to store the length of the boundary between two cells
int    **sticky;

/* to read from a parameter file */
int nrsimulations;     // the number of simulations to run
int duration;          //number of timesteps to take
int initialization;    //time for cells to blow up and come to first quasi steady state
int pictureinterval;   //interval at which to take pictures.
int printinfo;          //whether to print data to files
double divisionrate;    //rate of division in the growth zone
int cellnumber;         // number of cells to start with
int stiffcells;         //yieldvalue of cells in first row
int oscillationrate;    //period of switching between celltypes. when it is zero, assume there is no growth zone

int stripelength;    //how long a stripe should be (nr of cells)
int stripewidth;     //how many cell layers in the stripe
int nrstripes;       //nr of stripes when not doing growthzone
int staggered;       // flag whether to initialize rows in staggered configuration
int initcellwidth;   //width of initial cell square

int Jcell1med;
int Jcell1cell1;

double mm, hh;

/**********************************************************/
/******************* functions ****************************/
/**********************************************************/

void ReadParameters(char *file)
{
  int checkpars=0;

  if(ReadHeaderVars(file)){
    printf("there is a problem with reading the parameters. exiting...\n");
    exit(1);
  }
  checkpars+=ReadData(file,'i',"neighbourhoodorder",&order);
  checkpars+=ReadData(file,'i',"yield",&yieldzero);
  checkpars+=ReadData(file,'d',"alpha",&alpha);
  checkpars+=ReadData(file,'d',"labdasurf",&labdasurf);
  checkpars+=ReadData(file,'i',"targetvolume",&targetvolume);
  checkpars+=ReadData(file,'i',"targetsurface",&targetsurface);
  checkpars+=ReadData(file,'i',"surfaceconstraint",&surfaceconstraint);
  checkpars+=ReadData(file,'i',"pictureinterval",&pictureinterval);
  checkpars+=ReadData(file,'i',"duration",&duration);
  checkpars+=ReadData(file,'i',"nrsimulations",&nrsimulations);
  checkpars+=ReadData(file,'i',"initialization",&initialization);
  checkpars+=ReadData(file,'i',"cellnumber",&cellnumber);

  checkpars+=ReadData(file,'i',"printinfo",&printinfo);
  
  checkpars+=ReadData(file,'i',"stiffcells",&stiffcells);
  checkpars+=ReadData(file,'i',"oscillationrate",&oscillationrate);
  //checkpars+=ReadData(file,'d',"postdist",&postdist);
  checkpars+=ReadData(file,'d',"divisionrate",&divisionrate);
 
  checkpars+=ReadData(file,'i',"stripelength",&stripelength);
  checkpars+=ReadData(file,'i',"stripewidth",&stripewidth);
  checkpars+=ReadData(file,'i',"nrstripes",&nrstripes);
  checkpars+=ReadData(file,'i',"staggered",&staggered);
  checkpars+=ReadData(file,'i',"initcellwidth",&initcellwidth);

  checkpars+=ReadData(file,'i',"Jcell1med",&Jcell1med); 
  checkpars+=ReadData(file,'i',"Jcell1cell1",&Jcell1cell1);

  checkpars+=ReadData(file,'d',"mm",&mm);
  checkpars+=ReadData(file,'d',"hh",&hh);
 
  genenr=ADHESION;
 
  if (checkpars){
    exit(1);
  }

}


void updateBoundaries(int copied,int copying,int posi, int posj)//ids of cell copied into and copying, and central pos
{
  int neigh;
  for(int ii=-1;ii<=1;ii++)
	for(int jj=-1;jj<=1;jj++)
	  {
	    if((ii!=0||jj!=0) /*&& posi+ii>=0 && posi+ii<L && posj+jj>=0 && posj+jj<W*/ )
	      {
		neigh=CellIdGrid[posi+ii][posj+jj];
		if(neigh!=copied){
		  boundarylength[copied][neigh]-=1;
		  boundarylength[neigh][copied]-=1;
		}
		if(neigh!=copying){
		  boundarylength[copying][neigh]+=1;
		  boundarylength[neigh][copying]+=1;
		  
		}
	      }
	  }
  
  
  //  boundarylength[copied][copying]+=(after-before);
  //boundarylength[copying][copied]+=(after-before);
    
}

void checkBoundaryUpdate(void)
{
  const int size=anrcells_+1;
  int boundaries[size][size];
  int cellA, cellB;
  for(int i=0;i<anrcells_+1;i++)
    for(int j=0;j<anrcells_+1;j++)
      boundaries[i][j]=0;
  
  for(int i=5;i<L-5;i++)
    for(int j=5;j<W-5;j++)      
      if((cellA=CellIdGrid[i][j]))
	{
	  for(int ii=-1;ii<=1;ii++)
	    for(int jj=-1;jj<=1;jj++)//go through neighbours of cell
	      {
		cellB=CellIdGrid[i+ii][j+jj];
		if(cellA!=cellB && cellB)//automatically excludes counting the middle pixel
		  {
		    boundaries[cellA][cellB]++;
		  }
	      }
	}

  //check matrix
  for(int i=1;i<anrcells_+1;i++)
    for(int j=1;j<anrcells_+1;j++)
      if(boundarylength[i][j]!=boundaries[i][j])
	{
	  printf("Boundary check:measured boundary differs from stored: Cell A=%d, CellB=%d. Difference: %d\n",i,j,boundarylength[i][j]-boundaries[i][j]);
	  // exit(1);
	}
  
}

void initiateBoundaries(void)
{
  int cellA,cellB;

  //make sure array is fully 0
  for(int i=0;i<anrcells_+1;i++)
    for(int j=0;j<anrcells_+1;j++)
       boundarylength[i][j]=0;
  
  for(int i=5;i<L-5;i++)
    for(int j=5;j<W-5;j++)
      if((cellA=CellIdGrid[i][j]))
	{
	  for(int ii=-1;ii<=1;ii++)
	    for(int jj=-1;jj<=1;jj++)//go through neighbours of cell
	      {
		cellB=CellIdGrid[i+ii][j+jj];
		if((cellA!=cellB) && cellB)//automatically excludes counting the middle pixel
		  {
		    boundarylength[cellA][cellB]++;
		  }
	      }
	}
  
  //check matrix symmetry
  for(int i=1;i<anrcells_+1;i++)
    for(int j=1;j<anrcells_+1;j++)
      if(boundarylength[i][j]!=boundarylength[j][i])
	{
	  printf("InitiateBoundaries: error in initiating array. Non-symmetric...\n");
	  exit(1);
	}
  
}

void ClearMem(void)
{
  
  free(*CellIdGrid);
  free(CellIdGrid);
  CellIdGrid=NULL;
 
     
  free(CellArray);
  CellArray=NULL;
  
  delete(A);
  if(printinfo)
    CleanBlobVariables();

}

void EmptyVariables(void)
{
  anrcells_=0;
  
  for(int i=0;i<anrcells_+1;i++)
    for(int j=0;j<anrcells_+1;j++){
      boundarylength[i][j]=0;
      sticky[i][j]=0;
	
    }
}

int calculateAdhesion(int A, int B, int gene)
{
  double adhesion=0.;
  double adA,adB, ad;

  adA=CellArray[A].oldstates.at(gene);
  adB=CellArray[B].oldstates.at(gene);

  ad=((adA)<(adB) ? (adA) : (adB)); //take minimum of two values
  // printf("adA: %.2lf, adB: %.2lf, picked %.2lf\n",adA,adB,ad);
 
  adhesion=mm*ad/(ad+hh); //make m and h hill functions when multiple adhesion genes!
 

  return (int)adhesion;
}

void fillStickyArray(void)
{

  for(int i=1; i<=anrcells_;i++){
    for(int j=1; j<=anrcells_; j++){
      sticky[i][j]=calculateAdhesion(i,j,GENE);
      printf("i:%d, j:%d, sticky:%d\n",i,j,sticky[i][j]);
    }
    //  printf("\n");
  }
}


double Extra(int i,int j, int ii, int jj)//Extra function which contains all mechanisms (so they can be combined!)
{
  int idA, idB, idN;
  double DYield=0, DJ=0;
  double temp1,yield;
 
  idA=CellIdGrid[i][j];
  idB=CellIdGrid[ii][jj];

    
  if((idA && idA<=stripelength)||(idB && idB<=stripelength))
    DYield=stiffcells;

  //compute extra DeltaJ due to adhesion proteins. Note that also the yield needs to be adapted!

 

  for(int NB=0;NB<neighbourhoodsize;NB++)
    {    
      idN=CellIdGrid[i+wideneighbourhood[NB].yy][j+wideneighbourhood[NB].xx];
      if(idN!=idA)//now inside same cell, after change in different cell
	{
	  temp1-=sticky[idA][idN];
	  yield-=temp1;
	}
      if(idN!=idB)//now inside different cell, after change in same cell
	{
	  temp1+=sticky[idB][idN];
	  yield+=temp1;
	}
     
      DJ+=temp1;
      temp1=0;
    }

  return -DJ;
}

void After(int i,int j, int ii, int jj)//always use this when using Extra5!
{
  int idA, idB;
  idA=CellIdGrid[i][j];
  idB=CellIdGrid[ii][jj];


  if(idA){
    if(CellArray[idA].vol<1)
      printf("cell %d has no vol!\n",idA);
   CellArray[idA].UpdateCellPixel(i,j,-1);
  }
  if(idB){
    if(CellArray[idB].vol<1)
      printf("cell %d has no vol!\n",idB);
    CellArray[idB].UpdateCellPixel(i,j,1);
  }
  // UpdateBoundaries(idA,idB,i,j);
}

/* this function places 3 stripes, each of different celltype */
void PlaceStartStripes(int broadness,int stripelength)//width and length (in cells) of stripes
{
  // int stripelength=20; //nr of cells over the length of a stripe
  int i, j, ii, jj, m;
  int cellcounter=1;
  int type=1;
  
  CreateField(); //allocate memory for the field
  AllocateCellArray(); //allocate memory for the cellstruct
  
  if(CellIdGrid==NULL){
    printf("CreateAgentManyCells:error: field not yet made.\n");
    return;
  }
  
  if(CellArray==NULL){
    printf("CreateAgentManyCells:error: CellArray not yet made.\n");
    return;
  }

  while(W<stripelength*7){
    printf("PlaceStripes:error:going to place more cells than the field is wide. decreasing number of cells in width\n");
    stripelength--;
  }
 
  int startwidth=(W-7*stripelength)/2;

  while(3*7*broadness>L-20){
    printf("PlaceStripes:warning:stripes too big. width reduced to %d\n", broadness-1);
    broadness--;
  }
  if(broadness*3*stripelength>MaxNrCells){
    printf("PlaceStripes:error:going to place more cells than allocated. Increase MaxNrCells or reduce stripes\n");
    exit(1);
  }

  // int startlength=(L-7*broadness*number)/2;
  //we're starting at the border now: startlength is just that of the border
  int startlength=5;
  for(int p=0;p<4*broadness;p++){
   
    for(int q=0; q<stripelength; q++){

      ii=startlength+2+p*7;
      jj=startwidth+2+q*7;

     

      if(!(p%broadness)&& q==0 && type!=3){//switch celltype, except when we have celltype three: want that domain to be bigger.
	type=(int)(p/broadness)%3+1;
      }

      CellArray[cellcounter].CreateCell(cellcounter,type); //give the id
     
      for(i=-2;i<=2;i++)
	for(j=-2;j<=2;j++){
	  CellIdGrid[ii+i][jj+j]=cellcounter;
	  CellArray[cellcounter].UpdateMoments(ii+i,jj+j,1);
	  CellArray[cellcounter].vol++;
	}


	CellArray[cellcounter].tarvol=targetvolume;

	if(surfaceconstraint){
	  for(i=-2;i<=2;i++)
	    for(j=-2;j<=2;j++)
	      for(m=0;m<neighbourhoodsize;m++){
		if(CellIdGrid[ii+i+wideneighbourhood[m].yy][jj+j+wideneighbourhood[m].xx]!=CellIdGrid[ii+i][jj+j])
		  CellArray[cellcounter].surf++;
	      }
	  CellArray[cellcounter].tarsurf=targetsurface;
	}
	CellArray[cellcounter].UpdateCell();
	cellcounter++;
	anrcells_++;
    }
  }

  
}


void Start(int argc, char **argv)
{
  char command[]="mkdir ";
  char makedirname[50];
  int i,j;

  if(argc!=3){
    printf("usage: %s <dirname> <parfile>\n",argv[0]);
    exit(1);
  }
  else{
    strcpy(dirname,argv[1]);
  }

  /*** make directory if necessary ***/
  if (stat(dirname, &sb) == 0 && S_ISDIR(sb.st_mode)){//directory exists, just check if the cells directory exists.
    //do nothing; used to check the cells dir here
  }
    
  else if(stat(dirname, &sb) == 0 && !S_ISDIR(sb.st_mode)){//file exists but is not a directory
    printf("given file is not a directory. give another name\n");
    exit(1);
  }
  else{ //directory does not exist: make it.
    strcpy(makedirname,command);
    strcat(makedirname,dirname);
    if(system(makedirname)==-1){ //make directory for data
      printf("warning: could not make directory %s. Exiting...\n",dirname);
      exit(1);
    }
  }
  /***********************************/
  /***read the parameters ***/
  ReadParameters(argv[2]);
  CopyPars(argv[2],dirname);
 /*************************/
  if(MaxNrCells<cellnumber){
    printf("Error: too many cells, maxcellnr is too small. exiting..\n");
    exit(1);
  }

  /**allocate memory for arrays**/
  boundarylength=(int **)calloc((size_t)MaxNrCells,sizeof(int *));
  if (boundarylength==NULL){
    fprintf(stderr,"error in memory allocation boundarylength\n");
    exit(1);
  }
  
  boundarylength[0]=(int *)calloc((size_t)(MaxNrCells*MaxNrCells),sizeof(int));
  
  if (boundarylength[0]==NULL){
    fprintf(stderr,"error in memory allocation boundarylength \n");
    exit(1);
  }
  for(i=1,j=MaxNrCells; i<j; i++){
    boundarylength[i]=boundarylength[i-1]+MaxNrCells;
  }
 
  sticky=(int **)calloc((size_t)MaxNrCells,sizeof(int *));
  if (sticky==NULL){
    fprintf(stderr,"error in memory allocation sticky\n");
    exit(1);
  }
  
  sticky[0]=(int *)calloc((size_t)(MaxNrCells*MaxNrCells),sizeof(int));
  
  if (sticky[0]==NULL){
    fprintf(stderr,"error in memory allocation sticky \n");
    exit(1);
  }
  for(i=1,j=MaxNrCells; i<j; i++){
    sticky[i]=sticky[i-1]+MaxNrCells;
  }
 
   FindNeighbourhood(order);
 
  printf("neighbourhoodsize:%d\n",neighbourhoodsize);

  
  if(FillJTable(CELLTYPES,2,
		0,1,Jcell1med,
		1,1,Jcell1cell1
		))
    exit(1);
  
    
  int max=0;
  //determine a good value for yieldzero
  for(i=0;i<CELLTYPES; i++){
    for(j=0;j<CELLTYPES;j++){
      printf("%d\t",JTable[i][j]);
      if(JTable[i][j]>max){
	max=JTable[i][j];
      }
    }
    printf("\n");
  }
  yieldzero=0;//18*10;//max*neighbourhoodsize*2*alpha;

  
  /***********************/
 
}

void reverse(char s[])
{
  int c, i, j;
  for(i=0, j=strlen(s)-1; i<j; i++,j--){
    c=s[i];
    s[i]=s[j];
    s[j]=c;
  }
}

void SubStart(int run)
{
  int i,n;
  char rundir[50];
  char command[]="mkdir ";
  char makedirname[50];
  char a_string[50];
  char pictures[]="/Cells";

  //convert runnumber into char
  i=0;
  n=run;
  do{
    rundir[i++]=n%10+'0';
  }while((n/=10)>0);
  //rundir[i++]='/';
  rundir[i]='\0';
  reverse(rundir);
  
  //make the subdir string
  strcpy(subdirname,dirname);
  strcat(subdirname,"/");
  strcat(subdirname,rundir);
  strcpy(picturedirname,dirname); 
  strcat(picturedirname,pictures);
  strcat(picturedirname,rundir);

  
  /* make the picture directory if it does not exist */
  if (stat(picturedirname, &sb) == 0 && S_ISDIR(sb.st_mode)){//picturedirectory exists
    printf("directory %s exists\n",picturedirname);
  }
  else if(stat(picturedirname, &sb) == 0 && !S_ISDIR(sb.st_mode)){//file exists but is not a directory
    printf("in directory there already is another file called \"Cells%s\". give another name\n",rundir);
    exit(1);
  }

  else{ //make picturedir
    strcpy(makedirname,command);
    //strcat(makedirname,dirname);
    strcat(makedirname,picturedirname);
    printf("making directory as follows: %s\n",makedirname);
    if(system(makedirname)==-1){ //make subdirectory for pictures
      printf("warning: could not make directory %s. Exiting...\n",picturedirname);
      exit(1);
    }
  }
  strcpy(makedirname,"");//empty the string
  
  //create the subdirectory if necessary
  if (stat(subdirname, &sb) == 0 && S_ISDIR(sb.st_mode)){//directory exists
  }
  else if(stat(dirname, &sb) == 0 && !S_ISDIR(sb.st_mode)){//file exists but is not a directory
      printf("subdir name is already taken by another file. rename file and run again\n");
      exit(1);
  }
  else{ //directory does not exist: make it.
    strcpy(makedirname,command);
    strcat(makedirname,subdirname);
    if(system(makedirname)==-1){ //make directory for data
      printf("warning: could not make directory %s. Exiting...\n",dirname);
      exit(1);
    }
  }
 
  
  /*** open files to write data ***/
  if(printinfo){
    snprintf(a_string,50,"%s/locationdifferencesinfo.dat",subdirname);
    celldistdata=fopen(a_string,"w");
    if(celldistdata==NULL){
      fprintf(stderr,"could not open file %s. exiting...\n",a_string);
      exit(1);
    }
    
    snprintf(a_string,50,"%s/cellpos.dat",subdirname);
    cellposdata=fopen(a_string,"w");
    if( cellposdata==NULL){
      fprintf(stderr,"could not open file %s. exiting...\n",a_string);
      exit(1);
    }
    snprintf(a_string,50,"%s/genes.dat",subdirname);
    genedata=fopen(a_string,"w");
    if( genedata==NULL){
      fprintf(stderr,"could not open file %s. exiting...\n",a_string);
      exit(1);
    }
  }

    //to test if such a text file does not get too large...
  snprintf(a_string,50,"%s/field.dat",subdirname);
  fielddata=fopen(a_string,"w");
  if( fielddata==NULL){
    fprintf(stderr,"could not open file %s. exiting...\n",a_string);
    exit(1);
  }
  
   
}

int inArray(int* array, int val, int length)
{
  for(int i=0;i<length;i++)
    if(array[i]==val)
      return 1;
  return 0;
}


void PrintField(FILE *file)
{
  int i,j;

  for(i=0;i<L;i++){
    for(j=0;j<W;j++){
      fprintf(file,"%d ",CellIdGrid[i][j]);
    }
    fprintf(file,"\n");
  }
  fprintf(file,"\n");
}

void printCellStates(FILE *file)
{  
  for(int i=0;i<anrcells_+1;i++){
    CellArray[i].outputCellState(file);
  }
  fprintf(file,"\n\n");
}


int main(int argc, char **argv)
{
  int step;
  int runnr; //for the simulation loop
  //double postdist;
  /*** initialization stuff ***/
  Start(argc,argv);
  
  dsfmt_init_gen_rand(&dsfmt,seed);
  
  /***********************/
  /*** SIMULATION LOOP ***/   
  /***********************/
  for(runnr=0;runnr<nrsimulations;runnr++){//simulation loop
    
    //opens subdir, and files to write to
    SubStart(runnr);
    
    /*** Development of Agent ***/
    
    A=new Agent();
    
    //  PlaceStartStripes(stripewidth,stripelength);
    A->CreateAgentManyCells(cellnumber, CELL, 0,1); 
    /* for(int i=1; i<=cellnumber;i++){
       if(i<=cellnumber/2)
       CellArray[i].oldstates.at(0)=0.2;
       else
       CellArray[i].oldstates.at(0)=0.8;
       }	*/
    fillStickyArray();

    //initiateBoundaries();	
    // checkBoundaryUpdate();
    //CellArray[1].CheckDeque(1,1);
    printf("start time loop\n");

    /****************************/
    
    /*****************/
    /*** time loop ***/
    /*****************/
    for(step=0;step<duration;step++){//timeloop
      
      // printf("%d\n",step);
      
      if(anrcells_>MaxNrCells){
	printf("too many cells. exiting...\n");
	break;
      }
      
      /*** take a picture and print data ***/
      /**************************************/
      
      if(!(step%pictureinterval)){
	//printf("step: %d\n",step);
	A->DetermineCollectioninfo();
	if(printinfo && step>=initialization){
	  PrintBlobInfo(celldistdata,step-initialization);
	  PrintCellPositions(cellposdata,(int)((step-initialization)/pictureinterval));
	  printCellStates(genedata);
	}
	
	A->ShowExpression(1,picturedirname, 0, 1);

	A->Snapshot(1, picturedirname); 
	//	checkBoundaryUpdate();
	//	printf("time: %d\n",step);
      }
      /**************************************/
      /**************************************/
      
      //  LocatePosterior(0);
      
      
      
      if(step>=initialization){//after initialization

	/*Hamiltonian*/
	A->UpdateAgent_Random(&Extra,&After);
	//	printf("long: %4.3lf, short: %4.3lf, vol:%d\n",CellArray[1].ReturnLongval(),CellArray[1].ReturnMinorAxisLength(), CellArray[1].vol);  
      }
      
    
      else 
	A->UpdateAgent_Random(NULL,&After);

        

    }//time loop

    //functions to close files and free memory before starting a new simulation
    ClearMem();
    
    EmptyVariables();
    
    if(printinfo){
      fclose(celldistdata);
      fclose(cellposdata);
      fclose(genedata);
    }
    fclose(fielddata);
  }//simulation loop

  printf("Hello World!\n");
  
}
