/* this program incorporates (and can combine) two mechanisms of convergent extension: blokjes and targetvector. If so desired, a growth zone is implemented from a third celltype, which divides and differentiates.*/

#include "Header.hh"
#include "Agent.hh"
#include "IO.hh"
#include "Misc.hh"
#include <sys/stat.h>

/**********************************************************/
/**************** program variables ***********************/
/**********************************************************/



struct stat sb;        //allows to check for the existence of directories
int order;             //neighbourhoodorder
dsfmt_t dsfmt;
dsfmt_t dsfmt_output;

//char despath[]="testrun";
/**********************************************************/
/******************* functions ****************************/
/**********************************************************/

void Start(int argc, char **argv)
{
  char command[]="mkdir ";
  char makedirname[50];

 ReadPars(argc, argv);
  //cout<<"seedinitpop: "<<seedinitpop<<endl;
  dsfmt_init_gen_rand(&dsfmt, seed);
  //dsfmt_init_gen_rand(&dsfmt, seedmutations);
  
  /*** make directory if necessary ***/
  if (stat(despath, &sb) == 0 && S_ISDIR(sb.st_mode)){//directory exists, just check if the cells directory exists.
    //do nothing
  }
    
  else if(stat(despath, &sb) == 0 && !S_ISDIR(sb.st_mode)){//file exists but is not a directory
    printf("given file is not a directory. give another name\n");
    exit(1);
  }
  else{ //directory does not exist: make it.
    strcpy(makedirname,command);
    strcat(makedirname,despath);
    if(system(makedirname)==-1){ //make directory for data
      printf("warning: could not make directory %s. Exiting...\n",despath);
      exit(1);
    }
  }
 
   
}


int main(int argc, char **argv)
{
  int i;

   /*** initialization stuff ***/
  Start(argc,argv);

  Agent *A=new Agent();
  
  /** do stuff with agent **/
    
  A->DevelopAgent(despath);  
  
  delete A;
  
  return 0;
}

