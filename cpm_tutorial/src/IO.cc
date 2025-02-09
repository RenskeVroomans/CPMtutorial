#include "IO.hh"
#include "Agent.hh"
#include "Cell.hh"

#include <boost/program_options.hpp>
#include <iostream> 
#include <string>
#include <vector>
#include <fstream>

char despath[500];
int seed;

using namespace boost;
namespace po = boost::program_options;
using namespace std;


int ReadPars(int argc, char* argv[])
{
  ///temporary pars for filling in Jtable
  int J11,J22, J33, J12,J13, J23, Jcell1med, Jcell2med, Jcell3med;
  
  /** ******command line reading******** **/
  
  // Declare the supported options.
  po::options_description generic("command line options");
  //this function actually lists the options to read from command line
  generic.add_options()
  ("help", "produce help message")
  ("parfile",po::value< vector <string> >(),"file with the parameters")
  ("despath,d",po::value< vector <string> >()->required(),"destination path")
  ("seed,s", po::value<int>(&seed)->default_value(503))
  ;
    
  //to make sure that you don't need to say --parfile="name.cfg", but just name.cfg
  po::positional_options_description p;
  p.add("parfile", -1);
  
  po::variables_map vm;//stores values of options, and can store values of arbitrary types. 
  //store, parse_command_line and notify functions cause vm to contain all the options found on the command line.
  po::store(po::command_line_parser(argc, argv).options(generic).positional(p).run(), vm);
  
  if (vm.count("help")) { //checks if option was specified
    cout << generic << "\n";
    //return 1;
    exit(1);
  }
  
  //convert CC string DESPATH to char *
  vector <string> vs=vm["despath"].as< vector<string> >();
  strcpy(despath,&vs[0][0]);
   
  
  cout <<"put files in folder: "<<&vs[0][0]<<endl;
  
  if(!(vm.count("parfile"))){
    cout << "no parfile specified. Using default values at own risk..." <<endl;
    //exit(1);
  }
 
  /** ******parfile reading******** **/
  
  else{
    
    /**specify the parameters here!!**/
    po::options_description config("parameter values");
    config.add_options()
    ///Agent.cc variables
    ("labdavol", po::value<double>(&Agent::labdavol)->default_value(1.), "stiffness of cell volume")
    ("labdasurf", po::value<double>(&Agent::labdasurf)->default_value(0.), "stiffness of membrane surface")
    ("temperature", po::value<double>(&Agent::invT)->default_value(12))
    ("targetvolume", po::value<int>(&Agent::targetvolume)->default_value(100))
    ("targetsurface", po::value<int>(&Agent::targetsurface)->default_value(36))
    ("surfaceconstraint", po::value<int>(&Agent::surfaceconstraint)->default_value(0))
    ("L", po::value<int>(&Agent::L)->default_value(100), "field length")
    ("W", po::value<int>(&Agent::W)->default_value(100), "field width")
    ("zoom", po::value<int>(&Agent::zoom)->default_value(1), "scaling of png images")
    ("neighbourhoodsize", po::value<int>(&Agent::neighbourhoodsize)->default_value(2), "size of the neighbourhood for Hamiltonian")
    ("NrDevSteps", po::value<int>(&Agent::NrDevSteps)->default_value(10000), "nr of steps for which agent is run")
    ("NrCellTypes", po::value<int>(&Agent::NrCellTypes)->default_value(2), "nr of cell types")
    ("CellPlace", po::value<int>(&Agent::CellPlace)->default_value(0), "Whether to place cells randomly(0) or clustered (1)")
    ("InitNrCells",po::value<int>(&Agent::InitNrCells)->default_value(10), "nr of cells to start with")
    ("picinterval",po::value<int>(&Agent::picinterval)->default_value(10), "interval for making pictures")

    //T cell migration
    ("mu",po::value<double>(&Agent::mu)->default_value(0.0), "strength of cell migration")
    ("perstime",po::value<int>(&Agent::perstime)->default_value(10), "duration of persistent walk")
    
    //     ///local vars, only for setting up

    ("Jcell1med", po::value<int>(&Jcell1med)->default_value(10))
    ("Jcell2med", po::value<int>(&Jcell2med)->default_value(10))
    ("Jcell3med", po::value<int>(&Jcell3med)->default_value(10))
    ("Jcell1cell1", po::value<int>(&J11)->default_value(16))
    ("Jcell1cell2", po::value<int>(&J12)->default_value(12))
    ("Jcell2cell2", po::value<int>(&J22)->default_value(16))
    ("Jcell3cell3", po::value<int>(&J33)->default_value(16))
    ("Jcell1cell3", po::value<int>(&J13)->default_value(16))
    ("Jcell2cell3", po::value<int>(&J23)->default_value(16))
    ;
    
  
    /******required type conversions: from vector<string> to ifstream for parse_config_file()*******/
    vector <string> vs=vm["parfile"].as< vector<string> >();
    cout<<"parfile used: "<<&vs[0][0]<<endl;
    ifstream parf(&vs[0][0]);
    /*************/
  
    po::store(po::parse_config_file(parf, config), vm);
    
    po::notify(vm);    
    /*************************************/
    Agent::invT=1./Agent::invT;
//     printf("neighsize1: %d\n",Agent::neighbourhoodsize);
//     printf("temp1: %lf\n",Agent::invT);
    
    //fill the Jtable with appropriate J values. 
    //arguments: nr of celltypes other than medium, number of Jtable entries (depends on nr types), and then entries as follows: celltype 1, celltype 2, J value.
    //remember to halve J values with the medium.  1,2 3,5, 4,9
    if(Agent::FillJTable(4,9,
      0,1,Jcell1med,
      0,2,Jcell2med,
      0,3,Jcell3med,
      1,1,J11,
      1,2,J12,
      2,2,J22,
      3,3,J33,
      1,3,J13,
      2,3,J23
    ))
    {
      exit(1);
    }
    
    //create neighbourhood array
    if(Agent::neighbourhoodsize>0)
      Agent::FindNeighbourhood();
    else
    {
      printf("invalid neighbourhoodsize. Exiting...\n");
      exit(1);
    }
    
    
  }
  
  return 0;
  
}