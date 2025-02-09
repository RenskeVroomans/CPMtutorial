#include "Agent.hh"
#include <GL/glut.h>
#include <GL/glx.h>
#include <png.h>
#include <zlib.h>
#include <string>

FPosition *array;



static int typecolormap[11][3]=
  {{255, 255, 255},//white, reserved for medium in celltypes plane
   {0, 0, 0}, //black: reserved for boundaries in celltypes plane
   {255, 0, 0},//red 2
   {0, 255, 0},//green 3 
   {0, 0, 255},//blue 4
   {255, 255, 0}, //yellow 5
   {0, 255, 255},//cyan 6
   {255, 0, 255},//magenta 7
   {255, 153, 51}, //orange 8(not the 'perfect' orange but 1 shade lighter for better distinction)
   {153, 51, 255}, //purple 9 (same as orange)
   {128, 128, 128} //gray 10
  };



void Agent::Snapshot(int zoom, char *dirname)
{
  int i,j;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int WW=zoom*W;
  const int LL=zoom*L;
  int celltypes[LL][WW];
  int cellids[LL][WW];//plane holding the enlarged field
  int boundaries[LL][WW];//plane holding the boundary pixels
  unsigned char RGBdata[WW*3];
  
  for(i=0;i<500;i++)
    fname[i]='\0';


  //blow up data
  for(i=0;i<L;i++)
    for(j=0;j<W;j++)
      {
	for(ii=0;ii<zoom;ii++)
	  for(jj=0;jj<zoom;jj++)
	    {
	      if(CellIdGrid[i][j]!=0)//there is a cell
		{
		  cellids[zoom*i+ii][zoom*j+jj]=CellIdGrid[i][j]+30;//+1;
		  // if(CellIdGrid[i][j]%2==0)//should be outcommented when coloring celltypes
		  celltypes[zoom*i+ii][zoom*j+jj]=CellArray[CellIdGrid[i][j]].celltype+1;//30
		  // else//should be outcommented when coloring celltypes
		  //celltypes[zoom*i+ii][zoom*j+jj]=50;
		}
	      else//there is no cell
		{
		  cellids[zoom*i+ii][zoom*j+jj]=0;//254
		  celltypes[zoom*i+ii][zoom*j+jj]=0;//254
		}
	    }
      }
 //find boundaries
  for(i=1;i<LL-1;i++)
    for(j=1;j<WW-1;j++)
      {
	if(cellids[i][j]!=cellids[i-1][j] ||
	   cellids[i][j]!=cellids[i+1][j] ||
	   cellids[i][j]!=cellids[i][j-1] ||
	   cellids[i][j]!=cellids[i][j+1])
	  boundaries[i][j]=1;
	else
	  boundaries[i][j]=0;
      }
  //add boundaries
  
  
  for(i=1;i<LL-1;i++)
    for(j=1;j<WW-1;j++)
      {
	if(boundaries[i][j]==1)
	  {
	    cellids[i][j]=0;//black in the colormap array
	    celltypes[i][j]=1;
	  }
      }
  
  for(i=2;i<6;i++)
    celltypes[i][3]=1;

  
  
  sprintf(fname,"%s/Ag_%.6d.png",despath,piccounter);

  PNGFileP=fopen(fname, "wb");
  // printf("I'm here!\n");
  if(PNGFileP==NULL){
    fprintf(stderr,"Agent::Snapshot: error in opening the pngs. Does the directory exist?\n");
  }

  else{

    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,
				      (png_error_ptr) NULL, 
				      (png_error_ptr) NULL );
    if(!png_ptr)
      {
	printf("out of memory\n");
	exit(1);
      }
    info_ptr = png_create_info_struct ( png_ptr );
    if(!info_ptr)
      {
	png_destroy_write_struct(&png_ptr, NULL);
	printf("out of memory\n");
	exit(1);
      }
    png_init_io ( png_ptr, PNGFileP );
    png_set_IHDR(png_ptr, info_ptr,WW,LL,
		 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
   
    // write header
    png_write_info ( png_ptr, info_ptr );
    
    // write out image, one row at a time

   
    for(i=0;i<LL;i++){
      //transform info of one row into color
      for(j=0;j<WW;j++)
	{
	  
	  RGBdata[j*3+0]=typecolormap[celltypes[i][j]][0];
	  RGBdata[j*3+1]=typecolormap[celltypes[i][j]][1];
	  RGBdata[j*3+2]=typecolormap[celltypes[i][j]][2];
	  
	}
      row_pointer =( RGBdata);
      //write row to png
      png_write_rows ( png_ptr, &row_pointer, 1 );

    }
    
    // flush all info to file
    png_write_end ( png_ptr, info_ptr );
    fflush ( PNGFileP );
    png_destroy_write_struct ( &png_ptr,&info_ptr);
    fclose(PNGFileP);
   
    //free (datablock);
    
    piccounter++;
  
  }
   
}

static int doubletoint(double val, double maxval)
{
  // if(val/maxval>0.5)
  //printf("doubletoint: %d\n", (int)(val/maxval*255.));
  if(val/maxval>1)
    return 255;
  else
    return (int)((val/maxval)*255.);
}

int Agent::ColourGradient(double **field,int zoom,double maxval, char *dirname) //field, zoom, directory to put pictures in. Returns success
{
  int i,j;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int WW=zoom*W;
  const int LL=zoom*L;
  int colour[LL][WW];//plane holding the enlarged field
  unsigned char RGBdata[LL*WW*3];
  static int piccounter=0;

  
  for(i=0;i<L;i++)
    for(j=0;j<W;j++)
      {
	for(ii=0;ii<zoom;ii++)
	  for(jj=0;jj<zoom;jj++)
	    {
	      colour[zoom*i+ii][zoom*j+jj]=doubletoint(field[i][j],maxval);
	    }
      }
  
  //transform into color
  for(i=0;i<LL;i++)
    for(j=0;j<WW;j++)
      {
	RGBdata[j*3+i*WW*3+0]=max(0,255-colour[i][j]);//R: 
	RGBdata[j*3+i*WW*3+1]=max(0,255-2*colour[i][j]);//G 
	RGBdata[j*3+i*WW*3+2]=max(0,255-4*colour[i][j]);//B
      }

  sprintf(fname,"%s/field%.4d.png",dirname,piccounter);
  PNGFileP=fopen(fname, "wb");
  if(PNGFileP==NULL){
    fprintf(stderr,"Agent::ColourProperties: error in opening the pngs. Does the directory exist?\n");
  }

  else{
    // unsigned char *datablock;
    int row;
    //datablock = (unsigned char *) malloc ( 3*LL*WW*sizeof(unsigned char));
    /*for(i=0;i<LL;i++)
      for(j=0;j<WW;j++)
	{
	  
	  datablock[3*i*WW+3*j+0]=RGBdata[j*3+i*WW*3+0];
	  datablock[3*i*WW+3*j+1]=RGBdata[j*3+i*WW*3+1];
	  datablock[3*i*WW+3*j+2]=RGBdata[j*3+i*WW*3+2];
	  }*/
   
    
     /*allocate memory for png struct, specify writing*/
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,
				      (png_error_ptr) NULL, 
				      (png_error_ptr) NULL );
    if(!png_ptr)
      {
	printf("out of memory\n");
	exit(1);
      } 
    
    /*allocate memory for png info struct*/
    info_ptr = png_create_info_struct ( png_ptr );
    if(!info_ptr)
      {
	png_destroy_write_struct(&png_ptr, NULL);
	printf("out of memory\n");
	exit(1);
      }
    png_init_io ( png_ptr, PNGFileP );
    png_set_IHDR(png_ptr, info_ptr,WW,LL,
		 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    /* write header */
    png_write_info ( png_ptr, info_ptr );
    /* write out image, one row at a time */
    for ( row = 0; row <LL; row++ )
      {
	row_pointer =( RGBdata + WW*row*3);
	png_write_rows ( png_ptr, &row_pointer, 1 );
      }
   
    /* flush all info to file */
    png_write_end ( png_ptr, info_ptr );
    fflush ( PNGFileP );
    png_destroy_write_struct ( &png_ptr,&info_ptr);
    fclose(PNGFileP);
    
    //free (datablock);
    
    piccounter++;
  }
  return 0;
}

void Agent::printField(int time)
{
  ofstream myfile;
  char myname[200];
  sprintf(myname, "%s/time%05d.text",despath,time);
  myfile.open(myname);
 
  for(int i=0; i<L; i++)
  {
    for(int j=0; j<W; j++)
    {
      myfile << CellIdGrid[i][j] <<"\t";
    }
    myfile << endl;
  }
  myfile << endl;
  map<int, Cell>::iterator it;  
  map<int, int>::iterator neigh;
  for(it=CellArray.begin(); it!=CellArray.end(); ++it)
  {
    myfile <<"cell "<< it->first <<": ";
    for(neigh=it->second.neighbours.begin(); neigh!=it->second.neighbours.end(); ++neigh)
    {
      myfile <<neigh->first<<" ";
    }
    myfile <<endl;
  }
  
  myfile.close();
      
}
