#include "Agent.hh"
#include "stdio.h"
#include <stdlib.h> 
#include <png.h>
#include <zlib.h>
#include <string>

void Agent::WriteGenome(int c)
{
  FILE *f;
  char fname[800];
  Genome::iter it;
  Gene *gene;
  TFBS *tfbs;
  int teller;
  double YY;
  double step=0.75;

  
  if(c==0)
    sprintf(fname,"%s/CodedGenomeAgent%.10d%s",despath,agentid,"original");
  else 
    sprintf(fname,"%s/CodedGenomeAgent%.10d_%d",despath,agentid,c);


  
  f=fopen(fname,"w");
  for(it=G->ChromBBList->begin();it!=G->ChromBBList->end();it++)
    {
      if(G->IsGene(*it))
	{
	  gene=dynamic_cast<Gene *>(*it);
	  fprintf(f,"%i ",0);//to indicate that it is a gene
	  fprintf(f,"%i ",(gene->type));//type of gene
	}
      else if(G->IsTFBS(*it))
	{
	  tfbs=dynamic_cast<TFBS *>(*it);
	  if(tfbs->weight<0)
	    fprintf(f,"%i ",-1);//to indicate that it is a repressing tfbs
	  else
	    fprintf(f,"%i ",1);//to indicate that it is a activating tfbs
	  fprintf(f,"%i ",(tfbs->type));//type of tfbs
	}
    }
      
  fclose(f);

  if(c==0)
    sprintf(fname,"%s/GenomeAgent%.10d%s.dot",despath,agentid,"original");
  else 
    sprintf(fname,"%s/GenomeAgent%.10d_%d.dot",despath,agentid,c);
   
  f=fopen(fname,"w");
  fprintf(f,"graph genome_%i {\n",agentid);
  fprintf(f,"size=\"10,4\";\n");
  teller=0;
  for(it=G->ChromBBList->begin();it!=G->ChromBBList->end();it++)
    {
      if(G->IsGene(*it))
	{
	  gene=dynamic_cast<Gene *>(*it);
	  YY=step*teller;
	  fprintf(f,"\"%i\" [label=\"%i\",shape=rectangle,color=blue,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*gene).type,YY,0.);
	}
      else if(G->IsTFBS(*it))
	{
	  tfbs=dynamic_cast<TFBS *>(*it);
	  YY=step*teller;
	  if((*tfbs).type<NrMatGeneTypes)
	    {
	      if((*tfbs).weight==1)
		fprintf(f,"\"%i\" [label=\"%i\",shape=octagon,color=green,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	      else
		fprintf(f,"\"%i\" [label=\"%i\",shape=octagon,color=red,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	    }
	  else if((*tfbs).type<(NrMatGeneTypes+NrSignGeneTypes))
	    {
	      if((*tfbs).weight==1)
		fprintf(f,"\"%i\" [label=\"%i\",shape=diamond,color=green,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	      else
		fprintf(f,"\"%i\" [label=\"%i\",shape=diamond,color=red,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	    }
	  else
	    {
	      if((*tfbs).weight==1)
		fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=green,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	      else
		fprintf(f,"\"%i\" [label=\"%i\",shape=ellipse,color=red,fontsize=16., style=filled,pos=\"%f,\%f!\"];\n",teller,(*tfbs).type,YY,0.);
	    }
	    
	}
      teller++;
    }
  fprintf(f,"}\n");
  fclose(f);
}

void Agent::WriteBasicProperties(int aid)
{
  FILE *f;
  char fname[800];

  sprintf(fname,"%s/%s",despath,"BasicPropertiesOfAgents");
  f=fopen(fname,"a");
  fprintf(f,"%i\t",aid);//1
  fprintf(f,"%i\t",axislength);//2
  fprintf(f,"%i\t",CellArray.size());//3
  fprintf(f,"%i\t",deadcells);//4
  fprintf(f,"%f\t",fitness);//5
  fprintf(f,"%f\t",glpenalty);//6
  fprintf(f,"%i\t",loosepenalty);//7
  fprintf(f,"\n");
  fclose(f);
}


void Agent::DetermineGenomeAndNetworkProperties(char *filename)
{
  FILE *f;
  char fname[800];
  int i;

  Genome:: iter it;
  TFBS *tfbs;
  nrplustfbs=0;
  nrmintfbs=0;
  int indegree[200];
  int outdegree[200];
  int totdegree[200];
  int nrtfbsmatchinggenetype[NrGeneTypes];
  double indegreedistr[15];
  double outdegreedistr[15];
  double totdegreedistr[15];
  /*static double avgindegreedistr[15];
  static double avgoutdegreedistr[15];
  static double avgtotdegreedistr[15];*/ //seem to serve no purpose
  int deggenetypes[16][3]; //my addition

  for(i=0;i<200;i++)
    {
      indegree[i]=0;
      outdegree[i]=0;
      totdegree[i]=0;
     
    }
  for(i=0;i<NrGeneTypes;i++)
  {
    nrtfbsmatchinggenetype[i]=0; //outdegree
    deggenetypes[i][0]=0;//indegree (total not corrected for duplicates)
    deggenetypes[i][1]=0;//outdegree (total not corrected for duplicates)
    deggenetypes[i][2]=0;//nr of copies of this gene type
    
  }
  for(i=0;i<15;i++)
    {
      indegreedistr[i]=0;
      outdegreedistr[i]=0;
      totdegreedistr[i]=0;
    }

  int genecounter=0;  
  it=G->ChromBBList->begin();
 
  while(it!=G->ChromBBList->end())
    {
      if(G->IsTFBS(*it))//tfbs
	{
	  tfbs=dynamic_cast<TFBS *>(*it);
	  if(tfbs->weight==1)
	    nrplustfbs++;
	  else
	    nrmintfbs++;

	  indegree[genecounter]++; //don't take into account that there may be multiple genes of the same type impinging on this TFBS
	  totdegree[genecounter]++;
	  nrtfbsmatchinggenetype[tfbs->type]++;
	}
      else//gene
	genecounter++;
      
      it++;
    }
 
 Gene *gene;
  genecounter=0;
  it=G->ChromBBList->begin();
  while(it!=G->ChromBBList->end())
    {
      if(G->IsTFBS(*it))//tfbs
	;
      else
	{
	  gene=dynamic_cast<Gene *>(*it);
	  deggenetypes[gene->type][0]+=indegree[genecounter];
	  outdegree[genecounter]=nrtfbsmatchinggenetype[gene->type];
	  deggenetypes[gene->type][1]+=outdegree[genecounter];
	  totdegree[genecounter]+=nrtfbsmatchinggenetype[gene->type];
	  deggenetypes[gene->type][2]++;
	  genecounter++;
	}
      it++;
    }
    
  //print in and out degree per gene (when duplicated, it's an average) Also note that if a gene gets input from duplicated gene, this is only counted as 1! (one TFBS)
  sprintf(fname,"%s/%s%s",despath,"DegreesPerGeneType_",filename);
  f=fopen(fname,"w");
  for(i=0; i<NrGeneTypes; i++)
    fprintf(f, "%d\t%.2lf\t%.2lf\t%d\n",i,(double)deggenetypes[i][0]/(double)deggenetypes[i][2], (double) deggenetypes[i][1]/(double)deggenetypes[i][2],deggenetypes[i][2]);
 
  avgindegree=0;
  maxindegree=0;
  minindegree=1000;
  stdindegree=0;
  avgoutdegree=0;
  maxoutdegree=0;
  minoutdegree=1000;
  stdoutdegree=0;
  avgtotdegree=0;
  maxtotdegree=0;
  mintotdegree=1000;
  stdtotdegree=0;
  for(i=0;i<genecounter;i++)
    {
      if(indegree[i]<15)
	indegreedistr[indegree[i]]++;
      else
	indegreedistr[14]++;
      if(outdegree[i]<15)
	outdegreedistr[outdegree[i]]++;
      else
	outdegreedistr[14]++;
      if(totdegree[i]<15)
	totdegreedistr[totdegree[i]]++;
      else
	totdegreedistr[14]++;

      avgindegree+=indegree[i];
      avgoutdegree+=outdegree[i];
      avgtotdegree+=totdegree[i];
      
      stdindegree+=indegree[i]*indegree[i];
      stdoutdegree+=outdegree[i]*outdegree[i];
      stdtotdegree+=totdegree[i]*totdegree[i];

      if(indegree[i]>maxindegree)
	maxindegree=indegree[i];
      if(indegree[i]<minindegree)
	minindegree=indegree[i];
      
      if(outdegree[i]>maxoutdegree)
	maxoutdegree=outdegree[i];
      if(outdegree[i]<minoutdegree)
	minoutdegree=outdegree[i];

      if(totdegree[i]>maxtotdegree)
	maxtotdegree=totdegree[i];
      if(totdegree[i]<mintotdegree)
	mintotdegree=totdegree[i];
    }
  for(i=0;i<15;i++)
    {
      indegreedistr[i]/=(double)genecounter;
      outdegreedistr[i]/=(double)genecounter;
      totdegreedistr[i]/=(double)genecounter;
    }
 /* for(i=0;i<15;i++)
    {
      avgindegreedistr[i]+=indegreedistr[i];
      avgoutdegreedistr[i]+=outdegreedistr[i];
      avgtotdegreedistr[i]+=totdegreedistr[i];
    }
    for(i=0;i<15;i++)
    {
      avgindegreedistr[i]/=100;
      avgoutdegreedistr[i]/=100;
      avgtotdegreedistr[i]/=100;
    }*/



  avgindegree/=(double)genecounter;
  avgoutdegree/=(double)genecounter;
  avgtotdegree/=(double)genecounter;

  stdindegree/=(double)genecounter;
  stdindegree=stdindegree-avgindegree*avgindegree;
  stdindegree=sqrt(stdindegree);

  stdoutdegree/=(double)genecounter;
  stdoutdegree=stdoutdegree-avgoutdegree*avgoutdegree;
  stdoutdegree=sqrt(stdoutdegree);

  stdtotdegree/=(double)genecounter;
  stdtotdegree=stdtotdegree-avgtotdegree*avgtotdegree;
  stdtotdegree=sqrt(stdtotdegree);


  sprintf(fname,"%s/%s%s",despath,"DegreeDistr_",filename);
  f=fopen(fname,"w");
  for(i=0;i<15;i++)
  {
    fprintf(f,"%i\t%f\t%f\t%f\n",i,indegreedistr[i],outdegreedistr[i],totdegreedistr[i]);
    //fprintf(f,"%f\t%f\t%f\n",avgindegreedistr[i],avgoutdegreedistr[i],avgtotdegreedistr[i]);
  }
  fclose(f);

  /*for(i=0;i<15;i++)
  {
    avgindegreedistr[i]=0;
    avgoutdegreedistr[i]=0;
    avgtotdegreedistr[i]=0;
  }*/
  
}


void Agent::WriteGenomeAndNetworkProperties(int aid)
{
  FILE *f;
  char fname[800];

  sprintf(fname,"%s/%s",despath,"GenomeAndNetworkProperties");
  f=fopen(fname,"a");
  fprintf(f,"%i\t",aid);//1
  fprintf(f,"%i\t",G->glength_);//2
  fprintf(f,"%i\t",G->gnrgenes_);//3
  //fprintf(f,"%i\t",nridgenes);//4
  fprintf(f,"%i\t",G->gnrtfbs_);//5
  fprintf(f,"%i\t",nrplustfbs);//6
  fprintf(f,"%i\t",nrmintfbs);//7
  fprintf(f,"%f\t",avgindegree);//8
  fprintf(f,"%f\t",avgoutdegree);//9
  fprintf(f,"%f\t",avgtotdegree);//10
  fprintf(f,"%f\t",minindegree);//11
  fprintf(f,"%f\t",minoutdegree);//12
  fprintf(f,"%f\t",mintotdegree);//13
  fprintf(f,"%f\t",maxindegree);//14
  fprintf(f,"%f\t",maxoutdegree);//15
  fprintf(f,"%f\t",maxtotdegree);//16
  fprintf(f,"%f\t",stdindegree);//17
  fprintf(f,"%f\t",stdoutdegree);//18
  fprintf(f,"%f\t",stdtotdegree);//19
  fprintf(f,"\n");
  fclose(f);
}


bool compare_ints (int first,int second)
{
  return ( first < second);
}

void Agent::removeDoubles_type()
{
  nrpos=0;
  nrneg=0;
  
  vector<int> :: iterator il, el;
  vector<int>::iterator tal, ntal;

  list< vector<int> >::iterator walk, next, cwalk, cnext, pwalk, pnext; //loopstore, copy of loopstore, counterstore
  
  int ns=0;
  list< vector<int> > copy;
  
  copy=loopstore;
   
  for(walk=loopstore.begin(), tal=tallying.begin(), cwalk=copy.begin(), pwalk=counterstore.begin(); walk!=loopstore.end(); ++walk, ++tal, ++cwalk, ++pwalk)
  {
    next=walk;
    next++;
    cnext=cwalk;
    cnext++;
    pnext=pwalk;
    pnext++;
    ntal=tal;
    ntal++;
    
    while(next!=loopstore.end()){
      ns=0;
      
      if ((int)(*walk).size()==(int)(*next).size() && (*tal)==(*ntal) )
      {
	sort((*cwalk).begin(),(*cwalk).end());
	sort((*cnext).begin(),(*cnext).end());
	for(il=(*cwalk).begin(),el=(*cnext).begin(); il!=(*cwalk).end(); ++il, ++el)
	{
	  //printf("sorted: %d\n",(*il));
	  if((*il)!=(*el))
	  {
	    ns=1;
	    break;
	  }
	}
	if(!ns){ //the two vectors are the same
	  next=loopstore.erase(next);
	  cnext=copy.erase(cnext);
	  pnext=counterstore.erase(pnext);
	  ntal=tallying.erase(ntal);
	  continue; //don't have to move the iterator anymore
	}
      }
      ++next;
      ++cnext;
      ++ntal;
      ++pnext;
    }
    
    if((*tal)>0)
      nrpos++;
    else
      nrneg++;
  }//loop through loopstore
  
}

void Agent::removeDoubles_id(int finalcall)
{
//  removing double counted loops from loopstore, tallying the number of positive and negative loops 
  
  vector<int>::iterator tal, ntal, il, el;
  list< vector<int> >::iterator walk, next, cwalk, cnext;
  
  int ns=0;
  //list< vector<int> > copy;
  
  //copy=loopstore;
 
  for(walk=loopstore.begin(), tal=tallying.begin(), cwalk=counterstore.begin(); walk!=loopstore.end(); ++walk, ++tal, ++cwalk)
  {
    next=walk;
    next++;
    cnext=cwalk;
    cnext++;
    ntal=tal;
    ntal++;
    
    while(next!=loopstore.end()){
      ns=0;
      
      if ((int)(*walk).size()==(int)(*next).size() && (*tal)==(*ntal)) //same length and same sign
      {
	//printf("comparing...\n");
	sort((*cwalk).begin(),(*cwalk).end());
	sort((*cnext).begin(),(*cnext).end());
	
	for(il=(*cwalk).begin(),el=(*cnext).begin(); il!=(*cwalk).end(); ++il, ++el)
	{
	  //printf("sorted: %d\n",(*il));
	  if((*il)!=(*el))
	  {
	    ns=1;
	    break;
	  }
	}
	if(!ns){
	  next=loopstore.erase(next);
	  cnext=counterstore.erase(cnext);
	  ntal=tallying.erase(ntal);
	  continue;
	}
      }
      ++next;
      ++cnext;
      ++ntal;
    }
    if (finalcall){
      if((*tal)>0)
	nrpos++;
      else
      nrneg++;
    }
  }
}

void Agent::DetermineLoopAndMotifProperties()
{
  nrposauto=0;
  nrnegauto=0;
  nrpos=0;
  nrneg=0;
  int A[200][200];

  compreg=0; //number of complex (+-) regulations
  //printf("am here!\n");
  //int A[200][200];
  int genetypes[200];
  int i,j;
  for(i=0;i<200;i++)
    for(j=0;j<200;j++)
      A[i][j]=0;
  for(i=0;i<200;i++)
    genetypes[i]=-1;
  
  
  int genecounter=0;  
  Genome::iter it;
  TFBS *tfbs;
  Gene *gene;

  //which genes are located where in genome
  it=G->ChromBBList->begin();
  while(it!=G->ChromBBList->end())
    {
      if(G->IsGene(*it))
	{
	  gene=dynamic_cast<Gene *>(*it);
	  genetypes[genecounter]=gene->type;
	  genecounter++;
	}
      it++;
    }

  //which type of tfbs have genes upstream of them
  //multiple genes of the same type can bind to this tfbs
  int genecounter2=0;
  it=G->ChromBBList->begin();
  while(it!=G->ChromBBList->end())
  {
    if(G->IsTFBS(*it))//tfbs
    {
      tfbs=dynamic_cast<TFBS *>(*it);
      for(i=0;i<genecounter;i++)
      {
    
	if(tfbs->type==genetypes[i])
	{
	  if(genetypes[i]==genetypes[genecounter2]) //autoloop: don't put in matrix
            if(tfbs->weight<0)
	      nrnegauto++;
	    else
	      nrposauto++;

	  else if((tfbs->weight<0)==(A[i][genecounter2]<0) || !A[i][genecounter2]) //same sign or A is still empty
		A[i][genecounter2]+=tfbs->weight;
	  
	  else
	  {
	    A[i][genecounter2]=1; //a choice has to be made here: A is integer, so cannot account for having both neg and pos interaction. we'll count it as an "on" interaction for now.
	    compreg++;
	  }
	}
      }
    }
    else
      genecounter2++;
    it++;
  }
  
 
  
  //debugging
 /* for(i=0;i<200;i++)
  { 
    for(j=0;j<200;j++)
      if(A[i][j])
	printf("gene %d %d gene %d\n",genetypes[i], A[i][j], genetypes[j]);
  }*/
  
//    go through interaction matrix, find the loops 
  
  vector <int> nodes, counter; //temporary storage for the genetypes in the path followed, their genecounter (position in genome, unique identifier), and to keep score whether neg or pos loop.
  int genenr; 
  int tally=1; //to find whether loop is positive or negative
  
  int index=0; //to find position of iterator
  int found=0;
  
  int loopsize=0; //keep track of loop size: break free if it becomes too big... NEW!!
  
  vector <int> :: iterator il, el;
   
  for(int row=0; row<genecounter; row++)
  {
    //empty the temporary storage
    nodes.clear();
    counter.clear();
    //push first gene in list
    nodes.push_back(genetypes[row]);
    counter.push_back(row);
    tally=1;
    genenr=row;
    loopsize=1;
    //printf("start gene %d, ",nodes.back());
    //printf("row = %d\n", row);
    
    for(int col=0; col<genecounter; col++) //go through all genes that gene "genenr" regulates
    {
      if(A[genenr][col] && find(counter.begin(),counter.end(), col)==counter.end()){ //gene "genenr" regulates this gene, we have no loop yet
	//push gene at end of list
	nodes.push_back(genetypes[col]);
	counter.push_back(col);
	//tally*=A[genenr][col];
	genenr=col;
	col=-1;
	loopsize++;
	//printf("tally: %d\n",tally);
	//printf("next gene %d, counter=%d ", nodes.back(), counter.back());
      }
      else if (A[genenr][col] && (el=find(counter.begin(), counter.end(), col))!=counter.end()) //found a loooop
      {
	//to find the same position in the nodes vector
	index=el-counter.begin();
	il=nodes.begin()+index;
	//tally*=A[genenr][col];
	vector <int> vec (il, nodes.end());
	vector <int> cvec (el, counter.end());
	//store whether pos or neg fb loop, and the loop itself both as gene type and the unique identifier.
	//tallying.push_back(tally);
	counterstore.push_back(cvec);
	loopstore.push_back(vec);
	genenr=counter.back();
	found=1;
	
	for (int i=0; i<vec.size(); i++)
	{
	  //printf("gene %d, id %d\t",vec[i],cvec[i]);
	  if(i<vec.size()-1) tally*=A[cvec[i]][cvec[i+1]];
	  else tally*=A[cvec[i]][cvec[0]]; 
	  
	}

	tallying.push_back(tally);
	tally=1;
      }
      
            
      while( (found || loopsize >15 || (!A[row][col]  && col==genecounter-1) ) && nodes.size()>1 ) //iterate back to see if we can find more loops
      {
	//if(A[row][col]){
	  //printf("changing tally; was %d, is now %d; row=%d, col=%d\n",tally,tally/A[row][col], row, col );
	  //tally=tally/A[row][col];
	//}
	nodes.pop_back();
	//genenr=nodes.back();
	col=counter.back();
	counter.pop_back();
	genenr=counter.back();
	loopsize--;
	found=0; 
	
	//printf("back to gene %d, counter %d\n", nodes.back(), col);
      }
      
    }//col loop
    
    removeDoubles_id(0);
  } //row loop 
 
  removeDoubles_id(1);
}

void Agent::WriteLoopAndMotifProperties(int aid, char *name)
{
  FILE *f, *f2;
  char fname[800];
  int countarray[20];
    
  for(int i=0;i<20;i++)
    countarray[i]=0;
    
  sprintf(fname,"%s/%s_%s",despath,"LoopAndMotifProperties", name);
  f=fopen(fname,"w");
  
 list<vector<int> >::iterator walk; 
 vector<int>::iterator el; 
  

  /** some generic properties **/
  
  fprintf(f,"%i\t",aid);//1
  fprintf(f,"%i\t",nrposauto);//2
  fprintf(f,"%i\t",nrnegauto);//3
  fprintf(f,"%d\t", nrpos); //4
  fprintf(f,"%d\t", nrneg);//5
  removeDoubles_type();
  fprintf(f,"%d\t", nrpos); //6 //the non-redundant loops
  fprintf(f,"%d\n", nrneg);//7
  fprintf(f, "\n\n");
  //fprintf(f,"\n\nLoops of size:");
  /** print the loops themselves to a separate file **/
  for(walk=loopstore.begin(); walk!=loopstore.end(); ++walk)
  {
    if((*walk).size()<=15)
      countarray[(*walk).size()]++; //tally loop sizes
      else
	countarray[16]++;
      
    for(el=(*walk).begin(); el!=(*walk).end();++el) //print loop
   {
     fprintf(f,"%d\t",(*el));
   }
   fprintf(f,"\n\n");
  }
  
  fclose(f);
  
  /** ****** */
  sprintf(fname,"%s/%s_%s",despath,"LoopSizeFreqs", name); //corrected for doubles in gene type! 
  f2=fopen(fname,"w");
  
  for(int i=0; i<16;i++)
    fprintf(f2,"%d\t%d\n", i, countarray[i]); 
  
  fclose(f2);
  
}