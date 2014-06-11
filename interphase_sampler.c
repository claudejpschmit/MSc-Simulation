/* MC code for a random/self-avoiding walk with sticky sites*/
/* Compile line on LINUX
cc interphase_sampler.c -o interphase_sampler -O3 -lm -LNO
OR
gcc interphase_sampler.c -o interphase_sampler -O3 -lm -LNO
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*useful constants*/
#define Pi 3.141592653589793
#define TwoPi 6.283185307179586
#define sqr2 1.4142136

/*program parameters*/
#define L 151 // size of DNA
#define loop 10 // distance between successive polymerases
#define nstep 10000000 // number of steps

/*function and subroutine declarations*/
void inve(double costh, double sinth, double cospsi, double sinpsi, double cosom, double sinom, double rot[3][3]);
void dir(double costh, double sinth, double cospsi, double sinpsi, double cosom, double sinom, double rot[3][3]);
void streamfile(int step);
void pivot(int idum, double zn[L][3]);
void krank(int idum, double zn[L][3]);
void metropolis(double ti, int idum, double zn[L][3]);

/*global variables*/

int polymerase[L];
int interaction[L][L];

double z[L][3]; // polymer

double e0=1.0;
double kb=0.0; // stiffness, use for semiflexible polymers
double en,ev;

FILE *fp,*output4;
char name1[40];

  /*      open(unit=10,file='parameters.dat',status='unknown')
      write(10,*) 'energy between polymerases = ', e0
      write(10,*) 'loop length (nm / kbps) = ', 
     1     dfloat(loop)*8.5d0,' / ', dfloat(loop)*8.5d0*0.8d0/11.d0
     close(10)*/

int main(int argc, char** argv)
{
  int i,j,k,n,seed,nmoved;
  int poldist,index;
  int jmax,ind1,ind2,nr;
  int ct,ct1; // determine how often output gets written
  double z_new[L][3]; // polymer trial position
  double r,deltaE,enmax;
  double a1,a2;
  double t,ti;
  FILE *input1,*fp1,*fp2;

  seed=atof(argv[0]);
  poldist=loop;

  // L should be eg loop*15+1
  
  enmax=0.0;

  for(i=0;i<L;i++){
    polymerase[i]=0;
    if(i%poldist==0) polymerase[i]=1;
    for(j=0;j<L;j++){
      interaction[i][j]=0;
    }
  }
    
  if((input1 = fopen("interaction.dat","r"))==NULL)
    {
      printf("cannot open input file interaction.dat\n");
      exit(0);
    }

  
  fscanf(input1,"%d %lf",&jmax,&t);
      
  for(i=0;i<L;i++){
    if(polymerase[i]==1){
      for(j=0;j<jmax;j++){
	ind1=i;
	ind2=poldist*j;
	fscanf(input1,"%d",&interaction[ind1][ind2]); 
      }
    }
  }
    
  fclose(input1);
  
  fp1 = fopen("interactionparameters.dat","w");
  
  for(i=0;i<L;i++){
    for(j=0;j<L;j++){
      if(polymerase[i]+polymerase[j]==2) {
	fprintf(fp1,"%d %d %d\n",i,j,interaction[i][j]);
      }
    }
  }
    
  fclose(fp1);

  fp2= fopen("polymerasepositions.dat","w");
  
  for(i=0;i<L;i++){
    if(polymerase[i]!=0) fprintf(fp2,"%d\n",i);
  }
  
  fclose(fp2);

  if((output4 = fopen("energy.dat","w"))==NULL)
    {
      printf("cannot open energy output file\n");
      
      exit(0);
    }

  for(index=0;index<1;index++){

    for(i=0;i<L;i++){
      a1=drand48();
      a2=drand48();
      z[i][0]=0.0;
      z[i][1]=0.0;
      z[i][2]=i;
    }
    en=0;

    nr=0;
    kb=0.0;

    for(n=0;n<nstep;n++){

      if(n%1000==0) t=t*0.999;
  
      ti=1.0/t;

      ct=ct1=1000; 

      ev=en;
      a1=drand48();

      if(a1<0.5){
	pivot(seed,z_new); // move end
      }
      else{
	krank(seed,z_new); //move internal region
      }

      //     check self and mutual avoidance          
      metropolis(ti,seed,z_new);

      if ((n%ct1==0)||(n==0)) {
	streamfile(n);
	fprintf(output4,"%lf %lf\n",t,en);
	fflush(output4);

	}

    }
  }

  fflush(output4);
  	
}

/*----------------------- functions and subroutines ----------------------*/

double dist(double r1[3],double r2[3])
{

  double distance2,distance;

  distance2=pow(r1[0]-r2[0],2.0)+pow(r1[1]-r2[1],2.0)+pow(r1[2]-r2[2],2.0);
  distance=sqrt(distance2);

  return distance;

}

void pivot(int idum, double zn[L][3])
{

  int i,k,k1;
  int nm,nm1;
  double beta,teta,fi;
  double cospsi,sinpsi,costh,sinth,cosom,sinom;
  double a1,a2;
  double ver,sum;
  double rot[3][3];
  double pos[L][3],pos1[L][3];
        
  a1=drand48();
  a2=drand48();
  ver=1;
  if (a2<0.5) ver=-1;
  beta=a1*Pi/1.0*ver;        //max angle pi
  
  a1=drand48();
  ver=1;
  teta=acos(a1);        //max angle pi

  a1=drand48();
  a2=drand48();
  ver=1;
  if (a2<0.5) ver=-1;
  fi=a1*Pi/1.0*ver;        //max angle pi
  
  a1=drand48();
  nm=floor(a1*(L-1));
  nm1=L-1-nm;
  
  for(i=0;i<nm1+1;i++){
    for(k=0;k<3;k++){
      pos[i][k]=z[nm+i][k]-z[nm][k];
    }
  }
  
  costh=cos(teta);
  sinth=sin(teta);
  cospsi=cos(fi);
  sinpsi=sin(fi);
  cosom=cos(beta);
  sinom=sin(beta);
  
  dir(costh,sinth,cospsi,sinpsi,cosom,sinom,rot);
  
  for(i=0;i<nm1+1;i++){
    for(k=0;k<3;k++){
      sum=0;
      for(k1=0;k1<3;k1++){
	sum += rot[k][k1]*pos[i][k1];
      }
      pos1[i][k]=sum;
    }
  }
  
  inve(costh,sinth,cospsi,sinpsi,cosom,sinom,rot);
  
  for(i=0;i<nm1+1;i++){
    for(k=0;k<3;k++){
      sum=0;
      for(k1=0;k1<3;k1++){
	sum += rot[k][k1]*pos1[i][k1];
      }
      pos[i][k]=sum;
    }
  }
  
  for(i=0;i<nm1+1;i++){
    for(k=0;k<3;k++){
      zn[nm+i][k]=pos[i][k]+z[nm][k];
    }
  }
  
  if(nm!=0){
    for(i=0;i<nm;i++){
      for(k=0;k<3;k++){
	zn[i][k]=z[i][k];
      }
    }    
  }

}

void krank(int idum, double zn[L][3])
{

  int nm,nm1,nf;
  int i,k,k1;
  double a1,a2;
  double ver;
  double beta;
  double rag,sum;
  double costh,sinth,cospsi,sinpsi,cosom,sinom;
  double c[3],base2[3];
  double pos[L][3],pos1[L][3];
  double rot[3][3];
 
  Choose: a1=drand48(); 
  nm=floor(a1*L);
  if((nm>=L-2)||(nm==0)) goto Choose;

  a1=drand48();
  nm1=floor(a1*100+2);
  if((nm+nm1)>=L-1) nm1=L-1-nm;

  a1=drand48();
  a2=drand48();
  ver=1;
  if(a2<0.5) ver=-1;
  beta=a1*Pi/1.0*ver;           // max angle pi
  for(k=0;k<3;k++){
    base2[k]=z[nm+nm1][k];
  }
  nf=nm+nm1;

  for(k=0;k<3;k++){
    c[k]=z[nm][k];
  }
  for(i=0;i<nm1+1;i++){
    for(k=0;k<3;k++){
      pos[i][k]=z[nm+i][k]-c[k];
    }
  }

  rag=pow(pos[nm1][0],2.0)+pow(pos[nm1][1],2.0)+pow(pos[nm1][2],2.0);
  rag=sqrt(rag);
  costh=pos[nm1][2]/rag;
  sinth=sqrt(1.0-costh*costh);
  if (sinth<0.0001) goto ReAssembleChain;
  cospsi=pos[nm1][0]/(rag*sinth);
  sinpsi=pos[nm1][1]/(rag*sinth);
  cosom=cos(beta);
  sinom=sin(beta);
  dir(costh,sinth,cospsi,sinpsi,cosom,sinom,rot);  

  for(i=0;i<nm1+1;i++){
    for(k=0;k<3;k++){
      sum=0.0;
      for(k1=0;k1<3;k1++){
	sum += rot[k][k1]*pos[i][k1];
	  }
      pos1[i][k]=sum;
    }
  }

  inve(costh,sinth,cospsi,sinpsi,cosom,sinom,rot);
  
  for(i=0;i<nm1+1;i++){
    for(k=0;k<3;k++){
      sum=0.0;
      for(k1=0;k1<3;k1++){
	sum += rot[k][k1]*pos1[i][k1];
	  }
      pos[i][k]=sum;
    }
  }

 ReAssembleChain: for(i=0;i<nm1+1;i++){
    for(k=0;k<3;k++){
      zn[nm+i][k]=pos[i][k]+c[k];
    }  
  }



  if(nm!=0) {
    for(i=0;i<nm;i++){
      for(k=0;k<3;k++){
	zn[i][k]=z[i][k];
      }
    }
  }

  if((nm+nm1) != L-1){
    for(i=nm+nm1+1;i<=L-1;i++){
      for(k=0;k<3;k++){
	zn[i][k]=z[i][k];
      }
    }
  }

  for(k=0;k<3;k++){
    if (abs(base2[k]-zn[nm+nm1][k])> 0.01){
	fprintf(stderr,"aarrghhh!!!!\n");
	fprintf(stderr,"%d %lf\n",k,base2[k]-zn[nm+nm1][k]);
      }
  }   
      
}  

void inve(double costh, double sinth, double cospsi, double sinpsi, double cosom, double sinom, double rot[3][3])
{

  rot[0][0]=cospsi*costh;
  rot[0][1]=-sinpsi;
  rot[0][2]=cospsi*sinth;
  rot[1][0]=sinpsi*costh;
  rot[1][1]=cospsi;
  rot[1][2]=sinpsi*sinth;
  rot[2][0]=-sinth;
  rot[2][1]=0;
  rot[2][2]=costh;  

}


void dir(double costh, double sinth, double cospsi, double sinpsi, double cosom, double sinom, double rot[3][3])
{

  rot[0][0]=cosom*costh*cospsi+sinom*sinpsi;
  rot[0][1]=cosom*costh*sinpsi-sinom*cospsi;
  rot[0][2]=-sinth*cosom;
  rot[1][0]=sinom*cospsi*costh-sinpsi*cosom;
  rot[1][1]=sinpsi*sinom*costh+cosom*cospsi;
  rot[1][2]=-sinth*sinom;
  rot[2][0]=cospsi*sinth;
  rot[2][1]=sinpsi*sinth;
  rot[2][2]=costh;       
      
}

void metropolis(double ti, int idum, double zn[L][3])

{

  int i,i1,i2;
  double r0,r1,a1,asc;
  double rag1;
  double v1[3],v2[3];
  double tang[L][3];

  r0=0.95;
  r1=1.4;
  
  for(i1=0;i1<L-1;i1++){
    for(i2=0;i2<3;i2++){
      tang[i1][i2]=zn[i1+1][i2]-zn[i1][i2];
    }
  }
	 
  en=0.0;
  
  for(i1=0;i1<L-2; i1++){    
    en -= kb*(tang[i1][0]*tang[i1+1][0]+tang[i1][1]*tang[i1+1][1]+tang[i1][2]*tang[i1+1][2]);
  }

         
  for(i=0; i<L-2; i++){
    for(i2=i+1; i2<L-1; i2++) {

      if(i2!=i+1){
	
	if(polymerase[i]+polymerase[i2]==2){
	  for(i1=0;i1<3;i1++){
	    v1[i1]=zn[i2][i1];
	    v2[i1]=zn[i][i1];
	  }
	  rag1=dist(v1,v2);
	  if(rag1<r0){
	    en -= 100000.0;
	  }
	  if(rag1>=r0 && rag1<=r1){
	    if(interaction[i][i2]==1){
	      en += e0;
	    }
	  }
	}
	
      }
      
    }
  }


  asc=ti*(en-ev);
  a1=drand48();

  
  if(log(a1)<asc){
    for(i=0;i<L;i++){
      for(i1=0;i1<3;i1++){
	z[i][i1]=zn[i][i1];
      }
    }
  }

  if(log(a1)>=asc){
    en=ev;
  }

 
}


void streamfile(int step)
{
  int i;
  FILE *output,*output2,*output3;

if((output = fopen("configuration.dat","w"))==NULL)
{
 printf("cannot open configuration output file\n");
 exit(0);
}

if((output2 = fopen("polymerases.dat","w"))==NULL)
{
 printf("cannot open polymerase output file\n");
 exit(0);
}

if((output3 = fopen("distancecheck.dat","w"))==NULL)
{
 printf("cannot open distance check file\n");

 exit(0);
}

 for(i=0;i<L;i++){

   fprintf(output,"%lf %lf %lf\n",z[i][0],z[i][1],z[i][2]);
   
   if(polymerase[i]==1) fprintf(output2,"%lf %lf %lf\n",z[i][0],z[i][1],z[i][2]);   
   if(i!=L-1){
     fprintf(output3,"%d %lf\n",i,sqrt((z[i+1][0]-z[i][0])*(z[i+1][0]-z[i][0])+ (z[i+1][1]-z[i][1])*(z[i+1][1]-z[i][1])+ (z[i+1][2]-z[i][2])*(z[i+1][2]-z[i][2])));
   }	     

   fflush(output);
   fflush(output2);
   fflush(output3);

 }

 fclose(output);
 fclose(output2);
 fclose(output3);

}
