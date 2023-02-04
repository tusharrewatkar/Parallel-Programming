// waters.c
#include <stdio.h>
#include <math.h>
#include <time.h> // for CPU time
#include <sys/time.h> //for gettimeofday
#include <mpi.h>

#define LENGTH 80
// global variables
const int maxnum=100000;
double r[maxnum][3][3],rcutsq=1.44,L;
// r(number of molecule, atom 0=O,1=H,2=H, coordinate 0=x,1=y,2=z)

double sqr(double a){return a*a;}

double energy12(int i1,int i2){
// ============================
  int m,n,xyz;
  double shift[3],dr[3],mn[3],r6,distsq,dist,ene=0;
  const double sig=0.3166,eps=0.65,eps0=8.85e-12,e=1.602e-19,Na=6.022e23,q[3]={-0.8476,0.4238,0.4238};
  double elst,sig6;
  elst=e*e/(4*3.141593*eps0*1e-9)*Na/1e3,sig6=pow(sig,6);

  // periodic boundary conditions
  for(xyz=0;xyz<=2;xyz++){
   dr[xyz]=r[i1][0][xyz]-r[i2][0][xyz];shift[xyz]=-L*floor(dr[xyz]/L+.5); //round dr[xyz]/L to nearest integer
   dr[xyz]=dr[xyz]+shift[xyz];
  }
  distsq=sqr(dr[0])+sqr(dr[1])+sqr(dr[2]);
  if(distsq<rcutsq){ // calculate energy if within cutoff
    r6=sig6/pow(distsq,3);
    ene=4*eps*r6*(r6-1.); // LJ energy
    for(m=0;m<=2;m++){
      for(n=0;n<=2;n++){
        for(xyz=0;xyz<=2;xyz++) mn[xyz]=r[i1][m][xyz]-r[i2][n][xyz]+shift[xyz];
        dist=sqrt(sqr(mn[0])+sqr(mn[1])+sqr(mn[2]));
        ene=ene+elst*q[m]*q[n]/dist;
    } }
  }
  return ene;
}

main(int argc, char *argv[]){
  int i,j,natoms,nmol,nproc,me;
  double energy=0,dtime,total,count=0,total_count,minnpair,total_pair,maxnpair,load_imbalance,before_dtime,before1_dtime,wall_time_of_reading,cpu_time_of_reading,total_wall_job_time,total_cpu_job_time;
  double wall_time_max_reading,cpu_time_max_reading,wall_time_min_reading,cpu_time_min_reading,wall_time_max_calculation;
  double wall_time_max_total,cpu_time_max_total,wall_time_min_total,cpu_time_min_total,wall_time_of_calulation,cpu_time_of_calculation;
  double cpu_t1,cpu_t2,cpu_t3;
  FILE *fp;
  char line[LENGTH],nothing[LENGTH],name[20];
  clock_t cputime,t1,t2,t3; /* clock_t defined in <time.h> and <sys/types.h> as int */
  struct timeval start, end;

  
  MPI_Init(&argc,&argv); // initialise MPI
  MPI_Comm_size(MPI_COMM_WORLD,&nproc); // return total number of processors
  MPI_Comm_rank(MPI_COMM_WORLD,&me); // return number of this processor, me=0..nproc-1
  
t1= clock();      // Clock and wall time before Configuration
gettimeofday(&start, NULL); // returns structure with time in s and us (microseconds)
gettimeofday(&end, NULL);
before1_dtime = ((end.tv_sec  - start.tv_sec)+(end.tv_usec - start.tv_usec)/1e6);
cpu_t1=(float)t1/CLOCKS_PER_SEC;

 if (me==0)
 {
  printf("*************(Option 2)***************");
  printf("\n*********For %i number of threads********\n",nproc);
  printf("Input NAME of configuration file\n");
  scanf("%s",name); // reading of filename from keyboard 
  printf("Program to calculate energy of water\n");
 }
 MPI_Bcast(&name,20,MPI_CHAR,0,MPI_COMM_WORLD);
  
 for(int k=0;k<nproc;k++){
	if(me==k)
		{
   
			fp=fopen(name, "r"); //opening of file and beginning of reading from HDD
			fgets(line, LENGTH,fp); //skip first line
			fgets(line, LENGTH,fp); sscanf(line,"%i",&natoms);
			nmol=natoms/3; 
			printf("Number of molecules %i\n",nmol);
  
		for (i=0;i<nmol;i++){
			for(j=0;j<=2;j++){
				fgets(line, LENGTH,fp);
		sscanf(line, "%s %s %s %lf %lf %lf",nothing,nothing,nothing, &r[i][j][0],&r[i][j][1],&r[i][j][2]);
			} }
		printf("first line %lf %lf %lf\n",r[0][0][0],r[0][0][1],r[0][0][2]);
		fscanf(fp, "%lf",&L); // read box size
		printf("Box size %lf\n",L);
		fclose(fp);
 }
 MPI_Barrier(MPI_COMM_WORLD);

 }
  
 
//btime = clock();    // assign initial CPU time (IN CPU CLOCKS)
//gettimeofday(&start, NULL); // returns structure with time in s and us (microseconds)
 
//MPI_Bcast(&nmol,1,MPI_INT,0,MPI_COMM_WORLD); //broadcast the input from thread 0 to all threads
//MPI_Bcast(&r,nmol*9,MPI_DOUBLE,0,MPI_COMM_WORLD);
//MPI_Bcast(&L,1,MPI_DOUBLE,0,MPI_COMM_WORLD); //broadcast the input from thread 0 to all threads
  
  
t2= clock()-t1;      // Clock and wall time after configuration and before energy calculation
gettimeofday(&end, NULL);
before_dtime = ((end.tv_sec  - start.tv_sec)+(end.tv_usec - start.tv_usec)/1e6);
cpu_t2=(float)t2/CLOCKS_PER_SEC;


  for(i=me;i<nmol-1;i=nproc+i){ // calculate energy as sum over all pairs
    for(j=i+1;j<nmol;j++) {
		energy=energy+energy12(i,j);
		count++;
	}
	
  }

  MPI_Allreduce(&energy,&total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  
  t3= clock()-t1;      // clock and wall time after energy calculation
  gettimeofday(&end, NULL);
  dtime = ((end.tv_sec  - start.tv_sec)+(end.tv_usec - start.tv_usec)/1e6);
  cpu_t3=(float)t3/CLOCKS_PER_SEC;

  MPI_Allreduce(&count,&total_pair,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&count,&minnpair,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&count,&maxnpair,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD); 

  
  total_count=(nmol*(double)(nmol-1)/2);
  load_imbalance = ((maxnpair-minnpair)/((minnpair+maxnpair)/2));

  
  wall_time_of_reading=(before_dtime-before1_dtime);
  cpu_time_of_reading= (cpu_t2 - cpu_t1);
  
  wall_time_of_calulation=(dtime-before_dtime);
  cpu_time_of_calculation= (cpu_t3 - cpu_t2);
  
  total_wall_job_time=(dtime-before1_dtime);
  total_cpu_job_time= (cpu_t3 - cpu_t1);
  
  
  MPI_Allreduce(&wall_time_of_reading,&wall_time_max_reading,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&total_wall_job_time,&wall_time_max_total,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&wall_time_of_calulation,&wall_time_max_calculation,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  

  
 
  
 // MPI_Allreduce(&wall_time_of_reading,&wall_time_min_reading,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
 // MPI_Allreduce(&cpu_time_of_reading,&cpu_time_min_reading,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
 // MPI_Allreduce(&total_wall_job_time,&wall_time_min_total,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
 // MPI_Allreduce(&total_cpu_job_time,&cpu_time_min_total,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
 // MPI_Allreduce(&cpu_time_of_reading,&cpu_time_max_reading,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
 // MPI_Allreduce(&total_cpu_job_time,&cpu_time_max_total,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);


  
  
  if(me==0)
  {
	  
  printf("\nTotal energy: %lf\n ",total);
  printf("Energy per molecule: %lf\n ",total/nmol);
  printf("Load_imbalance:%lf\n",load_imbalance);
  
  
  printf("\nMax Wall Time of Reading:%lf\n", wall_time_max_reading);
  printf("Max Wall Time of Calculation:%lf\n", wall_time_max_calculation );
  printf("Max Wall Total Job Time: %lf\n", wall_time_max_total);
 
// printf("Max CPU Time of Reading: %lf\n", cpu_time_max_reading );
// printf("Max Total CPU Time since configuration: %lf\n", cpu_time_max_total );
// printf("Min Wall Time of Reading: %lf\n", wall_time_min_reading);
// printf("Min CPU Time of Reading: %lf\n", cpu_time_min_reading );
// printf("Min Total Wall Time since configuration: %lf\n", wall_time_min_total);
// printf("Min Total CPU Time since configuration: %lf\n", cpu_time_min_total );
  
  printf("\nTotal Number of pairs calculated by formula:%lf\n",total_count);
  printf("Actual Number of pairs calculated by Reduce:%lf\n", total_pair);

  printf("Minimum Number of pairs among all threads:%lf\n",minnpair);
  printf("Maximum Number of pairs among all threads:%lf\n",maxnpair);
  printf("\n");

  }
 
  printf("Number of pairs calculated by thread No %i is : %lf\n",me, count);
  
  MPI_Finalize();

}
