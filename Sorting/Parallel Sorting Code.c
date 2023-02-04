#include <time.h> 
#include <sys/time.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>


#define max_len 100000
#define LENGTH 40

//Function for Calculating Maxmimum Value of array

double max_value(double arr[], int n) { 

    double max = arr[1];     // Initialize maximum element
    for (int i = 1; i <= n; i++)
        if (arr[i] > max)
            max = arr[i];

    return max;
}

//Function for Calculating Minimum Value of array

double min_value(double arr[], int n) {  

    double min = arr[1];     // Initialize maximum element
    for (int i = 1; i <= n; i++)
        if (arr[i] < min)
            min = arr[i];

    return min;
}

main(int argc, char *argv[]) {

    int i=1,len,ind[max_len+1],j,cur,prev,k;
    double b[max_len+1],c[max_len+1],new,bnew,cnew,time;
    char name[LENGTH],line[LENGTH];
    FILE *fp;
    clock_t cpu0,cpu1,cpu2,cpu3, cpu4, cpu5, cpu_max, cpu2_max; // clock_t defined in <time.h> and <sys/types.h> as int
    struct timeval time0, time1,time2,time3; // for wall clock in s and us
    double dtime12,dtime03, dtime_max, dtime2_max; // for wall clock in s (real number)



    double highest_value, lowest_value, diff, sub_low, sub_high;
    int count = 1,maxnpair,minnpair;
    double b_sub[max_len+1], c_sub[max_len+1];\
    int load_imbalance;

    int me,nproc;

    //Clearing old values from sorted-p.txt OR trucating old values length to 0
    int fd = open("sorted-p.txt", O_WRONLY);
    ftruncate(fd, 0);
    close(fd);



    MPI_Init(&argc,&argv); // initialise MPI
    MPI_Comm_size(MPI_COMM_WORLD,&nproc); // return total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD,&me); // return number of this processor, me=0..nproc-1


    cpu0 = clock();    // assign initial CPU time (IN CPU CLOCKS)
    gettimeofday(&time0, NULL); // returns structure with time in s and us (microseconds)

  



    if(me==0) {
        fp=fopen("100k.txt","r");
        while(1) { //1 serves as true, i.e. condition which is always true
            if(fgets(line, LENGTH,fp)==NULL) 
                break; // finish reading when empty line is read
            if(sscanf(line, "%lf %lf",&b[i],&c[i])==-1) 
                break; // finish reading after error
    
            i++;
        }

        len=i-1;
        fclose(fp);

        highest_value = max_value(b, len);
        lowest_value = min_value(b, len);
        diff = highest_value - lowest_value;
        
        printf("\n*******************  < Parallel Sorting for %i threads> ******************\n",nproc);
        printf("Total Number of items to sort in text file: %i\n", len);
        printf("Highest value in text file  : %lf\n", highest_value);
        printf("Lowest number in text file : %lf\n", lowest_value);
        printf("*****************  < Dividing our array into %i parts >  *******************\n",nproc);


    }

    MPI_Bcast(&b,max_len+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&c,max_len+1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&len,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&highest_value,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&lowest_value,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&diff,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    cpu1 = clock();    // assign initial CPU time (IN CPU CLOCKS)
    gettimeofday(&time1, NULL); // returns structure with time in s and us (microseconds)


    sub_low = (me * (diff/nproc))  + lowest_value;
    sub_high = (me * (diff/nproc)) + diff/nproc + lowest_value;
    
    for(j = 1; j <= len; j = j + 1) {
        if ((b[j] >= sub_low) & (b[j] <= sub_high)) {

            b_sub[count] = b[j];
            c_sub[count] = c[j];
            count++;
        }
    }
    count = count - 1;

    printf("Thread %i has %i numbers calculates from low value = %lf to high value = %lf \n", me,count, sub_low, sub_high);

    //Sort the numbers for each thread

    ind[0] = 1;
    for(j = 2; j <= count; j++) { // start sorting with the second item

        new = b_sub[j];
        cnew = c_sub[j];
        cur = 0;
        
        for(i = 1; i < j; i++) {

            prev = cur; 
            cur = ind[cur];

            if(new == b_sub[cur]) {
                printf("Equal numbers are: %lf\n", new);
            }

            if( (new < b_sub[cur]) | ( (new == b_sub[cur]) & (cnew < c_sub[cur]) ) ) {
                ind[prev] = j;
                ind[j] = cur;
                goto loop;
            }
        }
        // new number is the largest so far
        ind[cur] = j;
        loop: ;
    }

    cpu2 = clock();    // assign CPU time (IN CPU CLOCKS)
    gettimeofday(&time2, NULL);
    dtime12 = ((time2.tv_sec  - time1.tv_sec)+(time2.tv_usec - time1.tv_usec)/1e6);



    // //Write sorted values of each thread to text file
    for(int c=0; c < nproc; c++) {
        if(me == c) { 
            cur = 0;
            fp = fopen("sorted-p.txt","a");
            for(i = 0; i < count; i++) {
                cur = ind[cur];
                fprintf(fp,"%lf %lf\n", b_sub[cur], c_sub[cur]);
            }
            
            fclose(fp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    

    cpu3 = clock();    // assign initial CPU time (IN CPU CLOCKS)
    gettimeofday(&time3, NULL); // returns structure with time in s and us (microseconds)
    dtime03 = ((time3.tv_sec  - time0.tv_sec)+(time3.tv_usec - time0.tv_usec)/1e6);



    //For load Imbalance calculation
    MPI_Allreduce(&count,&minnpair,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&count,&maxnpair,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD); 

    
    //Get the maximum of all thread times using mpi reduce
    MPI_Reduce(&dtime12,&dtime_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&dtime03,&dtime2_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);


    if(me==0) {

        printf("Load_imbalance:%lf\n",(double)(2*(maxnpair-minnpair)/(0.+maxnpair+maxnpair)));
        printf("Elapsed Wall time for Sorting : %f\n", dtime_max);
        printf("Elapsed Total Job time: %f\n", dtime2_max);
    }

    MPI_Finalize(); // finalize MPI peacefully (the system would kill the processes otherwise)

}