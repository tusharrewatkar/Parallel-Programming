// waters.c
#include <stdio.h>
#include <math.h>
#include <time.h>     // for CPU time
#include <sys/time.h> //for gettimeofday
#include <omp.h>

#define LENGTH 80
// global variables
const int maxnum = 1000000;
double r[maxnum][3][3], rcutsq = 1.44, L;
// r(number of molecule, atom 0=O,1=H,2=H, coordinate 0=x,1=y,2=z)

double sqr(double a) { return a * a; }

double energy12(int i1, int i2)
{
    // ============================
    int m, n, xyz;
    double shift[3], dr[3], mn[3], r6, distsq, dist, ene = 0;
    const double sig = 0.3166, eps = 0.65, eps0 = 8.85e-12, e = 1.602e-19, Na = 6.022e23, q[3] = {-0.8476, 0.4238, 0.4238};
    double elst, sig6;
    elst = e * e / (4 * 3.141593 * eps0 * 1e-9) * Na / 1e3, sig6 = pow(sig, 6);

    // periodic boundary conditions
    for (xyz = 0; xyz <= 2; xyz++)
    {
        dr[xyz] = r[i1][0][xyz] - r[i2][0][xyz];
        shift[xyz] = -L * floor(dr[xyz] / L + .5); // round dr[xyz]/L to nearest integer
        dr[xyz] = dr[xyz] + shift[xyz];
    }
    distsq = sqr(dr[0]) + sqr(dr[1]) + sqr(dr[2]);
    if (distsq < rcutsq)
    { // calculate energy if within cutoff
        r6 = sig6 / pow(distsq, 3);
        ene = 4 * eps * r6 * (r6 - 1.); // LJ energy
        for (m = 0; m <= 2; m++)
        {
            for (n = 0; n <= 2; n++)
            {
                for (xyz = 0; xyz <= 2; xyz++)
                    mn[xyz] = r[i1][m][xyz] - r[i2][n][xyz] + shift[xyz];
                dist = sqrt(sqr(mn[0]) + sqr(mn[1]) + sqr(mn[2]));
                ene = ene + elst * q[m] * q[n] / dist;
            }
        }
    }
    return ene;
}

// To calculate Max value of an array

double maxValue(double myArr[], int threads)

{
    double max = myArr[0];
    for (int m = 0; m < threads; m++)
    {
        if (myArr[m] > max)
        {
            max = myArr[m];
        }
    }
    return max;
}

// To calculate Min value of an array

double minValue(double myArr[], int threads)
{
    double min = myArr[0];
    for (int n = 0; n < threads; n++)
    {
        if (myArr[n] < min)
        {
            min = myArr[n];
        }
    }
    return min;
}

main()
{
    int i, j, natoms, nmol, threads = 12;
    double energy = 0, dtime, minnpair, maxnpair, count;
    double npair_shared[threads];
    FILE *fp;
    char line[LENGTH], nothing[LENGTH], name[20];
    clock_t cputime; /* clock_t defined in <time.h> and <sys/types.h> as int */
    struct timeval start, end;

    printf("\n**********Energy calculation OPEN MP***********\n");
    printf("Program to calculate energy of water\n");
    printf("Input NAME of configuration file\n");
    scanf("%s", name);       // reading of filename from keyboard
    fp = fopen(name, "r");   // opening of file and beginning of reading from HDD
    fgets(line, LENGTH, fp); // skip first line
    fgets(line, LENGTH, fp);
    sscanf(line, "%i", &natoms);
    nmol = natoms / 3;
    printf("Number of molecules %i\n", nmol);

    for (int i = 0; i < nmol; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            fgets(line, LENGTH, fp);
            sscanf(line, "%s %s %s %lf %lf %lf", nothing, nothing, nothing, &r[i][j][0], &r[i][j][1], &r[i][j][2]);
        }
    }
    printf("first line %lf %lf %lf\n", r[0][0][0], r[0][0][1], r[0][0][2]);
    fscanf(fp, "%lf", &L); // read box size
    printf("Box size %lf\n", L);

    cputime = clock();          // assign initial CPU time (IN CPU CLOCKS)
    gettimeofday(&start, NULL); // returns structure with time in s and us (microseconds)

    omp_set_num_threads(threads);

#pragma omp parallel
    {
        count = 0;
        int me = omp_get_thread_num();

        // npair_shared[me] = 0;

#pragma omp for private(i,count) reduction(+ : energy) schedule(static, 1)
        for (i = 0; i < nmol - 1; i++)
        { // calculate energy as sum over all pairs
            for (j = i + 1; j < nmol; j++)
            {
                energy = energy + energy12(i, j);
                count++;
            }
            npair_shared[me] = count;
        }

        printf("Number of pairs calculated by thread No %i is : %lf \n", me, npair_shared[me]);
    }

    cputime = clock() - cputime; // calculate  cpu clock time as difference of times after-before
    gettimeofday(&end, NULL);
    dtime = ((end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6);

    maxnpair = maxValue(npair_shared, threads);
    minnpair = minValue(npair_shared, threads);

    printf("Minimum Number of pairs among all threads:%lf \n", minnpair);
    printf("Maximum Number of pairs among all threads:%lf \n", maxnpair);

    double load_imbalance = (maxnpair - minnpair) / ((minnpair + maxnpair) / 2);

    printf("Load_Imbalance:%lf\n", load_imbalance);
    printf("Total energy %lf \n ", energy);
    printf("Energy per molecule %lf \n", energy / nmol);
    printf("Elapsed wall time: %f\n", dtime);
    printf("Elapsed CPU  time: %f\n", (float)cputime / CLOCKS_PER_SEC);
    fclose(fp);
}
