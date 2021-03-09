#define _GNU_SOURCE
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h> 
#include <pthread.h>
#include <errno.h>
#include <sys/time.h>
typedef int int32;

typedef char uint8;

int32 g_SIZE;
int32 g_THRESHOLD;
int32 g_SEED = 0;
char g_MULTITHREAD;
int32 g_PIECES = 10;
int32 g_MAXTHREADS = 4;

void _debug_print_array(char* name, int32* array, int32 n)
{
    int32 loop;

    printf("%s={ ", name);
    for(loop = 0; loop < n; loop++)
        printf("%d ", array[loop]);
    printf("}");
}

// swaps each value in the array with another value at a random location
void scramble_array(int32 *array, int32 n)
{
    if(n > 1) {
        for(int i = n - 1; i > 0; i-- ) {
            int rand_i = (int)(rand() % n);
            int tmp_j = array[rand_i]; // swap value at rand_i with value at i
            array[rand_i] = array[i];
            array[i] = tmp_j;
        }
    }
}

// checks if the array is sorted in ascending order
bool isSorted(int32* array, int32 n) {
    int32 i = 0;
    bool isSorted = true;
    
    while(i < n - 1) {
        if(array[i] >= array[i+1])
            isSorted = false;
        i++;
    }

    return isSorted;
}

typedef struct {
    int32 lo;
    int32 hi;
    int32 size;
} Segment; // a structure created for conveniently passing segment data around

// Partitions the array into two segments, first segment being [lo to j-1], second being [j+1 to hi]
// where j is the returned value
int Partition(int32* array, int32 lo, int32 hi) {
    
    int32 pivot = (rand() % (hi - lo + 1)) + lo; // generate a random number in the range [lo, hi]
    int32 i = lo - 1;
    int32 j = hi + 1;
    int32 temp;

    do {
        do i++; while(array[i] <= array[pivot]); // move i to the first number that's greater than pivot
        do j--; while(array[j] > array[pivot]); // move j to the first number that's less than or equal to pivot

        if(i < j) {
            temp = array[i];
            array[i] = array[j];
            array[j] = temp;

            if(pivot == j) { // if we swapped the pivot, keep track of its new position
                pivot = i;
            } else if(pivot == i) {
                pivot = j;
            }
        } else break; // stop when i and j have crossed
    } while(true);
    
    temp = array[pivot]; // swap pivot and j
    array[pivot] = array[j];
    array[j] = temp;

    return j;
}

// Built on top of Partition(), this function partitions the array, then builds
// a lower segment object and upper segment object from a partition from lo to hi
void GetPartitionSegments(int32* array, int32 lo, int32 hi, Segment* seg_lo, Segment* seg_hi) {
    int32 j = Partition(array, lo, hi);
    seg_lo->lo = lo; // lower segment starts at lo
    seg_lo->hi = (j-1); // higher segments ends at j-1
    seg_lo->size = j - lo; 

    seg_hi->lo = j+1; // upper segment starts at j+1
    seg_hi->hi = hi; // upper segment ends at hi
    seg_hi->size = hi - j -1 + 1;
}

// Sorts the array using a hybrid quicksort-shellsort algorithm approach- when the
// segment size is below the specified threshold, the remaining segment is shellsorted
void QuickSort(int32* array, int32 lo, int32 hi) {

    if (lo > hi) // nothing to do according to algorithm
        return;

    int32 segment_size = hi - lo + 1; // get the size of this segment
    int32* segment =  &array[lo]; // get a pointer to this segment of the array
    
    if(segment_size < 2) // the segment is 0 or 1 length- nothing to do
        return;

    if(segment_size == 2) { // only 2 values, swap if necessary
        if(segment[0] > segment[1]) {
            int32 tmp = segment[0];
            segment[0] = segment[1];
            segment[1] = tmp;
        }
        return;
    }

    if(segment_size <= g_THRESHOLD) {
        // if the segment size is at or below the threshold, SHELLSORT it
        int32 i;
        int32 j;
        int32 k = 1;
        while (k <= segment_size) k *= 2; // generate the sequence used to shellsort
        k = (k / 2) - 1;

        do {
            for(i = 0; i < (segment_size - k); i++) {
                for(j = i; j >= 0; j -= k) {
                    if(segment[j] <= segment[j + k]) break;
                    else {
                        int32 tmp = segment[j]; // swap j and j+k
                        segment[j] = segment[j + k];
                        segment[j + k] = tmp;
                    }
                }
            }

            k /= 2;
        } while (k > 0);
        return; // that's it, segment is sorted via shellsort, nothing else to do, so return
    }

    if(lo < hi) { // don't continue recursively quicksorting for segments where lo < hi, because...
                // ...according to the quicksort algorithm the array is trivially assorted
        int32 j = Partition(array, lo, hi);
        QuickSort(array, lo, j-1);
        QuickSort(array, j+1, hi);
    }
}

// Sorts segments by size in descending order
// using the bubble sort algorithm
void SortSegments(Segment array[], int32 n) {

    int32 i, j; 
    for (i = 0; i < n-1; i++) {
        for (j = 0; j < n-i-1; j++) {
            if (array[j].size < array[j+1].size) {
                Segment temp = array[j]; // swap j and j+1
                array[j] = array[j+1]; 
                array[j+1] = temp; 

            }            
        }
    }
}

typedef struct {
    int32* array;
    Segment seg;
} Parameters; // a struct used to pass the necessary parameter to a sorting thread

// function which threads will run
void* sortThread(void* ptr) {
    Parameters* params;
    params = (Parameters*)ptr; // get parameters passed into this thread
    QuickSort(params->array, 0, params->seg.size - 1);
}

// entry point to the program
int main(int argc, char *argv[])
{
	// Declare timing variables
	clock_t all_cpu_start, all_cpu_end,
            create_start, create_end,
        	init_start, init_end,
            scramble_start, scramble_end, 
            sortCPU_start, sortCPU_end,
			part_start, part_end;
	struct timeval 	all_wall_start, all_wall_end, sortWALL_start, sortWALL_end;
	double 	create_time, init_time, shuffle_time, part_time,
			sort_wall_total, sort_cpu_total, wall_total, cpu_total;

	all_cpu_start = clock(); // begin measuring total time elapsed for both CPU...
	gettimeofday(&all_wall_start, NULL); // ...and WALL clock time

    char* resultPtr; // An pointer that would be NULL for failed number conversions in `strtoull` (but we're not going to check)
    
    int num_args = argc - 1;
    if(num_args < 2) { // check that there are the minimum number of arguments needed to sort
        fprintf(stderr, "Error: not enough args\n");
        exit(1);
    }

    g_SIZE = atoi(argv[1]); // convert the string of the SIZE number to an int

    g_THRESHOLD = atoi(argv[2]); // convert the string of the THRESHOLD number to an int

    if(num_args >= 3 && num_args <= 6) {
        g_SEED = atoi(argv[3]); // convert the string of the SEED number to an int

        if(num_args >= 4) {
            g_MULTITHREAD = argv[4][0]; // get the character of the MULTITHREAD argument

            if(num_args >= 5) {
                g_PIECES = atoi(argv[5]); // convert the string of the PIECES number to an int


                if(num_args >= 6) {
                    g_MAXTHREADS = atoi(argv[6]); // convert the string of the MAXTHREADS number to an int

                }
            }
        }
    } else if(num_args >= 6) {
        fprintf(stderr, "Error: supplied more than 6 arguments\n");
        exit(1);
    }

    if(g_SEED > 0) // if seed arg is set to any positive integer, use it to seed random number generator
        srand(g_SEED);
    else if (g_SEED == -1) // else if seed is set to -1, seed random number generator with clock
        srand(clock());

    create_start = clock();
    int32* array = (int32*)malloc(g_SIZE * sizeof(int32)); // create array
    create_end = clock();

    init_start = clock();
    for(int32 i = 0; i < g_SIZE; i++) { // initialize array
        array[i] = i;
    }
    init_end = clock();

    scramble_start = clock();
    scramble_array(array, g_SIZE); // scramble array
    scramble_end = clock();

    if(g_MULTITHREAD == 'n' || g_MULTITHREAD == 'N') { // Do single-threaded sorting
        g_PIECES = 0; // don't need these...
        g_MAXTHREADS = 0; //...because we aren't using multi-threading

        sortCPU_start = clock(); // begin measuring..
        gettimeofday(&sortWALL_start, NULL); //..sorting time
        QuickSort(array, 0, g_SIZE - 1);
    } else { //  Do multi-threaded sorting
        int32 pieces_count = 0;

        Segment segments[g_PIECES]; // create a list of PIECES segments
        Segment seg_hi = {0}; 
        Segment seg_lo = {0};

        Segment next_biggest_seg = {0};
        next_biggest_seg.lo = 0;
        next_biggest_seg.hi = g_SIZE - 1;
        next_biggest_seg.size = g_SIZE;

        // Begin partitioning...
        part_start = clock();
        if(g_PIECES == 1) { // if the desired number of pieces is 1, just add one segment...
            //..which represents the entire array
            Segment single = {0};
            single.lo = 0;
            single.hi = g_SIZE - 1;
            single.size = g_SIZE;

            segments[pieces_count++] = single;
        } else { // else, operate normally:
            // Do as many partitions as necessary to produce PIECES pieces
            while(pieces_count < (g_PIECES-1) && next_biggest_seg.size > 0) {
            // printf("\nThis next size=%d", next_biggest_seg.size);
                GetPartitionSegments(array, next_biggest_seg.lo, next_biggest_seg.hi, &seg_lo, &seg_hi);
                if(seg_lo.size >= seg_hi.size) { // lower segment is bigger...
                    //... so partition the lower segment next
                    next_biggest_seg = seg_lo;
                    // and keep the hi (smaller) segment in the list
                    segments[pieces_count++] = seg_hi;

                    if(pieces_count == (g_PIECES - 1)) { // if we have 1 segment left, just add it to the list
                        // if we're here then we have the last remaining segment (and we don't want to partition it),
                        // so put it in the list
                        segments[pieces_count++] = seg_lo;
                    }
                }
                else {  // upper segment is bigger...
                    //... so partition the upper segment next
                    next_biggest_seg = seg_hi;
                    // and keep the lo (smaller) segment in the list
                    segments[pieces_count++] = seg_lo;

                    if(pieces_count == (g_PIECES - 1)) { // if we have 1 segment left, just add it to the list
                        // if we're here then we have the last remaining segment (and we don't want to partition it),
                        // so put it in the list
                        segments[pieces_count++] = seg_hi;
                    }
                }
            }
        }
        SortSegments(segments, pieces_count);
        part_end = clock();

        pthread_t thread_ids[g_MAXTHREADS];
        pthread_attr_t thread_attrs[g_MAXTHREADS];
        int32 threads_running = 0;
        for(int32 i = 0; i < pieces_count; i++) {

            Parameters *params = (Parameters *) malloc(sizeof(Parameters)); 
            params->array = &array[segments[i].lo];
            params->seg = segments[i];

            if(threads_running < g_MAXTHREADS) {
                if(!threads_running) { // start sorting wall clock on first thread creation
                    gettimeofday(&sortWALL_start, NULL);
                    sortCPU_start = clock();
                }

                pthread_attr_init(&thread_attrs[threads_running]);
                printf("(%d, %d, %d)\n", params->seg.lo, params->seg.hi, params->seg.size);
                pthread_create(&thread_ids[threads_running], NULL, sortThread, (void*)params);
                threads_running++; // another thread is now running
            }

        }

        int32 seg_count = threads_running; // used for accessing the next-up segment from the segments array
        while(threads_running) {
            for(int32 i = 0; i < g_MAXTHREADS; i++) {
                bool none_finished = true; // keep track of whether any threads have finished or not
                if(pthread_tryjoin_np(thread_ids[i], NULL) == 0) { // poll each thread
                    none_finished = false;
                    threads_running--; // Thread i finished, that's one less thread running
                    if(seg_count < g_PIECES) { // if we are here...
                        // ...then seg_count segments have already begun to be sorted in threads, but
                        // there is a total g_PIECES to sort: start another sort on an available thread
                        //printf("%d segment sorts started but %d total: start new sort on thread %d\n", seg_count, g_PIECES, i);

                        Parameters *params = (Parameters *) malloc(sizeof(Parameters)); 
                        params->array = &array[segments[seg_count].lo];
                        params->seg = segments[seg_count];
                        pthread_attr_init(&thread_attrs[i]);
                        printf("(%d, %d, %d)\n", params->seg.lo, params->seg.hi, params->seg.size);
                        pthread_create(&thread_ids[i], NULL, sortThread, (void*)params);

                        threads_running++; // another thread is now running
                        seg_count++; // another segment has begun to be sorted
                    }
                }
                if(none_finished) // if none have finished running, wait 50ms before polling again 
                    usleep(50000);
            }

        }

    }

    all_cpu_end = clock();
    sortCPU_end = clock();
    gettimeofday(&sortWALL_end, NULL);
    gettimeofday(&all_wall_end, NULL);

    // Get timing statistics
    create_time = ((double)create_end - (double)create_start)/CLOCKS_PER_SEC;
    init_time = ((double)init_end - (double)init_start)/CLOCKS_PER_SEC;
    shuffle_time = ((double)scramble_end - (double)scramble_start)/CLOCKS_PER_SEC;
    part_time = ((double)part_end - (double)part_start)/CLOCKS_PER_SEC;
    sort_wall_total = ((double)sortWALL_end.tv_sec-(double)sortWALL_start.tv_sec) + ((double)sortWALL_end.tv_usec-(double)sortWALL_start.tv_usec)/1000000;
    sort_cpu_total = ((double)sortCPU_end - (double)sortCPU_start)/CLOCKS_PER_SEC;
    wall_total = ((double)all_wall_end.tv_sec-(double)all_wall_start.tv_sec) + ((double)all_wall_end.tv_usec-(double)all_wall_start.tv_usec)/1000000;
    cpu_total = ((double)all_cpu_end - (double)all_cpu_start)/CLOCKS_PER_SEC;

    // Print nicely formatted timing statistics
    printf("    SIZE    THRESHOLD SD PC T CREATE   INIT  SHUFFLE   PART  SrtWall Srt CPU ALLWall ALL CPU\n");
	printf("  --------- --------- -- -- - ------ ------- ------- ------- ------- ------- ------- -------\n");
    if(g_SEED == 0) {
        // if seed was not specified, format it as '00'
        printf("F:%9d %9d %02d %2d %1d %0.3f  %0.3f   %0.3f   %0.3f   %0.3f   %0.3f   %0.3f   %0.3f\n",g_SIZE,g_THRESHOLD,g_SEED, g_PIECES,
            g_MAXTHREADS, create_time, init_time, shuffle_time, part_time, sort_wall_total, sort_cpu_total, wall_total, cpu_total);
    } else {
        // if seed was not specified, format without leading zeros
        printf("F:%9d %9d %2d %2d %1d %0.3f  %0.3f   %0.3f   %0.3f   %0.3f   %0.3f   %0.3f   %0.3f\n",g_SIZE,g_THRESHOLD,g_SEED, g_PIECES,
            g_MAXTHREADS, create_time, init_time, shuffle_time, part_time, sort_wall_total, sort_cpu_total, wall_total, cpu_total);
    }
	   
    return 0;
} // end main
