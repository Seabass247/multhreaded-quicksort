#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h> 
#include <pthread.h>
#include <errno.h>

typedef int int32;

typedef char uint8;

int32 g_SIZE;
int32 g_THRESHOLD;
int32 g_SEED = 0;
char g_MULTITHREAD;
int32 g_PIECES = 10;
int32 g_MAXTHREADS = 4;

int pthread_tryjoin_np(); // NOTE: for getting rid of compiler warning

void debug_print_array(char* name, int32* array, int32 n)
{
    int32 loop;

    printf("%s={ ", name);
    for(loop = 0; loop < n; loop++)
        printf("%d ", array[loop]);
    printf("}");
}

void scramble_array(int32 *array, int32 n)
{
    if(n > 1) {
        for(int i = n - 1; i > 0; i-- ) {
            int rand_i = (int)(rand() % n);
            int tmp_j = array[rand_i];
            array[rand_i] = array[i];
            array[i] = tmp_j;
        }
    }
}

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
} Segment;

int Partition(int32* array, int32 lo, int32 hi) {

    int32 pivot = (rand() % (hi - lo + 1)) + lo; // generate a random number in the range [lo, hi]
    int32 i = lo - 1;
    int32 j = hi + 1;
    int32 temp;

    do {
        do i++; while(array[i] <= array[pivot]);
        do j--; while(array[j] > array[pivot]);

        if(i < j) {
            temp = array[i];
            array[i] = array[j];
            array[j] = temp;

            if(pivot == j) { // if we swapped the pivot, keep track of its new position
                pivot = i;
            } else if(pivot == i) {
                pivot = j;
            }
        } else break;
    } while(true);
    
    temp = array[pivot];
    array[pivot] = array[j];
    array[j] = temp;

    return j;
}

void GetPartitionSegments(int32* array, int32 lo, int32 hi, Segment* seg_lo, Segment* seg_hi) {
    //printf("\nPartitioning %d -> %d (%d)...", lo, hi, hi - lo + 1);
    int32 j = Partition(array, lo, hi);
    seg_lo->lo = lo;
    seg_lo->hi = (j-1);
    seg_lo->size = j - lo;

    seg_hi->lo = j+1;
    seg_hi->hi = hi;
    seg_hi->size = hi - j -1 + 1;

    int32 total_size = seg_lo->size + seg_hi->size;
    double lo_percent =((double)seg_lo->size / (double)total_size) * 100.0;
    double hi_percent = ((double)seg_hi->size / (double)total_size) * 100.0;
    //printf("result: %d - %d (%f / %f)", seg_lo->size, seg_hi->size, lo_percent, hi_percent);
}

void QuickSort(int32* array, int32 lo, int32 hi) {         

    if (lo > hi) // nothing to do according to algorithm
        return;

    int32 segment_size = hi - lo + 1;
    int32* segment =  &array[lo];
    
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
        // if the segment size is at or below the threshold, shellsort it
        int32 i;
        int32 j;
        int32 k = 1;
        while (k <= segment_size) k *= 2;
        k = (k / 2) - 1;

        do {
            for(i = 0; i < (segment_size - k); i++) {
                for(j = i; j >= 0; j -= k) {
                    if(segment[j] <= segment[j + k]) break;
                    else {
                        int32 tmp = segment[j];
                        segment[j] = segment[j + k];
                        segment[j + k] = tmp;
                    }
                }
            }

            k /= 2;
        } while (k > 0);
        return;
    }

    if(lo < hi) {
        int32 j = Partition(array, lo, hi);
        QuickSort(array, lo, j-1);
        QuickSort(array, j+1, hi);
    }
}

void SortSegments(Segment array[], int32 n) {
    // sort segments by size in descending order
    // using the bubble sort algorithm
    int32 i, j; 
    for (i = 0; i < n-1; i++) {
        for (j = 0; j < n-i-1; j++) {
            if (array[j].size < array[j+1].size) {
                Segment temp = array[j]; 
                array[j] = array[j+1]; 
                array[j+1] = temp; 

            }            
        }
    }
}

typedef struct {
    int32* array;
    Segment seg;
} Parameters;

void* sortThread(void* ptr) {
    Parameters* params;
    params = (Parameters*)ptr;
    QuickSort(params->array, 0, params->seg.size - 1);
}

int main(int argc, char *argv[])
{
    char* resultPtr; // An pointer that would be NULL for failed number conversions in `strtoull` (but we're not going to check)
    
    int num_args = argc - 1;
    if(num_args < 2) { // check for the command-line parameter: n0
        fprintf(stderr, "Error: not enough args\n");
        exit(1);
    }

    g_SIZE = strtoull(argv[1], &resultPtr, 10); // convert the string of the input number to an int
    if(resultPtr == NULL) {
        fprintf(stderr, "Error: failed to convert argument to integer\n");
        exit(1);
    }
    //printf("SIZE=%d\n", g_arraySize);

    g_THRESHOLD = strtoull(argv[2], &resultPtr, 10); // convert the string of the input number to an int
    if(resultPtr == NULL) {
        fprintf(stderr, "Error: failed to convert argument to integer\n");
        exit(1);
    }
    //printf("THRESHOLD=%d\n", threshold_arg);

    if(num_args >= 3 && num_args <= 6) {
        g_SEED = strtoull(argv[3], &resultPtr, 10); // convert the string of the input number to an int
        if(resultPtr == NULL) {
            fprintf(stderr, "Error: failed to convert argument to integer\n");
            exit(1);
        }
        //printf("SEED=%d\n", g_SEED);

        if(num_args >= 4) {
            g_MULTITHREAD = argv[4][0]; // convert the string of the input number to an int
            if(resultPtr == NULL) {
                fprintf(stderr, "Error: failed to convert argument to integer\n");
                exit(1);
            }
            //printf("MULTITHREAD=%s\n", g_MULTITHREAD);

            if(num_args >= 5) {
                g_PIECES = strtoull(argv[5], &resultPtr, 10); // convert the string of the input number to an int
                if(resultPtr == NULL) {
                    fprintf(stderr, "Error: failed to convert argument to integer\n");
                    exit(1);
                }
                //printf("PIECES=%d\n", g_PIECES);

                if(num_args >= 6) {
                    g_MAXTHREADS = strtoull(argv[6], &resultPtr, 10); // convert the string of the input number to an int
                    if(resultPtr == NULL) {
                        fprintf(stderr, "Error: failed to convert argument to integer\n");
                        exit(1);
                    }

                    printf("\nParsed MAXTHREADS=%d\n", g_MAXTHREADS);

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

    int32* array = (int32*)malloc(g_SIZE * sizeof(int32)); // init array
    for(int32 i = 0; i < g_SIZE; i++) { // create array
        array[i] = i;
    }
    scramble_array(array, g_SIZE); // scramble array

    if(g_MULTITHREAD == 'n') { // Do single-threaded sorting
        //printf("\nDo single-threaded sort\n");

        //debug_print_array("unsorted_array", array, g_arraySize);

        QuickSort(array, 0, g_SIZE - 1);

        //debug_print_array("sorted_array", array, g_arraySize);

        if(isSorted(array, g_SIZE)) {
            printf("\nArray is sorted\n");
        } else {
            printf("\n\nArray is NOT sorted!\n");
        }
    } else { //  Do multi-threaded sorting
        //printf("\nDo multi-threaded sort\n");

        int32 pieces_count = 0;

        printf("\nSort %d pieces with %d threads\n", g_PIECES, g_MAXTHREADS);

        Segment segments[g_PIECES];
        Segment seg_hi = {0};
        Segment seg_lo = {0};

        Segment next_biggest_seg = {0};
        next_biggest_seg.lo = 0;
        next_biggest_seg.hi = g_SIZE - 1;
        next_biggest_seg.size = g_SIZE;

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
                    //printf("Lower segment is bigger, next[%d]: %d -> %d\n", pieces_count, seg_lo.lo, seg_lo.hi);
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
                    //printf("Upper segment is bigger, next[%d]: %d -> %d\n", pieces_count, seg_hi.lo, seg_hi.hi);
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

        // printf("\nTotal pieces=%d\n", pieces_count);
        // for(int32 i = 0; i < pieces_count; i++) {
        //     printf("%d -> %d, (%d)\n", segments[i].lo, segments[i].hi, segments[i].size);
        // }

        pthread_t thread_ids[g_MAXTHREADS];
        int32 threads_running = 0;
        // DEBUG: simulate threads iteratively (TODO: replace with multithreading)
        for(int32 i = 0; i < pieces_count; i++) {
            //QuickSort(&array[segments[i].lo], 0, segments[i].size - 1);

            Parameters *params = (Parameters *) malloc(sizeof(Parameters)); 
            params->array = &array[segments[i].lo];
            params->seg = segments[i];

            if(threads_running < g_MAXTHREADS)
                pthread_create(&thread_ids[threads_running++], NULL, sortThread, (void*)params);
            // if(isSorted(&array[segments[i].lo], segments[i].size)) {
            //     printf("\n(%d)Array is sorted\n", i);
            // } else {
            //     printf("\n\n(%d)Array is NOT sorted!\n", i);
            // }
            
        }

        int32 seg_count = threads_running;
        while(threads_running) {
            for(int32 i = 0; i < g_MAXTHREADS; i++) {
                if(pthread_tryjoin_np(thread_ids[i], NULL) == 0) {
                    threads_running--; // Thread i finished...
                    //printf("Thread %d finished...\n", i);
                    if(seg_count < g_PIECES) { // if we are here...
                        // ...then seg_count segments have already begun to be sorted in threads, but
                        // there is a total g_PIECES to sort: start another sort on an available thread
                        //printf("%d segment sorts started but %d total: start new sort on thread %d\n", seg_count, g_PIECES, i);

                        Parameters *params = (Parameters *) malloc(sizeof(Parameters)); 
                        params->array = &array[segments[seg_count].lo];
                        params->seg = segments[seg_count];
                        printf("(%d, %d, %d)\n", params->seg.lo, params->seg.hi, params->seg.size);
                        pthread_create(&thread_ids[i], NULL, sortThread, (void*)params);

                        threads_running++;
                        seg_count++;
                    }
                }
            }

        }

        if(isSorted(array, g_SIZE)) {
            printf("\nEntire array is sorted\n");
        } else {
            printf("\n\nEntire array is NOT sorted!\n");
        }
    }

    return 0;
} // end main
