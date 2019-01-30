/*
 * portable.h -- Isolate system dependencies for code portability
 */

/*
 * Copyright 1989 by John Kececioglu
 */


#ifndef PortableInclude
#define PortableInclude

#define Ansi
   /*
    * Should not be defined if the compiler does not support the ANSI standard
    */


/*
 * Argument lists for function prototypes
 */
#ifdef Ansi
#define Proto(arguments) arguments
#else
#define Proto(arguments) ()
#endif


/*
 * Void data type
 */
#ifdef Ansi
typedef void Void;
#else
typedef int Void;
#endif


/*
 * Generic pointer data type
 */
#ifdef Ansi
typedef void *Pointer;
#else
typedef char *Pointer;
#endif

/*
 * Dynamic memory allocation
 */
#ifdef Ansi
#include <stdlib.h>
#define MallocArgumentType size_t
#else
#define MallocArgumentType unsigned int
#endif

#define Allocate(bytes) malloc((MallocArgumentType) (bytes))
#define Free(memory)    free((Pointer) (memory))


/*
 * Halting the program
 */
#define Halt() exit(1)


/*
 * Machine constants
 */
#define LengthLongInteger  32
#define MaximumLongInteger ((long) 2147483647)


/*
 * Random number generation
 */
#define MaximumRandomInteger         MaximumLongInteger
#define GenerateRandomInteger()      ((long) lrand48())
#define SeedRandomIntegerGenerator() srand48((long) 1)


typedef float distancevalue; // the type of number returned by the distance function

#ifndef MAX
#define MAX(A,B) ( (A) > (B) ? (A) : (B) )
#endif


// this is what's kept in the heap-- both a lower bound distance between a point and a subtree
// and the actual distance between the point and that subtree's center
typedef struct {
    distancevalue lowerboundDistance;
    distancevalue centerdistance;
} CompoundKey;


#endif /* PortableInclude */
