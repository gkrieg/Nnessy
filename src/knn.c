#include "metricspaces.h"
#include "dispersionTree.h"
#include "partition.h"



/*
    This file is part of the dispersion tree code.

    The dispersion tree code is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Foobar is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the code for a dispersion tree.  If not, see <http://www.gnu.org/licenses/>.
*/


typedef struct {
    char *spacetype;
    FILE *data;
    FILE *queries;
    unsigned numNeighbors;
    unsigned numTreeLevels;
    char *binaryFileToWrite;
    char *binaryFileToRead;
} Args;

Args
parseArgs(int argc, char **argv) {

    Args args = {NULL, NULL, NULL, 1, 2, NULL, NULL};
    int i, help=0;
    for (i=1; i < argc; ++i) {
	if (argv[i][0] == '-') { 
	    if (argv[i][1] == 's') {
		args.spacetype = argv[i+1]; // e.g., PROTEIN_WORDS (see metricspaces.h)
		++i;
	    } else if (argv[i][1] == 'd') {
		args.data = fopen( argv[ i+1 ], "r" );
		++i;
	    } else if (argv[i][1] == 'q') {
		args.queries = fopen( argv[ i+1 ], "r" );
		++i;
	    } else if (argv[i][1] == 'n') {
		args.numNeighbors = atoi( argv[ i+1 ]);
		++i;
	    } else if (argv[i][1] == 'l') {
		args.numTreeLevels = atoi( argv[ i+1 ]);
		++i;
	    } else if (argv[i][1] == 'w') {
		args.binaryFileToWrite = argv[i+1];
		++i;
	    } else if (argv[i][1] == 'r') {
		args.binaryFileToRead = argv[i+1];
		++i;
	    } else if (argv[i][1] == 'h') {
		help=1;
	    } else {
		fprintf(stderr, "Unknown flag: %s\n", argv[i]);
		Halt();
	    }

	}
    }
    int errors=0;
    if (args.spacetype==NULL) {
	++errors;
	fprintf(stderr, "I need a space type to work!\n");
    }

    if (args.data==NULL) {
	++errors;
	fprintf(stderr, "I need a data file to work!\n");
    }

    if (args.queries==NULL) {
	++errors;
	fprintf(stderr, "I need queries to work!\n");
    }

    if (errors || help) {
	fprintf(stderr, "Correct usage: %s -s spaceType -d dataSetToIndex -q queries\n", argv[0]);
	fprintf(stderr, "Optional arguments:\n\t-l treeStages (default is %u)\n\t-n numNearestNeighbors (default is %u)\n", 
		args.numTreeLevels, args.numNeighbors);
	fprintf(stderr, "\t-w binaryFilename (filename used when writing a dtree to disk)\n\t-r binaryFilename (filename used when reading a dtree from disk)\n\t-h (help, prints this usage statement\n\n");
	Halt();
    }


    return args;
}


int
main(int argc, char **argv) {


    Args args = parseArgs(argc, argv);

    MetricSpace *space = CreateMetricSpace(args.spacetype, args.data);

    DispersionTree *tree;
    if (args.binaryFileToRead == NULL) { // either create a dtree from scratch
	tree= createDispersionTree(space->points, space->dist, space->numpoints, args.numTreeLevels, 0);
    } else { // or read it from disk
	tree= readDispersionTreeFromFile(space->points, space->dist, space->numpoints, args.binaryFileToRead);
    }

    unsigned i, j, k;
    Pointer neighbors[ args.numNeighbors ];
    char buffer[ 1024 ];

    printf("Number of distance computations needed to construct the tree: %lu\n", numDistanceComputations);
    MetricSpace *queries = CreateMetricSpace(args.spacetype, args.queries);

    for (i=0; i < queries->numpoints; ++i) {
	Pointer query = queries->points[ i ];
	queries->toString( query, buffer );
	printf("Nearest neighbors to %s:\n", buffer);
	k = knnSearch(tree, query, args.numNeighbors, neighbors); // KNN search
	for (j=0; j < k; ++j) {
	    queries->toString( neighbors[j], buffer );
	    printf("%u : %s --- %f\n", j, buffer,  tree->dist(query, neighbors[j]));      
	}
	printf("Number of distance computations: %lu\n", numDistanceComputations);

	float knnDistance = tree->dist(query, neighbors[k-1]);
	k = rangeSearch(tree, query, knnDistance-0.01, neighbors); // range search
	for (j=0; j < k; ++j) {
	    queries->toString( neighbors[j], buffer );
	    printf("%u : %s --- %f\n", j, buffer,  tree->dist(query, neighbors[j]));      
	}
	printf("Number of distance computations: %lu\n", numDistanceComputations);

	k = proximitySearch(tree, query, args.numNeighbors, neighbors, knnDistance); // hybrid range search
	for (j=0; j < k; ++j) {
	    queries->toString( neighbors[j], buffer );
	    printf("%u : %s --- %f\n", j, buffer,  tree->dist(query, neighbors[j]));      
	}
	printf("Number of distance computations: %lu\n\n", numDistanceComputations);
	// note that the hybrid range search set to the KNN distance 
	// will do the same number of distance computations as a knn search
	// this is a result of the range-optimality of doing a best-first search ordered by
	// the lowerbound distance (see Incremental Similarity Search in Multimedia Databases (2000) by Hjaltason and Samet )

    }


    if (args.binaryFileToWrite != NULL) {
	writeDispersionTreeToFile(args.binaryFileToWrite, tree);
    }


//    writeDispersionTreeToHeaderFile("ooga.h", "ooga", tree);

    
    destroyDispersionTree(tree);
    fclose(args.data);
    fclose(args.queries);
    DestroyMetricSpace(space);
    DestroyMetricSpace(queries);
    return 0;
}



