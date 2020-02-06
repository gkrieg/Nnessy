/* **********************************************************
 * Author: Benjamin B. Salazar
 * Date: June 30th, 2010
 *
 * For a given set of protein sequence/secondary structure
 * sequence pairs, this program generates three files:
 * each one consisting of amino acid words from either
 * alpha-helix, beta-strand, or coil secondary
 * structure types.
 * **********************************************************/

#include "filegen.h"

/* *************************************************************************************************************************** */
int main(int argc, char *argv[]) {
/* *************************************************************************************************************************** */
	ClassFiles files = {NULL, NULL, NULL, 0, 0, 0};
	OptionFlags flags = {0, 0, 0, 0, 0, 0};
	char alphaFileName[50];
	char betaFileName[50];
	char otherFileName[50];
	int optind;
	int i;

	/* ******************************* */
	/* handle the command-line options */
	/* ******************************* */
	optind = processOpts(argc, argv, &flags);
	progName = argv[0];
	argc -= optind;
	argv += optind;

	if (flags.unknownFlag || !flags.lengthFlag || flags.length < 3 || flags.length % 2 == 0 || argc < 1) usage(progName, stderr);
	if (flags.usageFlag) usage(progName, stdout);

	sprintf(alphaFileName, "%d", flags.length);
	sprintf(betaFileName, "%d", flags.length);
	sprintf(otherFileName, "%d", flags.length);

	debug = flags.debugFlag;

	/* ************************************ */
	/* initialize the class files to create */
	/* ************************************ */
	if (flags.pureFlag) {
		strcat(alphaFileName, "alphaPure.words");
		strcat(betaFileName, "betaPure.words");
		strcat(otherFileName, "otherPure.words");
	} else {
		strcat(alphaFileName, "alpha.words");
		strcat(betaFileName, "beta.words");
		strcat(otherFileName, "other.words");
	}

	files.alpha = fopen(alphaFileName, "wt");
	files.beta = fopen(betaFileName, "wt");
	files.other = fopen(otherFileName, "wt");
	files.numAlpha = 0;
	files.numBeta = 0;
	files.numOther = 0;

	/* ************************************************************* */
	/* iterate over each file that is specified via the command line */
	/* ************************************************************* */
	for (i = 0; i < argc; i++) {
		FILE *sequenceFile;

		if ((sequenceFile = fopen(argv[i], "rt")) == NULL) {
			fprintf(stderr, "%s: could not open '%s' for reading\n", progName, argv[i]);
			continue;
		}

		/* ****************************************************************************** */
		/* create amino acid words of length lengthOfWord for the specified sequence file */
		/* ****************************************************************************** */
		if (readSequenceFile(sequenceFile, flags.length, &files, flags.pureFlag)) {
			fprintf(stderr, "%s: invalid sequence file '%s'\n", progName, argv[i]);
			continue;
		}

		fclose(sequenceFile);
	}

	/* possibly, add number of words for each file? */

	/* *********************** */
	/* close the created files */
	/* *********************** */
	fclose(files.alpha);
	fclose(files.beta);
	fclose(files.other);
	return EXIT_SUCCESS;
}

/* *************************************************************************************************************************** */
int readSequenceFile(FILE* sequenceFile, int lengthOfWord, ClassFiles *files, int pure) {
/* *************************************************************************************************************************** */
	char proteinSequence[MAX_LENGTH_OF_LINE];
	char secondarySequence[MAX_LENGTH_OF_LINE];
	int flag;

	/* **************************** */
	/* zero-out the two char arrays */
	/* **************************** */
	int i;
	for (i = 0; i < MAX_LENGTH_OF_LINE; i++) {
		proteinSequence[i] = '\0';
		secondarySequence[i] = '\0';
	}

	/* ******************************************************************** */
	/* process each couple of sequences in the file (protein and secondary) */
	/* ******************************************************************** */
	flag = 0;
	while (!feof(sequenceFile)) {
		/* ********************** */
		/* read the protein label */
		/* ********************** */
		do {
			if (fgets(proteinSequence, sizeof(proteinSequence), sequenceFile) == NULL) {
				if (!feof(sequenceFile))
					return EXIT_FAILURE;
				else {
					flag = 1;
					break;
				}
			}
		} while (proteinSequence[0] != '>');

		if (flag == 1)
			break;

		/* ***************************** */
		/* retrieve the protein sequence */
		/* ***************************** */
		if (fgets(proteinSequence, sizeof(proteinSequence), sequenceFile) == NULL)
			return EXIT_FAILURE;

		/* ************************ */
		/* read the secondary label */
		/* ************************ */
		if (fgets(secondarySequence, sizeof(proteinSequence), sequenceFile) == NULL)
			return EXIT_FAILURE;

		secondarySequence[20] = '\0';

		if (strcmp(secondarySequence, ">Secondary Structure"))
			return EXIT_FAILURE;

		/* ******************************* */
		/* retrieve the secondary sequence */
		/* ******************************* */
		if (fgets(secondarySequence, sizeof(secondarySequence), sequenceFile) == NULL)
			return EXIT_FAILURE;

		/* ******************************************************** */
		/* process the protein and corresponding secondary sequence */
		/* ******************************************************** */
		if (processSequences(proteinSequence, secondarySequence, lengthOfWord, files, pure))
			return EXIT_FAILURE;

		/* ******************* */
		/* read the blank line */
		/* ******************* */
		fgets(proteinSequence, sizeof(proteinSequence), sequenceFile);
	}

	return EXIT_SUCCESS;
}

/* *************************************************************************************************************************** */
int processSequences(char *proteinSequence, char *secondarySequence, int lengthOfWord, ClassFiles *files, int pure) {
/* *************************************************************************************************************************** */
	int lengthOfSequence;
	int i;

	/* ***************************************** */
	/* reduce alphabet in the secondary sequence */
	/* ***************************************** */
	lengthOfSequence = 0;
	while (secondarySequence[lengthOfSequence] != '\n') {
		/* checks that the two sequences have the same length */
		if (proteinSequence[lengthOfSequence] == '\n')
			return EXIT_FAILURE;

		/* ***************************************************************** */
		/* institue mapping used by PSIPRED to reduce 8-letter alphabet to 3 */
		/* ***************************************************************** */
		switch (secondarySequence[lengthOfSequence]) {
		case 'H':
		case 'h':
		case 'g':
		case 'G':	secondarySequence[lengthOfSequence] = 'A';
				break;

		case 'e':
		case 'b':
		case 'E':
		case 'B':	secondarySequence[lengthOfSequence] = 'B';
				break;

		default:
				secondarySequence[lengthOfSequence] = 'C';
				break;
		}

		lengthOfSequence++;
	}

	/* checks that the two sequences have the same length */
	if (proteinSequence[lengthOfSequence] != '\n')
		return EXIT_FAILURE;

	/* ******************************************************************************* */
	/* alpha runs must be >= 5, beta runs must be >= 3 (within the secondary sequence) */
	/* ******************************************************************************* */
	int beginning;
	char currentType;
	int numberOfCharacters;

	beginning = 0;
	currentType = secondarySequence[0];
	numberOfCharacters = 1;
	for (i = 1; i < lengthOfSequence; i++) {
		/* run condition */
		if (secondarySequence[i] == currentType) {
			numberOfCharacters++;

			/* special case at end of sequence */
			if (i == lengthOfSequence - 1)
				if ((numberOfCharacters < 5 && currentType == 'A') || (numberOfCharacters < 3 && currentType == 'B'))
					for (; beginning <= i; beginning++)
						secondarySequence[beginning] = 'C';
		}
		/* change condition */
		else {
			if ((numberOfCharacters < 5 && currentType == 'A') || (numberOfCharacters < 3 && currentType == 'B'))
				for (; beginning < i; beginning++)
					secondarySequence[beginning] = 'C';

			/* special case at end of sequence */
			if (i == lengthOfSequence - 1) {
				switch (secondarySequence[i]) {
				case 'A':
				case 'B':	secondarySequence[i] = 'C'; break;
				default:	break;
				}
			} else {
				beginning = i;
				currentType = secondarySequence[i];
				numberOfCharacters = 1;
			}
		}
	}

	/* ********************************************************************************* */
	/* gather the protein words of length lengthOfWord into their respective class files */
	/* ********************************************************************************* */
	int currentResidue;
	int middleResidue;
	char word[lengthOfWord + 1];
	char firstSecondaryResidue;
	char currentSecondaryResidue;
	int flag;

	middleResidue = lengthOfWord / 2;

	for (currentResidue = 0; currentResidue <= (lengthOfSequence - lengthOfWord); currentResidue++) {
		flag = 0;		
		/* **************************************** */
		/* retrieve the word at the current residue */
		/* **************************************** */
		strncpy(word, (proteinSequence + currentResidue), lengthOfWord);
		word[lengthOfWord] = '\0';

		/* *************************************************** */
		/* examine the amino acid word and secondary structure */
		/* *************************************************** */
		firstSecondaryResidue = secondarySequence[currentResidue];
		for (i = 0; i < lengthOfWord; i++) {
			/* treat lowercase letters as uppercase */
			if (word[i] >= 'a' && word[i] <= 'z') {
				if (DEBUG && debug)
					printf("lower case letter %c found in word %s\n", word[i], word);

				word[i] = (char) ((int) word[i] - 32);
			}

			/* check for invalid residue */
			if (residueToInt[(int) word[i]] < 0) {
				if (DEBUG && debug)
					printf("invalid residue %c in word %s\n", word[i], word);

				if (word[i] >= 'A' && word[i] <= 'Z')
					word[i] = 'X';
				else {
					flag = 1;
					break;
				}
			}

			/* check for purity of secondary structure (each residue in amino acid word belongs to the same class) */
			if (pure && i) {
				currentSecondaryResidue = secondarySequence[currentResidue + i];

				if (firstSecondaryResidue != currentSecondaryResidue) {
					flag = 1;
					break;
				}
			}
		}

		if (flag)
			continue;

		/* ************************************************************************ */
		/* check if the word belongs in the alpha, beta, or other class by checking */
		/* the middle residue of the secondary sequence                             */
		/* ************************************************************************ */
		switch (secondarySequence[currentResidue + middleResidue]) {
		case 'A':	if (files->numAlpha != 0)
					fprintf(files->alpha, "\n");
				fprintf(files->alpha, "%s", word);
				files->numAlpha++;
				break;

		case 'B':	if (files->numBeta != 0)
					fprintf(files->beta, "\n");
				fprintf(files->beta, "%s", word);
				files->numBeta++;
				break;

		case 'C':	if (files->numOther != 0)
					fprintf(files->other, "\n");
				fprintf(files->other, "%s", word);
				files->numOther++;
				break;
		default:
				break;
		}
	}

	return EXIT_SUCCESS;
}

/* *********************************************************************************************************** */
void usage(char *progName, FILE *file) {
/* *********************************************************************************************************** */
	/* print the usage statement for mgenerator */
	fprintf(file, "\nUsage: %s [-dp] -l length file [file ...]\n", progName);
	fprintf(file, "  --usage           prints the usage information and exits\n");
	fprintf(file, "  -d --debug        activates debug messages\n");
	fprintf(file, "  -p --pure         create pure amino acid word files\n\n");
	fprintf(file, "length is the length of amino acid words that will be output (must be greater than 0).\n\n");
	fprintf(file, "file is a file containing training data that consists of protein and secondary structure sequences.\n");
	fprintf(file, "At least one file must be referenced, though an indefinite amount of files may be specified.\n\n");
	fprintf(file, "NOTES:\n");
	fprintf(file, "While the file(s) containing training data may consist of any number of comments not beginning\n");
	fprintf(file, "with the character '>', protein and secondary sequence couples must adhere to the structure demonstrated in\n");
	fprintf(file, "the following example:\n\n");
	fprintf(file, ">12AS\n");
	fprintf(file, "AYIAKQRQISFVKSH[...]\n");
	fprintf(file, ">Secondary Structure\n");
	fprintf(file, "LHHHHHHHHHHHHHH[...]\n\n");
	fprintf(file, "In order: >[protein name], [protein sequence], >Secondary Structure, [secondary structure sequence],\n");
	fprintf(file, "each of which exist on their own line.  Note that in the previous example, the four lines \n");
	fprintf(file, "will always be followed by one or more blank lines.\n\n");

	if (file == stdout)
		exit(EXIT_SUCCESS);
	else
		exit(EXIT_FAILURE);
}

int processOpts(int argc, char *argv[], OptionFlags *flags) {
	/* loop over every option from the command-line	and set the appropriate flags */
	char ch;
	while ((ch = getopt_long(argc, argv, "dpl:", longopts, NULL)) != -1) {
		switch (ch) {
		case 0: break;
		case 'd': flags->debugFlag = 1; break;
		case 'p': flags->pureFlag = 1; break;
		case 'l':
			if (*optarg == '=')
				optarg++;

			flags->lengthFlag = 1;
			flags->length = atoi(optarg);
			break;
		case 'u': flags->usageFlag = 1; break;
		case '?': flags->unknownFlag = 1; break;
		default: break;
		}
	}

	/* return the number of options that were processed */
	return optind;
}

