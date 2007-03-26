This is a beta of the idcoefs 2.0 software for computing generalized
kinship coefficients and condensed identity coefficients. This
software is copyright Mark Abney 2004, 2005 and licensed under the GPL
2.0. I may be contacted by email at abney@uchicago.edu. Following are
some instructions on compiling and running the software.

COMPILING:

All code is written in ANSI/ISO C and should compile on any platform
with a compliant compiler. The Makefile assumes gcc is installed, but
any ISO C compiler should work. To compile type "make" (without
quotes) from the command line. This will create a binary executable
file called "idcoefs". This file may be moved to anywhere you see fit.

NOTE: THERE IS AN OPTIMIZATION BUG IN GCC VERSIONS 3.0 - 3.3 THAT CAN
RESULT IN THIS PROGRAM GIVING INCORRECT RESULTS. To avoid this use a
version of gcc >= 3.3.1. If you must use a version of gcc affected by
this bug, type "make idcoefs33" at the command line. This will turn off
the buggy portion of the compiler optimization and give a correctly
functioning executable.

RUNNING:

The idcoefs program requires two input files and one output file,
which may be specified using command line arguments. Also, the program
can take an argument specifying the amount of RAM available for the
computations. The command line arguments are as follows:

-p pedfile    The name of the pedigree file. This defaults to
   	      "pedigree" if this flag is not used.

-s studyfile  The name of the file which includes the study sample
   	      individuals for whom identity coefficients will be 
	      computed. Defaults to "study.sample"

-o outfile    The name of the file which will contain the identity
   	      coefficients for all pairs specified in studyfile.
	      Defaults to "output"

-r N	      Allow the program to use approximately N megabytes of 
   	      RAM. For very large pedigrees N should be as large as
	      possible, but not so large as to cause the OS to start
	      using swap space on the hard disk. If N is too small,
	      the computations will take longer. Defaults to 100, but
	      you probably want to set it higher.

NOTE: There must be a space between the flag (e.g. -p) and the
argument for that flag!

FILE FORMATS:

Pedigree file:
Each line in the pedigree file should correspond to one individual in
the pedigree. Each line should be separated into at least three
columns with white space separating the columns. The first column is
the ID of the person. The next two columns are the IDs of the parents
(whether the mother or father comes first doesn't matter). Any
additional columns beyond these first three are ignored. The IDs must
all be integers. The ID 0 is reserved and means "unknown." If a person
is a founder, both parents must have ID 0. All non-founders must have
both parents specified. The pedigree must be ordered in such a way
that parents always appear earlier in the file than children.

Studyfile:
The study file may have either of two formats. The first is that each
line contains a pair of IDs separated by white space. The program will
then compute the condensed identity coefficients for each of these
pairs. The other format is to have one ID per line, in which case
condensed identity coefficients will be computed for every possible
pair of IDs.

Outfile:
The output file contains one line for each pair for which condensed
identity coefficients were computed. The first two columns are the IDs
for the pair and the next nine columns have coefficients Delta_1,...,
Delta_9, where Delta_x is Jacquard's xth condensed identity
coefficient. [Jacquard (1970 - French edition, 1974 - English) The
Genetic Structure of Populations; also in Lynch and Walsh (1998)
Genetics and Analysis of Quantitative Traits].

EXAMPLE: 
In the Example directory is an example pedigree file (ex.pedigree),
study sample file (ex.study) and the output (ex.output) that should be
obtained after running the program. To run the program copy the
idcoefs executable somewhere into your path and type from the command
line: 

idcoefs -p ex.pedigree -s ex.study -o new.output -r 200

Note that 100-200 MB of RAM seems to be about right for this size
pedigree. If you try it with 50 MB, say, the length of time to run
will probably go up. The contents of new.output should be the same as
in ex.output.
