ReformAlign (REFORMed ALIGNments)
=================================

ReformAlign is a profile-based meta-alignment approach that may improve on the 
quality of the deriving alignments from popular and widely used aligners. 
Currently, ReformAlign only supports nucleic (DNA/RNA) sequences that do not 
contain ambiguous characters. Gaps in the alignment should be represented 
using the dash (-) character.

ReformAlign was developed by

Dimitrios Lyras
e-mail: dimLyras@bio.lmu.de 

The homepage of ReformAlign is:
www: http://www.XXXXXXXXX.com
................................................................................

WHAT IS ReformAlign

The main idea behind ReformAlign is quite straightforward: at first, an existing 
alignment is used to construct a non-probabilistic (standard) profile which 
efficiently encapsulates all the alignment information and then all the 
examined sequences are re-aligned against the formed profile. Surprisingly, 
the employment of the ReformAlign logic may often result in alignments which are 
significantly better in terms of accuracy than the starting alignments.

ReformAlign does not come to replace other popular alignment techniques. Instead, 
users may continue to use their aligner(s) of preference and then 
complementarily employ ReformAlign as a final meta-alignment step to examine if 
the "reformed" alignment is more appropriate for their analyses.

Apparently, ReformAlign requires that a starting alignment is provided as
input to the program. Since no aligners are included in the distributed 
packages, you are advised to download and install your aligner of preference
before executing ReformAlign.

Also, please note that currently, ReformAlign only supports nucleic (DNA/RNA) 
sequences that do not contain ambiguous characters. Gaps in the alignment 
should be represented by the dash (-) character.
................................................................................

PROGRAM INSTALLATION AND REQUIREMENTS

ReformAlign is developed using the C++ programming language and also makes use 
of the Open Multi-Processing (OpenMP) parallelization API. It is also important 
to mention that the current implementation follows the C++11 Standard.

Consequently, a C++ compiler, the C++11 library as well as the OpenMP
parallelization library have to be installed to your machine prior to execution.
For a detailed list of compilers or open source communities that implement the 
OpenMP API as well as for compilation instructions please refer to the OpenMP 
website (http://openmp.org/wp/openmp-compilers/). 

The provided Makefile(s) also have to be updated in order to define the location
of the OpenMP library. Simply open the corresponding Makefile (for Windows or 
Unix) and change the LIB tag to point to the location of the:
a) libgomp-1.dll (for Windows)
b) -lgomp (for Unix)

Then, if you have the GNU compiler collection (see  
http://www.fsf.org/software/gcc/gcc.html ) and the "make" command is also 
available in your system, simply run the make command in the directory where 
the files:AlignmentMetrics.h, FNW.h, FSequence.h, Profile.h, 
ReformAlignInitializer.h, TypeDefinitions.h, Utils.h, AlignmentMetrics.cpp, 
FNW.cpp, FSequence.cpp, main.cpp, Profile.cpp, ReformAlignInitializer.cpp are 
located. If everything works, an executable file named ReformAlign will be 
created in the "bin\Release\" directory . 

Note: If you run ReformAlign in Windows, the following dlls:
1. libgcc_s_dw2-1.dll
2. libgomp-1.dll
3. libstdc++-6.dll
4. pthreadGC2.dll
should also be copied over to the folder where the ReformAlign executable is
compiled.
................................................................................

PROGRAM EXECUTION

After the program has been compiled successfully and the executable is created,
ReformAlign may be executed using one of the three supported calls:
   1.ReformAlign -i [InFile] -o [OutFile] -a [InitAlignFile]
   2.ReformAlign -i [InFile] -o [OutFile] -a [InitAlignFile] -p [ParamsFile] 
	  		     -m [SubWeightsFile]
   3.ReformAlign -i [InFile] -o [OutFile] -a [InitAlignFile] -q [HGOP] -w [HGEP] 
			     -e [VGOP] -r [VGEP] -f [FACTOR] -c [COEFFICIENT]

where:  [InFile] is the file containing the Unaligned Sequences in fasta format.
       [OutFile] is the name of the file where the Reformed alignment will be 
				 saved (in fasta format).
 [InitAlignFile] is the name of the file containing the initial alignment 
			    (in fasta format).
    [ParamsFile] is the name of the file containing the (user specified) Gap 
				 opening and Extension penalties.
[SubWeightsFile] is the name of the file containing the (user specified) 
				 Substitution Weights.
          [HGOP] is a (user-specified) value used for the Horizontal Gap Opening 
		         Penalty.
          [HGEP] is (user-specified) value used for the Horizontal Gap 
		         Extension Penalty.
          [VGOP] is (user-specified) value used for the Vertical Gap Opening 
		         Penalty.
          [VGEP] is a (user-specified) value used for the Vertical Gap Extension 
		         Penalty.
        [FACTOR] is a (user-specified) bonus value applied to the default HOXD 
		         Substitution Weights.
   [COEFFICIENT] is a (user-specified) value applied to the augmented HOXD 
		         Substitution Weights in order to further extrapolate the match 
				 and mismatch score differences.

Optionally, the switch -v may be added at the end of the call 
(e.g.ReformAlign -i [InFile] -o [OutFile] -a [InitAlignFile] -v)
to enable the "Verbose" function of the program. If the program is executed in 
Verbose mode then, comments regarding the alignment status will be printed on 
the standard output stream. By default (if the -v is not added), the Verbose 
option is turned off.

In a similar way the -t [ITERATIONS] switch may be optionally used to define the
(maximum) number of iterations that will be applied before ReformAlign outputs 
the reformed alignment. If this option is not specified, then the number 5 
will be employed by default as maximum number of iterations.
................................................................................

ReformAlign CALLS AND FILES DESCRIPTION

1.ReformAlign -i [InFile] -o [OutFile] -a [InitAlignFile]
---------------------------------------------------------
This is the weights auto tuning mode. In this mode, the gap penalties and the 
substitution weights are automatically decided based on the average pairwise 
sequence identity of the sequences in the starting alignment ([InitAlignFile]).
A detailed explanation of the required arguments and the data format is provided
below.

A. [InFile]
	The input file containing the Unaligned sequences should be formatted 
	according to the FASTA format. Each sequence should begin with a single
	line description starting with the ">" symbol followed by the sequence
	identifier. Then the sequence data should follow in the subsequent line(s).
	An example of three sequences in the desired input format is provided below:
	
	>Seq_A
	ACGGATCATTG
	>Seq_B
	GGGCTGAT
	>Seq_C
	GGAAAT
	
	The sequence data should contain nothing else but the letters A,C,G, and T 
	for DNA Sequences or the letters A,C,G, and U for RNA sequences. All 
	letters should be upper-cased.

B. [OutFile]
	The output file is the file where the the reformed alignment will be saved.
	As in the case of the input file, the sequences will be printed in the file
	according to the FASTA format. In this case, the aligned sequence data may
	be comprised of the letters A,C,G, and T for DNA Sequences or the letters 
	A,C,G, and U for RNA sequences, and also the dash character "-" may be used
	to represent gaps in the alignment. 
	
C. [InitAlignFile]
	The initial alignment should also be in FASTA format. Gaps in the alignment 
	should be represented using the dash "-" character. An example of the 
	contents of an initial alignment file is provided below:
	
	>Seq_A
	ACGGATCATTG
	>Seq_B
	--GGGCTGAT-
	>Seq_C
	----GGAAAT-
	

2.ReformAlign -i [InFile] -o [OutFile] -a [InitAlignFile] -p [ParamsFile] 
		   -m [SubWeightsFile]
-------------------------------------------------------------------------
In this mode the gap penalties and the substitution weights are explicitly
specified by the user and are received in the program from the corresponding
input files. For a detailed description of the [InFile], [OutFile] and
[InitAlignFile] arguments, please refer to the first function call above. 

A. [ParamsFile]
	This is the file where the Gap penalties should be specified. The format of
	the parameter file should be as follows:
	
	HGapOpeningPenalty=[Value]
	HGapExtensionPenalty=[Value]
	VGapOpeningPenalty=[Value]
	VGapExtensionPenalty=[Value]
	TerminalGapOpeningPenalty=0.0
	TerminalGapExtensionPenalty=0.0
	
	Please note that there must be not spaces between the gap tag, the equality
	operator and the input value.
	
B. [SubWeightsFile]
	This is the file where the Substitution weights should be specified. The 
	format of the Substitution weights file should be as follows:
	
	[A-A] [A-C] [A-G] [A-T] 0.0
	[C-A] [C-C] [C-G] [C-T] 0.0
	[G-A] [G-C] [G-G] [G-T] 0.0
	[T-A] [T-C] [T-G] [T-T] 0.0
	0.0   0.0   0.0   0.0   0.0
	
	where [X-Y] is the substitution value for matching letter X to Y.
	
	Please note that each value should be separated by the previous one by a 
	single tab.
	
3. ReformAlign -i [InFile] -o [OutFile] -a [InitAlignFile] -q [HGOP] -w [HGEP] 
	  	    -e [VGOP] -r [VGEP] -f [FACTOR] -c [COEFFICIENT]
------------------------------------------------------------------------------
In this mode the user may directly specify the desired values for the:
	   [HGOP]: Horizontal Gap Opening Penalty
	   [HGEP]: Horizontal Gap Extension Penalty
	   [VGOP]: Vertical Gap Opening Penalty
	   [VGEP]: Vertical gap Extension Penalty
	 [FACTOR]: the (positive bonus) value by which the default HOXD Substitution 
		       values will be increased
[COEFFICIENT]: the coefficient by which the augmented HOXD Substitution 
		       values will be increased in order to extrapolate the match - 
			   mismatch score differences
			   
	Please not that in this mode the HOXD Substitution Matrix (Chiaromonte F, 
	Yap VB, Miller W: Scoring pairwise genomic sequence alignments. Pac. Symp. 
	Biocomput. Pac. Symp. Biocomput. 2002:115–126) will be employed by default
	and all the matrix values will be increased by the provided [FACTOR] value 
	and then each value will be proportionally further increased by the 
	specified [COEFFICIENT] value.
	To specify a substitution weights matrix other than the HOXD, please use the
	second option (2) to execute ReformAlign.	
	
	For a detailed description of the [InFile], [OutFile] and [InitAlignFile] 
	arguments, please refer to the first function call above.

................................................................................

LICENSE

The software ReformAlign is distributed freely under the following license terms:

    Copyright (C) 2014  Dimitrios Lyras

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Please refer to the GPL-LICENSE.txt file for further licensing details.

If you use this software for scientific work, please cite the following 
article:

"ReformAlign: Improved multiple sequence alignments using a profile-based 
meta-alignment approach ",(Dimitrios P. Lyras, Dirk Metzler, BMC Bioinformatics,
2014, under review)
................................................................................

BUG REPORTS

For bugs or recommendations, please send an e-mail to:

Dimitrios Lyras
Email: dimLyras@bio.lmu.de


