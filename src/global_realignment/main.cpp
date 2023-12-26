/****************************************************************************
 *    This file is part of the program ReformAlign (Reformed Alignments)    *
 *               a profile-based meta-alignment approach                    *
 *                                                                          *
 *    Copyright (C) 2014  Dimitrios Lyras                                   *
 *    e-mail: dimLyras@bio.lmu.de                                           *
 *                                                                          *
 *  This program is free software; you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation; either version 3 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *  GNU General Public License for more details.                            *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.    *
 ****************************************************************************/

#include "ReformAlignInitializer.h"
#include "FSequence.h"

void ShowHelp(); //Function to display the Help options

int main(int argc, char* argv[]) {

    //0. To Make a Custom run with Command Line Arguments //{
    //char* dummy_args[] = {"ReformAlign", "-i", "./Sequences.fasta", "-o", "./ReformAlign_v10.fasta", "-a", "./KAlign.fasta", "-b" };
    //argv = dummy_args;
    //argc = 8;

    //}

    //1. Main Program Variables Definition//{
    string in           = "";  //The input  File name
    string out          = "";  //The output File name
    string param        = "";  //The Parameters File name
    string submatMotif  = "";  //The Substitution matrix File name for the Profile Construction step
    string submatAlign  = "";  //The Substitution matrix File name for the Final Alignment step
    string aligner      = "";  //The aligner name to construct the initial profile from (Prank or KAlign)

    double hgop   = -200.0;
    double hgep   = -20.0;
    double vgop   = -150.0;
    double vgep   = -10.0;
    double factor = 375.0;
    double coeff  = 0.05;

    bool   verbose = false;   //Switch to turn on outputting at stdout or not
    bool   recursive = true;  //Reform Alignments until two successive alignments are not longer different
    int    maxiterations = 5; //The maximum number of iterations to perform

    //Assisting variables to check if files have been defined by the user
    bool AlignerSet = false;
    bool ParamsFileSet = false;
    bool SubMatMotifFileSet = false;
    bool SubMatAlignFileSet = false;

    bool HGopSet = false;
    bool HGepSet = false;
    bool VGopSet = false;
    bool VGepSet = false;
    bool FactorSet = false;
    bool CoeffSet = false;
    //}

    //2. Initialize Program Parameters and Files from arguments //{
    int c;
    opterr = 0;

    while ((c = getopt (argc, argv, "bvhi:o:p:m:l:a:q:w:e:r:f:c:t:")) != -1) {
        switch (c) {
        case 'i':
            in      = string(optarg);
            break;
        case 'o':
            out     = string(optarg);
            break;
        case 'p':
            param   = string(optarg);
            ParamsFileSet = true;
            break;
        case 'm':
            submatMotif = string(optarg);
            SubMatMotifFileSet = true;
            break;
        case 'l':
            submatAlign = string(optarg);
            SubMatAlignFileSet = true;
            break;
        case 'a':
            aligner = string(optarg);
            AlignerSet = true;
            break;
        case 'q':
            hgop = atof(string(optarg).c_str());
            HGopSet = true;
            break;
        case 'w':
            hgep = atof(string(optarg).c_str());
            HGepSet = true;
            break;
        case 'e':
            vgop = atof(string(optarg).c_str());
            VGopSet = true;
            break;
        case 'r':
            vgep = atof(string(optarg).c_str());
            VGepSet = true;
            break;
        case 'f':
            factor = atof(string(optarg).c_str());
            FactorSet = true;
            break;
        case 'c':
            coeff = atof(string(optarg).c_str());
            CoeffSet = true;
            break;
        case 't':
            maxiterations = atoi(string(optarg).c_str());
            break;
        case 'v':
            verbose = true;
            break;
        case 'b':
            recursive = true;
            break;
        case 'h':
        case '?':
            ShowHelp();
            break;
        default:
            abort ();
        }
    }

    //Initialize with default values if arguments were not provided by the user
    if (!in.compare("")) {
        cerr<<"**Error: Input file was not defined.\nPlease type ReformAl -h for usage instructions\n******** ReformAlign will now terminate\n";
        exit(1);
    }

    if (!out.compare("")) {
        int found = in.find_last_of("/\\");
        if (found!=-1) {
            out  = in.substr(0, found) + "/ReformAlign.fasta";
        }else {
            out  = "ReformAlign.fasta";
        }

        cout<<"**Warning:The output file has not been specified. The default output filename: "<<out<<" will be used.\nNote:If the output file exists, it will be automatically overwritten.";
    }

    if (!aligner.compare("")) {
        cerr<<"**Error: Starting Alignmnet file was not defined.\nPlease type ReformAl -h for usage instructions\n******** ReformAlign will now terminate\n";
        exit(1);
    }

    //}

    //3. All Parameters Have been Set. Continue with program execution. //{
    if (ParamsFileSet && (SubMatAlignFileSet || SubMatMotifFileSet)) {
        if (!submatAlign.compare("")) {
            submatAlign = submatMotif;
        }

        if (!submatMotif.compare("")) {
            submatMotif = submatAlign;
        }
        ReformAlignInitializer RInit(in, out, param, submatMotif, submatAlign, aligner, verbose);
        RInit.AlignSequences();
    }else if (HGopSet && HGepSet && VGopSet && VGepSet && FactorSet && CoeffSet) {
        ReformAlignInitializer RInit(in, out, hgop, hgep, vgop, vgep, factor, coeff, aligner, verbose);
        RInit.AlignSequences();
    }else if (AlignerSet) {
        if (recursive) {
            vector<FSequence> SeqsPrev;  string AlphabetPrev;
            vector<FSequence> SeqsAfter; string AlphabetAfter;
            FSequence::ReadAlignmentFromFile(aligner, SeqsPrev, AlphabetPrev);

            int iterCounter = 0;
            while(iterCounter<maxiterations) {
                if (verbose) { cout<<"\n\t\t         ITERATION: "<<(iterCounter+1)<<"\n\t\t         ============"<<endl; }
                ReformAlignInitializer RInit(in, out, aligner, verbose);
                RInit.AlignSequences();
                SeqsAfter.clear();
                SeqsAfter = RInit.Get_Sequences();

                iterCounter++;
                if (FSequence::AreAlignmentsIdentical(SeqsPrev, SeqsAfter) ) {
                    break;
                }else {
                    SeqsPrev = SeqsAfter;
                }
            }
        }else {
            ReformAlignInitializer RInit(in, out, aligner, verbose);
            RInit.AlignSequences();
        }
    }else {
        ShowHelp();
    }
    //}

    return 0;
}

//Function to display the help options
void ShowHelp() {
    string EndDelimiter = "BUG REPORTS";
    string HelpContents = ReadFromFile("ReadMe.txt", true);
    size_t start_loc = HelpContents.find("PROGRAM EXECUTION",0);
    size_t end_loc = HelpContents.find(EndDelimiter,0);

    cout<<HelpContents.substr(start_loc, end_loc-start_loc)<<endl<<endl;
}


