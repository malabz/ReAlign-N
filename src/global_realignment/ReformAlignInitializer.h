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

#ifndef ReformAlignINITIALIZER_H
#define ReformAlignINITIALIZER_H

#include "Utils.h"
#include "TypeDefinitions.h"
#include "FSequence.h"
#include "Profile.h"
#include "FNW.h"
#include "AlignmentMetrics.h"

using namespace std;
using namespace Utils;
using namespace TypeDefinitions;

class ReformAlignInitializer {

    typedef vector<FSequence>  Sequences;
    typedef vector<FSequence>& SequencesRef;

    public:
        ReformAlignInitializer(string in, string out, string param, string submatMotif, string submatAlign, string aligner, bool Verbose);
        ReformAlignInitializer(string in, string out, double hgop, double hgep, double vgop, double vgep, double factor, double coeff, string aligner, bool Verbose);
        ReformAlignInitializer(string in, string out, string aligner, bool Verbose);
        virtual ~ReformAlignInitializer();

        //The Copy Constructor
        ReformAlignInitializer(const ReformAlignInitializer& other) {
            InputFile        = other.InputFile;
            OutputFile       = other.OutputFile;
            ParamsFile       = other.ParamsFile;
            SubMatsMotifFile = other.SubMatsMotifFile;
            SubMatsAlignFile = other.SubMatsAlignFile;
            BaseAligner      = other.BaseAligner;
            Alphabet         = other.Alphabet;

            HGapOpPenalty    = other.HGapOpPenalty;
            HGapExPenalty    = other.HGapExPenalty;
            VGapOpPenalty    = other.VGapOpPenalty;
            VGapExPenalty    = other.VGapExPenalty;
            TermGapOpPenalty = other.TermGapOpPenalty;
            TermGapExPenalty = other.TermGapExPenalty;
            BonusFactor      = other.BonusFactor;
            Coefficient      = other.Coefficient;

            SubMatMotif.resize(other.SubMatMotif.size(), DoubleVec(other.SubMatMotif[0].size(), 0.0));
            SubMatAlign.resize(other.SubMatAlign.size(), DoubleVec(other.SubMatAlign[0].size(), 0.0));
            for (unsigned int i=0;i<other.SubMatMotif.size();i++) {
                for (unsigned int j=0;j<other.SubMatMotif[i].size();j++) {
                    SubMatMotif[i][j] = other.SubMatMotif[i][j];
                    SubMatAlign[i][j] = other.SubMatAlign[i][j];
                }
            }

            for (unsigned int i=0;i<other.Seqs.size();i++) {
                FSequence seq(other.Seqs[i]);
                Seqs.push_back(seq);
            }

            for (unsigned int i=0;i<other.Aligners.size();i++) {
                FNW Align(other.Aligners[i]);
                Aligners.push_back(Align);
            }
        }

        //The Assignment Operator
        ReformAlignInitializer & operator= (const ReformAlignInitializer& other) {
            if (this != &other) { // protect against invalid self-assignment
                InputFile        = other.InputFile;
                OutputFile       = other.OutputFile;
                ParamsFile       = other.ParamsFile;
                SubMatsMotifFile = other.SubMatsMotifFile;
                SubMatsAlignFile = other.SubMatsAlignFile;
                BaseAligner      = other.BaseAligner;
                Alphabet         = other.Alphabet;

                HGapOpPenalty    = other.HGapOpPenalty;
                HGapExPenalty    = other.HGapExPenalty;
                VGapOpPenalty    = other.VGapOpPenalty;
                VGapExPenalty    = other.VGapExPenalty;
                TermGapOpPenalty = other.TermGapOpPenalty;
                TermGapExPenalty = other.TermGapExPenalty;
                BonusFactor      = other.BonusFactor;
                Coefficient      = other.Coefficient;

                SubMatMotif.resize(other.SubMatMotif.size(), DoubleVec(other.SubMatMotif[0].size(), 0.0));
                SubMatAlign.resize(other.SubMatAlign.size(), DoubleVec(other.SubMatAlign[0].size(), 0.0));
                for (unsigned int i=0;i<other.SubMatMotif.size();i++) {
                    for (unsigned int j=0;j<other.SubMatMotif[i].size();j++) {
                        SubMatMotif[i][j] = other.SubMatMotif[i][j];
                        SubMatAlign[i][j] = other.SubMatAlign[i][j];
                    }
                }

                for (unsigned int i=0;i<other.Seqs.size();i++) {
                    FSequence seq(other.Seqs[i]);
                    Seqs.push_back(seq);
                }

                for (unsigned int i=0;i<other.Aligners.size();i++) {
                    FNW Align(other.Aligners[i]);
                    Aligners.push_back(Align);
                }
            }
            return *this;
        }

        ///Setter and Getter Methods
        void Set_InputFile        (string val)  { InputFile         = val;  }
        void Set_OutputFile       (string val)  { OutputFile        = val;  }
        void Set_ParamsFile       (string val)  { ParamsFile        = val;  }
        void Set_SubMatsMotifFile (string val)  { SubMatsMotifFile  = val;  }
        void Set_SubMatsAlignFile (string val)  { SubMatsAlignFile  = val;  }
        void Set_BaseAligner      (string val)  { BaseAligner       = val;  }
        void Set_Alphabet         (string val)  { Alphabet          = val;  }

        void Set_HGapOpPenalty    (double val)  { HGapOpPenalty     = val;  }
        void Set_HGapExPenalty    (double val)  { HGapExPenalty     = val;  }
        void Set_VGapOpPenalty    (double val)  { VGapOpPenalty     = val;  }
        void Set_VGapExPenalty    (double val)  { VGapExPenalty     = val;  }
        void Set_TermGapOpPenalty (double val)  { TermGapOpPenalty  = val;  }
        void Set_TermGapExPenalty (double val)  { TermGapExPenalty  = val;  }
        void Set_BonusFactor      (double val)  { BonusFactor       = val;  }
        void Set_Coefficient      (double val)  { Coefficient       = val;  }
        void Set_IsVerbose        (bool   val)  { IsVerbose         = val;  }
        void Set_SubMatMotif(DoubleMatRef val)  {
            SubMatMotif.resize(val.size(), vector<double>(val[0].size(), 0.0));
            for (unsigned int i=0;i<val.size();i++) {
                for (unsigned int j=0;j<val[i].size();j++) {
                    SubMatMotif[i][j] = val[i][j];
                }
            }
        }
        void Set_SubMatAlign(DoubleMatRef val)  {
            SubMatAlign.resize(val.size(), vector<double>(val[0].size(), 0.0));
            for (unsigned int i=0;i<val.size();i++) {
                for (unsigned int j=0;j<val[i].size();j++) {
                    SubMatAlign[i][j] = val[i][j];
                }
            }
        }
        void Set_Sequences  (SequencesRef val)  {
            for (unsigned int i=0;i<val.size();i++) {
                FSequence seq(val[i]);
                Seqs.push_back(seq);
            }
        }


        string          Get_InputFile()               { return InputFile;         }
        string          Get_OutputFile()              { return OutputFile;        }
        string          Get_ParamsFile()              { return ParamsFile;        }
        string          Get_SubMatsMotifFile()        { return SubMatsMotifFile;  }
        string          Get_SubMatsAlignFile()        { return SubMatsAlignFile;  }
        string          Get_BaseAligner()             { return BaseAligner;       }
        string          Get_Alphabet()                { return Alphabet;          }

        double          Get_HGapOpPenalty()           { return HGapOpPenalty;     }
        double          Get_HGapExPenalty()           { return HGapExPenalty;     }
        double          Get_VGapOpPenalty()           { return VGapOpPenalty;     }
        double          Get_VGapExPenalty()           { return VGapExPenalty;     }
        double          Get_TermGapOpPenalty()        { return TermGapOpPenalty;  }
        double          Get_TermGapExPenalty()        { return TermGapExPenalty;  }
        double          Get_BonusFactor()             { return BonusFactor;       }
        double          Get_Coefficient()             { return Coefficient;       }
        DoubleMatRef    Get_SubMatMotif()             { return SubMatMotif;       }
        DoubleMatRef    Get_SubMatAlign()             { return SubMatAlign;       }
        SequencesRef    Get_Sequences()               { return Seqs;              }
        bool            Get_IsVerbose()               { return IsVerbose;         }

        void      MakeInitializations();
        void      MakeInitializations(double hgop, double hgep, double vgop, double vgep, double factor);
        void      MakeInitializations(bool autotune);
        void      AlignSequences();

        void AlignAgainstProfile(Profile& BasicProfile, vector<Profile>& SeqProfiles, int& steps); //Function to align (in Parallel) each sequence against the Profile (and perform Profile Fine-tuning if required)
        void CleanAlignment(); //Function to remove all Gapped Columns from the Alignment

        //Output Functions
        string    GetParamsToStr();
        string    GetSequencesToStr();
        string    GetSettingsToStr();
        string    GetSubMatToStr(DoubleMatRef mat);

    protected:

    private:
        /** Class Members **/
        string    InputFile;        //The input filename
        string    OutputFile;       //The output filename
        string    ParamsFile;       //The parameters filename (for the gap penalties)
        string    SubMatsMotifFile; //The substitution weights filename for the Profile Creation phase
        string    SubMatsAlignFile; //The substitution weights filename for the Sequence Alignment phase
        string    BaseAligner;      //The starting alignment file to create the starting profile from
        string    Alphabet;         //The employed alphabet for the residues(i.e. {A,C,G,T} for DNA, {A,C,G,U} for RNA)

        double    HGapOpPenalty;    //The Horizontal Gap Opening penalty
        double    HGapExPenalty;    //The Horizontal Gap Extension penalty
        double    VGapOpPenalty;    //The Vertical Gap Opening penalty
        double    VGapExPenalty;    //The Vertical Gap Extension penalty
        double    TermGapOpPenalty; //The Terminal Gap Opening penalty (set to 0.0 by default to not penalize termanal gaps)
        double    TermGapExPenalty; //The Terminal Gap Extension penalty (set to 0.0 by default to not penalize termanal gaps)
        double    BonusFactor;      //The (positive) bonus value by which all value in the HOXD Substitution Matrix will be increased
        double    Coefficient;      //The Coefficient by which the Substitution values are multiplied with to extrapolate differences between matches-mismatches
        DoubleMat SubMatMotif;      //The substitution weights matrix for the Profile Creation phase
        DoubleMat SubMatAlign;      //The substitution weights matrix for the Sequence Alignment phase
        Sequences Seqs;             //The input sequences

        vector<FNW> Aligners;       //A vector of aligners - each one having the appropriate adjusted weights according to the heuristic gap updates

        bool IsVerbose;               //If in verbose mode, the aligner status will be printed on the stdout
};

#endif // ReformAlignINITIALIZER_H
