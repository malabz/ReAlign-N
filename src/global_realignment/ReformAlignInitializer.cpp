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

/** Contructors and Decontructors **/ //{

ReformAlignInitializer::ReformAlignInitializer(string in, string out, string param, string submatMotif, string submatAlign, string aligner, bool Verbose) {

    //1. Init Members //{
    InputFile           = in;
    OutputFile          = out;
    ParamsFile          = param;
    SubMatsMotifFile    = submatMotif;
    SubMatsAlignFile    = submatAlign;
    BaseAligner         = aligner;
    Alphabet            = DNA_Order;

    BonusFactor         = 268.738034793977;
    Coefficient         = 0.2182;
    HGapOpPenalty       = -394.022939886315;
    HGapExPenalty       = -191.612476671029;
    VGapOpPenalty       = -276.993627921082;
    VGapExPenalty       = -126.071721004077;
    TermGapOpPenalty    = 0.0;
    TermGapExPenalty    = 0.0;

    IsVerbose           = Verbose;
    //}

    //2. Init Parameters From File //{
    ifstream inputstr(ParamsFile.c_str());
    if (!inputstr) {
        cerr << "Error: Cannot open Parameters file: " << ParamsFile.c_str() << endl;
        exit(1);
    }
    string hgo, hge, vgo, vge, tgo, tge;
    inputstr >> hgo >> hge >> vgo >> vge >> tgo >> tge;
    HGapOpPenalty = atof( hgo.substr(hgo.find("HGapOpeningPenalty=")+19).c_str());
    HGapExPenalty = atof(hge.substr(hge.find("HGapExtensionPenalty=")+21).c_str());
    VGapOpPenalty = atof( vgo.substr(vgo.find("VGapOpeningPenalty=")+19).c_str());
    VGapExPenalty = atof(vge.substr(vge.find("VGapExtensionPenalty=")+21).c_str());
    TermGapOpPenalty = atof( tgo.substr(tgo.find("TerminalGapOpeningPenalty=")+26).c_str());
    TermGapExPenalty = atof( tge.substr(tge.find("TerminalGapExtensionPenalty=")+28).c_str());

    inputstr.close();
    //}

    //3. Init Substitution Weights for the Profile Construction Phase from File//{
    SubMatMotif.clear();
    ifstream inputstrSb(SubMatsMotifFile.c_str());
    if (!inputstrSb) {
        cerr << "Error: Cannot open Substitution Matrix file: " << SubMatsMotifFile.c_str() << endl;
        exit(1);
    }
    int BaseCount;
    SubMatMotif.resize(5, vector<double>(5, 0.0));
    BaseCount = 5;

    int i=0;

    string line;
    while(std::getline(inputstrSb, line) && i<BaseCount) {
        SplitValues(line, '\t', SubMatMotif[i]);
        i++;
    }

    inputstrSb.close();
    //}

    //4. Init Substitution Weights for the Profile Construction Phase from File//{
    if (SubMatsAlignFile.compare(SubMatsMotifFile)!=0) {
        SubMatAlign.clear();
        ifstream inputstrSbAl(SubMatsAlignFile.c_str());
        if (!inputstrSbAl) {
            cerr << "Error: Cannot open Substitution Matrix file: " << SubMatsAlignFile.c_str() << endl;
            exit(1);
        }
        int BaseCountAl;
        SubMatAlign.resize(5, vector<double>(5, 0.0));
        BaseCountAl = 5;

        int iAl=0;

        string lineAl;
        while(std::getline(inputstrSb, lineAl) && iAl<BaseCountAl) {
            SplitValues(lineAl, '\t', SubMatAlign[iAl]);
            iAl++;
        }

        inputstrSbAl.close();
    }else {
        SubMatAlign.resize(SubMatMotif.size(), DoubleVec(SubMatMotif[0].size(), 0.0f));
        for (unsigned int i=0;i<SubMatMotif.size();i++) {
            for (unsigned int j=0;j<SubMatMotif[i].size();j++) {
                SubMatAlign[i][j] =  SubMatMotif[i][j];
            }
        }

    }
    //}

    //5. Init the Sequences //{
    Seqs.clear();
    FSequence::ReadSequencesFromFile(InputFile, Seqs, Alphabet);
    //}

    //6. Initialize the Aligners (one Aligner for each sequence - due to the different weights and penalties) //{
    for (unsigned int i=0;i<Seqs.size();i++) {
        FNW Aligner;
        double factor = (1/Seqs[i].Get_Ratio());

        Aligner.SetParameters(factor*HGapOpPenalty, factor*HGapExPenalty, factor*VGapOpPenalty, factor*VGapExPenalty, TermGapOpPenalty, TermGapExPenalty, SubMatMotif, SubMatAlign);
        Aligners.push_back(Aligner);
    }
    //}

}

ReformAlignInitializer::ReformAlignInitializer(string in, string out, double hgop, double hgep, double vgop, double vgep, double factor, double coeff, string aligner, bool Verbose) {

    //1. Init Members //{
    InputFile           = in;
    OutputFile          = out;
    ParamsFile          = "";
    SubMatsMotifFile    = "";
    SubMatsAlignFile    = "";
    BaseAligner         = aligner;
    Alphabet            = DNA_Order;

    HGapOpPenalty       = hgop;
    HGapExPenalty       = hgep;
    VGapOpPenalty       = vgop;
    VGapExPenalty       = vgep;
    TermGapOpPenalty    = 0.0;
    TermGapExPenalty    = 0.0;
    BonusFactor         = factor;
    Coefficient         = coeff;

    IsVerbose           = Verbose;
    //}

    //2. Init Substitution Weights for the Profile COnstruction Phase from Given Settings//{
    SubMatMotif.clear();
    SubMatMotif.resize(5, vector<double>(5, 0.0));
    SubMatMotif[0][0] =   91.0 + factor;
    SubMatMotif[0][1] = -114.0 + factor;
    SubMatMotif[0][2] =  -31.0 + factor;
    SubMatMotif[0][3] = -123.0 + factor;
    SubMatMotif[0][4] = 0.0;

    SubMatMotif[1][0] = -114.0 + factor;
    SubMatMotif[1][1] =  100.0 + factor;
    SubMatMotif[1][2] = -125.0 + factor;
    SubMatMotif[1][3] =  -31.0 + factor;
    SubMatMotif[1][4] = 0.0;

    SubMatMotif[2][0] =  -31.0 + factor;
    SubMatMotif[2][1] = -125.0 + factor;
    SubMatMotif[2][2] =  100.0 + factor;
    SubMatMotif[2][3] = -114.0 + factor;
    SubMatMotif[2][4] = 0.0;

    SubMatMotif[3][0] = -123.0 + factor;
    SubMatMotif[3][1] =  -31.0 + factor;
    SubMatMotif[3][2] = -114.0 + factor;
    SubMatMotif[3][3] =   91.0 + factor;
    SubMatMotif[3][4] = 0.0;

    SubMatMotif[4][0] = 0;
    SubMatMotif[4][1] = 0;
    SubMatMotif[4][2] = 0;
    SubMatMotif[4][3] = 0;
    SubMatMotif[4][4] = 0;

    //Extrapolate MATCH-MISMATCH Differences in the Substitution Weights
    for (unsigned int i=0;i<SubMatMotif.size();i++) {
        for (unsigned int j=0;j<SubMatMotif[i].size();j++) {
            SubMatMotif[i][j] += ( SubMatMotif[i][j] * Coefficient );
        }
    }
    //}

    //3. Init Substitution Weights for the Profile COnstruction Phase from File//{
    SubMatAlign.clear();
    SubMatAlign.resize(5, vector<double>(5, 0.0));
    SubMatAlign[0][0] =   91.0 + factor;
    SubMatAlign[0][1] = -114.0 + factor;
    SubMatAlign[0][2] =  -31.0 + factor;
    SubMatAlign[0][3] = -123.0 + factor;
    SubMatAlign[0][4] = 0.0;

    SubMatAlign[1][0] = -114.0 + factor;
    SubMatAlign[1][1] =  100.0 + factor;
    SubMatAlign[1][2] = -125.0 + factor;
    SubMatAlign[1][3] =  -31.0 + factor;
    SubMatAlign[1][4] = 0.0;

    SubMatAlign[2][0] =  -31.0 + factor;
    SubMatAlign[2][1] = -125.0 + factor;
    SubMatAlign[2][2] =  100.0 + factor;
    SubMatAlign[2][3] = -114.0 + factor;
    SubMatAlign[2][4] = 0.0;

    SubMatAlign[3][0] = -123.0 + factor;
    SubMatAlign[3][1] =  -31.0 + factor;
    SubMatAlign[3][2] = -114.0 + factor;
    SubMatAlign[3][3] =   91.0 + factor;
    SubMatAlign[3][4] = 0.0;

    SubMatAlign[4][0] = 0;
    SubMatAlign[4][1] = 0;
    SubMatAlign[4][2] = 0;
    SubMatAlign[4][3] = 0;
    SubMatAlign[4][4] = 0;

    //Extrapolate MATCH-MISMATCH Differences in the Substitution Weights
    for (unsigned int i=0;i<SubMatAlign.size();i++) {
        for (unsigned int j=0;j<SubMatAlign[i].size();j++) {
            SubMatAlign[i][j] += ( SubMatAlign[i][j] * Coefficient );
        }
    }
    //}

    //4. Init the Sequences //{
    Seqs.clear();
    FSequence::ReadSequencesFromFile(InputFile, Seqs, Alphabet);
    //}

    //5. Initialize the Aligners (one Aligner for each sequence - due to the different weights and penalties) //{
    for (unsigned int i=0;i<Seqs.size();i++) {
        FNW Aligner;
        double factor = (1/Seqs[i].Get_Ratio());

        Aligner.SetParameters(factor*HGapOpPenalty, factor*HGapExPenalty, factor*VGapOpPenalty, factor*VGapExPenalty, TermGapOpPenalty, TermGapExPenalty, SubMatMotif, SubMatAlign);
        Aligners.push_back(Aligner);
    }
    //}

}

ReformAlignInitializer::ReformAlignInitializer(string in, string out, string aligner, bool Verbose) {

    //1. Init Members //{
    InputFile           = in;
    OutputFile          = out;
    ParamsFile          = "";
    SubMatsMotifFile    = "";
    SubMatsAlignFile    = "";
    BaseAligner         = aligner;
    Alphabet            = DNA_Order;

    BonusFactor         = 268.738034793977;
    Coefficient         = 0.2182;
    HGapOpPenalty       = -394.022939886315;
    HGapExPenalty       = -191.612476671029;
    VGapOpPenalty       = -276.993627921082;
    VGapExPenalty       = -126.071721004077;
    TermGapOpPenalty    = 0.0;
    TermGapExPenalty    = 0.0;

    IsVerbose           = Verbose;
    //}

    //2. Init Parameters (automatically based on the APSI of the starting alignment //{
    double APSI  = 0.0;

    if (!BaseAligner.compare("") || !BaseAligner.compare("None") ) {
        cout<<"No starting alignment is defined and thus ReformAlign will now terminate!\n";
        exit(1);
    }

    //Calculate the APSI and Unary SP Score Values based on the starting alignmnet
    AlignmentMetrics metrics(BaseAligner);
    APSI = metrics.Get_APSI();

    Coefficient   = 0.224;
    BonusFactor   = (0.334 * APSI) + 317.206;
    HGapOpPenalty = (0.002 * APSI) - 435.566;
    HGapExPenalty = (0.026 * APSI) - 217.305;
    VGapOpPenalty = (0.1   * APSI) - 301.969;
    VGapExPenalty = (-0.14 * APSI) - 139.698;

    TermGapOpPenalty = 0.0;
    TermGapExPenalty = 0.0;
    //}

    //3. Init Substitution Weights for the Profile Construction Phase from Given Settings//{
    SubMatMotif.clear();
    SubMatMotif.resize(5, vector<double>(5, 0.0));
    SubMatMotif[0][0] =   91.0 + BonusFactor;
    SubMatMotif[0][1] = -114.0 + BonusFactor;
    SubMatMotif[0][2] =  -31.0 + BonusFactor;
    SubMatMotif[0][3] = -123.0 + BonusFactor;
    SubMatMotif[0][4] = 0.0;

    SubMatMotif[1][0] = -114.0 + BonusFactor;
    SubMatMotif[1][1] =  100.0 + BonusFactor;
    SubMatMotif[1][2] = -125.0 + BonusFactor;
    SubMatMotif[1][3] =  -31.0 + BonusFactor;
    SubMatMotif[1][4] = 0.0;

    SubMatMotif[2][0] =  -31.0 + BonusFactor;
    SubMatMotif[2][1] = -125.0 + BonusFactor;
    SubMatMotif[2][2] =  100.0 + BonusFactor;
    SubMatMotif[2][3] = -114.0 + BonusFactor;
    SubMatMotif[2][4] = 0.0;

    SubMatMotif[3][0] = -123.0 + BonusFactor;
    SubMatMotif[3][1] =  -31.0 + BonusFactor;
    SubMatMotif[3][2] = -114.0 + BonusFactor;
    SubMatMotif[3][3] =   91.0 + BonusFactor;
    SubMatMotif[3][4] = 0.0;

    SubMatMotif[4][0] = 0;
    SubMatMotif[4][1] = 0;
    SubMatMotif[4][2] = 0;
    SubMatMotif[4][3] = 0;
    SubMatMotif[4][4] = 0;

    //Extrapolate MATCH-MISMATCH Differences in the Substitution Weights
    for (unsigned int i=0;i<SubMatMotif.size();i++) {
        for (unsigned int j=0;j<SubMatMotif[i].size();j++) {
            SubMatMotif[i][j] += ( SubMatMotif[i][j] * Coefficient );
        }
    }

    //}

    //4. Init Substitution Weights for the Profile Construction Phase from Given Settings//{
    SubMatAlign.clear();
    SubMatAlign.resize(5, vector<double>(5, 0.0));
    SubMatAlign[0][0] =   91.0 + BonusFactor;
    SubMatAlign[0][1] = -114.0 + BonusFactor;
    SubMatAlign[0][2] =  -31.0 + BonusFactor;
    SubMatAlign[0][3] = -123.0 + BonusFactor;
    SubMatAlign[0][4] = 0.0;

    SubMatAlign[1][0] = -114.0 + BonusFactor;
    SubMatAlign[1][1] =  100.0 + BonusFactor;
    SubMatAlign[1][2] = -125.0 + BonusFactor;
    SubMatAlign[1][3] =  -31.0 + BonusFactor;
    SubMatAlign[1][4] = 0.0;

    SubMatAlign[2][0] =  -31.0 + BonusFactor;
    SubMatAlign[2][1] = -125.0 + BonusFactor;
    SubMatAlign[2][2] =  100.0 + BonusFactor;
    SubMatAlign[2][3] = -114.0 + BonusFactor;
    SubMatAlign[2][4] = 0.0;

    SubMatAlign[3][0] = -123.0 + BonusFactor;
    SubMatAlign[3][1] =  -31.0 + BonusFactor;
    SubMatAlign[3][2] = -114.0 + BonusFactor;
    SubMatAlign[3][3] =   91.0 + BonusFactor;
    SubMatAlign[3][4] = 0.0;

    SubMatAlign[4][0] = 0;
    SubMatAlign[4][1] = 0;
    SubMatAlign[4][2] = 0;
    SubMatAlign[4][3] = 0;
    SubMatAlign[4][4] = 0;

    //Extrapolate MATCH-MISMATCH Differences in the Substitution Weights
    for (unsigned int i=0;i<SubMatAlign.size();i++) {
        for (unsigned int j=0;j<SubMatAlign[i].size();j++) {
            SubMatAlign[i][j] += ( SubMatAlign[i][j] * Coefficient );
        }
    }

    //}

    //5. Init the Sequences //{
    Seqs.clear();
    FSequence::ReadSequencesFromFile(InputFile, Seqs, Alphabet);
    //}

    //6. Initialize the Aligners (one Aligner for each sequence - due to the different weights and penalties) //{
    for (unsigned int i=0;i<Seqs.size();i++) {
        FNW Aligner;
        double factor = (1/Seqs[i].Get_Ratio());

        Aligner.SetParameters(factor*HGapOpPenalty, factor*HGapExPenalty, factor*VGapOpPenalty, factor*VGapExPenalty, TermGapOpPenalty, TermGapExPenalty, SubMatMotif, SubMatAlign);
        Aligners.push_back(Aligner);
    }
    //}

}

ReformAlignInitializer::~ReformAlignInitializer() {
    SubMatMotif.clear();
    SubMatAlign.clear();
    Seqs.clear();
    Aligners.clear();
}

//}

/** Core Alignment Functions **/ //{

//Function to Align the Sequences
void ReformAlignInitializer::AlignSequences() {

    //1. Construct Initial Profile //{
    clock_t start, end;
    start = clock();

    //==========
    if (IsVerbose) { cout << "\n\t\t    ReformAlign Started.\n\t\t    --------------------" << "\n\n"; }
    vector<Profile> SeqProfiles; //The Profiles for the Sequences (one profile per Sequence)
    Profile         BasicProfile;//The Basic Profile against which all sequences will finally align against

    //Construct Profile from the Given Alignment File
    Profile::CreateProfileFromAlignmentFile(BasicProfile, BaseAligner, Alphabet);
    Profile::CreateProfilesFromSequences(SeqProfiles, Seqs);
    if (IsVerbose) { cout<<"\r"<<"  1. -> Estimating Initial Profile (Align. File)   - Completed: 100% "<<flush; }
    //cout<<"\n Initially Constructed Profile\n =============================\n"<<BasicProfile.GetProfileToStr(Alphabet, true)<<endl<<endl;
    //}

    //2. Align all Sequences against the Profile //{
    int steps = 1;
    if (IsVerbose) { cout<<"\r"<<"  Kband begin "<<flush; }
    AlignAgainstProfile(BasicProfile, SeqProfiles, steps);
    if (IsVerbose) { cout<<"\r"<<"  Kband sucessfully. "<<flush; }
    //}

    //3. Clean Alignment -  Delete all "Gapped" columns //{
    CleanAlignment();
    //}

    //4. Output new Alignment //{
    steps++;
    if (IsVerbose) { cout<<"  "<<steps<<". -> Saving the Alignment..."<<endl; }
    stringstream result;  result.str("");

    result<<FSequence::GetAlignmentsToStr(Seqs);
    WriteToFile(OutputFile, result.str(), false, true);
    //cout<<result.str()<<endl;

    end = clock();
    if (IsVerbose) { cout << "\n\t\tExecution Time: "<< (double)(end-start)/CLOCKS_PER_SEC<< " seconds." << "\n\n"; }
    if (IsVerbose) { cout<<"Sequences were aligned successfully at file: "<<OutputFile<<endl<<"******ReformAlign will now terminate******\n"; }
    //}

}

//Function to align (in Parallel) each sequence against the Profile (and perform Profile Fine-tuning if required)
void ReformAlignInitializer::AlignAgainstProfile(Profile& BasicProfile, vector<Profile>& SeqProfiles, int& steps) {
    bool ProfileChanged = false;
    IntVec TuningSeqs;
    do {
        if (IsVerbose) { cout<<endl; }
        int counter = 0;
        steps++;
        ProfileChanged = false;
        TuningSeqs.clear();

        //1. Align All Sequences Against the formed Profile //{
        if (IsVerbose) { cout<<"\r"<<"  Align All Sequences Against the formed Profile "<<flush; }
        // #pragma omp parallel for
        for (unsigned int i=0; i<Seqs.size(); i++) {
            if (TuningSeqs.size()==0) {
                // if (IsVerbose) { cout<<"\r"<<"  before getalignment "<<flush; }
                bool ProfileNeedsTuning = Aligners[i].GetAlignment(BasicProfile, SeqProfiles[i], Seqs[i].Get_Alignment(), Alphabet);
                if (IsVerbose) { cout<<"\r"<<"  getalignment suceffuly "<<flush; }

                if (ProfileNeedsTuning) {
                    TuningSeqs.push_back(i);
                }
                // #pragma omp atomic
                counter++;
                if (IsVerbose) {
                    // #pragma omp critical
                    cout<<"\r"<<"  "<<steps<<". -> Aligning against the Profile - Completed: " << (int) (100*(float)(counter)/(float)Seqs.size())<<"% "<<flush;
                }
            }
        }
        //}

        //2. Fine-Tune the Profile to accommodate for more Profile Insertions deriving from the final alignment step //{
        if (TuningSeqs.size()>0) {
            steps++;
            if (IsVerbose) { cout<<endl; }
            ProfileChanged = true;

            for (unsigned int i=0;i<Seqs.size();i++) {
                Aligners[i].GetMotif(BasicProfile, SeqProfiles[i], i);
                if (IsVerbose) { cout<<"\r"<<"  "<<steps<<". -> Fine-tuning the Profile      - Completed: " << (int) (100*(float)(i+1)/(float)Seqs.size())<<"% "<<flush; }
            }
        }
        //}

    } while(ProfileChanged);

    if (IsVerbose) { cout<<endl; }
}

//Function to remove all Gapped Columns from the Alignment
void ReformAlignInitializer::CleanAlignment() {
    int AlignLen = Seqs[0].Get_Alignment().size()-1;
    for (int j=AlignLen; j>=0 ;j--) {
        bool IsGapCol = true;
        for (unsigned int i=0; i<Seqs.size(); i++) {
            if (IsGapCol && Seqs[i].Get_Alignment()[j]!='-') {
                IsGapCol = false;
            }
        }

        if (IsGapCol){
            for (unsigned int i=0; i<Seqs.size(); i++) {
                Seqs[i].Get_Alignment().erase(j,1);
            }
        }
    }
}

//}

/** Output Functions **/ //{

string ReformAlignInitializer::GetParamsToStr() {
    stringstream ss; ss.str("");
    ss << endl;
    ss << "\t\tPARAMETERS\n\t\t==========\n\n";
    ss << "Residue Weight Matrix (Motifs):\n";

    for (unsigned int i=0; i<SubMatMotif.size(); i++) {
        if (i==0) {
            ss << setw(8) << " ";
        }

        ss << setw(8) <<Alphabet.at(i);
    }
    ss<<endl;

    for (unsigned int i=0; i<SubMatMotif.size(); i++) {
        for (unsigned int j=0; j<SubMatMotif[i].size(); j++) {
            if (j==0) {
                ss<<setw(8)<<Alphabet.at(i);
            }
            ss<<setw(8)<<SubMatMotif[i][j];
        }
        ss<<endl;
    }


    ss << "\nResidue Weight Matrix (Aligning):\n";

    for (unsigned int i=0; i<SubMatAlign.size(); i++) {
        if (i==0) {
            ss << setw(8) << " ";
        }
        ss << setw(8) <<Alphabet.at(i);
    }
    ss<<endl;

    for (unsigned int i=0; i<SubMatAlign.size(); i++) {
        for (unsigned int j=0; j<SubMatAlign[i].size(); j++) {
            if (j==0) {
                ss<<setw(8)<<Alphabet.at(i);
            }
            ss<<setw(8)<<SubMatAlign[i][j];
        }
        ss<<endl;
    }

    ss<<endl<<"  H. Gap Opening Penalty : "<< HGapOpPenalty<<endl;
    ss<<" H.Gap Extension Penalty : "<< HGapExPenalty<<endl;
    ss<<" V.Gap Opening   Penalty : "<< VGapOpPenalty<<endl;
    ss<<" V.Gap Extension Penalty : "<< VGapExPenalty<<endl;
    ss<<"Terminal Gap Op. Penalty : "<< TermGapOpPenalty<<endl;
    ss<<"Terminal Gap Ex. Penalty : "<< TermGapExPenalty<<endl;
    ss<<"             Coefficient : "<< Coefficient<<endl;
    ss<<"             BonusFactor : "<< BonusFactor<<endl;

    return ss.str();
}

string ReformAlignInitializer::GetSequencesToStr() {
    stringstream ss; ss.str("");
    ss << endl;
    ss << "\t\tSEQUENCES\n\t\t=========\n\n";
    ss << FSequence::GetSequencesToStr(Seqs);
    return ss.str();
}

string ReformAlignInitializer::GetSettingsToStr() {
    stringstream ss; ss.str("");
    ss << endl;
    ss << GetSequencesToStr();
	ss << GetParamsToStr();

	ss << "             Input File: " << InputFile<<endl;
	ss << "            Output File: " << OutputFile<<endl;
	ss << "        Parameters File: " << ParamsFile<<endl;
	ss << "  Sub. File for Profile: " << SubMatsMotifFile<<endl;
	ss << "Sub. File for Alignment: " << SubMatsAlignFile<<endl;
	ss << "           Base Aligner: " << BaseAligner<<endl;
	ss << "               Alphabet: " << Alphabet<<endl;
	ss << "                   HGOP: " << HGapOpPenalty<<endl;
	ss << "                   HGEP: " << HGapExPenalty<<endl;
	ss << "                   VGOP: " << VGapOpPenalty<<endl;
	ss << "                   VGEP: " << VGapExPenalty<<endl;
	ss << "                   TGOP: " << TermGapOpPenalty<<endl;
	ss << "                   TGEP: " << TermGapExPenalty<<endl;
	ss << "           Bonus Factor: " << BonusFactor<<endl;
	ss << "            Coefficient: " << Coefficient<<endl;
	ss << "     Number of Aligners: " << Aligners.size()<<endl;


    return ss.str();
}

string ReformAlignInitializer::GetSubMatToStr(DoubleMatRef mat) {
    stringstream ss; ss.str("");
    for (unsigned int i=0; i<mat.size(); i++) {
        if (i==0) {
            ss << setw(8) << " ";
        }

        ss << setw(8) <<Alphabet.at(i);
    }
    ss<<endl;

    for (unsigned int i=0; i<mat.size(); i++) {
        for (unsigned int j=0; j<mat[i].size(); j++) {
            if (j==0) {
                ss<<setw(8)<<Alphabet.at(i);
            }
            ss<<setw(8)<<mat[i][j];
        }
        ss<<endl;
    }

    return ss.str();
}

//}

