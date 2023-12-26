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

#include "Profile.h"


bool Profile::IsSiteSupportedBy(int PatternPos, int Letter, int SupporterID) {
    if (Patterns[PatternPos][Letter].Supporters.empty()) return false;

    IntMap::iterator it = Patterns[PatternPos][Letter].Supporters.find(SupporterID);
    if(it != Patterns[PatternPos][Letter].Supporters.end()){
        return true;
    }else{
        return false;
    }
}

bool Profile::PatternContainsSite(PatternRef Patt, int Letter) {
    if (Patt.empty() ) return false;

    Pattern::iterator it = Patt.find(Letter);
    if(it != Patt.end()){
        return true;
    }else{
        return false;
    }
}

bool Profile::PatternContainsSite(int PatternPos, int Letter) {
    if (Patterns[PatternPos].empty() ) return false;

    Pattern::iterator it = Patterns[PatternPos].find(Letter);
    if(it != Patterns[PatternPos].end()){
        return true;
    }else{
        return false;
    }
}

void Profile::IncreaseWeight(int PatternPos, int Letter, int IncrBy, int SupporterID) {
    Patterns[PatternPos][Letter].Weight += IncrBy;            //Increase Weight by the given value
    Patterns[PatternPos][Letter].Supporters[SupporterID] = 1; //Add this sequence as supporter to this site. The Supporters[] value (1) is chosen randmly and does not play any role.
}

void Profile::DecreaseWeight(int PatternPos, int Letter, int DecrBy, int SupporterID) {
    Patterns[PatternPos][Letter].Weight -= DecrBy;                //Increase Weight by the given value
    Patterns[PatternPos][Letter].Supporters.erase (SupporterID);  //Remove this Supporter from the Map of Supporters
}

int Profile::GetSumOfWeightsAt(int PatternPos) {
    int sum = 0;
    for ( auto it = Patterns[PatternPos].begin(); it != Patterns[PatternPos].end(); ++it ) {
        sum+=(*it).second.Weight;
    }
    return sum;
}

double Profile::GetMaxIdentityAt(int PatternPos) {
    int MaxWeight = 0;
    for ( auto it = Patterns[PatternPos].begin(); it != Patterns[PatternPos].end(); ++it ) {
        if ( (*it).second.Weight > MaxWeight  ) {
            MaxWeight = (*it).second.Weight;
        }
    }

    return (double)MaxWeight/GetSumOfWeightsAt(PatternPos);
}

double Profile::GetIdentityAt(int PatternPos, int Letter) {
    return (double)Patterns[PatternPos][Letter].Weight/GetSumOfWeightsAt(PatternPos);
}


void Profile::CreateProfileFromSequence(Profile& Prof, FSequence& Seq, int SupporterID) {
    for (unsigned int i=0;i<Seq.Get_SequenceInt().size();i++) {
        Pattern tempPattern;
        tempPattern[Seq.Get_SequenceInt()[i]].Weight = 1;
        tempPattern[Seq.Get_SequenceInt()[i]].Supporters[SupporterID] = 1;
        Prof.PushBackPattern(tempPattern);
    }
}

void Profile::CreateProfilesFromSequences(vector<Profile>& Profiles, vector<FSequence>& Seqs) {
    for (unsigned int i=0;i<Seqs.size();i++) {
        Profile tempProf;
        CreateProfileFromSequence(tempProf, Seqs[i], i);
        Profiles.push_back(tempProf);
    }
}

void Profile::CreateProfileFromAlignments(Profile& Prof, StringVecRef Alignments, string Alphabet) {
    for (unsigned int j=0;j<Alignments[0].size();j++) {
        Pattern tempPattern;
        for (unsigned int i=0;i<Alignments.size();i++) {
            if (Alignments[i][j]=='-') continue; //skip the gaps
            stringstream ss; ss<<Alignments[i].at(j);
            int Letter = FirstIndexOf(Alphabet, ss.str());
            if (Profile::PatternContainsSite(tempPattern, Letter)) {
                tempPattern[Letter].Weight++;
            }else {
                tempPattern[Letter].Weight=1;
            }
            tempPattern[Letter].Supporters[i] = 1;
        }
        Prof.PushBackPattern(tempPattern);
    }
}

void Profile::CreateProfileFromAlignmentFile(Profile& Prof, string AlignmentFileName, string Alphabet) {
    ClearInputStream();

    string line;
    ifstream myfile (AlignmentFileName.c_str());

    stringstream seqId;
    stringstream seqText;
    bool firstLine = true;

    StringVec Alignments;

    if (myfile.is_open()) {
        while ( myfile.good() ) {
            getline (myfile,line);

            if (line.compare("")==0 || line.compare("\n")==0 || line.compare(" ")==0) continue; //skip empty lines
            if (StartsWith(line, ">") || StartsWith(line,";")) {
                if (!firstLine) {
                    //Convert Sequences to Uppercase
                    string str = seqText.str();
                    std::transform(str.begin(), str.end(),str.begin(), ::toupper);

                    if (!Alphabet.compare("")) {
                        if (Contains(str, "U")) {
                            Alphabet = RNA_Order;
                        }else if (Contains(str, "R") || Contains(str, "N") || Contains(str, "D") ||
                                  Contains(str, "Q") || Contains(str, "E") || Contains(str, "H") ||
                                  Contains(str, "I") || Contains(str, "L") || Contains(str, "K") || Contains(str, "M") ||
                                  Contains(str, "F") || Contains(str, "P") || Contains(str, "S") ||
                                  Contains(str, "W") || Contains(str, "Y") || Contains(str, "V") || Contains(str, "B") ||
                                  Contains(str, "Z") || Contains(str, "X")
                                ) {
                            Alphabet = AA_Order;
                        }else {
                            Alphabet = DNA_Order;
                        }
                    }

                    Alignments.push_back(str);
                    seqId.str("");
                    seqText.str("");
                }
                seqId << line;
                firstLine = false;
            } else {
                seqText << line;
            }

        }
        //Convert Sequences to Uppercase
        string str = seqText.str();
        std::transform(str.begin(), str.end(),str.begin(), ::toupper);

        Alignments.push_back(str);
    } else {
        cerr << " Unable to open Sequences file: "<<AlignmentFileName<<" for reading.";
        exit(1);
    }
    myfile.close();

    CreateProfileFromAlignments(Prof, Alignments, Alphabet);
}

void Profile::CreatePatternCopy(PatternRef OriginalPattern, PatternRef CopyPattern) {
    CopyPattern.clear();
    for ( auto it = OriginalPattern.begin(); it != OriginalPattern.end(); ++it ) {
        Site tempSite;
        tempSite.Weight = (*it).second.Weight;
        for ( auto supit = (*it).second.Supporters.begin(); supit != (*it).second.Supporters.end(); ++supit ) {
            tempSite.Supporters[(*supit).first] = 1;
        }
        CopyPattern[(*it).first] = tempSite;
    }
}

void Profile::CreateProfileCopy(Profile& OriginalProfile, Profile& CopyProfile) {
    CopyProfile.Patterns.clear();
    for (int i=0;i<OriginalProfile.GetPatternsSize();i++) {
        Pattern PatternCopy;
        CreatePatternCopy(OriginalProfile.GetPatternAt(i), PatternCopy);
        CopyProfile.PushBackPattern(PatternCopy);
    }
}

void Profile::CreateSubProfile(Profile& OriginalProfile, Profile& SubProfile, int fromPos, int toPos) {
    SubProfile.Patterns.clear();
    for (int i=fromPos;i<=toPos;i++) {
        Pattern PatternCopy;
        CreatePatternCopy(OriginalProfile.GetPatternAt(i), PatternCopy);
        SubProfile.PushBackPattern(PatternCopy);
    }
}

void Profile::CreateSubProfile(Profile& OriginalProfile, Profile& SubProfile, int fromPos, int toPos, bool ResetWeights) {
    SubProfile.Patterns.clear();
    for (int i=fromPos;i<=toPos;i++) {
        Pattern PatternCopy;
        if (!ResetWeights) {
            CreatePatternCopy(OriginalProfile.GetPatternAt(i), PatternCopy);
        }else {
            PatternRef OriginalPattern = OriginalProfile.GetPatternAt(i);
            PatternRef CopyPattern = PatternCopy;
            CopyPattern.clear();
            for ( auto it = OriginalPattern.begin(); it != OriginalPattern.end(); ++it ) {
                Site tempSite;
                tempSite.Weight = 1;
                for ( auto supit = (*it).second.Supporters.begin(); supit != (*it).second.Supporters.end(); ++supit ) {
                    tempSite.Supporters[(*supit).first] = 1;
                }
                CopyPattern[(*it).first] = tempSite;
            }
        }
        SubProfile.PushBackPattern(PatternCopy);
    }
}

double Profile::GetScoreForAlignment(Profile& Prof, string& AlignmentStr, DoubleMatRef SubMat, string& Alphabet) {
    double score = 0.0f;
    for (unsigned int i=0;i<AlignmentStr.size();i++) {
        if (AlignmentStr.at(i)=='-') continue; //skip all gaps
        stringstream ss; ss<<AlignmentStr.at(i);
        PatternRef Pat = Prof.GetPatternAt(i);

        for ( auto it = Pat.begin(); it != Pat.end(); ++it ) {
            score += ((*it).second.Weight) * SubMat[(*it).first][FirstIndexOf(Alphabet, ss.str())];
        }
    }
    return score;
}


string Profile::CompareProfiles(Profile& Prof1, Profile& Prof2, string& Alphabet) {
    stringstream ss;
    ss.str("");
    int Prof1Size = Prof1.GetPatternsSize();
    int Prof2Size = Prof2.GetPatternsSize();
    if (Prof1Size!=Prof2Size) {
        ss<<"PROFILES DIFFER IN LENGTH: PROFILE 1 SIZE ["<<Prof1Size<<"] PROFILE 2 SIZE ["<<Prof2Size<<endl;
    }
    for (int i=0;i<max(Prof1Size, Prof2Size);i++) {
        string PatternAtProf1 = "", PatternAtProf2 = "";
        if (i<Prof1Size) {
            PatternAtProf1 = Prof1.GetSitesToStrAt(i, Alphabet);
        }

        if (i<Prof2Size) {
            PatternAtProf2 = Prof2.GetSitesToStrAt(i, Alphabet);
        }

        if(PatternAtProf1.compare(PatternAtProf2)) {
            //ss<<"PATTERN DIFFERENCES at Position: "<<(i)<<endl;
            //ss<<"\tPROFILE 1:\n"<<PatternAtProf1<<endl;
            //ss<<"\tPROFILE 2:\n"<<PatternAtProf2<<endl;
        }
    }
    if (!ss.str().compare("")) {
        ss<<"The Profiles are IDENTICAL!\n";
    }
    return ss.str();
}

string Profile::GetSupportersToStrAt(int PatternPos) {
    stringstream ss;
    ss.str("");
    for ( auto it = Patterns[PatternPos].begin(); it != Patterns[PatternPos].end(); ++it ) {
        if (it != Patterns[PatternPos].begin()) {
            ss<<"\t            ";
        }
        ss<<"Pattern: "<<(*it).first<<" -Supported By Seqs:";
        ss<<GetSupportersToStrAt(PatternPos, (*it).first);
        ss<<endl;
    }
    return ss.str();
}

string Profile::GetSupportersToStrAt(int PatternPos, int Letter) {
    stringstream ss;
    ss.str("");
    for ( auto it = Patterns[PatternPos][Letter].Supporters.begin(); it != Patterns[PatternPos][Letter].Supporters.end(); ++it ) {
        ss<<setw(5)<<(*it).first<<" | ";
    }
    return ss.str().substr(0, ss.str().size()-3);
}

string Profile::GetLettersToStrAt(int PatternPos) {
    stringstream ss;
    ss.str("");
    for ( auto it = Patterns[PatternPos].begin(); it != Patterns[PatternPos].end(); ++it ) {
        ss<<setw(5)<<(*it).first<<" | ";
    }
    return ss.str().substr(0, ss.str().size()-3);
}

string Profile::GetLettersToStrAt(int PatternPos, string Alphabet) {
    stringstream ss;
    ss.str("");
    for ( auto it = Patterns[PatternPos].begin(); it != Patterns[PatternPos].end(); ++it ) {
        ss<<setw(5)<<GetLetterAt((*it).first, Alphabet)<<" | ";
    }
    return ss.str().substr(0, ss.str().size()-3);
}

string Profile::GetWeightsToStrAt(int PatternPos) {
    stringstream ss;
    ss.str("");
    for ( auto it = Patterns[PatternPos].begin(); it != Patterns[PatternPos].end(); ++it ) {
        ss<<setw(5)<<(*it).second.Weight<<" | ";
    }
    return ss.str().substr(0, ss.str().size()-3);
}

string Profile::GetSitesToStrAt(int PatternPos) {
    stringstream ss1, ss2;
    ss1.str("");
    for ( auto it = Patterns[PatternPos].begin(); it != Patterns[PatternPos].end(); ++it ) {
        ss1<<setw(5)<<(*it).first<<" | ";
    }

    ss2<<ss1.str().substr(0, ss1.str().size()-3)<<endl;
    ss1.str("");
    for ( auto it = Patterns[PatternPos].begin(); it != Patterns[PatternPos].end(); ++it ) {
        ss1<<setw(5)<<(*it).second.Weight<<" | ";
    }
    ss2<<ss1.str().substr(0, ss1.str().size()-3)<<endl;
    return ss2.str();
}

string Profile::GetSitesToStrAt(int PatternPos, string Alphabet) {
    stringstream ss1, ss2;
    ss1.str("");
    for ( auto it = Patterns[PatternPos].begin(); it != Patterns[PatternPos].end(); ++it ) {
        ss1<<setw(5)<<GetLetterAt((*it).first, Alphabet)<<" | ";
    }

    ss2<<ss1.str().substr(0, ss1.str().size()-3)<<endl;
    ss1.str("");
    for ( auto it = Patterns[PatternPos].begin(); it != Patterns[PatternPos].end(); ++it ) {
        ss1<<setw(5)<<(*it).second.Weight<<" | ";
    }
    ss2<<ss1.str().substr(0, ss1.str().size()-3)<<endl;
    return ss2.str();
}

string Profile::GetProfileToStr(string Alphabet) {
    stringstream ss;
    ss.str("");
    for (unsigned int i=0;i<Patterns.size();i++) {
        ss<<GetLettersToStrAt(i, Alphabet)<<"\t";
    }
    ss<<endl;
    for (unsigned int i=0;i<Patterns.size();i++) {
        ss<<GetWeightsToStrAt(i)<<"\t";
    }
    return ss.str();
}

string Profile::GetProfileToStr(string Alphabet, bool ShowSupporters) {
    stringstream ss;
    ss.str("");
    for (unsigned int i=0;i<Patterns.size();i++) {
        ss<<"In Profile position ["<<setw(5)<<i<<"]"<<endl;
        ss<<"\t  Patterns: "<<GetLettersToStrAt(i, Alphabet)<<endl;
        ss<<"\t   Weights: "<<GetWeightsToStrAt(i)<<endl;
        ss<<"\tSupporters: "<<GetSupportersToStrAt(i)<<endl;
    }
    ss<<endl;
    return ss.str();
}


/** Definition of the toString Operation **/
ostream& operator<<(ostream& ostr, const Profile& prof) {
    ostr<<endl<<"Not Implemented"<<endl;
    return ostr;
}
