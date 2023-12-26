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

#ifndef PROFILE_H
#define PROFILE_H

#include "Utils.h"
#include "TypeDefinitions.h"
#include "FSequence.h"

using namespace Utils;
using namespace TypeDefinitions;

class Profile {
    class Site {
    public:
        int    Weight;      //The number of sequences that Aligned to this letter
        IntMap Supporters;  //A map<int, int> with the key being the i value of the sequence that supported this Pattern (used to decide on whether Weights should be updated during fine-tuning)

        //The Default Constructor
        Site() {  }

        //The Copy Constructor
        Site(const Site& other) {
            Weight = other.Weight;
            for ( auto it = other.Supporters.begin(); it != other.Supporters.end(); ++it ) {
                Supporters[(*it).first] = (*it).second;
            }
        }

        //The Assignment Operator
        Site & operator= (const Site& other) {
            if (this != &other) { // protect against invalid self-assignment
                Weight = other.Weight;
                for ( auto it = other.Supporters.begin(); it != other.Supporters.end(); ++it ) {
                    Supporters[(*it).first] = (*it).second;
                }
            }
            return *this;
        }

        //The Deconstructor
        virtual ~Site() {
            Supporters.clear();
        }
    };

public:
    typedef map<int, Site>  Pattern;
    typedef map<int, Site>& PatternRef;

    //The Default Constructor
    Profile() { };

    //The Copy Constructor
    Profile(const Profile& other) {
        for (unsigned int i=0;i<other.Patterns.size();i++) {
            Patterns.push_back(other.Patterns[i]);
        }
    }

    //The Assignment Operator
    Profile & operator= (const Profile& other) {
        if (this != &other) { // protect against invalid self-assignment
            for (unsigned int i=0;i<other.Patterns.size();i++) {
                Patterns.push_back(other.Patterns[i]);
            }
        }
        return *this;
    }

    //The Deconstructor
    virtual ~Profile() {
        Patterns.clear();
    }

    vector<Pattern>& GetPatterns()                 { return Patterns; }
    PatternRef       GetPatternAt(int PatternPos)  { return Patterns[PatternPos]; }
    int              GetPatternsSize()             { return Patterns.size(); }

    int       NumOfSitesAt(int PatternPos)  { return Patterns[PatternPos].size(); }
    string    GetLetterAt(int Letter, string Alphabet) { stringstream ss; ss<<(Alphabet.at(Letter)); return ss.str(); }
    int       GetWeightAt(int PatternPos, int Letter)  { return Patterns[PatternPos][Letter].Weight; }
    IntMapRef GetSupportersAt(int PatternPos, int Letter)  { return Patterns[PatternPos][Letter].Supporters; }
    int       GetNumOfSupportersAt(int PatternPos, int Letter) { return Patterns[PatternPos][Letter].Supporters.size(); }
    bool      IsSiteSupportedBy(int PatternPos, int Letter, int SupporterID);
    bool      PatternContainsSite(int PatternPos, int Letter);
    double    GetMaxIdentityAt(int PatternPos);
    double    GetIdentityAt(int PatternPos, int Letter);

    void   IncreaseWeight(int PatternPos, int Letter, int IncrBy, int SupporterID);
    void   DecreaseWeight(int PatternPos, int Letter, int DecrBy, int SupporterID);
    int    GetSumOfWeightsAt(int PatternPos);

    void   DeleteSiteAt(int PatternPos, int Letter) { Patterns[PatternPos].erase(Letter); }
    void   InsertSiteAt(int PatternPos, int Letter, int Weight, int SupporterID) { Site NewSite; NewSite.Weight = Weight; NewSite.Supporters[SupporterID] = 1; Patterns[PatternPos][Letter] = NewSite; }

    void   DeletePatternAt(int PatternPos) { Patterns.erase(Patterns.begin()+PatternPos); }
    void   InsertPatternAt(int PatternPos, Pattern NewPattern) { Patterns.insert(Patterns.begin()+PatternPos, NewPattern); }
    void   PushBackPattern(Pattern NewPattern)    { Patterns.push_back(NewPattern); }

    static void   CreateProfileFromSequence(Profile& Prof, FSequence& Seq, int SupporterID);
    static void   CreateProfilesFromSequences(vector<Profile>& Profiles, vector<FSequence>& Seqs);
    static void   CreateProfileFromAlignments(Profile& Prof, StringVecRef Alignments, string Alphabet);
    static void   CreateProfileFromAlignmentFile(Profile& Prof, string AlignmentFileName, string Alphabet);
    static bool   PatternContainsSite(PatternRef Patt, int Letter);
    static void   CreatePatternCopy(PatternRef OriginalPattern, PatternRef CopyPattern);
    static void   CreateProfileCopy(Profile& OriginalProfile, Profile& CopyProfile);
    static void   CreateSubProfile(Profile& OriginalProfile, Profile& SubProfile, int fromPos, int toPos);
    static void   CreateSubProfile(Profile& OriginalProfile, Profile& SubProfile, int fromPos, int toPos, bool ResetWeights);
    static double GetScoreForAlignment(Profile& Prof, string& AlignmentStr, DoubleMatRef SubMat, string& Alphabet);

    static string CompareProfiles(Profile& Prof1, Profile& Prof2, string& Alphabet);

    string GetSupportersToStrAt(int PatternPos);
    string GetSupportersToStrAt(int PatternPos, int Letter);
    string GetLettersToStrAt(int PatternPos);
    string GetLettersToStrAt(int PatternPos, string Alphabet);
    string GetWeightsToStrAt(int PatternPos);
    string GetSitesToStrAt(int PatternPos);
    string GetSitesToStrAt(int PatternPos, string Alphabet);
    string GetProfileToStr(string Alphabet);
    string GetProfileToStr(string Alphabet, bool ShowSupporters);

protected:

private:
    vector<Pattern> Patterns; //The key in the map is the Alphabet letter (e.g. A->0, C->1, G->2, T->3 for DNA alphabet) of this Pattern
};

/** Definition of the toString Operation **/
ostream& operator<<(ostream& ostr, const Profile& prof);

#endif // PROFILE_H
