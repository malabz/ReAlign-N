#ifndef FSEQUENCE_H
#define FSEQUENCE_H

#include "TypeDefinitions.h"
#include "Tokenizer.h"
#include "Utils.h"

using namespace std;
using namespace TypeDefinitions;
using namespace Utils;

class FSequence {
    typedef vector<FSequence>  Sequences;
    typedef vector<FSequence>& SequencesRef;

    public:
        FSequence(string id, string seq, string& SequencesType);

        //The Default Constructor
        FSequence() { };

        //The Copy Constructor
        FSequence(const FSequence& other) {
            SequenceID = other.SequenceID;
            Sequence  = other.Sequence;
            Alignment  = other.Alignment;
            Ratio  = other.Ratio;

            SequenceInt.resize(other.SequenceInt.size(), 0);
            for (unsigned int i=0;i<other.SequenceInt.size();i++) {
                SequenceInt[i] = other.SequenceInt[i];
            }
        }

        //The Assignment Operator
        FSequence & operator= (const FSequence& other) {
            if (this != &other) { // protect against invalid self-assignment
                SequenceID = other.SequenceID;
                Sequence  = other.Sequence;
                Alignment  = other.Alignment;
                Ratio  = other.Ratio;

                SequenceInt.resize(other.SequenceInt.size(), 0);
                for (unsigned int i=0;i<other.SequenceInt.size();i++) {
                    SequenceInt[i] = other.SequenceInt[i];
                }
            }
            return *this;
        }

        //The Deconstructor
        virtual ~FSequence() {
            SequenceInt.clear();
        }



        bool operator < (const FSequence& anotherSeq) const {
            return (Sequence.size() < anotherSeq.Sequence.size());
        }

        string&     Get_SequenceID()  { return SequenceID;  }
        string&     Get_Sequence()    { return Sequence;    }
        string&     Get_Alignment()   { return Alignment;   }
        double&     Get_Ratio()       { return Ratio;       }
        IntVecRef   Get_SequenceInt() { return SequenceInt; }

        static void   ReadSequencesFromFile(string SeqFile, SequencesRef Seqs, string& SequencesType);
        static void   ReadAlignmentFromFile(string SeqFile, SequencesRef Seqs, string& SequencesType);
        static void   UpdateSequencesRatio(SequencesRef Seqs);
        static double GetAverageRatio(SequencesRef Seqs);

        static int    GetPosOfSeqWithID(SequencesRef Seqs, string id);
        static bool   AreAlignmentsIdentical(SequencesRef Seqs1, SequencesRef Seqs2);

        static string GetSequencesToStr(SequencesRef Seqs);
        static string GetAlignmentsToStr(SequencesRef Seqs);
        static string GetSequencesIntToStr(SequencesRef Seqs);
        static string GetRatiosToStr(SequencesRef Seqs);

    protected:

    private:
        string      SequenceID;
        string      Sequence;
        string      Alignment;
        IntVec      SequenceInt;
        double      Ratio; //Representing how smaller each sequence is compared to the largest sequence (#of Residues in Seq/ Max # of Residues in all Sequences)
};

#endif // FSEQUENCE_H
