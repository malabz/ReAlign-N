#include "FSequence.h"

FSequence::FSequence(string id, string seq, string& Alphabet) {
    SequenceID = id;
    Sequence   = seq;
    Ratio      = 1.0;
    for (unsigned int i=0;i<Sequence.size();i++) {
        stringstream ss; ss<<Sequence.at(i);
        SequenceInt.push_back(FirstIndexOf(Alphabet, ss.str()));
    }
}

void FSequence::ReadSequencesFromFile(string SeqFile, SequencesRef Seqs, string& Alphabet) {
    ClearInputStream();

    string line;
    ifstream myfile (SeqFile.c_str());

    stringstream seqId;
    stringstream seqText;
    bool firstLine = true;

    Alphabet = "";

    if (myfile.is_open()) {
        while ( myfile.good() ) {
            getline (myfile,line);

            ReplaceAll(line, "\r", "");
            ReplaceAll(line, "\n", "");
            ReplaceAll(line, "\r\n", "");
            ReplaceAll(line, "\n\r", "");

            if (line.compare("")==0 || line.compare("\n")==0 || line.compare(" ")==0) continue; //skip empty lines
            if (StartsWith(line, ">") || StartsWith(line,";")) {
                if (!firstLine) {
                    //Convert Sequences to Uppercase
                    string str = seqText.str();
                    std::transform(str.begin(), str.end(),str.begin(), ::toupper);

                    if (!Alphabet.compare("")) {
                        if (Contains(str, "U")) {
                            stringstream st(RNA_Order);
                            Alphabet = st.str();
                        }else if (Contains(str, "R") || Contains(str, "N") || Contains(str, "D") ||
                                  Contains(str, "Q") || Contains(str, "E") || Contains(str, "H") ||
                                  Contains(str, "I") || Contains(str, "L") || Contains(str, "K") || Contains(str, "M") ||
                                  Contains(str, "F") || Contains(str, "P") || Contains(str, "S") ||
                                  Contains(str, "W") || Contains(str, "Y") || Contains(str, "V") || Contains(str, "B") ||
                                  Contains(str, "Z") || Contains(str, "X")
                                 ) {
                            stringstream st(AA_Order);
                            Alphabet = st.str();
                        }else {
                            stringstream st(DNA_Order);
                            Alphabet = st.str();
                        }
                    }
                    FSequence seq(seqId.str().substr(1), str, Alphabet);
                    Seqs.push_back(seq);
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

        FSequence seq(seqId.str().substr(1), str, Alphabet);
        Seqs.push_back(seq);
    } else {
        cerr << " Unable to open Sequences file: "<<SeqFile<<" for reading.";
        exit(1);
    }
    myfile.close();
    UpdateSequencesRatio(Seqs);
}

void FSequence::ReadAlignmentFromFile(string SeqFile, SequencesRef Seqs, string& Alphabet) {
    ClearInputStream();

    string line;
    ifstream myfile (SeqFile.c_str());

    stringstream seqId;
    stringstream seqText;
    bool firstLine = true;

    Alphabet = "";

    if (myfile.is_open()) {
        while ( myfile.good() ) {
            getline (myfile,line);

            ReplaceAll(line, "\r", "");
            ReplaceAll(line, "\n", "");
            ReplaceAll(line, "\r\n", "");
            ReplaceAll(line, "\n\r", "");

            if (line.compare("")==0 || line.compare("\n")==0 || line.compare(" ")==0) continue; //skip empty lines
            if (StartsWith(line, ">") || StartsWith(line,";")) {
                if (!firstLine) {
                    //Convert Sequences to Uppercase
                    string str = seqText.str();
                    std::transform(str.begin(), str.end(),str.begin(), ::toupper);

                    if (!Alphabet.compare("")) {
                        if (Contains(str, "U")) {
                            stringstream st(RNA_Order);
                            Alphabet = st.str();
                        }else if (Contains(str, "R") || Contains(str, "N") || Contains(str, "D") ||
                                  Contains(str, "Q") || Contains(str, "E") || Contains(str, "H") ||
                                  Contains(str, "I") || Contains(str, "L") || Contains(str, "K") || Contains(str, "M") ||
                                  Contains(str, "F") || Contains(str, "P") || Contains(str, "S") ||
                                  Contains(str, "W") || Contains(str, "Y") || Contains(str, "V") || Contains(str, "B") ||
                                  Contains(str, "Z") || Contains(str, "X")
                                 ) {
                            stringstream st(AA_Order);
                            Alphabet = st.str();
                        }else {
                            stringstream st(DNA_Order);
                            Alphabet = st.str();
                        }
                    }
                    string clearedStr(str);
                    ReplaceAll(clearedStr, "-", "");
                    FSequence seq(seqId.str().substr(1), clearedStr, Alphabet);
                    seq.Alignment = str;
                    Seqs.push_back(seq);
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

        string clearedStr(str);
        ReplaceAll(clearedStr, "-", "");
        FSequence seq(seqId.str().substr(1), clearedStr, Alphabet);
        seq.Alignment = str;
        Seqs.push_back(seq);
    } else {
        cerr << " Unable to open Alignment file: "<<SeqFile<<" for reading.";
        exit(1);
    }
    myfile.close();
}


void FSequence::UpdateSequencesRatio(SequencesRef Seqs) {
    //1. Find Max Size //{
    unsigned int MaxSize = Seqs[0].Get_Sequence().size();
    for (unsigned int i=0;i<Seqs.size();i++) {
        if (Seqs[i].Get_Sequence().size()>MaxSize) {
            MaxSize = Seqs[i].Get_Sequence().size();
        }
    }
    //}

    //2. Update the Ratio for all Sequences //{
    for (unsigned int i=0;i<Seqs.size();i++) {
        Seqs[i].Ratio = (double) Seqs[i].Get_Sequence().size() / MaxSize;
    }
    //}
}

double FSequence::GetAverageRatio(SequencesRef Seqs) {
    double AvgRatio = 0.0;
    for (unsigned int i=0;i<Seqs.size();i++) {
        AvgRatio += Seqs[i].Ratio;
    }
    AvgRatio /= Seqs.size();
    return AvgRatio;
}


int FSequence::GetPosOfSeqWithID(SequencesRef Seqs, string id) {
    for (unsigned int i=0;i<Seqs.size();i++) {
        if (Seqs[i].Get_SequenceID().compare(id)==0)
            return i;
    }
    return -1;
}

bool FSequence::AreAlignmentsIdentical(SequencesRef Seqs1, SequencesRef Seqs2) {
    for (unsigned int i=0;i<Seqs1.size();i++) {
        if (Seqs1[i].Get_Alignment().compare(Seqs2[i].Get_Alignment()))
            return false;
    }
    return true;
}



string FSequence::GetAlignmentsToStr(SequencesRef Seqs) {
    stringstream ss; ss.str("");
    for (unsigned int i=0;i<Seqs.size();i++) {
        ss << ">" << Seqs[i].Get_SequenceID() <<endl;
        ss << Seqs[i].Get_Alignment()<<endl;
    }
    return ss.str();
}

string FSequence::GetSequencesToStr(SequencesRef Seqs) {
    stringstream ss; ss.str("");
    for (unsigned int i=0;i<Seqs.size();i++) {
        ss << ">" << Seqs[i].Get_SequenceID() <<endl;
        ss << Seqs[i].Get_Sequence()<<endl<<endl;
    }
    return ss.str();
}

string FSequence::GetSequencesIntToStr(SequencesRef Seqs) {
    stringstream ss; ss.str("");
    for (unsigned int i=0;i<Seqs.size();i++) {
        ss << ">" << Seqs[i].Get_SequenceID() <<endl;
        for (unsigned int j=0;j<Seqs[i].Get_SequenceInt().size();j++) {
            ss <<setw(3) << Seqs[i].Get_SequenceInt()[j];
        }
        ss << endl;
    }
    ss<<endl;
    return ss.str();
}

string FSequence::GetRatiosToStr(SequencesRef Seqs) {
    stringstream ss; ss.str("");
    for (unsigned int i=0;i<Seqs.size();i++) {
        ss << ">" << setw(50)<< Seqs[i].Get_SequenceID() <<"\tRatio: "<<setw(12)<<Seqs[i].Get_Ratio()<<endl;
    }
    return ss.str();
}
