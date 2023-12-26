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

#include "AlignmentMetrics.h"

AlignmentMetrics::AlignmentMetrics(string SeqsFileName) {
    ClearInputStream();
    StringVec SequenceIds;
    StringVec Alignments;

    /**1. Read Starting Alignment and Sequence Ids from Input File**/ //{
    string line;
    ifstream myfile (SeqsFileName.c_str());
    string filename = SeqsFileName;

    stringstream seqId;
    stringstream seqText;
    bool firstLine = true;

    if (myfile.is_open()) {
        while ( myfile.good() ) {
            getline (myfile,line);

            ReplaceAll(line, "\r", "");
            ReplaceAll(line, "\n", "");
            ReplaceAll(line, "\r\n", "");
            ReplaceAll(line, "\n\r", "");
            if (line.compare("")==0 || line.compare("\r")==0 || line.compare("\r\n")==0   || line.compare("\n")==0 || line.compare(" ")==0) { continue; } //skip empty lines

            if (StartsWith(line, ">") || StartsWith(line,";")) {
                if (!firstLine) {
                    Alignments.push_back(seqText.str());
                    SequenceIds.push_back(seqId.str().substr(1));
                    seqId.str("");
                    seqText.str("");
                }
                seqId << line;
                firstLine = false;
            } else {
                seqText << line;
            }

        }
        Alignments.push_back(seqText.str());
        SequenceIds.push_back(seqId.str().substr(1));
    } else cerr << " Unable to open file: "<<filename<<" for reading.";
    myfile.close();
    //}

    /**2. Fix Alignment (trim blank lines, convert toupper etc) **/ //{
    for (unsigned int i=0; i<SequenceIds.size(); i++) {
        SequenceIds[i] = Trim(SequenceIds[i]);
        Alignments[i] = Trim(Alignments[i]);
        Alignments[i] = ConvertToUpperCase(Alignments[i]);

        for (unsigned int j=0; j<Alignments[i].size(); j++) {
            if (Alignments[i][j]!='A' && Alignments[i][j]!='C' && Alignments[i][j]!='G' &&
                    Alignments[i][j]!='T' && Alignments[i][j]!='U' && Alignments[i][j]!='-') {
            }
        }
    }
    //}

    /**3. Calculate the APSI Value **/ //{
    APSI = GetAPSI(SequenceIds, Alignments);
    //}

}

AlignmentMetrics::~AlignmentMetrics() { }

int AlignmentMetrics::min(int a, int b) {
    if (a<b) {
        return a;
    } else {
        return b;
    }
}

double AlignmentMetrics::average(int a, int b) {
    return (a+b)/2.0f;
}

double AlignmentMetrics::GetAPSI(StringVecRef SequenceIds, StringVecRef Alignments) {
    double ovpid1 = 0.0;
    double ovpid2 = 0.0;
    double ovpid3 = 0.0;

    int counter = 0;

    for (unsigned int k=0; k<SequenceIds.size()-1; k++) {
        string fSeq    = Alignments[k];
        string fSeqCpy = Alignments[k];
        ReplaceAll(fSeqCpy, "-", "");
        int fSeqLen  = fSeqCpy.size();
        int alignLen = fSeq.size();

        for (unsigned int l=k+1; l<SequenceIds.size(); l++) {
            string sSeq    = Alignments[l];
            string sSeqCpy = Alignments[l];

            ReplaceAll(sSeqCpy, "-", "");
            int sSeqLen = sSeqCpy.size();
            int minLen = min(fSeqLen, sSeqLen);

            double avgLen = average(fSeqLen, sSeqLen);

            int Pid1Count = 0, Pid3Count = 0;
            for (int i=0; i<alignLen; i++) {
                if (fSeq.at(i)!='-' && sSeq.at(i)!='-') {
                    if (fSeq.at(i)==sSeq.at(i)) {
                        Pid1Count++;
                    }
                }

                if (fSeq.at(i)==sSeq.at(i)) {
                    Pid3Count++;
                }
            }
            double pid1 = 100.0*(Pid1Count/(double)minLen);
            double pid2 = 100.0*(Pid1Count/(double)avgLen);
            double pid3 = 100.0*(Pid3Count/(double)alignLen);

            ovpid1 += pid1;
            ovpid2 += pid2;
            ovpid3 += pid3;
            counter++;
        }
    }
    ovpid1 /= counter;
    ovpid2 /= counter;
    ovpid3 /= counter;

    return (ovpid1+ovpid2+ovpid3)/3.0;
}
