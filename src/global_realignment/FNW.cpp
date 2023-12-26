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


#include "FNW.h"

// K-band's utils begin {
inline double FNW::maxi(double a, double b) {
    if (a > b)return a;
    else return b;
}
inline double FNW::maxi(double a, double b, double c) {
    double max; if (a > b)max = a; else max = b;
    if (max > c)return max; else return c;
}


inline bool FNW::InsideStrip(int i, int j, int k, int diff) {
    if( j - i >= -k and j - i <= k + diff ) {
        return true;
    } else {
        return false;
    }
}

void FNW::Init(int m, int k, int diff, CharMatRef bt, DoubleMatRef pm){
    bt.clear();
    pm.clear();

    bt.resize((m + 1), std::vector<unsigned char>((diff + 2 * k + 1), '\0'));
    pm.resize(3, DoubleVec((diff + 2 * k + 1), MY_INF_MIN));

    pm[0][k] = TermGapOpPenalty;
    bt[0][k] = (char)16;

    for (int j = 1; j < (diff + k + 1); j++) {
        pm[0][j + k] = TermGapOpPenalty + TermGapExPenalty * j;
        pm[1][j + k] = TermGapOpPenalty + TermGapExPenalty * j;
        pm[2][j + k] = TermGapOpPenalty + TermGapExPenalty * j;
        bt[0][j + k] = (char)8;
    }
    for (int i = 1; i < (k + 1); i++)
        bt[i][k - i] = (char)3;

}

void FNW::InitTwo(int ii, int k, int diff, DoubleMatRef pm2) {
    pm2.clear();
    pm2.resize(3, DoubleVec((diff + 2 * k + 1), MY_INF_MIN));
    if (ii < k + 1) {
        pm2[0][k - ii] = TermGapOpPenalty + TermGapExPenalty * (ii - 1);
        pm2[1][k - ii] = TermGapOpPenalty + TermGapExPenalty * (ii - 1);
        pm2[2][k - ii] = TermGapOpPenalty + TermGapExPenalty * (ii - 1);
    }
}

int FNW::ChooseWay(int p0, int p1, int p2, bool state) {
    if (p0 >= p1)
    {
        if (p0 >= p2)
            return state ? 16 : 0;
        else
            return state ? 48 : 2;
    }
    else if (p1 >= p2)
        return state ? 32 : 1;
    else
        return state ? 48 : 2;
}

inline int FNW::parse(int b, int s) {
    b = (b >> (4 - s * 2));
    return (b & 3) - 1;
}
// }




void FNW::GetMotif(Profile& V_Seq, Profile& H_Seq, int SupporterID) {

    //1. Calculate the Distances Matrix //{
    CharMat      bt;
    DoubleMat    pm;
    DoubleMat    pm2;
    int dir;
    dir = GetDistanceMatrix(V_Seq, H_Seq, bt, pm, pm2, SubMatMotif);
    //}

    //2. (Optional for Debugging) Print Matrices //{
    //cout<<GetAllDoubleMatsToStr(XMat, YMat, MMat)<<endl;
    //cout<<GetAllDirMatsToStr(XtrMat, YtrMat, MtrMat)<<endl;
    //}

    //3. Traverse Backwards //{
    int long_len = V_Seq.GetPatternsSize();
    int short_len = H_Seq.GetPatternsSize();

    int k = KOfKband;
    int i = short_len;
    int b_j = long_len;
    int j = long_len - short_len + k;

    map<int, Profile::Pattern> Insertions;

    while((i > 0)||(j > k)) {
        if ((i == 0) && (dir != 1)) {
            dir = 1;
        }
        if ((b_j == 0) && (dir != 2)) {
            dir = 2;
        }
        switch (dir) {
            case 0: {
                if ((i >= 1) && (b_j >= 1)) {
                    for ( auto H_it = H_Seq.GetPatternAt(i-1).begin(); H_it != H_Seq.GetPatternAt(i-1).end(); ++H_it ) {
                        if ( V_Seq.PatternContainsSite(b_j-1, (*H_it).first ) ) {
                            if ( !V_Seq.IsSiteSupportedBy(b_j-1, (*H_it).first , SupporterID) ) {
                                V_Seq.IncreaseWeight(b_j-1, (*H_it).first, 1, SupporterID);
                            }
                        }else {
                            V_Seq.InsertSiteAt(b_j-1, (*H_it).first, 1, SupporterID);
                        }
                    }
                }
                dir = parse(bt[i][j], 0);
                i -= 1;
                b_j -= 1;
                break;
            }
            case 1: {
                dir = parse(bt[i][j], 1);
                b_j -= 1;
                j -= 1;
                break;
            }
            case 2: {
                Profile::Pattern H_Pat;
                Profile::CreatePatternCopy(H_Seq.GetPatternAt(i-1), H_Pat);
                Insertions[b_j] = H_Pat;
                dir = parse(bt[i][j], 2);
                i -= 1;
                j += 1;
                break;
            }
            default: {
                cout<<"Unknown Direction!"<<endl;
                exit(1);
            }
        }
    }

    //Now Insert all New Positions to the V_Seq Profile
    for (auto rit=Insertions.rbegin(); rit!=Insertions.rend(); ++rit) {
        V_Seq.InsertPatternAt((*rit).first, (*rit).second);
    }
    //}

}

bool FNW::GetAlignment(Profile& V_Seq, Profile& H_Seq, string& SeqAlignment, string& Alphabet) {

    //1. Calculate the Distances Matrix //{
    CharMat      bt;
    DoubleMat    pm;
    DoubleMat    pm2;
    int dir;
    dir = GetDistanceMatrix(V_Seq, H_Seq, bt, pm, pm2, SubMatMotif);
    //}

    //2. (Optional for Debugging) Print Matrices //{
    //cout<<GetAllDoubleMatsToStr(XMat, YMat, MMat)<<endl;
    //cout<<GetAllDirMatsToStr(XtrMat, YtrMat, MtrMat)<<endl;
    //}

    //3. Traverse Backwards //{
    int long_len = V_Seq.GetPatternsSize();
    int short_len = H_Seq.GetPatternsSize();

    int k = KOfKband;
    int i = short_len;
    int b_j = long_len;
    int j = long_len - short_len + k;

    stringstream ss;
    ss.str();

    while((i > 0)||(j > k)) {
        if ((i == 0) && (dir != 1)) {
            dir = 1;
        }
        if ((b_j == 0) && (dir != 2)) {
            dir = 2;
        }
        switch (dir) {
            case 0: { //GO DIAGONALLY
                if ((i >= 1) && (b_j >= 1)) {
                    auto H_it = H_Seq.GetPatternAt(i-1).begin();
                    ss<<Alphabet.at((*H_it).first);
                }
                dir = parse(bt[i][j], 0);
                i -= 1;
                b_j -= 1;
                break;
            }
            case 1: { //GO LEFT
                ss<<"-";
                dir = parse(bt[i][j], 1);
                b_j -= 1;
                j -= 1;
                break;
            }
            case 2: { //GO UP
                dir = parse(bt[i][j], 2);
                return true;
                i -= 1;
                j += 1;
                break;
            }
            default: {
                cout<<"Unknown Direction!"<<endl;
                exit(1);
            }
        }
    }
    //}

    //4. Assign the "Reversed" Alignment to the given argument and return false (no Fine tuning is required) //{
    string TempAlignStr = ss.str();
    std::reverse(TempAlignStr.begin(), TempAlignStr.end());
    SeqAlignment = std::move(TempAlignStr);
    return false;
    //}

}

int FNW::GetDistanceMatrix(Profile& V_Seq, Profile& H_Seq, CharMatRef bt, DoubleMatRef pm, DoubleMatRef pm2, DoubleMatRef SubMat) {

    int longer_len, shorter_len, diff, k = 1;

    int profile_len = V_Seq.GetPatternsSize();
    int sequence_len = H_Seq.GetPatternsSize();

    if (profile_len >= sequence_len) {
        longer_len = profile_len;
        shorter_len = sequence_len;
    } else {
        longer_len = sequence_len;
        shorter_len = profile_len;
    }

    diff = longer_len - shorter_len;
    double new_score, current_score = MY_INF_MIN;

//    cout << "long length:" << longer_len << "\n";
//    cout << "short length:" << shorter_len << "\n";

    while (k <= shorter_len) {
        Init(shorter_len, k, diff, bt, pm);
        for (int i = 1; i < (shorter_len + 1); i++) {
            InitTwo(i, k, diff, pm2);
            for (int z = -k; z < (diff + k + 1); z++) {
                int j = z;
                if ((1 <= (j + i)) && ((j + i) <= longer_len)) {
                    j += k;
                    int bt1 = 0;
                    int bt2 = 0;
                    int bt3 = 0;

                    double TotSubScore = 0.0f;
                    Profile::PatternRef H_Pat = H_Seq.GetPatternAt(i - 1);
                    Profile::PatternRef V_Pat = V_Seq.GetPatternAt(j + i - k - 1);

                    for ( auto H_it = H_Pat.begin(); H_it != H_Pat.end(); ++H_it ) {
                        for ( auto V_it = V_Pat.begin(); V_it != V_Pat.end(); ++V_it ) {
                            TotSubScore += ((*H_it).second.Weight * (*V_it).second.Weight) * SubMat[(*H_it).first][(*V_it).first];
                        }
                    }
                    double Pattern_Uniformity = V_Seq.GetMaxIdentityAt(j + i - k - 1);
                    // t : A[i] ~ B[j]
                    bt1 = ChooseWay(pm[0][j], pm[1][j], pm[2][j], true);
                    pm2[0][j] = maxi(pm[0][j], pm[1][j], pm[2][j]) + TotSubScore;

                    if (InsideStrip(i, j+i-k-1, k, diff)) {
                        double val1 = pm2[0][j-1] + VGapOpPenalty - (TotSubScore * Pattern_Uniformity);
                        double val2 = pm2[1][j-1] + VGapExPenalty;
                        pm2[1][j] = maxi(val1, val2);
                        if (val1 >= val2) {
                            bt2 = 4; // GO DIAGONALLY
                        } else {
                            bt2 = 8; // GO LEFT
                        }
                    }

                    if (InsideStrip(i-1, j+i-k, k, diff)) {
                        double val1 = pm[0][j+1] + HGapOpPenalty - (TotSubScore * Pattern_Uniformity);
                        double val2 = pm[2][j+1] + HGapExPenalty;
                        pm2[2][j] = maxi(val1, val2);
                        if (val1 >= val2) {
                            bt3 = 1;  // GO DIAGONALLY
                        } else {
                            bt3 = 3;  // GO UP
                        }
                    }
                    bt[i][j] = (char)(bt1 + bt2 + bt3);
                    // cout << (int)bt[i][j] << "\n";
                }
            }
            pm = pm2;
        }
        new_score = maxi(pm[0][diff + k], pm[1][diff + k], pm[2][diff + k]);
        // cout << new_score <<"\n";
        if (current_score == new_score || (k * 2) > shorter_len) {
            break;
        } else {
            current_score = new_score;
            k *= 2;
        }
    }
    KOfKband = k;
    return ChooseWay(pm[0][diff + k], pm[1][diff + k], pm[2][diff + k], false);

}
// }

string FNW::GetDoubleMatToStr(DoubleMatRef Mat) {
    stringstream ss;
    ss.str("");

    ss<<"\n\tVALUES MATRIX\n\t=============\n";
    ss<<endl;
    ss.setf(std::ios::fixed);
    ss.precision(2);
    for (unsigned int i=0; i<Mat.size(); i++) {
        for (unsigned int j=0; j<Mat[i].size(); j++) {
            ss<<left<<setw(PRINTWIDTH)<<Mat[i][j];
        }
        ss<<endl;
    }

    return ss.str();
}

string FNW::GetDirMatToStr(DirMatRef Mat) {
    stringstream ss;
    ss.str("");

    ss<<"\n\tDIRECTIONS MATRIX\n\t=================\n";
    ss<<endl;
    ss.setf(std::ios::fixed);
    ss.precision(2);
    for (unsigned int i=0; i<Mat.size(); i++) {
        for (unsigned int j=0; j<Mat[i].size(); j++) {
            ss<<left<<setw(PRINTWIDTH)<<Mat[i][j];
        }
        ss<<endl;
    }

    return ss.str();
}

string FNW::GetAllDoubleMatsToStr(DoubleMatRef XMat, DoubleMatRef YMat, DoubleMatRef MMat) {
    stringstream ss;
    ss.str("");

    ss<<"\n\tALL VALUE MATRICES\n\t==================\n";
    ss<<endl;
    ss.setf(std::ios::fixed);
    ss.precision(2);

    for (unsigned int i=0; i<MMat.size(); i++) {
        for (unsigned int j=0; j<MMat[i].size(); j++) {
            ss<<left<<setw(PRINTWIDTH)<<MMat[i][j]<<left<<setw(PRINTWIDTH+2)<<XMat[i][j];
        }
        ss<<endl;
        for (unsigned int j=0; j<MMat[i].size(); j++) {
            ss<<left<<setw(PRINTWIDTH)<<YMat[i][j]<<setw(PRINTWIDTH+2)<<"";
        }
        ss<<endl<<endl;
    }
    return ss.str();
}

string FNW::GetAllDirMatsToStr(DirMatRef XtrMat, DirMatRef YtrMat, DirMatRef MtrMat) {
    stringstream ss;
    ss.str("");

    ss<<"\n\tALL DIRECTION MATRICES\n\t======================\n";
    ss<<endl;
    ss.setf(std::ios::fixed);
    ss.precision(2);

    for (unsigned int i=0; i<MtrMat.size(); i++) {
        for (unsigned int j=0; j<MtrMat[i].size(); j++) {
            ss<<left<<setw(PRINTWIDTH)<<MtrMat[i][j]<<left<<setw(PRINTWIDTH+2)<<XtrMat[i][j];
        }
        ss<<endl;
        for (unsigned int j=0; j<MtrMat[i].size(); j++) {
            ss<<left<<setw(PRINTWIDTH)<<YtrMat[i][j]<<setw(PRINTWIDTH+2)<<"";
        }
        ss<<endl<<endl;
    }

    return ss.str();
}

string FNW::GetParamsToStr(string& Alphabet) {
    stringstream ss; ss.str("");
    ss << endl;
    ss << "\t\tPARAMETERS\n\t\t==========\n\n";
    ss << "Residue Weight Matrix (Motifs):\n";

    for (unsigned int i=0; i<SubMatMotif.size(); i++) {
        if (i==0) {
            ss << setw(4) << " ";
        }

        ss << setw(4) <<Alphabet.at(i);
    }
    ss<<endl;

    for (unsigned int i=0; i<SubMatMotif.size(); i++) {
        for (unsigned int j=0; j<SubMatMotif[i].size(); j++) {
            if (j==0) {
                ss<<setw(4)<<Alphabet.at(i);
            }
            ss<<setw(4)<<SubMatMotif[i][j];
        }
        ss<<endl;
    }


    ss << "\nResidue Weight Matrix (Aligning):\n";

    for (unsigned int i=0; i<SubMatAlign.size(); i++) {
        if (i==0) {
            ss << setw(4) << " ";
        }
        ss << setw(4) <<Alphabet.at(i);
    }
    ss<<endl;

    for (unsigned int i=0; i<SubMatAlign.size(); i++) {
        for (unsigned int j=0; j<SubMatAlign[i].size(); j++) {
            if (j==0) {
                ss<<setw(4)<<Alphabet.at(i);
            }
            ss<<setw(4)<<SubMatAlign[i][j];
        }
        ss<<endl;
    }

    ss<<endl<<"  H. Gap Opening Penalty : "<< HGapOpPenalty<<endl;
    ss<<" H.Gap Extension Penalty : "<< HGapExPenalty<<endl;
    ss<<" V.Gap Opening   Penalty : "<< VGapOpPenalty<<endl;
    ss<<" V.Gap Extension Penalty : "<< VGapExPenalty<<endl;
    ss<<"Terminal Gap Op. Penalty : "<< TermGapOpPenalty<<endl;
    ss<<"Terminal Gap Ex. Penalty : "<< TermGapExPenalty<<endl;

    return ss.str();
}

