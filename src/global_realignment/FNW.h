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

#ifndef FNW_H
#define FNW_H

#include "Utils.h"
#include "TypeDefinitions.h"
#include "Profile.h"
#include "limits"

#define DIAGONALLY 1
#define UPWARDS 2
#define LEFTWARDS 3

#define PRINTWIDTH 4

using namespace std;
using namespace Utils;
using namespace TypeDefinitions;

class FNW {
public:
    //The Default Constructor
    FNW() { };

    //The Copy Constructor
    FNW(const FNW& other) {
        KOfKband = other.KOfKband;
        HGapOpPenalty = other.HGapOpPenalty;
        HGapExPenalty = other.HGapExPenalty;
        VGapOpPenalty = other.VGapOpPenalty;
        VGapExPenalty = other.VGapExPenalty;
        TermGapOpPenalty = other.TermGapOpPenalty;
        TermGapExPenalty = other.TermGapExPenalty;
        SubMatMotif.resize(other.SubMatMotif.size(), DoubleVec(other.SubMatMotif[0].size(), 0.0));
        SubMatAlign.resize(other.SubMatAlign.size(), DoubleVec(other.SubMatAlign[0].size(), 0.0));
        for (unsigned int i=0;i<other.SubMatMotif.size();i++) {
            for (unsigned int j=0;j<other.SubMatMotif[i].size();j++) {
                SubMatMotif[i][j] = other.SubMatMotif[i][j];
                SubMatAlign[i][j] = other.SubMatAlign[i][j];
            }
        }
    }

    //The Assignment Operator
    FNW & operator= (const FNW& other) {
        if (this != &other) { // protect against invalid self-assignment
            KOfKband = other.KOfKband;
            HGapOpPenalty = other.HGapOpPenalty;
            HGapExPenalty = other.HGapExPenalty;
            VGapOpPenalty = other.VGapOpPenalty;
            VGapExPenalty = other.VGapExPenalty;
            TermGapOpPenalty = other.TermGapOpPenalty;
            TermGapExPenalty = other.TermGapExPenalty;
            SubMatMotif.resize(other.SubMatMotif.size(), DoubleVec(other.SubMatMotif[0].size(), 0.0));
            SubMatAlign.resize(other.SubMatAlign.size(), DoubleVec(other.SubMatAlign[0].size(), 0.0));
            for (unsigned int i=0;i<other.SubMatMotif.size();i++) {
                for (unsigned int j=0;j<other.SubMatMotif[i].size();j++) {
                    SubMatMotif[i][j] = other.SubMatMotif[i][j];
                    SubMatAlign[i][j] = other.SubMatAlign[i][j];
                }
            }
        }
        return *this;
    }

    //The Deconstructor
    virtual ~FNW() {
        SubMatMotif.clear();
        SubMatAlign.clear();
    }


    void    SetParameters(double hgo, double hge, double vgo, double vge, double tgo, double tge, DoubleMatRef smatmot, DoubleMatRef smatalign) { HGapOpPenalty = hgo; HGapExPenalty = hge; VGapOpPenalty = vgo; VGapExPenalty = vge; TermGapOpPenalty = tgo; TermGapExPenalty = tge; Set_SubMatMotif(smatmot); Set_SubMatAlign(smatalign); }
    void    SetParameters(double hgo, double hge, double vgo, double vge, double tgo, double tge)                                               { HGapOpPenalty = hgo; HGapExPenalty = hge; VGapOpPenalty = vgo; VGapExPenalty = vge; TermGapOpPenalty = tgo; TermGapExPenalty = tge; }
    void    GetMotif(Profile& V_Seq, Profile& H_Seq, int SupporterID);
    bool    GetAlignment(Profile& V_Seq, Profile& H_Seq, string& SeqAlignment, string& Alphabet);
    // double  GetDistanceMatrix(Profile& V_Seq, Profile& H_Seq, DoubleMatRef MMat, DoubleMatRef XMat, DoubleMatRef YMat, DirMatRef MtrMat, DirMatRef XtrMat, DirMatRef YtrMat, DoubleMatRef SubMat);
    // updated first time{
    // double  GetDistanceMatrix(Profile& V_Seq, Profile& H_Seq, DoubleMatRef MMat, DoubleMatRef XMat, DoubleMatRef YMat, DoubleMat SMat, DoubleMatRef SubMat);
    // }
    int GetDistanceMatrix(Profile& V_Seq, Profile& H_Seq, CharMatRef bt, DoubleMatRef pm, DoubleMatRef pm2, DoubleMatRef SubMat);
    void    Set_SubMatMotif(DoubleMatRef val)  {
        SubMatMotif.resize(val.size(), vector<double>(val[0].size(), 0.0));
        for (unsigned int i=0;i<val.size();i++) {
            for (unsigned int j=0;j<val[i].size();j++) {
                SubMatMotif[i][j] = val[i][j];
            }
        }
    }
    void    Set_SubMatAlign(DoubleMatRef val)  {
        SubMatAlign.resize(val.size(), vector<double>(val[0].size(), 0.0));
        for (unsigned int i=0;i<val.size();i++) {
            for (unsigned int j=0;j<val[i].size();j++) {
                SubMatAlign[i][j] = val[i][j];
            }
        }
    }

    string  GetDoubleMatToStr(DoubleMatRef Mat);
    string  GetDirMatToStr(DirMatRef Mat);
    string  GetAllDoubleMatsToStr(DoubleMatRef XMat, DoubleMatRef YMat, DoubleMatRef MMat);
    string  GetAllDirMatsToStr(DirMatRef XtrMat, DirMatRef YtrMat, DirMatRef MtrMat);
    string  GetParamsToStr(string& Alphabet);

    // Kband {
    inline double maxi(double a, double b);
    inline double maxi(double a, double b, double c);
    inline int parse(int b, int s);
    inline bool InsideStrip(int i, int j, int k, int diff);
    void Init(int m, int k, int diff, CharMatRef bt, DoubleMatRef pm);
    void InitTwo(int ii, int k, int diff, DoubleMatRef pm2);
    int ChooseWay(int p0, int p1, int p2, bool state);
    // }
protected:

private:
    int          KOfKband;
    double       MY_INF_MIN = -std::numeric_limits<double>::infinity();
    double       HGapOpPenalty;
    double       HGapExPenalty;
    double       VGapOpPenalty;
    double       VGapExPenalty;
    double       TermGapOpPenalty;
    double       TermGapExPenalty;
    DoubleMat    SubMatMotif;
    DoubleMat    SubMatAlign;
};

#endif // FNW_H
