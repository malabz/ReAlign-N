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

#ifndef TYPEDEFINITIONS_H
#define TYPEDEFINITIONS_H

#include <iostream>
#include <vector>
#include <map>

using namespace std;

namespace TypeDefinitions {
    const string DNA_Order = "ACGT*";
    const string RNA_Order = "ACGU*";
    const string  AA_Order = "ARNDCQEGHILKMFPSTWYVBZX*";

    typedef vector<string>   StringVec;
    typedef vector<string>&  StringVecRef;
    typedef vector<string>*  StringVecPtr;

    typedef const vector<string>   StringCVec;
    typedef const vector<string>&  StringCVecRef;
    typedef const vector<string>*  StringCVecPtr;

    typedef vector<bool>   BoolVec;
    typedef vector<bool>&  BoolVecRef;
    typedef vector<bool>*  BoolVecPtr;

    typedef vector<double>   DoubleVec;
    typedef vector<double>&  DoubleVecRef;
    typedef vector<double>*  DoubleVecPtr;

    typedef const vector<double>   DoubleCVec;
    typedef const vector<double>&  DoubleCVecRef;
    typedef const vector<double>*  DoubleCVecPtr;

    typedef vector<int>      IntVec;
    typedef vector<int>&     IntVecRef;
    typedef vector<int>*     IntVecPtr;

    typedef vector<size_t>      DirVec;
    typedef vector<size_t>&     DirVecRef;
    typedef vector<size_t>*     DirVecPtr;

    typedef const vector<int>      IntCVec;
    typedef const vector<int>&     IntCVecRef;
    typedef const vector<int>*     IntCVecPtr;

    typedef vector <vector< vector<int > > >   IntHyp;
    typedef vector <vector< vector<int > > >&  IntHypRef;
    typedef vector <vector< vector<int > > >*  IntHypPtr;

    typedef const vector <vector< vector<int > > >   IntCHyp;
    typedef const vector <vector< vector<int > > >&  IntCHypRef;
    typedef const vector <vector< vector<int > > >*  IntCHypPtr;

    typedef vector< vector<string> >  StringMat;
    typedef vector< vector<string> >& StringMatRef;
    typedef vector< vector<string> >* StringMatPtr;

    typedef const vector< vector<string> >  StringCMat;
    typedef const vector< vector<string> >& StringCMatRef;
    typedef const vector< vector<string> >* StringCMatPtr;

    typedef vector< vector<double> >  DoubleMat;
    typedef vector< vector<double> >& DoubleMatRef;
    typedef vector< vector<double> >** DoubleMatPtr;

    typedef const vector< vector<double> >  DoubleCMat;
    typedef const vector< vector<double> >& DoubleCMatRef;
    typedef const vector< vector<double> >* DoubleCMatPtr;

    typedef vector <vector< vector<double > > >   DoubleHyp;
    typedef vector <vector< vector<double > > >&  DoubleHypRef;
    typedef vector <vector< vector<double > > >*  DoubleHypPtr;

    typedef const vector <vector< vector<double > > >   DoubleCHyp;
    typedef const vector <vector< vector<double > > >&  DoubleCHypRef;
    typedef const vector <vector< vector<double > > >*  DoubleCHypPtr;

    typedef vector< vector<int> >     IntMat;
    typedef vector< vector<int> >&    IntMatRef;
    typedef vector< vector<int> >*    IntMatPtr;

    typedef vector< vector<size_t> >     DirMat;
    typedef vector< vector<size_t> >&    DirMatRef;
    typedef vector< vector<size_t> >*    DirMatPtr;

    typedef const vector< vector<int> >     IntCMat;
    typedef const vector< vector<int> >&    IntCMatRef;
    typedef const vector< vector<int> >*    IntCMatPtr;

    typedef map<int, int>  IntMap;
    typedef map<int, int>& IntMapRef;
    typedef map<int, int>* IntMapPtr;

    typedef const map<int, int>  IntCMap;
    typedef const map<int, int>& IntCMapRef;
    typedef const map<int, int>* IntCMapPtr;

    typedef vector< map<int, int> > IntMapVec;
    typedef vector< map<int, int> >& IntMapVecRef;
    typedef vector< map<int, int> >* IntMapVecPtr;

    typedef const vector< map<int, int> > IntMapCVec;
    typedef const vector< map<int, int> >& IntMapCVecRef;
    typedef const vector< map<int, int> >* IntCMapCVecPtr;

    typedef long long int LongInt;

    typedef vector< vector<unsigned char> >  CharMat;
    typedef vector< vector<unsigned char> >& CharMatRef;
    typedef vector< vector<unsigned char> >* CharMatPtr;
};

#endif // TYPEDEFINITIONS_H
