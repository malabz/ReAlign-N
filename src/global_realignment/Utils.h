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


#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

/** Header file Inclusions **/
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <map>
#include <cstdlib>
#include <ctime>
#include <assert.h>
#include <unistd.h>
#include <algorithm>
#include <functional>
#include <regex>
#include <unordered_map>
#include <map>
#include <cfloat>
#include <omp.h>    //For parallelizing the program

using namespace std;


namespace Utils {

/*******************************/
/** File Management Functions **/
/*******************************/ //{
//Function to Save a string (with contents) to a File
inline bool WriteToFile(string filename, string contents, bool append, bool clearContents) {
    ofstream myfile;

    if (clearContents) {
        myfile.open(filename.c_str(),ios::out|ios::trunc);
        myfile.close();
    }

    if (append) {
        myfile.open(filename.c_str(), fstream::out | fstream::app);
    } else {
        myfile.open(filename.c_str(), fstream::out);
    }

    if (myfile.is_open()) {
        myfile << contents;
        myfile.close();

        return true;
    } else cerr << "Unable to open file: "<<filename<<" for writing.";

    return false;
}

//Function to Read Contents from File and get the returned result as a string
inline string ReadFromFile(string filename, bool AddNewLine) {
    string contents = "";
    string line;
    ifstream myfile (filename.c_str());
    if (myfile.is_open()) {
        while ( myfile.good() ) {
            getline (myfile,line);
            contents+=line;
            if (AddNewLine) {
                contents+="\n";
            }
        }
        myfile.close();
    } else cerr << "Unable to open file: "<<filename<<" for reading.";

    return contents;
}

//}

/************************************/
/** Runtime Intervention Functions **/
/************************************/ //{

//Function that Pauses the Command-line output Until a key is pressed by the user
inline void Pause() {
    fflush(stdin);
    cout<<endl<<"Press enter to continue..."<<endl;
    getchar();
}

//Function to flush stdin
inline void ClearInputStream() {
    rewind(stdin);
}

//}

/***********************************/
/** String Manipulation Functions **/
/***********************************/ //{

//Function that replaces every appearance of part of a string with another string
inline void ReplaceAll(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
}

//Function that Returns whether a string starts with the given preffix of not
inline bool StartsWith (std::string const &fullString, std::string const &prefix) {
    return prefix.length() <= fullString.length() && equal(prefix.begin(), prefix.end(), fullString.begin());
}

//Trim a string from the start
inline std::string &TrimLeadingSpaces(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

//Trim a string from the end
inline std::string &TrimTrailingSpace(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

//Trim from both ends
inline std::string &Trim(std::string &s) {
    return TrimLeadingSpaces(TrimTrailingSpace(s));
}

//Convert To UpperCase
inline std::string ConvertToUpperCase(std::string s) {
    string upS(s);
    std::transform(upS.begin(), upS.end(),upS.begin(), ::toupper);
    return upS;
}

//Convert To LowerCase
inline std::string ConvertToLowerCase(std::string s) {
    string upS(s);
    std::transform(upS.begin(), upS.end(),upS.begin(), ::tolower);
    return upS;
}

//Function that Returns the first position (index) where a string (target) is found in another string (source) or returns -1 if not found
inline int FirstIndexOf(string source, string target) {
    string::size_type location = source.find( target, 0 );
    if( location != string::npos ) {
        return location;
    } else {
        return -1;
    }
}

//Overloaded version of the previous Function. Acts as previously but starting after the given starting position (stratPos)
inline int FirstIndexOf(string source, string target, int startPos) {
    string::size_type location = source.find( target, startPos );
    if( location != string::npos ) {
        return location;
    } else {
        return -1;
    }
}

//Function that returns the reversed version of a given string
inline string GetReversedString(string source) {
    return string ( source.rbegin(), source.rend() );
}

//Function that Returns the last position (index) where a string (target) is found in another string (source) or returns -1 if not found
inline int LastIndexOf(string source, string target) {
    int pos = FirstIndexOf(GetReversedString(source), target);
    if (pos==-1) return -1;
    else return source.size() -1 - pos;
}

//Function that Returns whether a given string contains another string or not
inline bool Contains(string source, string target) {
    size_t found = source.find(target);
    if (found!=string::npos) return true;
    else                     return false;
}

//Split a string based on a given delimiter and return a vector of strings
inline void SplitStr(std::string const& input, char delimiter, vector<std::string>& results) {
    stringstream str(input);
    string s;
    results.clear();
    while(getline(str, s, ',')) {
        results.push_back(s);
    }
}

//Function that Splits a Given String According to a delimiter and returns an already constructed container also given as argument
inline void SplitValues(const std::string &s, char delim, std::vector<double> &elems) {
    std::stringstream ss(s);
    std::string item;
    int i=0;
    while(std::getline(ss, item, delim)) {
        elems[i++] = (atof(item.c_str()));
    }
}

//}

}



#endif // UTILS_H_INCLUDED
