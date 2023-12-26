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

#ifndef ALIGNMENTMETRICS_H
#define ALIGNMENTMETRICS_H

#include "TypeDefinitions.h"
#include "Utils.h"

using namespace TypeDefinitions;
using namespace Utils;


class AlignmentMetrics {
    public:
        AlignmentMetrics(string SeqsFileName); //The Starting Alignment Filename
        virtual ~AlignmentMetrics();

        double Get_APSI()         { return APSI;         }

    protected:

    private:
        double APSI;

        //Assisting Functions for Metrics Calculation
        int min(int a, int b);
        double average(int a, int b);
        double GetAPSI(StringVecRef SequenceIds, StringVecRef Alignments);
};

#endif // ALIGNMENTMETRICS_H
