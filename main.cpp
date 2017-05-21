/*
    EDSO: Elastic Degenerate Sequence Outputter

    Copyright (C) 2017 Ahmad Retha

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <Variant.h>

#define BUFFERSIZE 1000000

using namespace std;
using namespace vcflib;
typedef vector<string> Segment;

char BUFF[BUFFERSIZE];
int BUFFLIMIT = 0;
int POS = 0;

char getNextChar(ifstream & f)
{
    if (BUFFLIMIT == 0) {
        f.read(BUFF, BUFFERSIZE);
        BUFFLIMIT = f.gcount();
        POS = 0;
        if (BUFFLIMIT == 0) {
            return '\0';
        }
    }
    char c = BUFF[POS];
    if (++POS == BUFFLIMIT) {
        BUFFLIMIT = 0;
    }
    return c;
}

void output(ofstream & f, Segment & segment)
{
    if (segment.size() == 1)
    {
        f << segment[0];
    }
    else
    {
        int j = 1;
        int l = segment.size();
        f << "{";
        for (const string & s : segment) {
            f << s;
            if (j != l) {
                f << ",";
            }
            j++;
        }
        f << "}";
    }
}

int main(int argc, char * argv[])
{
    string help = "EDSO (Elastic Degenerate Sequence Outputter) takes a reference \
    fasta file and VCF file and produces an EDS format file.\n\n \
    \tUsage example: ./edso reference.fasta variants.vcf.gz [outfile.eds]";

    if (argc == 1 || (argc == 2 && (strcmp("--help", argv[1]) == 0 || strcmp("-h", argv[1]) == 0))) {
        cout << help << endl;
        return 0;
    }

    if (!(argc == 3 || argc == 4)) {
        cerr << "Invalid number of arguments!" << endl;
        cout << help << endl;
        return 1;
    }

    //ref file
    string refName = argv[1];
    ifstream rf(refName.c_str(), ios::in);
    if (!rf.good()) {
        cerr << "Error: Failed to open reference file!" << endl;
        return 1;
    }

    //vcf file
    string vcfName = argv[2];
    VariantCallFile vf;
    vf.open(vcfName);
    if (!vf.is_open()) {
        cerr << "Error: Failed to open variants file!" << endl;
        return 1;
    }

    //out file
    string ofName;
    if (argc == 4)
    {
        ofName = argv[3];
    }
    else
    {
        size_t pos = refName.find_last_of("\\/");
        if (pos == string::npos) {pos = 0;} else {pos += 1;}
        ofName = refName.substr(pos) + ".eds";
    }
    ofstream of(ofName.c_str(), ios::out);
    if (!of.good()) {
        cerr << "Error: Failed to open output file!" << endl;
        return 1;
    }

    cout << "EDSO processing..." << endl;

    string tBuff = "";
    tBuff.reserve(BUFFERSIZE);
    char c;
    unsigned int rfIdx = 1, vfIdx = 0, i = 0;
    Segment segment;
    getline(rf, tBuff); //skip first line of fasta file
    tBuff = "";

    //create variables for reading through vcf records and looking for duplicates
    Variant var(vf), vBuffer(vf), var2(vf);
    bool hasMoreVariants = true;
    Segment vAlleles;

    // read first variant and possibly successive duplicates for the same position, removing duplicate alleles
    hasMoreVariants = vf.getNextVariant(var);
    if (hasMoreVariants)
    {
        vfIdx = (unsigned int) var.position;
        for (const auto & a : var.alleles) {
            if (a[0] != '<') {
                vAlleles.push_back(a);
            }
        }
        while (true) {
            hasMoreVariants = vf.getNextVariant(var2);
            if (hasMoreVariants)
            {
                if (var2.position == var.position) {
                    for (const auto & a : var2.alt) {
                        if (a[0] != '<') {
                            vAlleles.push_back(a);
                        }
                    }
                } else {
                    vBuffer = var2;
                    break;
                }
            }
            else
            {
                break;
            }
        }
    }

    // go through the reference sequence
    while ((c = getNextChar(rf)) != '\0')
    {
        if (!(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N')) {
            continue;
        }

        if (rfIdx != vfIdx)
        {
            tBuff += c;
            i++;
            if (i >= BUFFERSIZE) {
                segment.clear();
                segment.push_back(tBuff);
                output(of, segment);
                tBuff = "";
                i = 0;
            }
        }
        else
        {
            segment.clear();
            if (tBuff.length() > 0) {
                segment.push_back(tBuff);
                output(of, segment);
                tBuff = "";
            }

            //then search current variant
            if (vAlleles.size() > 0) {
                output(of, vAlleles);
                vAlleles.clear();
            }

            //fetch the next variant to be searched when its position comes up
            if (vBuffer.alleles.size() > 0)
            {
                vfIdx = (unsigned int) vBuffer.position;
                for (const auto & a : vBuffer.alleles) {
                    if (a[0] != '<') {
                        vAlleles.push_back(a);
                    }
                }

                vBuffer.alleles.clear();

                while (true) {
                    hasMoreVariants = vf.getNextVariant(var2);
                    if (hasMoreVariants)
                    {
                        if (vfIdx == (unsigned int) var2.position) {
                            for (const auto & a : var2.alt) {
                                if (a[0] != '<') {
                                    vAlleles.push_back(a);
                                }
                            }
                        } else {
                            vBuffer = var2;
                            break;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
        rfIdx++;
    }
    if (tBuff.length() > 0)
    {
        segment.clear();
        segment.push_back(tBuff);
        output(of, segment);
        tBuff = "";
    }

    rf.close();
    of.close();

    cout << "Done! Output to " << ofName << endl;

    return 0;
}
