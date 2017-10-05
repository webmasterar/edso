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

using namespace std;
using namespace vcflib;

typedef vector<string> Segment;

struct VarItem
{
    unsigned int pos;
    unsigned int skip;
    Segment seg;
};

typedef vector<struct VarItem> VarItemArray;

#define BUFFERSIZE 1000000

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

    if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') {
        return c;
    } else {
        return getNextChar(f);
    }
}

bool populateVarItemArray(string vcfName, VarItemArray & variantItems)
{
    VariantCallFile vf;
    vf.open(vcfName);
    if (!vf.is_open()) {
        cerr << "Error: Failed to open variants file!" << endl;
        return false;
    }
    Variant v(vf);
    unsigned int prevPos = 0, prevRefLen = 0;
    unsigned int currPos, currRefLen;

    while (vf.getNextVariant(v))
    {
        currPos = v.position;
        currRefLen = v.ref.length();
        //handle duplicates
        if (currPos == prevPos)
        {
            //check duplicate doesn't have a longer ref -- if it does merge new ref into segment
            if (currRefLen > prevRefLen && !(v.ref[0] == '<' || v.ref[0] == '.'))
            {
                Segment updSeg;
                for (const string & t : variantItems.back().seg)
                {
                    string r = t;
                    if (r.length() < currRefLen) {
                        r += v.ref.substr(prevRefLen);
                    }
                    updSeg.push_back(r);
                }
                updSeg[0] = v.ref;
                variantItems.back().seg = updSeg;
                variantItems.back().skip = currRefLen;
            }
            //add new alts
            for (const string & t : v.alt)
            {
                if (!(v.ref[0] == '<' || v.ref[0] == '.')) {
                    variantItems.back().seg.push_back(t);
                }
            }
        }
        //handle nested variants
        else if (currPos < (prevPos + prevRefLen))
        {
            Segment lastSeg;
            for (const string & s : variantItems.back().seg)
            {
                for (const string & t : v.alt)
                {
                    if ((prevPos + s.length()) > currPos && !(t[0] == '<' || t[0] == '.')) {
                        string r;
                        r = s.substr(0, currPos - prevPos);
                        r += t;
                        r += s.substr(currPos - prevPos + 1);
                        lastSeg.push_back(r);
                    }
                }
            }
            for (const string & r : lastSeg) {
                variantItems.back().seg.push_back(r);
            }
        }
        //handle regular variants
        else
        {
            Segment currSeg;
            for (const string & t : v.alleles) {
                if (!(t[0] == '<' || t[0] == '.')) {
                    currSeg.push_back(t);
                }
            }
            struct VarItem varItem = {currPos, currRefLen, currSeg};
            variantItems.push_back(varItem);
            prevPos = currPos;
            prevRefLen = currRefLen;
        }
    }
    return true;
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
    unsigned int rfIdx = 1, i = 0, j;
    Segment segment;
    getline(rf, tBuff); //skip first line of fasta file
    tBuff = "";

    //create variables for reading through vcf records and looking for duplicates
    VarItemArray variantItems;
    bool parsedVCF = populateVarItemArray(string(argv[2]), variantItems);
    if (!parsedVCF) {
        cerr << "Exiting!" << endl;
        return 1;
    }
    VarItemArray::iterator vit = variantItems.begin();
    unsigned int vPos = vit->pos;
    unsigned int vSkip = vit->skip;
    Segment vSeg = vit->seg;

    // go through the reference sequence
    while ((c = getNextChar(rf)) != '\0')
    {
        if (rfIdx != vPos)
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
            rfIdx++;
        }
        else //rfIdx == vPos
        {
            //output buffered text
            segment.clear();
            if (tBuff.length() > 0) {
                segment.push_back(tBuff);
                output(of, segment);
                tBuff = "";
            }

            //then output current variant segment and skip required number of characters
            output(of, vSeg);
            for (j = 1; j < vSkip; j++) {
                getNextChar(rf);
            }
            rfIdx += vSkip;

            //get the next variant position, skipping repeats or nested variants
            do {
                if (++vit == variantItems.end()) {
                    break;
                }
                vPos = vit->pos;
                vSkip = vit->skip;
                vSeg = vit->seg;
            }
            while (vPos < rfIdx);
        }
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
