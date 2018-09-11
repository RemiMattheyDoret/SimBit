/*

 Author: Remi Matthey-Doret

    MIT License

    Copyright (c) 2017 Remi Matthey-Doret

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

 */

T1_locusDescription::T1_locusDescription(int byteIndex, int bitIndex, int Locus)
: byte_index(byteIndex), bit_index(bitIndex), locus(Locus)
{
    assert(byte_index * EIGHT + bit_index == locus);
}

T1_locusDescription::T1_locusDescription(int Locus) // minimalist class used only for creating object to search in 'isLocusInT1_outputs_sequenceList'
: byte_index(0), bit_index(0), locus(Locus)
{}

T1_locusDescription::T1_locusDescription(){}

bool T1_locusDescription::operator==(const T1_locusDescription& other) const
{
    // useless to compare byte_index or bit_index
    return this->locus == other.locus;
}

bool T1_locusDescription::operator>(const T1_locusDescription& other) const
{
    // useless to compare byte_index or bit_index
    return this->locus > other.locus;
}

bool T1_locusDescription::operator<(const T1_locusDescription& other) const
{
    // useless to compare byte_index or bit_index
    return this->locus < other.locus;
}

bool T1_locusDescription::operator>(int l) const
{
    // useless to compare byte_index or bit_index
    return this->locus > l;
}

bool T1_locusDescription::operator<(int l) const
{
    // useless to compare byte_index or bit_index
    return this->locus < l;
}

bool T1_locusDescription::operator==(std::pair<int,int>& pa) const
{
    // useless to compare byte_index or bit_index
    return this->byte_index == pa.first && this->bit_index == pa.second;
}

std::string T1_locusDescription::toString()
{
    return std::to_string(this->locus);
}
