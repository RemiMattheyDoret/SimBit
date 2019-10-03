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



class T1_locusDescription
{
public:
    int byte_index;
    int bit_index;
    int locus; // locus = byte_index * EIGHT + bit_index

    T1_locusDescription(int byteIndex, int bitIndex, int Locus);
    T1_locusDescription(int Locus);
    T1_locusDescription();

    bool operator==(const T1_locusDescription& other) const;
    bool operator>(const T1_locusDescription& other) const;
    bool operator<(const T1_locusDescription& other) const;
    bool operator==(int l) const;
    bool operator>(int l) const;
    bool operator<(int l) const;
    bool operator==(std::pair<int,int>& pa) const;
    std::string toString();
};

