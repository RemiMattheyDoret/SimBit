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
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRAreadLineNTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

 */

void BinaryFileToRead::testReachedEnd()
{
    if (strm.eof())
    {
        std::cout << "There seems to have not enough data in the binary file '" << path << "', given the input parameters.\n";
        abort();   
    }

    char nothing;
    read(nothing); // sets EOF to true if there was nothing left to read

    if (!strm.eof())
    {
        std::cout << "There seems to have more data in the binary file '" << path << "', than expected from the input parameters.\n";
        int nbExtraBytes = 1;
        while (strm >> nothing)
        {
            nbExtraBytes++;
            if (nbExtraBytes > 99999)
            {
                std::cout << "There were more than " << nbExtraBytes << " extra bytes!\n";
                abort();   
            }
        }
        std::cout << "There were " << nbExtraBytes << " extra bytes!\n";
        abort();
    }
}

void BinaryFileToRead::open(std::string& p)
{
    path=p;
    strm.open(path.c_str(), std::ios::in | std::ios::binary);

    if (!strm.is_open())
    {
        std::cout << "GP->BinaryFileToRead (" << path << ") failed to open! Please note that 'readPopFromBinary' should NOT be indicated relative to 'GeneralPath' but should be a root path.\n";
        abort();
    }
}

void BinaryFileToRead::close()
{
    strm.close();
}

template<typename T> void BinaryFileToRead::read(T& x)
{
    strm.read(reinterpret_cast<char *>(&x), sizeof(T));
    //std::cout << "read: " << x << " of type " << typeid(T).name() << "\n";
}

template<typename T> void BinaryFileToRead::readVector(std::vector<T>& x)
{
    size_t size;
    read(size);
    assert(size >= 0);

    x.resize(size);
    strm.read(reinterpret_cast<char*>(&x[0]), size*sizeof(T)); // Read up until vector is full

    /*
    std::cout <<"read: ";
    for (size_t i = 0 ; i < x.size() ; ++i)
    {
        std::cout << x[i] << " ";
    }
    std::cout << " of type " << typeid(T).name()<<"\n";
    */
}

