template<typename vectorType, typename iteratorType>
class ZipIterator
{
private:

    iteratorType haploP;
    iteratorType flippedP;
    iteratorType haploEndP;
    iteratorType flippedEndP;

    void skipPotentialValueFoundInBothHaploAndFlipped();
    
public:
    ZipIterator<vectorType, iteratorType>(iteratorType a, iteratorType b, iteratorType c, iteratorType d);
    ZipIterator<vectorType, iteratorType>(const ZipIterator<vectorType, iteratorType>& other);
    ZipIterator<vectorType, iteratorType>();
    ZipIterator<vectorType, iteratorType>& operator=(const ZipIterator<vectorType, iteratorType>&& other);
    ZipIterator<vectorType, iteratorType>& operator=(const ZipIterator<vectorType, iteratorType>& other);

    ZipIterator<vectorType, iteratorType> operator++();
    ZipIterator<vectorType, iteratorType> operator++(int) = delete;
    unsigned int operator*();
    template<typename vType, typename iType>
    friend bool operator==(ZipIterator<vType, iType>& lhs, ZipIterator<vType, iType>& rhs);
    template<typename vType, typename iType>
    friend bool operator!=(ZipIterator<vType, iType>& lhs, ZipIterator<vType, iType>& rhs);
    template<typename vType, typename iType>
    friend bool operator>(ZipIterator<vType, iType>& lhs, ZipIterator<vType, iType>& rhs);
    template<typename vType, typename iType>
    friend bool operator<(ZipIterator<vType, iType>& lhs, ZipIterator<vType, iType>& rhs);
    /*friend bool operator==(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType> rhs);
    friend bool operator!=(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType> rhs);
    friend bool operator>(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType> rhs);
    friend bool operator<(ZipIterator<vectorType, iteratorType>& lhs, ZipIterator<vectorType, iteratorType> rhs);*/
    bool isMore(); // deprecated
    iteratorType get_haploP();
    iteratorType get_flippedP();

    void lower_bound(unsigned int target);
    void upper_bound(unsigned int target);
};
