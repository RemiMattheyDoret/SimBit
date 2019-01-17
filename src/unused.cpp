

template <typename T, typename _ = void>
struct is_vector { 
    static const bool value = false;
};

template <typename T>
struct is_vector< T,
                  typename std::enable_if<
                      std::is_same<T,
                              std::vector< typename T::value_type,
                                           typename T::allocator_type >
                             >::value
                  >::type
                >
{
    static const bool value = true;
};




