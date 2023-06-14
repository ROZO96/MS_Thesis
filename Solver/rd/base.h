#ifndef RD_BASE_H
#define RD_BASE_H

//
//double max_val(double, double);
//
//double min_val(double, double);

template <typename T>
inline T const& max_val (T const& a, T const& b)
{
    return a > b ? a:b;
}

template <typename T>
inline T const& min_val (T const& a, T const& b)
{
    return a < b ? a:b;
}

#endif //RD_BASE_H
