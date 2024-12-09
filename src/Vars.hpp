#ifndef VARS_HPP
#define VARS_HPP

#include <vector>
#include <memory>
#include <cmath>
#include <sstream>
#include <variant>

template <int N>
class Vars
{
    public:

        Vars() : data() {}
        Vars(const std::array<double, N>& in) : data(in) {}
        Vars(double in) : data() { for (size_t i = 0; i < N; i++) data[i] = in; }

        virtual ~Vars() {}

        inline const double& operator[](int i) const
        {
            return data[i];
        }

        inline double& operator[] (int i)
        {
            return data[i];
        }

        void operator+=(const Vars<N>& v);
        void operator-=(const Vars<N>& v);

        Vars<N> operator-();

        const std::array<double, N>& getData() const
        {
            return data;
        }

        std::array<double, N>& getData()
        {
            return data;
        }
        
    protected:
        std::array<double, N> data;
};

template <int N>
inline void Vars<N>::operator+= (const Vars<N>& v)
{
    for (int i = 0; i < N; i++)
    {
        data[i] += v[i];
    }
}

template <int N>
inline void Vars<N>::operator-= (const Vars<N>& v)
{
    for (int i = 0; i < N; i++)
    {
        data[i] -= v[i];
    }
}

template <int N>
inline Vars<N> Vars<N>::operator-()
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = -data[i];
    }
    
    return out;
}

//////////////Non member operators///////////////////

template <int N>
inline bool operator== (const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        if(u[i] != v[i])
            return false;
    }
    return true;
}

template <int N>
inline Vars<N> operator+ (const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i] + v[i];
    }
    return out;
}

// u - v
template <int N>
inline Vars<N> operator- (const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i] - v[i];
    }
    return out;
}

// w * u
template <int N>
inline Vars<N> operator* (const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i]*v[i];
    }
    return out;
}

// w / u (po slozkach)
template <int N>
inline Vars<N> operator/ (const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i]/v[i];
    }
    return out;
}

// a * u
template <int N>
inline Vars<N> operator* (double a, const Vars<N>& u)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i]*a;
    }
    return out;
}

// u * a
template <int N>
inline Vars<N> operator* (const Vars<N>& u, double a)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i]*a;
    }
    return out;
}

// u / a
template <int N>
inline Vars<N> operator/ (const Vars<N>& u, double a)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = u[i]/a;
    }
    return out;
}


//////////////Non member function///////////////////

template <int N>
inline double norm2(const Vars<N>& u)
{
    double out = 0;
    for (int i = 0; i < N; i++)
    {
        out += u[i]*u[i];
    }
    return std::sqrt(out);
}

template <int N>
inline Vars<N> unit(const Vars<N>& u)
{
    return u/norm2(u);
}

template <int N>
inline Vars<N> sqrt(const Vars<N>& u)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = sqrt(u[i]);
    }
    return out;
}

template <int N>
inline Vars<N> abs(const Vars<N>& u)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = fabs(u[i]);
    }
    return out;
}

template <int N>
inline double dot(const Vars<N>& u, const Vars<N>& v)
{
    double out = 0;
    for (int i = 0; i < N; i++)
    {
        out += u[i]*v[i];
    }
    return out;
}

template <int N>
inline double sum(const Vars<N>& u)
{
    double out;
    for (int i = 0; i < N; i++)
    {
        out += u[i];
    }
    return out;
}

template <int N>
inline Vars<N> max(const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = std::max(u[i], v[i]);
    }
    return out;
}

template <int N>
inline Vars<N> min(const Vars<N>& u, const Vars<N>& v)
{
    Vars<N> out;
    for (int i = 0; i < N; i++)
    {
        out[i] = std::min(u[i], v[i]);
    }
    return out;
}

template <int N>
inline Vars<N> zeroSmallNumbers(const Vars<N>& u)
{
    constexpr double SMALLNUM = 1e-12;

    Vars<N> out = u;
    for (int i = 0; i < N; i++)
    {
        if(std::fabs(out[i]) < SMALLNUM)
        {
            out[i] = 0.0;
        }
    }
    return out;
}

template <int N>
inline std::ostream& operator<<(std::ostream& os, const Vars<N>& u)
{
    for (int i = 0; i < N; i++)
    {
        os << u[i] << " ";
    }

    return os;
}


Vars<3> cross(const Vars<3>& u, const Vars<3>& v);

Vars<3> angleAngleToUnit(double xy, double xz);


#endif // VARS_HPP