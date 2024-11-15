#ifndef THERMOVAR_HPP
#define THERMOVAR_HPP

#include <vector>
#include <memory>

#include "Vars.hpp"

class ThermoVar : public Vars<3>
{
    public:
        enum {T, P, A};

        ThermoVar() : Vars<3>() {}
        ThermoVar(const Vars<3>& in) : Vars<3>(in) {}
        ThermoVar(const std::array<double, 3>& in) : Vars<3>(in) {}

        virtual ~ThermoVar() {}

        //void operator+=(const ThermoVar& v);
        //void operator-=(const ThermoVar& v);


        double temperature() const;
        double pressure() const;        
        double soundSpeed() const;

    private:
};



//////////////Non member operators///////////////////

// u == v
//bool operator== (const ThermoVar& u, const ThermoVar& v);

// u + v
//ThermoVar operator+ (const ThermoVar& u, const ThermoVar& v);

// u - v
//ThermoVar operator- (const ThermoVar& u, const ThermoVar& v);

// u * v
//ThermoVar operator* (const ThermoVar& u, const ThermoVar& v);

// u / v
//ThermoVar operator/ (const ThermoVar& u, const ThermoVar& v);

// a * u
//ThermoVar operator* (double a, const ThermoVar& u);

// u * a
//ThermoVar operator* (const ThermoVar& u, double a);

// u / a
//ThermoVar operator/ (const ThermoVar& u, double a);

// ThermoVar, Vars

// u + v
//ThermoVar operator+ (const ThermoVar& u, const Vars<5>& v);

// u - v
//ThermoVar operator- (const ThermoVar& u, const Vars<5>& v);

// w * u
//ThermoVar operator* (const ThermoVar& u, const Vars<5>& v);



//////////////Non member function///////////////////

//ThermoVar sqrt(const ThermoVar& u);

#endif // THERMOVAR_HPP