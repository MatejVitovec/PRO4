#ifndef LIMITER_HPP
#define LIMITER_HPP

#include "../Mesh/Mesh.hpp"
#include "../Field.hpp"
#include "../VolField.hpp"
#include "../Compressible.hpp"
#include "../Condensation/CompressibleMixture.hpp"
#include "../Mat.hpp"

class Limiter
{
    public:

        Limiter() {}

        virtual ~Limiter() {}

        virtual Field<Vars<5>> calculateLimiter(const VolField<Compressible>& w, Field<Mat<5,3>>& grad, const Mesh& mesh) const;
        virtual Field<Vars<9>> calculateLimiter(const VolField<CompressibleMixture>& w, Field<Mat<9,3>>& grad, const Mesh& mesh) const;

    protected:
        virtual double limiterFunction(double y) const;

};

#endif // LIMITER_HPP