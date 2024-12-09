#ifndef EXPLICITEULERCONDENSATIONMIXTURE_HPP
#define EXPLICITEULERCONDENSATIONMIXTURE_HPP

#include "FVMSchemeCondensationMixture.hpp"

class ExplicitEulerCondensationMixture : public FVMSchemeCondensationMixture
{
    public:

        ExplicitEulerCondensationMixture(Mesh&& mesh_, std::unique_ptr<FluxSolver> fluxSolver_, std::unique_ptr<Thermo> thermo_) : FVMSchemeCondensationMixture(std::move(mesh_), std::move(fluxSolver_), std::move(thermo_)) {}

        virtual ~ExplicitEulerCondensationMixture() {}

        void solve();        
        

    private:

};



#endif // EXPLICITEULERCONDENSATIONMIXTURE_HPP