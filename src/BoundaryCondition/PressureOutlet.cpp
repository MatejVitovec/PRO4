#include "PressureOutlet.hpp"


Compressible PressureOutlet::calculateState(const Compressible& w, const ThermoVar& thermoVar, const Face& f, const Thermo * const thermoModel) const
{
    if(w.normalVelocity(f.normalVector)/thermoVar.soundSpeed() >= 1.0)
    //if(w.absVelocity()/w.soundSpeed() >= 1.0)
    {
        return w;
    }

    return thermoModel->primitiveToConservative(Vars<5>({w.density(),
                                                         w.velocityU(),
                                                         w.velocityV(),
                                                         w.velocityW(),
                                                         pressure}));
}