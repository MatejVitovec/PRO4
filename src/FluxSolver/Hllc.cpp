#include "Hllc.hpp"

Vars<5> Hllc::claculateFlux(const Compressible& wl, const Compressible& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const
{
    double rhoL = wl.density();
    double pL = thermoL.pressure();
    double nuL = wl.normalVelocity(normalVector);
    double aL = thermoL.soundSpeed();

    double rhoR = wr.density();
    double pR = thermoR.pressure();
    double nuR = wr.normalVelocity(normalVector);
    double aR = thermoR.soundSpeed();
    
    //PVRS
    double pm = std::fmax(0, 0.5*(pL + pR) - 0.5*(nuR - nuL)*0.5*(rhoL + rhoR)*0.5*(aL + aR));

    double sl;
    double sr;
    double sm;

    if (pm <= pL)
    {
        sl = nuL - aL;
    }        
    else
    {
        sl = 0.5*(nuL + nuR) - 0.5*(aL + aR);
    }
        
    if (pm <= pR)
    {
        sr = nuR + aR;
    }        
    else
    {
        sr = 0.5*(nuL + nuR) + 0.5*(aL + aR);
    }
        
    //contact wave speed
    sm = pR - pL + rhoL*nuL*(sl - nuL) - rhoR*nuR*(sr - nuR);
    sm = sm / (rhoL*(sl - nuL) - rhoR*(sr - nuR));


    /*Vars<3> wSpeed = waveSpeedsEstimate2(wl, wr, normalVector);
    double sl = wSpeed[0];
    double sr = wSpeed[2];
    double sm = wSpeed[1];
    double pm = 0.0;*/

    //HLLC scheme
    if (sl >= 0)
    {
        //left state
        return wl.flux(thermoL, normalVector);
    }
    else if (sr <= 0)
    {
        //right state
        return wr.flux(thermoR, normalVector);
    }
    else if (sm >= 0)
    {
        //middle-left state

        double uL = wl.velocityU();
        double vL = wl.velocityV();
        double wL = wl.velocityW();
        double EL = wl.totalEnergy();

        pm = pL + rhoL*(sl - nuL)*(sm - nuL);
        double rhoM = rhoL*(sl - nuL)/(sl - sm);

        return Vars<5>({rhoL*nuL                         + sl*(rhoM - rhoL),
                        rhoL*nuL*uL + pL*normalVector[0] + sl*((rhoM - rhoL)*uL + (pm - pL)/(sl - sm)*normalVector[0]),
                        rhoL*nuL*vL + pL*normalVector[1] + sl*((rhoM - rhoL)*vL + (pm - pL)/(sl - sm)*normalVector[1]),
                        rhoL*nuL*wL + pL*normalVector[2] + sl*((rhoM - rhoL)*wL + (pm - pL)/(sl - sm)*normalVector[2]),
                        rhoL*nuL*(EL + pL/rhoL)          + sl*((rhoM - rhoL)*EL + (pm*sm - pL*nuL)/(sl - sm))});
    }
    else
    {
        //middle-right state

        double uR = wr.velocityU();
        double vR = wr.velocityV();
        double wR = wr.velocityW();
        double ER = wr.totalEnergy();

        pm = pR + rhoR*(sr - nuR)*(sm - nuR);
        double rhoM = rhoR*(sr - nuR)/(sr - sm);
        return Vars<5>({rhoR*nuR + sr*(rhoM - rhoR),
                        rhoR*nuR*uR + pR*normalVector[0] + sr*((rhoM - rhoR)*uR + (pm - pR)/(sr - sm)*normalVector[0]),
                        rhoR*nuR*vR + pR*normalVector[1] + sr*((rhoM - rhoR)*vR + (pm - pR)/(sr - sm)*normalVector[1]),
                        rhoR*nuR*wR + pR*normalVector[2] + sr*((rhoM - rhoR)*wR + (pm - pR)/(sr - sm)*normalVector[2]),
                        rhoR*nuR*(ER + pR/rhoR) + sr*((rhoM - rhoR)*ER + (pm*sm - pR*nuR)/(sr - sm))});
    }
}

Vars<9> Hllc::claculateFlux(const CompressibleMixture& wl, const CompressibleMixture& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const
{
    double rhoL = wl.density();
    double pL = thermoL.pressure();
    double nuL = wl.normalVelocity(normalVector);
    double aL = thermoL.soundSpeed();

    double rhoR = wr.density();
    double pR = thermoR.pressure();
    double nuR = wr.normalVelocity(normalVector);
    double aR = thermoR.soundSpeed();
    
    //PVRS
    double pm = std::fmax(0, 0.5*(pL + pR) - 0.5*(nuR - nuL)*0.5*(rhoL + rhoR)*0.5*(aL + aR));

    double sl;
    double sr;
    double sm;

    if (pm <= pL)
    {
        sl = nuL - aL;
    }        
    else
    {
        sl = 0.5*(nuL + nuR) - 0.5*(aL + aR);
    }
        
    if (pm <= pR)
    {
        sr = nuR + aR;
    }        
    else
    {
        sr = 0.5*(nuL + nuR) + 0.5*(aL + aR);
    }
        
    //contact wave speed
    sm = pR - pL + rhoL*nuL*(sl - nuL) - rhoR*nuR*(sr - nuR);
    sm = sm / (rhoL*(sl - nuL) - rhoR*(sr - nuR));

    //HLLC scheme
    if (sl >= 0)
    {
        //left state
        return wl.flux(thermoL, normalVector);
    }
    else if (sr <= 0)
    {
        //right state
        return wr.flux(thermoR, normalVector);
    }
    else if (sm >= 0)
    {
        //middle-left state

        double uL = wl.velocityU();
        double vL = wl.velocityV();
        double wL = wl.velocityW();
        double EL = wl.totalEnergy();

        double Q0 = wl.Q0();
        double Q1 = wl.Q1();
        double Q2 = wl.Q2();
        double Q3 = wl.Q3();

        pm = pL + rhoL*(sl - nuL)*(sm - nuL);
        double rhoM = rhoL*(sl - nuL)/(sl - sm);

        return Vars<9>({rhoL*nuL                         + sl*(rhoM - rhoL),
                        rhoL*nuL*uL + pL*normalVector[0] + sl*((rhoM - rhoL)*uL + (pm - pL)/(sl - sm)*normalVector[0]),
                        rhoL*nuL*vL + pL*normalVector[1] + sl*((rhoM - rhoL)*vL + (pm - pL)/(sl - sm)*normalVector[1]),
                        rhoL*nuL*wL + pL*normalVector[2] + sl*((rhoM - rhoL)*wL + (pm - pL)/(sl - sm)*normalVector[2]),
                        rhoL*nuL*(EL + pL/rhoL)          + sl*((rhoM - rhoL)*EL + (pm*sm - pL*nuL)/(sl - sm)),
                        rhoL*nuL*Q0                      + sl*((rhoM - rhoL)*Q0),
                        rhoL*nuL*Q1                      + sl*((rhoM - rhoL)*Q1),
                        rhoL*nuL*Q2                      + sl*((rhoM - rhoL)*Q2),
                        rhoL*nuL*Q3                      + sl*((rhoM - rhoL)*Q3)});
    }
    else
    {
        //middle-right state

        double uR = wr.velocityU();
        double vR = wr.velocityV();
        double wR = wr.velocityW();
        double ER = wr.totalEnergy();

        double Q0 = wr.Q0();
        double Q1 = wr.Q1();
        double Q2 = wr.Q2();
        double Q3 = wr.Q3();

        pm = pR + rhoR*(sr - nuR)*(sm - nuR);
        double rhoM = rhoR*(sr - nuR)/(sr - sm);
        return Vars<9>({rhoR*nuR + sr*(rhoM - rhoR),
                        rhoR*nuR*uR + pR*normalVector[0] + sr*((rhoM - rhoR)*uR + (pm - pR)/(sr - sm)*normalVector[0]),
                        rhoR*nuR*vR + pR*normalVector[1] + sr*((rhoM - rhoR)*vR + (pm - pR)/(sr - sm)*normalVector[1]),
                        rhoR*nuR*wR + pR*normalVector[2] + sr*((rhoM - rhoR)*wR + (pm - pR)/(sr - sm)*normalVector[2]),
                        rhoR*nuR*(ER + pR/rhoR) + sr*((rhoM - rhoR)*ER + (pm*sm - pR*nuR)/(sr - sm)),
                        rhoR*nuR*Q0                      + sr*((rhoM - rhoR)*Q0),
                        rhoR*nuR*Q1                      + sr*((rhoM - rhoR)*Q1),
                        rhoR*nuR*Q2                      + sr*((rhoM - rhoR)*Q2),
                        rhoR*nuR*Q3                      + sr*((rhoM - rhoR)*Q3)});
    }
}



Vars<3> Hllc::waveSpeedsEstimate2(const Compressible& wl, const Compressible& wr, const ThermoVar& thermoL, const ThermoVar& thermoR, const Vars<3>& normalVector) const
{
    double al = thermoL.soundSpeed();
    double ar = thermoR.soundSpeed();
    double rhol = wl.density();
    double rhor = wr.density();
    double ul = wl.normalVelocity(normalVector);
    double ur = wr.normalVelocity(normalVector);
    double rholsqrt = sqrt(rhol);
    double rhorsqrt = sqrt(rhor);

    double n2 = 0.5*((rholsqrt*rhorsqrt)/std::pow(rholsqrt + rhorsqrt, 2));
    double d = sqrt((rholsqrt*al*al + rhorsqrt*ar*ar)/(rholsqrt + rhorsqrt) + n2*std::pow(ur - ul, 2));

    double uAvg = (rholsqrt*ul + rhorsqrt*ur)/(rholsqrt + rhorsqrt);

    double sl = uAvg - d;
    double sr = uAvg + d;
    double ss = (thermoR.pressure() - thermoL.pressure() + rhol*ul*(sl - ul) - rhor*ur*(sr - ur))/(rhol*sl - rhol*ul - rhor*sr + rhor*ur);

    return Vars<3>({sl, ss, sr});
}