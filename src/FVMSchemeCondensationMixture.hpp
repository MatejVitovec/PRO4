#ifndef FVMSCHEMECONDENSATIONMIXTURE_HPP
#define FVMSCHEMECONDENSATIONMIXTURE_HPP

#include <vector>
#include <memory>

#include "Condensation/CompressibleMixture.hpp"
#include "ThermoVar.hpp"
#include "ComponentThermoVar.hpp"
#include "Field.hpp"
#include "VolField.hpp"
#include "Mesh/Mesh.hpp"
#include "FluxSolver/FluxSolver.hpp"
#include "Thermo/Thermo.hpp"
#include "Thermo/IapwsLiquidThermo.hpp"
#include "GradientScheme/LeastSquare.hpp"
#include "GradientScheme/LeastSquaresPS.hpp"
#include "Limiter/Limiter.hpp"
#include "Limiter/BarthJespersen.hpp"
#include "Limiter/Venkatakrishnan.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"

#include "Condensation/Mixture.hpp"

class FVMSchemeCondensationMixture
{
    public:

        //vznikne obraz meshe - neprekopiruje se
        //FVMSchemeCondensationMixture() : mesh(Mesh()), w(Field<CompressibleMixture>()), wl(Field<CompressibleMixture>()), wr(Field<CompressibleMixture>()) {}
        FVMSchemeCondensationMixture(Mesh&& mesh_, std::unique_ptr<FluxSolver> fluxSolver_, std::unique_ptr<Thermo> thermo_) : mesh(std::move(mesh_)),
                                                                                                            fluxSolver(std::move(fluxSolver_)),
                                                                                                            thermo(std::move(thermo_)),
                                                                                                            w(VolField<CompressibleMixture>(mesh)),
                                                                                                            mixtureThermoField(VolField<ThermoVar>(mesh)),
                                                                                                            vaporThermoField(VolField<ComponentThermoVar>(mesh)),
                                                                                                            liquidThermoField(VolField<ComponentThermoVar>(mesh)),
                                                                                                            wl(Field<CompressibleMixture>()),
                                                                                                            wr(Field<CompressibleMixture>()),
                                                                                                            mixtureThermoFieldL(Field<ThermoVar>()),
                                                                                                            mixtureThermoFieldR(Field<ThermoVar>()),
                                                                                                            cfl(0.8),
                                                                                                            maxIter(10000000),
                                                                                                            targetError(0000005),
                                                                                                            localTimeStep(false),
                                                                                                            reconstruction(false),
                                                                                                            time(0.0),
                                                                                                            gradientScheme(std::make_unique<GradientScheme>()),
                                                                                                            limiter(std::make_unique<Limiter>()) {}


        virtual ~FVMSchemeCondensationMixture() {}

        void setReconstructionGradient(std::unique_ptr<GradientScheme> gradScheme_);
        void setReconstructionLimiter(std::unique_ptr<Limiter> limiter_);

        void setCfl(double cfl_);
        void setMaxIter(int maxIter_);
        void setSaveEveryIter(int saveEveryIter_);
        void setTargetError(double targetError_);
        void setLocalTimeStep(bool localTimeStep_);
        void setReconstructionSettings(bool reconstruction_);
        void setSavePath(std::string path) {savePath = path;}

        void setThermoModel(std::unique_ptr<Thermo> thermo_) {thermo = std::move(thermo_);}
        
        double getCfl() const;
        int getMaxIter() const;
        double getTargetError() const;
        bool getTimeStepsettings() const;
        bool getReconstructionSettings() const;
        std::string getSavePath() const { return savePath; }
        const Mesh& getMesh() const;
        const Thermo* getThermoRef();

        void setInitialConditions(Compressible initialCondition);
        void setInitialConditionsPrimitive(Vars<5> initialCondition);

        void setBoundaryConditions(std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditions);
        void applyFreeBoundaryCondition();

        void init();
        virtual void solve() = 0;

        Field<CompressibleMixture> getResults() const;
        Field<ThermoVar> getResultsThermo() const;
        
        
    protected:
        std::unique_ptr<FluxSolver> fluxSolver;

        Mixture mixture;
        std::unique_ptr<Thermo> thermo;
        IapwsLiquidThermo liquidThermo; //TODO na smart pointer

        std::unique_ptr<GradientScheme> gradientScheme;
        std::unique_ptr<Limiter> limiter;

        

        Mesh mesh;
        std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditionList;

        VolField<CompressibleMixture> w; //cell size
        VolField<ThermoVar> mixtureThermoField;
        VolField<ComponentThermoVar> vaporThermoField;
        VolField<ComponentThermoVar> liquidThermoField;
        //std::vector<std::vector<CompressibleMixture>> boundaryFields;

        Field<CompressibleMixture> wl; //faces size
        Field<CompressibleMixture> wr;

        Field<ThermoVar> mixtureThermoFieldL; //faces size
        Field<ThermoVar> mixtureThermoFieldR;

        Field<Vars<9>> fluxes;

        Field<Mat<9,3>> grad;
        Field<Vars<9>> phi;

        double cfl;
        int maxIter;
        int saveEveryIter;
        double targetError;
        bool localTimeStep;
        bool reconstruction;
        std::string savePath;

        int fixGradient;
        int fixedGradStep;

        Field<double> timeSteps;

        double time;
        int iter;

        void updateTimeStep();
        void calculateWlWr();
        void interpolateToFaces();
        void calcBoundaryConditionFields();
        void calculateFluxes();
        Field<Vars<9>> calculateResidual();
        void boundField();

        //std::vector<std::vector<CompressibleMixture>> calcBoundaryConditionsToBoundaryFields();

    private:

};



#endif // FVMSCHEMECONDENSATIONMIXTURE_HPP