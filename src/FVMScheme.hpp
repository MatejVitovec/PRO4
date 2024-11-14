#ifndef FVMSCHEME_HPP
#define FVMSCHEME_HPP

#include <vector>
#include <memory>

#include "Compressible.hpp"
#include "Field.hpp"
#include "VolField.hpp"
#include "Mesh/Mesh.hpp"
#include "FluxSolver/FluxSolver.hpp"
#include "Thermo/Thermo.hpp"
#include "GradientScheme/LeastSquare.hpp"
#include "GradientScheme/LeastSquaresPS.hpp"
#include "Limiter/Limiter.hpp"
#include "Limiter/BarthJespersen.hpp"
#include "Limiter/Venkatakrishnan.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"

class FVMScheme
{
    public:

        //vznikne obraz meshe - neprekopiruje se
        //FVMScheme() : mesh(Mesh()), w(Field<Compressible>()), wl(Field<Compressible>()), wr(Field<Compressible>()) {}
        FVMScheme(Mesh&& mesh_, std::unique_ptr<FluxSolver> fluxSolver_, std::unique_ptr<Thermo> thermo_) : mesh(std::move(mesh_)),
                                                                                                            fluxSolver(std::move(fluxSolver_)),
                                                                                                            thermo(std::move(thermo_)),
                                                                                                            w(Field<Compressible>()),
                                                                                                            wl(Field<Compressible>()),
                                                                                                            wr(Field<Compressible>()),
                                                                                                            boundaryFields(),
                                                                                                            cfl(0.8), maxIter(10000000),
                                                                                                            targetError(0000005),
                                                                                                            localTimeStep(false),
                                                                                                            reconstruction(false),
                                                                                                            time(0.0),
                                                                                                            gradientScheme(std::make_unique<GradientScheme>()),
                                                                                                            limiter(std::make_unique<Limiter>()) {}


        virtual ~FVMScheme() {}

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

        Field<Compressible> getResults() const;
        
        
    protected:
        std::unique_ptr<FluxSolver> fluxSolver;
        std::unique_ptr<Thermo> thermo;
        std::unique_ptr<GradientScheme> gradientScheme;
        std::unique_ptr<Limiter> limiter;

        Mesh mesh;
        std::vector<std::shared_ptr<BoundaryCondition>> boundaryConditionList;

        VolField<Compressible> w; //cell size
        std::vector<std::vector<Compressible>> boundaryFields;

        Field<Compressible> wl; //faces size
        Field<Compressible> wr;
        Field<Vars<5>> fluxes;

        Field<Mat<5,3>> grad;
        Field<Vars<5>> phi;

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
        void applyBoundaryConditions();
        void calculateWlWr();
        void interpolateToFaces();
        void calcBoundaryConditionFields();
        void calculateFluxes();
        Field<Vars<5>> calculateResidual();
        void reconstruct();
        void boundField();

        std::vector<std::vector<Compressible>> calcBoundaryConditionsToBoundaryFields();

    private:

};



#endif // FVMSCHEME_HPP