#include <cxxtest/TestSuite.h>
#include "CardiacSimulationArchiver.hpp"
#include "BidomainProblem.hpp"
#include "LuoRudy1991.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PlaneStimulusCellFactory.hpp"

class TestCheckpointingAndRestarting_RA : public CxxTest::TestSuite
{
public:
    void TestCheckpointing()
    {
        HeartConfig::Instance()->Reset();

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML,2> cell_factory(-2000000);
        HeartConfig::Instance()->SetSimulationDuration(15.0); //ms
        HeartConfig::Instance()->SetOutputDirectory("TestBidomainCheckpointing");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_800_elements", cp::media_type::Orthotropic);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        
        double scale = 2;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75*scale, 0.19*scale));
        BidomainProblem<2> bidomain_problem( &cell_factory );
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
        CardiacSimulationArchiver<BidomainProblem<2> >::Save(bidomain_problem, "TestBidomainCheckpointing/saved_simulation");
    }

    void TestRestarting()
    {
        BidomainProblem<2>* p_bidomain_problem = CardiacSimulationArchiver<BidomainProblem<2> >::Load("TestBidomainCheckpointing/saved_simulation");

        HeartConfig::Instance()->SetSimulationDuration(30); //ms

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(3.0, 0.3));

        p_bidomain_problem->Solve();

        delete p_bidomain_problem;
    }
};
