
#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "AbstractElement.hpp"
#include "Node.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "ChastePoint.hpp"
#include "Debug.hpp"
#include "AbstractConductivityModifier.hpp"

#include "../src/NeuralComponents.hpp"
#include "NeuralComponents.hpp"
//#include "ChasteCuboid.hpp" //uncommend if using a circle shaped activation
#include "ChasteEllipsoid.hpp"

#include <fstream>
#include <vector>
#include <set>
#include <cassert> // standard debugging tool (evaluation assertion)
#include <cmath> // for sqrt
#include "../src/CellICCBioPhy.hpp"
#include "../src/Imtiaz2002_neural.hpp"
#include "../src/DummyCell.hpp"

// *************************** CELL FACTORY ************************************* //
class ICCCellFactory : public AbstractCardiacCellFactory<3>
{
private:
	std::set<unsigned> setICCNode;
    HistogramData& excitatory;
    HistogramData& inhibitory;
public:
    ICCCellFactory(std::set<unsigned> iccNodes, HistogramData& excitData, HistogramData& inhibData) //read list of ICC from the test
        : AbstractCardiacCellFactory<3>(), setICCNode(iccNodes), excitatory(excitData), inhibitory(inhibData)
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)

    {
        ChastePoint<3> centre(-6.0, -11.0, -31.0);
        ChastePoint<3> radii(3.0, 3.0, 3.0);

        ChasteEllipsoid<3> pacemaker(centre, radii);

        // // Define pacemaker region
        // double x = pNode->rGetLocation()[0];
        // double y = pNode->rGetLocation()[1];
        // double z = pNode->rGetLocation()[2];

        double r = 200e-4; // set size of the radius
				double scale = 2e-4; // Set the same like for mesh.scale

        // Find all nodes which are ICC by reading the list
        //".find" tries to locate an element in a set. Returns iterator point to the sought-after element
        // otherwise .end() if not found
        unsigned index = pNode->GetIndex();
        if (setICCNode.find(index) != setICCNode.end())
        {

            SingleModParam E_K(new ModifiableParams("E_K", -70.0));
            SingleModParam eta(new ModifiableParams("eta", 0.0389));
            SingleModParam excitatory_neural(new ModifiableParams("excitatory_neural", 0.000975));
            SingleModParam inhibitory_neural(new ModifiableParams("inhibitory_neural", 1.2));

            ModParamHolder modParams = {E_K, eta, excitatory_neural, inhibitory_neural} ;


            // ICC Cells NOT IN pacemaker region
            Cellimtiaz_2002d_noTstart_CORFromCellML* cell = new Cellimtiaz_2002d_noTstart_CORFromCellML(mpSolver, mpZeroStimulus, modParams);

            cell->SetParameter("cor", 1.0);

            // Active ICC Cells inside the pacemaker region (circle shaped over the whole z-depth)
            if  (pacemaker.DoesContain(pNode->GetPoint()))
            //  for a quadratic pacemaker region: if (x>=(306.968*scale-100e-4) && y>=(307.1210*scale-100e-4))
            {
                cell->SetParameter("cor", 1.2);
            }
            //Sets ICC
            return cell;
         }
        // All other cells which are not ICC or Bath
         else
         {
             CellDummyCellFromCellML* i_cell = new CellDummyCellFromCellML(mpSolver, mpZeroStimulus);
             return i_cell;
         }
    }
};

// *************************** CONDUCTIVITY MODIFIER ************************************* //
// class ICCConductivityModifier : public AbstractConductivityModifier<3,3>
// {
// private:
//     std::set<unsigned> setICCElement; //Copy list of Indexes with ICC attributes to this class
//     c_matrix<double, 3,3> mSpecialMatrix;
//     c_matrix<double, 3,3> mSpecialMatrix2;
//     c_matrix<double, 3,3> mSpecialMatrix3;
//     c_matrix<double, 3,3> mSpecialMatrix4;

// public:
//     ICCConductivityModifier(std::set<unsigned> elementIndexesICC)
//         : AbstractConductivityModifier<3,3>(),
//           setICCElement(elementIndexesICC),
//           mSpecialMatrix( zero_matrix<double>(3,3) ),
//           mSpecialMatrix2( zero_matrix<double>(3,3) ),
//           mSpecialMatrix3( zero_matrix<double>(3,3) ),
//           mSpecialMatrix4( zero_matrix<double>(3,3) )
//           {
//             //Conductivities for ICC and Dummy (Bath conductivity is set at HeartConfig::Instance())
//             double intraICC = 0.024; // Intracellular ICC [mS/cm] = [S/m]*e-2
//             double extraICC = 0.036; // Extracellular ICC
//             double intraDummy = 0.02; // Intracellular Dummy
//             double extraDummy = 0.02; // Extracellular Dummy

// 			// Introducing faster propagation in x direction set mSpecialMatrix(0,0) higher
// 			// e.g. intraICC + 0.05 and extraICC + 0.05
//               mSpecialMatrix(0,0) = intraICC; // Intracellular ICC Conductivities
//               mSpecialMatrix(1,1) = intraICC;
//               mSpecialMatrix(2,2) = intraICC;

//               mSpecialMatrix2(0,0) = extraICC; // Extracellular ICC Conductivities
//               mSpecialMatrix2(1,1) = extraICC;
//               mSpecialMatrix2(2,2) = extraICC;

//               mSpecialMatrix3(0,0) = intraDummy; // Intracellular Dummy Conductivities
//               mSpecialMatrix3(1,1) = intraDummy;
//               mSpecialMatrix3(2,2) = intraDummy;

//               mSpecialMatrix4(0,0) = extraDummy; // Extracellular Dummy Conductivities
//               mSpecialMatrix4(1,1) = extraDummy;
//               mSpecialMatrix4(2,2) = extraDummy;

//           }

//     c_matrix<double,3,3>& rCalculateModifiedConductivityTensor(unsigned elementIndex,
//                                                                const c_matrix<double,3,3>& rOriginalConductivity,
//                                                                unsigned domainIndex)
//     {
//         // Conductivities for ICC
//         // depending on the `domainIndex` (intra/extracellular).
//         if (setICCElement.find(elementIndex) != setICCElement.end())
//         {
//             if (domainIndex == 0) // domainIndex==0 implies intracellular
//             {
//                 return mSpecialMatrix;
//             }
//             else //domain Index==1 implies extracellular
//             {
//                 return mSpecialMatrix2;
//             }

//         }
//         // Conductivities for Dummy cells
//         else
//         {
//             if (domainIndex == 0) // domainIndex==0 implies intracellular
//             {
//                 return mSpecialMatrix3;
//             }
//             else //domain Index==1 implies extracellular
//             {
//                 return mSpecialMatrix4;
//             }

//         }
//     }
// };

// *************************** TEST ************************************* //
class TestICC3D : public CxxTest::TestSuite
{
public:
    void TestMesh3D() //throw(Exception)
    {
        // Read mesh created by TetGen
	    TrianglesMeshReader<3,3> reader("projects/mesh/ICC3D/sel_elem_cm_layer_ventCorpus");

        DistributedTetrahedralMesh<3,3> mesh; // Data shared among processes if run in parallel
        mesh.ConstructFromMeshReader(reader);

        // Scale mesh from 512x512x22pix to 303.64x303.64x22um by 1e-4: not needed for rat stomach
        // mesh.Scale(2e-4, 2e-4, 2e-4); // double of original image size (because of the Jacobian and converging errors)
        // mesh.SetMeshHasChangedSinceLoading();

        // Create list with ICC indexes
        std::set<unsigned> iccNodes;
        // std::set<unsigned> elementIndexesICC;


        TRACE("3");
        // Define boundary nodes
        double eleIdentify = 0;
        for (DistributedTetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin(); iter != mesh.GetElementIteratorEnd();
            ++iter)
        {
            eleIdentify = iter->GetAttribute();
            if (eleIdentify == 1) // ICC=1 and Bath=0
            {
                if(!iter->GetNode(0)->IsBoundaryNode())
                {
                    iccNodes.insert(iter->GetNodeGlobalIndex(0));
                }
                if(!iter->GetNode(1)->IsBoundaryNode())
                {
                    iccNodes.insert(iter->GetNodeGlobalIndex(1));
                }
                if(!iter->GetNode(2)->IsBoundaryNode())
                {
                    iccNodes.insert(iter->GetNodeGlobalIndex(2));
                }
            }
        }


        // // Iterating trough all Elements in the mesh and assigning attributes, conductivities and saving all ICC nodes
        // for (DistributedTetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
        //                 iter != mesh.GetElementIteratorEnd(); ++iter)
        // {
        //     // Read Attributes
        //     double attribute = iter->GetAttribute();

        //     // Set all small islands in the mesh to Dummy cells
        //     if (attribute != 1 && attribute != 2 && attribute != 3) // ICC=1(red), Dummy=2(green), Bath=3(blue)
        //     {
        //         iter->SetAttribute(2);
        //         attribute = iter->GetAttribute();
        //     }

        //     // Copy all nodes of the element to the elementIndexesICC list
        //     if (attribute == 1) // Check if ICC node
        //     {
        //         elementIndexesICC.insert(iter->GetIndex());
        //         for(int j = 0; j<=3; ++j)
        //         {
        //             iccNodes.insert(iter->GetNodeGlobalIndex(j));
        //         }
        //     }
        // }

        // Check that if we're in parallel no single process owns every element (to ensure that the conductivities
        // really are distributed).
        if (PetscTools::IsParallel())
               {
                   TS_ASSERT_DIFFERS( mesh.GetNumElements(), mesh.GetNumLocalElements() );
               }

 // Setting Attributes
        std::set<unsigned> background_ids;
        static unsigned background_id1 = 0; // Set Bath
        background_ids.insert(background_id1);
        // Uncomment for only ICC and Bath:
        //static unsigned background_id2 = 2; // Set Dummy as Bath
        //background_ids.insert(background_id2);

		// Set Attributes for ICC and Dummy
        std::set<unsigned> ICC_ids;
        static unsigned ICC_id1 = 1; // Set ICCs
        ICC_ids.insert(ICC_id1);
        // static unsigned ICC_id2 = 2; // Set Dummy (Comment for ICC-Bath simulation without Dummy)
        // ICC_ids.insert(ICC_id2);

        HeartConfig::Instance()->SetTissueAndBathIdentifiers(ICC_ids, background_ids); // tissue and bath ids

    	// Set Information for simulation
        HeartConfig::Instance()->SetSimulationDuration(20000); //ms for one cycle 10,000
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1, 1, 10); //timesteps: ode, pde, printing
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000); // Ratio for each cell
        HeartConfig::Instance()->SetUseAbsoluteTolerance(2e-4); //Changed to get around the DIVERGED_ITS error default:2e-4
        HeartConfig::Instance()->SetCapacitance(3); // Membrane Capacitance
        HeartConfig::Instance()->SetBathConductivity(0.02); // Bath capacitance

		// Set outputfile name
        HeartConfig::Instance()->SetOutputDirectory("RatMesh3D_Imtiaz");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true); // Set for visualizing with Meshlab
				//HeartConfig::Instance()->SetVisualizeWithCmgui(true);

        // Load neural info
        HistogramData eData = HistogramData("projects/NeuralData/hist_ex_mean.txt", 20, 30, 1, 0.0639, 0.0639);
        HistogramData iData = HistogramData("projects/NeuralData/hist_in_mean.txt", 20, 30, 1, 0.0639, 0.0639);

        // Initialize cell_factory
        ICCCellFactory cell_factory(iccNodes, eData, iData);

        // Declare the problem class, `BidomainProblem<3>`
        BidomainProblem<3> bidomain_problem( &cell_factory, true); // true indicates we are solving a bath problem

        // When not used 'HeartConfig' for reading the mesh. Has to be called before Initialise
        bidomain_problem.SetMesh(&mesh);

        //  min/max voltage is printed as the simulation runs, useful for verifying that cells are stimulated
        bidomain_problem.SetWriteInfo();

        //Initialise
        bidomain_problem.Initialise(); // Has to be initialised before setting the conductivities

        // Set conductivities with the Conductivity Modifier for every element
        BidomainTissue<3>* p_bidomain_tissue = bidomain_problem.GetBidomainTissue();
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.12, 0.12)); // these are quite smaller than cm values
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.2, 0.2)); // these are quite smaller than cm values
        // ICCConductivityModifier modifier(elementIndexesICC); // Initialise Conductivity Modifier
        // p_bidomain_tissue->SetConductivityModifier( &modifier );

        //Solve
        bidomain_problem.Solve();

        // Switch on the output into log.
        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};
