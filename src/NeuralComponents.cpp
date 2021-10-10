#include <cmath>
#include <fstream>
#include "NeuralComponents.hpp"
#include "Exception.hpp"

double CalibrationFunctions::All_FromData(double data)
{
    return data;
}

double CalibrationFunctions::Beta_Zhang2011(double f_EFS)
{
    double a = 0.000756965785697425;
    double b = 0.0001;
    
    double toHz = 500.0 * 0.001;
    double endings_per_ICC = 0.04;

    double eval_val = a * f_EFS * toHz * endings_per_ICC + b;

    if (eval_val > 0.0067)
    {
        return 0.0067;
    } else if (eval_val < 0.0001)
    {
        return 0.0001;
    } else
    {
        return eval_val;
    }
}

double CalibrationFunctions::GBKmax_Kim2003(double f_EFS)
{
    double a = 0.147525397773565;
    double b = -0.175087725001323;
    double c = 1.00016775704983;
    double d = 0.0285035244071441;

    double toHz = 500.0 * 0.01;
    double endings_per_ICC = 0.04;

    double eval_val = a*exp(b * f_EFS * toHz * endings_per_ICC) + c*exp(d * f_EFS * toHz * endings_per_ICC);

    if (eval_val > 2.5)
    {
        return 2.5;
    } else if (eval_val < 1.15)
    {
        return 1.15;
    } else
    {
        return eval_val;
    }
}

ModifiableParams::ModifiableParams(const std::string& name, const double init, double timeStep, double timeMax, const std::vector<double>& timeDep, const std::string& fName)
{
    this->pName = name;
    this->pInit = init;
    this->tStep = timeStep;
    this->tMax = timeMax;
    this->pTimeDep = timeDep;
    this->isTimeVarying = true;
    this->funcName = fName;
    this->calibFunc = calibFuncMap[funcName];
}

ModifiableParams::ModifiableParams(const std::string &name, const double init)
{
    this->pName = name;
    this->pInit = init;
    this->isTimeVarying = false;
}

double ModifiableParams::GetValue(int index) const
{
    if (this->isTimeVarying)
    {
        return this->calibFunc(pTimeDep[index]);
    }
    else
    {
        return this->pInit;
    }
}

double ModifiableParams::GetValue(double time) const
{

    int index = (int) ((fmod(time, tMax))/(tStep));
    if (this->isTimeVarying) 
    {
        return this->calibFunc(pTimeDep[index]);
    }
    else
    {
        return this->pInit;
    }
}

double ModifiableParams::GetValue() const
{
    return this->pInit;
}

std::vector<double> HistogramData::GetValueOverTime(double xCoord, double yCoord, int numT)
{
    int xInd = (xCoord/(xLen/xDivs)); // Get2DIndex(xCoord, 20.0, 2.0);
    int yInd = (yCoord/(yLen/yDivs)); // Get2DIndex(xCoord, 20.0, 2.0);
    typedef boost::multi_array_types::index_range range;

    std::vector<double> vect;

    for (double v : A[boost::indices[xInd][yInd][range(0,numT)]])
    {
        vect.push_back(v);
    }

    return vect;
}

HistogramData::HistogramData(const std::string fName, int X, int Y, int T, double xL, double yL):
A(boost::extents[X][Y][T]), xLen(xL), yLen(yL), xDivs(X), yDivs(Y)
{

    // Open file
    std::ifstream inFile;
    inFile.open(fName); 

    // Create a 3D array

    double input1;
    // Assign values to the elements
    for(indexHist k = 0; k != T; k++) 
    {
        for(indexHist j = 0; j != Y; j++)
        {
            for(indexHist i = 0; i != X; i++)
            {
                inFile >> input1;
                A[i][j][k] = input1;
            }
        }
    }

    //Close file
    inFile.close();
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ModifiableParams)