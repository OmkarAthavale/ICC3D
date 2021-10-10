#ifndef NEURALCOMPONENTS_HPP_
#define NEURALCOMPONENTS_HPP_

#include <boost/multi_array.hpp>
#include <unordered_map>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>

#include "ChasteSerialization.hpp"


class CalibrationFunctions
{
    public:
    static double All_FromData(double data);
    static double Beta_Zhang2011(double f_EFS);
    static double GBKmax_Kim2003(double f_EFS);
};

class HistogramData
{
    private:
    typedef boost::multi_array<double, 3> Hist3D;
    typedef Hist3D::index indexHist;
    Hist3D A;
    double xLen;
    double yLen;
    int xDivs;
    int yDivs;
    // int Get2DIndex(double val, double dimLength, double nLength);

    public:
    std::vector<double> GetValueOverTime(double xCoord, double yCoord, int numT);
    HistogramData(const std::string fName, int X, int Y, int T, double xL, double yL);

};
class ModifiableParams
{

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        //archive & mLocation; //earlier versions of boost are unable to do this. See #1709
    }

    private:
    std::string pName;
    double pInit;
    std::vector<double> pTimeDep;
    double tStep;
    double tMax;
    bool isTimeVarying;
    std::string funcName;
    typedef double (*CFunc)(double);
    CFunc calibFunc;

    std::unordered_map<std::string, CFunc> calibFuncMap {  {"All_FromData", &CalibrationFunctions::All_FromData}, 
                                                            {"Beta_Zhang2011", &CalibrationFunctions::Beta_Zhang2011},  
                                                            {"GBKmax_Kim2003", &CalibrationFunctions::GBKmax_Kim2003}
                                                        };

    public:
    ModifiableParams(const std::string &name, const double init, double timeStep, double timeMax, const std::vector<double> &pTimeDep, const std::string& fName);
    ModifiableParams(const std::string &name, const double init);
    double GetValue(double time) const;
    double GetValue(int index) const;
    double GetValue() const;

    // pure getters
    const std::string& GetName() const {return pName;};
    const std::string& GetFuncName() const {return funcName;};
    double GetStep() const {return tStep;};
    double GetMax() const {return tMax;};
    const std::vector<double>& GetVals() const {return pTimeDep;};
    bool GetTimeDepBool() const {return isTimeVarying;};

};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(ModifiableParams)

namespace boost
{
namespace serialization
{

template<class Archive>
inline void save_construct_data(
    Archive & ar, const ModifiableParams * t, const unsigned int file_version)
{

    ar & t->GetName();
    ar & t->GetValue(); //pInit
    ar & t->GetStep();
    ar & t->GetMax();
    ar & t->GetVals();
    ar & t->GetTimeDepBool();
    ar & t->GetFuncName();
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a ChastePoint instance (using existing constructor)
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, ModifiableParams * t, const unsigned int file_version)
{

    std::string pName;
    double pInit;
    std::vector<double> pTimeDep;
    double tStep;
    double tMax;
    bool isTimeVarying;
    std::string funcName;

    ar & pName;
    ar & pInit;
    ar & tStep;
    ar & tMax;
    ar & pTimeDep;
    ar & isTimeVarying;
    ar & funcName;

    if (isTimeVarying)
    {
        ::new(t)ModifiableParams(pName, pInit, tStep, tMax, pTimeDep, funcName);

    } else 
    {
        ::new(t)ModifiableParams(pName, pInit);
    }
}
}
} // namespace ...



#endif // NEURALCOMPONENTS_HPP_