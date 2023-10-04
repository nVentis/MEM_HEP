#include "MEMTransferFunction.hpp"
#include "cmath"

/**
 * params[0]: xo
 * params[1]: gamma
*/
class MEMTransferFunctionLorentz: public MEMTransferFunction
{
    private:
        /* data */
        
    public:
        using MEMTransferFunction::MEMTransferFunction;

        double calc_tf(double reco_val, double true_val);
};

double MEMTransferFunctionLorentz::calc_tf(double reco_val, double true_val) {
    return 1/(M_PI*params[1]*(1 + std::pow( ((reco_val - true_val) - params[0])/params[1], 2.) ));
}