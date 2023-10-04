#include "MEElement.hpp"
#include "CPPProcess.h"
#include "rambo.h"

class MEElementMG5 : public MEElement
{
    private:
        /* data */
        CPPProcess *_cppp;

    protected:
        PhysConstants constants;

    public:
        MEElementMG5(/* args */);
        ~MEElementMG5();

        double* multi(double momenta[], int n_elements, double buffer[]);
        double single(std::vector<double*> momenta);
        double rambo();
};

MEElementMG5::MEElementMG5(/* args */)
{
}

MEElementMG5::~MEElementMG5()
{
}
