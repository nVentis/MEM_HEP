#include "CalcME.hpp"

//----------- Implementation of class calc_me -----------//

calc_me::calc_me()
{
    _cppp = new CPPProcess();

    std::cout << _cppp->name() << std::endl;

    _cppp->initProc("mg5/zhh/Cards/param_card.dat");

    std::cout << " [TRACE] CPPP created Ok." << std::endl;
}

calc_me::~calc_me()
{
    std::cout << " [TRACE] CPPP deleted OK" << std::endl;
}

double calc_me::calc_me_calc()
{
    return 1.;
}

// Px,Py,Pz,E
double calc_me::calc_me_rambo()
{
    double energy = 1500;
    double weight;
    
    // Get phase space point
    std::vector<double*> p = get_momenta(_cppp->ninitial, energy,
                    _cppp->getMasses(), weight);

    // Set momenta for this event
    _cppp->setMomenta(p);

    // Evaluate matrix element
    _cppp->sigmaKin();

    const double* matrix_elements = _cppp->getMatrixElements();

    /*

    std::cout << "Momenta:" << std::endl;
    for(int i=0;i < _cppp->nexternal; i++)
    std::cout << std::setw(4) << i+1 
        << std::setiosflags(ios::scientific) << std::setw(14) << p[i][0]
        << std::setiosflags(ios::scientific) << std::setw(14) << p[i][1]
        << std::setiosflags(ios::scientific) << std::setw(14) << p[i][2]
        << std::setiosflags(ios::scientific) << std::setw(14) << p[i][3] << std::endl;
    std::cout << " -----------------------------------------------------------------------------" << std::endl;

    */

    return matrix_elements[0];
}

//---------- C-Interface for class lcmezhh ---------------------//

pStat calc_me_new()
{
    return new (std::nothrow) calc_me();
}

void calc_me_del(pStat self)
{
    delete reinterpret_cast<calc_me*>(self);
}

double calc_me_calc(pStat self)
{
    auto p = reinterpret_cast<calc_me*>(self);
    return p->calc_me_calc();
}

double calc_me_rambo(pStat self)
{
    auto p = reinterpret_cast<calc_me*>(self);
    return p->calc_me_rambo();
}
