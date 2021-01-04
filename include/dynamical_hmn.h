#include "in.h"
#ifndef DYNAMICAL_HMN_H
#define DYNAMICAL_HMN_H


class dynamical_HMN
{
public:
    dynamical_HMN(int inm0, int inl);
    int outll,outm0;//to print number of modules and nodules sizes
    void intraModule(iMatrix& aa);//building hierarchical modules
    void intraModuleLine(iMatrix& aa);//building linear modules
    void intraModuleRing(iMatrix& aa);//building ring modules
    void uniformHIMC(iMatrix& aa, int iLink);//uniform and hemogenius increasing modular connectivity
    void uniformNHIMC(iMatrix& aa, int iLink);//uniform and non hemogenius increasing modular connectivity
    void interModule(iMatrix &aa);//connecting modules together in hierarchical structure
private:
    int m0, l, ll, iniLinks;//initial conditions
    std::vector<int> Links;
    std::vector<int> PossibleLinks;
    void Cal_inverseM(iMatrix& aa, iMatrix& inv);
    //void Cal_modules_links(iMatrix& aa);
};

#endif // DYNAMICAL_HMN_H
