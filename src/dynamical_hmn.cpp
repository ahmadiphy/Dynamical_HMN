#include "dynamical_hmn.h"
using namespace std;

dynamical_HMN::dynamical_HMN(int inm0,int inl)
{
    m0=inm0;
    l=inl;
    ll=pow(2,l);
    outll=ll;
    outm0=m0;
    iniLinks=2;
    for(int i=0;i<l;++i)
    {
        Links.push_back(iniLinks);
        PossibleLinks.push_back(pow(4,i)*pow(m0,2));
    }
    Links.shrink_to_fit();
    PossibleLinks.shrink_to_fit();
}

void dynamical_HMN::Cal_inverseM(iMatrix &aa, iMatrix &inv)
{
    int modulB,modulE;
    for(int i=0;i<ll;++i)
    {
        modulB=i*m0;
        modulE=(i+1)*m0-1;
        for(int j=modulB;j<=modulE;++j)
        {
            for(int jk=modulB;jk<=modulE;++jk)
            {
                bool state=0;
                for(int k=0;k<aa[j].size();++k)
                {
                    if(aa[j][k]==jk || j==jk)
                    {
                        state=1;
                    }
                }
                if(state==0)
                {
                    inv[j].push_back(jk);
                }
            }
        }

    }
}

/*
void dynamical_HMN::Cal_modules_links(iMatrix& aa)
{
    int modulB,modulE;
    for(int i=0;i<ll;++i)
    {
        modulB=i*m0;
        modulE=(i+1)*m0;
        Links[i]=0;
        for(int j=modulB;j<modulE;++j)
        {
            Links[i]+=aa[j].size();
        }
        PossibleLinks[i]=m0-1-Links[i];

    }
}
*/

void dynamical_HMN::intraModule(iMatrix& aa)
{
    int nn=m0;
    int nl=0;
    int modulB,modulE;
    for(int i=0;i<ll;++i)
    {
        modulB=i*nn;
        modulE=(i+1)*nn;
        for(int ii=0;ii<nn;++ii)
        {
            for(int j=modulB;j<nl;++j)
            {
                aa[nl].push_back(j);
            }
            for(int j=nl+1;j<modulE;++j)
            {
                aa[nl].push_back(j);
            }
            nl++;
        }

    }
}

void dynamical_HMN::intraModuleLine(iMatrix& aa)
{
    int modulB,modulE;
    for(int i=0;i<ll;++i)
    {
        modulB=i*m0;
        modulE=(i+1)*m0;
        for(int j=modulB;j<modulE-1;++j)
        {
            aa[j].push_back(j+1);
            aa[j+1].push_back(j);
        }

    }
}

void dynamical_HMN::intraModuleRing(iMatrix& aa)
{
    int modulB,modulE;
    for(int i=0;i<ll;++i)
    {
        modulB=i*m0;
        modulE=(i+1)*m0;
        for(int j=modulB;j<modulE;++j)
        {
            aa[j].push_back(j+1);
            aa[j+1].push_back(j);
        }

    }
}

void dynamical_HMN::interModule(iMatrix &aa)
{
    int nn=m0;
    int alpha=iniLinks;
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist(0,1);
    int ip=1;
    for(int ii=0;ii<l;ii++)
    {
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            uniform_int_distribution<> dis1(dis1B, dis1E);
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            uniform_int_distribution<> dis2(dis2B, dis2E);
            int tch=0;
            while(tch<alpha){
                int r1=dis1(gen);
                int r2=dis2(gen);
                int rcheck=0;
                for(int ri=0;ri<aa[r1].size();++ri)
                {
                    if(aa[r1][ri]==r2)
                        rcheck=1;
                }
                if(rcheck==0)
                {
                    aa[r1].push_back(r2);
                    aa[r2].push_back(r1);
                    tch++;
                }
            }
        }

    }
}

void dynamical_HMN::uniformHIMC(iMatrix &aa, int iLink)
{
    random_device rd;//random device to randomize tha seed
    mt19937 gen(rd());  // to seed mersenne twister.
    iMatrix inverseM(ll*m0,iRow());
    Matrixf inmf;
    Cal_inverseM(aa,inverseM);
    inmf.printMatrixForm2File(inverseM,"inmf");
    for(int i=0;i<ll;++i)
    {
        int modulB,modulE;
        modulB=i*m0;
        modulE=(i+1)*m0-1;
        vector<int> vec1;
        for(int j=modulB;j<=modulE;++j)
        {
            if(inverseM[j].size()!=0)
                vec1.push_back(j);
        }
        vec1.shrink_to_fit();
        for(int ln=0;ln<iLink && vec1.size()!=0;++ln)
        {

            uniform_int_distribution<> dist1(0, vec1.size()-1);
            int r1=dist1(gen);
            int targetI=vec1[r1];
            int invS=inverseM[targetI].size();
            uniform_int_distribution<> dist2(0, invS-1);
            int r2=dist2(gen);
            int targetJ=inverseM[targetI][r2];
            aa[targetI].push_back(targetJ);
            aa[targetJ].push_back(targetI);
            vec1.erase(vec1.begin() + r1);
            inverseM[targetI].erase(inverseM[targetI].begin()+r2);
            for(int id=0;id<inverseM[targetJ].size();++id)
            {
                if(inverseM[targetJ][id]==targetI)
                    inverseM[targetJ].erase(inverseM[targetJ].begin()+id);
            }

        }
    }
    inmf.printMatrixForm2File(inverseM,"_final_inverce");
}

void dynamical_HMN::uniformNHIMC(iMatrix &aa, int iLink)
{

}
