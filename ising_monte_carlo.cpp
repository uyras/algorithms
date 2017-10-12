#include <iostream>
#include <vector>
#include <random>

using namespace std;

#define SEED 1000
#define L 4
#define N L*L*L
#define NEIGH 6
#define J -1

int E=0;
double T=0.01;
vector<int> spins(N),neigh(N*NEIGH);
default_random_engine generator;

void drop();
void eCalc();
void mc(unsigned long steps);
void dbg();

inline int nnum (int a, int b, int c){return a*L*L+b*L+c;}

int main()
{
    generator.seed(SEED);
    drop();
    eCalc();
    mc(1000000);
    dbg();

    cout << "E=" <<E<< endl;
    return 0;
}

void drop(){
    int num=0;
    for (int i=0;i<L;++i){
        for (int j=0;j<L;++j){
            for (int k=0;k<L;++k){
                //set all up
                //spins[num] = (i%2)?1:-1;
                //spins[num] = (i%2^j%2^k%2)?1:-1;
                //spins[num] = (j%2^k%2)?1:-1;
                spins[num] = 1; //all up


                //left neigh
                if (k==0)
                    neigh[num*NEIGH+0]=nnum(i,j,L-1);
                else
                    neigh[num*NEIGH+0]=nnum(i,j,k-1);

                //right neigh
                if (k==L-1)
                    neigh[num*NEIGH+1]=nnum(i,j,0);
                else
                    neigh[num*NEIGH+1]=nnum(i,j,k+1);

                //up neigh
                if (j==0)
                    neigh[num*NEIGH+2]=nnum(i,L-1,k);
                else
                    neigh[num*NEIGH+2]=nnum(i,j-1,k);

                //down neigh
                if (j==L-1)
                    neigh[num*NEIGH+3]=nnum(i,0,k);
                else
                    neigh[num*NEIGH+3]=nnum(i,j+1,k);

                //near neigh
                if (i==0)
                    neigh[num*NEIGH+4]=nnum(L-1,j,k);
                else
                    neigh[num*NEIGH+4]=nnum(i-1,j,k);

                //far neigh
                if (i==L-1)
                    neigh[num*NEIGH+5]=nnum(0,j,k);
                else
                    neigh[num*NEIGH+5]=nnum(i+1,j,k);

                ++num;
            }
        }
    }
}

void eCalc(){
    E=0;
    for (int i=0;i<N;++i){
        E -= spins[i] * (spins[neigh[i*NEIGH+0]] + spins[neigh[i*NEIGH+2]] + spins[neigh[i*NEIGH+4]]);
    }
    E*=J;
}

void mc(unsigned long steps){
    unsigned assump;
    int de;

    uniform_int_distribution<int> distr1(0,N-1);
    uniform_real_distribution<double> distr2(0.,1.);
    do{
        assump=distr1(generator);
        de = -2 * J * -spins[assump] *
                (spins[neigh[assump*NEIGH+0]] +
                spins[neigh[assump*NEIGH+1]] +
                spins[neigh[assump*NEIGH+2]] +
                spins[neigh[assump*NEIGH+3]] +
                spins[neigh[assump*NEIGH+4]] +
                spins[neigh[assump*NEIGH+5]]);

        if (exp(-de/T)>distr2(generator)){
            spins[assump]=-spins[assump];
            E+=de;
        }

        --steps;
    } while (steps>0);
}

void dbg(){
    for (int i=0;i<L;++i){
        for (int j=0;j<L;++j){
            for (int k=0;k<L;++k){
                cout<<i<<"\t"<<j<<"\t"<<k<<"\t"<<spins[nnum(i,j,k)];
                cout<<endl;
            }
        }
    }
}
