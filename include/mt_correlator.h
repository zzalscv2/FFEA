#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <vector>
#include <fstream>
#include <string>
#include <random>
#include<sstream>
#include<unistd.h>
#include<cmath>
#include "mat_vec_types.h"
#include "mat_vec_fns_II.h"

using namespace std;
#ifndef MT_CORRELATOR_H
#define MT_CORRELATOR_H


const int ntpoints = 1000;

class TCorrelator{
    private:
        scalar *aa, *cor;
        long long *ncor;
        int numcorr, pcor,length,p2;
        void resizeCorrelator(int k);
    public:
        scalar *t,*f,*fav,*fsqav,*tav;
        int npcorr,npcorrmax,nexp;
        TCorrelator(int numcorrin = 35, int pcorin = 16);
        ~TCorrelator();
        TCorrelator(const TCorrelator &obj);
        int resizecount=0, addcount=0;
        int get_numcorr();
        void add(scalar w, int k=1);
        void init();
        void evaluate();
        void toaverage();
        void clear();
        void save(std::string savename);
        void save_out(FILE *fout);
        void read(std::string readname);
        void save_ffea(FILE *fout);
        void read_ffea(FILE *fout);
};

class TCorrelatorDiffusion{
    private:
        scalar *aa, *cor;
        long long *ncor;
        int numcorr, pcor,length,p2;
        void resizeCorrelator(int k);
    public:
        scalar *t,*f,*fav,*fsqav,*tav;
        int npcorr,npcorrmax,nexp;
        TCorrelatorDiffusion(int numcorrin = 35, int pcorin = 16);
        ~TCorrelatorDiffusion();
        TCorrelatorDiffusion(const TCorrelatorDiffusion &obj);
        int resizecount=0, addcount=0;
        void add(scalar w, int k=1);
        void init();
        void evaluate();
        void toaverage();
        int getlength();
        void clear();
        int get_numcorr();
        void save(std::string savename);
        void save_out(std::string savename);
        void save_ffea(FILE *fout);
        void read_ffea(FILE *fout);
        void read(std::string readname);
};

class TCrossCorrelator{
    private:
        scalar *aa,*aa2, *cor;
        long long *ncor;
        int numcorr, pcor,length,p2;
        void resizeCorrelator(int k);
    public:
        scalar *t,*f,*fav,*fsqav,*tav;
        int npcorr,npcorrmax,nexp;
        TCrossCorrelator(int numcorrin = 35, int pcorin = 16);
        ~TCrossCorrelator();
        int resizecount=0, addcount=0;
        void add(scalar w,scalar w2, int k=1);
        void init();
        void evaluate();
        void toaverage();
        void clear();
        void save(std::string savename);
        void read(std::string readname);
};

class TCorrelatorVector{
    private:
        arr3 *aa;
        scalar *cor;
        long long *ncor;
        int numcorr, pcor,length,p2;
        void resizeCorrelator(int k);
    public:
        scalar *t,*f,*fav,*fsqav,*tav;
        int npcorr,npcorrmax,nexp;
        TCorrelatorVector(int numcorrin = 35, int pcorin = 16);
        ~TCorrelatorVector();
        int resizecount=0, addcount=0;
        void add(arr3 w, int k=1);
        void init();
        void evaluate();
        void toaverage();
        void clear();
        void save(std::string savename);
        void read(std::string readname);
};

class TCorrelatorDiffusionVector{
    private:
        arr3 *aa;
        scalar *cor;
        long long *ncor;
        int numcorr, pcor,length,p2;
        void resizeCorrelator(int k);
    public:
        scalar *t,*f,*fav,*fsqav,*tav;
        int npcorr,npcorrmax,nexp;
        TCorrelatorDiffusionVector(int numcorrin = 35, int pcorin = 16);
        ~TCorrelatorDiffusionVector();
        TCorrelatorDiffusionVector(const TCorrelatorDiffusionVector &obj);
        int resizecount=0, addcount=0;
        void add(arr3 w, int k=1);
        void init();
        void evaluate();
        void toaverage();
        int getlength();
        void clear();
        void save(std::string savename);
        void save_out(std::string savename);
        void save_ffea(FILE *fout);
        void read_ffea(FILE *fout);
        void read(std::string readname);
};

struct cossindata{
    arr3 cos;
    arr3 sin;
};

class TCorrelatorSq{
    private:
        cossindata *aa;
        scalar *cor;
        long long *ncor;
        int numcorr, pcor,length,p2;
        void resizeCorrelator(int k);
    public:
        scalar *t,*f,*fav,*fsqav,*tav;
        int npcorr,npcorrmax,nexp;
        TCorrelatorSq(int numcorrin = 35, int pcorin = 16);
        ~TCorrelatorSq();
        int resizecount=0, addcount=0;
        void add(cossindata w, int k=1);
        void init();
        void evaluate();
        void toaverage();
        void clear();
        void save(std::string savename);
        void read(std::string readname);
};
/*
class TCorrelatorStress{
    private:
        TCorrelator *c0;
        TCorrelator *c1;
        TCorrelator *c2;
        TCorrelator *c3;
        TCorrelator *c4;
        TCorrelator *c5;
        TCorrelator *c6;
    public:
        double t[ntpoints],f[ntpoints];
        TCorrelatorStress(int numcorrin = 35,int pcorin=16);
        ~TCorrelatorStress();
        void add(matrix3 w, int k=1);
        void evaluate();
        void clear();
        void save_ffea(FILE *fout);
        void read_ffea(FILE *fout);
        void save(std::string savename);
        void save_out(FILE *fout);
        int read(std::string readname);
};

*/
class Fmm_blob{
public:
    Fmm_blob();
    ~Fmm_blob();
    arr3 pos;
    int PBC_Count[3];
    void init();
};

#endif // MT_CORRELATOR_H
