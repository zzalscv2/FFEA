#include "mt_correlator.h"
TCorrelator::TCorrelator(int numcorrin, int pcorin){
    numcorr = numcorrin;
    pcor = pcorin;
    p2 = pcor/2;
    length = (numcorr+1)*p2;
    aa = new scalar[(numcorr+1)*(pcor+3)];
    cor = new scalar[(numcorr+1)*(pcor+3)];
    ncor = new long long[((numcorr+1)*(pcor+3))];
    t = new scalar[length];
    f = new scalar[length];
    fav = new scalar[length];
    fsqav = new scalar[length];
    tav = new scalar[length];
    init();
}

TCorrelator::~TCorrelator(){
    numcorr = 0;
    pcor = 0;
    p2 = 0;
    length = 0;
    delete aa;
    aa = 0;
    delete cor;
    cor = 0;
    delete ncor;
    ncor = 0;
    delete t;
    t = 0;
    delete f;
    f = 0;
    delete fav;
    fav = 0;
    delete fsqav;
    fsqav = 0;
    delete tav;
    tav = 0;
}

int TCorrelator::get_numcorr(){
    return numcorr;
}

void TCorrelator::resizeCorrelator(int k){
    int oldlength;

    numcorr = k;
    oldlength = length;
    length = (numcorr + 1)*p2;


    scalar *new_aa = new scalar[(numcorr+1)*(pcor+3)];
    scalar *new_cor = new scalar[(numcorr+1)*(pcor+3)];
    long long *new_ncor = new long long[(numcorr+1)*(pcor+3)];


    for(int i=0;i<numcorr;i++){
        for(int j=0;j<pcor+3;j++){
            new_aa[i*(pcor+3) + j]=aa[i*(pcor+3) + j];
            new_cor[i*(pcor+3) + j]=cor[i*(pcor+3) + j];
            new_ncor[i*(pcor+3) + j]=ncor[i*(pcor+3) + j];
        }
    }

    for(int j=0;j<pcor+3;j++){
    new_aa[numcorr*(pcor+3) + j]= -2E10;
    new_cor[numcorr*(pcor+3) + j]=0;
    new_ncor[numcorr*(pcor+3) + j]=0;
    }
    delete[] aa;
    aa = new_aa;
    delete[] cor;
    cor = new_cor;
    delete[] ncor;
    ncor = new_ncor;

    scalar *newt = new scalar[length];
    scalar *newf = new scalar[length];
    scalar *newfav = new scalar[length];
    scalar *newfsqav = new scalar[length];
    scalar *newtav = new scalar[length];



    if(oldlength<=length)
    for(int i=0; i<oldlength;i++){
        newt[i] = t[i];
        newf[i] = f[i];
        newfav[i] = fav[i];
        newfsqav[i] = fsqav[i];
        newtav[i] = tav[i];
    }
    for(int i=oldlength; i<length;i++){
        newt[i] = 0;
        newf[i] = 0;
        newfav[i] = 0;
        newfsqav[i] = 0;
        newtav[i] = 0;
    }

    delete[] t;
    t = newt;
    delete[] f;
    f = newf;
    delete[] fav;
    fav = newfav;
    delete[] fsqav;
    fsqav = newfsqav;
    delete[] tav;
    tav = newtav;


}


void TCorrelator::init(){

    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){
            aa[i*(pcor+3) + j]= -2E10;
            cor[i*(pcor+3) + j]= 0;
            ncor[i*(pcor+3) + j]= 0;
        }
    }


    for(int i=0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
        fav[i] = 0;
        fsqav[i] = 0;
        tav[i] = 0;
    }

    npcorr = 0;
    npcorrmax = 0;
    nexp = 0;
}

void TCorrelator::add(scalar w,int k){

    if(k>numcorr){
        resizecount++;
        resizeCorrelator(k);

    }

    for(int i = pcor;i>1;i--){
        aa[k*(pcor+3) + i]= aa[k*(pcor+3) + i-1];
    }

    aa[k*(pcor+3) + 1] = w;

    aa[k*(pcor+3) + pcor + 2] = aa[k*(pcor+3) +  pcor + 1];
    aa[k*(pcor+3) + pcor + 1] = w;

    for(int m =1;m<=p2;m++){
        if(aa[k*(pcor+3) + p2+m]>-1E10){
            cor[k*(pcor+3) +  m] = cor[k*(pcor+3) +  m] + aa[k*(pcor+3) +  1]* aa[k*(pcor+3) + p2+m];
            ncor[k*(pcor+3)+m] +=1;
        }
    }

    if(k==1){
        for(int m=0;m<p2;m++){
            if(aa[k*(pcor+3) +  m+1]>-1E10){
                cor[m] = cor[m]+aa[k*(pcor+3) +  1]*aa[k*(pcor+3) + 1+m];
                ncor[m]++;
            }
        }
    }


    if(aa[k*(pcor+3) +  pcor+2]>-1E10){
        add((aa[k*(pcor+3) + pcor+2]+aa[k*(pcor+3) + pcor+1])/2, k+1);
        aa[k*(pcor+3) + pcor+2]=-2E10;
        aa[k*(pcor+3) + pcor+1]=-2E10;
    }
}

void TCorrelator::evaluate(){
    long im;
    im = -1;

    for(int m =0;m<p2;m++){
        if(ncor[m]>0){
            im++;
            t[im]=m;
            f[im] = cor[m] / ncor[m];
        }
    }
    for(int k=1;k<=numcorr;k++){
        for(int m=1;m<=p2;m++){
            if(ncor[k*(pcor+3) + m]>0){
                im++;
                t[im] = (p2 - 1 + m) * pow(2, k - 1);
                f[im] = cor[k*(pcor+3) + m] / ncor[k*(pcor+3) + m];
            }
        }
    }
    npcorr = im;
}


void TCorrelator::toaverage(){
    for(int i=0;i<=npcorr;i++){
        fav[i] = fav[i] + f[i];
        fsqav[i] = fsqav[i] + f[i]*f[i];
        tav[i] = t[i];
    }
    if(npcorr>npcorrmax){
        npcorrmax = npcorr;
    }
    nexp++;
}

void TCorrelator::clear(){
    for(int i=0;i<=pcor+2;i++){
        for(int j=0;j<=numcorr;j++){
            aa[j*(pcor+3) +  i] = -2E10;
            cor[j*(pcor+3) +  i] = 0;
            ncor[j*(pcor+3) +  i] = 0;
        }
    }

    for(int i = 0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
    }
    npcorr = 0;
}

void TCorrelator::save(std::string savename){

    std::ofstream fsave;

    fsave.precision(14);

    fsave.open(savename);

    fsave<<numcorr<<"\t"<<pcor<<"\t"<<length<<endl;
    fsave<<npcorr<<"\t"<<npcorrmax<<"\t"<<nexp<<endl;
    for(int i = 0;i<=numcorr;i++){
        for(int j = 0;j<pcor+3;j++){
            fsave<<std::scientific<<aa[i*(pcor+3) + j]<<"\t"<<cor[i*(pcor+3) + j]<<"\t"<<ncor[i*(pcor+3) + j]<<endl;
        }
    }

    fsave.close();
}
/*
void TCorrelator::save_out(std::string savename){
    std::ofstream fsave;

    fsave.precision(14);

    fsave.open(savename);

    for(int i = 0;i<=length;i++){
            fsave<<std::scientific<<t[i]<<"\t"<<f[i]<<"\t"<<fav[i]<<"\t"<<fsqav[i]<<"\t"<<tav[i]<<endl;
    }

    fsave.close();
}*/

void TCorrelator::save_out(FILE *fout){
    for(int i = 0;i<=length;i++){
        fprintf(fout,"%f\t%f\t%f\t%f\t%f\n",t[i],f[i],fav[i],fsqav[i],tav[i]);
    }

}



void TCorrelator::save_ffea(FILE *fout){
    
    fprintf(fout,"%d\t%d\t%d\n",numcorr,pcor,length);
    fprintf(fout,"%d\t%d\t%d\n",npcorr,npcorrmax,nexp);
    for(int i = 0;i<=numcorr;i++){
        for(int j = 0;j<pcor+3;j++){
            fprintf(fout,"%.15e\t%.15e\t%lld\n",aa[i*(pcor+3) + j],cor[i*(pcor+3) + j],ncor[i*(pcor+3) + j]);
        }
    }
}

void TCorrelator::read_ffea(FILE *fout){
   
    if(fscanf(fout,"%d\t%d\t%d\n",&numcorr,&pcor,&length)!=3){cout<<"Reading Problem"<<endl;};
    if(fscanf(fout,"%d\t%d\t%d\n",&npcorr,&npcorrmax,&nexp)!=3){cout<<"Reading Problem"<<endl;};
    for(int i = 0;i<=numcorr;i++){
        for(int j = 0;j<pcor+3;j++){
            if(fscanf(fout,"%lf\t%lf\t%lld\n",&aa[i*(pcor+3) + j],&cor[i*(pcor+3) + j],&ncor[i*(pcor+3) + j])!=3){cout<<"Reading Problem"<<endl;};
        }
    }
}

void TCorrelator::read(std::string readname){
    cout<<"initiating the read"<<endl;

    std::string storage;

    std::ifstream fread;
    fread.open(readname);


    std::getline(fread,storage,'\t');
    numcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    pcor = stoi(storage);
    std::getline(fread,storage,'\n');
    length = stoi(storage);

    init();

    std::getline(fread,storage,'\t');
    npcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    npcorrmax = stoi(storage);
    std::getline(fread,storage,'\n');
    nexp = stoi(storage);


    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){

            storage.clear();


            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j] = stod(storage);
            std::getline(fread,storage,'\t');
            cor[i*(pcor+3) + j] = stod(storage);
            std::getline(fread,storage,'\n');
            ncor[i*(pcor+3) + j]= stoll(storage);
        }
    }
    fread.close();

}

TCrossCorrelator::TCrossCorrelator(int numcorrin, int pcorin){
    numcorr = numcorrin;
    pcor = pcorin;
    p2 = pcor/2;
    length = (numcorr+1)*p2;
    aa = new scalar[(numcorr+1)*(pcor+3)];
    aa2 = new scalar[(numcorr+1)*(pcor+3)];
    cor = new scalar[(numcorr+1)*(pcor+3)];
    ncor = new long long[((numcorr+1)*(pcor+3))];
    t = new scalar[length];
    f = new scalar[length];
    fav = new scalar[length];
    fsqav = new scalar[length];
    tav = new scalar[length];
    init();
}

TCrossCorrelator::~TCrossCorrelator(){
    numcorr = 0;
    pcor = 0;
    p2 = 0;
    length = 0;
    delete aa;
    delete cor;
    delete ncor;
    delete t;
    delete f;
    delete fav;
    delete fsqav;
    delete tav;
}

void TCrossCorrelator::resizeCorrelator(int k){
    int oldlength;

    numcorr = k;
    oldlength = length;
    length = (numcorr + 1)*p2;

    scalar *new_aa = new scalar[(numcorr+1)*(pcor+3)];
    scalar *new_aa2 = new scalar[(numcorr+1)*(pcor+3)];
    scalar *new_cor = new scalar[(numcorr+1)*(pcor+3)];
    long long *new_ncor = new long long[(numcorr+1)*(pcor+3)];


    for(int i=0;i<numcorr;i++){
        for(int j=0;j<pcor+3;j++){
            new_aa[i*(pcor+3) + j]=aa[i*(pcor+3) + j];
            new_aa2[i*(pcor+3) + j]=aa2[i*(pcor+3) + j];
            new_cor[i*(pcor+3) + j]=cor[i*(pcor+3) + j];
            new_ncor[i*(pcor+3) + j]=ncor[i*(pcor+3) + j];
        }
    }

    for(int j=0;j<pcor+3;j++){
    new_aa[numcorr*(pcor+3) + j]= -2E10;
    new_aa2[numcorr*(pcor+3) + j]= -2E10;
    new_cor[numcorr*(pcor+3) + j]=0;
    new_ncor[numcorr*(pcor+3) + j]=0;
    }
    delete[] aa;
    aa = new_aa;
    delete[] aa2;
    aa2 = new_aa2;
    delete[] cor;
    cor = new_cor;
    delete[] ncor;
    ncor = new_ncor;

    scalar *newt = new scalar[length];
    scalar *newf = new scalar[length];
    scalar *newfav = new scalar[length];
    scalar *newfsqav = new scalar[length];
    scalar *newtav = new scalar[length];



    if(oldlength<=length)
    for(int i=0; i<oldlength;i++){
        newt[i] = t[i];
        newf[i] = f[i];
        newfav[i] = fav[i];
        newfsqav[i] = fsqav[i];
        newtav[i] = tav[i];
    }

    for(int i=oldlength; i<length;i++){
        newt[i] = 0;
        newf[i] = 0;
        newfav[i] = 0;
        newfsqav[i] = 0;
        newtav[i] = 0;
    }

    delete[] t;
    t = newt;
    delete[] f;
    f = newf;
    delete[] fav;
    fav = newfav;
    delete[] fsqav;
    fsqav = newfsqav;
    delete[] tav;
    tav = newtav;


}

void TCrossCorrelator::init(){

    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){
            aa[i*(pcor+3) + j]= -2E10;
            aa2[i*(pcor+3) + j]= -2E10;
            cor[i*(pcor+3) + j]= 0;
            ncor[i*(pcor+3) + j]= 0;
        }
    }


    for(int i=0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
        fav[i] = 0;
        fsqav[i] = 0;
        tav[i] = 0;
    }

    npcorr = 0;
    npcorrmax = 0;
    nexp = 0;
}

void TCrossCorrelator::add(scalar w,scalar w2,int k){

    if(k>numcorr){
        resizecount++;
        resizeCorrelator(k);

    }

    for(int i = pcor;i>1;i--){
        aa[k*(pcor+3) + i]= aa[k*(pcor+3) + i-1];
        aa2[k*(pcor+3) + i]= aa2[k*(pcor+3) + i-1];
    }

    aa[k*(pcor+3) + 1] = w;
    aa2[k*(pcor+3) + 1] = w2;

    aa[k*(pcor+3) + pcor + 2] = aa[k*(pcor+3) +  pcor + 1];
    aa2[k*(pcor+3) + pcor + 2] = aa2[k*(pcor+3) +  pcor + 1];
    aa[k*(pcor+3) + pcor + 1] = w;
    aa2[k*(pcor+3) + pcor + 1] = w2;

    for(int m =1;m<=p2;m++){
        if(aa[k*(pcor+3) + p2+m]>-1E10){
            cor[k*(pcor+3) +  m] = cor[k*(pcor+3) +  m] + aa2[k*(pcor+3) +  1]* aa[k*(pcor+3) + p2+m];
            ncor[k*(pcor+3)+m] +=1;
        }
    }

    if(k==1){
        for(int m=0;m<p2;m++){
            if(aa[k*(pcor+3) +  m+1]>-1E10){
                cor[m] = cor[m]+aa[k*(pcor+3) +  1]*aa2[k*(pcor+3) + 1+m];
                ncor[m]++;
            }
        }
    }


    if(aa[k*(pcor+3) +  pcor+2]>-1E10){
        add((aa[k*(pcor+3) + pcor+2]+aa[k*(pcor+3) + pcor+1])/2, (aa2[k*(pcor+3) + pcor+2]+aa2[k*(pcor+3) + pcor+1])/2, k+1);
        aa[k*(pcor+3) + pcor+2]=-2E10;
        aa[k*(pcor+3) + pcor+1]=-2E10;
    }
}

void TCrossCorrelator::evaluate(){
    long im;
    im = -1;

    for(int m =0;m<p2;m++){
        if(ncor[m]>0){
            im++;
            t[im]=m;
            f[im] = cor[m] / ncor[m];
        }
    }
    for(int k=1;k<=numcorr;k++){
        for(int m=1;m<=p2;m++){
            if(ncor[k*(pcor+3) + m]>0){
                im++;
                t[im] = (p2 - 1 + m) * pow(2, k - 1);
                f[im] = cor[k*(pcor+3) + m] / ncor[k*(pcor+3) + m];
            }
        }
    }
    npcorr = im;
}

void TCrossCorrelator::toaverage(){
    for(int i=0;i<=npcorr;i++){
        fav[i] = fav[i] + f[i];
        fsqav[i] = fsqav[i] + f[i]*f[i];
        tav[i] = t[i];
    }
    if(npcorr>npcorrmax){
        npcorrmax = npcorr;
    }
    nexp++;
}

void TCrossCorrelator::clear(){
    for(int i=0;i<=pcor+2;i++){
        for(int j=0;j<=numcorr;j++){
            aa[j*(pcor+3) +  i] = -2E10;
            aa2[j*(pcor+3) +  i] = -2E10;
            cor[j*(pcor+3) +  i] = 0;
            ncor[j*(pcor+3) +  i] = 0;
        }
    }

    for(int i = 0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
    }
    npcorr = 0;
}



void TCrossCorrelator::save(std::string savename){

    std::ofstream fsave;

    fsave.precision(14);

    fsave.open(savename);

    fsave<<numcorr<<"\t"<<pcor<<"\t"<<length<<"\n";
    fsave<<npcorr<<"\t"<<npcorrmax<<"\t"<<nexp<<"\n";
    for(int i = 0;i<=numcorr;i++){
        for(int j = 0;j<pcor+3;j++){
            fsave<<std::scientific<<aa[i*(pcor+3) + j]<<"\t"<<aa2[i*(pcor+3) + j]<<"\t"<<cor[i*(pcor+3) + j]<<"\t"<<ncor[i*(pcor+3) + j]<<endl;
        }
    }



    fsave.close();
}

void TCrossCorrelator::read(std::string readname){
    cout<<"initiating the read"<<endl;

    std::string storage;

    std::ifstream fread;
    fread.open(readname);


    std::getline(fread,storage,'\t');
    numcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    pcor = stoi(storage);
    std::getline(fread,storage,'\n');
    length = stoi(storage);

    init();

    std::getline(fread,storage,'\t');
    npcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    npcorrmax = stoi(storage);
    std::getline(fread,storage,'\n');
    nexp = stoi(storage);


    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){

            storage.clear();


            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j] = stod(storage);
            std::getline(fread,storage,'\t');
            aa2[i*(pcor+3) + j] = stod(storage);
            std::getline(fread,storage,'\t');
            cor[i*(pcor+3) + j] = stod(storage);
            std::getline(fread,storage,'\n');
            ncor[i*(pcor+3) + j]= stoll(storage);
        }
    }
    fread.close();

}

TCorrelatorVector::TCorrelatorVector(int numcorrin, int pcorin){
    numcorr = numcorrin;
    pcor = pcorin;
    p2 = pcor/2;
    length = (numcorr+1)*p2;
    aa = new arr3[(numcorr+1)*(pcor+3)];
    cor = new scalar[(numcorr+1)*(pcor+3)];
    ncor = new long long[((numcorr+1)*(pcor+3))];
    t = new scalar[length];
    f = new scalar[length];
    fav = new scalar[length];
    fsqav = new scalar[length];
    tav = new scalar[length];
    init();
}

TCorrelatorVector::~TCorrelatorVector(){
    numcorr = 0;
    pcor = 0;
    p2 = 0;
    length = 0;
    delete aa;
    delete cor;
    delete ncor;
    delete t;
    delete f;
    delete fav;
    delete fsqav;
    delete tav;
}

void TCorrelatorVector::init(){

    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){
            for(int m=0; m<3;m++){
                aa[i*(pcor+3) + j][m]= -2E10;
            }
            cor[i*(pcor+3) + j]= 0;
            ncor[i*(pcor+3) + j]= 0;
        }
    }



    for(int i=0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
        fav[i] = 0;
        fsqav[i] = 0;
        tav[i] = 0;
    }

    npcorr = 0;
    npcorrmax = 0;
    nexp = 0;
}

void TCorrelatorVector::resizeCorrelator(int k){
    int oldlength;

    numcorr = k;
    oldlength = length;
    length = (numcorr + 1)*p2;

    arr3 *new_aa = new arr3[(numcorr+1)*(pcor+3)];
    scalar *new_cor = new scalar[(numcorr+1)*(pcor+3)];
    long long *new_ncor = new long long[(numcorr+1)*(pcor+3)];


    for(int i=0;i<numcorr;i++){
        for(int j=0;j<pcor+3;j++){
            for(int m=0;m<3;m++){
                new_aa[i*(pcor+3) + j][m]=aa[i*(pcor+3) + j][m];
            }
            new_cor[i*(pcor+3) + j]=cor[i*(pcor+3) + j];
            new_ncor[i*(pcor+3) + j]=ncor[i*(pcor+3) + j];
        }
    }

    for(int j=0;j<pcor+3;j++){
    for(int m=0;m<3;m++){
        new_aa[numcorr*(pcor+3) + j][m]= -2E10;
    }
    new_cor[numcorr*(pcor+3) + j]=0;
    new_ncor[numcorr*(pcor+3) + j]=0;
    }
    delete[] aa;
    aa = new_aa;
    delete[] cor;
    cor = new_cor;
    delete[] ncor;
    ncor = new_ncor;

    scalar *newt = new scalar[length];
    scalar *newf = new scalar[length];
    scalar *newfav = new scalar[length];
    scalar *newfsqav = new scalar[length];
    scalar *newtav = new scalar[length];



    if(oldlength<=length)
    for(int i=0; i<oldlength;i++){
        newt[i] = t[i];
        newf[i] = f[i];
        newfav[i] = fav[i];
        newfsqav[i] = fsqav[i];
        newtav[i] = tav[i];
    }
    for(int i=oldlength; i<length;i++){
        newt[i] = 0;
        newf[i] = 0;
        newfav[i] = 0;
        newfsqav[i] = 0;
        newtav[i] = 0;
    }

    delete[] t;
    t = newt;
    delete[] f;
    f = newf;
    delete[] fav;
    fav = newfav;
    delete[] fsqav;
    fsqav = newfsqav;
    delete[] tav;
    tav = newtav;


}

void TCorrelatorVector::add(arr3 w,int k){

    arr3 f;

    if(k>numcorr){
        resizecount++;
        resizeCorrelator(k);
    }

    for(int i = pcor;i>1;i--){
        for(int m=0;m<3;m++){
            aa[k*(pcor+3) + i][m]= aa[k*(pcor+3) + i-1][m];
        }
    }
    for(int m=0;m<3;m++){
    aa[k*(pcor+3) + 1][m] = w[m];

    aa[k*(pcor+3) + pcor + 2][m] = aa[k*(pcor+3) +  pcor + 1][m];
    aa[k*(pcor+3) + pcor + 1][m] = w[m];
    }
    for(int m =1;m<=p2;m++){
        if(aa[k*(pcor+3) + p2+m][0]>-1E10){
            cor[k*(pcor+3) +  m] = cor[k*(pcor+3) +  m] + arr3arr3DotProduct<scalar,arr3>(aa[k*(pcor+3) +  1], aa[k*(pcor+3) + p2+m]);
            ncor[k*(pcor+3)+m] +=1;
        }
    }

    if(k==1){
        for(int m=0;m<p2;m++){
            if(aa[k*(pcor+3) +  m+1][0]>-1E10){
                cor[m] = cor[m]+arr3arr3DotProduct<scalar,arr3>(aa[k*(pcor+3) +  1],aa[k*(pcor+3) + 1+m]);
                ncor[m]++;
            }
        }
    }


    if(aa[k*(pcor+3) +  pcor+2][0]>-1E10){

        arr3arr3Add<scalar,arr3>(aa[k*(pcor+3) + pcor+2],aa[k*(pcor+3) + pcor+1],f);
        arr3Resize<scalar,arr3>(0.5,f);
        add(f, k+1);
        aa[k*(pcor+3) + pcor+2][0]=-2E10;
        aa[k*(pcor+3) + pcor+1][0]=-2E10;
    }
}

void TCorrelatorVector::save(std::string savename){

    std::ofstream fsave;

    fsave.precision(14);

    fsave.open(savename);

    fsave<<numcorr<<"\t"<<pcor<<"\t"<<length<<endl;
    fsave<<npcorr<<"\t"<<npcorrmax<<"\t"<<nexp<<endl;
    for(int i = 0;i<=numcorr;i++){
        for(int j = 0;j<pcor+3;j++){
            fsave<<std::scientific<<aa[i*(pcor+3) + j][0]<<"\t"<<aa[i*(pcor+3) + j][1]<<"\t"<<aa[i*(pcor+3) + j][2]<<"\t"<<cor[i*(pcor+3) + j]<<"\t"<<ncor[i*(pcor+3) + j]<<endl;
        }
    }



    fsave.close();
}

void TCorrelatorVector::read(std::string readname){
    cout<<"initiating the read"<<endl;

    std::string storage;

    std::ifstream fread;
    fread.open(readname);


    std::getline(fread,storage,'\t');
    numcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    pcor = stoi(storage);
    std::getline(fread,storage,'\n');
    length = stoi(storage);

    init();

    std::getline(fread,storage,'\t');
    npcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    npcorrmax = stoi(storage);
    std::getline(fread,storage,'\n');
    nexp = stoi(storage);


    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){

            storage.clear();


            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j][0] = stod(storage);
            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j][1] = stod(storage);
            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j][2] = stod(storage);
            std::getline(fread,storage,'\t');
            cor[i*(pcor+3) + j] = stod(storage);
            std::getline(fread,storage,'\n');
            ncor[i*(pcor+3) + j]= stoll(storage);
        }
    }
    fread.close();

}

void TCorrelatorVector::evaluate(){
    long im;
    im = -1;

    for(int m =0;m<p2;m++){
        if(ncor[m]>0){
            im++;
            t[im]=m;
            f[im] = cor[m] / ncor[m];
        }
    }
    for(int k=1;k<=numcorr;k++){
        for(int m=1;m<=p2;m++){
            if(ncor[k*(pcor+3) + m]>0){
                im++;
                t[im] = (p2 - 1 + m) * pow(2, k - 1);
                f[im] = cor[k*(pcor+3) + m] / ncor[k*(pcor+3) + m];
            }
        }
    }
    npcorr = im;
}

void TCorrelatorVector::toaverage(){
    for(int i=0;i<=npcorr;i++){
        fav[i] = fav[i] + f[i];
        fsqav[i] = fsqav[i] + f[i]*f[i];
        tav[i] = t[i];
    }
    if(npcorr>npcorrmax){
        npcorrmax = npcorr;
    }
    nexp++;
}

void TCorrelatorVector::clear(){
    for(int i=0;i<=pcor+2;i++){
        for(int j=0;j<=numcorr;j++){
            for(int m=0;m<3;m++){
                aa[j*(pcor+3) +  i][m] = -2E10;
            }
            cor[j*(pcor+3) +  i] = 0;
            ncor[j*(pcor+3) +  i] = 0;
        }
    }

    for(int i = 0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
    }
    npcorr = 0;
}

TCorrelatorDiffusionVector::TCorrelatorDiffusionVector(int numcorrin, int pcorin){
    numcorr = numcorrin;
    pcor = pcorin;
    p2 = pcor/2;
    length = (numcorr+1)*p2;
    aa = new arr3[(numcorr+1)*(pcor+3)];
    cor = new scalar[(numcorr+1)*(pcor+3)];
    ncor = new long long[((numcorr+1)*(pcor+3))];
    t = new scalar[length];
    f = new scalar[length];
    fav = new scalar[length];
    fsqav = new scalar[length];
    tav = new scalar[length];
    init();
}

TCorrelatorDiffusionVector::~TCorrelatorDiffusionVector(){
    numcorr = 0;
    pcor = 0;
    p2 = 0;
    length = 0;
    delete[] aa;
    delete[] cor;
    delete[] ncor;
    delete[] t;
    delete[] f;
    delete[] fav;
    delete[] fsqav;
    delete[] tav;
}

TCorrelatorDiffusionVector::TCorrelatorDiffusionVector(const TCorrelatorDiffusionVector &obj){
    numcorr = obj.numcorr;
    pcor = obj.pcor;
    p2 = obj.pcor/2;
    length = (numcorr+1)*p2;
    aa = new arr3[(numcorr+1)*(pcor+3)];
    cor = new scalar[(numcorr+1)*(pcor+3)];
    ncor = new long long[((numcorr+1)*(pcor+3))];
    t = new scalar[length];
    f = new scalar[length];
    fav = new scalar[length];
    fsqav = new scalar[length];
    tav = new scalar[length];

    for(int i=0;i<numcorr;i++){
        for(int j=0;j<pcor+3;j++){
            for(int m=0;m<3;m++){
                aa[i*(pcor+3) + j][m]=obj.aa[i*(pcor+3) + j][m];
            }
            cor[i*(pcor+3) + j]=obj.cor[i*(pcor+3) + j];
            ncor[i*(pcor+3) + j]=obj.ncor[i*(pcor+3) + j];
        }
    }

    for(int i=0; i<length;i++){
        t[i] = obj.t[i];
        f[i] = obj.f[i];
        fav[i] = obj.fav[i];
        fsqav[i] = obj.fsqav[i];
        tav[i] = obj.tav[i];
    }

}

int TCorrelatorDiffusionVector::getlength(){
return length;
}


void TCorrelatorDiffusionVector::init(){

    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){
            for(int m=0; m<3;m++){
                aa[i*(pcor+3) + j][m]= -2E10;
            }
            cor[i*(pcor+3) + j]= 0;
            ncor[i*(pcor+3) + j]= 0;
        }
    }

    for(int i=0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
        fav[i] = 0;
        fsqav[i] = 0;
        tav[i] = 0;
    }

    npcorr = 0;
    npcorrmax = 0;
    nexp = 0;
}

void TCorrelatorDiffusionVector::resizeCorrelator(int k){
    int oldlength;

    numcorr = k;
    oldlength = length;
    length = (numcorr + 1)*p2;

    arr3 *new_aa = new arr3[(numcorr+1)*(pcor+3)];
    scalar *new_cor = new scalar[(numcorr+1)*(pcor+3)];
    long long *new_ncor = new long long[(numcorr+1)*(pcor+3)];


    for(int i=0;i<numcorr;i++){
        for(int j=0;j<pcor+3;j++){
            for(int m=0;m<3;m++){
                new_aa[i*(pcor+3) + j][m]=aa[i*(pcor+3) + j][m];
            }
            new_cor[i*(pcor+3) + j]=cor[i*(pcor+3) + j];
            new_ncor[i*(pcor+3) + j]=ncor[i*(pcor+3) + j];
        }
    }

    for(int j=0;j<pcor+3;j++){
    for(int m=0;m<3;m++){
        new_aa[numcorr*(pcor+3) + j][m]= -2E10;
    }
    new_cor[numcorr*(pcor+3) + j]=0;
    new_ncor[numcorr*(pcor+3) + j]=0;
    }
    delete[] aa;
    aa = new_aa;
    delete[] cor;
    cor = new_cor;
    delete[] ncor;
    ncor = new_ncor;

    scalar *newt = new scalar[length];
    scalar *newf = new scalar[length];
    scalar *newfav = new scalar[length];
    scalar *newfsqav = new scalar[length];
    scalar *newtav = new scalar[length];



    if(oldlength<=length)
    for(int i=0; i<oldlength;i++){
        newt[i] = t[i];
        newf[i] = f[i];
        newfav[i] = fav[i];
        newfsqav[i] = fsqav[i];
        newtav[i] = tav[i];
    }
    for(int i=oldlength; i<length;i++){
        newt[i] = 0;
        newf[i] = 0;
        newfav[i] = 0;
        newfsqav[i] = 0;
        newtav[i] = 0;
    }

    delete[] t;
    t = newt;
    delete[] f;
    f = newf;
    delete[] fav;
    fav = newfav;
    delete[] fsqav;
    fsqav = newfsqav;
    delete[] tav;
    tav = newtav;

}

void TCorrelatorDiffusionVector::add(arr3 w,int k){

    //cout<<"added x for w is "<<w[0]<<" added y for w is "<<w[1]<<" added z for w is "<<w[2]<<endl;
    arr3 f;

    if(k>numcorr){
        resizecount++;
        resizeCorrelator(k);
    }

    for(int i = pcor;i>1;i--){
        for(int m=0;m<3;m++){
            aa[k*(pcor+3) + i][m]= aa[k*(pcor+3) + i-1][m];
        }
    }
    for(int m=0;m<3;m++){
    aa[k*(pcor+3) + 1][m] = w[m];

    aa[k*(pcor+3) + pcor + 2][m] = aa[k*(pcor+3) +  pcor + 1][m];
    aa[k*(pcor+3) + pcor + 1][m] = w[m];
    }
    for(int m =1;m<=p2;m++){
        if(aa[k*(pcor+3) + p2+m][0]>-1E10){
            cor[k*(pcor+3) +  m] = cor[k*(pcor+3) +  m] + pow(arr3arr3Distance<scalar,arr3>(aa[k*(pcor+3) +  1], aa[k*(pcor+3) + p2+m]),2);
            ncor[k*(pcor+3)+m] +=1;
        }
    }

    if(k==1){
        for(int m=0;m<p2;m++){
            if(aa[k*(pcor+3) +  m+1][0]>-1E10){
                cor[m] = cor[m]+pow(arr3arr3Distance<scalar,arr3>(aa[k*(pcor+3) +  1],aa[k*(pcor+3) + 1+m]),2);
                ncor[m]++;
            }
        }
    }


    if(aa[k*(pcor+3) +  pcor+2][0]>-1E10){
        arr3arr3Add<scalar,arr3>(aa[k*(pcor+3) + pcor+2],aa[k*(pcor+3) + pcor+1],f);
        arr3Resize<scalar,arr3>(0.5,f);
        add(f, k+1);
        aa[k*(pcor+3) + pcor+2][0]=-2E10;
        aa[k*(pcor+3) + pcor+1][0]=-2E10;
    }
}

void TCorrelatorDiffusionVector::save(std::string savename){

    std::ofstream fsave;

    fsave.precision(14);

    fsave.open(savename);

    fsave<<numcorr<<"\t"<<pcor<<"\t"<<length<<endl;
    fsave<<npcorr<<"\t"<<npcorrmax<<"\t"<<nexp<<endl;
    for(int i = 0;i<=numcorr;i++){
        for(int j = 0;j<pcor+3;j++){
            fsave<<std::scientific<<aa[i*(pcor+3) + j][0]<<"\t"<<aa[i*(pcor+3) + j][1]<<"\t"<<aa[i*(pcor+3) + j][2]<<"\t"<<cor[i*(pcor+3) + j]<<"\t"<<ncor[i*(pcor+3) + j]<<endl;
        }
    }

    fsave.close();
}

void TCorrelatorDiffusionVector::save_out(std::string savename){
    std::ofstream fsave;

    fsave.precision(14);

    fsave.open(savename);

    for(int i = 0;i<=length;i++){
            fsave<<std::scientific<<t[i]<<"\t"<<f[i]<<"\t"<<fav[i]<<"\t"<<fsqav[i]<<"\t"<<tav[i]<<endl;
    }

    fsave.close();
}


void TCorrelatorDiffusionVector::read(std::string readname){
    cout<<"initiating the read"<<endl;

    std::string storage;

    std::ifstream fread;
    fread.open(readname);


    std::getline(fread,storage,'\t');
    numcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    pcor = stoi(storage);
    std::getline(fread,storage,'\n');
    length = stoi(storage);

    init();

    std::getline(fread,storage,'\t');
    npcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    npcorrmax = stoi(storage);
    std::getline(fread,storage,'\n');
    nexp = stoi(storage);


    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){

            storage.clear();


            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j][0] = stod(storage);
            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j][1] = stod(storage);
            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j][2] = stod(storage);
            std::getline(fread,storage,'\t');
            cor[i*(pcor+3) + j] = stod(storage);
            std::getline(fread,storage,'\n');
            ncor[i*(pcor+3) + j]= stoll(storage);
        }
    }
    fread.close();

}

void TCorrelatorDiffusionVector::save_ffea(FILE *fout){
    //cout<<"IT'S SAVING TIME"<<endl;
        
        fprintf(fout,"%d\t%d\t%d\n",numcorr,pcor,length);
        
        fprintf(fout,"%d\t%d\t%d\n",npcorr,npcorrmax,nexp);
        
        for(int i = 0;i<=numcorr;i++){
            for(int j = 0;j<pcor+3;j++){
                //cout<<std::scientific<<aa[i*(pcor+3) + j][0]<<"\t"<<aa[i*(pcor+3) + j][1]<<"\t"<<aa[i*(pcor+3) + j][2]<<"\t"<<cor[i*(pcor+3) + j]<<"\t"<<ncor[i*(pcor+3) + j]<<endl;
                fprintf(fout,"%.15e\t%.15e\t%.15e\t%.15e\t%lld\n",aa[i*(pcor+3) + j][0],aa[i*(pcor+3) + j][1],aa[i*(pcor+3) + j][2],cor[i*(pcor+3) + j],ncor[i*(pcor+3) + j]);
            }
        }
        //cout<<"THE PLANET IS SAVED"<<endl;
}

void TCorrelatorDiffusionVector::read_ffea(FILE *fout){
    cout<<"started read_ffea"<<endl;
    if (fscanf(fout,"%d\t%d\t%d\n",&numcorr,&pcor,&length)!=3){cout<<"Reading PROBLEM";};
    cout<<"read first line"<<endl;
    if (fscanf(fout,"%d\t%d\t%d\n",&npcorr,&npcorrmax,&nexp)!=3){cout<<"Reading PROBLEM";};
    cout<<"read second line"<<endl;
    for(int i = 0;i<=numcorr;i++){
        for(int j = 0;j<pcor+3;j++){
            if (fscanf(fout,"%lf\t%lf\t%lf\t%lf\t%lld\n",&aa[i*(pcor+3) + j][0],&aa[i*(pcor+3) + j][1],&aa[i*(pcor+3) + j][2],&cor[i*(pcor+3) + j],&ncor[i*(pcor+3) + j])!=5){cout<<"SAVING PROBLEM";};
        }
    }
    cout<<"Read whole correlator!"<<endl;
}

void TCorrelatorDiffusionVector::evaluate(){
    long im;
    double delta;
    im = -1;

    for(int m =0;m<p2;m++){
        if(ncor[m]>0){
            im++;
            t[im]=m;
            f[im] = cor[m] / ncor[m];
        }
    }
    for(int k=1;k<=numcorr;k++){
        for(int m=1;m<=p2;m++){
            if(ncor[k*(pcor+3) + m]>0){
                im++;
                t[im] = (p2 - 1 + m) * pow(2, k - 1);
                f[im] = cor[k*(pcor+3) + m] / ncor[k*(pcor+3) + m];

                if(m==1&&k>1){
                    delta = f[im - 2] + (t[im] - t[im - 2]) / (t[im - 1] - t[im - 2]) * (f[im - 1] - f[im - 2]) - f[im];
                }
                if(k>1){
                    f[im]+=delta;
                }
            }
        }
    }
    npcorr = im;
}

void TCorrelatorDiffusionVector::toaverage(){
    for(int i=0;i<=npcorr;i++){
        fav[i] = fav[i] + f[i];
        fsqav[i] = fsqav[i] + f[i]*f[i];
        tav[i] = t[i];
    }
    if(npcorr>npcorrmax){
        npcorrmax = npcorr;
    }
    nexp++;
}

void TCorrelatorDiffusionVector::clear(){
    for(int i=0;i<=pcor+2;i++){
        for(int j=0;j<=numcorr;j++){
            for(int m=0;m<3;m++){
                aa[j*(pcor+3) +  i][m] = -2E10;
            }
            cor[j*(pcor+3) +  i] = 0;
            ncor[j*(pcor+3) +  i] = 0;
        }
    }

    for(int i = 0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
    }
    npcorr = 0;
}

TCorrelatorSq::TCorrelatorSq(int numcorrin, int pcorin){
    numcorr = numcorrin;
    pcor = pcorin;
    p2 = pcor/2;
    length = (numcorr+1)*p2;
    aa = new cossindata[(numcorr+1)*(pcor+3)];
    cor = new scalar[(numcorr+1)*(pcor+3)];
    ncor = new long long[((numcorr+1)*(pcor+3))];
    t = new scalar[length];
    f = new scalar[length];
    fav = new scalar[length];
    fsqav = new scalar[length];
    tav = new scalar[length];
    init();
}

TCorrelatorSq::~TCorrelatorSq(){
    numcorr = 0;
    pcor = 0;
    p2 = 0;
    length = 0;
    delete aa;
    delete cor;
    delete ncor;
    delete t;
    delete f;
    delete fav;
    delete fsqav;
    delete tav;
}

void TCorrelatorSq::init(){

    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){
            for(int m=0; m<3;m++){
                aa[i*(pcor+3) + j].cos[m]= -2E10;
                aa[i*(pcor+3) + j].sin[m]= -2E10;
            }
            cor[i*(pcor+3) + j]= 0;
            ncor[i*(pcor+3) + j]= 0;
        }
    }

    for(int i=0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
        fav[i] = 0;
        fsqav[i] = 0;
        tav[i] = 0;
    }

    npcorr = 0;
    npcorrmax = 0;
    nexp = 0;
}

void TCorrelatorSq::resizeCorrelator(int k){
    int oldlength;

    numcorr = k;
    oldlength = length;
    length = (numcorr + 1)*p2;

    cossindata *new_aa = new cossindata[(numcorr+1)*(pcor+3)];
    scalar *new_cor = new scalar[(numcorr+1)*(pcor+3)];
    long long *new_ncor = new long long[(numcorr+1)*(pcor+3)];


    for(int i=0;i<numcorr;i++){
        for(int j=0;j<pcor+3;j++){
            for(int m=0;m<3;m++){
                new_aa[i*(pcor+3) + j].cos[m]=aa[i*(pcor+3) + j].cos[m];
                new_aa[i*(pcor+3) + j].sin[m]=aa[i*(pcor+3) + j].sin[m];
            }
            new_cor[i*(pcor+3) + j]=cor[i*(pcor+3) + j];
            new_ncor[i*(pcor+3) + j]=ncor[i*(pcor+3) + j];
        }
    }

    for(int j=0;j<pcor+3;j++){
    for(int m=0;m<3;m++){
        new_aa[numcorr*(pcor+3) + j].cos[m]= -2E10;
        new_aa[numcorr*(pcor+3) + j].sin[m]= -2E10;
    }
    new_cor[numcorr*(pcor+3) + j]=0;
    new_ncor[numcorr*(pcor+3) + j]=0;
    }
    delete[] aa;
    aa = new_aa;
    delete[] cor;
    cor = new_cor;
    delete[] ncor;
    ncor = new_ncor;

    scalar *newt = new scalar[length];
    scalar *newf = new scalar[length];
    scalar *newfav = new scalar[length];
    scalar *newfsqav = new scalar[length];
    scalar *newtav = new scalar[length];


    if(oldlength<=length)
    for(int i=0; i<oldlength;i++){
        newt[i] = t[i];
        newf[i] = f[i];
        newfav[i] = fav[i];
        newfsqav[i] = fsqav[i];
        newtav[i] = tav[i];
    }
    for(int i=oldlength; i<length;i++){
        newt[i] = 0;
        newf[i] = 0;
        newfav[i] = 0;
        newfsqav[i] = 0;
        newtav[i] = 0;
    }

    delete[] t;
    t = newt;
    delete[] f;
    f = newf;
    delete[] fav;
    fav = newfav;
    delete[] fsqav;
    fsqav = newfsqav;
    delete[] tav;
    tav = newtav;


}

void TCorrelatorSq::add(cossindata w,int k){

    cossindata f;

    if(k>numcorr){
        resizecount++;
        resizeCorrelator(k);
    }

    for(int i = pcor;i>1;i--){
        for(int m=0;m<3;m++){
            aa[k*(pcor+3) + i].cos[m]= aa[k*(pcor+3) + i-1].cos[m];
            aa[k*(pcor+3) + i].sin[m]= aa[k*(pcor+3) + i-1].sin[m];
        }
    }
    for(int m=0;m<3;m++){
    aa[k*(pcor+3) + 1].cos[m] = w.cos[m];
    aa[k*(pcor+3) + 1].sin[m] = w.sin[m];

    aa[k*(pcor+3) + pcor + 2].cos[m] = aa[k*(pcor+3) +  pcor + 1].cos[m];
    aa[k*(pcor+3) + pcor + 2].sin[m] = aa[k*(pcor+3) +  pcor + 1].sin[m];
    aa[k*(pcor+3) + pcor + 1].cos[m] = w.cos[m];
    aa[k*(pcor+3) + pcor + 1].sin[m] = w.sin[m];
    }
    for(int m =1;m<=p2;m++){
        if(aa[k*(pcor+3) + p2+m].cos[0]>-1E10){
            cor[k*(pcor+3) +  m] = cor[k*(pcor+3) +  m] + (arr3arr3DotProduct<scalar,arr3>(aa[k*(pcor+3) +  1].cos, aa[k*(pcor+3) + p2+m].cos)+arr3arr3DotProduct<scalar,arr3>(aa[k*(pcor+3) +  1].sin, aa[k*(pcor+3) + p2+m].sin))/3;
            ncor[k*(pcor+3)+m] +=1;
        }
    }

    if(k==1){
        for(int m=0;m<p2;m++){
            if(aa[k*(pcor+3) +  m+1].cos[0]>-1E10){
                cor[m] = cor[m]+(arr3arr3DotProduct<scalar,arr3>(aa[k*(pcor+3) +  1].cos,aa[k*(pcor+3) + 1+m].cos)+arr3arr3DotProduct<scalar,arr3>(aa[k*(pcor+3) +  1].sin,aa[k*(pcor+3) + 1+m].sin))/3;
                ncor[m]++;
            }
        }
    }


    if(aa[k*(pcor+3) +  pcor+2].cos[0]>-1E10){
        arr3arr3Add<scalar,arr3>(aa[k*(pcor+3) + pcor+2].cos,aa[k*(pcor+3) + pcor+1].cos,f.cos);
        arr3arr3Add<scalar,arr3>(aa[k*(pcor+3) + pcor+2].sin,aa[k*(pcor+3) + pcor+1].sin,f.sin);
        arr3Resize<scalar,arr3>(0.5,f.cos);
        arr3Resize<scalar,arr3>(0.5,f.sin);
        add(f, k+1);
        aa[k*(pcor+3) + pcor+2].cos[0]=-2E10;
        aa[k*(pcor+3) + pcor+2].sin[0]=-2E10;
        aa[k*(pcor+3) + pcor+1].cos[0]=-2E10;
        aa[k*(pcor+3) + pcor+1].sin[0]=-2E10;
    }
}

void TCorrelatorSq::save(std::string savename){

    std::ofstream fsave;

    fsave.precision(14);

    fsave.open(savename);

    fsave<<numcorr<<"\t"<<pcor<<"\t"<<length<<endl;
    fsave<<npcorr<<"\t"<<npcorrmax<<"\t"<<nexp<<endl;
    for(int i = 0;i<=numcorr;i++){
        for(int j = 0;j<pcor+3;j++){
            fsave<<std::scientific<<aa[i*(pcor+3) + j].cos[0]<<"\t"<<aa[i*(pcor+3) + j].cos[1]<<"\t"<<aa[i*(pcor+3) + j].cos[2]<<"\t"<<aa[i*(pcor+3) + j].sin[0]<<"\t"<<aa[i*(pcor+3) + j].sin[1]<<"\t"<<aa[i*(pcor+3) + j].sin[2]<<"\t"<<cor[i*(pcor+3) + j]<<"\t"<<ncor[i*(pcor+3) + j]<<endl;
        }
    }



    fsave.close();
}

void TCorrelatorSq::read(std::string readname){
    cout<<"initiating the read"<<endl;

    std::string storage;

    std::ifstream fread;
    fread.open(readname);


    std::getline(fread,storage,'\t');
    numcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    pcor = stoi(storage);
    std::getline(fread,storage,'\n');
    length = stoi(storage);

    init();

    std::getline(fread,storage,'\t');
    npcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    npcorrmax = stoi(storage);
    std::getline(fread,storage,'\n');
    nexp = stoi(storage);


    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){

            storage.clear();


            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j].cos[0] = stod(storage);
            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j].cos[1] = stod(storage);
            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j].cos[2] = stod(storage);
            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j].sin[0] = stod(storage);
            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j].sin[1] = stod(storage);
            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j].sin[2] = stod(storage);
            std::getline(fread,storage,'\t');
            cor[i*(pcor+3) + j] = stod(storage);
            std::getline(fread,storage,'\n');
            ncor[i*(pcor+3) + j]= stoll(storage);
        }
    }
    fread.close();

}

void TCorrelatorSq::evaluate(){
    long im;
    im = -1;

    for(int m =0;m<p2;m++){
        if(ncor[m]>0){
            im++;
            t[im]=m;
            f[im] = cor[m] / ncor[m];
        }
    }
    for(int k=1;k<=numcorr;k++){
        for(int m=1;m<=p2;m++){
            if(ncor[k*(pcor+3) + m]>0){
                im++;
                t[im] = (p2 - 1 + m) * pow(2, k - 1);
                f[im] = cor[k*(pcor+3) + m] / ncor[k*(pcor+3) + m];
            }
        }
    }
    npcorr = im;
}

void TCorrelatorSq::toaverage(){
    for(int i=0;i<=npcorr;i++){
        fav[i] = fav[i] + f[i];
        fsqav[i] = fsqav[i] + f[i]*f[i];
        tav[i] = t[i];
    }
    if(npcorr>npcorrmax){
        npcorrmax = npcorr;
    }
    nexp++;
}

void TCorrelatorSq::clear(){
    for(int i=0;i<=pcor+2;i++){
        for(int j=0;j<=numcorr;j++){
                aa[j*(pcor+3) +  i].cos[0] = -2E10;
            cor[j*(pcor+3) +  i] = 0;
            ncor[j*(pcor+3) +  i] = 0;
        }
    }

    for(int i = 0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
    }
    npcorr = 0;
}

TCorrelatorStress::TCorrelatorStress(int numcorrin,int pcorin){
    TCorrelator c0(numcorrin,pcorin);
    TCorrelator c1(numcorrin,pcorin);
    TCorrelator c2(numcorrin,pcorin);
    TCorrelator c3(numcorrin,pcorin);
    TCorrelator c4(numcorrin,pcorin);
    TCorrelator c5(numcorrin,pcorin);
    TCorrelator c6(numcorrin,pcorin);
    
}

TCorrelatorStress::~TCorrelatorStress(){

}

void TCorrelatorStress::add(matrix3 w,int k){
    c0.add(w[0][1]);
    c1.add(w[1][2]);
    c2.add(w[2][0]);
    c3.add(w[1][1]-w[2][2]);
    c4.add(w[0][0]-w[1][1]);
    c5.add(w[0][0]-w[2][2]);
    c6.add(w[0][0]+w[1][1]+w[2][2]);
}

void TCorrelatorStress::evaluate(){
    //for(int i = 0;i<5;i++){
        c0.evaluate();
        c1.evaluate();
        c2.evaluate();
        c3.evaluate();
        c4.evaluate();
        c5.evaluate();
        c6.evaluate();
    //}
/*
    for(int i=0;i<ntpoints;i++){
        f[i] = 0;

            f[i] = f[i]+0.2 * c0.f[i];
            f[i] = f[i]+0.2 * c1.f[i];
            f[i] = f[i]+0.2 * c2.f[i];
            f[i] = f[i]+0.2 * c3.f[i];
            f[i] = f[i]+0.2 * c4.f[i];
            t[i] = c0.t[i];
     }
*/
}



void TCorrelatorStress::clear(){
    c0.clear();
    c1.clear();
    c2.clear();
    c3.clear();
    c4.clear();
    c5.clear();
    c6.clear();

    for(int i = 0;i<ntpoints;i++){
        t[i] = 0;
        f[i] = 0;
    }
}
/*
this won't work as is, and I don't use it

void TCorrelatorStress::save(std::string savename){

        c0.save(savename+std::to_string(0)+".txt");
        c1.save(savename+std::to_string(1)+".txt");
        c2.save(savename+std::to_string(2)+".txt");
        c3.save(savename+std::to_string(3)+".txt");
        c4.save(savename+std::to_string(4)+".txt");
        c5.save(savename+std::to_string(5)+".txt");
        c6.save(savename+std::to_string(6)+".txt");

  std::ofstream fsave;

  fsave.open(savename+"_base.txt");

  fsave<<c0.get_numcorr()<<"\t"<<c1.get_numcorr()<<"\t"<<c2.get_numcorr()<<"\t"<<c3.get_numcorr()<<"\t"<<c4.get_numcorr()<<endl;
  for(int i=0;i<ntpoints;i++){
      fsave<<t[i]<<"\t"<<f[i]<<endl;
  }
    fsave.close();

}
*/
/*
int TCorrelatorStress::read(std::string readname){
    std::string storage;

    int ncor0,ncor1,ncor2, ncor3, ncor4;

    std::ifstream fread;
    fread.open(readname+"_base.txt");
    std::getline(fread,storage,'\t');
    ncor0 = stoi(storage);
    std::getline(fread,storage,'\t');
    ncor1 = stoi(storage);
    std::getline(fread,storage,'\t');
    ncor2 = stoi(storage);
    std::getline(fread,storage,'\t');
    ncor3 = stoi(storage);
    std::getline(fread,storage,'\t');
    ncor4 = stoi(storage);

    if(!(ncor0==ncor1&&ncor1==ncor2&&ncor2==ncor3&&ncor3==ncor4)){
        cout<<"Correlator's not of equal size. Aborting..."<<endl;
        return 0;
    }

        c0.read(readname+std::to_string(0)+".txt");
        c1.read(readname+std::to_string(1)+".txt");
        c2.read(readname+std::to_string(2)+".txt");
        c3.read(readname+std::to_string(3)+".txt");
        c4.read(readname+std::to_string(4)+".txt");
}
*/

void TCorrelatorStress::save_ffea(FILE *fout){
   
   c0.save_ffea(fout);
   c1.save_ffea(fout);
   c2.save_ffea(fout);
   c3.save_ffea(fout);
   c4.save_ffea(fout);
   c5.save_ffea(fout);
   c6.save_ffea(fout);
   
}
    

void TCorrelatorStress::read_ffea(FILE *fout){
   
    c0.read_ffea(fout);
    c1.read_ffea(fout);
    c2.read_ffea(fout);
    c3.read_ffea(fout);
    c4.read_ffea(fout);
    c5.read_ffea(fout);
    c6.read_ffea(fout);
       
}

void TCorrelatorStress::save_out(FILE *fout){
    c0.save_out(fout);
    c1.save_out(fout);
    c2.save_out(fout);
    c3.save_out(fout);
    c4.save_out(fout);
    c5.save_out(fout);
    c6.save_out(fout);
}

Fmm_blob::Fmm_blob(){
    init();
}

void Fmm_blob::init(){
    for(int i=0;i<3;i++){
        pos[i]=0;
        PBC_Count[i]=0;
    }
   /* for(int i=0;i<6;i++){
        FFt[i]=0;
    }*/

}


Fmm_blob::~Fmm_blob(){
    for(int i=0;i<3;i++){
        pos[i]=0;
        PBC_Count[i]=0;
    }
    /*
    for(int i=0;i<6;i++){
        FFt[i]=0;
    }*/

}

TCorrelatorDiffusion::TCorrelatorDiffusion(int numcorrin, int pcorin){
    numcorr = numcorrin;
    pcor = pcorin;
    p2 = pcor/2;
    length = (numcorr+1)*p2;
    aa = new scalar[(numcorr+1)*(pcor+3)];
    cor = new scalar[(numcorr+1)*(pcor+3)];
    ncor = new long long[((numcorr+1)*(pcor+3))];
    t = new scalar[length];
    f = new scalar[length];
    fav = new scalar[length];
    fsqav = new scalar[length];
    tav = new scalar[length];
    init();
}

TCorrelatorDiffusion::~TCorrelatorDiffusion(){
    numcorr = 0;
    pcor = 0;
    p2 = 0;
    length = 0;
    delete[] aa;
    delete[] cor;
    delete[] ncor;
    delete[] t;
    delete[] f;
    delete[] fav;
    delete[] fsqav;
    delete[] tav;
}

TCorrelatorDiffusion::TCorrelatorDiffusion(const TCorrelatorDiffusion &obj){
    numcorr = obj.numcorr;
    pcor = obj.pcor;
    p2 = obj.pcor/2;
    length = (numcorr+1)*p2;
    aa = new scalar[(numcorr+1)*(pcor+3)];
    cor = new scalar[(numcorr+1)*(pcor+3)];
    ncor = new long long[((numcorr+1)*(pcor+3))];
    t = new scalar[length];
    f = new scalar[length];
    fav = new scalar[length];
    fsqav = new scalar[length];
    tav = new scalar[length];

    for(int i=0;i<numcorr;i++){
        for(int j=0;j<pcor+3;j++){
            aa[i*(pcor+3) + j]=obj.aa[i*(pcor+3) + j];
            cor[i*(pcor+3) + j]=obj.cor[i*(pcor+3) + j];
            ncor[i*(pcor+3) + j]=obj.ncor[i*(pcor+3) + j];
        }
    }

    for(int i=0; i<length;i++){
        t[i] = obj.t[i];
        f[i] = obj.f[i];
        fav[i] = obj.fav[i];
        fsqav[i] = obj.fsqav[i];
        tav[i] = obj.tav[i];
    }

}

void TCorrelatorDiffusion::resizeCorrelator(int k){
    int oldlength;

    numcorr = k;
    oldlength = length;
    length = (numcorr + 1)*p2;


    scalar *new_aa = new scalar[(numcorr+1)*(pcor+3)];
    scalar *new_cor = new scalar[(numcorr+1)*(pcor+3)];
    long long *new_ncor = new long long[(numcorr+1)*(pcor+3)];


    for(int i=0;i<numcorr;i++){
        for(int j=0;j<pcor+3;j++){
            new_aa[i*(pcor+3) + j]=aa[i*(pcor+3) + j];
            new_cor[i*(pcor+3) + j]=cor[i*(pcor+3) + j];
            new_ncor[i*(pcor+3) + j]=ncor[i*(pcor+3) + j];
        }
    }

    for(int j=0;j<pcor+3;j++){
    new_aa[numcorr*(pcor+3) + j]= -2E10;
    new_cor[numcorr*(pcor+3) + j]=0;
    new_ncor[numcorr*(pcor+3) + j]=0;
    }
    delete[] aa;
    aa = new_aa;
    delete[] cor;
    cor = new_cor;
    delete[] ncor;
    ncor = new_ncor;

    scalar *newt = new scalar[length];
    scalar *newf = new scalar[length];
    scalar *newfav = new scalar[length];
    scalar *newfsqav = new scalar[length];
    scalar *newtav = new scalar[length];



    if(oldlength<=length)
    for(int i=0; i<oldlength;i++){
        newt[i] = t[i];
        newf[i] = f[i];
        newfav[i] = fav[i];
        newfsqav[i] = fsqav[i];
        newtav[i] = tav[i];
    }
    for(int i=oldlength; i<length;i++){
        newt[i] = 0;
        newf[i] = 0;
        newfav[i] = 0;
        newfsqav[i] = 0;
        newtav[i] = 0;
    }

    delete[] t;
    t = newt;
    delete[] f;
    f = newf;
    delete[] fav;
    fav = newfav;
    delete[] fsqav;
    fsqav = newfsqav;
    delete[] tav;
    tav = newtav;


}


void TCorrelatorDiffusion::init(){

    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){
            aa[i*(pcor+3) + j]= -2E10;
            cor[i*(pcor+3) + j]= 0;
            ncor[i*(pcor+3) + j]= 0;
        }
    }


    for(int i=0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
        fav[i] = 0;
        fsqav[i] = 0;
        tav[i] = 0;
    }

    npcorr = 0;
    npcorrmax = 0;
    nexp = 0;
}

void TCorrelatorDiffusion::add(scalar w,int k){

    if(k>numcorr){
        resizecount++;
        resizeCorrelator(k);

    }

    for(int i = pcor;i>1;i--){
        aa[k*(pcor+3) + i]= aa[k*(pcor+3) + i-1];
    }

    aa[k*(pcor+3) + 1] = w;

    aa[k*(pcor+3) + pcor + 2] = aa[k*(pcor+3) +  pcor + 1];
    aa[k*(pcor+3) + pcor + 1] = w;

    for(int m =1;m<=p2;m++){
        if(aa[k*(pcor+3) + p2+m]>-1E10){
            cor[k*(pcor+3) +  m] = cor[k*(pcor+3) +  m] + (aa[k*(pcor+3) +  1]- aa[k*(pcor+3) + p2+m])*(aa[k*(pcor+3) +  1]- aa[k*(pcor+3) + p2+m]);
            ncor[k*(pcor+3)+m] +=1;
        }
    }

    if(k==1){
        for(int m=0;m<p2;m++){
            if(aa[k*(pcor+3) +  m+1]>-1E10){
                cor[m] = cor[m]+(aa[k*(pcor+3) +  1]-aa[k*(pcor+3) + 1+m])*(aa[k*(pcor+3) +  1]-aa[k*(pcor+3) + 1+m]);
                ncor[m]++;
            }
        }
    }


    if(aa[k*(pcor+3) +  pcor+2]>-1E10){
        add((aa[k*(pcor+3) + pcor+2]+aa[k*(pcor+3) + pcor+1])/2, k+1);
        aa[k*(pcor+3) + pcor+2]=-2E10;
        aa[k*(pcor+3) + pcor+1]=-2E10;
    }
}

void TCorrelatorDiffusion::evaluate(){
    long im;
    double delta;
    im = -1;

    for(int m =0;m<p2;m++){
        if(ncor[m]>0){
            im++;
            t[im]=m;
            f[im] = cor[m] / ncor[m];
        }
    }
    for(int k=1;k<=numcorr;k++){
        for(int m=1;m<=p2;m++){
            if(ncor[k*(pcor+3) + m]>0){
                im++;
                t[im] = (p2 - 1 + m) * pow(2, k - 1);
                f[im] = cor[k*(pcor+3) + m] / ncor[k*(pcor+3) + m];

                if(m==1&&k>1){
                    delta = f[im - 2] + (t[im] - t[im - 2]) / (t[im - 1] - t[im - 2]) * (f[im - 1] - f[im - 2]) - f[im];
                }
                if(k>1){
                    f[im]+=delta;
                }
            }
        }
    }
    npcorr = im;
}


void TCorrelatorDiffusion::toaverage(){
    for(int i=0;i<=npcorr;i++){
        fav[i] = fav[i] + f[i];
        fsqav[i] = fsqav[i] + f[i]*f[i];
        tav[i] = t[i];
    }
    if(npcorr>npcorrmax){
        npcorrmax = npcorr;
    }
    nexp++;
}

void TCorrelatorDiffusion::clear(){
    for(int i=0;i<=pcor+2;i++){
        for(int j=0;j<=numcorr;j++){
            aa[j*(pcor+3) +  i] = -2E10;
            cor[j*(pcor+3) +  i] = 0;
            ncor[j*(pcor+3) +  i] = 0;
        }
    }

    for(int i = 0;i<length;i++){
        t[i] = 0;
        f[i] = 0;
    }
    npcorr = 0;
}

void TCorrelatorDiffusion::save(std::string savename){

    std::ofstream fsave;

    fsave.precision(14);

    fsave.open(savename);

    fsave<<numcorr<<"\t"<<pcor<<"\t"<<length<<endl;
    fsave<<npcorr<<"\t"<<npcorrmax<<"\t"<<nexp<<endl;
    for(int i = 0;i<=numcorr;i++){
        for(int j = 0;j<pcor+3;j++){
            fsave<<std::scientific<<aa[i*(pcor+3) + j]<<"\t"<<cor[i*(pcor+3) + j]<<"\t"<<ncor[i*(pcor+3) + j]<<endl;
        }
    }

    fsave.close();
}

void TCorrelatorDiffusion::save_out(std::string savename){
    std::ofstream fsave;

    fsave.precision(14);

    fsave.open(savename);

    for(int i = 0;i<=length;i++){
            fsave<<std::scientific<<t[i]<<"\t"<<f[i]<<"\t"<<fav[i]<<"\t"<<fsqav[i]<<"\t"<<tav[i]<<endl;
    }

    fsave.close();
}

void TCorrelatorDiffusion::save_ffea(FILE *fout){
   
    fprintf(fout,"%d\t%d\t%d\n",numcorr,pcor,length);
    fprintf(fout,"%d\t%d\t%d\n",npcorr,npcorrmax,nexp);
    for(int i = 0;i<=numcorr;i++){
        for(int j = 0;j<pcor+3;j++){
            fprintf(fout,"%.15e\t%.15e\t%lld\n",aa[i*(pcor+3) + j],cor[i*(pcor+3) + j],ncor[i*(pcor+3) + j]);
        }
    }
       
}


void TCorrelatorDiffusion::read_ffea(FILE *fout){
   
    if(fscanf(fout,"%d\t%d\t%d\n",&numcorr,&pcor,&length)!=3){cout<<"Reading Problem"<<endl;};
    if(fscanf(fout,"%d\t%d\t%d\n",&npcorr,&npcorrmax,&nexp)!=3){cout<<"Reading Problem"<<endl;};
    for(int i = 0;i<=numcorr;i++){
        for(int j = 0;j<pcor+3;j++){
            if(fscanf(fout,"%lf\t%lf\t%lld\n",&aa[i*(pcor+3) + j],&cor[i*(pcor+3) + j],&ncor[i*(pcor+3) + j])!=3){cout<<"Reading Problem"<<endl;};
        }
    }
       
}



void TCorrelatorDiffusion::read(std::string readname){
    cout<<"initiating the read"<<endl;

    std::string storage;

    std::ifstream fread;
    fread.open(readname);


    std::getline(fread,storage,'\t');
    numcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    pcor = stoi(storage);
    std::getline(fread,storage,'\n');
    length = stoi(storage);

    init();

    std::getline(fread,storage,'\t');
    npcorr = stoi(storage);
    std::getline(fread,storage,'\t');
    npcorrmax = stoi(storage);
    std::getline(fread,storage,'\n');
    nexp = stoi(storage);


    for(int i=0;i<numcorr+1;i++){
        for(int j=0;j<pcor+3;j++){

            storage.clear();


            std::getline(fread,storage,'\t');
            aa[i*(pcor+3) + j] = stod(storage);
            std::getline(fread,storage,'\t');
            cor[i*(pcor+3) + j] = stod(storage);
            std::getline(fread,storage,'\n');
            ncor[i*(pcor+3) + j]= stoll(storage);
        }
    }
    fread.close();

}
