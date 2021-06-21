#ifndef _inc_lattice_reduction_wrapper
#define _inc_lattice_reduction_wrapper

#include "bkzlibwrapper.hpp"

#define gsupdate 0x00
#define nogsupdate 0x00

//Subroutine for LLL reduction

//return Euclidean norm
template<typename T=double> T LengthOf(vec_ZZ& b) {
    //approximated length
    ZZ ip;
    InnerProduct(ip, b,b);
    if (ip==0) return 0;

    //To prevent overflow
    double lip = to_double(log(ip));
    T r = lip;
    return exp(r*0.5);  
}

#define define_normal
#include "LLLcore.cpp"
#undef define_normal

#define define_small
#include "LLLcore.cpp"
#undef define_small

template <typename T> vec_ZZ NearestPlane(vec_ZZ& vin,LatticeBasis<T>& B) {
    
    LatticeBasis<mpfr_float> B2;
    B2 = B;

    std::vector<mpfr_float> gsrept;
    vec_ZZ v = vin;
    
    while (1) {
        GSrepresentation(gsrept,v,B2);
        int zeroflag = 0;
        for (int i=B.dim;i>=1;i--) {
            mpfr_float bias;
            bias = round(gsrept[i]);
            if (bias!=0) {
                gsrept[i] -= bias;
                for (int k=1;k<i;k++) gsrept[k] -= B2.gs.mu[i][k] * bias;
                v -= to_ZZ(bias.str(1000).c_str()) * B2.L[i-1];
                zeroflag = 1;
            }
        } 
        if (zeroflag == 0) break;
    }
    return v;
}

mat_ZZ BUBuffer;

template <typename CFLOAT,typename T> int SmallLocalLLL(LatticeBasis<T>& B,mat_ZZ*U,CFLOAT delta,int istart,int iend,int vl,int &swap,int& maxshift) {
    
    //Allocate mu-buffer 
    int cols = iend-istart+1;
    int rows = cols;

    SmallLatticeBasis<CFLOAT,CFLOAT> LB;   //local basis
    pbkzmatrix<int> LU;
    LB.L.SetDims(rows,cols);
    LU.SetDims(rows,cols);

    //diagonal
    B.updateGSBasis(istart,iend);

    //Find factor
    T mings = B.gs.c[1];
    for (int i=1;i<=cols;i++) {
        mings = min(mings,B.gs.c[i]);
    }    
    mings = mings / 100.0;
    if (mings < 1.0) mings = 1.0;

    for (int i=0;i<cols;i++) {
        for (int j=0;j<cols;j++) {
            if (j<i) {
                T prod = B.gs.c[j+istart] * B.gs.mu[i+istart][j+istart] / mings;
                LB.L[i][j] = (CFLOAT)prod;
            }
            if (j==i) {
                T prod = B.gs.c[i+istart] / mings;
                LB.L[i][i] = (CFLOAT)prod;  //unconditional jump error??
            }
            if (j>i) LB.L[i][j] = 0;
        }
    }
    for (int i=0;i<rows;i++) {
        for (int j=0;j<cols;j++) {
            LU[i][j] = 0;
            if (i==j) LU[i][i] = 1;
        }
    }
    int final_iend = local_LLL_small(LB,&LU,delta,1,rows,stnormal,vl,swap);
    
    //Reconstruct original matrix
    int n = B.L.NumCols();
   BUBuffer.SetDims(cols,n);

    for (int i=0;i<cols;i++) {
        clear(BUBuffer[i]);
        for (int j=0;j<cols;j++) {
            if (LU[i][j]!=0) RowTransformwrap(BUBuffer[i],  B.L[istart + j-1], to_ZZ(-LU[i][j]));
        }
    }
    for (int i=0;i<cols;i++) {
        B.L[istart+i-1] = BUBuffer[i];
    }    

    if (U!=0) {
        BUBuffer.SetDims(cols,(*U).NumCols());
        for (int i=0;i<cols;i++) {
            clear(BUBuffer[i]);
            for (int j=0;j<cols;j++) {
                if (LU[i][j]!=0) RowTransformwrap(BUBuffer[i],  (*U)[istart + j-1], to_ZZ(-LU[i][j]));
            }
        }
        for (int i=0;i<cols;i++) {
            (*U)[istart+i-1] = BUBuffer[i];
        }    
    }

    B.GScomputed = istart-1;
    B.updateGSBasis(1,istart+cols-1);


    for (int i=0;i<cols;i++) {
        bool flag_unitvector = true;
        for (int j=0;j<cols;j++) {
            if ((i!=j) && (LU[i][j]!=0)) {
                flag_unitvector = false;
                break;
            }
            if ((i==j) && (LU[i][j]!=1)) {
                flag_unitvector = false;
                break;
            }
        }
        maxshift = istart + i;
        if (flag_unitvector==false) {
            break;
        }
    }

    return final_iend;
}

template <typename T> int local_LLL(LatticeBasis<T>& B, mat_ZZ* U, double delta,T zerolimit=1) {
    int swapcount;
    return local_LLL(B,U,delta,1,B.dim,stnormal,VL0,swapcount,zerolimit);
}

template <typename T> void BigLLL(LatticeBasis<T>& B, mat_ZZ* U, double delta,int istart,int iend,int vl=0) {
    //For lattices with large elements
    
    int s = 1;
    int swap=0;
    int llsize = 160;
    double ss = gettimeofday_sec();
    for (int j=istart+s;j<=iend;j+=s) {

        B.updateGSBasis(1,j);
        
        int localswap=0;
        int bstart = istart;
        if (j > istart+s) {
            //cout << B.gs.c[j-1] << endl;
            T gap = log(B.gs.c[j-1]) / 2.30258509299;
            if (B.gs.c[j]<1) gap -= log(B.gs.c[j]) / 2.30258509299;  
            //cout << gap << endl;
            int llen = 1+(int)boost::lexical_cast<double>(gap);
            bstart = max(istart,j-llen);
        }
        
        local_LLL(B,U,0.5,bstart,j,stnormal,VL0,localswap);
        if (bstart>1) {
           for (int i=bstart;i<=j;i++) {
                SizeReduce(B,U,i);
           }
        }

        swap += localswap;
        if (vl>=1) {
            cout << "BigLLL: time=" << gettimeofday_sec() -ss << " " << bstart << ":" << j << "]/" << iend << " swap=" << swap << " \r";
            if (vl==1) cout.flush();
            if (vl>=2) cout << endl;;
        }            

        bstart = j - llsize + 1;
        if (bstart < istart) bstart = istart;
        int bend;
        int maxshift;
        while (1) {
            bend = min(bstart+llsize -1,j);
            SmallLocalLLL<long double>(B,U,(long double)0.999,bstart,bend,VL0,localswap,maxshift);
            if (vl>=1) {
                cout << "small LLL[" << bstart << ":" << bend << "] swap=" << swap << " maxshift=" << maxshift << "   \r";
                if (vl==1) cout.flush();
                if (vl>=2) cout << endl;;
            }
            swap += localswap;
            if (bstart>1) {
                for (int i=maxshift;i<=bend;i++) {
                    SizeReduce(B,U,i);
                }
            }
            if (localswap > 0) {
                if ((bstart==istart) && (bend==j)) break;
                bstart = maxshift - llsize/2;
                if (bstart < istart) bstart = istart;
            } else {
                if (bend == j) break;
                bstart += llsize;
            }
            bend = min(bstart+llsize -1,j);
        }
        if (vl>=1) {
            cout << "time=" << gettimeofday_sec() -ss << " Small-j=" << j << "/" << iend << " swap=" << swap <<  "    \r";;
            if (vl==1) cout.flush();
            if (vl>=2) cout << endl;;
        }            
        
    }
}

void BigLLL(mat_ZZ& L,mat_ZZ*U,double delta,int vl) {
    LatticeBasis<float50> B;
    B = L;
    //B.updateGSBasis();
    BigLLL(B,U,delta,1,B.dim,vl);
    L = B.L;
}

#endif