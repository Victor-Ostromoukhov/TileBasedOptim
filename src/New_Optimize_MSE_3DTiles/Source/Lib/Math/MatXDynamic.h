//
// Created by Loïs Paulin on 2019-05-20.

#ifndef STATTOOLSSAMPLING_MATXDYNAMIC_H
#define STATTOOLSSAMPLING_MATXDYNAMIC_H

#include <vector>
#include "VecXDynamic.h"
#include <string>

class MatXDynamic {

public:
    int width, height;
    std::vector<double> data;
    explicit MatXDynamic(int w = 0, int h = 0);
    explicit MatXDynamic(const std::vector<double>& dat, int w = 0, int h = 0);

    static MatXDynamic Identity(int d);

    double& operator()(int i, int j);
    double operator()(int i, int j) const;
    double& access(int i, int j);
    double access(int i, int j) const;

    MatXDynamic& operator+=(const MatXDynamic& m);
    MatXDynamic& operator-=(const MatXDynamic& m);
    explicit operator double() const;

    MatXDynamic transpose() const;

    bool LUPDecompose(MatXDynamic& L, MatXDynamic& U, std::vector<int>& P) const;
    double det() const;
    double highestEigenPair(VecXDynamic& v) const;
    MatXDynamic inverse() const;
    VecXDynamic solve(const VecXDynamic& b) const;

    friend MatXDynamic operator* (const MatXDynamic& m1, const MatXDynamic& m2);
    friend std::ostream& operator<<(std::ostream& out, const MatXDynamic& m);
};

MatXDynamic operator+ (const MatXDynamic& m1, const MatXDynamic& m2);
MatXDynamic operator- (const MatXDynamic& m1, const MatXDynamic& m2);
MatXDynamic operator* (const MatXDynamic& m1, const MatXDynamic& m2);
MatXDynamic operator* (double m1, const MatXDynamic& m2);
MatXDynamic operator* (const MatXDynamic& m1, double m2);

MatXDynamic transpose(const MatXDynamic& m);

template<class VECTYPE>
MatXDynamic operator* (const VECTYPE& v1, const MatXDynamic& m2){

    if (v1.dim() != m2.width && m2.height != 1){
        std::cerr << "Error: Incompatible matrix sizes (1, " + std::to_string(v1.dim()) + ") x ("
                     + std::to_string(m2.width) + ", " + std::to_string(m2.height) + ") for * operator" << std::endl;
        exit(1);
    }

    MatXDynamic m(v1.dim(), v1.dim());

    for (int j = 0; j < m.height; ++j){
        for (int i = 0; i < m.width; ++i){
            m(i, j) += v1[j] * m2(i,0);
        }
    }
    return m;
}

template<class VECTYPE>
VECTYPE operator* (const MatXDynamic& m1, const VECTYPE& v2){

    if (v2.dim() != m1.width){
        std::cerr << "Error: Incompatible matrix sizes ("
                     + std::to_string(m1.width) + ", " + std::to_string(m1.height) + ") x ("
                     + std::to_string(v2.dim()) + ") for * operator" << std::endl;
        exit(1);
    }

    VecXDynamic v(m1.height);

    for (int i = 0; i < m1.height; ++i){
        for (int k = 0; k < m1.width; ++k){
            v[i] += v2[k] * m1(k,i);
        }
    }
    return v;
}

template<class VECTYPE>
MatXDynamic transpose(const VECTYPE& v){
    MatXDynamic m(v.dim(), 1);
    for (int i = 0; i < v.dim(); ++i){
        m(i,0) = v[i];
    }
    return m;
}


#endif //STATTOOLSSAMPLING_MATXDYNAMIC_H
