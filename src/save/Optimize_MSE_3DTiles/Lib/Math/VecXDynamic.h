
//
// Created by lpaulin on 25/03/19.
//

#ifndef SLICEDOPTIM_VECXDYNAMIC_H
#define SLICEDOPTIM_VECXDYNAMIC_H

#include <vector>
#include <iostream>
#include <cmath>
#include "VecX.h"

class VecXDynamic {
private:
    int N;
    std::vector<double> coefs;
public:

    int dim() const {
        return N;
    }

    inline VecXDynamic(int n, const std::vector<double>& values):N(n), coefs(values){
        if (values.size() != size_t(N)){
            std::cerr << "input vector must be of size " << N << " when creating VecXDynamic<" << N << std::endl;
            exit(1);
        }
    }

    const double* data() const{
        return coefs.data();
    }

    double* data(){
        return coefs.data();
    }

    explicit inline VecXDynamic(int n=0):N(n), coefs(n){
    }

    explicit inline operator double() const{
        if (N != 1){
            std::cerr << "Error: Can't convert Vector to double if dimension is higher than 1." << std::endl;
            std::cerr << "\tCurrent dimension is " << N << std::endl;
            exit(1);
        }
        return coefs[0];
    }

    inline VecXDynamic(const VecXDynamic& b):N(b.N), coefs(b.coefs){}

    template <int n>
    inline VecXDynamic(const VecX<n>& v):N(n), coefs(n){
        for (int i = 0; i < n; ++i){
            coefs[i] = v[i];
        }
    }

    inline VecXDynamic(double x, double y): N(2), coefs(2){
        coefs[0] = x;
        coefs[1] = y;
    }

    inline VecXDynamic(double x, double y, double z): N(3), coefs(3){
        coefs[0] = x;
        coefs[1] = y;
        coefs[2] = z;
    }

    inline VecXDynamic(double x, double y, double z, double w): N(4), coefs(4){
        coefs[0] = x;
        coefs[1] = y;
        coefs[2] = z;
        coefs[3] = w;
    }


    inline double norm() const{
        return std::sqrt(norm2());
    }
    inline double norm2() const{
        double n = 0.;
        for (const double& v : coefs){
            n += v*v;
        }
        return n;
    }
    inline void normalize(){
        double n = norm();
        for (double& v : coefs){
            v /= n;
        }
    }

    inline double& operator[] (int i){
        return coefs[i];
    }

    inline const double& operator[] (int i) const{
        return coefs[i];
    }

    inline VecXDynamic& operator+= (const VecXDynamic& b){
        if (b.N != N){
            std::cerr << "Operator +=  should be used between same sized VecXDynamic" << std::endl;
            exit(1);
        }
        for (size_t i = 0; i < coefs.size(); ++i){
            coefs[i] += b.coefs[i];
        }
        return *this;
    }
    inline VecXDynamic& operator-= (const VecXDynamic b){
        if (b.N != N){
            std::cerr << "Operator -= should be used between same sized VecXDynamic" << std::endl;
            exit(1);
        }
        for (size_t i = 0; i < coefs.size(); ++i){
            coefs[i] -= b.coefs[i];
        }
        return *this;

    }
    inline VecXDynamic& operator/= (double s) {
        for (size_t i = 0; i < coefs.size(); ++i) {
            coefs[i] /= s;
        }
        return *this;
    }
    inline VecXDynamic& operator*= (double s) {
        for (size_t i = 0; i < coefs.size(); ++i) {
            coefs[i] *= s;
        }
        return *this;
    }
    inline VecXDynamic& operator+= (double s) {
        for (size_t i = 0; i < coefs.size(); ++i) {
            coefs[i] += s;
        }
        return *this;
    }
    inline VecXDynamic& operator-= (double s) {
        for (size_t i = 0; i < coefs.size(); ++i) {
            coefs[i] -= s;
        }
        return *this;
    }

    inline VecXDynamic operator- (const VecXDynamic& b) const{
        if (b.N != N){
            std::cerr << "Operator - should be used between same sized VecXDynamic" << std::endl;
            exit(1);
        }
        VecXDynamic a = *this;
        for (size_t i = 0; i < coefs.size(); ++i) {
            a.coefs[i] -= b.coefs[i];
        }
        return a;
    }
    inline VecXDynamic operator+ (const VecXDynamic& b) const{
        if (b.N != N){
            std::cerr << "Operator + should be used between same sized VecXDynamic" << std::endl;
            exit(1);
        }
        VecXDynamic a = *this;
        for (size_t i = 0; i < coefs.size(); ++i) {
            a.coefs[i] += b.coefs[i];
        }
        return a;
    }
    inline VecXDynamic operator/ (const VecXDynamic& b) const{
        if (b.N != N){
            std::cerr << "Operator / should be used between same sized VecXDynamic" << std::endl;
            exit(1);
        }
        VecXDynamic a = *this;
        for (size_t i = 0; i < coefs.size(); ++i) {
            a.coefs[i] /= b.coefs[i];
        }
        return a;
    }
    inline VecXDynamic operator- (double s) const{
        VecXDynamic a = *this;
        for (size_t i = 0; i < coefs.size(); ++i) {
            a.coefs[i] -= s;
        }
        return a;
    }
    inline VecXDynamic operator+ (double s) const{
        VecXDynamic a = *this;
        for (size_t i = 0; i < coefs.size(); ++i) {
            a.coefs[i] += s;
        }
        return a;
    }
    inline VecXDynamic operator* (double s) const{
        VecXDynamic a = *this;
        for (size_t i = 0; i < coefs.size(); ++i) {
            a.coefs[i] *= s;
        }
        return a;
    }
    inline VecXDynamic operator/ (double s) const{
        VecXDynamic a = *this;
        for (size_t i = 0; i < coefs.size(); ++i) {
            a.coefs[i] /= s;
        }
        return a;
    }

    inline double operator* (const VecXDynamic& b) const {
        if (b.N != N){
            std::cerr << "Operator * should be used between same sized VecXDynamic" << std::endl;
            exit(1);
        }
        double res = 0;
        for (size_t i = 0; i < coefs.size(); ++i) {
            res += coefs[i] * b.coefs[i];
        }
        return res;
    }

    inline VecXDynamic operator- () const{
        VecXDynamic ret(*this);
        for (double& v : ret.coefs){
            v = -v;
        }
        return ret;
    }

};

inline VecXDynamic operator- (double s, const VecXDynamic& v){
    VecXDynamic a(v.dim());
    for (int i = 0; i < a.dim(); ++i) {
        a[i] = s - v[i];
    }
    return a;
}

inline VecXDynamic operator+ (double s, const VecXDynamic& v){
    VecXDynamic a(v.dim());
    for (int i = 0; i < a.dim(); ++i) {
        a[i] = s + v[i];
    }
    return a;
}

inline VecXDynamic operator* (double s, const VecXDynamic& v){
    VecXDynamic a(v.dim());
    for (int i = 0; i < a.dim(); ++i) {
        a[i] = s * v[i];
    }
    return a;
}

inline VecXDynamic operator/ (double s, const VecXDynamic& v){
    VecXDynamic a(v.dim());
    for (int i = 0; i < a.dim(); ++i) {
        a[i] = s / v[i];
    }
    return a;
}


inline bool operator< (const VecXDynamic& a, const VecXDynamic& b){
    if (a.dim() != b.dim()){
        std::cerr << "Operator < should be used between same sized VecXDynamic" << std::endl;
        exit(1);
    }
    for (int i = 0; i < a.dim(); ++i) {
        if (a[i] < b[i])
            return true;
        else if (a [i] > b[i])
            return false;
    }
}


inline double scalar(const VecXDynamic& a, const VecXDynamic& b){
    return a*b;
}


inline std::ostream& operator<<(std::ostream& out, const VecXDynamic& v){
    for (int i = 0; i < v.dim()-1; ++i){
        out << v[i] << " ";
    }
    out << v[v.dim()-1];
    return out;
}

#endif //SLICEDOPTIM_VECXDYNAMIC_H
