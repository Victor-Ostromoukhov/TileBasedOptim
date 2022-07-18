//
// Created by lpaulin on 21/03/19.
//

#ifndef SLICEDOPTIM_VECN_H
#define SLICEDOPTIM_VECN_H
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
template <int N>
class VecX {


public:
		double coefs[N];
    inline VecX(const std::vector<double>& values){
        if (values.size() != N){
            std::cerr << "input vector must be of size " << N << " when creating VecX<" << N << std::endl;
        }
		memcpy(coefs, &values[0], N * sizeof(double));
    }

    inline VecX(){}

	explicit inline VecX(int dim) {
		if (dim != N) {
			std::cerr << "input vector must be of size " << N << " when creating VecX<" << N << "> but received " << dim << std::endl;
		}
	}

    inline VecX(double x, double y){
        if (N < 2){
            std::cerr << "Creating a Vec1 from 2 values" << std::endl;
            exit(1);
        }
        coefs[0] = x;
        coefs[1] = y;
    }

    //inline VecX(const VecX<N>& b){}

    const double* data() const{
        return coefs;
    }

    double* data(){
        return coefs;
    }

	inline int dim() const { return (int) N; }

    inline double norm() const{
        return std::sqrt(norm2());
    }
    inline double norm2() const{
        double n = 0.;
		for (int i = 0; i < N; i++) {
			double v = coefs[i];
            n += v*v;
        }
        return n;
    }
    inline void normalize(){
        double n = norm();
		for (int i = 0; i < N; i++) {
			coefs[i] /= n;
        }
    }

    inline double& operator[] (int i){
        return coefs[i];
    }

    inline const double& operator[] (int i) const{
        return coefs[i];
    }

    inline VecX& operator+= (const VecX<N>& b){
        for (int i = 0; i < N; ++i){
            coefs[i] += b.coefs[i];
        }
        return *this;
    }
    inline VecX& operator-= (const VecX b){
        for (int i = 0; i < N; ++i){
            coefs[i] -= b.coefs[i];
        }
        return *this;

    }
    inline VecX& operator/= (double s) {
        for (int i = 0; i < N; ++i) {
            coefs[i] /= s;
        }
        return *this;
    }
    inline VecX& operator*= (double s) {
        for (int i = 0; i < N; ++i) {
            coefs[i] *= s;
        }
        return *this;
    }
    inline VecX& operator+= (double s) {
        for (int i = 0; i < N; ++i) {
            coefs[i] += s;
        }
        return *this;
    }
    inline VecX& operator-= (double s) {
        for (int i = 0; i < N; ++i) {
            coefs[i] -= s;
        }
        return *this;
    }

    inline VecX operator- (const VecX& b) const {
        VecX a = *this;
        for (int i = 0; i < N; ++i) {
            a.coefs[i] -= b.coefs[i];
        }
        return a;
    }
    inline VecX operator+ (const VecX& b) const {
        VecX a = *this;
        for (int i = 0; i < N; ++i) {
            a.coefs[i] += b.coefs[i];
        }
        return a;
    }
    inline VecX operator/ (const VecX& b) const {
        VecX a = *this;
        for (int i = 0; i < N; ++i) {
            a.coefs[i] /= b.coefs[i];
        }
        return a;
    }
    inline VecX operator- (double s) const {
        VecX a = *this;
        for (int i = 0; i < N; ++i) {
            a.coefs[i] -= s;
        }
        return a;
    }
    inline VecX operator+ (double s) const {
        VecX a = *this;
        for (int i = 0; i < N; ++i) {
            a.coefs[i] += s;
        }
        return a;
    }
    inline VecX operator* (double s) const {
        VecX a = *this;
        for (int i = 0; i < N; ++i) {
            a.coefs[i] *= s;
        }
        return a;
    }
    inline VecX operator/ (double s) const {
        VecX a = *this;
        for (int i = 0; i < N; ++i) {
            a.coefs[i] /= s;
        }
        return a;
    }

    inline double operator* (const VecX& b) const {
        double res = 0;
        for (int i = 0; i < N; ++i) {
            res += coefs[i] * b.coefs[i];
        }
        return res;
    }

};

template <int N>
inline VecX<N> operator-(const VecX<N>& v){
    VecX<N> v1;
    for (int i = 0; i < N; ++i){
        v1[i] = - v[i];
    }
    return v1;
}

template <int N>
inline VecX<N> operator- (double s, const VecX<N>& v){
    VecX<N> a;
    for (int i = 0; i < N; ++i) {
        a[i] = s - v[i];
    }
    return a;
}
template <int N>
inline VecX<N> operator+ (double s, const VecX<N>& v){
    VecX<N> a;
    for (int i = 0; i < N; ++i) {
        a[i] = s + v[i];
    }
    return a;
}
template <int N>
inline VecX<N> operator* (double s, const VecX<N>& v){
    VecX<N> a;
    for (int i = 0; i < N; ++i) {
        a[i] = s * v[i];
    }
    return a;
}
template <int N>
inline VecX<N> operator/ (double s, const VecX<N>& v){
    VecX<N> a;
    for (int i = 0; i < N; ++i) {
        a[i] = s / v[i];
    }
    return a;
}

template <int N>
inline VecX<N> operator/ (const VecX<N>& v, double s) {
	VecX<N> a;
	for (int i = 0; i < N; ++i) {
		a[i] =  v[i] / s;
	}
	return a;
}

template <int N>
inline bool operator< (const VecX<N>& a, const VecX<N>& b){
    for (int i = 0; i < N; ++i) {
        if (a[i] < b[i])
            return true;
        else if (a [i] > b[i])
            return false;
    }
}

template <int N>
inline double scalar(const VecX<N>& a, const VecX<N>& b){
    return a*b;
}

template <int N>
inline std::ostream& operator<<(std::ostream& out, const VecX<N>& v){
    for (int i = 0; i < N-1; ++i){
        out << std::setprecision(20) << v[i] << "\t";
    }
    out << v[N-1];
    return out;
}

template<int N>
inline VecX<N> mod(const VecX<N>& v, double m){
    VecX<N> res;
    for (int i = 0; i < N; ++i){
        res[i] = std::fmod(v[i], m);
    }
    return res;
}



#endif //SLICEDOPTIM_VECN_H
