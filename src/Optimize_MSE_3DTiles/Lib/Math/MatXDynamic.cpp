//

#include "MatXDynamic.h"
#include <string>

MatXDynamic::MatXDynamic(int w, int h) : width(w), height(h), data(w*h, 0.){}

MatXDynamic::MatXDynamic(const std::vector<double>& dat, int w, int h) : width(w), height(h), data(dat){}

MatXDynamic::operator double() const{
    if (width != 1 || height != 1){
        std::cerr << "Error: Can't convert Vector to double if dimension is higher than 1." << std::endl;
        std::cerr << "\tCurrent dimension is (" << width << ", " << height << ")" << std::endl;
        exit(1);
    }
    return data[0];
}

double& MatXDynamic::operator()(int i, int j){
    return data[i + j * width];
}

double MatXDynamic::operator()(int i, int j) const{
    return data[i + j * width];
}

double& MatXDynamic::access(int i, int j){
    return data[i + j * width];

}
double MatXDynamic::access(int i, int j) const {
    return data[i + j * width];
}

MatXDynamic MatXDynamic::transpose() const{
    MatXDynamic m(height, width);

    for (int j = 0; j < height; ++j){
        for (int i = 0; i < width; ++i){
            m(j, i) = access(i, j);
        }
    }

    return m;
}

MatXDynamic operator* (const MatXDynamic& m1, const MatXDynamic& m2){
    if (m1.width != m2.height){
        std::cerr << "Error: Incompatible matrix sizes ("
                    + std::to_string(m1.width) + ", " + std::to_string(m1.height) + ") x ("
                    + std::to_string(m2.width) + ", " + std::to_string(m2.height) + ") for * operator" << std::endl;
        exit(1);
    }

    MatXDynamic m(m2.width, m1.height);


    for (int j = 0; j < m.height; ++j){
        for (int i = 0; i < m.width; ++i){

            for (int k = 0; k < m1.width; ++k){
                m(i,j) += m1(k, j) * m2(i, k);
            }

        }
    }

    return m;
}

VecXDynamic operator* (const VecXDynamic& v1, const MatXDynamic& m2){

    if (v1.dim() != int(m2.height)){
        std::cerr << "Error: Incompatible matrix sizes ("
                     + std::to_string(v1.dim()) + ") x ("
                     + std::to_string(m2.width) + ", " + std::to_string(m2.height) + ") for * operator" << std::endl;
        exit(1);
    }

    VecXDynamic v(v1.dim());

    for (int i = 0; i < v.dim(); ++i){
        for (int k = 0; k < v.dim(); ++k){
            v[i] += v1[k] * m2(i,k);
        }
    }
    return v;
}

VecXDynamic operator* (const MatXDynamic& m1, const VecXDynamic& v2){

    if (v2.dim() != m1.width){
        std::cerr << "Error: Incompatible matrix sizes ("
                     + std::to_string(m1.width) + ", " + std::to_string(m1.height) + ") x ("
                     + std::to_string(v2.dim()) + ") for * operator" << std::endl;
        exit(1);
    }

    VecXDynamic v(v2.dim());

    for (int i = 0; i < v.dim(); ++i){
        for (int k = 0; k < v.dim(); ++k){
            v[i] += v2[k] * m1(k,i);
        }
    }
    return v;
}

std::ostream& operator<<(std::ostream& out, const MatXDynamic& m){
    for (int j = 0; j < m.height; ++j){
        for (int i = 0; i < m.width; ++i){
            out << m(i,j);
            if (i != m.width - 1)
                out << " ";
        }
        if (j != m.height - 1)
            out << std::endl;
    }
    return out;
}

/**
 * Computes LUP decomposition of the matrix such that L has 1 on its diagonal.
 * The last cell of P contains the number of permutations
 *
 * @param L Lower triangular matrix
 * @param U Upper triangular matrix
 * @param P Permutations between lines of the original matrix
 * @return
 */
bool MatXDynamic::LUPDecompose(MatXDynamic& L, MatXDynamic& U, std::vector<int>& P) const{

    if (width != height){
        std::cerr << "Error: LUDecompose must be called on square matrices" << std::endl;
        exit(1);
    }

    //We compute L and U in a single matrix A such that A = (L - I) + U
    //Since L has 1 on the diagonal both half of this sum do not interfere with each other
    MatXDynamic A(width, height);
    for (int j = 0; j < height; ++j){
        for (int i = 0; i < width; ++i){
            A(i,j) = access(i,j);
        }
    }

    static const double eps = std::pow(10, -15);
    int N = width;
    P.resize(N + 1);

    //P starts with no permutation
    for (int i = 0; i < N; i++) {
        P[i] = i;
    }
    //P[N] will contain the number of permutations
    P[N] = 0;

    for (int i = 0; i < N; i++) {
        double maxA = 0.0;
        int imax = i;
        double absA;

        //Find line with highest first factor and move it to the top
        for (int k = i; k < N; k++) {
            if ((absA = fabs(A(i, P[k]))) > maxA) {
                maxA = absA;
                imax = k;
            }
        }

        //If a whole column is 0 then the matrix is degenerate and cannot be factorised
        if (maxA < eps) return false;

        //If best line is not the first one then swap them
        if (imax != i) {
            int j = P[i];
            P[i] = P[imax];
            P[imax] = j;
            //counting pivots (for determinant)
            P[N]++;
        }


        for (int j = i + 1; j < N; j++) {
            //We divide the lower triangle to make it so that L has 1 on the diagonal
            A(i, P[j]) /= A(i, P[i]);
            //Sum from Doolittle method
            for (int k = i + 1; k < N; k++) {
                A(k, P[j]) -= A(i, P[j]) * A(k, P[i]);
            }
        }
    }
    for (int j = 0; j < height; ++j){
        for (int i = 0; i < j; ++i){
            L(i,j) = A(i,P[j]);
        }
    }
    for (int i = 0; i < width; ++i){
        L(i,i) = 1.;
    }
    for (int j = 0; j < height; ++j){
        for (int i = j; i < width; ++i){
            U(i,j) = A(i,P[j]);
        }
    }

    return true;  //decomposition done

}

double MatXDynamic::det() const{
    MatXDynamic L(width, height), U(width, height);
    std::vector<int> P(width + 1);
    //If matrix is degenerate det = 0
    if (!LUPDecompose(L, U, P)){
        return 0.;
    }
    double res = 1.;
    for (int i = 0; i < width; ++i){
        res *= U(i, i);
    }
    return (P[width]%2 ? -1 : 1) * res;
}

double MatXDynamic::highestEigenPair(VecXDynamic& v) const{
    std::vector<VecXDynamic> vs(2, VecXDynamic(width));
    static const double eps = pow(10, -10);

    double l = 0.;

    for (int i = 0; i < width; ++i){
        vs[0][i] = 1.;
    }
    vs[0].normalize();

    int act = 0;
    VecXDynamic diff;

    do{
        vs[(1 + act)%2] = *this * vs[act];
        act = (act + 1) % 2;
        l = vs[act].norm();
        vs[act] /= l;
        diff = vs[0] - vs[1];
    }while (diff.norm() > eps);

    v = vs[act];
    return l;
}

/**
 * Finds x such that Ax = b with A this matrix
 * Only works with square matrix
 *
 * @param b
 * @return
 */
VecXDynamic MatXDynamic::solve(const VecXDynamic& b) const{

    MatXDynamic L(width, width), U(width, width);
    std::vector<int> P(width + 1);

    if (!LUPDecompose(L, U, P)){
        std::cerr << "Error: Underdefined system does not have a unique solution" << std::endl;
        exit(1);
    }

    VecXDynamic x(b.dim());

    for (int i = 0; i < b.dim(); ++i){
        x[i] = b[P[i]];
    }

    //Solve L*y = P*b
    //L is lower triangular so we can iterate on the rows starting from the first one
    for (int j = 0; j < height; ++j){
        x[j] /= L(j, j);
        //Cancel y_j * L(j, i) contribution for each row i
        for (int i = j + 1; i < height; ++i){
            x[i] -= x[j] * L(j, i);
        }
    }

    //Solve U*x = y
    //U is upper triangular so we can iterate on the rows starting from the last one
    for (int j = height - 1; j >= 0; --j){
        x[j] /= U(j, j);
        //Cancel x_j * U(j, i) contribution for each row i
        for (int i = j - 1; i >= 0; --i){
            x[i] -= x[j] * U(j, i);
        }
    }

    return x;
}


MatXDynamic MatXDynamic::inverse() const{

    MatXDynamic L(width, width), U(width, width);
    std::vector<int> P(width + 1);

    if(!LUPDecompose(L, U, P)){
        std::cerr << "Error: Matrix cannot be inverted" << std::endl;
    }

    //Construct idententity with row swaps from LUP
    MatXDynamic inv(height, width);
    for (int i = 0; i < width; ++i){
        inv(P[i], i) = 1.;
    }

    //Solve for each column of identity
    for (int col = 0; col < width; ++col) {
        //Solve L*y = P*b
        //L is lower triangular so we can iterate on the rows starting from the first one
        for (int j = 0; j < height; ++j) {
            inv(col, j) /= L(j, j);
            //Cancel y_j * L(j, i) contribution for each row i
            for (int i = j + 1; i < height; ++i) {
                inv(col, i) -= inv(col, j) * L(j, i);
            }
        }

        //Solve U*x = y
        //U is upper triangular so we can iterate on the rows starting from the last one
        for (int j = height - 1; j >= 0; --j) {
            inv(col, j) /= U(j, j);
            //Cancel x_j * U(j, i) contribution for each row i
            for (int i = j - 1; i >= 0; --i) {
                inv(col, i) -= inv(col, j) * U(j, i);
            }
        }
    }

    return inv;

}


MatXDynamic& MatXDynamic::operator+=(const MatXDynamic& m){
    if (width != m.width || height != m.height){
        std::cerr << "Error: Incompatible matrix sizes ("
                     + std::to_string(width) + ", " + std::to_string(height) + ") x ("
                     + std::to_string(m.width) + ", " + std::to_string(m.height) + ") for += operator" << std::endl;
        exit(1);
    }
    for (int j = 0; j < height; ++j){
        for (int i = 0; i < width; ++i) {
            access(i, j) += m(i,j);
        }
    }
    return *this;
}

MatXDynamic operator+ (const MatXDynamic& m1, const MatXDynamic& m2){
    if (m1.width != m2.width || m1.height != m2.height){
        std::cerr << "Error: Incompatible matrix sizes ("
                     + std::to_string(m1.width) + ", " + std::to_string(m2.height) + ") x ("
                     + std::to_string(m2.width) + ", " + std::to_string(m2.height) + ") for + operator" << std::endl;
        exit(1);
    }
    MatXDynamic m(m1.width, m2.width);
    for (int j = 0; j < m1.height; ++j){
        for (int i = 0; i < m1.width; ++i) {
            m(i, j) = m1(i,j) + m2(i,j);
        }
    }
    return m;
}

MatXDynamic& MatXDynamic::operator-=(const MatXDynamic& m){
    if (width != m.width || height != m.height){
        std::cerr << "Error: Incompatible matrix sizes ("
                     + std::to_string(width) + ", " + std::to_string(height) + ") x ("
                     + std::to_string(m.width) + ", " + std::to_string(m.height) + ") for -= operator" << std::endl;
        exit(1);
    }
    for (int j = 0; j < height; ++j){
        for (int i = 0; i < width; ++i) {
            access(i, j) -= m(i,j);
        }
    }
    return *this;
}

MatXDynamic operator-(const MatXDynamic& m1, const MatXDynamic& m2){
    if (m1.width != m2.width || m1.height != m2.height){
        std::cerr << "Error: Incompatible matrix sizes ("
                     + std::to_string(m1.width) + ", " + std::to_string(m2.height) + ") x ("
                     + std::to_string(m2.width) + ", " + std::to_string(m2.height) + ") for - operator" << std::endl;
        exit(1);
    }
    MatXDynamic m(m1.width, m2.width);
    for (int j = 0; j < m1.height; ++j) {
        for (int i = 0; i < m1.width; ++i) {
            m(i, j) = m1(i,j) - m2(i,j);
        }
    }
    return m;
}


MatXDynamic operator* (double m1, const MatXDynamic& m2){
    MatXDynamic m(m2);
    for (int j = 0; j < m.height; ++j) {
        for (int i = 0; i < m.width; ++i) {
            m(i,j) *= m1;
        }
    }
    return m;
}
MatXDynamic operator* (const MatXDynamic& m1, double m2){
    return m2 * m1;
}

MatXDynamic transpose(const MatXDynamic& m){
    return m.transpose();
}

MatXDynamic MatXDynamic::Identity(int d){
    MatXDynamic matXDynamic(d, d);

    for(int i = 0; i < d; ++i){
        matXDynamic(i,i) = 1.;
    }

    return matXDynamic;
}
