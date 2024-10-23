#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <cstdlib>
#include <iostream>

class Vec3;
static inline Vec3 operator + (Vec3 const & a , Vec3 const & b);
static inline Vec3 operator * (float a , Vec3 const & b);

class Vec3 {
private:
    float mVals[3];
public:
    Vec3() {mVals[0] = mVals[1] = mVals[2] = 0.f;}
    Vec3( float x , float y , float z ) {
       mVals[0] = x; mVals[1] = y; mVals[2] = z;
    }
    float & operator [] (unsigned int c) { return mVals[c]; }
    float operator [] (unsigned int c) const { return mVals[c]; }
    Vec3 operator = (Vec3 const & other) {
       mVals[0] = other[0] ; mVals[1] = other[1]; mVals[2] = other[2];
       return *this;
    }
    float squareLength() const {
       return mVals[0]*mVals[0] + mVals[1]*mVals[1] + mVals[2]*mVals[2];
    }
    float length() const { return sqrt( squareLength() ); }
    inline
    float norm() const { return length(); }
    inline
    float squareNorm() const { return squareLength(); }
    void normalize() { float L = length(); mVals[0] /= L; mVals[1] /= L; mVals[2] /= L; }
    static float dot( Vec3 const & a , Vec3 const & b ) {
       return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    static Vec3 cross( Vec3 const & a , Vec3 const & b ) {
       return Vec3( a[1]*b[2] - a[2]*b[1] ,
                    a[2]*b[0] - a[0]*b[2] ,
                    a[0]*b[1] - a[1]*b[0] );
    }
    void operator += (Vec3 const & other) {
        mVals[0] += other[0];
        mVals[1] += other[1];
        mVals[2] += other[2];
    }
    void operator -= (Vec3 const & other) {
        mVals[0] -= other[0];
        mVals[1] -= other[1];
        mVals[2] -= other[2];
    }
    void operator *= (float s) {
        mVals[0] *= s;
        mVals[1] *= s;
        mVals[2] *= s;
    }
    void operator /= (float s) {
        mVals[0] /= s;
        mVals[1] /= s;
        mVals[2] /= s;
    }
    static Vec3 compProduct(Vec3 const & a , Vec3 const & b) {
        return Vec3(a[0]*b[0] , a[1]*b[1] , a[2]*b[2]);
    }

    unsigned int getMaxAbsoluteComponent() const {
        if( fabs(mVals[0]) > fabs(mVals[1]) ) {
            if( fabs(mVals[0]) > fabs(mVals[2]) ) {
                return 0;
            }
            return 2;
        }
        if( fabs(mVals[1]) > fabs(mVals[2]) ) {
            return 1;
        }
        return 2;
    }
    Vec3 getOrthogonal() const {
        unsigned int c1 = getMaxAbsoluteComponent();
        unsigned int c2 = (c1+1)%3;
        Vec3 res( 0,0,0 );
        res[c1] = mVals[c2];
        res[c2] = -mVals[c1];
        return res;
    }
};

static inline Vec3 operator + (Vec3 const & a , Vec3 const & b) {
   return Vec3(a[0]+b[0] , a[1]+b[1] , a[2]+b[2]);
}
static inline Vec3 operator - (Vec3 const & a , Vec3 const & b) {
   return Vec3(a[0]-b[0] , a[1]-b[1] , a[2]-b[2]);
}
static inline Vec3 operator * (float a , Vec3 const & b) {
   return Vec3(a*b[0] , a*b[1] , a*b[2]);
}
static inline Vec3 operator * (Vec3 const & b , float a ) {
   return Vec3(a*b[0] , a*b[1] , a*b[2]);
}
static inline Vec3 operator / (Vec3 const &  a , float b) {
   return Vec3(a[0]/b , a[1]/b , a[2]/b);
}
static inline std::ostream & operator << (std::ostream & s , Vec3 const & p) {
    s << p[0] << " " << p[1] << " " << p[2];
    return s;
}
static inline std::istream & operator >> (std::istream & s , Vec3 & p) {
    s >> p[0] >> p[1] >> p[2];
    return s;
}

class Mat3 {
public:
    ////////////         CONSTRUCTORS          //////////////
    Mat3() {
        vals[0] = 0;
        vals[1] = 0;
        vals[2] = 0;
        vals[3] = 0;
        vals[4] = 0;
        vals[5] = 0;
        vals[6] = 0;
        vals[7] = 0;
        vals[8] = 0;
    }

    Mat3(float v1, float v2, float v3, float v4, float v5, float v6, float v7, float v8, float v9) {
        vals[0] = v1;
        vals[1] = v2;
        vals[2] = v3;
        vals[3] = v4;
        vals[4] = v5;
        vals[5] = v6;
        vals[6] = v7;
        vals[7] = v8;
        vals[8] = v9;
    }

    Mat3(const Mat3 &m) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                (*this)(i, j) = m(i, j);
    }

    // Multiplication de matrice avec un Vec3 : m.p
    //--> application d'un matrice de rotation à un point ou un vecteur
    Vec3 operator*(const Vec3 &p) {
        //Pour acceder a un element de la matrice (*this)(i,j) et du point p[i]
        Vec3 res = Vec3(
                    (*this)(0, 0) * p[0] + (*this)(0, 1) * p[1] + (*this)(0, 2) * p[2],
                (*this)(1, 0) * p[0] + (*this)(1, 1) * p[1] + (*this)(1, 2) * p[2],
                (*this)(2, 0) * p[0] + (*this)(2, 1) * p[1] + (*this)(2, 2) * p[2]);
        return res;
    }

    Mat3 operator*(const Mat3 &m2) { // calcul du produit matriciel m1.m2
        //Pour acceder a un element de la premiere matrice (*this)(i,j) et de la deuxième m2(k,l)
        Mat3 res = Mat3(
                    (*this)(0, 0) * m2(0, 0) + (*this)(0, 1) * m2(1, 0) + (*this)(0, 2) * m2(2, 0),
                    (*this)(0, 0) * m2(0, 1) + (*this)(0, 1) * m2(1, 1) + (*this)(0, 2) * m2(2, 1),
                    (*this)(0, 0) * m2(0, 2) + (*this)(0, 1) * m2(1, 2) + (*this)(0, 2) * m2(2, 2),
                    (*this)(1, 0) * m2(0, 0) + (*this)(1, 1) * m2(1, 0) + (*this)(1, 2) * m2(2, 0),
                    (*this)(1, 0) * m2(0, 1) + (*this)(1, 1) * m2(1, 1) + (*this)(1, 2) * m2(2, 1),
                    (*this)(1, 0) * m2(0, 2) + (*this)(1, 1) * m2(1, 2) + (*this)(1, 2) * m2(2, 2),
                    (*this)(2, 0) * m2(0, 0) + (*this)(2, 1) * m2(1, 0) + (*this)(2, 2) * m2(2, 0),
                    (*this)(2, 0) * m2(0, 1) + (*this)(2, 1) * m2(1, 1) + (*this)(2, 2) * m2(2, 1),
                    (*this)(2, 0) * m2(0, 2) + (*this)(2, 1) * m2(1, 2) + (*this)(2, 2) * m2(2, 2)
                    );
        return res;
    }

    bool isnan() const {
        return std::isnan(vals[0]) || std::isnan(vals[1]) || std::isnan(vals[2])
                || std::isnan(vals[3]) || std::isnan(vals[4]) || std::isnan(vals[5])
                || std::isnan(vals[6]) || std::isnan(vals[7]) || std::isnan(vals[8]);
    }

    void operator=(const Mat3 &m) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                (*this)(i, j) = m(i, j);
    }

    void operator+=(const Mat3 &m) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                (*this)(i, j) += m(i, j);
    }

    void operator-=(const Mat3 &m) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                (*this)(i, j) -= m(i, j);
    }

    void operator/=(double s) {
        for (unsigned int c = 0; c < 9; ++c)
            vals[c] /= s;
    }

    Mat3 operator-(const Mat3 &m2) {
        return Mat3((*this)(0, 0) - m2(0, 0), (*this)(0, 1) - m2(0, 1), (*this)(0, 2) - m2(0, 2),
                    (*this)(1, 0) - m2(1, 0), (*this)(1, 1) - m2(1, 1), (*this)(1, 2) - m2(1, 2),
                    (*this)(2, 0) - m2(2, 0), (*this)(2, 1) - m2(2, 1), (*this)(2, 2) - m2(2, 2));
    }

    Mat3 operator+(const Mat3 &m2) {
        return Mat3((*this)(0, 0) + m2(0, 0), (*this)(0, 1) + m2(0, 1), (*this)(0, 2) + m2(0, 2),
                    (*this)(1, 0) + m2(1, 0), (*this)(1, 1) + m2(1, 1), (*this)(1, 2) + m2(1, 2),
                    (*this)(2, 0) + m2(2, 0), (*this)(2, 1) + m2(2, 1), (*this)(2, 2) + m2(2, 2));
    }

    Mat3 operator/(float s) {
        return Mat3((*this)(0, 0) / s, (*this)(0, 1) / s, (*this)(0, 2) / s, (*this)(1, 0) / s, (*this)(1, 1) / s,
                    (*this)(1, 2) / s, (*this)(2, 0) / s, (*this)(2, 1) / s, (*this)(2, 2) / s);
    }

    Mat3 operator*(float s) {
        return Mat3((*this)(0, 0) * s, (*this)(0, 1) * s, (*this)(0, 2) * s, (*this)(1, 0) * s, (*this)(1, 1) * s,
                    (*this)(1, 2) * s, (*this)(2, 0) * s, (*this)(2, 1) * s, (*this)(2, 2) * s);
    }

    ////////        ACCESS TO COORDINATES      /////////
    float operator()(unsigned int i, unsigned int j) const { return vals[3 * i + j]; }

    float &operator()(unsigned int i, unsigned int j) { return vals[3 * i + j]; }

    ////////        BASICS       /////////
    inline float sqrnorm() {
        return vals[0] * vals[0] + vals[1] * vals[1] + vals[2] * vals[2]
                + vals[3] * vals[3] + vals[4] * vals[4] + vals[5] * vals[5]
                + vals[6] * vals[6] + vals[7] * vals[7] + vals[8] * vals[8];
    }

    inline float norm() { return sqrt(sqrnorm()); }

    inline float determinant() const {
        return vals[0] * (vals[4] * vals[8] - vals[7] * vals[5])
                - vals[1] * (vals[3] * vals[8] - vals[6] * vals[5])
                + vals[2] * (vals[3] * vals[7] - vals[6] * vals[4]);
    }

    inline float trace() const { return vals[0] + vals[4] + vals[8]; }

    ////////        TRANSPOSE       /////////
    inline
    void transpose() {
        float xy = vals[1], xz = vals[2], yz = vals[5];
        vals[1] = vals[3];
        vals[3] = xy;
        vals[2] = vals[6];
        vals[6] = xz;
        vals[5] = vals[7];
        vals[7] = yz;
    }

    Mat3 getTranspose() const {
        return Mat3(vals[0], vals[3], vals[6], vals[1], vals[4], vals[7], vals[2], vals[5], vals[8]);
    }

    Mat3 operator-() const {
        return Mat3(-vals[0], -vals[1], -vals[2], -vals[3], -vals[4], -vals[5], -vals[6], -vals[7], -vals[8]);
    }


private:
    float vals[9];
    // will be noted as :
    // 0 1 2
    // 3 4 5
    // 6 7 8
};


inline static
Mat3 operator*(float s, const Mat3 &m) {
    return Mat3(m(0, 0) * s, m(0, 1) * s, m(0, 2) * s, m(1, 0) * s, m(1, 1) * s, m(1, 2) * s, m(2, 0) * s, m(2, 1) * s,
                m(2, 2) * s);
}


inline static std::ostream &operator<<(std::ostream &s, Mat3 const &m) {
    s << m(0, 0) << " \t" << m(0, 1) << " \t" << m(0, 2) << std::endl << m(1, 0) << " \t" << m(1, 1) << " \t" << m(1, 2)
      << std::endl << m(2, 0) << " \t" << m(2, 1) << " \t" << m(2, 2) << std::endl;
    return s;
}


#endif
