#ifndef VEC3_H
#define VEC3_H

#include <algorithm>
#include <cmath>
#include <iostream>

class Vec3;
static inline Vec3 operator + (Vec3 const & a , Vec3 const & b);
static inline Vec3 operator * (float a , Vec3 const & b);

class Vec3 {
private:
    float mVals[3];

public:
    Vec3() : mVals{0.f, 0.f, 0.f} {}
    Vec3(const float x, const float y, const float z) : mVals{x, y, z} {}
    explicit Vec3(const float x) : mVals{x, x, x} {}

    float &operator[](const unsigned int c) { return mVals[c]; }
    float operator[](const unsigned int c) const { return mVals[c]; }

    Vec3 &operator=(Vec3 const &other) {
        std::copy_n(other.mVals, 3, mVals);
        return *this;
    }

    Vec3 clamp(const float min, const float max) {
        for (float & mVal : mVals) {
            mVal = std::clamp(mVal, min, max);
        }
        return *this;
    }

    [[nodiscard]] float squareLength() const {
        return mVals[0] * mVals[0] + mVals[1] * mVals[1] + mVals[2] * mVals[2];
    }

    [[nodiscard]] float distanceSquared(Vec3 const &other) const {
        return (mVals[0] - other[0]) * (mVals[0] - other[0]) +
               (mVals[1] - other[1]) * (mVals[1] - other[1]) +
               (mVals[2] - other[2]) * (mVals[2] - other[2]);
    }

    [[nodiscard]] float length() const { return std::sqrt(squareLength()); }
    [[nodiscard]] float norm() const { return length(); }
    [[nodiscard]] float squareNorm() const { return squareLength(); }

    Vec3 &normalize() {
        const float L = length();
        if (L > 0) {
            mVals[0] /= L;
            mVals[1] /= L;
            mVals[2] /= L;
        }
        return *this;
    }

    [[nodiscard]] Vec3 normalized() const {
        const float L = length();
        return L > 0 ? Vec3(mVals[0] / L, mVals[1] / L, mVals[2] / L) : Vec3(0.f);
    }

    bool operator==(const Vec3 &vec3) const {
        return mVals[0] == vec3[0] && mVals[1] == vec3[1] && mVals[2] == vec3[2];
    }

    Vec3 &min(const Vec3 &vec3) {
        for (int i = 0; i < 3; ++i) {
            mVals[i] = std::min(mVals[i], vec3[i]);
        }
        return *this;
    }

    static float dot(Vec3 const &a, Vec3 const &b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    Vec3 operator-() const {
        return {-mVals[0], -mVals[1], -mVals[2]};
    }

    static Vec3 cross(Vec3 const &a, Vec3 const &b) {
        return {a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]};
    }

    Vec3 &operator+=(Vec3 const &other) {
        for (int i = 0; i < 3; ++i) {
            mVals[i] += other[i];
        }
        return *this;
    }

    Vec3 &operator-=(Vec3 const &other) {
        for (int i = 0; i < 3; ++i) {
            mVals[i] -= other[i];
        }
        return *this;
    }

    Vec3 &operator*=(float s) {
        for (float & mVal : mVals) {
            mVal *= s;
        }
        return *this;
    }

    Vec3 &operator/=(const float s) {
        for (float & mVal : mVals) {
            mVal /= s;
        }
        return *this;
    }

    Vec3 operator/(const Vec3 &other) const {
        return {mVals[0] / other[0], mVals[1] / other[1], mVals[2] / other[2]};
    }

    static Vec3 compProduct(Vec3 const &a, Vec3 const &b) {
        return {a[0] * b[0], a[1] * b[1], a[2] * b[2]};
    }

    [[nodiscard]] unsigned int getMaxAbsoluteComponent() const {
        unsigned int maxIndex = 0;
        for (unsigned int i = 1; i < 3; ++i) {
            if (std::fabs(mVals[i]) > std::fabs(mVals[maxIndex])) {
                maxIndex = i;
            }
        }
        return maxIndex;
    }

    [[nodiscard]] Vec3 getOrthogonal() const {
        const unsigned int c1 = getMaxAbsoluteComponent();
        const unsigned int c2 = (c1 + 1) % 3;
        Vec3 res(0.f, 0.f, 0.f);
        res[c1] = mVals[c2];
        res[c2] = -mVals[c1];
        return res;
    }
};

static inline Vec3 operator+(Vec3 const &a, Vec3 const &b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

static inline Vec3 operator-(Vec3 const &a, Vec3 const &b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

static inline Vec3 operator*(const float a, Vec3 const &b) {
    return {a * b[0], a * b[1], a * b[2]};
}

static inline Vec3 operator*(Vec3 const &b, const float a) {
    return a * b;
}

static inline Vec3 operator/(Vec3 const &a, const float b) {
    return {a[0] / b, a[1] / b, a[2] / b};
}

static inline std::ostream &operator<<(std::ostream &s, Vec3 const &p) {
    return s << p[0] << " " << p[1] << " " << p[2];
}

static inline std::istream &operator>>(std::istream &s, Vec3 &p) {
    return s >> p[0] >> p[1] >> p[2];
}

class Mat3 {
private:
    float vals[9]{};

public:
    Mat3() : vals{0.f} {}

    Mat3(const float v1, const float v2, const float v3,
         const float v4, const float v5, const float v6,
         const float v7, const float v8, const float v9) : vals{v1, v2, v3, v4, v5, v6, v7, v8, v9} {}

    Mat3(const Mat3 &m) {
        std::copy_n(m.vals, 9, vals);
    }

    Vec3 operator*(const Vec3 &p) const {
        return {
            vals[0] * p[0] + vals[1] * p[1] + vals[2] * p[2],
            vals[3] * p[0] + vals[4] * p[1] + vals[5] * p[2],
            vals[6] * p[0] + vals[7] * p[1] + vals[8] * p[2]
        };
    }

    Mat3 operator*(const Mat3 &m2) const {
        Mat3 res;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                res(i, j) = (*this)(i, 0) * m2(0, j) +
                            (*this)(i, 1) * m2(1, j) +
                            (*this)(i, 2) * m2(2, j);
            }
        }
        return res;
    }

    [[nodiscard]] bool isnan() const {
        for (const auto &val : vals) {
            if (std::isnan(val)) return true;
        }
        return false;
    }

    Mat3 &operator=(const Mat3 &m) {
        std::copy_n(m.vals, 9, vals);
        return *this;
    }

    Mat3 &operator+=(const Mat3 &m) {
        for (unsigned int i = 0; i < 9; ++i) {
            vals[i] += m.vals[i];
        }
        return *this;
    }

    Mat3 &operator-=(const Mat3 &m) {
        for (unsigned int i = 0; i < 9; ++i) {
            vals[i] -= m.vals[i];
        }
        return *this;
    }

    Mat3 &operator/=(const double s) {
        for (float & val : vals) {
            val /= static_cast<float>(s);
        }
        return *this;
    }

    Mat3 &operator*=(const double s) {
        for (float & val : vals) {
            val *= static_cast<float>(s);
        }
        return *this;
    }

    float &operator()(const unsigned int r, const unsigned int c) {
        return vals[r * 3 + c];
    }

    float operator()(const unsigned int r, const unsigned int c) const {
        return vals[r * 3 + c];
    }

    static Mat3 identity() {
        return {
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        };
    }

    static Mat3 transpose(const Mat3 &m) {
        return {
            m(0, 0), m(1, 0), m(2, 0),
            m(0, 1), m(1, 1), m(2, 1),
            m(0, 2), m(1, 2), m(2, 2)
        };
    }
};

static inline Mat3 operator+(Mat3 const &a, Mat3 const &b) {
    Mat3 res = a;
    res += b;
    return res;
}

static inline Mat3 operator-(Mat3 const &a, Mat3 const &b) {
    Mat3 res = a;
    res -= b;
    return res;
}

static Mat3 operator*(Mat3 const &a, const float s) {
    Mat3 res = a;
    res *= s;
    return res;
}

[[maybe_unused]] static Mat3 operator*(const float s, Mat3 const &a) {
    return a * s;
}


[[maybe_unused]] static std::ostream &operator<<(std::ostream &s, Mat3 const &m) {
    s << m(0, 0) << " \t" << m(0, 1) << " \t" << m(0, 2) << std::endl << m(1, 0) << " \t" << m(1, 1) << " \t" << m(1, 2)
      << std::endl << m(2, 0) << " \t" << m(2, 1) << " \t" << m(2, 2) << std::endl;
    return s;
}

#endif // VEC3_H