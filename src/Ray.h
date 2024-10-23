#ifndef RAY_H
#define RAY_H
#include "Line.h"
class Ray : public Line {
public:
    Ray() : Line() {}
    Ray( Vec3 const & o , Vec3 const & d ) : Line(o,d) {}
};
#endif
