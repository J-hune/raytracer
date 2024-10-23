#ifndef matrixUtilities_H
#define matrixUtilities_H

#include "Vec3.h"

template< class T >
bool gluInvertMatrix(const T m[16], T invOut[16])
{
    T inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] +
              m[4]  * m[11] * m[14] +
              m[8]  * m[6]  * m[15] -
              m[8]  * m[7]  * m[14] -
              m[12] * m[6]  * m[11] +
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] +
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] -
               m[8]  * m[6] * m[13] -
               m[12] * m[5] * m[10] +
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] +
              m[1]  * m[11] * m[14] +
              m[9]  * m[2] * m[15] -
              m[9]  * m[3] * m[14] -
              m[13] * m[2] * m[11] +
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] +
              m[0]  * m[11] * m[13] +
              m[8]  * m[1] * m[15] -
              m[8]  * m[3] * m[13] -
              m[12] * m[1] * m[11] +
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] +
              m[0]  * m[7] * m[14] +
              m[4]  * m[2] * m[15] -
              m[4]  * m[3] * m[14] -
              m[12] * m[2] * m[7] +
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] +
               m[0]  * m[6] * m[13] +
               m[4]  * m[1] * m[14] -
               m[4]  * m[2] * m[13] -
               m[12] * m[1] * m[6] +
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
              m[1] * m[7] * m[10] +
              m[5] * m[2] * m[11] -
              m[5] * m[3] * m[10] -
              m[9] * m[2] * m[7] +
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
               m[0] * m[7] * m[9] +
               m[4] * m[1] * m[11] -
               m[4] * m[3] * m[9] -
               m[8] * m[1] * m[7] +
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return false;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return true;
}






template< class T  >
void mult2( const T m[16] , T x , T y , T z , T w , T & resX , T & resY , T & resZ , T & resW ) {
    resX = m[0] * x + m[1] * y + m[2] * z + m[3] * w;
    resY = m[4] * x + m[5] * y + m[6] * z + m[7] * w;
    resZ = m[8] * x + m[9] * y + m[10] * z + m[11] * w;
    resW = m[12] * x + m[13] * y + m[14] * z + m[15] * w;
}

template< class T  >
void mult2( const T m[16] , T x[4] , T res[4] ) {
    res[0] = m[0] * x[0] + m[1] * x[1] + m[2] * x[2] + m[3] * x[3];
    res[1] = m[4] * x[0] + m[5] * x[1] + m[6] * x[2] + m[7] * x[3];
    res[2] = m[8] * x[0] + m[9] * x[1] + m[10] * x[2] + m[11] * x[3];
    res[3] = m[12] * x[0] + m[13] * x[1] + m[14] * x[2] + m[15] * x[3];
}







template< class T  >
void mult( const T m[16] , T x , T y , T z , T w , T & resX , T & resY , T & resZ , T & resW ) {
    resX = m[0] * x + m[4] * y + m[8] * z + m[12] * w;
    resY = m[1] * x + m[5] * y + m[9] * z + m[13] * w;
    resZ = m[2] * x + m[6] * y + m[10] * z + m[14] * w;
    resW = m[3] * x + m[7] * y + m[11] * z + m[15] * w;
}


template< class T  >
void mult( const T m[16] , T x[4] , T res[4] ) {
    res[0] = m[0] * x[0] + m[4] * x[1] + m[8] * x[2] + m[12] * x[3];
    res[1] = m[1] * x[0] + m[5] * x[1] + m[9] * x[2] + m[13] * x[3];
    res[2] = m[2] * x[0] + m[6] * x[1] + m[10] * x[2] + m[14] * x[3];
    res[3] = m[3] * x[0] + m[7] * x[1] + m[11] * x[2] + m[15] * x[3];
}







// These functions are not optimized, because you probably don't want to invert the camera matrices every time!
Vec3 cameraSpaceToWorldSpace(Vec3 const & pCS) { // pCS : p in Camera Space
    GLdouble modelview[16];
    GLdouble modelviewInverse[16];
    glMatrixMode (GL_MODELVIEW);
    glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    gluInvertMatrix( modelview , modelviewInverse );
    GLdouble res[4];
    mult(modelviewInverse , (GLdouble)pCS[0] , (GLdouble)pCS[1] , (GLdouble)pCS[2] , (GLdouble)1.0 , res[0] , res[1] , res[2] , res[3]);
    return Vec3( res[0] / res[3] , res[1] / res[3] , res[2] / res[3] );
}
Vec3 screen_space_to_worldSpace( float u , float v ) {
    // u et v sont entre 0 et 1 (0,0 est en haut a gauche de l'ecran)
    GLdouble projection[16];
    GLdouble projectionInverse[16];
    glMatrixMode (GL_PROJECTION);
    glGetDoublev( GL_PROJECTION_MATRIX, projection );
    gluInvertMatrix( projection , projectionInverse );
    GLdouble nearAndFarPlanes[2];
    glGetDoublev( GL_DEPTH_RANGE , nearAndFarPlanes );
    GLdouble modelview[16];
    GLdouble modelviewInverse[16];
    glMatrixMode (GL_MODELVIEW);
    glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    gluInvertMatrix( modelview , modelviewInverse );
    GLdouble resInt[4];
    mult(projectionInverse , (GLdouble)2.f*u - 1.f , -((GLdouble)2.f*v - 1.f) , nearAndFarPlanes[0] , (GLdouble)1.0 , resInt[0] , resInt[1] , resInt[2] , resInt[3]);
    GLdouble res[4];
    mult(modelviewInverse , resInt[0] , resInt[1] , resInt[2] , resInt[3] , res[0] , res[1] , res[2] , res[3]);
    return Vec3( res[0] / res[3] , res[1] / res[3] , res[2] / res[3] );
}
void screen_space_to_world_space_ray(float u , float v , Vec3 & position , Vec3 & direction) {
    position = cameraSpaceToWorldSpace( Vec3(0,0,0) );
    direction = screen_space_to_worldSpace(u,v) - position;
    direction.normalize();
}

#endif // matrixUtilities_H
