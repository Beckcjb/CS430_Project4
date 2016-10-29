#ifndef VECTORMATH_H_INCLUDED
#define VECTORMATH_H_INCLUDED
typedef double* v3;

static inline void v3_add(v3 x, v3 y, v3 c){
    c[0] = x[0] + y[0];
    c[1] = x[1] + y[1];
    c[2] = x[2] + y[2];
}

static inline double v3_dot(v3 a, v3 b){
    double c = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    return c;
}

static inline void v3_cross(v3 a, v3 b, v3 c){
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}
static inline void  v3_subtract(v3 x, v3 y, v3 c){
    c[0] = x[0] - y[0];
    c[1] = x[1] - y[1];
    c[2] = x[2] - y[2];
}
static inline void v3_scale(v3 x, double s, v3 z){ 
    z[0] = x[0] * s;
    z[1] = x[1] * s;
    z[2] = x[2] * s;
}



#endif // VECTORMATH_H_INCLUDED