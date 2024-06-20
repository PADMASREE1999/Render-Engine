#ifndef PP_TRANSFORM_H_INCLUDED
#define PP_TRANSFORM_H_INCLUDED

#include <string>
#include <stack>
#include "attr_point.h"
using std::string;

struct point{
    float x;
    float y;
    float z;
};

struct pointh{
    float x;
    float y;
    float z;
    float w;
};

struct Vector{
    float x;
    float y;
    float z;
};

struct xform{
    float matrix[4][4];
};

struct Edge{
    int yLast;
    attr_point p;
    attr_point inc;

    Edge * next;
};

enum Boundary{
    LEFT,
    RIGHT,
    TOP,
    BOTTOM,
    FRONT,
    BACK
};

extern std::stack<xform> transforms;

point multiply(float scalar, const point& p);

Vector multiply(float scalar, const Vector& v);

void copy(point& dest, const point& src);
void copy(Vector& dest,  const Vector& src);

float dot(const Vector& v1, const Vector& v2);

Vector cross(const Vector& v1, const Vector& v2);

float mag2(Vector v);

Vector normalize(const Vector& v);


Vector add(const Vector& v1, const Vector& v2);

Vector subtract(const Vector& v1, const Vector& v2);

Vector subtract(const point& p1, const point& p2);

pointh convert(const point& p);

pointh convert(const Vector& v);

//transformations

xform identity();

xform multiply(const xform& m1, const xform m2);
pointh multiply(const xform& m1, const pointh& p1);

void push(const xform& m1);

xform pop();


xform translate(xform& m1, double tx, double ty, double tz);

xform scale(xform& m1, double sx, double sy, double sz);

void rotate_xy(xform& m1, double theta);

void rotate_yz(xform& m1, double theta);

void rotate_zx(xform& m1, double theta);


xform world_to_camera(const point& eye, const point& at, const Vector& up);


xform camera_to_clip(double fov, double near, double far, double aspect);

xform clip_to_device(int width, int height);


#endif // PP_TRANSFORM_H_INCLUDED
