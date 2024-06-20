#include <iostream>
#include <math.h>
#include <string>
#include <cmath>
#include <stack>
#include "pp_transform.h"
#include "rd_error.h"

using namespace::std;
using std::string;

stack<xform> transforms;
point multiply(float scalar, const point& p){

    point product;
     product.x = p.x * scalar;
     product.y = p.y * scalar;
     product.z = p.z * scalar;

     return product;
}

Vector multiply(float scalar, const Vector& v){

    Vector product;
    product.x = scalar * v.x;
    product.y = scalar * v.y;
    product.z = scalar * v.z;

    return product;
}

void copy(point& dest, const point& src){

    dest.x = src.x;
    dest.y = src.y;
    dest.z = src.z;
}

void copy(Vector& dest,  const Vector& src){

    dest.x = src.x;
    dest.y = src.y;
    dest.z = src.z;
}

float dot(const Vector& v1, const Vector& v2){

    return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

Vector cross(const Vector& v1, const Vector& v2){

    Vector result;
    result.x = (v1.y * v2.z) - (v1.z * v2.y);
    result.y = (v1.z * v2.x) - (v1.x * v2.z);
    result.z = (v1.x * v2.y) - (v1.y * v2.x);

    return result;
}

float mag2(Vector v){

    return (v.x * v.x) + (v.y * v.y) + (v.z * v.z);
}

Vector normalize(const Vector& v){

    Vector v1;
    float mag = sqrt(mag2(v));

    v1.x = v.x / mag;
    v1.y = v.y/mag;
    v1.z = v.z/mag;

    return v1;
}

Vector add(const Vector& v1, const Vector& v2){

    Vector result;

    result.x = v1.x + v2.x;
    result.y = v1.y + v2.y;
    result.z = v1.z + v2.z;

    return result;
}

Vector subtract(const Vector& v1, const Vector& v2){

    Vector result;

    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    result.z = v1.z - v2.z;

    return result;
}

Vector subtract(const point& p1, const point& p2){

    Vector result;

    result.x = p2.x - p1.x;
    result.y = p2.y - p1.y;
    result.z = p2.z - p1.z;

    return result;
}

pointh convert(const point& p){

    pointh result;

    result.x = p.x;
    result.y = p.y;
    result.z = p.z;
    result.w = 1.0f;

    return result;
}

pointh convert(const Vector& v){

    pointh result;

    result.x = v.x;
    result.y = v.y;
    result.z = v.z;
    result.w = 1.0f;
    return result;
}

//transformations

xform identity(){

    xform m1;
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            if(i==j){
                m1.matrix[i][j] = 1.0f;
            }
            else{
                m1.matrix[i][j] = 0.0f;
            }
        }
    }
    return m1;
}

xform multiply(const xform& m1, const xform m2){

    xform result;

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            result.matrix[i][j] =0;

            for(int k=0;k<4;k++){
                result.matrix[i][j] += m1.matrix[i][k] * m2.matrix[k][j];
            }
        }

    }

    return result;
}

pointh multiply(const xform& m1, const pointh& p1){

    pointh result;
    result.x = m1.matrix[0][0] * p1.x + m1.matrix[0][1] * p1.y + m1.matrix[0][2] * p1.z + m1.matrix[0][3] * p1.w;
    result.y = m1.matrix[1][0] * p1.x + m1.matrix[1][1] * p1.y + m1.matrix[1][2] * p1.z + m1.matrix[1][3] * p1.w;
    result.z = m1.matrix[2][0] * p1.x + m1.matrix[2][1] * p1.y + m1.matrix[2][2] * p1.z + m1.matrix[2][3] * p1.w;
    result.w = m1.matrix[3][0] * p1.x + m1.matrix[3][1] * p1.y + m1.matrix[3][2] * p1.z + m1.matrix[3][3] * p1.w;

    return result;
}

void push(const xform& m1){

    transforms.push(m1);
}

xform pop(){

    //change identity matrix.
    if(transforms.empty()){
       return identity();
    }

    xform transformation = transforms.top();

    transforms.pop();

    return transformation;

}


xform translate(xform& m1, double tx, double ty, double tz){

    xform translation_matrix = identity();

    translation_matrix.matrix[0][3] = tx;
    translation_matrix.matrix[1][3] = ty;
    translation_matrix.matrix[2][3] = tz;

    m1 = multiply(m1, translation_matrix);

    return m1;

}

xform scale(xform& m1, double sx, double sy, double sz){

    xform scale_matrix = identity();

    scale_matrix.matrix[0][0] = sx;
    scale_matrix.matrix[1][1] = sy;
    scale_matrix.matrix[2][2] = sz;

    m1 = multiply(m1, scale_matrix);

    return m1;
}

void rotate_xy(xform& m1, double theta){

    double radian = theta * M_PI/180;

    xform rotation_matrix = identity();

    rotation_matrix.matrix[0][0] = cos(radian);
    rotation_matrix.matrix[0][1] = -(sin(radian));
    rotation_matrix.matrix[1][0] = sin(radian);
    rotation_matrix.matrix[1][1] = cos(radian);

    m1 = multiply(m1, rotation_matrix);
}

void rotate_yz(xform& m1, double theta){

    double radian = theta * M_PI/180;

    xform rotation_matrix = identity();

    rotation_matrix.matrix[1][1] = cos(radian);
    rotation_matrix.matrix[1][2] = -(sin(radian));
    rotation_matrix.matrix[2][1] = sin(radian);
    rotation_matrix.matrix[2][2] = cos(radian);

    m1 = multiply(m1, rotation_matrix);
}

void rotate_zx(xform& m1, double theta){

    double radian = theta * M_PI/180;

    xform rotation_matrix = identity();

    rotation_matrix.matrix[0][0] = cos(radian);
    rotation_matrix.matrix[0][2] = sin(radian);
    rotation_matrix.matrix[2][0] = -(sin(radian));
    rotation_matrix.matrix[2][2] = cos(radian);

    m1 = multiply(m1, rotation_matrix);
}

xform world_to_camera(const point& eye, const point& at, const Vector& up){

    xform world_to_camera_matrix;
    Vector A = normalize(subtract(eye, at));

    Vector U = normalize(cross(A, up));

    Vector V = cross(U,A);

    xform translation_matrix;
    translation_matrix = identity();

    translation_matrix.matrix[0][3]=-eye.x;
    translation_matrix.matrix[1][3]=-eye.y;
    translation_matrix.matrix[2][3]=-eye.z;

    xform rotation_matrix;
    rotation_matrix.matrix[0][0]=U.x;
    rotation_matrix.matrix[0][1]=U.y;
    rotation_matrix.matrix[0][2]=U.z;
    rotation_matrix.matrix[0][3]=0.0f;

    rotation_matrix.matrix[1][0]=V.x;
    rotation_matrix.matrix[1][1]=V.y;
    rotation_matrix.matrix[1][2]=V.z;
    rotation_matrix.matrix[1][3]=0.0f;

    rotation_matrix.matrix[2][0]=A.x;
    rotation_matrix.matrix[2][1]=A.y;
    rotation_matrix.matrix[2][2]=A.z;
    rotation_matrix.matrix[2][3]=0.0f;

    rotation_matrix.matrix[3][0]=0.0f;
    rotation_matrix.matrix[3][1]=0.0f;
    rotation_matrix.matrix[3][2]=0.0f;
    rotation_matrix.matrix[3][3]=1.0f;

    world_to_camera_matrix = multiply(rotation_matrix, translation_matrix);

    return world_to_camera_matrix;
}


xform camera_to_clip(double fov, double near, double far, double aspect){

    xform camera_to_clip_matrix;
    xform projection_matrix;
    double theta = fov * M_PI / 180.0;


    projection_matrix = identity();

    projection_matrix.matrix[0][0] = 1/( aspect * tan(theta/2.0));
    projection_matrix.matrix[1][1] = 1/tan(theta/2.0);
    projection_matrix.matrix[2][2] = far/(far-near);
    projection_matrix.matrix[2][3] = -(far*near)/(far-near);
    projection_matrix.matrix[3][2] = 1.0f;
    projection_matrix.matrix[3][3] = 0.0f;

    xform transformation_matrix = identity();
    transformation_matrix.matrix[0][0] = 0.5f;
    transformation_matrix.matrix[1][1] = 0.5f;
    transformation_matrix.matrix[0][3] = 0.5f;
    transformation_matrix.matrix[1][3] = 0.5f;

    camera_to_clip_matrix = multiply(transformation_matrix, projection_matrix);

    return camera_to_clip_matrix;

}

xform clip_to_device(int width, int height){

    xform clip_to_device_matrix = identity();
    float epsilon = 1/10;

    clip_to_device_matrix.matrix[0][0] = width - epsilon;
    clip_to_device_matrix.matrix[1][1] = -(height-epsilon);
    clip_to_device_matrix.matrix[1][3] = (height - epsilon);
    clip_to_device_matrix.matrix[0][3] = 0.5;

    return clip_to_device_matrix;

}

