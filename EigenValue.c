/*
 *
 *  EigenValue.c
 *  cMat
 *
 *  Created by 吕文韬 on 2021/3/9.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "cmat.h"

//This function is defined to initial a vector
void vInitial(vector *T, int m)
{
    (*T).row = m;
    (*T).v = (double _Complex *)malloc(m * sizeof(double _Complex));
}


//This function is used to make zero vector.
void vZero(vector *T, int m)
{
    vInitial(T, m);
    for (int i = 0; i < m; ++i)
    {
        (*T).v[i] = 0;
    }
}


//This function is used to free vector memory
void vFree(vector *T)
{
    free((*T).v);
}


//This function is used to reInitial the vector memory
void vReInitial(vector *T, int n)
{
    vFree(T);
    vInitial(T, n);
}


//This function is defined to get the i-th column vector of matrix A and save to result
void colVector(vector *result, mat *T, int i)
{
    vZero(result, (*T).row);
    for (int k = 0; k < (*T).row; ++k)
    {
        (*result).v[k] = (*T).m[k][i - 1];
    }
}


//This function is defined to get a vector product by a number
void vNumProduct(vector *a, double num, vector *b){
    vZero(a, (*b).row);
    for (int i = 0; i < (*b).row; ++i){
        (*a).v[i] = num*(*b).v[i];
    }
}

//This function is defined to get the inner product of two vectors
double _Complex vInnerProduct(vector *a, vector *b)
{
    double _Complex product = 0;
    for (int i = 0; i < (*a).row; ++i){
        product += conj((*a).v[i])*(*b).v[i];
    }
    return product;
}


//This function is defined to get the norm of a vector
double vNorm(vector *a)
{
    double norm = cabs(vInnerProduct(a, a));
    norm = pow(norm, 0.5);
    return norm;
}


//This function is defined to get vectors joint
void vJoint(int column, mat * result, vector *s[column])
{
    int row = (*s[0]).row;
    _Bool flag = 1;
    for (int i = 0; i < row; ++i){
        if (s[i]->row != row){
            reInitial(result, 1, 1);
            (*result).m[0] = NULL;
            flag = 0;
        }
    }
    
    if (flag == 0){
        
    }
    else{
        reInitial(result, row, row);
        for (int i = 0; i < row; ++i){
            for (int j = 0; i < column; ++j){
                result->m[i][j] = s[j]->v[i];
            }
        }
    }
}


//This function is defined to carry Gram-Schmidt procedure on vectors.
void vGS(mat *result, mat *input){
    int size = input ->row;
    vector *v[size];
    
    for (int i = 0; i < size; ++i){
        colVector(v[i], input, i+1);
    }
    
    zeros(result, size, size);
    
    vector temp;
    for (int i = 0; i < size; ++i){
        result ->m[i][i] = vNorm(v[i]);
        vNumProduct(&temp, 1/(result->m[i][i]), v[i]);
        for (int j = i+1; j < size; ++j){
            (*result).m[i][j] = vInnerProduct(v[j], &temp);
            for (int k = 0; k < size; ++k){
                v[j] -> v[k] = v[j] -> v[k] - (*result).m[i][j]*temp.v[i];
            }
        }
    }
}


//This function is defined to carry Householder Transformation on matrix

