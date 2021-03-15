//
//  cvector.c
//  cMat
//
//  Created by 吕文韬 on 2021/3/13.
//

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include "cmat.h"

/* This function is defined to initial a vector */
void vInitial( vector *T, int m )
{
    (*T).v        = (double _Complex *) malloc( m * sizeof(double _Complex) );
    (*T).row    = m;
}


/* This function is used to make zero vector. */
void vZero( vector *T, int m )
{
    vInitial( T, m );
    for ( int i = 0; i < m; ++i )
    {
        (*T).v[i] = 0;
    }
}


void vStaticZero( vector *T ){
    for ( int i = 0; i < (*T).row; ++i )
    {
        (*T).v[i] = 0;
    }
}


/* This function is used to free vector memory */
void vFree( vector *T )
{
    free( (*T).v );
    (*T).v = NULL;
}


/* This function is used to reInitial the vector memory */
void vReInitial( vector *T, int n )
{
    vFree( T );
    vInitial( T, n );
}


/* This function is defined to get the i-th column vector of matrix A and save to result */
void colVector( vector *result, mat *T, int i )
{
    vStaticZero( result );
    for ( int k = 0; k < (*T).row; ++k )
    {
        (*result).v[k] = (*T).m[k][i - 1];
    }
}


/* This function is defined to get a vector product by a number */
void vNumProduct( vector *a, double num, vector *b )
{
    for ( int i = 0; i < (*b).row; ++i )
    {
        (*a).v[i] = num * (*b).v[i];
    }
}


/* This function is defined to get the inner product of two vectors */
double _Complex vInnerProduct( vector *a, vector *b )
{
    double _Complex product = 0;
    for ( int i = 0; i < (*a).row; ++i )
    {
        product += conj( (*a).v[i] ) * (*b).v[i];
    }
    return(product);
}


/* This function is defined to get the norm of a vector */
double vNorm( vector *a )
{
    double norm = cabs( vInnerProduct( a, a ) );
    norm = cpow( norm, 0.5 );
    return(norm);
}


/* This function is defined to get vectors joint */
void vJoint( int column, mat * result, vector *s[column] )
{
    int    row    = (*s[0]).row;
    _Bool    flag    = 1;
    for ( int i = 0; i < row; ++i )
    {
        if ( s[i]->row != row )
        {
            reInitial( result, 1, 1 );
            (*result).m[0]    = NULL;
            flag        = 0;
        }
    }
    
    if ( flag == 0 )
    {
    }else  {
        reInitial( result, row, row );
        for ( int i = 0; i < row; ++i )
        {
            for ( int j = 0; i < column; ++j )
            {
                result->m[i][j] = s[j]->v[i];
            }
        }
    }
}

//This function is defined to get the the Housholder matrix from a to b
void vHousholder(mat *H /*The Householder Matrix*/, vector *a, vector *b)
{
    vector vtemp;
    vInitial(&vtemp, (*a).row);
    
    for (int i = 0; i < (*a).row; ++i)
    {
        vtemp.v[i] = (*a).v[i] - (*b).v[i];
    }
    
    double norm = vNorm(&vtemp);
    
    vNumProduct(&vtemp, 1/norm, &vtemp);
    
    for (int i = 0; i < (*a).row; ++i)
    {
        for (int j = 0; j < (*a).row; ++j){
            if (i == j){
                (*H).m[i][j] = 1 - 2*vtemp.v[i]*conj(vtemp.v[j]);
            }
            else
                (*H).m[i][j] = -2*vtemp.v[i]*conj(vtemp.v[j]);
        }
    }
    
    vFree(&vtemp);
}

