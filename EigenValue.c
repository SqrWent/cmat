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


/* This function is defined to carry Gram-Schmidt procedure on vectors. */
void MGS( mat *Qresult, mat *Rresult, mat *input )
{
    int    size = input->row;
    vector    *v;
    v = (vector *) malloc( size * size * sizeof(double _Complex) );
    for ( int i = 0; i < size; ++i )
    {
        vInitial( &v[i], size );
        colVector( &v[i], input, i + 1 );
    }

    zeros( Qresult, size, size );
    zeros( Rresult, size, size );

    vector temp;

    for ( int i = 0; i < size; ++i )
    {
        Rresult->m[i][i] = vNorm( &v[i] );

        vNumProduct( &temp, 1 / (Rresult->m[i][i]), &v[i] );

        for ( int j = i + 1; j < size; ++j )
        {
            (*Rresult).m[i][j] = vInnerProduct( &v[j], &temp );
            for ( int k = 0; k < size; ++k )
            {
                v[j].v[k] = v[j].v[k] - (*Rresult).m[i][j] * temp.v[k];
            }
        }
        for ( int k = 0; k < size; ++k )
        {
            Qresult->m[k][i] = temp.v[k];
        }
    }

    /* Free temporary vectors */
    vFree( &temp );
    for ( int i = 0; i < size; ++i )
    {
        vFree( &v[i] );
    }
    free( v );
}


/* This function is defined get the Eigen Value of matrix via MGS method. */
void eigValueMGS( double _Complex * eigvalues, mat *T, int times )
{
    mat Q, R, temp;
    mEqual( &temp, T );
    int size = T->row;

    for ( int i = 0; i < times; ++i )
    {
        MGS( &Q, &R, &temp );
        mProduct( &temp, &Q, &R );
    }

    for ( int i = 0; i < size; ++i )
    {
        eigvalues[i] = temp.m[i][i];
    }
}


/* This function is defined to do Hessenburg transformation to matrix inputed. */
