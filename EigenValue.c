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

    mStaticZero( Qresult );
    mStaticZero( Rresult );

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


/* This function is defined to do Housholder transformation to matrix inputed. */
void mHousholder( mat * Qresult /*The unitary matrix*/,
          mat * Rresult /* The upper triangle matrix*/,
          mat * input /*The input matrix*/ )
{
    int size = (*input).row;


    /*
     * Initial the Q to be an identity matrix
     * Initial R to be equal to input
     */
    for ( int i = 0; i < size; ++i )
    {
        for ( int j = 0; j < size; ++j )
        {
            if ( i != j )
            {
                (*Qresult).m[i][j] = 0;
            }else
                (*Qresult).m[i][j] = 1;

            (*Rresult).m[i][j] = (*input).m[i][j];
        }
    }

    vector    temp;
    vector    e;
    mat    mTemp;
    mat    H;
    initial( &H, size, size );
    double norm;
    for ( int i = 0; i < size; ++i )
    {
        vInitial( &temp, size - i );
        vZero( &e, size - i );
        initial( &mTemp, size - i, size - i );

        /* get the value of temp and e */
        for ( int k = i; k < size; ++k )
        {
            temp.v[k] = (*Rresult).m[k][i];
        }
        norm    = vNorm( &temp );
        e.v[0]    = norm;


        vHousholder( &mTemp, &temp, &e );

        for ( int k = 0; k < size; ++k )
        {
            for ( int j = 0; j < size; ++j )
            {
                if ( k >= i && j >= i )
                {
                    H.m[k][j] = mTemp.m[k - i][j - i];
                }else if ( k == j )
                {
                    H.m[k][j] = 1;
                }else
                    H.m[k][j] = 0;
            }
        }

        mProduct( Rresult, &H, Rresult );
        mProduct( Qresult, Qresult, &H );

        mFree( &mTemp );
        vFree( &temp );
        vFree( &e );
    }
    mFree( &H );
}


/* This function is defined to get the Eigen Value of matrix via Householder method. */
void eigValueHS( double _Complex * eigvalues/* The pointer to eigenvalues*/,
         mat *T/*The target matrix*/,
         int times/*The iteration times*/ )
{
    mat Q, R, temp;
    mEqual( &temp, T );
    int size = T->row;

    for ( int i = 0; i < times; ++i )
    {
        mHousholder( &Q, &R, &temp );
        mProduct( &temp, &Q, &R );
    }

    for ( int i = 0; i < size; ++i )
    {
        eigvalues[i] = temp.m[i][i];
    }
}


/*This function is defined to do QR decomposition on a matrix via Givens means*/
void mGivens(mat * Qresult/*The unitary matrix*/,
             mat * Rresult/*The tpper triangle matrix*/,
             mat * input/*The input matrix*/)
{
    int size = (*input).row;
    
    /*
     Initial Qresult and Rresult
     */
    for ( int i = 0; i < size; ++i )
    {
        for ( int j = 0; j < size; ++j )
        {
            if ( i != j )
            {
                (*Qresult).m[i][j] = 0;
            }else
                (*Qresult).m[i][j] = 1;

            (*Rresult).m[i][j] = (*input).m[i][j];
        }
    }
    
    mat Givens;
    zeros(&Givens, size, size);
    double _Complex c,s;
    c = 0;
    s = 0;
    double norm1,norm2;
    double _Complex temp;
    temp = 0;
    norm1 = 0;
    norm2 = 0;
    
    for (int i = 0; i < size; ++i)
    {
        
        mStaticZero(&Givens);
        for (int j = 0; j < size; ++j)
        {
            Givens.m[j][j] = 1;
        }
        
        for (int j = i + 1; j < size; ++j)
        {
            
            // Check if denominator == 0
            if (cabs((*Rresult).m[j][i]) == 0 && cabs((*Rresult).m[i][i]) == 0)
            {
                c = 1;
                s = 0;
            }
            else
            {
                norm1 = cabs((*Rresult).m[i][i]);
                norm2 = cabs((*Rresult).m[j][i]);
                
                c = (*Rresult).m[i][i]/sqrt(pow(norm1, 1)+pow(norm2, 2));
                s = (*Rresult).m[j][i]/sqrt(pow(norm1, 1)+pow(norm2, 2));
            }
            
            // Set Value for Givens Matrix
            if (j == i + 1)
            {
                Givens.m[i][i] = conj(c);
                Givens.m[j][j] = c;
                Givens.m[i][j] = conj(s);
                Givens.m[j][i] = -s;
            }
            else
            {
                Givens.m[i][i] = conj(c);
                Givens.m[j][j] = c;
                Givens.m[i][j] = conj(s);
                Givens.m[j][i] = -s;
                Givens.m[j-1][j-1] = 1;
                Givens.m[i][j-1] = 0;
                Givens.m[j-1][i] = 0;
            }
            
            mProduct(Rresult, &Givens, Rresult);
            
            //Transpose Givens Matrix.
            Givens.m[i][i] = c;
            Givens.m[j][j] = conj(c);
            Givens.m[i][j] = -conj(s);
            Givens.m[j][i] = s;
            
            mProduct(Qresult, Qresult, &Givens);
            
        }
    }
    
}


/* This function is defined to get the Eigen Value of matrix via Given method. */
void eigValueGVS( double _Complex * eigvalues/* The pointer to eigenvalues*/,
         mat *T/*The target matrix*/,
         int times/*The iteration times*/ )
{
    mat Q, R, temp;
    mEqual( &temp, T );
    int size = T->row;

    for ( int i = 0; i < times; ++i )
    {
        mGivens( &Q, &R, &temp );
        mProduct( &temp, &Q, &R );
    }

    for ( int i = 0; i < size; ++i )
    {
        eigvalues[i] = temp.m[i][i];
    }
}
