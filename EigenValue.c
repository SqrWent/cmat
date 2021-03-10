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
#include "cmat.c"
#include "cmat.h"


/*
 * The function was defined to get the rank of a matrix input.
 */
int mRank( mat *T )
{
    int    rank;
    int    row    = (*T).row;
    int    column    = (*T).col;
    mat    Temp;
    mTranspose( &Temp, T );
    if ( column >= row )
    {
        rank = row;
        upTra( &Temp );
        for ( int i = 0; i < row; ++i )
        {
            if ( (Temp).m[i][i] == 0 )
            {
                --rank;
            }
        }
    }
    else
    {
        mEqual( &Temp, T );
        rank = column;
        upTra( &Temp );
        for ( int i = 0; i < row; ++i )
        {
            if ( (Temp).m[i][i] == 0 )
            {
                --rank;
            }
        }
    }
    return(rank);
}


//This function is defined to initial a vector
void vInitial(vector *T, int m)
{
    (*T).row = m;
    (*T).v = (double *)malloc(m * sizeof(double));
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
    vReInitial(result, (*T).row);
    for (int k = 0; k < (*T).row; ++k)
    {
        (*result).v[k] = (*T).m[k][i];
    }
}




