//
//  cmat.c
//  cmat
//
//  Created by 吕文韬 on 2021/3/9.
//

#include "cmat.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>


/*Initialize a matrix*/
void initial( mat *T, int row, int column )
{
    (*T).m = (double _Complex * *) malloc( row * sizeof(double _Complex *) );
    for ( int i = 0; i < row; ++i )
    {
        ( (*T).m)[i] = (double _Complex *) malloc( column * sizeof(double _Complex) );
    }
    (*T).row    = row;
    (*T).col    = column;
}


/*
 * Generalize zero matrix
 */
void zeros( mat *T, int row, int column )
{
    initial( T, row, column );
    for ( int i = 0; i < row; ++i )
    {
        for ( int j = 0; j < column; ++j )
        {
            (*T).m[i][j] = 0;
        }
    }
}


/*
 * Free matrix memory
 */
void mFree( mat *T )
{
    for ( int i = 0; i < (*T).row; ++i )
    {
        free( (*T).m[i] );
        (*T).m[i] = NULL;
    }
    free( (*T).m );
    (*T).m = NULL;
}


/*
 * Read matrix
 */
/*
 Function temporarily unavailabe since compex number added.
void mScanf( mat *T )
{
    for ( int i = 0; i < (*T).row; ++i )
    {
        for ( int j = 0; j < (*T).col; ++j )
        {
            scanf( "%lf", &(*T).m[i][j] );
        }
    }
}
*/

/*
 * This function was defined to re-initial the matrix
 */
void reInitial( mat *T, int row, int column )
{
    mFree( T );
    initial( T, row, column );
}


/*
 * The function was defined to get the sum of matrix inputed
 */
void mSum( int n, mat *result, mat *A, ... )
{
    va_list arg_ptr;
    va_start( arg_ptr, A );

    int    row    = (*A).row;
    int    column    = (*A).col;

    reInitial( result, row, column );
    int i, j;


    /*
     * The code was used to check the dimension of matrixs.
     */
    if ( (*A).row != row || (*A).col != column )
    {
        printf( "Matrixs don't have the same dimension.\n" );
        char *p;
        p    = NULL;
        p[0]    = 0;
    }
    for ( i = 1; i < n; ++i )
    {
        mat *M;
        M = va_arg( arg_ptr, mat * );
        if ( (*M).row != row || (*M).col != column )
        {
            printf( "Matrixs don't have the same dimension.\n" );
            char *p;
            p    = NULL;
            p[0]    = 0;
        }
    }


    /*
     * Sum all the matrix and change value of result.
     */
    for ( i = 0; i < row; ++i )
    {
        for ( j = 0; j < row; ++j )
        {
            (*result).m[i][j] += (*A).m[i][j];
        }
    }
    va_start( arg_ptr, A );
    mat *M;
    for ( int k = 1; k < n; ++k )
    {
        M = va_arg( arg_ptr, mat * );
        for ( i = 0; i < row; ++i )
        {
            for ( j = 0; j < column; ++j )
            {
                (*result).m[i][j] += (*M).m[i][j];
            }
        }
    }
}


/*
 * This function was designed to get the product of the input matrix
 */
void mProduct( mat *result, mat *A, mat *B )
{
    if ( (*A).col != (*B).row )
    {
        printf( "Matrixs don't have the required dimension.\n" );
        char *p;
        p    = NULL;
        p[0]    = 0;
    }else {
        int m_len = (*A).col;
        mFree( result );
        zeros( result, (*A).row, (*B).col );
        for ( int i = 0; i < (*A).row; ++i )
        {
            for ( int j = 0; j < (*B).col; ++j )
            {
                for ( int k = 0; k < m_len; ++k )
                {
                    (*result).m[i][j] += (*A).m[i][k] * (*B).m[k][j];
                }
            }
        }
    }
}


/*
 * This function was used to get the upper triangle form of a matrix.
 */
int upTra( mat *T )
{
    double _Complex    temp;
    int    exp = 1;
    for ( int i = 0; i < (*T).row; ++i )
    {
        if ( (*T).m[i][i] != 0 )
        {
            for ( int j = i + 1; j < (*T).row; ++j )
            {
                if ( (*T).m[j][i] != 0 )
                {
                    double _Complex num = (*T).m[j][i] / (*T).m[i][i];
                    for ( int k = i; k < (*T).col; ++k )
                    {
                        (*T).m[j][k] -= num * (*T).m[i][k];
                    }
                }
            }
        }else {
            _Bool flag = 0;
            for ( int check_1 = i + 1; check_1 < (*T).row && flag != 1; ++check_1 )
            {
                if ( (*T).m[check_1][i] != 0 )
                {
                    flag    = 1;
                    exp    = exp * (-1);
                    for ( int check_2 = i; check_2 < (*T).col; ++check_2 )
                    {
                        temp                = (*T).m[check_1][check_2];
                        (*T).m[check_1][check_2]    = (*T).m[i][check_2];
                        (*T).m[i][check_2]        = temp;
                    }
                    for ( int j = i + 1; j < (*T).row; ++j )
                    {
                        if ( (*T).m[j][i] != 0 )
                        {
                            double _Complex num = (*T).m[j][i] / (*T).m[i][i];
                            for ( int k = i; k < (*T).col; ++k )
                            {
                                (*T).m[j][k] -= num * (*T).m[i][k];
                            }
                        }
                    }
                }
            }
        }
    }
    return(exp);
}


/*
 * This function was used to get the determinant of a matrix
 */
double _Complex det( mat *T )
{
    if ( (*T).row != (*T).col )
    {
        printf( "The matrix can't calculate determinant, check the dimension!" );
        char *p;
        p    = NULL;
        p[0]    = 0;
        return(0);
    }else {
        mat Temp;
        initial( &Temp, (*T).row, (*T).col );
        for ( int i = 0; i < (*T).row; ++i )
        {
            for ( int j = 0; j < (*T).col; ++j )
            {
                Temp.m[i][j] = (*T).m[i][j];
            }
        }
        int    exp    = upTra( &Temp );
        double _Complex    det    = 1;
        for ( int i = 0; i < (*T).row; ++i )
        {
            det = det * Temp.m[i][i];
        }
        mFree( &Temp );
        return(exp * det);
    }
}


/*
 * print a matrix
 */
/*
 Function temporarily unavilavle since complex number added.
void mPrint( mat *T )
{
    int    row    = (*T).row;
    int    column    = (*T).col;
    for ( int i = 0; i < row; ++i )
    {
        for ( int j = 0; j < column; ++j )
        {
            if ( j == column - 1 )
            {
                printf( "%8.5lf\n", (*T).m[i][j] );
            }else {
                printf( "%8.5lf\t", (*T).m[i][j] );
            }
        }
    }
    printf( "\n" );
}
 */


/*
 * This function was designed to joint matrix, the input is an array consists of pointers to mat.
 */
void mJoint( int row, int column, mat *result, mat *T[row][column] )
{
    /*
     * 'row' is the row dimension of partitioned matrix
     * 'column' is the column of partioned matrix
     */
    int    row_result    = 0;
    int    column_result    = 0;
    for ( int i = 0; i < column; ++i )
    {
        int k = (T[0][i])->col;
        for ( int j = 0; j < row; ++j )
        {
            if ( (T[j][i])->col != k )
            {
                printf( "Matrixs can't joint!" );
                char *p = NULL;
                p[0] = 0;
            }
        }
        column_result += k;
    }
    for ( int i = 0; i < row; ++i )
    {
        int k = (T[i][0])->row;
        for ( int j = 0; j < column; ++j )
        {
            if ( (T[i][j])->row != k )
            {
                printf( "Matrixs can't joint!" );
                char *p = NULL;
                p[0] = 0;
            }
        }
        row_result += k;
    }

    reInitial( result, row_result, column_result );


    /*
     * Join matrixs into result
     */
    int col_result = 0;
    for ( int j = 0; j < column; ++j )
    {
        int ro_result = 0;
        for ( int i = 0; i < row; ++i )
        {
            for ( int k = 0; k < (T[i][j])->row; ++k )
            {
                for ( int m = 0; m < T[i][j]->col; ++m )
                {
                    (*result).m[k + ro_result][m + col_result] = (*T[i][j]).m[k][m];
                }
            }
            ro_result += (T[i][j])->row;
        }
        col_result += (T[0][j])->col;
    }
}


/*
 * Find the submatrix of *input from row 'row1' to `row2`, `column1` to `column2`, and save to *result
 */
void subMat( mat *result, mat *input, int row1, int row2, int column1, int column2 )
{
    reInitial( result, row2 - row1 + 1, column2 - column1 + 1 );
    for ( int i = 0; i < row2 - row1 + 1; ++i )
    {
        for ( int j = 0; j < column2 - column1 + 1; ++j )
        {
            (*result).m[i][j] = (*input).m[i + row1 - 1][j + column1 - 1];
        }
    }
}


/*
 * The function was designed to get the inverse of the matrix.
 */
char mInverse( mat *result, mat * input1 )
{
    double _Complex tras; /* An variant saving temporary numbers*/
    if ( (*input1).row != (*input1).col )
    {
        printf( "The matrix don't have inverse!" );
        char *p = NULL;
        p[0] = 0;
    }
    char    ifInvert    = 1;
    int    row        = (*input1).row;
    int    column        = (*input1).col;
    mat    Id;              /*Identity Matrix*/
    zeros( &Id, row, row );  /*Initial the identity matrix and set zero*/

    for ( int i = 0; i < row; ++i )
    {
        Id.m[i][i] = 1;
    }

    mat temp;
    initial( &temp, row, 2 * row );
    mat *arr[1][2];
    arr[0][0]    = input1;
    arr[0][1]    = &Id;

    mJoint( 1, 2, &temp, arr );


    upTra( &temp );
    for ( int i = 0; i < row; ++i )
    {
        if ( temp.m[i][i] == 0 )
        {
            ifInvert = 0;
            break;
        }else  {
            tras = temp.m[i][i];
            for ( int k = 0; k < 2 * row; ++k )
            {
                temp.m[i][k] = temp.m[i][k] / tras;
            }
        }
    }
    if ( ifInvert == 0 )
    {
        mFree( result );
        zeros( result, row, row );
    }else  {
        for ( int i = column - 1; i > 0; --i )
        {
            for ( int j = 0; j < i; ++j )
            {
                tras = temp.m[j][i];
                for ( int k = i; k < 2 * row; ++k )
                {
                    temp.m[j][k] -= temp.m[i][k] * tras;
                }
            }
        }
        reInitial( result, row, row ); /*Initial the result matrix*/
        subMat( result, &temp, 1, row, column + 1, 2 * column );
    }

    mFree( &temp );
    mFree( &Id );
    return(ifInvert);
}


/*
 This function is used to get the transpose of T and save the result to *result
 */
void mTranspose(mat *result,mat *T){
    reInitial(result, (*T).col, (*T).row);
    for (int i = 0; i < (*T).col; ++i){
        for (int j = 0; j < (*T).row; ++j){
            (*result).m[j][i] = (*T).m[i][j];
        }
    }
}


/*
 This function is used to equal two matrixs
 */
void mEqual(mat * result,mat * T){
    reInitial(result, (*T).row,(*T).col);
    for (int i = 0; i < (*T).row; ++i){
        for (int j = 0; j < (*T).col; ++j){
            (*result).m[i][j] = (*T).m[i][j];
        }
    }
}
