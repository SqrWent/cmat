/*
 *
 *  cmat.h
 *  cmat
 *
 *  Created by 吕文韬 on 2021/3/9.
 *
 */

#ifndef cmat_h
#define cmat_h

/*Declare matrix structure*/
typedef struct {
    double _Complex **m;
    int        row, col;
}
mat;

/*
 * Declare vector structure
 * By default, it's a column vector.
 */
typedef struct {
    double _Complex *v;
    int        row;
}vector;


/*Initialize a matrix*/
void initial( mat *T, int row, int column );


/*
 * Generalize zero matrix
 */
void zeros( mat *T, int row, int column );


/*
 * Free matrix memory
 */
void mFree( mat *T );


/*
 * Read matrix
 */
/* void mScanf( mat *T ); */


/*
 * This function was defined to re-initial the matrix
 */
void reInitial( mat *T, int row, int column );


/*
 * The fmProductfined to get the sum of matrix inputed
 */
void mSum( int n, mat *result, mat *A, ... );


/*
 * This function was designed to get the product of the input matrix
 */
void mProduct( mat *result, mat *A, mat *B );


/*
 * This function was used to get the upper triangle form of a matrix.
 */
int upTra( mat *T );


/*
 * This function was used to get the determinant of a matrix
 */
double _Complex det( mat *T );


/*
 * print a matrix
 */
void mPrint( mat *T );


/*
 * This function was designed to joint matrix, the input is an array consists of pointers to mat.
 */
void mJoint( int row, int column, mat *result, mat *T[row][column] );


/*
 * Find the submatrix of *input from row 'row1' to `row2`, `column1` to `column2`, and save to *result
 */
void subMat( mat *result, mat *input, int row1, int row2, int column1, int column2 );


/*
 * The function was designed to get the inverse of the matrix.
 */
char mInverse( mat *result, mat * input1 );


/*
 * This function is used to get the transpose of T and save the result to *result
 */
void mTranspose( mat *result, mat *T );


/*
 * This function is used to equal two matrixs
 */
void mEqual( mat * result, mat * T );


/* This function is defined to get the i-th column vector of matrix A and save to result */
void colVector( vector * result, mat * T, int i );


/*
 * The function was defined to get the rank of a matrix input.
 */
int mRank( mat *T );


/* This function is defined to initial a vector */
void vInitial( vector *T, int m );


/* This function is used to make zero vector. */
void vZero( vector *T, int m );

/*Should only cited by other functions within cmat*/
void vStaticZero( vector *T );

/*Should only cited by other functions within cmat*/
void mStaticZero(mat *T/* The input matrix*/);


/* This function is used to free vector memory */
void vFree( vector *T );


/* This function is used to reInitial the vector memory */
void vReInitial( vector *T, int n );


/* This function is defined to get the i-th column vector of matrix A and save to result */
void colVector( vector *result, mat *T, int i );


/* This function is defined to get the inner product of two vectors */
double _Complex vInnerProduct( vector *a, vector *b );


/* This function is defined to get the norm of a vector */
double vNorm( vector *a );


/* This function is defined to get vectors joint */
void vJoint( int column, mat * result, vector *s[column] );

/* This function is defined to get a vector product by a number */
void vNumProduct( vector *a, double num, vector *b );


/* This function is defined to carry Gram-Schmidt procedure on vectors. */
void MGS( mat *Qresult, mat *Rresult, mat *input );

#ifndef _DEFVALUE
#define _DEFVALUE(arg,defvalue) ( (#arg[0]) ? (arg + 0) : defvalue)//The macro defined to set defvalue for args.

#ifndef eigValueMGS
#define eigValueMGS(arg1, arg2) _eigValueMGS(arg1, _DEFVALUE(arg2, 1000))

/* This function is defined get the Eigen Value of matrix via MGS method. */
void _eigValueMGS( double _Complex * eigvalues, mat *T, int times );

#endif /* eigValueMGS */


//This function is defined to get the the Housholder matrix from a to b
void vHousholder(mat *H /*The Householder Matrix*/, vector *a, vector *b);

/* This function is defined to do Housholder transformation to matrix inputed. */
void mHousholder(mat * Qresult /*The unitary matrix*/,
                 mat * Rresult /* The upper triangle matrix*/,
                 mat * input /*The input matrix*/ );

#ifndef eigValueHS
#define eigValueHS(arg1, arg2) _eigValueHS(arg1, _DEFVALUE(arg2, 1000))

/* This function is defined to get the Eigen Value of matrix via Householder method. */
void _eigValueHS(double _Complex * eigvalues/* The pointer to eigenvalues*/,
                mat *T/*The target matrix*/,
                int times/*The iteration times*/ );

#endif /* eigValueHS */

/*This function is defined to do QR decomposition on a matrix via Givens means*/
void mGivens(mat * Qresult/*The unitary matrix*/,
             mat * Rresult/*The tpper triangle matrix*/,
             mat * input/*The input matrix*/);

#ifndef eigValueGVS
#define eigValueGVS(arg1, arg2) _eigValueGVS(arg1, _DEFVALUE(arg2, 1000))

/* This function is defined to get the Eigen Value of matrix via Given method. */
void _eigValueGVS(double _Complex * eigvalues/* The pointer to eigenvalues*/,
                 mat *T/*The target matrix*/,
                 int times/*The iteration times*/ );
#endif /* eigValueGVS */

#endif /* _DEFVALUE */

#endif /* cmat_h */
