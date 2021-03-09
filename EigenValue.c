//
//  EigenValue.c
//  cMat
//
//  Created by 吕文韬 on 2021/3/9.
//

#include <stdio.h>
#include <stdlib.h>
#include "cmat.c"

/*
 The function was defined to get the rank of a matrix input.
 */
int mRank(mat *T){
    int rank;
    int row = (*T).row;
    int column = (*T).col;
    mat Temp;
    mTranspose(&Temp, T);
    if (column >= row){
        rank = row;
        upTra(&Temp);
        for (int i = 0; i < row; ++i){
            if ((Temp).m[i][i] == 0){
                --rank;
            }
        }
    }
    else{
        mEqual(&Temp, T);
        rank = column;
        upTra(&Temp);
        for (int i = 0; i < row; ++i){
            if ((Temp).m[i][i] == 0){
                --rank;
            }
        }
    }
    return rank;
}


/*
 This function is defined to get the eigenvalues of matrix input.
 */


