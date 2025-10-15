#include <iostream>
#include <cmath>
#include "lib.h"
using namespace std;

void printlxn(const double *a, int size, int l, int n, int r)
{
    size = size;
    if(n>=r && l>=r) //l<r -> n<r ->n<r and l<r
    {
        
        for(int i =0 ; i < r ; i++){
        cout<<endl;
        for(int j = 0; j < r; j++)
            {

                
                    printf("%10.3e ",a[i*size+j]);
                
            }
            
        }

    }else if (l>=r && n<=r)
    {
        for(int i =0 ; i < r ; i++){
        cout<<endl;
        for(int j = 0; j < n; j++)
        {
             printf("%10.3e ",a[i*size+j]);
        }   
        
        }
    }else if (l<=r && n>=r)
    {
        for(int i =0 ; i < l ; i++){
        cout<<endl;
        for(int j = 0; j < r ; j++)
        {
            printf("%10.3e ",a[i*size+j]);
        }   
        
        }
    }
    
    else
    {
     for(int i =0 ; i < l ; i++){
        cout<<endl;
        for(int j = 0; j < n; j++)
        {
            printf("%10.3e ",a[i*size+j]);
        }   
        
        }
    }

    printf("\n");
}

double f (int s , int n , int i , int j)
{
    if(s == 1) 
        return n - max(i+1,j+1) +1;
    else if (s == 2)
        return max(i+1,j+1);
    else if (s == 3)
        return fabs(i-j);
    else if (s == 4)
        return static_cast<double>(1) / (static_cast<double>(i+1) + static_cast<double>(j+1) - static_cast<double>(1));
    else
        return -1;
}

void init(double *a,  double (*f)(int,int,int,int),  int n, int s )
{
    if(!a || !f)
    {
        printf("a or f = nullptr\n");
        return;
    }

    for(int i = 0; i < n ; i++)
    {
        for(int j=0 ; j < n; j++)
        {
            a[i*n+j] = f(s,n,i,j);
        }
    }
}

double vectornorm(double *a , int n)
{
    double res = 0;

    if(!a)
    {
        printf("nullptr in vector norm\n");
        return -1;
    }
    for(int i = 0; i<n; i++)
    {
        res += fabs(a[i]);
    }

    return res;
}

void mat_x_vector(double *res,double *a, double *b, int n)
{
    

    if(!a || !b)
    {
        return;
    }

    for(int i = 0; i< n; i++)
    {
        double s = 0;
        for(int j = 0; j < n; j++)
        {
            s += a[i*n +j] * b[j];
        }
        res[i] = s;
    }

  
}

double* vectorsub(double *a,double *b, int n)
{
    double *res = new double[n];

    if(!a || !b)
    {
        printf("nullptr in vector subtract\n");
        return 0;
    }

    for(int i = 0 ; i< n ; i++)
    {
        
        res[i] = a[i] - b[i];
    }

    return res;
}



void residuals(double &r1,double &r2,double (*vectornorm)(double* , int ), double* (*mat_x_vector)(double *, double *, int ),double* (*vectorsub)(double *,double *, int ),double *a,double *b,double *x,double *realx,int n)
{
    double *Ax = mat_x_vector(a,x,n);
    double *Ax_b = vectorsub( Ax , b, n);
    double *x_realx = vectorsub( x , realx, n);

    r1 = vectornorm(Ax_b ,  n) / vectornorm(b,n);
    r2 = vectornorm(x_realx ,  n) / vectornorm(realx,n);

    delete []Ax;
    delete []Ax_b;
    delete []x_realx;
    
}






void report(char *title, int task, double r1, double r2 ,double t1,  double t2 ,int s, int n , int m )
{
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2lf T2 = %.2lf S = %d N = %d M = %d\n",
title, task, r1, r2, t1, t2, s, n, m);
}


void matmult(double *res,double *a, double *b, int n, int m,int l) //a^{n*m} * b^{m*l} = res^{n*l}
{
    int i = 0,j=0;
    
   
    if(!a || !b)
    {
        return;
    }

    for(i = 0 ; i < n ; i++)
    {
        for(j = 0; j< l; j++)
        {
            double s = 0;
            for(int k = 0 ; k < m ; k++)
            {
                s += a[i*m+k]*b[k*l+j];
            }
            res[i*l+j] = s;
        }
    }

}

 void get_block(double *a, double *b, int n, int m, int i, int j)
{
    int i1=0, j1=0, k, l, r, h;
    if(m == 0) {
        printf("m == 0 in i = %d , j = %d\n",i,j);
        return;
    }
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    if(j < k) h = m;
    else h = l;


    double *bl = a + i*n*m + j*m; //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < r ; i1++)
    {
        for (j1=0 ; j1 < h ; j1++)
        {
            b[i1*h + j1] = bl[i1*n + j1];
            

        }
    }
    

}



void set_block(double *a, double *b, int n, int m, int i, int j)
{
    int i1=0, j1=0, k, l, r, h;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    if(j < k) h = m;
    else h = l;


    double *bl = a + i*n*m + j*m; //start of block

    for (i1 =0 ; i1 < r ; i1++)
    {
        for (j1=0 ; j1 < h ; j1++)
        {
            bl[i1*n + j1] = b[i1*h + j1];
            

        }
    }
    

}

double normofmatrix(double *a , int size)
{
    if (!a)
    { 
        printf("nullptr in norm of matrix\n");
        return -1;
    }

    double mm = 0;

    for(int j = 0; j< size ; j++)
    {
        double s = 0;
        for(int i = 0 ;i <size; i++)
        {
            s+=fabs(a[i*size+j]);
        }
        if(mm < s) mm = s;
    }

    return mm;

}

double* inverse(double *result,const double* A, int size,double eps)
{
    if(!A)
    {
        printf("nullptr in inverse\n");
        return nullptr;
    }

    // for(int i = 0; i<size*size; i++)
    // {
    //     if(A[i])
    //     {
    //         printf("cant inverse non square matrix\n");
    //         return nullptr;
    //     }
    // }

    double* a = new double[size*size];
    memcpy(a, A, size*size * sizeof(double));


    double* E = new double[size*size];
    int* colsw = new int[size];
    
    // init E
    for(int i = 0; i<size; i++)
    {
        for(int j = 0; j< size; j++)
        {
            if (i!=j) E[i*size+j] = 0;
            else E[i*size +j] = 1;
        }
    }

    //init colsw
    for(int i = 0; i <size ;i++) colsw[i] = i;


    for(int i = 0 ; i< size ; i++)
    {
        double maxEl = fabs(a[i*size+i]);
        int pivot = i;
        for(int j = i+1; j<size;j++)
        {
            if(fabs(a[i*size +j]) > maxEl)
            {
                maxEl = fabs(a[i*size +j]);
                pivot = j;
            }
        }

        
        //swap columns
        if(pivot!=i)
        {
            for(int k = 0; k<size ; k++)
            {
                swap(a[k*size+i],a[k*size+pivot]);
                swap(E[k*size+i],E[k*size+pivot]);
            }
            swap(colsw[i],colsw[pivot]);
        }

        if(fabs(a[i*size+i]) < eps)
        {
            // printf("matrix has no inverse\n");
            delete []E;
            delete []colsw;
            delete []a;
            return nullptr;
        }

        //devide row i 
        // cout<<"A matr before devideing"<<endl;
        // printlxn(a,size,size,size,size);

        double mainEl = a[i*size+i];
        for(int k = 0 ; k<size; k++)
        {
            a[i*size+k] /= mainEl;
            E[i*size+k] /= mainEl;
        }
        // cout<<"A matr after devideing"<<endl;
        // printlxn(a,size,size,size,size);


        for(int k = 0; k< size;k++)
        {
            if(k!=i)
            {
                double factor = a[k*size+i];
                for(int j = 0; j <size;j++)
                {
                    a[k*size+j] -= factor * a[i*size+j];
                    E[k*size+j] -= factor * E[i*size+j];
                }
            }
        }

    }

// cout<<"Last a presentation before swap"<<endl;
// printlxn(a,size,size,size,size);


    // swap back columns
    

    for (int i = 0; i < size*size; ++i) result[i] = 0.0;

    
    
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            result[ colsw[i]*size + colsw[j] ] = E[i*size + j];
// cout<<"Last a presentation after swap"<<endl;

// printlxn(a,size,size,size,size);

    delete []a;
    delete [] colsw;
    delete []E;
    return result;
    
    
}

void blocksize(int i, int j,int n,int m, int &r, int &h)
{
    int k,l;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
        else r = l;

    if(j < k) h = m;
    else h = l;
}


void swap_block_columns(double *a, int n,int m, int i, int j)
{
    for(int p =0 ; p < m ; p++)
                {
                    for(int c = 0; c<n ; c++)
                    {
                        //должны свапнуть столбцы блоков 
                        swap(a[c*n+i*m+p],a[c*n+j*m+p]);
                        
                    }
                }
                
                // printf("swapped blocks %d and %d\n",i,j);
}

void get_vec_block(double *b,double *block,int n, int m, int i)
{
    int i1=0, k, l, r;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    for (i1 =0 ; i1 < r ; i1++)
    {
            block[i1] = b[i*m+i1];
    }
    

}


void mat_mult_sub(double* Result, double* Block_A, double* Block_B, const int row_A, const int col_B, const int col_row) {
    int row_l = row_A % 3;
    int col_l = col_B % 3;
    int row_k = (row_A - row_l) / 3;
    int col_k = (col_B - col_l) / 3;
    double res_00 = 0, res_01 = 0, res_02 = 0;
    double res_10 = 0, res_11 = 0, res_12 = 0;
    double res_20 = 0, res_21 = 0, res_22 = 0;

    for (int b_i = 0; b_i < row_k; b_i++) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            res_20 = 0, res_21 = 0, res_22 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];

                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_12 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_22 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];
            }

            Result[b_i * 3 * col_B + b_j * 3] -= res_00; Result[b_i * 3 * col_B + b_j * 3 + 1] -= res_01; Result[b_i * 3 * col_B + b_j * 3 + 2] -= res_02;
            Result[(b_i * 3 + 1) * col_B + b_j * 3] -= res_10; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 1] -= res_11; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 2] -= res_12;
            Result[(b_i * 3 + 2) * col_B + b_j * 3] -= res_20; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 1] -= res_21; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 2] -= res_22;
        }

        if (col_l != 0) {
            res_00 = 0, res_01 = 0, res_10 = 0;
            res_11 = 0, res_20 = 0, res_21 = 0;

            for (int s = 0; s < col_row; s++) {
                if(col_l > 1) {
                    res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }

                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3];
            }

            Result[b_i * 3 * col_B + col_k * 3] -= res_00; 
            Result[(b_i * 3 + 1) * col_B + col_k * 3] -= res_10;
            Result[(b_i * 3 + 2) * col_B + col_k * 3] -= res_20;

            if(col_l > 1) {
                Result[b_i * 3 * col_B + col_k * 3 + 1] -= res_01; 
                Result[(b_i * 3 + 1) * col_B + col_k * 3 + 1] -= res_11;
                Result[(b_i * 3 + 2) * col_B + col_k * 3 + 1] -= res_21;
            }
        }
            
    }

    if(row_l != 0) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                    res_11 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                    res_12 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 
                }
            }

            Result[row_k * 3 * col_B + b_j * 3] -= res_00; Result[row_k * 3 * col_B + b_j * 3 + 1] -= res_01; Result[row_k * 3 * col_B + b_j * 3 + 2] -= res_02;

            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + b_j * 3] -= res_10; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 1] -= res_11; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 2] -= res_12;
            }
        }

        if(col_l != 0) {
            res_00 = 0, res_01 = 0;
            res_10 = 0, res_11 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                if (col_l > 1) {
                    res_01 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                }
                if (col_l > 1 && row_l > 1) {
                    res_11 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
            }

            Result[row_k * 3 * col_B + col_k * 3] -= res_00;

            if (col_l > 1) {
                Result[row_k * 3 * col_B + col_k * 3 + 1] -= res_01;
            }
            
            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3] -= res_10;
            }

            if (row_l > 1 && col_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3 + 1] -= res_11;
            }
        }     
    }
}

void set_vec_block(double *b,double *block,int n, int m, int i)
{
    int i1=0, k, l, r;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    


    //start of block 

    for (i1 =0 ; i1 < r ; i1++)
    {
            b[i*m+i1]=block[i1] ;
    }
    

}

void get_block_ml(double *a, double *b, int n, int m,int l, int i)
{
    
    int i1=0, j1=0;
   

    double *bl = a + i*n*m + (n-l); //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < m ; i1++)
    {
        for (j1=0 ; j1 < l ; j1++)
        {
            b[i1*l + j1] = bl[i1*n + j1];
            

        }
    }
    

    
}

void set_block_ml(double *a, double *b, int n, int m,int l, int i)
{
    
    int i1=0, j1=0;
   

    double *bl = a + i*n*m + (n-l); //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < m ; i1++)
    {
        for (j1=0 ; j1 < l ; j1++)
        {
             bl[i1*n + j1]=b[i1*l + j1] ;
            

        }
    }
    

    
}

void get_block_lm(double *a, double *b, int n, int m,int l, int j)
{
int i1=0, j1=0;
   

    double *bl = a + (n-l)*n + j*m; //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < l ; i1++)
    {
        for (j1=0 ; j1 < m ; j1++)
        {
            b[i1*m + j1] = bl[i1*n + j1];
            

        }
    }
    
}

void set_block_lm(double *a, double *b, int n, int m,int l, int j)
{
int i1=0, j1=0;
   

    double *bl = a + (n-l)*n + j*m; //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < l ; i1++)
    {
        for (j1=0 ; j1 < m ; j1++)
        {
            bl[i1*n + j1]= b[i1*m + j1];
            

        }
    }
    
}

void matsub(double *res,double *a, double *b, int m,int l)
{
    for(int i = 0 ;i<m;i++)
    {
        for(int j = 0; j<l; j++)
        {
            res[i*l+j] = a[i*l+j] - b[i*l+j];
        }
    }
}

void vec_mult_sub(double* Result, double* A, double* vec, int m) {
    double* temp = new double[m];
    for (int i = 0; i < m; i++)
    {
        double sum = 0;
        for (int j = 0; j < m; j++) {
            
            
                sum += A[i*m + j]*vec[j];

        }
        temp[i] += sum;
            // printf("temp[%d] = %lf\n",i,temp[i]);
        
    }

    for (int i = 0; i < m; i++)
        Result[i] -= temp[i];

    // cout<<"vec_mult:\n";
    // for(int i = 0 ; i < m ; i++)
    // {
    //     cout<<Result[i]<<" ";
    // }
    delete[] temp;
}

void vec_mult_sub_lm(double* Result, double* A, double* vec, int l,int m) {
    double* temp = new double[l];
    for (int i = 0; i < l; i++)
    {
        double sum = 0;
        for (int j = 0; j < m; j++) {
            
            
                sum += A[i*m + j]*vec[j];

        }
        temp[i] += sum;
            // printf("temp[%d] = %lf\n",i,temp[i]);
        
    }

    for (int i = 0; i < l; i++)
        Result[i] -= temp[i];

    // cout<<"vec_mult:\n";
    // for(int i = 0 ; i < m ; i++)
    // {
    //     cout<<Result[i]<<" ";
    // }
    delete[] temp;
}





void undo_block_column_permutation_and_build_x(int n, int m, const int* colsw, const double* x_cur, double* x_out)
{
    int k = n / m;
    int l = n - k * m;
    int block_count = k + (l > 0 ? 1 : 0);

    // Проверка безопасности (по желанию)
    if (!colsw || !x_cur || !x_out) return;

    // sizes per block (block i has size block_sz[i])
    std::vector<int> block_sz(block_count);
    for (int i = 0; i < block_count; ++i)
        block_sz[i] = (i < k ? m : l); // если l == 0, последний блок_sz будет 0 — но тогда block_count == k

    // offsets: offs[0] = 0, offs[1] = block_sz[0], ...
    std::vector<int> offs(block_count + 1, 0);
    for (int i = 0; i < block_count; ++i)
        offs[i + 1] = offs[i] + block_sz[i];

    // Проверка: суммарный размер должен равняться n
    if (offs[block_count] != n) {
        // неконсистентность — сигнал к отладке
        return;
    }

    // Инициализируем выходной вектор (опционально нулями)
    // for (int i = 0; i < n; ++i) x_out[i] = 0.0;

    // Для каждой текущей позиции p: поместимBlock x_cur[p] в x_out в позицию colsw[p]
    for (int p = 0; p < block_count; ++p) {
        int src_off = offs[p];
        int src_len = block_sz[p];

        int dst_block = colsw[p]; // original block index for current position p
        if (dst_block < 0 || dst_block >= block_count) {
            // некорректный индекс — сигнал к отладке
            return;
        }
        int dst_off = offs[dst_block];
        int dst_len = block_sz[dst_block];

        // В идеале src_len == dst_len. Если они отличаются, нужно разбираться:
        if (src_len != dst_len) {
            // Перестановки блоков с разными размерами — опасная ситуация,
            // но её можно обработать как минимум до min_len:
            int min_len = (src_len < dst_len) ? src_len : dst_len;
            for (int t = 0; t < min_len; ++t)
                x_out[dst_off + t] = x_cur[src_off + t];
            // если хотим, можно заполнить оставшиеся элементы нулями либо сигнализировать об ошибке
        } else {
            // Копируем блок
            // memcpy(&x_out[dst_off], &x_cur[src_off], src_len * sizeof(double));
            for (int t = 0; t < src_len; ++t)
                x_out[dst_off + t] = x_cur[src_off + t];
        }
    }
}


void back_substitution(double* A, double* b, double* x, int n, int m)
{
    int k = n / m;
    int l = n - k * m;
    int block_count = k + (l > 0 ? 1 : 0);

    // размер каждого блока
    std::vector<int> block_sz(block_count);
    for (int i = 0; i < block_count; ++i)
        block_sz[i] = (i < k ? m : l);

    // смещения
    std::vector<int> offset(block_count + 1, 0);
    for (int i = 0; i < block_count; ++i)
        offset[i + 1] = offset[i] + block_sz[i];

    // временные вектора
    std::vector<double> tmp(m * m, 0.0);
    std::vector<double> rhs_vec(m, 0.0);

    // обратный ход (снизу вверх)
    for (int bi = block_count - 1; bi >= 0; --bi)
    {
        int bs = block_sz[bi]; // размер блока
        int row_off = offset[bi];

        // b_i' = b_i - sum_{j>i} A_ij * x_j
        std::vector<double> rhs(bs, 0.0);
        for (int t = 0; t < bs; ++t)
            rhs[t] = b[row_off + t];

        for (int bj = bi + 1; bj < block_count; ++bj)
        {
            int cs = block_sz[bj];
            int col_off = offset[bj];

            // A_ij * x_j
            for (int r = 0; r < bs; ++r)
            {
                double s = 0.0;
                for (int c = 0; c < cs; ++c)
                    s += A[(row_off + r) * n + (col_off + c)] * x[col_off + c];
                rhs[r] -= s;
            }
        }

        // Теперь решаем A_ii * x_i = rhs
        int diag_off = row_off;
        for (int i = 0; i < bs; ++i)
        {
            double s = rhs[i];
            for (int j = 0; j < i; ++j)
                s -= A[(diag_off + i) * n + (diag_off + j)] * x[diag_off + j];
            x[diag_off + i] = s / A[(diag_off + i) * n + (diag_off + i)];
        }
    }
}




void multiplication(double* Result, double* Block_A, double* Block_B, const int row_A, const int col_row,const int col_B)
{
    int row_l = row_A % 3;
    int col_l = col_B % 3;
    int row_k = (row_A - row_l) / 3;
    int col_k = (col_B - col_l) / 3;
    double res_00 = 0, res_01 = 0, res_02 = 0;
    double res_10 = 0, res_11 = 0, res_12 = 0;
    double res_20 = 0, res_21 = 0, res_22 = 0;

    for (int b_i = 0; b_i < row_k; b_i++) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            res_20 = 0, res_21 = 0, res_22 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];

                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_12 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_22 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];
            }

            Result[b_i * 3 * col_B + b_j * 3] = res_00; Result[b_i * 3 * col_B + b_j * 3 + 1] = res_01; Result[b_i * 3 * col_B + b_j * 3 + 2] = res_02;
            Result[(b_i * 3 + 1) * col_B + b_j * 3] = res_10; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 1] = res_11; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 2] = res_12;
            Result[(b_i * 3 + 2) * col_B + b_j * 3] = res_20; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 1] = res_21; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 2] = res_22;
        }

        if (col_l != 0) {
            res_00 = 0, res_01 = 0, res_10 = 0;
            res_11 = 0, res_20 = 0, res_21 = 0;

            for (int s = 0; s < col_row; s++) {
                if(col_l > 1) {
                    res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }

                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3];
            }

            Result[b_i * 3 * col_B + col_k * 3] = res_00; 
            Result[(b_i * 3 + 1) * col_B + col_k * 3] = res_10;
            Result[(b_i * 3 + 2) * col_B + col_k * 3] = res_20;

            if(col_l > 1) {
                Result[b_i * 3 * col_B + col_k * 3 + 1] = res_01; 
                Result[(b_i * 3 + 1) * col_B + col_k * 3 + 1] = res_11;
                Result[(b_i * 3 + 2) * col_B + col_k * 3 + 1] = res_21;
            }
        }
            
    }

    if(row_l != 0) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                    res_11 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                    res_12 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 
                }
            }

            Result[row_k * 3 * col_B + b_j * 3] = res_00; Result[row_k * 3 * col_B + b_j * 3 + 1] = res_01; Result[row_k * 3 * col_B + b_j * 3 + 2] = res_02;

            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + b_j * 3] = res_10; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 1] = res_11; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 2] = res_12;
            }
        }

        if(col_l != 0) {
            res_00 = 0, res_01 = 0;
            res_10 = 0, res_11 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                if (col_l > 1) {
                    res_01 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                }
                if (col_l > 1 && row_l > 1) {
                    res_11 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
            }

            Result[row_k * 3 * col_B + col_k * 3] = res_00;

            if (col_l > 1) {
                Result[row_k * 3 * col_B + col_k * 3 + 1] = res_01;
            }
            
            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3] = res_10;
            }

            if (row_l > 1 && col_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3 + 1] = res_11;
            }
        }     
    }      
}



int solution(int n, int m, double *a, double *b, double *x,
    double *block_mm, double *block_ml, double *block_ll,double *invblock_mm, double *diaginvblock_mm, 
    double *invblock_ll,double *diagblock_mm,
    int *colsw,double *vecb_m,double *vecb_l,
    double *tmpblock_mm,double *tmpblock_ml,double *tmpblock_ml1,double *tmpblock_ll,double *tmpvecb_m,double *tmpvecb_l,const double eps)
{
    
    m=m;
    a=a;
    b=b;
    block_ml=block_ml;
    diaginvblock_mm=diaginvblock_mm;
    diagblock_mm=diagblock_mm;
    vecb_l=vecb_l;
    tmpblock_ml1 = tmpblock_ml1;

    

    

    int  k, l/*, r, h*/;
    // int i = 1,
    // j=0;

    k = n/m; l = n - m*k ;
    
    int is_l = (l == 0) ? 0:1; 
    
    
    for(int c =0; c < k  ;c++) colsw[c]=c;

    for(int i = 0 ; i < k + is_l; i++)
    {   
        double minNorm = 1e16;
        int mainBlock = i;

        if(i != k)
        {
            for(int j = i ; j< k ; j++)
            {

            get_block(a,block_mm,n,m,i,j);

            // printf("Block[%d,%d]\n",i,j);
            // printlxn(block_mm,m,m,m,m);

            inverse(invblock_mm,block_mm,m,eps);
            // printlxn(invblock_mm,m,m,m,m);

            if(invblock_mm)
                {
                //     cout<<"inverse "<<i<<" "<<j<<" with norm = "<<normofmatrix(invblock_mm,m)<<endl;

                // printlxn(invblock_mm,m,m,m,m);

                if(normofmatrix(invblock_mm,m) < minNorm) 
                    {
                        minNorm = normofmatrix(invblock_mm,m);
                        mainBlock = j;
                    }
                }

            }

        }else{
            get_block(a,block_ll,n,m,k,k);
            inverse(invblock_ll,block_ll,l,eps);
            if(!invblock_ll)
            {
                printf("Block [%d,%d] has no inverse\n",k,k);
                return -1;
            }
            minNorm = normofmatrix(invblock_ll,l);
        }

        if(fabs(minNorm - 1e16) < eps)
        {
            printf("No inverse matrix in row %d\n",i);
            return -1;
        }

        if(mainBlock != i)
            {
                swap_block_columns(a,n,m,i,mainBlock);
                swap(colsw[i],colsw[mainBlock]);
                cout<<"swapped "<< i<<" "<<mainBlock<<"in row "<<i<<endl;
            }
        
        
        // printlxn(a,n,n,n,n);
        // cout<<"TEST1"<<endl;
        if(i<k)
        {
            get_block(a,diagblock_mm,n,m,i,i);
            
            if(!(inverse(diaginvblock_mm,diagblock_mm,m,eps)))
            {
                        printf("block has no inverse\n");
                        return -1;
                    }

            get_vec_block(b,vecb_m,n,m,i);
            mat_x_vector(tmpvecb_m,diaginvblock_mm,vecb_m,m);// double *resvec = mat_x_vector(diaginvblock_mm,vecb_m,m);
            // cout<<"tmpvecb_m : "<<endl;
            // printlxn(tmpvecb_m,m,1,m,m);    
            set_vec_block(b,tmpvecb_m,n,m,i);

            for(int j = i ; j < k ; j++) //mb try j = i
            {
                get_block(a,block_mm,n,m,i,j);
                
               multiplication(tmpblock_mm,diaginvblock_mm,block_mm,m,m,m);// matmult(tmpblock_mm,diaginvblock_mm,block_mm,m,m,m);// double *resmult = matmult(diaginvblock_mm,block_mm,m,m,m)

                set_block(a,tmpblock_mm,n,m,i,j);
                
                            if (!block_mm || !vecb_m || !invblock_mm) {
                fprintf(stderr, "Error: temporary buffers not initialized!\n");
                return -1;
            }
            }
            // printlxn(a,n,n,n,n);
            if(is_l != 0)
            {
                get_block_ml(a,block_ml,n,m,l,i);
                multiplication(tmpblock_ml,diaginvblock_mm,block_ml,m,m,l);// matmult(tmpblock_ml,diaginvblock_mm,block_ml,m,m,l);
                set_block_ml(a,tmpblock_ml,n,m,l,i);
            }
            
        }else
            {   
                // printlxn(a,n,n,n,n);
                // printlxn(b,n,1,n,n);
                get_block(a,block_ll,n,m,i,i);
                get_vec_block(b,vecb_l,n,m,i);

                // printlxn(block_ll,l,l,l,n);

                // cout<<"vecb_l:"<<endl;
                // printlxn(vecb_l,l,1,l,n);

                if(!(inverse(invblock_ll,block_ll,l,eps)))
                    {
                        printf("block has no inverse\n");
                        return -1;
                    }

                
                // printlxn(invblock_ll,l,l,l,n);

                multiplication(tmpblock_ll,invblock_ll,block_ll,l,l,l);// matmult(tmpblock_ll,invblock_ll,block_ll,l,l,l);

                mat_x_vector(tmpvecb_l,invblock_ll,vecb_l,l);
                
                // printlxn(tmpvecb_l,l,1,l,n);

                set_block(a,tmpblock_ll,n,m,i,i);
                set_vec_block(b,tmpvecb_l,n,m,i);
                // cout<<"WE ARE IN i = k"<<endl;
            }

            // printlxn(a,n,n,n,n);
            // printlxn(b,n,1,n,n);
            //начинаем обнулять столбцы
            
        for(int r = i+1 ; r < k + is_l ; r++)
        {
            // cout<<"TEST "<<r<<endl;
            if(r < k)
            {
                get_block(a,block_mm,n,m,r,i);
                get_block(a,tmpblock_mm,n,m,r,i);
                memset(tmpblock_mm,0, m*m*sizeof(double));
                set_block(a,tmpblock_mm,n,m,r,i);

                // not in i for
                get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
                get_vec_block(b,tmpvecb_m,n,m,r);
                vec_mult_sub(tmpvecb_m,block_mm,vecb_m,m);
                set_vec_block(b,tmpvecb_m,n,m,r);

                // cout<<"tmpvecb_m in subtract i= "<<i<<" r="<<r<<endl;
                // printlxn(tmpvecb_m,m,1,m,m);
                // printlxn(b,n,1,n,n);

                for (int j = i + 1; j < k; j++) {
                    get_block(a,invblock_mm,n,m,i,j);
                    get_block(a,diagblock_mm,n,m,r,j);
                    mat_mult_sub(diagblock_mm,block_mm,invblock_mm,m,m,m);
                    // get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
                    // get_vec_block(b,tmpvecb_m,n,m,r);
                    // vec_mult_sub(tmpvecb_m,block_mm,vecb_m,m);
                    // set_vec_block(b,tmpvecb_m,n,m,r);
                    set_block(a,diagblock_mm,n,m,r,j);
                }

                if (is_l!= 0) {
                get_block_ml(a,tmpblock_ml,n,m,l,i);
                get_block_ml(a,tmpblock_ml1,n,m,l,r);
                mat_mult_sub(tmpblock_ml1,block_mm,tmpblock_ml,m,l,m);
                set_block_ml(a,tmpblock_ml1,n,m,l,r);
                }
            }else
            {
                // printlxn(a,n,n,n,n);
               get_block_lm(a, block_ml, n, m, l, i);

            //    printf("block_lm in col %d\n",i);
            //    printlxn(block_ml,m,l,m,n);

               get_block_lm(a, tmpblock_ml, n, m, l, i);
               memset(tmpblock_ml,0,m*l*sizeof(double));
               set_block_lm(a, tmpblock_ml, n, m, l, i);

               get_vec_block(b,vecb_m,n,m,i);// get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
            //    cout<<"vecb_m in subtract i= "<<i<<" r="<<r<<endl;
            //     printlxn(vecb_m,m,1,m,m);
               get_vec_block(b,tmpvecb_l,n,m,r);// get_vec_block(b,tmpvecb_m,n,m,r);
               vec_mult_sub_lm(tmpvecb_l,block_ml,vecb_m,l,m);// vec_mult_sub(tmpvecb_m,block_mm,vecb_m,m);
               set_vec_block(b,tmpvecb_l,n,m,r);  // set_vec_block(b,tmpvecb_m,n,m,r);


                for(int j = i + 1; j < k; j++) {
                get_block(a,tmpblock_mm,n,m,i,j);
                get_block_lm(a, tmpblock_ml, n, m, l, j);

                // cout<<"tmpblock_ml in col "<<j<<endl;
                // printlxn(tmpblock_ml,m,l,m,m);

                mat_mult_sub(tmpblock_ml,block_ml,tmpblock_mm,l,m,m);

                
                // get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
                // get_vec_block(b,tmpvecb_m,n,m,r);//
                vec_mult_sub(tmpvecb_m,block_ml,vecb_m,m);//
                set_block_lm(a, tmpblock_ml, n, m, l, j);
                // set_vec_block(b,tmpvecb_m,n,m,r);//set_vec...
                }

                if (is_l != 0) {
                    get_block_ml(a,tmpblock_ml,n,m,l,i);
                    get_block(a,tmpblock_ll,n,m,k,k);
                    mat_mult_sub(tmpblock_ll,block_ml,tmpblock_ml,l,l,m);
                    set_block(a,tmpblock_ll,n,m,k,k);
                }

            }
            
        }

        // cout<<"LAST PRINT"<<endl;
        // printlxn(a,n,n,n,n);
        // printlxn(b,n,1,n,n);
         
        
    }
       
    
    //начало обратного хода
        

    // Выполняем обратный ход для получения решения в текущем порядке столбцов
    back_substitution(a, b, x, n, m);

    // Формируем полный массив перестановок блок‑столбцов (включая последний l-блок)
    int block_count = k + is_l;
    std::vector<int> colsw_full(block_count, 0);
    for (int i = 0; i < k; ++i) colsw_full[i] = colsw[i];
    if (is_l) colsw_full[k] = k;

    // Переносим решение в исходный порядок неизвестных
    std::vector<double> x_out(n, 0.0);
    undo_block_column_permutation_and_build_x(n, m, colsw_full.data(), x, x_out.data());
    for (int i = 0; i < n; ++i) x[i] = x_out[i];

    // printf("a[%d,%d] = %lf, b[%d] = %lf\n",k,k,a[k*n+k],k,b[k]);
    // for(int i = n-1 ; i >= 0 ;i--)
    // {
    //     for(int j = 0; j < i ; j++)
    //     {
    //         // printf("a[%d,%d] = %lf, b[%d] = %lf,a[%d,%d] = %lf\n",j,i,a[j*n+i],j,b[j],i,i,a[i*n+i]);
    //         double dd = a[i*n+i]*a[j*n+i];
    //         a[j*n+i] -= a[i*n+i]*dd;
    //         b[j] -= dd*b[i];
    //         // printf("a[%d,%d] = %lf, b[%d] = %lf,a[%d,%d] = %lf\n",j,i,a[j*n+i],j,b[j],i,i,a[i*n+i]);
    //     }
    //     // break;
    //     // return 0;
    // }



    return 0;

}


