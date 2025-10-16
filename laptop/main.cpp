#include <iostream>
#include "lib.h"
#include <chrono>


using namespace std;

int main(int argc, char *argv[])
{

    int n, m , r, s, task = 10,_c=0;
    char *filename=nullptr;
    double r1=0, r2=0,el;
    double t1=0, t2=0;


    
    filename=filename;

    if(!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n)==1 && sscanf(argv[2], "%d", &m)==1 && sscanf(argv[3], "%d", &r)==1 && sscanf(argv[4], "%d", &s)==1)) 
    {
        cout<<"Usage : "<<argv[0]<<" <n> "<<" <m> "<<" <r> "<<" <s>"<<endl;
      
        return 0;
    }

    double *a = new double[n*n]; //create matrix a
    double *b = new double[n];  // create vector b
    double *x = new double[n];  // create vector x
    double *realx = new double[n];  // create vector real x

    if(argc == 5 && s>=1 && s<=4)
    {
        init(a,&f,n,s);
    }else if(argc == 5 && (s<1 || s>4))
    {
        printf("bad argument for s, s=1,2,3,4 \n");
        report(argv[0],task,r1,r2,t1,t2,s,n,m); 
        delete []a;
        delete []b;
        delete []x;
        delete []realx;
        return 0;
    }

    if(argc == 6 && s == 0 ) 
    {
        filename = argv[5];
        FILE *file = fopen(filename,"r");
        if(!file){
            
            printf("File %s doesnt exist or wrong file name!\n",filename);
            report(argv[0],task,r1,r2,t1,t2,s,n,m); 
            delete []a;
            delete []b;
            delete []x;
            delete []realx;
            
            return 0;
        }

        while(fscanf(file,"%lf",&el)==1)
        {
            if(_c<n*n) 
            {
                a[_c] = el;
                _c++;
            }else
            {
                printf("Bad scan from file %s\n",filename);
                report(argv[0],task,r1,r2,t1,t2,s,n,m); 
                delete []a;
                delete []b;
                delete []x;
                delete []realx;
                fclose(file);
                return 0;
            }

        
        }
        if(!feof(file) || _c!=n*n){
            printf("Bad file %s\n",filename);
            report(argv[0],task,r1,r2,t1,t2,s,n,m); 
            fclose(file);
            delete []a;
            delete []b;
            delete []x;
            delete []realx;
            return 0;
        }
        fclose(file);
    }
    else if(argc == 6 && s!=0)
    {
        printf("Wrong usage! If s!=0 dont use initialization from file or s == 0 and file name not specified \n");
        report(argv[0],task,r1,r2,t1,t2,s,n,m); 
        delete []a;
        delete []b;
        delete []x;
        delete []realx;
        return 0;
    }

    if(normofmatrix(a,n) < EPS64)
    {
        printf("Norm of matrix A < 1e-64 \n");
        report(argv[0],task,r1,r2,t1,t2,s,n,m); 
        delete []a;
        delete []b;
        delete []x;
        delete []realx;
        return 0;
    }

    //print matrix a

    cout<<"\n MATRIX A :\n";
    printlxn(a,n,n,n,r);
       
    
    //init for b
    cout<<"\nVector b : ";
    for (int i = 0; i < n; i++)
    {   
        double sumbi = 0;
        for(int k = 0; k <(n-1)/2+1 ; k++)
        {
            sumbi+= a[i*n+2*k];
            
            
        }
        
        b[i] = sumbi;
        
    }
    printlxn(b,n,1,n,r);

   

    
    //init for x

    
    for(int i =0 ; i< n; i++)
    {
        realx[i] = (i+1)%2 ;
        
    }

    
    
    

    int k,l;

    k = n/m; l = n - m*k ;





    double *block_mm = new double[m*m];
    double *block_ml = new double[m*l];
    double *block_ll = new double[l*l];
    double *tmpblock_mm = new double[m*m];
    double *tmpblock_ml = new double[m*l];
    double *tmpblock_ml1 = new double[m*l];
    double *tmpblock_ll = new double[l*l];
    double *invblock_mm = new double[m*m];
    double *invblock_ll = new double[l*l];
    double *diagblock_mm = new double[m*m];
    double *diaginvblock_mm = new double[m*m];
    double *vecb_m = new double[m];
    double *vecb_l = new double[l];
    double *tmpvecb_m = new double[m];
    double *tmpvecb_l = new double[l]; 
    int *colsw = new int[k];
    double eps = 1e-15*normofmatrix(a,n);

    auto start_sol= std::chrono::high_resolution_clock::now();

    if(solution(n,m,a,b,x,
        block_mm,block_ml,block_ll,invblock_mm,diaginvblock_mm,
        invblock_ll,diagblock_mm,colsw,vecb_m,vecb_l,
        tmpblock_mm,tmpblock_ml,tmpblock_ml1,tmpblock_ll,tmpvecb_m,tmpvecb_l,eps)<0)
    {
        printf("Method stopped\n");
        r1=-1;r2=-1;
        report(argv[0],task,r1,r2,t1,t2,s,n,m); 

        delete []a;
        delete []b;
        delete []x;
        delete []realx;
        delete []block_mm ;
        delete []block_ml ;
        delete []block_ll ;
        delete []tmpblock_mm ;
        delete []tmpblock_ml ;
        delete []tmpblock_ml1 ;
        delete []tmpblock_ll ;
        delete []invblock_mm ;
        delete []invblock_ll ;
        delete []diagblock_mm ;
        delete []diaginvblock_mm ;
        delete []vecb_m ;
        delete []vecb_l ;
        delete []tmpvecb_m ;
        delete []tmpvecb_l ; 
        delete []colsw ;
        
        return 0;

    }


    

    auto end_sol = std::chrono::high_resolution_clock::now();

    //print vector x
    cout<<"\nSolution vector x : ";
    printlxn(x,n,1,n,r);
    

    auto start_res= std::chrono::high_resolution_clock::now();

    

    //reinit


    if(argc == 5 && s>=1 && s<=4)
    {
        init(a,&f,n,s);
    }

    if(argc == 6 && s == 0 ) 
    {
        filename = argv[5];
        FILE *file = fopen(filename,"r");
        

        while(fscanf(file,"%lf",&el)==1)
        {
            if(_c<n*n) 
            {
                a[_c] = el;
                _c++;
            }

        
        }
        if(!feof(file) || _c!=n*n){
            printf("Bad file\n");
            report(argv[0],task,r1,r2,t1,t2,s,n,m); 
            fclose(file);
            delete []a;
            delete []b;
            delete []x;
            delete []realx;
            return 0;
        }
        fclose(file);
    }
    else if(argc == 6 && s!=0)
    {
        printf("Wrong usage! If s!=0 dont use initialization from file or s == 0 and file name not specified \n");
        report(argv[0],task,r1,r2,t1,t2,s,n,m); 
        delete []a;
        delete []b;
        delete []x;
        delete []realx;
        return 0;
    }

  

    //reinit vector b

    for (int i = 0; i < n; i++)
    {   
        double sumbi = 0;
        for(int k = 0; k <(n-1)/2+1 ; k++)
        {
            sumbi+= a[i*n+2*k];
            
            
        }
        
        b[i] = sumbi;
        
    }

    double *Ax = new double[n];//mat_x_vector(a,x,n);
    double *Ax_b = new double[n];//vectorsub( Ax , b, n);
    double *x_realx = new double[n];//vectorsub( x , realx, n);

    residuals(r1,r2,a,b,x,realx,n,Ax,Ax_b,x_realx);

    delete []Ax;
    delete []Ax_b;
    delete []x_realx;
    

    auto end_res= std::chrono::high_resolution_clock::now();

     

    t1 = chrono::duration<double>(end_sol - start_sol ).count();
    t2 = chrono::duration<double>(end_res - start_res).count();



    report(argv[0],task,r1,r2,t1,t2,s,n,m); 

    delete []a;
        delete []b;
        delete []x;
        delete []realx;
        delete []block_mm ;
        delete []block_ml ;
        delete []block_ll ;
        delete []tmpblock_mm ;
        delete []tmpblock_ml ;
        delete []tmpblock_ml1 ;
        delete []tmpblock_ll ;
        delete []invblock_mm ;
        delete []invblock_ll ;
        delete []diagblock_mm ;
        delete []diaginvblock_mm ;
        delete []vecb_m ;
        delete []vecb_l ;
        delete []tmpvecb_m ;
        delete []tmpvecb_l ; 
        delete []colsw ;
    return 0;
}