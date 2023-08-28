
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
   int j,m=101;
   double Re,U,H,dt,dy,u[m],un[m];
    Re = 100;
    U = 1.0;
    H = 1.0;
    dt = 0.01;
    dy = H/(m-1);
    
   double aj,bj,cj,gamma,d[m],P[m],Q[m];
    
    gamma =  dt/(Re*dy*dy);
    
    aj=-(1+gamma);
    bj=gamma/2;
    cj=gamma/2;
    
    //Boundary conditions

    un[1] = 0.0;
    un[m] = U;
    for(j=2;j<m;j++)  
    {  
        un[j] = 0.0;
    }
    u[m]=U;
     for(j=1;j<m;j++)  
    {  
        u[j] = 0.0;
    }
    
    FILE *fp = fopen("cn_steady.dat","w");

    fprintf(fp,"TITLE = VELOCITY\nVARIABLES = \"u\", \"y\"\nZone T = \"BLOCK1\", j=%d, F=POINT\n",m);

    //CRANC NICOLSON TDMA
    d[1]=-u[1];
    d[m]=-u[m];
    P[1]=-bj/aj;
    P[m]=0.0;
    Q[1]=d[1]/aj;
    Q[m]=u[m];

    double error = 0.0;
    int iteration=0;

   do
    {
    for(j=2; j<m; j++)
    {
        d[j]=(gamma-1)*u[j]-(gamma/2)*(u[j+1]+u[j-1]);

        P[j]=-(bj/(aj+(cj*P[j-1])));
            
        Q[j]=(d[j]-(cj*Q[j-1]))/(aj+(cj*P[j-1]));
    }
    
    for(j=(m-1);j>1;j--)
    {
        un[j] = P[j] * un[j+1] + Q[j];
    } 

    error =0.0;
    for(j=1;j<=m;j++)
    {
        error = error + pow((un[j]-u[j]),2);
    }
    error = sqrt(error/m);
    for(j=1;j<=m;j++)
    {
        u[j] = un[j];
    }
    iteration++;
    }while(error > 0.000001);
    
    for(j=1;j<=m;j++)
    {
        fprintf(fp,"%lf\t\t%.3lf\n",un[j],(j-1)*dy);
    }
    
    printf("Number of iterations = %d\n",iteration);
    printf("Data stored in cn_steady.dat");
    return 0;
}

