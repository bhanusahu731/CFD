#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main()
{
   int j,m=101;
   float time;
    printf("Enter the value of t\n");
    scanf("%f",&time);
   double Re,U,H,dt,dy,u[m+1];
    Re = 100;
    U = 1.0;
    H = 1.0;
    dt = 0.01;
    dy = H/(m-1);
    
   double aj,bj,cj,gamma,d[m],P[m],Q[m];
    
    gamma =  dt/(Re*dy*dy);
    
    aj=-(1+2*gamma);
    bj=gamma;
    cj=gamma;

    FILE *fp = fopen("btcs.dat","w");

    fprintf(fp,"TITLE = VELOCITY\nVARIABLES = \"u\", \"y\"\nZone T = \"BLOCK1\", j=%d, F=POINT\n",m);
    
    //Boundary conditions

    u[1] = 0.0;
    u[m] = U;
    for(j=2;j<m;j++)  
    {  
        u[j] = 0.0;
    }
    
    //BTCS TDMA
    double t;
    int iteration=0;

    for(t=0.0;t<=time;t=t+dt)
    {
    
    d[1]=-u[1];
    d[m]=-u[m];
    P[1]=-bj/aj;
    P[m]=0.0;
    Q[1]=d[1]/aj;
    Q[m]=u[m];

    for(j=2; j<m; j++)
    {
        d[j]=-u[j];

        P[j]=-(bj/(aj+cj*P[j-1]));
            
        Q[j]=(d[j]-cj*Q[j-1])/(aj+cj*P[j-1]);
    }
    
    for(j=(m-1);j>1;j--)
    {
        u[j] = P[j] * u[j+1] + Q[j];
    } 
    iteration++;

    }
    
    for(j=1;j<=m;j++)
    {
        fprintf(fp,"%lf\t\t%.3lf\n",u[j],(j-1)*dy);
    }
    
    printf("Number of iterations = %d\n",iteration);
    printf("Data stored in btcs.dat");
    return 0;
}