#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{  
    int j,m=101;
    float time;
    
    printf("Enter the value of t\n");
    scanf("%f",&time);
    
    double Re,U,H,dt,dy,t,gamma,u[m+1],un[m+1];
    Re = 100;
    U = 1.0;
    H = 1.0;
    dt = 0.005;
    dy = H/(m-1);
    gamma = dt/(Re*dy*dy) ;
    
    FILE *fp = fopen("ftcs.dat","w");

    fprintf(fp,"TITLE = VELOCITY\nVARIABLES = \"u\", \"y\"\nZone T = \"BLOCK1\", j=%d, F=POINT\n",m);
    
    //Boundary condition
    un[1] = 0.0;
    un[m] = U;
    for(j=2;j<m;j++)  
    {  
        un[j] = 0.0;
    }
    //FTCS
    int iteration =0;
    
        for(t=0.0;t<=time;t=t+dt)
        {
   
        for(j=1;j<=m;j++)
        {
            u[j] = un[j] ;      //updating u 
        }
        
        for(j=2;j<m;j++)
        {
            un[j] = (1-2*gamma)*u[j] + gamma*(u[j+1] + u[j-1]) ;
        }
       
        iteration++;
        }
         for(j=1;j<=m;j++)
        {
            fprintf(fp,"%lf\t\t%.3lf\n",un[j],(j-1)*dy);
        }
        fclose(fp);
    printf("Number of itertion = %d\n",iteration);
    printf("Data stored in ftcs.dat");
    return 0;
}

