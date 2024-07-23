#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#define G (6.67430*pow(10,-11))
#define d0 (149597870700)
#define M0 (1.99*pow(10,30))
#define t0 (sqrt(pow(d0,3)/(G*M0)))
#define N  780
double rnorm(double r);
double mnorm(double m);
double tnorm(double t);
double vnorm(double v);

double m[4];


double F(double x,double y,double xj, double yj, double mj);



double A[4][4][N+1];


int main(){
    int t=84600; // t serà els segons de cada interval de temps (es a dir tindrem N intervals de t segons des del dia 0 fins el 780)
    double dt=tnorm(t);
    double L[4][N+1];

    
    m[0]=mnorm(1.99*pow(10,30));
    m[1]=mnorm(5.97*pow(10,24));
    m[2]=mnorm(6.42*pow(10,23));
    m[3]=mnorm(1.90*pow(10,27));

    A[0][0][0]=rnorm(-1.218*pow(10,9));
    A[0][1][0]=rnorm(6.361*pow(10,8));
    A[0][2][0]=vnorm(-8.131);
    A[0][3][0]=vnorm(-13.68);
    
    A[1][0][0]=rnorm(1.471*pow(10,11));
    A[1][1][0]=rnorm(-2.514*pow(10,10));
    A[1][2][0]=vnorm(4.598*pow(10,3));
    A[1][3][0]=vnorm(2.922*pow(10,4));

    A[2][0][0]=rnorm(-2.471*pow(10,11));
    A[2][1][0]=rnorm(-1.422*pow(10,10));
    A[2][2][0]=vnorm(2.354*pow(10,3));
    A[2][3][0]=vnorm(-2.212*pow(10,4));

    A[3][0][0]=rnorm(6.427*pow(10,11));
    A[3][1][0]=rnorm(-3.853*pow(10,11));
    A[3][2][0]=vnorm(6.560*pow(10,3));
    A[3][3][0]=vnorm(1.182*pow(10,4));
  

    double K1[4][4]={0};
    double K2[4][4]={0};
    double K3[4][4]={0};
    double K4[4][4]={0};

    for (int j = 0; j < N; j++)// j es el bucle de dies, la k calcula per cada cos i la i fa el sumatori per trobar les K de les velocitats
    {   
        for (int k = 0; k < 4; k++)
        {
            for (int i = 0; i < 4; i++)
            {
                K1[k][i]=0;
                K2[k][i]=0;
                K3[k][i]=0;
                K4[k][i]=0;                
            }
        
        }    
        
        for (int k = 0; k < 4; k++)
        {
            K1[k][0]=A[k][2][j];
            K1[k][1]=A[k][3][j];
            for (int i = 0; i < 4; i++)
            {
                if (i !=k)
                {
                    K1[k][2]=K1[k][2]+F(A[k][0][j],A[k][1][j],A[i][0][j],A[i][1][j],m[i]);
                    K1[k][3]=K1[k][3]+F(A[k][1][j],A[k][0][j],A[i][1][j],A[i][0][j],m[i]);
                }
                
            }
        }//Aqui hem calculat TOTES les K1 (K1 de x de y de vx i de vy ) de TOTS els cossos
        for (int d = 0; d < 4; d++)
        {
        
            K2[d][0]=A[d][2][j]+dt/2*K1[d][2];
            K2[d][1]=A[d][3][j]+dt/2*K1[d][3];
            for (int l = 0; l < 4; l++)
            {
                if (l!=d)
                {
                    K2[d][2]=K2[d][2]+F(A[d][0][j]+dt/2*K1[d][0],A[d][1][j]+dt/2*K1[d][1],A[l][0][j]+dt/2*K1[l][0],A[l][1][j]+dt/2*K1[l][1],m[l]);
                    K2[d][3]=K2[d][3]+F(A[d][1][j]+dt/2*K1[d][1],A[d][0][j]+dt/2*K1[d][0],A[l][1][j]+dt/2*K1[l][1],A[l][0][j]+dt/2*K1[l][0],m[l]);
                }
                
                
            }
        } //Aqui hem calculat TOTES les K2 (K2 de x de y de vx i de vy ) de TOTS els cossos a partir de les K1 que ja teniem
        for (int o = 0; o < 4; o++)
        {
            K3[o][0]=A[o][2][j]+dt/2*K2[o][2];
            K3[o][1]=A[o][3][j]+dt/2*K2[o][3];
            for (int n = 0; n < 4; n++)
            {
                if (n!=o)
                {
                    K3[o][2]=K3[o][2]+F(A[o][0][j]+dt/2*K2[o][0],A[o][1][j]+dt/2*K2[o][1],A[n][0][j]+dt/2*K2[n][0],A[n][1][j]+dt/2*K2[n][1],m[n]);
                    K3[o][3]=K3[o][3]+F(A[o][1][j]+dt/2*K2[o][1],A[o][0][j]+dt/2*K2[o][0],A[n][1][j]+dt/2*K2[n][1],A[n][0][j]+dt/2*K2[n][0],m[n]);
                }
                
                
            }
        }//Aqui hem calculat TOTES les K3 (K3 de x de y de vx i de vy ) de TOTS els cossos a partir de les K2 que ja teniem
        for (int g = 0; g < 4; g++)
        {
            K4[g][0]=A[g][2][j]+dt*K3[g][2];
            K4[g][1]=A[g][3][j]+dt*K3[g][3];
            for (int p = 0; p < 4; p++)
            {
                if (p!=g)
                {
                    K4[g][2]=K4[g][2]+F(A[g][0][j]+dt*K3[g][0],A[g][1][j]+dt*K3[g][1],A[p][0][j]+dt*K3[p][0],A[p][1][j]+dt*K3[p][1],m[p]);
                    K4[g][3]=K4[g][3]+F(A[g][1][j]+dt*K3[g][1],A[g][0][j]+dt*K3[g][0],A[p][1][j]+dt*K3[p][1],A[p][0][j]+dt*K3[p][0],m[p]);
                }
                
                
            }
        }//Aqui hem calculat TOTES les K4 (K4 de x de y de vx i de vy ) de TOTS els cossos a partir de les K3 que ja teniem
        // Ara ja tenim TOTES les K (K1,K2,K3,K4) de TOTES les variables(x,y,vx,vy) de TOTS els cossos. S'identifiquen com K1[0][0] la K1 de la x del sol
        // Ara per actualitzar la x,y,vx,vy de cada cos fem un bucle en f per cada cos
         for (int f = 0; f < 4; f++)
        {

            A[f][0][j+1]=A[f][0][j]+dt/6*(K1[f][0]+2*K2[f][0]+2*K3[f][0]+K4[f][0]);
            A[f][1][j+1]=A[f][1][j]+dt/6*(K1[f][1]+2*K2[f][1]+2*K3[f][1]+K4[f][1]);
            A[f][2][j+1]=A[f][2][j]+dt/6*(K1[f][2]+2*K2[f][2]+2*K3[f][2]+K4[f][2]);
            A[f][3][j+1]=A[f][3][j]+dt/6*(K1[f][3]+2*K2[f][3]+2*K3[f][3]+K4[f][3]);
        }      
        
        
    }

    /* Sabem que el moment angluar es una quantitat conservada, 
    i per tant si el calculem a l'inici i al final hauria de ser el mateix.
    El càlcul és el següent */   
    for (int iii=0; iii<4; iii++)
    {
        for (int jjj = 0; jjj < N+1; jjj++)
        {
        L[iii][jjj]=m[iii]*(A[iii][0][jjj]*A[iii][3][jjj]-A[iii][1][jjj]*A[iii][2][jjj]);
        
        }
    }
    
    
    FILE* output; //generem un fitxer que direm output 
    
    
    output=fopen("Carpeta\\traj(1dia-mom-ang).txt","w"); 
    // On posa Carpeta cal canviar la direcció on vulguis que et guardi el fitxer dins l'ordinador perquè el codi funcioni. Gracies.
    
        for (int jjj = 0; jjj < N+1; jjj++)
    {
      fprintf(output, "%lf   " "%lf   " "%lf   " "%lf   " "%lf\n", jjj*dt , fabs((L[0][jjj]-L[0][0])*100/L[0][0]),fabs((L[1][jjj]-L[1][0])*100/L[1][0]),fabs((L[2][jjj]-L[2][0])*100/L[2][0]),fabs((L[3][jjj]-L[3][0])*100/L[3][0]));  
    }
    
    fclose(output);
    
    return 0;
}

double rnorm(double r){
    double rnorm=r/d0;
    return rnorm;
}
double mnorm(double m){
    double mnorm=m/M0;
    return mnorm;
}
double tnorm(double t){
    double tnorm=t/t0;
    return tnorm;
}
double vnorm(double v){
    double vnorm=v*t0/d0;
    return vnorm;
}

double F(double x,double y,double xj, double yj, double mj){
    
    return -mj*(x-xj)/(pow((pow((x-xj),2)+pow((y-yj),2)),1.5));
}

