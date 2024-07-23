#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#define G (6.67430*pow(10,-11))
#define d0 (149597870700)
#define M0 (1.99*pow(10,30))
#define t0 (sqrt(pow(d0,3)/(G*M0)))
#define N 780
#define K 4  //Aquest K representarà el nombre de cossos, si es vol afegir més s'haurà de canviar el número
double rnorm(double r);
double mnorm(double m);
double tnorm(double t);
double vnorm(double v);




double F(double x,double y,double z,double xj, double yj, double zj, double mj);



double A[K][6][N+1];

int main(){
    int t=86400; // t serà els segons de cada interval de temps (es a dir tindrem N intervals de t segons des del dia 0 fins el 780)    
    double dt=tnorm(t);
    double m[K];
    
    m[0]=mnorm(1.99*pow(10,30));
    m[1]=mnorm(5.97*pow(10,24));
    m[2]=mnorm(6.42*pow(10,23));
    m[3]=mnorm(1.90*pow(10,27));
    //m[4]=mnorm(); Si s'afageix un planeta cal especificar la massa dins la funció mnorm

    A[0][0][0]=rnorm(-1.218*pow(10,9));
    A[0][1][0]=rnorm(6.361*pow(10,8));
    A[0][2][0]=vnorm(-8.131);
    A[0][3][0]=vnorm(-13.68);
    A[0][4][0]=rnorm(2.326*pow(10,7));
    A[0][5][0]=vnorm(3.036*pow(10,-1));
    
    A[1][0][0]=rnorm(1.471*pow(10,11));
    A[1][1][0]=rnorm(-2.514*pow(10,10));
    A[1][2][0]=vnorm(4.598*pow(10,3));
    A[1][3][0]=vnorm(2.922*pow(10,4));
    A[1][4][0]=rnorm(2.398*pow(10,7));
    A[1][5][0]=vnorm(1.65*pow(10,-2));

    A[2][0][0]=rnorm(-2.471*pow(10,11));
    A[2][1][0]=rnorm(-1.422*pow(10,10));
    A[2][2][0]=vnorm(2.354*pow(10,3));
    A[2][3][0]=vnorm(-2.212*pow(10,4));
    A[2][4][0]=rnorm(5.744*pow(10,9));
    A[2][5][0]=vnorm(-5.211*pow(10,2));

    A[3][0][0]=rnorm(6.427*pow(10,11));
    A[3][1][0]=rnorm(-3.853*pow(10,11));
    A[3][2][0]=vnorm(6.560*pow(10,3));
    A[3][3][0]=vnorm(1.182*pow(10,4));
    A[3][4][0]=rnorm(-1.278*pow(10,10));
    A[3][5][0]=vnorm(-1.958*pow(10,2));

   /* A[4][0][0]=rnorm();
    A[4][1][0]=rnorm();
    A[4][2][0]=vnorm();
    A[4][3][0]=vnorm();
    A[4][4][0]=rnorm();
    A[4][5][0]=vnorm();*/  //Si s'afegeix un planeta cal posar les seves condicions inicials 
  

    for (int j = 0; j < N; j++)// j es el bucle de dies, la k calcula per cada cos i la i fa el sumatori per trobar les K de les velocitats
    {

        double K1[K][6]={0};
        double K2[K][6]={0};
        double K3[K][6]={0};
        double K4[K][6]={0};
        for (int k = 0; k < K; k++)
        {
            K1[k][0]=A[k][2][j];
            K1[k][1]=A[k][3][j];
            K1[k][4]=A[k][5][j];// K1 de la Z es la Vz que es el A[k][5][j]
            for (int i = 0; i < K; i++)
            {
                if (i !=k)
                {
                    K1[k][2]=K1[k][2]+F(A[k][0][j],A[k][1][j],A[k][4][j],A[i][0][j],A[i][1][j],A[i][4][j],m[i]);
                    K1[k][3]=K1[k][3]+F(A[k][1][j],A[k][0][j],A[k][4][j],A[i][1][j],A[i][0][j],A[i][4][j],m[i]);
                    K1[k][5]=K1[k][5]+F(A[k][4][j],A[k][0][j],A[k][1][j],A[i][4][j],A[i][0][j],A[i][1][j],m[i]);
                }
                
            }
        }//Aqui hem calculat TOTES les K1 de TOTS els cossos
        for (int k = 0; k < K; k++)
        {
        
            K2[k][0]=A[k][2][j]+dt/2*K1[k][2];
            K2[k][1]=A[k][3][j]+dt/2*K1[k][3];
            K2[k][4]=A[k][5][j]+dt/2*K1[k][5];
            for (int i = 0; i < K; i++)
            {
                if (i!=k)
                {
                    K2[k][2]=K2[k][2]+F(A[k][0][j]+dt/2*K1[k][0],A[k][1][j]+dt/2*K1[k][1],A[k][4][j]+dt/2*K1[k][4],A[i][0][j]+dt/2*K1[i][0],A[i][1][j]+dt/2*K1[i][1],A[i][4][j]+dt/2*K1[i][4],m[i]);
                    K2[k][3]=K2[k][3]+F(A[k][1][j]+dt/2*K1[k][1],A[k][0][j]+dt/2*K1[k][0],A[k][4][j]+dt/2*K1[k][4],A[i][1][j]+dt/2*K1[i][1],A[i][0][j]+dt/2*K1[i][0],A[i][4][j]+dt/2*K1[i][4],m[i]);
                    K2[k][5]=K2[k][5]+F(A[k][4][j]+dt/2*K1[k][4],A[k][0][j]+dt/2*K1[k][0],A[k][1][j]+dt/2*K1[k][1],A[i][4][j]+dt/2*K1[i][4],A[i][0][j]+dt/2*K1[i][0],A[i][1][j]+dt/2*K1[i][1],m[i]);

                }
                
                
            }
        } //Aqui hem calculat TOTES les K2 de TOTS els cossos a partir de les K1 que ja teniem
        for (int k = 0; k < K; k++)
        {
            K3[k][0]=A[k][2][j]+dt/2*K2[k][2];
            K3[k][1]=A[k][3][j]+dt/2*K2[k][3];
            K3[k][4]=A[k][5][j]+dt/2*K2[k][5];
            for (int i = 0; i < K; i++)
            {
                if (i!=k)
                {
                    K3[k][2]=K3[k][2]+F(A[k][0][j]+dt/2*K2[k][0],A[k][1][j]+dt/2*K2[k][1],A[k][4][j]+dt/2*K2[k][4],A[i][0][j]+dt/2*K2[i][0],A[i][1][j]+dt/2*K2[i][1],A[i][4][j]+dt/2*K2[i][4],m[i]);
                    K3[k][3]=K3[k][3]+F(A[k][1][j]+dt/2*K2[k][1],A[k][0][j]+dt/2*K2[k][0],A[k][4][j]+dt/2*K2[k][4],A[i][1][j]+dt/2*K2[i][1],A[i][0][j]+dt/2*K2[i][0],A[i][4][j]+dt/2*K2[i][4],m[i]);
                    K3[k][5]=K3[k][5]+F(A[k][4][j]+dt/2*K2[k][4],A[k][0][j]+dt/2*K2[k][0],A[k][1][j]+dt/2*K2[k][1],A[i][4][j]+dt/2*K2[i][4],A[i][0][j]+dt/2*K2[i][0],A[i][1][j]+dt/2*K2[i][1],m[i]);

                }
                
                
            }
        }//Aqui hem calculat TOTES les K3 de TOTS els cossos a partir de les K2 que ja teniem
        for (int k = 0; k < K; k++)
        {
            K4[k][0]=A[k][2][j]+dt*K3[k][2];
            K4[k][1]=A[k][3][j]+dt*K3[k][3];
            K4[k][4]=A[k][5][j]+dt*K3[k][5];
            for (int i = 0; i < K; i++)
            {
                if (i!=k)
                {
                    K4[k][2]=K4[k][2]+F(A[k][0][j]+dt*K3[k][0],A[k][1][j]+dt*K3[k][1],A[k][4][j]+dt*K3[k][4],A[i][0][j]+dt*K3[i][0],A[i][1][j]+dt*K3[i][1],A[i][4][j]+dt*K3[i][4],m[i]);
                    K4[k][3]=K4[k][3]+F(A[k][1][j]+dt*K3[k][1],A[k][0][j]+dt*K3[k][0],A[k][4][j]+dt*K3[k][4],A[i][1][j]+dt*K3[i][1],A[i][0][j]+dt*K3[i][0],A[i][4][j]+dt*K3[i][4],m[i]);
                    K4[k][5]=K4[k][5]+F(A[k][4][j]+dt*K3[k][4],A[k][0][j]+dt*K3[k][0],A[k][1][j]+dt*K3[k][1],A[i][4][j]+dt*K3[i][4],A[i][0][j]+dt*K3[i][0],A[i][1][j]+dt*K3[i][1],m[i]);
                }
                
                
            }
        }//Aqui hem calculat TOTES les K4 de TOTS els cossos a partir de les K3 que ja teniem
        // Ara ja tenim TOTES les K (K1,K2,K3,K4) de TOTES les variables(x,y,vx,vy,z,vz) de TOTS els cossos. S'identifiquen com K1[0][0] la K1 de la x del sol
        // Ara per actualitzar la x,y,vx,vy,z,vz de cada cos fem un bucle en f per cada cos
         for (int k = 0; k < K; k++)
        {

            A[k][0][j+1]=A[k][0][j]+dt/6*(K1[k][0]+2*K2[k][0]+2*K3[k][0]+K4[k][0]);
            A[k][1][j+1]=A[k][1][j]+dt/6*(K1[k][1]+2*K2[k][1]+2*K3[k][1]+K4[k][1]);
            A[k][2][j+1]=A[k][2][j]+dt/6*(K1[k][2]+2*K2[k][2]+2*K3[k][2]+K4[k][2]);
            A[k][3][j+1]=A[k][3][j]+dt/6*(K1[k][3]+2*K2[k][3]+2*K3[k][3]+K4[k][3]);
            A[k][4][j+1]=A[k][4][j]+dt/6*(K1[k][4]+2*K2[k][4]+2*K3[k][4]+K4[k][4]);
            A[k][5][j+1]=A[k][5][j]+dt/6*(K1[k][5]+2*K2[k][5]+2*K3[k][5]+K4[k][5]);
        }      
        
        
    }
    
    FILE* output; //generem un fitxer que direm output
    
    output=fopen("Carpeta\\traj(3D_4C_1dia).txt","w");  //en aquest cas 4C perquè hi ha 4 cossos
    // // On posa Carpeta cal canviar la direcció on vulguis que et guardi el fitxer dins l'ordinador perquè el codi funcioni. Gracies.
    
    //fprintf(output,"t xsol ysol zsol vxsol vysol vzsol \n"...);  Aquest serà l'ordre, primer les posicions i despres les velocitats en columnes

    for (int j = 0; j < N+1; j++)
    {
        fprintf(output, "%lf  ", dt*j); 
        for (int k = 0; k <K ; k++)
        {
            fprintf(output, "%lf  " "%lf  " "%lf  " "%lf  " "%lf  " "%lf  ", A[k][0][j], A[k][1][j], A[k][4][j],A[k][2][j],A[k][3][j],A[k][5][j]);
        }
        fprintf(output, "\n");
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

double F(double x,double y,double z,double xj, double yj, double zj, double mj){
    
    return -mj*(x-xj)/(pow(pow((x-xj),2)+pow((y-yj),2)+pow((z-zj),2),1.5));
}




