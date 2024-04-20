#include <math.h>
#include <iostream>
# include <fstream>
# include <stdlib.h>

using namespace std;
const int N=20;

void iniciar(double u[N][N], double v[N][N], int N){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
        u[i][j]=0;
        v[i][j]=0;
        }
    }
}

void copiar(double a[N][N], double b[N][N], int N){
    for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
            a[i][j]=b[i][j];
            }
        }
}

void Thomas(double x[N], double a[N], double b[N],double c[N], double d[N])
{
    double w;

    for(int i=1; i<N; i++){
        w=a[i]/b[i-1];
        b[i]=b[i]-w*c[i-1];
        d[i]=d[i]-w*d[i-1];
    }
    
    x[N-1]=d[N-1]/b[N-1];
    for(int i=N-2;i>-1;i--){
        x[i]=(d[i]-c[i]*x[i+1])/b[i];
    }
}

void Sist_eq(double x[N], double a, double b, double d[N])
{
    double w;

    w=a/b;
    b=b-w*a;
    for(int i=1; i<N; i++){
        d[i]=d[i]-w*d[i-1];
    }
    
    x[N-1]=d[N-1]/b;
    for(int i=N-2;i>-1;i--){
        x[i]=(d[i]-a*x[i+1])/b;
    }
}

void ADI(double u[N][N], double v[N][N],double t,double l,double D,double C1,double C2){

    double x1[N][N],x2[N][N],a,b,d[N],aux[N];

    a=-D/pow(l,2);
    b=1/(t/2)+2*D/pow(l,2);

    //Calculo u_n+1
    for(int i=0;i<N;i++){
        d[0]=u[i][0]/(t/2.0)+(C1-(C2+1)*u[i][0]+pow(u[i][0],2)*v[i][0])/pow(l,2);
        d[N-1]=u[i][N-1]/(t/2.0)+(C1-(C2+1)*u[i][N-1]+pow(u[i][0],2)*v[i][N-1])/pow(l,2);
        for(int j=1;j<N-1;j++){
            d[j]=u[i][j]/(t/2.0)+D*(u[i][j+1]-2*u[i][j]+u[i][j-1])/pow(l,2)+(C1-(C2+1)*u[i][j]+pow(u[i][j],2)*v[i][j])/pow(l,2);
        }
        Sist_eq(aux,a,b,d);
        for(int j=0;j<N;j++){
            x1[i][j]=aux[j];
        }
    }

    for(int j=0;j<N;j++){
        d[0]=u[0][j]/(t/2.0)+(C1-(C2+1)*u[0][j]+pow(u[0][j],2)*v[0][j])/pow(l,2);
        d[N-1]=u[N-1][j]/(t/2.0)+(C1-(C2+1)*u[N-1][j]+pow(u[N-1][j],2)*v[N-1][j])/pow(l,2);
        for(int i=1;i<N-1;i++){
            d[i]=x1[i][j]/(t/2.0)+D*(x1[i+1][j]-2*x1[i][j]+x1[i-1][j])/pow(l,2)+(C1-(C2+1)*u[i][j]+pow(u[i][j],2)*v[i][j])/pow(l,2);
        }
        Sist_eq(aux,a,b,d);
        for(int i=0;i<N;i++){
            x1[i][j]=aux[i];
        }
    }

    //Calculo v_n+1
    for(int i=0;i<N;i++){
        d[0]=u[i][0]/(t/2.0)+(C2*u[i][0]-pow(u[i][0],2)*v[i][0])/pow(l,2);
        d[N-1]=u[i][N-1]/(t/2.0)+(C2*u[i][N-1]-pow(u[i][N-1],2)*v[i][N-1])/pow(l,2);
        for(int j=1;j<N-1;j++){
            d[j]=u[i][j]/(t/2.0)+D*(u[i][j+1]-2*u[i][j]+u[i][j-1])/pow(l,2)+(C2*u[i][j]-pow(u[i][j],2)*v[i][j])/pow(l,2);
        }
        Sist_eq(aux,a,b,d);
        for(int j=0;j<N;j++){
            x2[i][j]=aux[j];
        }
    }

    for(int j=0;j<N;j++){
        d[0]=u[0][j]/(t/2.0)+(C2*u[0][j]-pow(u[0][j],2)*v[0][j])/pow(l,2);
        d[N-1]=u[N-1][j]/(t/2.0)+(C2*u[N-1][j]-pow(u[N-1][j],2)*v[N-1][j])/pow(l,2);
        for(int i=1;i<N-1;i++){
            d[i]=x2[i][j]/(t/2.0)+D*(x2[i+1][j]-2*x2[i][j]+x2[i-1][j])/pow(l,2)+(C2*u[i][j]-pow(u[i][j],2)*v[i][j])/pow(l,2);
        }
        Sist_eq(aux,a,b,d);
        for(int i=0;i<N;i++){
            x2[i][j]=aux[i];
        }
    }

        copiar(u,x1,N);
        copiar(v,x2,N);

}

void crear_fichero(string nombre)
{
    ofstream salida(nombre);
    salida << "";
    salida.close();
}

void escribir_datos(double u[][N],double v[][N], int N)
{
    ofstream salida("Chemical_oscillations.dat", ios::app);
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N-1; j++)
        {
            salida << u[i][j] << ",";
        }
        salida << u[i][N-1] <<"\n";
    }
    salida << "\n";
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N-1; j++)
        {
            salida << u[i][j] << ",";
        }
        salida << u[i][N-1] <<"\n";
    }

    salida << "\n";
    salida.close();
}

int main()
{
    double u[N][N],v[N][N];
    double t,l,D,C1,C2;

    t=0.01;
    l=0.01;
    D=1;
    C1=1;
    C2=3;
    crear_fichero("Chemical_oscillations.dat");
    iniciar(u,v,N);

    for(int n=0;n<100;n++){
        escribir_datos(u,v,N);
        ADI(u,v,t,l,D,C1,C2);
    }
    escribir_datos(u,v,N);


    return 0;
}
