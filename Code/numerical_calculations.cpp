#include <math.h>
#include <iostream>
# include <fstream>
# include <stdlib.h>

using namespace std;
const int N=20;
const int M=N-2;

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

void Sist_eq(double x[M], double a, double b, double d[M])
{
    double w,baux[M];
    for(int i=0;i<M;i++){
        baux[i]=b;
    }

    for(int i=1; i<M; i++){
        w=a/baux[i-1];
        baux[i]=baux[i]-w*a;
        d[i]=d[i]-w*d[i-1];

    }
    
    x[M-1]=d[M-1]/baux[M-1];
    for(int i=M-2;i>-1;i--){
        x[i]=(d[i]-a*x[i+1])/baux[i];
    }
}

void ADI(double u[N][N], double v[N][N],double t,double l,double D1,double D2,double C1,double C2){

    double x1[N][N],x2[N][N],y1[N][N],y2[N][N],a,b,d[N-2],aux[N-2];

   

    //Calculo u_n+1

    a=-D1*t/pow(l,2);
    b=(1-2*a);

    for(int i=0;i<N;i++){
        x1[0][i]=0;
        x1[N-1][i]=0;
        x1[i][0]=0;
        x1[i][N-1]=0;

        x2[0][i]=0;
        x2[N-1][i]=0;
        x2[i][0]=0;
        x2[i][N-1]=0;
    }

    for(int j=1;j<N-1;j++){
        for(int i=1;i<N-1;i++){
            d[i-1]=C1*t-a*(u[i][j+1]+u[i][j-1])+(1+2*a-t*(C2+1))*u[i][j]+t*pow(u[i][j],2)*v[i][j];
        }
        Sist_eq(aux,a,b,d);
        for(int i=1;i<N-1;i++){
            x1[i][j]=aux[i-1];
        }
    }

    for(int i=1;i<N-1;i++){
        for(int j=1;j<N-1;j++){
            d[j-1]=C1*t-a*(x1[i+1][j]+x1[i-1][j])+(1+2*a)*x1[i][j]-t*(C2+1)*u[i][j]+t*pow(u[i][j],2)*v[i][j];


        }
        Sist_eq(aux,a,b,d);
        for(int j=1;j<N-1;j++){
            x2[i][j]=aux[i-1];
        }
    }    
   

    //Calculo v_n+1
    a=-D2*t/pow(l,2);
    b=(1-2*a);

    for(int i=0;i<N;i++){
        y1[0][i]=0;
        y1[N-1][i]=0;
        y1[i][0]=0;
        y1[i][N-1]=0;

        y2[0][i]=0;
        y2[N-1][i]=0;
        y2[i][0]=0;
        y2[i][N-1]=0;
    }

    for(int j=1;j<N-1;j++){
        for(int i=1;i<N-1;i++){
            d[i-1]=-a*(v[i][j+1]+v[i][j-1])+(1+2*a)*v[i][j]+t*C2*u[i][j]-t*pow(u[i][j],2)*v[i][j];
        }
        Sist_eq(aux,a,b,d);
        for(int i=1;i<N-1;i++){
            y1[i][j]=aux[i-1];
        }
    }

    for(int i=1;i<N-1;i++){
        for(int j=1;j<N-1;j++){
            d[j-1]=-a*(y1[i+1][j]+y1[i-1][j])+(1+2*a)*y1[i][j]+t*C2*u[i][j]-t*pow(u[i][j],2)*v[i][j];
        }
        Sist_eq(aux,a,b,d);
        for(int j=1;j<N-1;j++){
            y2[i][j]=aux[i-1];
        }
    }


        copiar(u,x2,N);
        copiar(v,y2,N);

}

void crear_fichero(string nombre)
{
    ofstream salida(nombre);
    salida << "";
    salida.close();
}

void escribir_datos(double u[][N],double v[][N], int N)
{
    ofstream salida1("Chemical_oscillations_u.dat", ios::app), salida2("Chemical_oscillations_v.dat", ios::app);
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N-1; j++)
        {
            salida1 << u[i][j] << ",";
        }
        salida1 << u[i][N-1] <<"\n";
    }
    salida1 << "\n";
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N-1; j++)
        {
            salida2 << v[i][j] << ",";
        }
        salida2 << v[i][N-1] <<"\n";
    }

    salida2 << "\n";
    salida1.close();
    salida2.close();
}

int main()
{
    double u[N][N],v[N][N];
    double t,l,D1,D2,C1,C2;

    t=0.01;
    l=0.01;
    D1=1.0/4.0;
    D2=1.0/4.0;
    C1=1;
    C2=1;
    crear_fichero("Chemical_oscillations_u.dat");
    crear_fichero("Chemical_oscillations_v.dat");
    iniciar(u,v,N);

    for(int n=0;n<100;n++){
        escribir_datos(u,v,N);
        ADI(u,v,t,l,D1,D2,C1,C2);
    }
    escribir_datos(u,v,N);
        return 0;



}
