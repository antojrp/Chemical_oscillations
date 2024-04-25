#include <math.h>
#include <iostream>
# include <fstream>
# include <stdlib.h>

using namespace std;
const int N=100;
const int M=N-2;
const double PI = 3.14159265358979323846264;

void iniciar(double u[N][N], double v[N][N], int N, string mode){
    if(mode == "zeros")
    {
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++)
            {
                u[i][j]=0;
                v[i][j]=0;
            }
        }
    }

    if(mode == "centro")
    {
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++)
            {
                u[i][j]=0;
                v[i][j]=0;
            }
        }

        for(int i=N/2 - 30; i<N/2 + 30; i++){
            for(int j=N/2 - 30; j<N/2 + 30; j++)
            {
                u[i][j]=1;
                v[i][j]=0;
            }
        }
    }

    if(mode == "mitad")
    {
        for(int i=0; i<N/2.0; i++){
            for(int j=0; j<N; j++)
            {
                u[i][j]=1;
                u[N-i-1][j]=0;
                v[i][j]=0;
                v[N-i-1][j]=1;
            }
        }
    }

    if(mode == "seno")
    {
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++)
            {
                u[i][j]= pow(sin( ((i + j) * 2.0 * PI) / N), 2);
                v[i][j]= pow(cos( ((i + j) * 2.0 * PI) / N), 2);
            }
        }
    }
    
    if(mode == "random")
    {
        for (int i = 0; i < N; ++i) 
        {
            for (int j = 0; j < N; ++j) 
            {
                u[i][j] = rand() / (double)(RAND_MAX + 1.0);
                v[i][j] = rand() / (double)(RAND_MAX + 1.0);
            }
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

void cyclic_thomas(double x[N], double a, double b) {
    
    double alpha = a;
    double beta = a;
    double cmod[N], u[N];
    double gamma = -b;
    double m;
    double fact;

    cmod[0] = alpha / (b - gamma);
    u[0] = gamma / (b - gamma);
    x[0] /= (b - gamma);


    for (int i = 1; i + 1 < N; i++) {
        m = 1.0 / (b - a * cmod[i - 1]);
        cmod[i] = a * m;
        u[i] = (0.0f  - a * u[i - 1]) * m;
        x[i] = (x[i] - a * x[i - 1]) * m;
    }

    m = 1.0 / (b - alpha * beta / gamma - beta * cmod[N - 2]);
    u[N - 1] = (alpha    - a * u[N - 2]) * m;
    x[N - 1] = (x[N - 1] - a* x[N - 2]) * m;

    for (int i = N - 2; i >= 0; i--) {
        u[i] -= cmod[i] * u[i + 1];
        x[i] -= cmod[i] * x[i + 1];
    }

    fact = (x[0] + x[N - 1] * beta / gamma) / (1.0 + u[0] + u[N - 1] * beta / gamma);

    for (int i = 0; i < N; i++)
        x[i] -= fact * u[i];
}

void ADI(double u[N][N], double v[N][N],double t,double l,double D1,double D2,double C1,double C2){

    double x1[N][N],x2[N][N],y1[N][N],y2[N][N],a,b,d[N-2],aux[N-2];


    //Calculo u_n+1

    a=-D1*t/pow(l,2);
    b=(1-2*a);

    for(int j=0;j<N;j++){
        for(int i=0;i<N;i++){
            d[i]=C1*t-a*(u[i][(j+1)%N]+u[i][(N+j-1)%N])+(1+2*a-t*(C2+1))*u[i][j]+t*pow(u[i][j],2)*v[i][j];
        }
        cyclic_thomas(d,a,b);
        for(int i=0;i<N;i++){
            x1[i][j]=d[i];
        }
    }

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            d[j]=C1*t-a*(x1[(i+1)%N][j]+x1[(N+i-1)%N][j])+(1+2*a)*x1[i][j]-t*(C2+1)*u[i][j]+t*pow(u[i][j],2)*v[i][j];

        }
        cyclic_thomas(d,a,b);
        for(int j=0;j<N;j++){
            x2[i][j]=d[j];
        }
    }    
   

    //Calculo v_n+1
    a=-D2*t/pow(l,2);
    b=(1-2*a);


    for(int j=0;j<N;j++){
        for(int i=0;i<N;i++){
            d[i]=-a*(v[i][(j+1)%N]+v[i][(N+j-1)%N])+(1+2*a)*v[i][j]+t*C2*u[i][j]-t*pow(u[i][j],2)*v[i][j];
        }
        cyclic_thomas(d,a,b);
        for(int i=0;i<N;i++){
            y1[i][j]=d[i];
        }
    }

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            d[j-1]=-a*(y1[(i+1)%N][j]+y1[(N+i-1)%N][j])+(1+2*a)*y1[i][j]+t*C2*u[i][j]-t*pow(u[i][j],2)*v[i][j];
        }
        cyclic_thomas(d,a,b);
        for(int j=0;j<N;j++){
            y2[i][j]=d[j];
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
    ofstream salida1("Chemical_oscillations_u.txt", ios::app), salida2("Chemical_oscillations_v.txt", ios::app);
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

    // mode = "zeros", "random", "centro", "mitad", "seno"
    string mode = "random";

    t=0.001;
    l=0.01  ;
    D1=0.0/4.0;
    D2=0.0/4.0; 
    C1=30;
    C2=C1 * pow(30, 0.5);
    crear_fichero("Chemical_oscillations_u.txt");
    crear_fichero("Chemical_oscillations_v.txt");
    iniciar(u,v,N,mode);

    for(int n=0;n<500;n++){
        escribir_datos(u,v,N);
        ADI(u,v,t,l,D1,D2,C1,C2);
    }
    escribir_datos(u,v,N);

    return 0;

}
