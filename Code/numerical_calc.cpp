#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <time.h>

using namespace std;
const int N=250;
const int M=N-2;
const float PI = 3.14159265358979323846264;

void iniciar(float u[N][N], float v[N][N], int N, string mode){
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

    if(mode == "esquina")
    {
        for(int i=0; i<N/4; i++){
            for(int j=0; j<N/4; j++)
            {
                u[i][j]=1;
                v[i][j]=0;
            }
        }

        for(int i=N/4; i<N; i++){
            for(int  j=N/4; j<N; j++)
            {
                u[i][j]=0;
                v[i][j]=1;
            }
        }

        for(int i=3*N/4; i<N; i++){
            for(int  j=3*N/4; j<N; j++)
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
    
    if (mode == "random")
    {
        srand(time(NULL));  // Initialize the random seed with the current time
        for (int i = 0; i < N; ++i) 
        {
            for (int j = 0; j < N; ++j) 
            {
                u[i][j] = 1.0 * rand() / (float)RAND_MAX;
                v[i][j] = 1.0 * rand() / (float)RAND_MAX;
            }
        }
    }

    
}

void copiar(float a[N][N], float b[N][N], int N){
    for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
            a[i][j]=b[i][j];
            }
        }
}

void cyclic_thomas(float x[N], float a, float b) {
    
    float alpha = a;
    float beta = a;
    float cmod[N], u[N];
    float gamma = -b;
    float m;
    float fact;

    cmod[0] = alpha / (b - gamma);
    u[0] = gamma / (b - gamma);
    x[0] = x[0] / (b - gamma);


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

    for (int i = 0; i < N; i++){
        x[i] -= fact * u[i];
    }
}

void ADI(float u[N][N], float v[N][N],float t,float l,float D1,float D2,float C1,float C2){

    float x1[N][N],x2[N][N],y1[N][N],y2[N][N],a,b,d[N];


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
            d[j]=-a*(y1[(i+1)%N][j]+y1[(N+i-1)%N][j])+(1+2*a)*y1[i][j]+t*C2*u[i][j]-t*pow(u[i][j],2)*v[i][j];
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

void comparativa_rgb(float u[][N], float v[][N], float p_u[][N], float p_v[][N])
{
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            if(u[i][j] - v[i][j] >= 0)
            {
                p_u[i][j] = (u[i][j] - v[i][j])*(1);
                p_v[i][j] = 0; // (u[i][j] - v[i][j]);
            }
            else
            {
                p_v[i][j] = (v[i][j] - u[i][j])*(1);
                p_u[i][j] = 0; // (v[i][j] - u[i][j]);
            }
        }
    }
}

void escribir_datos(float u[][N],float v[][N], int N, int i)
{
    ofstream salida1("Chemical_oscillations_u"+to_string(i)+".txt", ios::app), salida2("Chemical_oscillations_v"+to_string(i)+".txt", ios::app);
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
    float u[N][N],v[N][N], p_u[N][N],p_v[N][N];
    float t,l,D1,D2,a,b,t0,t1,iteraciones,eta,bc,mu;

    // mode = "zeros", "random", "centro", "mitad", "seno", "esquina"

    for(int i=5; i<7;i++){
        string mode = "random";

        t0=clock();

        iteraciones = 100000;

        t=0.8*pow(10.0, -3);

        D1=5.0;
        D2=40.0;
        eta = pow(D1 / D2, 0.5);
        a=5.0; // rojo
        bc=pow(1+eta*a,2);
        mu=-0.5+i/2.0;
        b=bc*mu+bc; // azul

        l = 3.5 * (5.0 * a * eta + 7*pow(a*eta, 2) - 3 - 3*pow(a*eta, 3) ) / (pow(a,3)* eta * (1 + a *eta));

        cout << "The value of i is: " << i << endl;
        cout << "The value of mu is: " << mu << endl;
        
        crear_fichero("Chemical_oscillations_u"+to_string(i)+".txt");
        crear_fichero("Chemical_oscillations_v"+to_string(i)+".txt");
        iniciar(u,v,N,mode);

        cout << "Progress: " << endl;
        cout << "░░░░░░░░░░" << endl;
        for(int n=0;n<iteraciones;n++){
            if( n%400 == 0)
                {
                escribir_datos(u,v,N,i);
                }
            if (( n % 10000 == 0 ) && ( n != 0))
                {
                    cout << "▐";
                }
            ADI(u,v,t,l,D1,D2,a,b);
        }
        escribir_datos(u,v,N,i);;

        t1=clock();

        cout << "El programa ha hecho " << iteraciones << " iteraciones para una red " << N << "x" << N << endl;
        cout << "El programa ha tardado " << (t1-t0) / CLOCKS_PER_SEC << " segundos" << endl;

    }
    

    return 0;

}