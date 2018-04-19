#ifndef COMPENSATION_HPP
#define COMPENSATION_HPP
#include <cmath>
#include <fftw3.h>
class complex
{
public:
    double re;
    double im;
    complex() { re = 0.0; im = 0.0; }
    complex(double x, double y) { re = x; im = y; }
    ~complex(){}
    double Abs() {return (sqrt(pow(re,2)+pow(im,2)));}
    //complex & operator =(complex &c){return (complex(re+c.re, im+c.im));}
    complex operator +(complex c){return (complex(re+c.re, im+c.im));}
    complex operator +(double x){return (complex(1.0*x+re, im));}//complex + int тогда все работает
    complex operator -(complex c){return (complex(re-c.re, im-c.im));}
    complex operator *(complex c){return (complex(re * c.re - im * c.im, re * c.im + im * c.re));}
    complex operator *(double x){return (complex(re*x,im*x));}
    complex inv(){return complex(re/pow(Abs(),2),-im/pow(Abs(),2));}

};

class Compensator
{
public:
    
    explicit Compensator(int Num)
    {
        P = new double[Num];
        t = new double[Num];
        f_ = new double[Num];
        P_sp = new complex[Num];
        P2 = new complex[Num];
        P3 = new complex[Num];
        P_rec = new complex[Num];
        T = dt * Num; 
        in  = new fftwf_complex[Num];
        in2  = new fftwf_complex[Num];
        out  = new fftwf_complex[Num];
        out2  = new fftwf_complex[Num];
    };
    ~Compensator();


    void recovery_signal();

    complex recovery_one_f(complex P0,double f);

    complex S(complex P0,double f);

    complex S0(complex P0,double f);

    double L=50;//длина объекта, мм
    double a=4 ;// начальный радиус пучка, мм
    double r=0; // поперечное расстояние, мм
    double c=5860; // м/с
    double f=5; // Мгц
    double r0=a; // радиус приемника, мм
    double Norma=0;
    double NormD=0;
    double s=1;
    double f_d= c*1000*L/(pow(a,2)*M_PI);
    double F_d_const= c*1000*L/(pow(a,2)*M_PI);
    double* P=nullptr;//начальный сигнал
    double* t=nullptr;//ось времени
    double*  f_=nullptr;//ось частот
    double dt=1.0/100000000;//шаг времени
    double F_N=1.0/dt/2.0;//частота Найквиста
    double T{0};//длина импульса
    //double T=dt*Num;//длина импульса
    double df=1.0/T;//шаг частот
    double tau_0=100;// коэффициент : ширина импульса: T/tau_0
    double t1=19.5,t2=20.5;
    complex*  P_sp=nullptr;//спектр  сигнала
    complex*  P2=nullptr;// обработанный сигнал // тут и спектр и сигнал по ситуации
    complex*  P3=nullptr;//
    complex*  P_rec=nullptr; // восстановленный сигнал
    fftwf_complex* in;
    fftwf_complex* in2;
    fftwf_complex* out;
    fftwf_complex* out2;


};

class Sleeper: public QThread
{
    public:
        static void msleep(int ms)
        {
            QThread::msleep(ms);
        }
}
;

#endif
