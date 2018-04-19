#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include<QThread>
#include <fftw3.h>
#define Num 4000
namespace Ui {
class MainWindow;
}
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
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void Graph();

    void Init_signal_and_spectrum();

    void Graph3();

    void Graph4();

    void draw_recovery_signal();

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
    double P[Num];//начальный сигнал
    double t[Num];//ось времени
    double f_[Num];//ось частот
    double dt=1.0/100000000;//шаг времени
    double F_N=1.0/dt/2.0;//частота Найквиста
    double T=dt*Num;//длина импульса
    double df=1.0/T;//шаг частот
    double tau_0=100;// коэффициент : ширина импульса: T/tau_0
    double t1=19.5,t2=20.5;
    complex P_sp[Num];//спектр  сигнала
    complex P2[Num];// обработанный сигнал // тут и спектр и сигнал по ситуации
    complex P3[Num];//
    complex P_rec[Num]; // восстановленный сигнал
private slots:

    void on_horizontalSlider_L_sliderMoved(int position);

    void on_horizontalSlider_c_sliderMoved(int position);

    void on_horizontalSlider_f_sliderMoved(int position);

    void on_horizontalSlider_a_sliderMoved(int position);

    void on_horizontalSlider_d_sliderMoved(int position);

    void on_doubleSpinBox_L_valueChanged(double arg1);

    void on_doubleSpinBox_c_valueChanged(double arg1);

    void on_doubleSpinBox_f_valueChanged(double arg1);

    void on_doubleSpinBox_a_valueChanged(double arg1);

    void on_doubleSpinBox_d_valueChanged(double arg1);

    void on_horizontalSlider_valueChanged(int value);

    void on_horizontalSlider_tau_0_valueChanged(int value);

    void on_init_signal_and_spectrum_clicked();

    void on_transformation_clicked();



    void on_horizontalSlider_2_valueChanged(int value);

    void on_doubleSpinBox_tau_0_valueChanged(double arg1);

    void on_pushButton_clicked();



    void on_horizontalSlider_3_sliderMoved(int position);

    void on_doubleSpinBox_3_valueChanged(double arg1);

    void on_radioButton_clicked();


private:
    Ui::MainWindow *ui;


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

#endif // MAINWINDOW_H
