#include "compensation.h"
#include <cmath>
#include <complex>
#include <algorithm>
using namespace QtCharts;
#define N 4000
/*
fftwf_complex* in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
fftwf_comple* out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
*/
double cmp2 (double a,double c ,double f, double z, double f_d )//на вход частота в герцах ,f_d / f, для воввтановления передаем f_d
{
    return f_d/f;
}
double cmp (double a,double c ,double f, double z )//на вход частота в герцах ,f_d / f
{
    return c*1000*z/(pow(a,2)*M_PI*f);
}
double SQR(double a,double c ,double f, double z )
{

    return sqrt(1+pow(cmp(a,c,f,z),2));
    //return sqrt(1+pow(c*1000*z/(pow(a,2)*M_PI*f),2));
}
//компенсация с учетом размера:
complex Compensator::recovery_one_f(complex P0,double f0)//на вход частота в мегагерцах
{
    bool rigth_compensation =ui->radioButton->isChecked();

    if(rigth_compensation)
    {
        if (f0>F_N) {f0 = f0 - 2*F_N;}
        double B = pow(r0,2)/(a*a*(1+pow(cmp2(a,c,f0,L, f_d),2)));
        double ex1 = exp(-B);
        double module = M_PI * pow(a,2) * (pow(1.0 - ex1* cos(cmp2(a,c,f0,L, f_d)*B),2) + pow(ex1* sin(cmp2(a,c,f0,L, f_d)*B),2));
        complex q=complex((1.0 - ex1* cos(cmp2(a,c,f0,L, f_d)*B))/module,
                      -ex1* sin(cmp2(a,c,f0,L, f_d)*B)/module);
        return P0*q;
    }
    else
    {
        if (f0>F_N) {f0 = f0 - 2*F_N;}//когда этого нет а есть f-F_N в теле то работает

        complex q = complex(1.0  , -cmp2(a,c,f0,L, f_d));//получается что тут  1 - i*f_d/f
        return P0*q;

    }
}
//с учетом размера приемника:
complex Compensator::S(complex P0,double f0)//на вход частота в мегагерцах
{

    if (f0>F_N) {f0 = f0 - 2*F_N;}
    double B = pow(r0,2)/(a*a*(1+pow(cmp(a,c,f0,L),2)));
    double ex1 = exp(-B);
    complex q=complex(M_PI * pow(a,2) *  (1.0 - ex1* cos(cmp(a,c,f0,L)*B)), M_PI * pow(a,2) * ex1* sin(cmp(a,c,f0,L)*B));
    return P0*q;
}
//просто дифракция:
complex Compensator::S0(complex P0, double f0)//на вход частота в герцах
{
    if (f0>F_N) {f0 = f0 - 2*F_N;}//когда этого нет а есть f-F_N в теле то работает
    double b = 1 + pow(cmp(a,c,f0,L),2);
    double ksi = cmp(a,c,f0,L);//знак зависит от того че за фурье
    //double B = pow(0,2)/(a*a*b);// r0 тут стоял размерность!!
    //double ex1 = exp(-B);
    //complex q=complex(1.0 / b * (cos(ksi*B)+ksi*sin(ksi*B)), -1.0 / b * (ksi*cos(ksi*B) - sin(ksi*B)));
    complex q = complex(1.0 / b , 1.0 / b * ksi);//получается что тут  1 - i*f_d/f
    return P0*q;
}


Compensator::~Compensator()
{
    delete ui;
}
void Compensator::Init_signal_and_spectrum()
{

    for(int i=0; i<N; i++)
    {
        t[i]=i*dt;//в секундах
        f_[i]=i*df;//в герцах
        P[i]=exp(-pow(t[i]-T/2,2)/pow(tau_0/1000000000,2));//tau_0 в наносекундах
    }
    fftwf_plan plan;
    for (int i=0;i<N;i++)
    {
        in[i][0]=P[i];
        in[i][1]=0;
    }
    plan = fftwf_plan_dft_1d(N, in, out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    for(int i=0;i<N;i++)
    {
        P_sp[i].re=out[i][0];
        P_sp[i].im=out[i][1];
    }
    QLineSeries *series1 = new QLineSeries();

    for (int i =0; i<N;i++)
    {
        series1->append(f_[i]/1000000,P_sp[i].Abs());//выводим начальный спектр
    }

    QPen pen;
    pen.setColor(QColor(255,0,0));
    pen.setWidth(5);
    QChart *chart1 = new QChart();
    chart1->addSeries(series1);
    chart1->setTitle("Начальный спектр");
    chart1->createDefaultAxes();
    chart1->axisX()->setTitleText("f, MHz");
    chart1->axisY()->setTitleText("P(f)");
    QChartView * chartView1 = new QChartView(chart1,ui->widget_3);
    chartView1->setRenderHint(QPainter::Antialiasing);
    chartView1->resize(350,280);
    chartView1->show();
    QLineSeries *series2 = new QLineSeries();
    QPen pen2;
    pen2.setColor(QColor(0,255,0));
    series2->setPen(pen2);
    for (int i=0;i<N;i++)
    {
        series2->append(t[i]*1000000,P[i]);
    }
    QChart *chart2 = new QChart();
    chart2->addSeries(series2);
    chart2->setTitle(" Начальный импульс");
    chart2->createDefaultAxes();
    chart2->axisX()->setTitleText("t, микросекунды");
    chart2->axisY()->setTitleText("P(t)");
    series2->setName("Границы пучка");
    QChartView * chartView2 = new QChartView(chart2,ui->widget);
    chartView2->setRenderHint(QPainter::Antialiasing);
    chartView2->resize(350,280);
    chartView2->show();
}


void Compensator::Graph3()
{
    for (int i=1; i<N;i++)//домножаем спектр нулевую частоту пропустил

    {
        P2[i]=S0(P_sp[i],f_[i]);
        P3[i]=S(P_sp[i],f_[i]);
        in[i][0]=P2[i].re;
        in[i][1]=P2[i].im;
        in2[i][0]=P3[i].re;
        in2[i][1]=P3[i].im;
    }
    fftwf_plan plan1;
    plan1 = fftwf_plan_dft_1d(N, in, out,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwf_execute(plan1);
    fftwf_destroy_plan(plan1);
    fftwf_plan plan2;
    plan2 = fftwf_plan_dft_1d(N, in2, out2,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwf_execute(plan2);
    fftwf_destroy_plan(plan2);
    for(int i=0;i<N;i++)
    {
        P2[i].re=out[i][0]/(N);
        P2[i].im=out[i][1]/(N);
        P3[i].re=out2[i][0]/(N);
        P3[i].im=out2[i][1]/(N);
    }
    Graph4();
}

void Compensator::recovery_signal()
{
    //P3 - дифрагированный сигнал с учетом размера применика
    for (int i=0; i < N; i++)
    {
        in[i][0]=P3[i].re;
        in[i][1] =0;
    }
    fftwf_plan plan1;
    plan1 = fftwf_plan_dft_1d(N, in, out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftwf_execute(plan1);
    fftwf_destroy_plan(plan1);
    for (int i=0; i < N; i++)
    {
        P_sp[i].re = out[i][0];
        P_sp[i].im = out[i][1];
    }
    for (int i=1; i<N;i++)//домножаем спектр нулевую частоту пропустил
    {

        P_rec[i]=recovery_one_f(P_sp[i],f_[i]);
        in[i][0]=P_rec[i].re;
        in[i][1]=P_rec[i].im;
    }
    fftwf_plan plan2;
    plan2 = fftwf_plan_dft_1d(N, in, out,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwf_execute(plan2);
    fftwf_destroy_plan(plan2);
    for(int i=0;i<N;i++)
    {
        P_rec[i].re=out[i][0]/(N);
        P_rec[i].im=out[i][1]/(N);
    }

}
void Compensator::draw_recovery_signal()
{
    QLineSeries *series1 = new QLineSeries();
    QLineSeries *series2 = new QLineSeries();
    double max = P_rec[0].re;
    int max_n = 0;
    for (int i =0; i<N;i++)
    {
        series1->append(t[i]*1000000,P_rec[i].re);//!!!!!! rec  vs sp
        if (P_rec[i].re >= max)
        {max = P_rec[i].re; max_n =i;}
    }
    series2->append(t[max_n]*1000000, max);
    series2->append(t[max_n]*1000000, 0);
    ui->doubleSpinBox_5->setValue(t[max_n]*1000000);
    QPen pen;
    pen.setColor(QColor(255,0,0));
        pen.setWidth(1);
    series2->setPen(pen);
    QChart *chart1 = new QChart();
    chart1->addSeries(series1);
    chart1->addSeries(series2);
    chart1->setTitle("Востановленный сигнал ");
    chart1->createDefaultAxes();
    chart1->axisX()->setTitleText("t, микросекунд");
    chart1->axisY()->setTitleText("Р_rec(t)");
    chart1->axisX()->setRange(t1,t2);
    QChartView * chartView1 = new QChartView(chart1,ui->widget_6);
    chartView1->setRenderHint(QPainter::Antialiasing);
    chartView1->resize(350,280);
    chartView1->show();

}
void Compensator::Graph4()
{
    QLineSeries *series1 = new QLineSeries();
    QLineSeries *series2 = new QLineSeries();
    for (int i =0; i<N;i++)
    {
        series1->append(t[i]*1000000,P3[i].re);
        series2->append(t[i]*1000000,P2[i].re);
    }
    QPen pen;
    pen.setColor(QColor(255,0,0));
        pen.setWidth(5);

    QChart *chart1 = new QChart();
    QChart *chart2 = new QChart();
    chart1->addSeries(series1);
    chart1->setTitle("Сигнал на приенике ");
    chart1->createDefaultAxes();
    chart1->axisX()->setTitleText("t, микросекунд");
    chart1->axisY()->setTitleText("Р2(t)-P(t)");
    chart1->axisX()->setRange(t1,t2);
    QChartView * chartView1 = new QChartView(chart1,ui->widget_2);
    chartView1->setRenderHint(QPainter::Antialiasing);
    chartView1->resize(350,280);
    chartView1->show();

    /*chart2->addSeries(series2);
    chart2->setTitle("Сигнал на оси");
    chart2->createDefaultAxes();
    chart2->axisX()->setTitleText("t, микросекунд ");
    chart2->axisY()->setTitleText("P(t)");
    chart2->axisX()->setRange(t1,t2);
    QChartView * chartView2 = new QChartView(chart2,ui->widget_5);
    chartView2->setRenderHint(QPainter::Antialiasing);
    chartView2->resize(350,280);
    chartView2->show();*/
}

void Compensator::Graph()//рисует профиль распространения, частота в герцах , все в герцах сделал
{
    QLineSeries *series1 = new QLineSeries();
    QLineSeries *series2 = new QLineSeries();
    QLineSeries *series3 = new QLineSeries();
    QPen pen, pen2;
    pen.setColor(QColor(0,255,0));
    pen2.setColor(QColor(255,0,0));
    pen2.setWidth(5);
    series1->setPen(pen);
    series2->setPen(pen);
    series3->setPen(pen2);
    for (double z=0;z<L;z=z+0.1)
    {
        series1->append(z,a*SQR(a,c,f*1000000,z));
        series2->append(z,-a*SQR(a,c,f*1000000,z));
    }
    series3->append(L,r0);
    series3->append(L,-r0);
    QChart *chart1 = new QChart();
    chart1->addSeries(series1);
    chart1->addSeries(series2);
    chart1->addSeries(series3);
    chart1->setTitle(QString("Поперечное распределение"));
    chart1->createDefaultAxes();
    chart1->axisX()->setTitleText("Z, mm");
    chart1->axisY()->setTitleText("r, mm");
    series1->setName("Границы пучка");
    series3->setName("Приемник");
    QChartView * chartView1 = new QChartView(chart1,ui->widget_4);
    chartView1->setRenderHint(QPainter::Antialiasing);
    chartView1->resize(350,280);
    chartView1->show();
}


void Compensator::on_init_signal_and_spectrum_clicked()
{
    Graph();
    Init_signal_and_spectrum();
    Graph3();
    F_d_const = c*1000*L/(pow(a,2)*M_PI);
}





