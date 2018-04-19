#include <cstdlib>
#include <cstdio>
#include <vector>
#include <unistd.h>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <cmath>
// size of double - 8 bytes always
// long - mb 4 mb 8
//fseek(File, 9 bute,SEEK_CUR(SET,END))
size_t LONG_SIZE = 4;
size_t DOUBLE_SIZE = 8;
struct InitParameters
{
    int version, wigth, height, sample_lenght;
    double horz_step, vert_step, freq,zero_pos_1{1},zero_pos_2{1};
    InitParameters();
};
InitParameters::InitParameters()
{
    fread(&version, LONG_SIZE, 1, stdin);
    fread(&wigth, LONG_SIZE, 1, stdin);
    fread(&height, LONG_SIZE, 1, stdin);
    fread(&horz_step, DOUBLE_SIZE, 1, stdin);
    fread(&vert_step, DOUBLE_SIZE, 1, stdin);
    fread(&sample_lenght,LONG_SIZE, 1, stdin);
    fread(&freq, DOUBLE_SIZE, 1, stdin);
}
double contrast(std::vector<double>& payload, double max, int left, int right)
{
    double contrast = 0;
    for(auto it = payload.begin()+left; it <= payload.begin() + right; it++)
        contrast += (*it) * (*it);
    contrast /= (right - left +1) * (right - left );//n = r-l +1 (?)
    contrast = sqrt(contrast);
    contrast = max / contrast;
    return contrast;
}


double by_3_points(std::vector<double>& y, int n)
{
    return n + (y[n-1] - y[n+1])/(2 *(y[n-1] - 2  * y[n] + y[n+1]));
}

int main(int argc, char* argv[])
{

    std::vector<double> payload;
    
    FILE * file = fopen("matrix.txt", "w");
    InitParameters init_parameters;
     
    fprintf(stderr,"wigth %d\n", init_parameters.wigth);
    fprintf(stderr,"height %d\n",init_parameters.height);
    int current_wigth =1;
    int current_height =1;
    int bottom_time = atoi(argv[1]);
    fprintf(stderr,"botom tim %d\n", bottom_time);
    double coeff_1 = 0.8;
    double coeff_2 = 1.1;
//    double h = i
    int left = bottom_time * coeff_1;
    int right = bottom_time * coeff_2;
    double max =0;
    /*double** matrix;
    matrix = new double* [init_parameters.height];
    for(int i=0; i<init_parameters.height;i++)
            matrix[i] = new double[init_parameters.wigth];
    */
    for(current_height = 1; current_height <= init_parameters.height; current_height++)
    {       
        for(current_wigth = 1; current_wigth <= init_parameters.wigth; current_wigth++)
        {
            payload.resize(4096);
            fread(&init_parameters.zero_pos_1, DOUBLE_SIZE, 1, stdin);
            if(init_parameters.zero_pos_1)
                    fprintf(stderr, "ALARM %f \n", init_parameters.zero_pos_1);
            fread(&payload[0], DOUBLE_SIZE, 4096, stdin);
            std::vector<double>::iterator max_n = (std::min_element(payload.begin() + left ,payload.begin() + right));
            max = *max_n;
            int n = max_n - payload.begin();
            //fprintf(stderr,"%f \n",contrast(payload, max, left, right));
            //if(contrast(payload, max, left, right) < 10)
            //       fprintf(stderr,"contrast < 10");
            fprintf(file,"%f ", by_3_points(payload,n));
            
            //if (n<2970)
              //  fprintf(stderr,"h %d w %d max %d\n",current_height,current_wigth, n);

            //fprintf(stderr,"%f \n", by_3_points(payload,n));
        }
        fprintf(file,"\n");
    }
    fclose(file);
    return 0;
}




