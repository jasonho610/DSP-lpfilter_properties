//
//  main.cpp
//  DSP HW3
//
//  Created by 何冠勳 on 2020/11/4.
//  Copyright © 2020 何冠勳. All rights reserved.
//
#define f_s 16000       // sample freq = sample rate (Hz)
#define MP 2           // music period
#define NoS f_s*MP      // number of samples

#include <iostream>
#include <complex>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
using namespace std;

typedef struct WaveHeader
{
    // riff wave header
    char ChunkID[4] = {'R','I','F','F'};
    unsigned int ChunkSize;        // 0 ~ FFFF,FFFF
    char Format[4] = {'W','A','V','E'};

    // fmt subchunk
    char SubChunk1ID[4] = {'f','m','t',' '};
    unsigned int SubChunk1Size;    // 0 ~ FFFF,FFFF
    unsigned short AudioFormat;    // 0 ~ FFFF
    unsigned short NumChannels;    // 0 ~ FFFF
    unsigned int SampleRate;       // 0 ~ FFFF,FFFF
    unsigned int ByteRate;         // 0 ~ FFFF,FFFF
    unsigned short BlockAlign;     // 0 ~ FFFF
    unsigned short BitsPerSample;  // 0 ~ FFFF

    // data subchunk
    char SubChunk2ID[4] = {'d','a','t','a'};
    unsigned int SubChunk2Size;    // 0 ~ FFFF,FFFF
} WaveHeader;

WaveHeader make_WaveHeader(unsigned short const NumChannels, unsigned int const SampleRate, unsigned short const BitsPerSample)
{
    WaveHeader myWH;

    myWH.AudioFormat = 1;                  // 1 for PCM...
    myWH.SampleRate = SampleRate;
    myWH.NumChannels = NumChannels;        // 1 for Mono, 2 for Stereo
    myWH.BitsPerSample = BitsPerSample;
    myWH.ByteRate = (myWH.SampleRate * myWH.NumChannels * myWH.BitsPerSample)/8;
    myWH.BlockAlign = (myWH.NumChannels * myWH.BitsPerSample)/8;
    myWH.SubChunk2Size = (NoS * myWH.NumChannels * myWH.BitsPerSample)/8;
    myWH.SubChunk1Size = 16;               // 16 for PCM
    myWH.ChunkSize = 4+(8+myWH.SubChunk1Size)+(8+myWH.SubChunk2Size);

    return myWH;
}

bool step(float time)
{
    if(time<0.5) return 1;
    else return 0;
}

double sinc(double x)
{
    if(x == 0) return 1;
    else return sin(x)/x;
}

double impulse_response(int n, int M, double w_c, double pi)
{
    double centr = (double)(n - M);
    return sinc(centr*w_c)*(w_c/pi);
}

vector<double> linspace(double start, double end, int num)
{
    vector<double> linspaced;
    if(num==1)
    {
        linspaced.push_back(start);
        return linspaced;
    }
    else
    {
        double delta = (end-start)/(num-1);
        for(int i=0;i<num-1;i++)
            linspaced.push_back(start+delta*i);
        linspaced.push_back(end);
        return linspaced;
    }
}

complex<double> DTFT(double M, double w, vector<double> h_n)     // H = sig_k (h[k])*exp(-j*w*k))
{
    complex<double> H_w(0,0);
    double k = h_n.size();
    for(int i=0;i<k;i++)
    {
        complex<double> c(0, -(w*i));
        c = exp(c);
        H_w = H_w + c*h_n[i];
    }
    return H_w;
}

vector<short int> conv(vector<short int> signal, vector<double> imp_resp, string opt)
{
    // this function is to pattern on python model numpy.convolve()
    int M = (int)signal.size();
    int L = (int)imp_resp.size();
    int N = M+L-1;
    
    cout << "Calculating convolution, please wait";
    vector<short int> con;
    for(int n=0;n<N;n++)
    {
        short int temp = 0;
        for(int k=0;k<L;k++)
        {
            double input = 0.;
            if((n-k>=0)&&(n-k<M))
                input = signal[n-k];
            temp = temp+(short int)(input*imp_resp[k]);
        }
        con.push_back(temp);
        if(n%10666==0) cout << ".";
        //cout << "Calculating Convolution - " << round((float)n/N*100.0) << "% \r";
    }

    if(opt=="same")
    {
        int max_len = max(M,L);
        vector<short int> con_cut(con.begin()+(N-max_len)/2,con.begin()+(N-max_len)/2+max_len);
        return con_cut;
    }
    else return con;
}

int main()
{
    /* section of signal */
    vector<vector<short int>> x;     // [freq][n]
    double t[NoS] = {0};
    double const pi = acos(-1);
    double f[8] = {3520, 3729.310, 3951.066, 4186.009, 4434.922, 4698.636, 4978.032, 5274.041};
    string pitch[8] = {"A7", "A#7", "B7", "C8", "C#8", "D8", "D#8", "E8"};
    WaveHeader myWH = make_WaveHeader(1,f_s,16);
    int f_c = 4000;       // Cutoff freq is 4000 Hz.
    double T_s = (double)1/f_s;
    double w_c = 2*pi*f_c*T_s;      // Return to time manner (normalized cutoff freq).
    int M[8] = {1, 4, 16, 64, 256, 512, 1024, 2048};

    for(int i=0;i<NoS;i++)
        t[i]=(1.0/f_s)*i;

    for(int i=0;i<(int)sizeof(f)/sizeof(double);i++)
    {
        vector<short int> temp;
        for(int j=0;j<NoS;j++)
            temp.push_back(0.1 * cos(2 * pi * f[i] * t[j]) * step(t[j]) * pow(2,15));
        
        string wav_filename = "wave_" + pitch[i] + "Hz.wav";
        ofstream wavefile;
        cout << wav_filename << endl;
        wavefile.open(wav_filename, ofstream::binary|ofstream::out);
        wavefile.write((char*)&myWH, sizeof(myWH));
        for(int j=0;j<temp.size();j++)
            wavefile.write((char*) &temp[j], sizeof(temp[j]));
        wavefile.close();
        
        string x_filename = "wave_" + pitch[i] + "Hz.csv";
        cout << x_filename << endl;
        ofstream x_file(x_filename, ios::out);
        x_file << "Output_Signal" << endl;
        for(int j=0;j<temp.size();j++)
            x_file << temp[j] << endl;
        x_file.close();
        
        x.push_back(temp);
    }
    /* end */
    
    for(int i=0;i<sizeof(M)/sizeof(int);i++)
    //for(int i=0;i<4;i++)
    {
        /* Responses */
        vector<double> h_n;
        vector<complex<double>> H_w;
        vector<double> H_w_mag;
        vector<double> H_w_phase;

        for(int j=0;j<2*M[i]+1;j++)
        {
            h_n.push_back(impulse_response(j, M[i], w_c, pi));
            // uncentralized -> index: 0 ~ 2M (Ts)
            // centralized -> index: -M ~ M (Ts)
        }
        cout << "Calculating M = " << M[i] << ", h_n[" << h_n.size() << "] done" << endl;

        vector<double> unit_omega = linspace(-pi, pi, 10001);
        for(int j=0;j<unit_omega.size();j++)
        {
            complex<double> c = DTFT(M[i], unit_omega[j], h_n);
            H_w.push_back(c);
            H_w_mag.push_back(abs(c));
            H_w_phase.push_back(atan2(c.imag(), c.real())*180/pi);
            // uncentralized -> index: 0 ~ 1001 (unit_omega = 2pi/1000 = pi/500)
            // centralized -> index: (-pi), (-pi + pi/500), ~ (+pi)
        }
        cout << "Calculating DTFT - H_w[" << H_w.size() << "] done" << endl;
        /* end */

        /* Write FR info for ploting */
        string H_filename = "ImpR_M" + to_string(M[i]) + ".csv";
        ofstream H_file(H_filename, ios::out);
        if (!H_file)
        {
           cout << "Open file failed!" << endl;
           exit(1);
        }
        H_file << "Impulse_Response" << endl;
        for(int j=0;j<h_n.size();j++)
        {
            H_file << h_n[j] << endl;
        }
        H_file.close();

        string H_mag_filename = "FR_Mag_M" + to_string(M[i]) + ".csv";
        ofstream H_mag_file(H_mag_filename, ios::out);
        if (!H_mag_file)
        {
           cout << "Open file failed!" << endl;
           exit(1);
        }
        H_mag_file << "Magnitude" << endl;
        for(int j=0;j<H_w_mag.size();j++)
            H_mag_file << H_w_mag[j] << endl;
        H_mag_file.close();

        string H_phase_filename = "FR_Phase_M" + to_string(M[i]) + ".csv";
        ofstream H_phase_file(H_phase_filename, ios::out);
        if (!H_phase_file)
        {
           cout << "Open file failed!" << endl;
           exit(1);
        }
        H_phase_file << "Phase(deg)" << endl;
        for(int j=0;j<H_w_phase.size();j++)
            H_phase_file << H_w_phase[j] << endl;
        H_phase_file.close();
        /* end */

        /* Convolution */
        for(int j=0;j<x.size();j++)        // [freq]
        {
            vector<short int> y = conv(x[j], h_n, "same");
            cout << " done" << endl;
            
            /* Write output wav file*/
            string wav_filename = "filter_M" + to_string(M[i]) + "_wave_" + pitch[j] + "Hz.wav";
            ofstream wavefile;
            
            cout << wav_filename << endl;
            wavefile.open(wav_filename, ofstream::binary|ofstream::out);
            wavefile.write((char*)&myWH, sizeof(myWH));
            for(int k=0;k<y.size();k++)
                wavefile.write((char*) &y[k], sizeof(y[k]));
            wavefile.close();
            
            string y_filename = "filter_M" + to_string(M[i]) + "_wave_" + pitch[j] + "Hz.csv";
            cout << y_filename << endl;
            ofstream y_file(y_filename, ios::out);
            y_file << "Output_Signal" << endl;
            for(int k=0;k<y.size();k++)
                y_file << y[k] << endl;
            y_file.close();
        }
        /* end */
        cout << endl << endl;
    }
    return 0;
}

/*
 // Example program for linspace
 #include <iostream>
 #include <string>
 #include <vector>

 using namespace std;

 vector<double> linspace(double start, double end, int num)
 {
     vector<double> linspaced;
     if(num==1)
     {
         linspaced.push_back(start);
         return linspaced;
     }
     else
     {
         double delta = (end-start)/(num-1);
         for(int i=0;i<num-1;i++)
             linspaced.push_back(start+delta*i);
         linspaced.push_back(end);
         return linspaced;
     }
 }

 int main()
 {
     vector<double> centr = linspace(-10, 10, 21);
     for(int i=0;i<centr.size();i++)
         cout<<centr[i]<<endl;
 }
 */
