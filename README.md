# DSP HW3
<center> <font size=2> Jason < 2020/11 > </font> </center>

original hackmd : https://hackmd.io/@jasonho610/S1ISowT_D

## Problem
In this homework, there're two things to be accomplished.

1. Design several ideal low pass filter, and probe into each in terms of impulse response $h[n]$ , magnitude $\left\lvert H_{lp}(e^{j\omega}) \right \rvert$ and phase $\measuredangle H(e^{j\omega })$ of freqency reponse.
2. Implement the filter with input $x[n]$ , and yet produced output $y[n]$. 

But instead of using wavefile produced in HW1, I decided to substitude inputs $x[n]$ as:

|          | Pitch |Frequency   |
| -------- | ----  | --------    |
| $x_1[n]$ | A7    | 3520Hz      |
| $x_2[n]$ | A#7   | $\approx$ 3729.310Hz  |
| $x_3[n]$ | B7    | $\approx$ 3951.066Hz  |
| $x_4[n]$ | C8    | $\approx$ 4186.009Hz  |
| $x_5[n]$ | C#8   | $\approx$ 4434.922Hz  |
| $x_6[n]$ | D8    | $\approx$ 4698.636Hz  |
| $x_7[n]$ | D#8   | $\approx$ 4978.032Hz  |
| $x_8[n]$ | E8    | $\approx$ 5274.041Hz  |

<font size=2> Here's a little reminder, my applications are contructed through mac-Xcode and anaconda-jupyter_notebook which is safari-based. Therefore,
- Couldn't write *mkdir* on Xcode, haven't solved it yet.
- Couldn't download .py on safari, it kept turn into .html file automatically.
- Some uncompitable errors might occur if changing to other envoirment from time to time.

Sorry about these problems. </font>

## Code Demonstration
<font size=2>Divided into two parts, c++ for all the computation, and python for graphing to ensure the results are correct.</font>
### Cpp:

#### 1. Libraries, Wave Structure and some Side Functions<font size=2> (mostly same as HW1)</font>
```cpp=
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
```
```cpp=
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
```
```cpp=
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
```
#### 2. Variables
Why $w_c$ isn't $2\pi f_c$ but $2\pi f_cT_s$?
>Normalized frequency is a unit of measurement of frequency equivalent to cycles/sample. In digital signal processing (DSP), the continuous time variable, $t$, with units of seconds, is replaced by the discrete integer variable, $n$, with units of samples. More precisely, the time variable, in seconds, has been normalized (divided) by the sampling interval, $T_s$ (seconds/sample), which causes time to have convenient integer values at the moments of sampling. This practice is analogous to the concept of natural units, meaning that the natural unit of time in a DSP system is samples. The normalized value of a frequency variable, $f$, is $f/f_s$, where $f_s = 1/T_s$ is the sampling rate in samples/sec.

> In the "analog case" (more formally, the "continuous-time case") it's $\Omega = 2\pi f$, and in the "digital case" (more formally the "discrete-time case") it's $\omega = 2\pi f/f_s = 2\pi f_cT_s$

```cpp=
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
```
#### 3. Writing input $x[n]$ in .wav and .csv
```cpp=
for(int i=0;i<(int)sizeof(f)/sizeof(double);i++)
{
    vector<short int> temp;
    for(int j=0;j<NoS;j++)
        temp.push_back(0.1 * cos(2 * pi * f[i] * t[j]) * step(t[j]) * pow(2,15));
    x.push_back(temp)
...
```
#### 4. Impulse Response
The IR despicted in assignment goes with

$$h[n]=\begin{cases} \begin{eqnarray}&\frac{sin(\omega_c n)}{n\pi}*\delta(n-M)&, 0\leq n\leq 2M \\ \\ &0&, otherwise\end{eqnarray}\end{cases}$$

Our goal is to simplify $h[n]$ to be more visual.

First, We can rewrite $\frac{sin(\omega_c n)}{n\pi}$ as

$$\frac{sin(\omega_c n)}{n\pi}=\frac{sin(\omega_c n)}{\omega_c n}\cdot \frac{\omega_c}{\pi}=sinc(\omega_c n)\cdot \frac{\omega_c}{\pi}$$

to cope with the limitation problem, as we already knew $\displaystyle\lim_{x\to 0}sinc(x)=1$. Second, continue on our simplification.

$$
\begin{eqnarray}h[n]&=&
\begin{cases} 
&\left(sinc(\omega_c n)\cdot \frac{\omega_c}{\pi}\right)*\delta(n-M)&, 0\leq n\leq 2M
\\ \\
&0&, otherwise
\end{cases}
\\ \\
&=&\begin{cases} 
&\displaystyle\frac{\omega_c}{\pi}\sum_{k=-\infty}^{\infty}sinc(\omega_c k)\cdot \delta(n-k-M)&, 0\leq n\leq 2M
\\ \\
&0&, otherwise
\end{cases}
\\ \\
&=&\begin{cases} 
&\displaystyle\frac{\omega_c}{\pi}\sum_{k=-\infty}^{\infty}sinc(\omega_c k)\cdot \delta(-k-M)&, n=0
\\ \\
&\displaystyle\frac{\omega_c}{\pi}\sum_{k=-\infty}^{\infty}sinc(\omega_c k)\cdot \delta(1-k-M)&, n=1
\\ \\
&...&
\\ \\
&\displaystyle\frac{\omega_c}{\pi}\sum_{k=-\infty}^{\infty}sinc(\omega_c k)\cdot \delta(M-k)&, n=2M
\\ \\
&0&, otherwise
\end{cases}
\\ \\
&=&\begin{cases} 
&\frac{\omega_c}{\pi}sinc(\omega_c (-M))&, n=0
\\ \\
&\frac{\omega_c}{\pi}sinc(\omega_c (1-M))&, n=1
\\ \\
&...&
\\ \\
&\frac{\omega_c}{\pi}\cdot 1&, n=M
\\ \\
&...&
\\ \\
&\frac{\omega_c}{\pi}sinc(\omega_c (M))&, n=2M
\\ \\
&0&, otherwise
\end{cases}
\\ \\
&=&\begin{cases} 
&\frac{\omega_c}{\pi}sinc(\omega_c (n-M))&, 0\leq n\leq 2M
\\ \\
&0&, otherwise
\end{cases}
\end{eqnarray}
$$

We end up with the simple conclution above, yet can be implemented easily, and also we find out $h[n]$ is a even function.

```cpp=
double impulse_response(int n, int M, double w_c, double pi)
{
    double centr = (double)(n - M);
    return sinc(centr*w_c)*(w_c/pi);
}
```
```cpp=
for(int j=0;j<2*M[i]+1;j++)
{
    h_n.push_back(impulse_response(j, M[i], w_c, pi));
    // uncentralized -> index: 0 ~ 2M (Ts)
    // centralized -> index: -M ~ M (Ts)
}
cout << "Calculating M = " << M[i] << ", h_n[" << h_n.size() << "] done" << endl;
```
#### 5. DTFT

The definition of DTFT:
$$H(\omega)=\displaystyle \sum_{n=-\infty}^{\infty}h[n]e^{−(j\omega n)}$$

```cpp=
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
```
The units of $\omega$  is set to 10,000 points.
```cpp=
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
```
<font size=2> Be cautious of the domain of phase. </font>

> <font size=2> atan returns the arctangent of x in the range -π/2 to π/2 radians. atan2 returns the arctangent of y/x in the range -π to π radians. If x is 0, atan returns 0. If both parameters of atan2 are 0, the function returns 0. All results are in radians. </font>

#### 6. Convolution
The definition of Convolution:

$$\begin{eqnarray}
(f*g)[n]&=&\sum_{m=-\infty}^{\infty}f[m]g[n-m]
\\
&=&\sum_{m=-\infty}^{\infty}f[n-m]g[m]
\end{eqnarray}$$

```cpp=
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
```
```cpp=
for(int j=0;j<x.size();j++)        // [freq]
{
    vector<short int> y = conv(x[j], h_n, "same");
    cout << " done" << endl;
...
```
<font size=2> This function is written to pattern on np.convolve, especially refercing on the *opt* part. 
[https://numpy.org/doc/stable/reference/generated/numpy.convolve.html#numpy.convolve](https://)</font>

### Python:
#### 1. Graphing $h[n]$, $\left\lvert H_{lp}(\omega) \right \rvert$, $\measuredangle H(\omega)$

```python=
def plot_ir(h_n, M, filename):
    n = np.arange(-M, M+1)
    plt.figure(figsize = (30, 10))
    #plt.xticks(n)
    plt.stem(n, h_n)                       # plot h[n]
    plt.plot(n, h_n, 'r-')     # plot h[n]'s envolop
    plt.title("FIR Impulse Response h[n], M= {:d}, ({:d} <= n <= {:d})".format(M, -M, M))
    plt.savefig(filename + ".svg", format="svg", facecolor='white')
    print("Generate SVG File : " + filename)
    if run_in_nb:
        plt.show()
    else:
        plt.close()  
```

```python=
def plot_mag(w, w_c, H_mag, M, dB, filename):
    plt.figure(figsize = (30, 10))
    unit = ""
    if dB == True:
        H_mag = 20*np.log10(H_mag)
        unit = "(dB)"
    plt.plot(w, H_mag)            # plot |H(ejw)|  Magnitude Response 
    plt.plot(-w, H_mag, color='b')            # plot |H(ejw)|  Magnitude Response 
    plt.title("Magnitude, {:d} Lobes (M = {:d})".format(2*M+1, M))
    plt.xlabel("$\omega\ /\ sample$ (-$\pi$ to $\pi$)")
    plt.ylabel("$|H(e^{j\omega })|$ "+unit)
    ymin = plt.ylim()[0]
    xmin = plt.xlim()[0]
    xmax = plt.xlim()[1]
    plt.hlines(y = 0, color = 'y', xmin=-w_c, xmax=w_c)
    plt.hlines(y = ymin, color = 'y', xmin=w_c, xmax=xmax)
    plt.hlines(y = ymin, color = 'y', xmin=xmin, xmax=-w_c)
    plt.vlines(x = w_c, color = 'y', ymin=ymin, ymax=0)
    plt.vlines(x = -w_c, color = 'y', ymin=ymin, ymax=0)
    plt.annotate(s = "$w_c$", xy = (w_c, ymin), xytext = (w_c-0.05, ymin-1.5), xycoords = 'data', color = "b")
    plt.annotate(s = "$-w_c$", xy = (-w_c, ymin), xytext = (-w_c-0.05, ymin-1.5), xycoords = 'data', color = "b")
    plt.annotate(s = "Ideal FIR", xy = (np.pi,ymin), xytext = (np.pi, ymin+1), xycoords = 'data', color = "b")
    plt.grid(True)
    print("Generate SVG File :" + filename)
    plt.savefig(filename + ".svg", format = "svg", facecolor = 'white')
    if run_in_nb:
        plt.show()
    else:
        plt.close()  
```

```python=
def plot_phase(w, w_c, phase, filename):
    plt.figure(figsize = (30, 10))
    plt.plot(w, phase)           # plot H(ejw) - Phase Response y unit: degree
    plt.title("Phase, M = {:d}".format(M))
    plt.xlabel("$\omega\ /\ sample$ (-$\pi$ to $\pi$)")
    plt.ylabel("$ \measuredangle H(e^{j\omega}) (\degree) $")
    plt.axvline(x = w_c, color = 'y')
    plt.axvline(x = -w_c, color = 'y')
    plt.axhline(y = -180 , color = 'y')
    plt.axhline(y = 180 , color = 'y')
    ymax = 180
    ymin = -180
    plt.annotate(s = "$w_c$", xy = (w_c, ymin), xytext = (w_c-0.05, ymin-6), xycoords = 'data', color = "b")
    plt.annotate(s = "$-w_c$", xy = (-w_c, ymin), xytext = (-w_c-0.05, ymin-6), xycoords = 'data', color = "b")
    plt.annotate(s = "$180 \degree $", xy = (np.pi, ymax), xytext = (np.pi, ymax), xycoords = 'data', color = "b")
    plt.annotate(s = "$-180 \degree $", xy = (np.pi, ymin), xytext = (np.pi, ymin), xycoords = 'data', color = "b")
    plt.grid(True)
    print("Generate SVG File : " + filename)
    plt.savefig(filename + ".svg", format = "svg", facecolor = 'white')
    if run_in_nb:
        plt.show()
    else:
        plt.close()
```

```python=
for M in M_list:
    h_n = pd.read_csv("ImpR_M{:d}.csv".format(M))
    h_n = h_n.values.tolist()
    
    mag = pd.read_csv("FR_Mag_M{:d}.csv".format(M))
    mag = mag.values.tolist()
    
    phase = pd.read_csv("FR_Phase_M{:d}.csv".format(M))
    phase = phase.values.tolist()
    
    plot_ir(h_n, M, "ImpR_M{:d}".format(M))
    plot_mag(unit_omega, w_c, mag, M, True, "FR_Mag_M{:d}".format(M))
    plot_phase(unit_omega, w_c, phase, "FR_Phase_M{:d}".format(M))
```
#### 2. Graphing $x[n]$ and $y[n]$

<font size=2> Only plotting the signals when $0\leq n\leq 100$ is enough to observe the transient/steady state. </font>
```python=
def plot_IO_signal(sig_in, sig_out, N, filename):
    n = np.arange(0, N)
    #print(n.shape)
    plt.figure(figsize = (30, 10))
    plt.plot(n, sig_in, 'deeppink', label = "x[n]")
    plt.plot(n, sig_out, 'deepskyblue', label = "y[n]")
    plt.title("Input x[n] and Output y[n], n shown in [0, {:d}]".format(N))
    plt.legend(fontsize = 20)
    plt.savefig(filename + ".svg", format = "svg", facecolor = 'white')
    print("Generate SVG File : " + filename)
    if run_in_nb:
        plt.show()
    else:
        plt.close()  
```
```python=
for pitch in pitch_list:
    x = pd.read_csv("wave_{}Hz.csv".format(pitch))
    x = x.values.tolist()
    bound = 100
    
    for M in M_list:
        cut_x = x[0:bound]
        #print(x)
        y = pd.read_csv("filter_M{:d}_wave_{}Hz.csv".format(M, pitch))
        y = y.values.tolist()
        cut_y = y[0:bound]
        #print(y)

        plot_IO_signal(cut_x, cut_y, bound, "Comparison_filter_M{:d}_wave_{}Hz".format(M, pitch))
```
## Results
All the outputs include:
- 8 ImpR .svg and .csv
- 8 input $x[n]$ .csv and .wav
- 8 FR magnitude .svg and .csv
- 8 FR phase .svg and .csv
- 64 output $y[n]$ .csv and .wav
- 64 comparison .svg between $x[n]$ and $y[n]$

where .csv are just for plotting.

Image Examples:
<font size=2>---- ImpR_M16.svg </font>
![](https://i.imgur.com/fn5W91a.png)

<font size=2>---- FR_Mag_M16.svg </font>
![](https://i.imgur.com/pzjs72J.png)

<font size=2>---- FR_Phase_M64.svg </font>
![](https://i.imgur.com/pywAnX1.png)

<font size=2>---- Comparison_filter_M2048_wave_A#7Hz.svg </font>
![](https://i.imgur.com/7q0BmKc.png)

<font size=2>---- Comparison_filter_M16_wave_C8Hz.svg </font>
![](https://i.imgur.com/CTgw7v0.png)

<font size=2>---- Comparison_filter_M256_wave_D8Hz.svg </font>
![](https://i.imgur.com/ORSQGIy.png)

## Discussion
We can observe several characteristics,

1. Ideal lowpass may be characterized by a gain of 1 for all frequencies below some cut-off frequency $f_c$ in Hz, and a gain of 0 for all higher frequencies. Our cutoff frequency is 4000Hz, meaning that pitches higher than B7 should be eliminated. (C8~E8 are all muted.) Although in some cases, lp filter soothes the signal but doen't eliminate it, when M = 1, 4 specifically. (lower the loudness)

> The low-pass filter is sometimes called a high-cut filter, or treble-cut filter in audio applications.

2. The longer the IR is (the larger M is), the filter gives the better peformance on cutoff; however, the later the output gets out the transient/gets in to the steady state.

3. All phases in $[-\omega_c,\omega_c]$ should be linear.

These phenomena can be proved either on theory, mathmetics, and the images.

## Proof

How we derive $h_{lp}$ from $H_{lp}$, and correspond to the form depicted in assignment?

$$\begin{eqnarray}
h_{lp}[n]&=&DTFT^{-1}(H_{lp}(\omega))
\\
&=&\frac{1}{2\pi}\int_{-\pi}^{\pi}d\omega e^{j\omega n}
\begin{cases}1, \quad\lvert\omega\rvert\leq\omega_c
\\
0,\quad otherwise
\end{cases}
\\
&=&\frac{1}{2\pi}\int_{-\omega_c}^{\omega_c}e^{j\omega n}d\omega
\\
&=&\frac{1}{2\pi jn}[e^{j\omega_c n}-e^{-j\omega_c n}]
\\
&=&\frac{sin(\omega_cn)}{\pi n}\quad ,n\in \mathbb{Z}
\end{eqnarray}$$


