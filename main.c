
#include <stdlib.h>
#include <math.h>
#include <ti/dsplib/dsplib.h>

// бібліотека з коефіцієнтами фільтру з МатЛабу
#include "fdacoefs.h"


/* Global definitions */
//#define _LITTLE_ENDIAN
/* Число вибірок для сигналу*/
#define N 128
/* Число вибірок для яких необхідно виконати ШПФ*/
#define N_FFT 256

/* вирівнювання даних має бути подвійне слово  8 */
#define AL 8

#pragma DATA_ALIGN(brev, AL);
unsigned char brev[64] = {
    0x0, 0x20, 0x10, 0x30, 0x8, 0x28, 0x18, 0x38,
    0x4, 0x24, 0x14, 0x34, 0xc, 0x2c, 0x1c, 0x3c,
    0x2, 0x22, 0x12, 0x32, 0xa, 0x2a, 0x1a, 0x3a,
    0x6, 0x26, 0x16, 0x36, 0xe, 0x2e, 0x1e, 0x3e,
    0x1, 0x21, 0x11, 0x31, 0x9, 0x29, 0x19, 0x39,
    0x5, 0x25, 0x15, 0x35, 0xd, 0x2d, 0x1d, 0x3d,
    0x3, 0x23, 0x13, 0x33, 0xb, 0x2b, 0x1b, 0x3b,
    0x7, 0x27, 0x17, 0x37, 0xf, 0x2f, 0x1f, 0x3f
};


/* масив відліків вхідного сигналу та вихідного */
#pragma DATA_ALIGN(in_sig,AL) // вирівнювання
float   in_sig [2*N]; // для заповнення сигналом
// 2*N - для REAL i IMAGIN складових  комплексного сигналу

/* масив відліків вхідного сигналу та вихідного для ШПФ */
#pragma DATA_ALIGN(in_sig_fft,AL) // вирівнювання
float   in_sig_fft [2*N_FFT]; // для заповнення сигналом
// після ШПФ цей масив зміниться

#pragma DATA_ALIGN(out_signal,AL) // вирівнювання
float   out_signal [2*N_FFT]; // для спектру сигналу

/* масив відліків імпульсної характеристики фільтра   для ШПФ */
#pragma DATA_ALIGN(in_inpulse_fft,AL) // вирівнювання
float   in_inpulse_fft [2*N_FFT]; // для заповнення сигналом
// після ШПФ цей масив зміниться

#pragma DATA_ALIGN(out_impulse,AL) // вирівнювання
float   out_impulse[2*N_FFT]; // для спектру сигналу

#pragma DATA_ALIGN(w_sp, AL)
float   w_sp [2*N_FFT]; // для таблиці cos | sin - Twiddle factor
#pragma DATA_ALIGN(w_sp_ifft, AL)
float   w_sp_ifft [2*N_FFT]; // для таблиці cos | sin - Twiddle factor

#pragma DATA_ALIGN(a4h, AL)
float a4h[N_FFT];


#pragma DATA_ALIGN(a4h, AL)
float out_sig_fin[N+BL-1];


// прототипи ф-ій
void generateInput();
void seperateRealImg (float *real, float * img, float * cmplx, int size_cmplx);
void tw_gen_fft (float *w, int n);
void tw_gen_ifft (float *w, int n);
void magnitude(float *in_real, float * in_img, float * out_arr, int size);
void magnitude_cplx(float *cmplx_in, float * out_arr, int size);
void vecmul_cplx(float *a_vec_in, float * b_vec_in,  float * vec_out, int size);

int main(void) {
	int i;
	/* Генерування сигналу */
	generateInput(in_sig, N);

	// копія масиву в in_sig_fft
	for (i = 0; i < N; i++) {
		in_sig_fft[2*i] = in_sig[2*i];
		in_sig_fft[2*i + 1]  = in_sig[2*i+1];
	}


	// copy impulse response filter to real place 0, 2, 4
	for(i=0; i< BL; ++i){
		in_inpulse_fft[i*2] = B[i];
	}

	/* Функція генерує масив (таблицю) косинусів та сінусів*/
	tw_gen_fft(w_sp, N_FFT);
	tw_gen_ifft(w_sp_ifft, N_FFT);
	// реалізація ШПФ
	DSPF_sp_fftSPxSP(N_FFT, in_sig_fft, w_sp, out_signal, brev, 4, 0, N_FFT);
	//DSPF_sp_fftSPxSP(N, in_sig_fft, w_sp, out_signal, brev, 4, 0, N);

	// реалізація ШПФ impusle resp
	DSPF_sp_fftSPxSP(N_FFT, in_inpulse_fft, w_sp, out_impulse, brev, 4, 0, N_FFT);

	// multiply complex vectors
	vecmul_cplx(out_signal, out_impulse, out_impulse, 2*N_FFT);

	// Зворотнє ШПФ
	DSPF_sp_ifftSPxSP(N_FFT, out_impulse, w_sp_ifft, out_signal, brev, 4, 0, N_FFT);

	// for test point
	//	magnitude_cplx(out_impulse, a4h, N_FFT);

	// forming final out signal array - first N samples, like signal size
	for (i=0; i < N; ++i){
		out_sig_fin[i] =out_signal[2*i] ;
	}

	i=0;// for breakpoint

	return 0;
}



/*
    Функція генерує сигнал певної частоти і записує відліки в масив
*/
void generateInput (float * in_signal, int iN ) {
    int   i;
    float FreqSignal1, sinWaveMag1 ;
    float FreqSignal2, sinWaveMag2 ;
    float real_s, img_s ; // для реальної та уявної складової

    float FreqSample;

    /* частота  дискретизації Hz*/
    FreqSample = 48000;

    /* Амплітуда сигналу */
    sinWaveMag1 = 5;
    sinWaveMag2 = 10;

    /* частота  сигналу Hz*/
    FreqSignal1 = 2000;
    FreqSignal2 = 5550;

    /* генерування сигналу дійсних та уявних складових сигналу */
	for (i = 0; i < iN; i++) {
		real_s  = (float) (
				sinWaveMag1 * cos(2*3.14*i * FreqSignal1 /FreqSample)+
				sinWaveMag2 * cos(2*3.14*i* FreqSignal2 /FreqSample)
		);

		img_s = (float)0.0;

		in_signal[2*i] = real_s;
		in_signal[2*i + 1]  = img_s;
	}

}

/*
    Функція генерує масив (таблицю) косинусів та сінусів  для ШПФ
*/

void tw_gen_fft (float *w, int n)
{
    int i, j, k;
    const double PI = 3.141592654;

    for (j = 1, k = 0; j <= n >> 2; j = j << 2)
    {
        for (i = 0; i < n >> 2; i += j)
        {
#ifdef _LITTLE_ENDIAN
            w[k]     = (float) sin (2 * PI * i / n);
            w[k + 1] = (float) cos (2 * PI * i / n);
            w[k + 2] = (float) sin (4 * PI * i / n);
            w[k + 3] = (float) cos (4 * PI * i / n);
            w[k + 4] = (float) sin (6 * PI * i / n);
            w[k + 5] = (float) cos (6 * PI * i / n);
#else
            w[k]     = (float)  cos (2 * PI * i / n);
            w[k + 1] = (float) -sin (2 * PI * i / n);
            w[k + 2] = (float)  cos (4 * PI * i / n);
            w[k + 3] = (float) -sin (4 * PI * i / n);
            w[k + 4] = (float)  cos (6 * PI * i / n);
            w[k + 5] = (float) -sin (6 * PI * i / n);
#endif
            k += 6;
        }
    }
}

/* Function for generating Specialized sequence of twiddle factors */

void tw_gen_ifft (float *w, int n)
{
    int i, j, k;
    const double PI = 3.141592654;

    for (j = 1, k = 0; j <= n >> 2; j = j << 2)
    {
        for (i = 0; i < n >> 2; i += j)
        {
#ifdef _LITTLE_ENDIAN
            w[k]     = (float) -sin (2 * PI * i / n);
            w[k + 1] = (float)  cos (2 * PI * i / n);
            w[k + 2] = (float) -sin (4 * PI * i / n);
            w[k + 3] = (float)  cos (4 * PI * i / n);
            w[k + 4] = (float) -sin (6 * PI * i / n);
            w[k + 5] = (float)  cos (6 * PI * i / n);
#else
            w[k]     = (float) cos (2 * PI * i / n);
            w[k + 1] = (float) sin (2 * PI * i / n);
            w[k + 2] = (float) cos (4 * PI * i / n);
            w[k + 3] = (float) sin (4 * PI * i / n);
            w[k + 4] = (float) cos (6 * PI * i / n);
            w[k + 5] = (float) sin (6 * PI * i / n);
#endif
            k += 6;
        }
    }
}

/*
    Функція seperateRealImg
    для розділення real and imaginary даних після ШПФ
    Вона потрібна для того, щоб відобразити графічно дані,
    використовуючи CCS graph
*/

void seperateRealImg (float *real, float * img, float * cmplx, int size_cmplx) {
    int i, j;

    for (i = 0, j = 0; j < size_cmplx; i+=2, j++) {
    	real[j] = cmplx[i];
    	img[j] = cmplx[i + 1];
    }
}

/*
 Функція  обчислення модуля комплексного числа
 */
void magnitude(float *in_real, float * in_img, float * out_arr, int size){
	int i;
	 for (i = 0; i < size; i++) {
		 out_arr[i] = (float)sqrt(in_real[i]* in_real[i] + in_img[i]*in_img[i]);
	    }
}
// vers 2

void magnitude_cplx(float *cmplx_in, float * out_arr, int size){
	int i;
	 for (i = 0; i < size; i++) {
		 out_arr[i] = (float)sqrt(cmplx_in[i*2]* cmplx_in[i*2] + cmplx_in[i*2+1]*cmplx_in[i*2+1]);
	    }
}


void vecmul_cplx(float *a_vec_in, float * b_vec_in,  float * vec_out, int size){
	int i;
	float a1, a2, b1, b2;

	for (i = 0; i < size; i+=2) {
		 a1 = a_vec_in[i]; // real a_vec
		 b1 = a_vec_in[i+1]; // img a_vec
		 a2 = b_vec_in[i]; // real b_vec
		 b2 = b_vec_in[i+1]; // img b_vec

		 vec_out[i]  = a1 * a2 - b1 * b2;
		 vec_out[i+1] =  a1 * b2 + a2 * b1;
	 }
}
