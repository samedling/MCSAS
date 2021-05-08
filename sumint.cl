// This file contains the OpenCL code for much faster calculations of the intensity.
// Separate kernels exist both with and without symmetry and the small angle approximation, as well as those without any points masked.
// The first number in the title indicates if it is symmetric (1=symmetric,0=asymmetric).
// The second number in the title indicates if the small-angle approximation is used.

#define PYOPENCL_DEFINE_CDOUBLE

//This line must be commented out to work on HD4000 graphics.
//#pragma OPENCL EXTENSION cl_khr_fp64: enable


inline float my_dot (float* a, float*b, int length) {
   float c = 0;
   for (int i = 0; i < length; i++) {
      c += a[i] * b[i];
   }
   return c;
}

inline float square (float x) {
    float y = x*x;
    return y;
}


__kernel void sumint_asym (
   const float qsize, const float ehc, const int x_pixels, const int y_pixels,
   __global const float4* points, const int npts, __global float* intensity)
{
// For no symmetry and no small angle approximation
   float temp_intensity,temp_intensity_2,QdotR;
   int n = get_global_id(0);
   int i = n/x_pixels;
   int j = n%x_pixels;
   float max_pixels = max(x_pixels,y_pixels);
   float Q[3] = { (i-0.5*y_pixels)*qsize/max_pixels, (j-0.5*x_pixels)*qsize/max_pixels,
        2*ehc*pow(sin(sqrt(square(i-0.5*y_pixels)+square(j-0.5*x_pixels))*qsize/(max_pixels*2*ehc)),2) };
   temp_intensity = 0;
   temp_intensity_2 = 0;
   for ( int p = 0; p < npts; p++) {
       float R[3] = {points[p][0],points[p][1],points[p][2]};
       QdotR = my_dot(Q,R,3);
       temp_intensity += points[p][3]*cos(QdotR);
       temp_intensity_2 += points[p][3]*sin(QdotR);
   }
   intensity[n] = pow(temp_intensity,2)+pow(temp_intensity_2,2);
}

__kernel void sumint_sym (
   const float qsize, const float ehc, const int x_pixels, const int y_pixels,
   __global const float4* points, const int npts, __global float* intensity)
{
// For symmetry but no small angle approximation
   float temp_intensity;
   int n = get_global_id(0);
   int i = n/x_pixels;
   int j = n%x_pixels;
   float max_pixels = max(x_pixels,y_pixels);
   float Q[3] = { (i-0.5*y_pixels)*qsize/max_pixels, (j-0.5*x_pixels)*qsize/max_pixels,
        2*ehc*pow(sin(sqrt(square(i-0.5*y_pixels)+square(j-0.5*x_pixels))*qsize/(max_pixels*2*ehc)),2) };
   temp_intensity = 0;
   for ( int p = 0; p < npts; p++) {
       float R[3] = {points[p][0],points[p][1],points[p][2]};
       temp_intensity += points[p][3]*cos(my_dot(Q,R,3));
   }
   intensity[n] = pow(temp_intensity,2);
}


//Pass in separate x and y arrays.//
__kernel void sumint_asym_mask (
   const float qsize, const float ehc, const int x_pixels, const int y_pixels, __global const int* xval, __global const int* yval,
   __global const float4* points, const int npts, __global float* intensity)
{
// For no symmetry and no small angle approximation
   float temp_intensity,temp_intensity_2,QdotR;
   int n = get_global_id(0);
   int i = xval[n];
   int j = yval[n];
   float max_pixels = max(x_pixels,y_pixels);
   float Q[3] = { (i-0.5*y_pixels)*qsize/max_pixels, (j-0.5*x_pixels)*qsize/max_pixels,
        2*ehc*pow(sin(sqrt(square(i-0.5*x_pixels)+square(j-0.5*y_pixels))*qsize/(max_pixels*2*ehc)),2) };
   temp_intensity = 0;
   temp_intensity_2 = 0;
   for ( int p = 0; p < npts; p++) {
       float R[3] = {points[p][0],points[p][1],points[p][2]};
       QdotR = my_dot(Q,R,3);
       temp_intensity += points[p][3]*cos(QdotR);
       temp_intensity_2 += points[p][3]*sin(QdotR);
   }
   intensity[n] = pow(temp_intensity,2)+pow(temp_intensity_2,2);
}


__kernel void sumint_sym_mask (
   const float qsize, const float ehc, const int x_pixels, const int y_pixels, __global const int* xval, __global const int* yval,
   __global const float4* points, const int npts, __global float* intensity)
{
// For symmetry but no small angle approximation
   float temp_intensity;
   int n = get_global_id(0);
   int i = xval[n];
   int j = yval[n];
   float max_pixels = max(x_pixels,y_pixels);
   float Q[3] = { (i-0.5*y_pixels)*qsize/max_pixels, (j-0.5*x_pixels)*qsize/max_pixels,
        2*ehc*pow(sin(sqrt(square(i-0.5*x_pixels)+square(j-0.5*y_pixels))*qsize/(max_pixels*2*ehc)),2) };
   temp_intensity = 0;
   for ( int p = 0; p < npts; p++) {
       float R[3] = {points[p][0],points[p][1],points[p][2]};
       temp_intensity += points[p][3]*cos(my_dot(Q,R,3));
   }
   intensity[n] = pow(temp_intensity,2);
}

__kernel void sumint_long_mask (
    const float qsize, const float ehc, const float coherence_length, const int x_pixels, const int y_pixels, __global const int* xval, __global const int* yval,
    __global const float4* points, const int npts, __global float* intensity)
{
	float temp_intensity;
	int n = get_global_id(0);
	int i = xval[n];
	int j = yval[n];
   float max_pixels = max(x_pixels,y_pixels);
   float Q[3] = { (i-0.5*y_pixels)*qsize/max_pixels, (j-0.5*x_pixels)*qsize/max_pixels,
		2*ehc*pow(sin(sqrt(square(i-0.5*x_pixels)+square(j-0.5*y_pixels))*qsize/(max_pixels*2*ehc)),2) };
	temp_intensity = 0;
	for ( int p1 = 0; p1 < npts; p1++) {
		for ( int p2 = 0; p2 < npts; p2++) {
			float R[3] = {points[p1][0]-points[p2][0],points[p1][1]-points[p2][1],points[p1][2]-points[p2][2]};
			if (pow(R[0],2)+pow(R[1],2)+pow(R[2],2) < pow(coherence_length,2)) {
				temp_intensity += points[p1][3]*points[p2][3]*cos(my_dot(Q,R,3));
			}
		}
	}
	intensity[n] = temp_intensity;
}
	 
__kernel void sumint_long (
   const float qsize, const float ehc, const float coherence_length, const int x_pixels, const int y_pixels,
   __global const float4* points, const int npts, __global float* intensity)
{
   float temp_intensity;
   int n = get_global_id(0);
   int i = n/x_pixels;
   int j = n%x_pixels;
   float max_pixels = max(x_pixels,y_pixels);
   float Q[3] = { (i-0.5*y_pixels)*qsize/max_pixels, (j-0.5*x_pixels)*qsize/max_pixels,
   		2*ehc*pow(sin(sqrt(square(i-0.5*x_pixels)+square(j-0.5*y_pixels))*qsize/(max_pixels*2*ehc)),2) };
	temp_intensity = 0;
	for ( int p1 = 0; p1 < npts; p1++) {
		for ( int p2 = 0; p2 < npts; p2++) {
   		float R[3] = {points[p1][0]-points[p2][0],points[p1][1]-points[p2][1],points[p1][2]-points[p2][2]};
			if (pow(R[0],2)+pow(R[1],2)+pow(R[2],2) < pow(coherence_length,2)) {
				temp_intensity += points[p1][3]*points[p2][3]*cos(my_dot(Q,R,3));
			}
   	}
   }
	intensity[n] = temp_intensity;
}

__kernel void sumint_sym_small (
   const float qsize, const int x_pixels, const int y_pixels,
   __global const float4* points, const int npts, __global float* intensity)
{
// For symmetry and small angle approximation
   float temp_intensity = 0;
   int n = get_global_id(0);
   int i = n/x_pixels;
   int j = n%x_pixels;
   float max_pixels = max(x_pixels,y_pixels);
   float Q[2] = { (i-0.5*y_pixels)*qsize/max_pixels, (j-0.5*x_pixels)*qsize/max_pixels};
   for ( int p = 0; p < npts; p++) {
       float R[2] = {points[p][0],points[p][1]};
       temp_intensity += points[p][3]*cos(my_dot(Q,R,2));
   }
   intensity[n] = pow(temp_intensity,2);
}
	 


__kernel void sumint_sym_small_mask (
   const float qsize, const int x_pixels, const int y_pixels, __global const int* xval, __global const int* yval,
   __global const float4* points, const int npts, __global float* intensity)
{
// For symmetry and small angle approximation
   float temp_intensity = 0;
   int n = get_global_id(0);
   int i = xval[n];
   int j = yval[n];
   float max_pixels = max(x_pixels,y_pixels);
   float Q[2] = { (i-0.5*y_pixels)*qsize/max_pixels, (j-0.5*x_pixels)*qsize/max_pixels};
   for ( int p = 0; p < npts; p++) {
       float R[2] = {points[p][0],points[p][1]};
       temp_intensity += points[p][3]*cos(my_dot(Q,R,2));
   }
   intensity[n] = pow(temp_intensity,2);
}




__kernel void sumint_asym_small (
   const float qsize, const float ehc, const int x_pixels, const int y_pixels,
   __global const float4* points, const int npts, __global float* intensity)
{
// For no symmetry and small angle approximation
   float temp_intensity,temp_intensity_2,QdotR,QPdotR;
   int n = get_global_id(0);
   int i = n/x_pixels;
   int j = n%x_pixels;
   float max_pixels = max(x_pixels,y_pixels);
   float Q[2] = { (i-0.5*y_pixels)*qsize/max_pixels, (j-0.5*x_pixels)*qsize/max_pixels };
   temp_intensity = 0;
   temp_intensity_2 = 0;
   for ( int p = 0; p < npts; p++) {
       float R[2] = {points[p][0],points[p][1]};
       QdotR = my_dot(Q,R,2);
       temp_intensity += points[p][3]*cos(QdotR);
       temp_intensity_2 += points[p][3]*sin(QdotR);
   }
   intensity[n] = pow(temp_intensity,2)+pow(temp_intensity_2,2);
}

__kernel void sumint_asym_small_mask (
   const float qsize, const float ehc, const int x_pixels, const int y_pixels, __global const int* xval, __global const int* yval,
   __global const float4* points, const int npts, __global float* intensity)
{
// For no symmetry and small angle approximation
   float temp_intensity,temp_intensity_2,QdotR,QPdotR;
   int n = get_global_id(0);
   int i = xval[n];
   int j = yval[n];
   float max_pixels = max(x_pixels,y_pixels);
   float Q[2] = { (i-0.5*y_pixels)*qsize/max_pixels, (j-0.5*x_pixels)*qsize/max_pixels };
   temp_intensity = 0;
   temp_intensity_2 = 0;
   for ( int p = 0; p < npts; p++) {
       float R[2] = {points[p][0],points[p][1]};
       QdotR = my_dot(Q,R,2);
       temp_intensity += points[p][3]*cos(QdotR);
       temp_intensity_2 += points[p][3]*sin(QdotR);
   }
   intensity[n] = pow(temp_intensity,2)+pow(temp_intensity_2,2);
}
