#define PYOPENCL_DEFINE_CDOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64: enable

inline float my_dot (float* a, float*b, int length) {
   float c = 0;
   for (int i = 0; i < length; i++) {
      c += a[i] * b[i];
   }
   return c;
}

__kernel void sumint00 (
   const float qsize, const float ehc, const int pixels,
   __global const float4* points, const int npts, __global float* intensity)
{
// For no symmetry and no small angle approximation //
   float temp_intensity,temp_intensity_2,QdotR;
   int n = get_global_id(0);
   int i = n/pixels;
   int j = n%pixels;
   float Q[3] = { i*qsize/pixels-0.5*qsize, j*qsize/pixels-0.5*qsize,
        2*ehc*pow(sin(sqrt(pow((i-0.5*pixels),2)+pow(j-0.5*pixels,2))*qsize/(pixels*2*ehc)),2) };
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

__kernel void sumint10 (
   const float qsize, const float ehc, const int pixels,
   __global const float4* points, const int npts, __global float* intensity)
{
// For symmetry but no small angle approximation //
   float temp_intensity;
   int n = get_global_id(0);
   int i = n/pixels;
   int j = n%pixels;
   float Q[3] = { i*qsize/pixels-0.5*qsize, j*qsize/pixels-0.5*qsize,
        2*ehc*pow(sin(sqrt(pow((i-0.5*pixels),2)+pow(j-0.5*pixels,2))*qsize/(pixels*2*ehc)),2) };
   temp_intensity = 0;
   for ( int p = 0; p < npts; p++) {
       float R[3] = {points[p][0],points[p][1],points[p][2]};
       temp_intensity += points[p][3]*cos(my_dot(Q,R,3));
   }
   intensity[n] = pow(temp_intensity,2);
}

__kernel void sumint11 (
   const float qsize, const int pixels,
   __global const float4* points, const int npts, __global float* intensity)
{
// For symmetry and small angle approximation //
   float temp_intensity = 0;
   int n = get_global_id(0);
   int i = n/pixels;
   int j = n%pixels;
   float Q[2] = { i*qsize/pixels-0.5*qsize, j*qsize/pixels-0.5*qsize};
   for ( int p = 0; p < npts; p++) {
       float R[2] = {points[p][0],points[p][1]};
       temp_intensity += points[p][3]*cos(my_dot(Q,R,2));
   }
   intensity[n] = pow(temp_intensity,2);
}


