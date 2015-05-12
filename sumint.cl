__kernel void sumint00 (
   __global real qsize, __global real ehc, __global int pixels,
   __global float* points, __global int npts, __global float* intensity)
{
// For no symmetry and no small angle approximation //
   int n = get_global_id(0);
   i = n%pixels;
   j = n/pixels;
   Q = { i*qsize/pixels-0.5*qsize, j*qsize/pixels-0.5*qsize,
        2*ehc*pow(sin(sqrt(pow((i-0.5*pixels),2)+pow(j-0.5*pixels,2))*qsize/(pixels*2*ehc)),2) };
   temp_intensity = 0;
   temp_intensity2 = 0;
   for ( int p = 0; p < npts; p++) {
       // QdotR = inner_product(Q,{points[p][0],points[p][1],points[p][2]}); //
       QdotR = inner_product(Q,points[p][0..2]);
       temp_intensity += points[p][3]*cos(QdotR);
       temp_intensity_2 += points[p][3]*sin(QdotR);
   }
   intensity[n] = pow(temp_intensity,2)+pow(temp_intensity_2,2)
}

__kernel void sumint10 (
   __global real qsize, __global real ehc, __global int pixels,
   __global float* points, __global int npts, __global float* intensity)
{
// For symmetry but no small angle approximation //
   int n = get_global_id(0);
   i = n%pixels;
   j = n/pixels;
   Q = { i*qsize/pixels-0.5*qsize, j*qsize/pixels-0.5*qsize,
        2*ehc*pow(sin(sqrt(pow((i-0.5*pixels),2)+pow(j-0.5*pixels,2))*qsize/(pixels*2*ehc)),2) };
   temp_intensity = 0;
   for ( int p = 0; p < npts; p++) {
       temp_intensity += points[p][3]*cos(inner_product(Q,points[p][0..2]));
   }
   intensity[n] = pow(temp_intensity,2)
}

__kernel void sumint11 (
   __global real qsize, __global int pixels,
   __global float* points, __global int npts, __global float* intensity)
{
// For symmetry and small angle approximation //
   int n = get_global_id(0);
   i = n%pixels;
   j = n/pixels;
   Q = { i*qsize/pixels-0.5*qsize, j*qsize/pixels-0.5*qsize};
   temp_intensity = 0;
   temp_intensity2 = 0;
   for ( int p = 0; p < npts; p++) {
       temp_intensity += points[p][3]*cos(inner_product(Q,points[p][0..1]));
   }
   intensity[n] = pow(temp_intensity,2)+pow(temp_intensity_2,2)
}
