#!/usr/bin/python

import numpy as np
import pyopencl as cl

class OpenCL:
   def __init__(self):
      '''Initilize OpenCL.'''
      self.context=cl.create_some_context()
      self.queue=cl.CommandQueue(self.context)

   def load_program(self,filename):
      '''Loads OpenCL kernel from file.'''
      f = open(filename,'r')
      fstr = "".join(f.readlines())
      self.program=cl.Program(self.context,fstr).build()
      
   def sumint(self,qsize,ehc,pixels,points,sym=0,small=0):
      '''Returns the normalized intensity.'''
      #Copy arguments into memory.
      #mf = cl.mem_flags
      npts=points.shape()
      buf_qsize = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.float32(qsize))
      if sym == 0 or small == 0:
         buf_ehc = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.float32(ehc))
      buf_pixels = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.int32(pixels))
      buf_points = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=points)
      buf_npts = cl.Buffer(self.context,cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR, hostbuf=np.int32(npts))
      out=np.empty(pixels**2)
      #out=np.empty(pixels**2).reshape(pixels,pixels)
      out_buffer = cl.Buffer(self.context,cl.mem_flags.WRITE_ONLY,out.nbytes)   #for the output

      if sym == 0:
         self.program.sumint00(self.queue,out.shape,None,buf_qsize,buf_ehc,buf_pixels,buf_points,buf_npts,out_buffer)
      elif small == 0:
         self.program.sumint10(self.queue,out.shape,None,buf_qsize,buf_ehc,buf_pixels,buf_points,buf_npts,out_buffer)
      else:
         self.program.sumint11(self.queue,out.shape,None,buf_qsize,buf_pixels,buf_points,buf_npts,out_buffer)
   def exe(self):
      cl.enqueue_read_buffer(self.queue,out_buffer,out).wait()

      return out/np.sum(out).reshape(pixels,-1)    #normalizes and converts to square

#no mask; implement mask so x,y coordinates of relevant pixels are passed in.
#no separate x_,y_pixels

if __name__ == "__main__":
   instance = OpenCL()
   instance.load_program('sumint.cl')
   intensity = instance.sumint(qsize,ehc,pixels,points)
