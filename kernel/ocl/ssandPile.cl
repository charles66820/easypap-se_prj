#include "kernel/ocl/common.cl"

__kernel void ssandPile_ocl(__global unsigned *in, __global unsigned *out)
{
  int x = get_global_id(0);
  int y = get_global_id(1);

  //traitement
  int posX = x; //-xloc;
  int posY = y; //-yloc;

  int myPos = posY * DIM + posX;

  if (!(x == 0 || y == 0 || x == DIM - 1 || y == DIM - 1))
  {
    unsigned result = in[myPos] % 4;
    result += in[myPos + DIM] / 4;
    result += in[myPos - DIM] / 4;
    result += in[myPos + 1] / 4;
    result += in[myPos - 1] / 4;
    out[myPos] = result;
  }
}

static void ssandPile(__global unsigned *in, __global unsigned *out, __global unsigned *buffer, unsigned offset)
{
  int x = get_global_id(0);
  int y = get_global_id(1) + offset;


  int myPos = y * DIM + x;

  if (!(x == 0 || y == 0 || x == DIM - 1 || y == DIM - 1))
  {
    unsigned result = in[myPos] % 4;
    result += in[myPos + DIM] / 4;
    result += in[myPos - DIM] / 4;
    result += in[myPos + 1] / 4;
    result += in[myPos - 1] / 4;
    out[myPos] = result;
    buffer[myPos] = result != in[myPos];
  }
}

__kernel void ssandPile_ocl_term(__global unsigned *in, __global unsigned *out, __global unsigned *buffer)
{
  ssandPile(in, out, buffer, 0);
}

__kernel void ssandPile_ocl_omp(__global unsigned *in, __global unsigned *out, __global unsigned *buffer, unsigned offset)
{
  ssandPile(in, out, buffer, offset);
}

// DO NOT MODIFY: this kernel updates the OpenGL texture buffer
// This is a ssandPile-specific version (generic version is defined in common.cl)
__kernel void ssandPile_update_texture(__global unsigned *cur, __write_only image2d_t tex)
{
  int y      = get_global_id(1);
  int x      = get_global_id(0);
  int2 pos   = (int2) (x, y);
  unsigned c = cur[y * DIM + x];
  unsigned r = 0, v = 0, b = 0;

  if (c == 1)
    v = 255;
  else if (c == 2)
    b = 255;
  else if (c == 3)
    r = 255;
  else if (c == 4)
    r = v = b = 255;
  else if (c > 4)
    r = v = b = (2 * c);

  c = rgba(r, v, b, 0xFF);
  write_imagef(tex, pos, color_scatter(c));
}
