__kernel void ssandPile_ocl_term(__global unsigned *in, __global unsigned *out, __global unsigned *buffer)
{
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

static uint countIter = 0;
static uint checkTermIterm = 20;

err |= clSetKernelArg(compute_kernel, 2, sizeof(cl_mem), &term_buffer);

if (countIter > checkTermIterm) {
  err = clEnqueueReadBuffer(queue, term_buffer, CL_TRUE, 0, DIM * DIM * sizeof(TYPE), tmpTab, 0, NULL, NULL);
  check(err, "Failed to read buffer from GPU");

  bool notChange = true;
  for (size_t i = 0; i < DIM * DIM; i++)
    if (tmpTab[i] != 0) {
      notChange = false;
      break;
    }

  if (notChange) {
    ret = 1;
    break;
  }
  countIter = 1;
} else countIter++;

