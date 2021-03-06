#include "easypap.h"

#include <omp.h>
#include <stdbool.h>
#include <sys/mman.h>
#include <unistd.h>

#pragma region init

typedef unsigned int TYPE;

static TYPE *restrict TABLE = NULL;
static TYPE *TILED_TABLE1   = NULL;
static TYPE *TILED_TABLE2   = NULL;

static inline TYPE *atable_cell(TYPE *restrict i, int y, int x)
{
  return i + y * DIM + x;
}

#define atable(y, x) (*atable_cell(TABLE, (y), (x)))

static inline TYPE *table_cell(TYPE *restrict i, int step, int y, int x)
{
  return DIM * DIM * step + i + y * DIM + x;
}

#define table(step, y, x) (*table_cell(TABLE, (step), (y), (x)))

static int in  = 0;
static int out = 1;

static inline void swap_tables()
{
  int tmp = in;
  in      = out;
  out     = tmp;
}

#define RGB(r, g, b) rgba(r, g, b, 0xFF)

static TYPE max_grains;

void asandPile_refresh_img()
{
  unsigned long int max = 0;
  for (int i = 1; i < DIM - 1; i++)
    for (int j = 1; j < DIM - 1; j++)
    {
      int g = table(in, i, j);
      int r, v, b;
      r = v = b = 0;
      if (g == 1)
        v = 255;
      else if (g == 2)
        b = 255;
      else if (g == 3)
        r = 255;
      else if (g == 4)
        r = v = b = 255;
      else if (g > 4)
        r = b = 255 - (240 * ((double) g) / (double) max_grains);

      cur_img(i, j) = RGB(r, v, b);
      if (g > max)
        max = g;
    }
  max_grains = max;
}

/////////////////////////////  Initial Configurations

static inline void set_cell(int y, int x, unsigned v)
{
  atable(y, x) = v;
  if (opencl_used)
    cur_img(y, x) = v;
}

void asandPile_draw_4partout(void);

void asandPile_draw(char *param)
{
  // Call function ${kernel}_draw_${param}, or default function (second
  // parameter) if symbol not found
  hooks_draw_helper(param, asandPile_draw_4partout);
}

void ssandPile_draw(char *param)
{
  hooks_draw_helper(param, asandPile_draw_4partout);
}

void asandPile_draw_4partout(void)
{
  max_grains = 8;
  for (int i = 1; i < DIM - 1; i++)
    for (int j = 1; j < DIM - 1; j++)
      set_cell(i, j, 4);
}

void asandPile_draw_DIM(void)
{
  max_grains = DIM;
  for (int i = DIM / 4; i < DIM - 1; i += DIM / 4)
    for (int j = DIM / 4; j < DIM - 1; j += DIM / 4)
      set_cell(i, j, i * j / 4);
}

void asandPile_draw_alea(void)
{
  max_grains = 5000;
  for (int i = 0; i < DIM >> 3; i++)
  {
    set_cell(1 + random() % (DIM - 2), 1 + random() % (DIM - 2), 1000 + (random() % (4000)));
  }
}

void asandPile_draw_big(void)
{
  const int i = DIM / 2;
  set_cell(i, i, 100000);
}

static void one_spiral(int x, int y, int step, int turns)
{
  int i = x, j = y, t;

  for (t = 1; t <= turns; t++)
  {
    for (; i < x + t * step; i++)
      set_cell(i, j, 3);
    for (; j < y + t * step + 1; j++)
      set_cell(i, j, 3);
    for (; i > x - t * step - 1; i--)
      set_cell(i, j, 3);
    for (; j > y - t * step - 1; j--)
      set_cell(i, j, 3);
  }
  set_cell(i, j, 4);

  for (int i = -2; i < 3; i++)
    for (int j = -2; j < 3; j++)
      set_cell(i + x, j + y, 3);
}

static void many_spirals(int xdebut, int xfin, int ydebut, int yfin, int step, int turns)
{
  int i, j;
  int size = turns * step + 2;

  for (i = xdebut + size; i < xfin - size; i += 2 * size)
    for (j = ydebut + size; j < yfin - size; j += 2 * size)
      one_spiral(i, j, step, turns);
}

static void spiral(unsigned twists)
{
  many_spirals(1, DIM - 2, 1, DIM - 2, 2, twists);
}

void asandPile_draw_spirals(void)
{
  spiral(DIM / 32);
}

// shared functions

#define ALIAS(fun)       \
  void ssandPile_##fun() \
  {                      \
    asandPile_##fun();   \
  }

ALIAS(refresh_img);
ALIAS(draw_4partout);
ALIAS(draw_DIM);
ALIAS(draw_alea);
ALIAS(draw_big);
ALIAS(draw_spirals);

#pragma endregion init

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Synchronous Kernel
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

#pragma region synchronousKernel

#pragma region ssandInit

void ssandPile_init()
{
  TABLE = calloc(2 * DIM * DIM, sizeof(TYPE));
}

void ssandPile_finalize()
{
  free(TABLE);
}

int ssandPile_do_tile_default(int x, int y, int width, int height)
{
  int diff = 0;

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
    {
      table(out, i, j) = table(in, i, j) % 4;
      table(out, i, j) += table(in, i + 1, j) / 4;
      table(out, i, j) += table(in, i - 1, j) / 4;
      table(out, i, j) += table(in, i, j + 1) / 4;
      table(out, i, j) += table(in, i, j - 1) / 4;
      if (table(out, i, j) != table(in, i, j))
        diff = 1;
    }

  return diff;
}

int ssandPile_do_tile_opt(int x, int y, int width, int height)
{
  int diff = 0;

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
    {
      int result = table(in, i, j) % 4;
      result += table(in, i + 1, j) / 4;
      result += table(in, i - 1, j) / 4;
      result += table(in, i, j + 1) / 4;
      result += table(in, i, j - 1) / 4;
      table(out, i, j) = result;
      diff |= result != table(in, i, j);
    }

  return diff;
}

#pragma endregion ssandInit

#pragma region ssandSeq

/////////////////////////////  Sequential version (seq)
// Suggested cmdline(s):
// ./run -k ssandPile -v seq -s 512 -m
//
// ./run -k ssandPile -v seq -wt opt -s 512 -m
//
// Renvoie le nombre d'it??rations effectu??es avant stabilisation, ou 0
unsigned ssandPile_compute_seq(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = do_tile(1, 1, DIM - 2, DIM - 2, 0);
    swap_tables();
    if (change == 0)
      return it;
  }
  return 0;
}

/////////////////////////////  Tiled sequential version (tiled)
// Suggested cmdline(s):
// ./run -k ssandPile -v tiled -s 512 -m
//
// ./run -k ssandPile -v tiled -wt opt -s 512 -m
//
unsigned ssandPile_compute_tiled(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = 0;

    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
        change |= do_tile(x + (x == 0),
                          y + (y == 0),
                          TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                          TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                          0 /* CPU id */);
    swap_tables();
    if (change == 0)
      return it;
  }

  return 0;
}

#pragma endregion ssandSeq

#pragma region ssandOmp

/////////////////////////////  Sequential version (omp)
// Suggested cmdline(s):
// ./run -k ssandPile -v omp -s 512 -m
//
// ./run -k ssandPile -v omp -wt opt -s 512 -m
//
// Renvoie le nombre d'it??rations effectu??es avant stabilisation, ou 0
unsigned ssandPile_compute_omp(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = 0;
#pragma omp parallel for schedule(runtime) reduction(| : change)
    for (int y = 1; y < DIM - 1; y += 1)
      for (int x = 1; x < DIM - 1; x += 1)
        change |= do_tile(x, y, 1, 1, omp_get_thread_num());
    swap_tables();
    if (change == 0)
      return it;
  }

  return 0;
}

/////////////////////////////  Tiled parallale version (omp_tiled)
// Suggested cmdline(s):
// ./run -k ssandPile -v omp_tiled -s 512 -m
//
// ./run -k ssandPile -v omp_tiled -wt opt -s 512 -m
//
unsigned ssandPile_compute_omp_tiled(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = 0;

#pragma omp parallel for schedule(runtime) collapse(2) reduction(| : change)
    for (int y = 0; y < DIM - 1; y += TILE_H)
      for (int x = 0; x < DIM - 1; x += TILE_W)
        change |= do_tile(x + (x == 0),
                          y + (y == 0),
                          TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                          TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                          omp_get_thread_num());
    swap_tables();
    if (change == 0)
      return it;
  }

  return 0;
}

/////////////////////////////  Tiled parallale version (omp_tiled)
// Suggested cmdline(s):
// ./run -k ssandPile -v omp_taskloop -s 512 -m
//
// ./run -k ssandPile -v omp_taskloop -wt opt -s 512 -m
//
unsigned ssandPile_compute_omp_taskloop(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = 0;

#pragma omp parallel
#pragma omp single
    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
#pragma omp task shared(change)
      {
        int localChange = do_tile(x + (x == 0),
                                  y + (y == 0),
                                  TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                  TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                  omp_get_thread_num());
        if (change == 0 && localChange != 0)
          change = 1;
      }
    swap_tables();
    if (change == 0)
      return it;
  }

  return 0;
}

#pragma endregion ssandOmp

#pragma region ssandLazy

static int tt = 0;

static inline void swap_tt()
{
  tt = tt ^ 1; // tt = !tt;
}

void ssandPile_init_lazy()
{
  ssandPile_init();
  int nbCell   = NB_TILES_Y * (NB_TILES_X + 1);
  TILED_TABLE1 = malloc(nbCell * sizeof(TYPE));
  for (TYPE i = 0; i < nbCell; i++)
    TILED_TABLE1[i] = 1;

  TILED_TABLE2 = calloc(nbCell, sizeof(TYPE));
}

void ssandPile_finalize_lazy()
{
  ssandPile_finalize();
  free(TILED_TABLE1);
  free(TILED_TABLE2);
}

#define tiled_table1(y, x) (TILED_TABLE1[((y) *NB_TILES_X) + (x)])
#define tiled_table2(y, x) (TILED_TABLE2[((y) *NB_TILES_X) + (x)])

/////////////////////////////  Tiled sequential version (tiled lazy)
// Suggested cmdline(s):
// ./run -k ssandPile -v lazy -s 512 -m
//
// ./run -k ssandPile -v lazy -wt opt -s 512 -m
//
unsigned ssandPile_compute_lazy(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = 0;

    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
      {
        int localChange = 0;
        int ty          = y / TILE_H;
        int tx          = x / TILE_W;

        if ((tt == 0 &&
             ((ty != 0 && tiled_table1(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table1(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table1(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table1(ty, tx + 1) == 1) ||
              tiled_table1(ty, tx) == 1)) ||
            (tt == 1 &&
             ((ty != 0 && tiled_table2(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table2(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table2(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table2(ty, tx + 1) == 1) ||
              tiled_table2(ty, tx) == 1)))
        {
          localChange = do_tile(x + (x == 0),
                                y + (y == 0),
                                TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                0 /* CPU id */);
        }

        if (tt == 0)
          tiled_table2(ty, tx) = localChange;
        else
          tiled_table1(ty, tx) = localChange;
        change |= localChange;
      }
    swap_tables();
    swap_tt();
    if (change == 0)
      return it;
  }

  return 0;
}

void ssandPile_init_omp_lazy()
{
  ssandPile_init_lazy();
}
void ssandPile_finalize_omp_lazy()
{
  ssandPile_finalize_lazy();
}
/////////////////////////////  Tiled sequential version (tiled lazy)
// Suggested cmdline(s):
// ./run -k ssandPile -v omp_lazy -s 512 -m
//
// ./run -k ssandPile -v omp_lazy -wt opt -s 512 -m
//
unsigned ssandPile_compute_omp_lazy(unsigned nb_iter)
{
  int res = 0;

  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = 0;

#pragma omp parallel for schedule(runtime) collapse(2) reduction(| : change)
    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
      {
        int localChange = 0;
        int ty          = y / TILE_H;
        int tx          = x / TILE_W;

        if ((tt == 0 &&
             ((ty != 0 && tiled_table1(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table1(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table1(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table1(ty, tx + 1) == 1) ||
              tiled_table1(ty, tx) == 1)) ||
            (tt == 1 &&
             ((ty != 0 && tiled_table2(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table2(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table2(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table2(ty, tx + 1) == 1) ||
              tiled_table2(ty, tx) == 1)))
        {
          localChange = do_tile(x + (x == 0),
                                y + (y == 0),
                                TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                omp_get_thread_num());
        }

        if (tt == 0)
          tiled_table2(ty, tx) = localChange;
        else
          tiled_table1(ty, tx) = localChange;
        change |= localChange;
      }
    swap_tables();
    swap_tt();
    if (change == 0)
    {
      res = it;
      break;
    }
  }

  return res;
}
#pragma endregion ssandLazy

#pragma region ssandAVX
// Intrinsics functions
#ifdef ENABLE_VECTO

#if __AVX2__ == 1

#include <immintrin.h>

void ssandPile_tile_check_avx (void)
{
  // Tile width must be larger than AVX vector size
  easypap_vec_check (AVX_VEC_SIZE_INT, DIR_HORIZONTAL);
}

int ssandPile_do_tile_avx(int x, int y, int width, int height)
{
  const __m256i vec3_i = _mm256_set1_epi32(3);

  // Outer tiles are computed the usual way
  // if (x == 1 || x == (DIM - 1) - width || y == 1 || y == (DIM - 1) - height)
  //   return ssandPile_do_tile_opt(x, y, width, height);
  if (y == 1 || y == (DIM - 1) - height)
    return ssandPile_do_tile_opt(x, y, width, height);
  if (x == (DIM - 1) - width)
    x -= 1;

  // Inner tiles involve no border test
  int diff = 0;
  // We travel AVX_VEC_SIZE_INT and not 1 by on
  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j += AVX_VEC_SIZE_INT)
    {
      __m256i result_i;
      __m256i tmp;

      // load table(in, i, j)
      __m256i currentPixelsRow_i = _mm256_loadu_si256((__m256i *) &table(in, i, j));
      // load table(in, i + 1, j)
      __m256i bottomPixelsRow_i = _mm256_loadu_si256((__m256i *) &table(in, i + 1, j));
      // load table(in, i - 1, j)
      __m256i topPixelsRow_i = _mm256_loadu_si256((__m256i *) &table(in, i - 1, j));
      // load table(in, i, j + 1)
      __m256i rightPixelsRow_i = _mm256_loadu_si256((__m256i *) &table(in, i, j + 1));
      // load table(in, i, j - 1)
      __m256i leftPixelsRow_i = _mm256_loadu_si256((__m256i *) &table(in, i, j - 1));

      // result = currentPixelsRow % 4;
      result_i = _mm256_and_si256(currentPixelsRow_i, vec3_i); // currentPixelsRow & (4 - 1)

      // result += topPixelsRow / 4;
      tmp = _mm256_srli_epi32(topPixelsRow_i, 2);
      result_i = _mm256_add_epi32(result_i, tmp);

      // result += bottomPixelsRow / 4;
      tmp = _mm256_srli_epi32(bottomPixelsRow_i, 2);
      result_i = _mm256_add_epi32(result_i, tmp);

      // result += rightPixelsRow / 4;
      tmp = _mm256_srli_epi32(rightPixelsRow_i, 2);
      result_i = _mm256_add_epi32(result_i, tmp);

      // result += leftPixelsRow / 4;
      tmp = _mm256_srli_epi32(leftPixelsRow_i, 2);
      result_i = _mm256_add_epi32(result_i, tmp);

      // table(out, i, j) = result;
      _mm256_storeu_si256((__m256i *) &table(out, i, j), result_i);

      // diff |= result != currentPixelsRow;
      __m256i mask = _mm256_xor_si256(result_i, currentPixelsRow_i);

      // all 0 => res = 0 : diff no change
      // exists 1 => res = 1 : diff change to 1
      // All the mask is not equal to zero
      if (_mm256_testz_si256(mask, mask) != 1)
        diff = 1;
    }

  return diff;
}

#endif
#endif
#pragma endregion ssandAVX

#pragma region ssandOpenCL

// Only called when --dump or --thumbnails is used
void ssandPile_refresh_img_ocl()
{
  cl_int err;

  err = clEnqueueReadBuffer(queue, cur_buffer, CL_TRUE, 0, sizeof(unsigned) * DIM * DIM, TABLE, 0, NULL, NULL);
  check(err, "Failed to read buffer from GPU");

  ssandPile_refresh_img();
}


void ssandPile_refresh_img_ocl_term()
{
  ssandPile_refresh_img_ocl();
}

static cl_mem term_buffer;
void ssandPile_init_ocl_term(void)
{
  ssandPile_init();

  const int size = DIM * DIM * sizeof(TYPE);

  term_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, NULL);
  if (!term_buffer)
    exit_with_error ("Failed to allocate termination buffer");
}

static uint countIter = 0;
static uint checkTermIterm = 20;
unsigned ssandPile_invoke_ocl_term(unsigned nb_iter)
{
  TYPE * tmpTab = calloc(DIM * DIM, sizeof(TYPE));
  uint ret = 0;

  size_t global[2] = {GPU_SIZE_X, GPU_SIZE_Y}; // global domain size for our calculation
  size_t local[2]  = {GPU_TILE_W, GPU_TILE_H}; // local domain size for our calculation
  cl_int err;

  monitoring_start_tile(easypap_gpu_lane(TASK_TYPE_COMPUTE));

  for (unsigned it = 1; it <= nb_iter; it++) {

    // Set kernel arguments
    //
    err = 0;
    err |= clSetKernelArg(compute_kernel, 0, sizeof(cl_mem), &cur_buffer);
    err |= clSetKernelArg(compute_kernel, 1, sizeof(cl_mem), &next_buffer);
    err |= clSetKernelArg(compute_kernel, 2, sizeof(cl_mem), &term_buffer);
    check(err, "Failed to set kernel arguments");

    err = clEnqueueNDRangeKernel(queue, compute_kernel, 2, NULL, global, local,
                                  0, NULL, NULL);
    check(err, "Failed to execute kernel");


    // Swap buffers
    {
      cl_mem tmp  = cur_buffer;
      cur_buffer  = next_buffer;
      next_buffer = tmp;
    }

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
  }

  clFinish(queue);

  free(tmpTab);

  monitoring_end_tile(0, 0, DIM, DIM, easypap_gpu_lane(TASK_TYPE_COMPUTE));

  return ret;
}

#pragma endregion ssandOpenCL


#pragma region ssandOpenCLOpenMp
void ssandPile_refresh_img_ocl_omp()
{
  ssandPile_refresh_img_ocl();
}

// Threashold = 10%
#define THRESHOLD 10

static unsigned cpu_y_part;
static unsigned gpu_y_part;

// static cl_mem term_buffer;
void ssandPile_init_ocl_omp(void)
{
  ssandPile_init();

  const int size = DIM * DIM * sizeof(TYPE);

  term_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, NULL);
  if (!term_buffer)
    exit_with_error ("Failed to allocate termination buffer");

  if (GPU_TILE_H != TILE_H)
    exit_with_error ("CPU and GPU Tiles should have the same height (%d != %d)",
                     GPU_TILE_H, TILE_H);

  cpu_y_part = (NB_TILES_Y / 2) * GPU_TILE_H; // Start with fifty-fifty
  gpu_y_part = DIM - cpu_y_part;
}

static long gpu_duration = 0, cpu_duration = 0;

static int much_greater_than (long t1, long t2)
{
  return (t1 > t2) && ((t1 - t2) * 100 / t1 > THRESHOLD);
}

// static uint countIter = 0;
// static uint checkTermIterm = 1;

#define DIMENSION 2
#define LOAD_BALANCING false

// Suggested cmdline:
// ./run -k ssandPile -o -v ocl_omp -m -ts 16
//
unsigned ssandPile_invoke_ocl_omp(unsigned nb_iter)
{
  // Vars
  TYPE * tmpTab = calloc(DIM * DIM, sizeof(TYPE));
  uint ret = 0;

  size_t global[DIMENSION] = {GPU_SIZE_X, gpu_y_part}; // global domain size for our calculation
  size_t local[DIMENSION]  = {GPU_TILE_W, GPU_TILE_H}; // local domain size for our calculation
  cl_int err;

  cl_event kernel_event;
  long startTime, endTime;
  int gpu_accumulated_lines = 0;

  // Iters
  for (unsigned it = 1; it <= nb_iter; it++) {

    // Load balancing
    if (LOAD_BALANCING && gpu_duration != 0) {
      if (much_greater_than (gpu_duration, cpu_duration) &&
          gpu_y_part > GPU_TILE_H) {
        gpu_y_part -= GPU_TILE_H;
        cpu_y_part += GPU_TILE_H;
        global[1] = gpu_y_part;
      } else if (much_greater_than (cpu_duration, gpu_duration) &&
                 cpu_y_part > GPU_TILE_H) {
        gpu_y_part += GPU_TILE_H;
        cpu_y_part -= GPU_TILE_H;
        global[1] = gpu_y_part;
      }
    }

    // Set kernel arguments
    //
    err = 0;
    err |= clSetKernelArg(compute_kernel, 0, sizeof(cl_mem), &cur_buffer);
    err |= clSetKernelArg(compute_kernel, 1, sizeof(cl_mem), &next_buffer);
    err |= clSetKernelArg(compute_kernel, 2, sizeof(cl_mem), &term_buffer);
    err |= clSetKernelArg (compute_kernel, 3, sizeof (unsigned), &cpu_y_part);
    check(err, "Failed to set kernel arguments");

    // Launch GPU kernel
    err = clEnqueueNDRangeKernel(queue, compute_kernel, DIMENSION, NULL, global, local,
                                  0, NULL, &kernel_event);
    check(err, "Failed to execute kernel");
    clFlush(queue); // submit all commands

    // Compute CPU part
    startTime = what_time_is_it();

    {
      int change = 0;

#pragma omp parallel for collapse(2) schedule(runtime) reduction(| : change)
        for (int y = 0; y < cpu_y_part; y += TILE_H)
          for (int x = 0; x < DIM - 1; x += TILE_W)
            change |= do_tile(x + (x == 0), y + (y == 0),
                          TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                          TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                          omp_get_thread_num());
      ret = !change;
    }

    endTime = what_time_is_it ();
    cpu_duration = endTime - startTime;

    // Wait all GPU command finish
    clFinish(queue);

    // Monitor
    gpu_duration = ocl_monitor(kernel_event, 0, cpu_y_part, global[0],
                                global[1], TASK_TYPE_COMPUTE);
    clReleaseEvent(kernel_event);

    gpu_accumulated_lines += gpu_y_part;

    // Swap buffers
    {
      cl_mem tmp  = cur_buffer;
      cur_buffer  = next_buffer;
      next_buffer = tmp;
    }
    swap_tables();

    // termination
    if (countIter > checkTermIterm)
    {
      err = clEnqueueReadBuffer(queue, term_buffer, CL_TRUE, 0, DIM * DIM * sizeof(TYPE), tmpTab, 0, NULL, NULL);
      check(err, "Failed to read buffer from GPU");

      bool notChange = true;
      for (size_t i = DIM * cpu_y_part; i < DIM * DIM; i++)
        if (tmpTab[i] != 0)
        {
          notChange = false;
          break;
        }

      if (notChange)
      {
        ret = 1;
        break;
      }
      countIter = 1;
    }
    else countIter++;
  }

  // clFinish(queue);

  if (do_display) {
    // Send CPU contribution to GPU memory
    err = clEnqueueWriteBuffer(queue, cur_buffer, CL_TRUE, 0,
                                DIM * cpu_y_part * sizeof (TYPE), &table(in,0,0), 0,
                                NULL, NULL);
    check(err, "Failed to write to buffer");
  } else
    PRINT_DEBUG ('u', "In average, GPU took %.1f%% of the lines\n",
                 (float)gpu_accumulated_lines * 100 / (DIM * nb_iter));

  free(tmpTab);

  return ret;
}

void printTableDebug() {
  for (int i = 0; i < DIM; i++)
  {
    for (int j = 0; j < DIM; j++)
      printf("%d", TABLE[i * DIM + j]);
    printf("\n");
  }
  printf("\033[0;0H");

  for (int i = 0; i < DIM; i++)
  {
    for (int j = 0; j < DIM; j++)
      printf("\033[%d;%dH%d\n", i + 5, DIM + 4 + j, TABLE[i * DIM + j]);
    printf("\n");
  }
}
#pragma endregion ssandOpenCLOpenMp

#pragma endregion synchronousKernel

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Asynchronous Kernel
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

#pragma region asynchronousKernel

#pragma region asandInit

void asandPile_init()
{
  in = out = 0;
  if (TABLE == NULL)
  {
    const unsigned size = DIM * DIM * sizeof(TYPE);

    PRINT_DEBUG('u', "Memory footprint = 2 x %d bytes\n", size);

    TABLE = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  }
}

void asandPile_finalize()
{
  const unsigned size = DIM * DIM * sizeof(TYPE);

  munmap(TABLE, size);
}

///////////////////////////// Version s??quentielle simple (seq)
// Renvoie le nombre d'it??rations effectu??es avant stabilisation, ou 0

int asandPile_do_tile_default(int x, int y, int width, int height)
{
  int change = 0;

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
      if (atable(i, j) >= 4)
      {
        atable(i, j - 1) += atable(i, j) / 4;
        atable(i, j + 1) += atable(i, j) / 4;
        atable(i - 1, j) += atable(i, j) / 4;
        atable(i + 1, j) += atable(i, j) / 4;
        atable(i, j) %= 4;
        change = 1;
      }
  return change;
}

int asandPile_do_tile_opt(int x, int y, int width, int height)
{
  int change = 0;

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
    {
      int result = atable(i, j);
      if (result >= 4)
      {
        result /= 4;
        atable(i, j - 1) += result;
        atable(i, j + 1) += result;
        atable(i - 1, j) += result;
        atable(i + 1, j) += result;
        atable(i, j) %= 4;
        change = 1;
      }
    }
  return change;
}

#pragma endregion asandInit

#pragma region asandSeq

/////////////////////////////  Sequential version (seq)
// Suggested cmdline(s):
// ./run -k asandPile -v seq -s 512 -m
//
// ./run -k asandPile -v seq -wt opt -s 512 -m
//
unsigned asandPile_compute_seq(unsigned nb_iter)
{
  int change = 0;
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    // On traite toute l'image en un coup (oui, c'est une grosse tuile)
    change = do_tile(1, 1, DIM - 2, DIM - 2, 0);

    if (change == 0)
      return it;
  }
  return 0;
}

/////////////////////////////  Tiled sequential version (tiled)
// Suggested cmdline(s):
// ./run -k asandPile -v tiled -s 512 -m
//
// ./run -k asandPile -v tiled -wt opt -s 512 -m
//
unsigned asandPile_compute_tiled(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = 0;

    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
        change |= do_tile(x + (x == 0),
                          y + (y == 0),
                          TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                          TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                          0 /* CPU id */);
    if (change == 0)
      return it;
  }

  return 0;
}

#pragma endregion asandSeq

#pragma region asandOmp

/////////////////////////////  Tiled sequential version (tiled)
// Suggested cmdline(s):
// ./run -k asandPile -v omp_tiled -s 512 -m
//
// ./run -k asandPile -v omp_tiled -wt opt -s 512 -m
//
unsigned asandPile_compute_omp_tiled(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = 0;

    //BLEU
#pragma omp parallel for schedule(runtime) reduction(| : change)
    for (int y = 0; y < DIM; y += 2 * TILE_H)
    {
      for (int x = y % (TILE_H * 2); x < DIM; x += TILE_W * 2)
      {
        change |= do_tile(x + (x == 0),
                          y + (y == 0),
                          TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                          TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                          omp_get_thread_num());
      }
    }

    //ROUGE
#pragma omp parallel for schedule(runtime) reduction(| : change)
    for (int y = 0; y < DIM; y += 2 * TILE_H)
    {
      for (int x = (y + TILE_H) % (TILE_H * 2); x < DIM; x += TILE_W * 2)
      {
        change |= do_tile(x + (x == 0),
                          y + (y == 0),
                          TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                          TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                          omp_get_thread_num());
      }
    }

    //VERT
#pragma omp parallel for schedule(runtime) reduction(| : change)
    for (int y = TILE_H; y < DIM; y += 2 * TILE_H)
    {
      for (int x = y % (TILE_H * 2); x < DIM; x += TILE_W * 2)
      {
        change |= do_tile(x + (x == 0),
                          y + (y == 0),
                          TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                          TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                          omp_get_thread_num());
      }
    }

    //NOIR
#pragma omp parallel for schedule(runtime) reduction(| : change)
    for (int y = TILE_H; y < DIM; y += 2 * TILE_H)
    {
      for (int x = (y + TILE_H) % (TILE_H * 2); x < DIM; x += TILE_W * 2)
      {
        change |= do_tile(x + (x == 0),
                          y + (y == 0),
                          TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                          TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                          omp_get_thread_num());
      }
    }

    if (change == 0)
      return it;
  }
  return 0;
}

/////////////////////////////  Tiled sequential version (tiled)
// Suggested cmdline(s):
// ./run -k asandPile -v omp -s 512 -m
//
// ./run -k asandPile -v omp_task -wt opt -s 512 -m
//
unsigned asandPile_compute_omp_task(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = 0;
    int tuile[NB_TILES_Y][NB_TILES_X + 1] __attribute__((unused));

#pragma omp parallel
#pragma omp single
    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
#pragma omp task depend(in                                                                                   \
                        : tuile[x / TILE_W][(y / TILE_H) + 1]) depend(in                                     \
                                                                      : tuile[(x / TILE_W) - 1][y / TILE_H]) \
    depend(in                                                                                                \
           : tuile[(x / TILE_W) - 1][(y / TILE_H) + 1]) depend(out                                           \
                                                               : tuile[x / TILE_W][y / TILE_H]) shared(change)
      {
        int localChange = do_tile(x + (x == 0),
                                  y + (y == 0),
                                  TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                  TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                  omp_get_thread_num());
        if (change == 0 && localChange != 0)
          change = 1;
      }
#pragma omp taskwait
    if (!change)
      return it;
  }
  return 0;
}

#pragma endregion asandOmp

#pragma region asandLazy

void asandPile_init_lazy()
{
  ssandPile_init_lazy();
}
void asandPile_finalize_lazy()
{
  ssandPile_finalize_lazy();
}
/////////////////////////////  Lazy version (tiled)
// Suggested cmdline(s):
// ./run -k asandPile -v lazy -s 512 -m
//
// ./run -k asandPile -v lazy -wt opt -s 512 -m
//
unsigned asandPile_compute_lazy(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = 0;

    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
      {
        int localChange = 0;
        int ty          = y / TILE_H;
        int tx          = x / TILE_W;

        if ((tt == 0 &&
             ((ty != 0 && tiled_table1(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table1(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table1(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table1(ty, tx + 1) == 1) ||
              tiled_table1(ty, tx) == 1)) ||
            (tt == 1 &&
             ((ty != 0 && tiled_table2(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table2(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table2(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table2(ty, tx + 1) == 1) ||
              tiled_table2(ty, tx) == 1)))
        // if ((tt == 0 && tiled_table1(ty, tx)) ||
        //   (tt == 1 && tiled_table2(ty, tx)))
        {
          localChange = do_tile(x + (x == 0),
                                y + (y == 0),
                                TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                0 /* CPU id */);
        }

        if (tt == 0)
          tiled_table2(ty, tx) = localChange;
        else
          tiled_table1(ty, tx) = localChange;
        change |= localChange;
      }

    swap_tt();
    if (change == 0)
      return it;
  }

  return 0;
}

void asandPile_init_omp_lazy()
{
  asandPile_init_lazy();
}
void asandPile_finalize_omp_lazy()
{
  asandPile_finalize_lazy();
}
/////////////////////////////  Lazy version (tiled)
// Suggested cmdline(s):
// ./run -k asandPile -v omp_lazy -s 512 -m
//
// ./run -k asandPile -v omp_lazy -wt opt -s 512 -m
//
unsigned asandPile_compute_omp_lazy(unsigned nb_iter)
{
  int res = 0;

  for (unsigned it = 1; it <= nb_iter; it++)
  {
    int change = 0;

    //BLEU
#pragma omp parallel for schedule(runtime) reduction(| : change)
    for (int y = 0; y < DIM; y += 2 * TILE_H)
    {
      for (int x = y % (TILE_H * 2); x < DIM; x += TILE_W * 2)
      {
        int localChange = 0;
        int ty          = y / TILE_H;
        int tx          = x / TILE_W;

        if ((tt == 0 &&
             ((ty != 0 && tiled_table1(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table1(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table1(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table1(ty, tx + 1) == 1) ||
              tiled_table1(ty, tx) == 1

              )) ||
            (tt == 1 &&
             ((ty != 0 && tiled_table2(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table2(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table2(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table2(ty, tx + 1) == 1) ||
              tiled_table2(ty, tx) == 1

              )))
        {
          localChange = do_tile(x + (x == 0),
                                y + (y == 0),
                                TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                omp_get_thread_num());
        }

        if (tt == 0)
          tiled_table2(ty, tx) = localChange;
        else
          tiled_table1(ty, tx) = localChange;
        change |= localChange;
      }
    }

    //ROUGE
#pragma omp parallel for schedule(runtime) reduction(| : change)
    for (int y = 0; y < DIM; y += 2 * TILE_H)
    {
      for (int x = (y + TILE_H) % (TILE_H * 2); x < DIM; x += TILE_W * 2)
      {
        int localChange = 0;
        int ty          = y / TILE_H;
        int tx          = x / TILE_W;

        if ((tt == 0 &&
             ((ty != 0 && tiled_table1(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table1(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table1(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table1(ty, tx + 1) == 1) ||
              tiled_table1(ty, tx) == 1

              )) ||
            (tt == 1 &&
             ((ty != 0 && tiled_table2(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table2(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table2(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table2(ty, tx + 1) == 1) ||
              tiled_table2(ty, tx) == 1

              )))
        {
          localChange = do_tile(x + (x == 0),
                                y + (y == 0),
                                TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                omp_get_thread_num());
        }

        if (tt == 0)
          tiled_table2(ty, tx) = localChange;
        else
          tiled_table1(ty, tx) = localChange;
        change |= localChange;
      }
    }

    //VERT
#pragma omp parallel for schedule(runtime) reduction(| : change)
    for (int y = TILE_H; y < DIM; y += 2 * TILE_H)
    {
      for (int x = y % (TILE_H * 2); x < DIM; x += TILE_W * 2)
      {
        int localChange = 0;
        int ty          = y / TILE_H;
        int tx          = x / TILE_W;

        if ((tt == 0 &&
             ((ty != 0 && tiled_table1(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table1(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table1(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table1(ty, tx + 1) == 1) ||
              tiled_table1(ty, tx) == 1

              )) ||
            (tt == 1 &&
             ((ty != 0 && tiled_table2(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table2(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table2(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table2(ty, tx + 1) == 1) ||
              tiled_table2(ty, tx) == 1

              )))
        {
          localChange = do_tile(x + (x == 0),
                                y + (y == 0),
                                TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                omp_get_thread_num());
        }

        if (tt == 0)
          tiled_table2(ty, tx) = localChange;
        else
          tiled_table1(ty, tx) = localChange;
        change |= localChange;
      }
    }

    //NOIR
#pragma omp parallel for schedule(runtime) reduction(| : change)
    for (int y = TILE_H; y < DIM; y += 2 * TILE_H)
    {
      for (int x = (y + TILE_H) % (TILE_H * 2); x < DIM; x += TILE_W * 2)
      {
        int localChange = 0;
        int ty          = y / TILE_H;
        int tx          = x / TILE_W;

        if ((tt == 0 &&
             ((ty != 0 && tiled_table1(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table1(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table1(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table1(ty, tx + 1) == 1) ||
              tiled_table1(ty, tx) == 1

              )) ||
            (tt == 1 &&
             ((ty != 0 && tiled_table2(ty - 1, tx) == 1) || (ty != NB_TILES_Y - 1 && tiled_table2(ty + 1, tx) == 1) ||
              (tx != 0 && tiled_table2(ty, tx - 1) == 1) || (tx != NB_TILES_X - 1 && tiled_table2(ty, tx + 1) == 1) ||
              tiled_table2(ty, tx) == 1

              )))
        {
          localChange = do_tile(x + (x == 0),
                                y + (y == 0),
                                TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                omp_get_thread_num());
        }

        if (tt == 0)
          tiled_table2(ty, tx) = localChange;
        else
          tiled_table1(ty, tx) = localChange;
        change |= localChange;
      }
    }

    swap_tt();
    if (change == 0)
    {
      res = it;
      break;
    }
  }

  return res;
}

#pragma endregion asandLazy

#pragma region asandAVX
// Intrinsics functions
#ifdef ENABLE_VECTO

#if __AVX2__ == 1 && __AVX512F__ == 1

#include <immintrin.h>

#define DEBUG 0

int asandPile_do_tile_avx(int x, int y, int width, int height)
{
#if DEBUG == 1
  bool debug = true;
  if (debug) printf("\n");
#endif
  if (x == (DIM - 1) - width)
    x -= 1;

  // $$\overrightarrow{X} == vecX$$
  const __m256i vec3_i = _mm256_set1_epi32(3);
  const __m256i vec0_i = _mm256_set1_epi32(0);

  int diff = 0;
  for (int j = y; j < y + height; j++)
    for (int i = x; i < x + width; i += AVX_VEC_SIZE_INT)
    {
      // vecT_{j-1,i} <-- (t_{j-1,i+k}, ..., t_{j-1,i})
      __m256i topVec_i = _mm256_loadu_si256((__m256i *) &table(in, j - 1, i));
      // vecT_{j,i} <-- (t_{j,i+k}, ..., t_{j,i})
      __m256i vec_i = _mm256_loadu_si256((__m256i *) &table(in, j, i)); // load?
      // vecT_{j+1,i} <-- (t_{j+1,i+k}, ..., t_{j+1,i})
      __m256i bottomVec_i = _mm256_loadu_si256((__m256i *) &table(in, j + 1, i));
#if DEBUG == 1
      if (debug) {
        int *ptr = (int*)&topVec_i;
        printf("top %d %d %d %d %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
      if (  debug) {
        int *ptr = (int*)&vec_i;
        printf("cur %d %d %d %d %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
      if (debug) {
        int *ptr = (int*)&bottomVec_i;
        printf("bot %d %d %d %d %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
#endif
      // vecD <-- vec_i / 4
      __m256i vecD = _mm256_srli_epi32(vec_i, 2);
#if DEBUG == 1
      if (debug) {
        int *ptr = (int*)&vecD;
        printf("\nvecD %d %d %d %d %d %d %d %d    ", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
#endif
      // (vecD << 1)
      __m256i vecDShiftLeft = _mm256_alignr_epi32(vec0_i, vecD, 1);
#if DEBUG == 1
      if (debug) {
        int *ptr = (int*)&vecDShiftLeft;
        printf("vecD SL %d %d %d %d %d %d %d %d    ", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
#endif
      // (vecD >> 1)
      __m256i vecDShiftRight = _mm256_alignr_epi32(vecD, vec0_i, 7);
#if DEBUG == 1
      if (debug) {
        int *ptr = (int*)&vecDShiftRight;
        printf("vecD SR %d %d %d %d %d %d %d %d    ", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
#endif
      // vec_i <-- vec_i % 4 + vecDShiftLeft + vecDShiftRight
      __m256i res_vec_i = _mm256_add_epi32(_mm256_and_si256(vec_i, vec3_i),
                               _mm256_add_epi32(vecDShiftLeft, vecDShiftRight));
#if DEBUG == 1
      if (debug) {
        int *ptr = (int*)&res_vec_i;
        printf("res_rec_i %d %d %d %d %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
#endif
      // topVec_i <-- topVec_i + vecD
      topVec_i = _mm256_add_epi32(topVec_i, vecD);
#if DEBUG == 1
      if (debug) {
        int *ptr = (int*)&topVec_i;
        printf("new top %d %d %d %d %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
#endif
      // bottomVec_i <-- bottomVec_i + vecD
      bottomVec_i = _mm256_add_epi32(bottomVec_i, vecD);
#if DEBUG == 1
      if (debug) {
        int *ptr = (int*)&bottomVec_i;
        printf("new bot %d %d %d %d %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
#endif
      // t_{j,i-1} <-- t_{j,i-1} + vecD[0]
      __m256i leftVec_i = _mm256_loadu_si256((__m256i *) &table(in, j, i - 1));
#if DEBUG == 1
      if (debug) {
        int *ptr = (int*)&leftVec_i;
        printf("leftVec %d %d %d %d %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
#endif
      leftVec_i         = _mm256_add_epi32(leftVec_i, vecD);
#if DEBUG == 1
      if (debug) {
        int *ptr = (int*)&leftVec_i;
        printf("new leftVec %d %d %d %d %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
#endif
      _mm256_storeu_si256((__m256i *) &table(out, j, i - 1), leftVec_i);

      // t_{j,i+k+1} <-- t_{j,i+k+1} + vecD[k]  :  k = AVX_VEC_SIZE_INT-1
      __m256i rightVec_i = _mm256_loadu_si256((__m256i *) &table(in, j, i + 1));
#if DEBUG == 1
      if (debug) {
        int *ptr = (int*)&rightVec_i;
        printf("rightVec %d %d %d %d %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
#endif
      rightVec_i         = _mm256_add_epi32(rightVec_i, vecD);
#if DEBUG == 1
      if (debug) {
        int *ptr = (int*)&rightVec_i;
        printf("new rightVec %d %d %d %d %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7]);
      }
#endif
      _mm256_storeu_si256((__m256i *) &table(out, j, i + 1), rightVec_i);
#if DEBUG == 1
      debug = false;
#endif
      // (t_{j-1,i+k}, ..., t_{j-1,i}) <-- vecT_{j-1,i}
      _mm256_storeu_si256((__m256i *) &table(out, j - 1, i), topVec_i);
      // (t_{j,i+k}, ..., t_{j,i}) <-- vecT_{j,i}
      _mm256_storeu_si256((__m256i *) &table(out, j, i), res_vec_i);
      // (t_{j+1,i+k}, ..., t_{j+1,i}) <-- vecT_{j+1,i}
      _mm256_storeu_si256((__m256i *) &table(out, j + 1, i), bottomVec_i);


      __m256i mask = _mm256_xor_si256(res_vec_i, vec_i);
      if (_mm256_testz_si256(mask, mask) != 1)
        diff = 1;
    }

  return diff;
}

#endif
#endif
#pragma endregion asandAVX

#pragma endregion asynchronousKernel

#ifdef ENABLE_MPI
#include <mpi.h>
#pragma region MPI

static int rank, size;

void ssandPile_init_mpi()
{
    easypap_check_mpi(); // check if MPI was correctly configured

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ssandPile_init();
}

static int rankTop(int rank)
{
  return rank*(DIM/size);
}

static int rankDown(int rank)
{
  return rankTop(rank+1) - DIM;
}

static int* leftBorder(int rank, int borderSize)
{
  int* border = (int*)malloc(sizeof(int)*borderSize);
  if(border==NULL)
  {
    fprintf(stderr, "malloc failure in leftBorder at rank %d!\n", rank);
    exit(EXIT_FAILURE);
  }
  
  for(int i=0; i<borderSize; i++)
  {
    border[i] = rankTop(rank)+i*DIM;
  }
  return border;
}

static int* rightBorder(int rank, int borderSize)
{
  int* border = (int*)malloc(sizeof(int)*borderSize);
  if(border==NULL)
  {
    fprintf(stderr, "malloc failure in rightBorder at rank %d!\n", rank);
    exit(EXIT_FAILURE);
  }
  
  for(int i=0; i<borderSize; i++)
  {
    border[i] = (rankTop(rank)+DIM-1)*(i+1);
  }
  return border;
}

//PERHAPS ???
static int rankSize(int rank)
{
  if(rank==size-1)
  {
    return DIM-rankTop(rank);
  }
  return DIM/size;
}

void ssandPile_refresh_img_mpi()
{
  MPI_Status status;
  //processus Ma??tre = reception
  if(rank==0)
  {
    for(int i=1; i<size; i++)
    {
      MPI_Recv(&cur_img(rankTop(i), 0), rankSize(i)*DIM, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
    }
  }
  //processus ouvrier
  else
  {
    // int borderSize = DIM/size;
    if(rank%2==0)
    {
      //processus pair --> send ; compute? ; receive
      int upperZone[DIM];
      int lowerZone[DIM];

      MPI_Send(&table(in, 0, rankTop(rank)), sizeof(int)*DIM, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
      if(rank!=size-1)
      {
        MPI_Send(&table(in, 0, rankDown(rank)), sizeof(int)*DIM, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
      }

      //compute ?

      if(rank!=size-1)
        MPI_Recv(&lowerZone, sizeof(int)*DIM, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&upperZone, sizeof(int)*DIM, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status);
    }
    else
    {
      //processus impair --> receive ; compute? ; send
      int upperZone[DIM];
      int lowerZone[DIM];

      if(rank!=size-1)
        MPI_Recv(&lowerZone, sizeof(int)*DIM, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &status);
      if(rank!=1)
        MPI_Recv(&upperZone, sizeof(int)*DIM, MPI_INT, rank+1, 0, MPI_COMM_WORLD, &status);

      //compute ?

      if(rank!=1)
        MPI_Send(&table(in, 0, rankTop(rank)), sizeof(int)*DIM, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
      if(rank!=size-1)
        MPI_Send(&table(in, 0, rankDown(rank)), sizeof(int)*DIM, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }
    MPI_Send(&table(in, 0, rankTop(rank)), rankSize(rank), MPI_INT, 0, 0, MPI_COMM_WORLD);
  }
}

int ssandPile_do_tile_mpi(int x, int y, int width, int height)
{
  //same as refresh_img_mpi ????
  return 0;
}

unsigned ssandPile_compute_mpi(unsigned nb_iter)
{
    for (unsigned it = 1; it <= nb_iter; it++)
    {
      do_tile(0, rankTop(rank), DIM, rankSize(rank), 0);
      //zoom();
    }
    return 0;
}
#endif