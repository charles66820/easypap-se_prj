#include "easypap.h"

#include <omp.h>
#include <stdbool.h>
#include <sys/mman.h>
#include <unistd.h>

typedef unsigned int TYPE;

static TYPE *TABLE = NULL;

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

static int in = 0;
static int out = 1;

static inline void swap_tables()
{
  int tmp = in;
  in = out;
  out = tmp;
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
        r = b = 255 - (240 * ((double)g) / (double)max_grains);

      cur_img(i, j) = RGB(r, v, b);
      if (g > max)
        max = g;
    }
  max_grains = max;
}

/////////////////////////////  Initial Configurations

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
      atable(i, j) = 4;
}

void asandPile_draw_DIM(void)
{
  max_grains = DIM;
  for (int i = DIM / 4; i < DIM - 1; i += DIM / 4)
    for (int j = DIM / 4; j < DIM - 1; j += DIM / 4)
      atable(i, j) = i * j / 4;
}

void asandPile_draw_alea(void)
{
  max_grains = 5000;
  for (int i = 0; i < DIM >> 3; i++)
  {
    atable(1 + random() % (DIM - 2), 1 + random() % (DIM - 2)) =
        1000 + (random() % (4000));
  }
}

void asandPile_draw_big(void)
{
  const int i = DIM / 2;
  atable(i, i) = 100000;
}

static void one_spiral(int x, int y, int step, int turns)
{
  int i = x, j = y, t;

  for (t = 1; t <= turns; t++)
  {
    for (; i < x + t * step; i++)
      atable(i, j) = 3;
    for (; j < y + t * step + 1; j++)
      atable(i, j) = 3;
    for (; i > x - t * step - 1; i--)
      atable(i, j) = 3;
    for (; j > y - t * step - 1; j--)
      atable(i, j) = 3;
  }
  atable(i, j) = 4;

  for (int i = -2; i < 3; i++)
    for (int j = -2; j < 3; j++)
      atable(i + x, j + y) = 3;
}

static void many_spirals(int xdebut, int xfin, int ydebut, int yfin, int step,
                         int turns)
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

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Synchronous Kernel
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

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
      if (table(out, i, j) >= 4)
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
      if (result >= 4)
        diff = 1;
    }

  return diff;
}

/////////////////////////////  Sequential version (seq)
// Suggested cmdline(s):
// ./run -k ssandPile -v seq -s 512 -m
//
// ./run -k ssandPile -v seq -wt opt -s 512 -m
//
// Renvoie le nombre d'itérations effectuées avant stabilisation, ou 0
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

/////////////////////////////  Sequential version (omp)
// Suggested cmdline(s):
// ./run -k ssandPile -v omp -s 512 -m
//
// ./run -k ssandPile -v omp -wt opt -s 512 -m
//
// Renvoie le nombre d'itérations effectuées avant stabilisation, ou 0
unsigned ssandPile_compute_omp(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {

    int change = 0;
#pragma omp parallel for schedule(runtime) shared(change)
      for (int y = 1; y < DIM-1; y += 1)
        for (int x = 1; x < DIM-1; x += 1)
        {
          int localChange = do_tile(x, y, 1, 1, omp_get_thread_num());
          if (change == 0 && localChange != 0)
          {
  #pragma omp critical
            change |= localChange;
          }
        }
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
        change |=
            do_tile(x + (x == 0), y + (y == 0),
                    TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                    TILE_H - ((y + TILE_H == DIM) + (y == 0)), 0 /* CPU id */);
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

#pragma omp parallel for schedule(runtime) collapse(2) shared(change)
    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
      {
        int localChange =
            do_tile(x + (x == 0), y + (y == 0),
                    TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                    TILE_H - ((y + TILE_H == DIM) + (y == 0)), omp_get_thread_num());
        if (change == 0 && localChange != 0)
        {
#pragma omp critical
          change |= localChange;
        }
      }
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
        int localChange =
            do_tile(x + (x == 0), y + (y == 0),
                    TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                    TILE_H - ((y + TILE_H == DIM) + (y == 0)), omp_get_thread_num());
        if (change == 0 && localChange != 0)
        {
#pragma omp critical
          change |= localChange;
        }
      }
    swap_tables();
    if (change == 0)
      return it;
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
///////////////////////////// Asynchronous Kernel
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void asandPile_init()
{
  in = out = 0;
  if (TABLE == NULL)
  {
    const unsigned size = DIM * DIM * sizeof(TYPE);

    PRINT_DEBUG('u', "Memory footprint = 2 x %d bytes\n", size);

    TABLE = mmap(NULL, size, PROT_READ | PROT_WRITE,
                 MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
  }
}

void asandPile_finalize()
{
  const unsigned size = DIM * DIM * sizeof(TYPE);

  munmap(TABLE, size);
}

///////////////////////////// Version séquentielle simple (seq)
// Renvoie le nombre d'itérations effectuées avant stabilisation, ou 0

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
        result/=4;
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
        change |=
            do_tile(x + (x == 0), y + (y == 0),
                    TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                    TILE_H - ((y + TILE_H == DIM) + (y == 0)), 0 /* CPU id */);
    if (change == 0)
      return it;
  }

  return 0;
}

/////////////////////////////  Tiled sequential version (tiled)
// Suggested cmdline(s):
// ./run -k asandPile -v omp -s 512 -m
//
// ./run -k asandPile -v omp -wt opt -s 512 -m
//
unsigned asandPile_compute_omp(unsigned nb_iter)
{
    for (unsigned it = 1; it <= nb_iter; it++)
    {
      int change = 0;

      for(int i=0; i<2; i++)
      {
#pragma omp parallel for schedule(runtime) shared(change)
        for (int y = 0; y < DIM; y += TILE_H)
          for (int x = ((y - (TILE_H *i))%(TILE_H*2)); x < DIM; x += TILE_W*2) //==TILE_H*i) ? TILE_W : 0
          {
            int localChange =
                do_tile(x + (x == 0), y + (y == 0),
                        TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                        TILE_H - ((y + TILE_H == DIM) + (y == 0)), omp_get_thread_num());
            if (change == 0 && localChange != 0)
            {
#pragma omp critical
              change |= localChange;
            }
          }
        if (change == 0)
          return it;
      }
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
  unsigned res = 0;
#pragma omp parallel
#pragma omp master
  {
    for (unsigned it = 1; it <= nb_iter; it++) {
      int change = 0;
      int tuile[NB_TILES_Y][NB_TILES_X + 1] __attribute__ ((unused));

      for (int y = 0; y < DIM; y += TILE_H)
        for (int x = 0; x < DIM; x += TILE_W)
        {
#pragma omp task depend(in:tuile[x/TILE_W][(y/TILE_H)-1]) depend(in:tuile[(x/TILE_W)-1][y/TILE_H]) depend(out:tuile[x/TILE_W][y/TILE_H])
          change |= do_tile(x + (x == 0), y + (y == 0),
                            TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                            TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                            omp_get_thread_num());
        }
#pragma omp taskwait
      if (!change) {
        res = it;
        break;
      }
    }
  }
  return res;
}