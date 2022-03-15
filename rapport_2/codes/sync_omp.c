unsigned ssandPile_compute_omp(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {

    int change = 0;
#pragma omp parallel for schedule(runtime) reduction(|: change)
      for (int y = 1; y < DIM-1; y += 1)
        for (int x = 1; x < DIM-1; x += 1)
          change |= do_tile(x, y, 1, 1, omp_get_thread_num());
      swap_tables();
      if (change == 0)
        return it;
  }

  return 0;
}