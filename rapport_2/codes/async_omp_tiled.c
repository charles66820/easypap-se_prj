unsigned asandPile_compute_omp_tiled(unsigned nb_iter)
{
    for (unsigned it = 1; it <= nb_iter; it++)
    {
        int change = 0;

        //BLEU
        #pragma omp parallel for schedule(runtime) shared(change)  
        for(int y=0; y<DIM; y+=2*TILE_H)
        {
            for (int x = y%(TILE_H*2); x < DIM; x += TILE_W*2)
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
        }

        //ROUGE
        #pragma omp parallel for schedule(runtime) shared(change)  
        for(int y=0; y<DIM; y+=2*TILE_H)
        {
            for (int x = (y+TILE_H)%(TILE_H*2); x < DIM; x += TILE_W*2)
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
        }
        
        //VERT
        #pragma omp parallel for schedule(runtime) shared(change)  
        for(int y=TILE_H; y<DIM; y+=2*TILE_H)
        {
            for (int x = y%(TILE_H*2); x < DIM; x += TILE_W*2)
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
        }

        //NOIR
        #pragma omp parallel for schedule(runtime) shared(change)  
        for(int y=TILE_H; y<DIM; y+=2*TILE_H)
        {
            for (int x = (y+TILE_H)%(TILE_H*2); x < DIM; x += TILE_W*2)
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
        }

        if (change == 0)
            return it;
    }
  return 0;
}