int asandPile_do_tile_avx(int x, int y, int width, int height)
{
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

      // vecD <-- vec_i / 4
      __m256i vecD = _mm256_srli_epi32(vec_i, 2);

      // (vecD << 1)
      __m256i vecDShiftLeft = _mm256_alignr_epi32(vec0_i, vecD, 1);

      // (vecD >> 1)
      __m256i vecDShiftRight = _mm256_alignr_epi32(vecD, vec0_i, 7);

      // vec_i <-- vec_i % 4 + vecDShiftLeft + vecDShiftRight
      __m256i res_vec_i = _mm256_add_epi32(_mm256_and_si256(vec_i, vec3_i),
                               _mm256_add_epi32(vecDShiftLeft, vecDShiftRight));

      // topVec_i <-- topVec_i + vecD
      topVec_i = _mm256_add_epi32(topVec_i, vecD);

      // bottomVec_i <-- bottomVec_i + vecD
      bottomVec_i = _mm256_add_epi32(bottomVec_i, vecD);

      // t_{j,i-1} <-- t_{j,i-1} + vecD[0]
      __m256i leftVec_i = _mm256_loadu_si256((__m256i *) &table(in, j, i - 1));

      leftVec_i         = _mm256_add_epi32(leftVec_i, vecD);

      _mm256_storeu_si256((__m256i *) &table(out, j, i - 1), leftVec_i);

      // t_{j,i+k+1} <-- t_{j,i+k+1} + vecD[k]  :  k = AVX_VEC_SIZE_INT-1
      __m256i rightVec_i = _mm256_loadu_si256((__m256i *) &table(in, j, i + 1));

      rightVec_i         = _mm256_add_epi32(rightVec_i, vecD);

      _mm256_storeu_si256((__m256i *) &table(out, j, i + 1), rightVec_i);

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