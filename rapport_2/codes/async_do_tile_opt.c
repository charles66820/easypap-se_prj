int asandPile_do_tile_default(int x, int y, int width, int height)
{
  int change = 0;

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
-      if (atable(i, j) >= 4)
+      int result = atable(i, j);
+      if (result >= 4)
      {
+        result/=4;
-        atable(i, j - 1) += atable(i, j) / 4;
+        atable(i, j - 1) += result;
-        atable(i, j + 1) += atable(i, j) / 4;
+        atable(i, j + 1) += result;
-        atable(i - 1, j) += atable(i, j) / 4;
+        atable(i - 1, j) += result;
-        atable(i + 1, j) += atable(i, j) / 4;
+        atable(i + 1, j) += result;
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