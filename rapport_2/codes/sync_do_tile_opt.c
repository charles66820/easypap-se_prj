int ssandPile_do_tile_opt(int x, int y, int width, int height)
{
  int diff = 0;

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
    {
-      table(out, i, j) = table(in, i, j) % 4;
+      int result = table(in, i, j) % 4;
-      table(out, i, j) += table(in, i + 1, j) / 4;
+      result += table(in, i + 1, j) / 4;
-      table(out, i, j) += table(in, i - 1, j) / 4;
+      result += table(in, i - 1, j) / 4;
-      table(out, i, j) += table(in, i, j + 1) / 4;
+      result += table(in, i, j + 1) / 4;
-      table(out, i, j) += table(in, i, j - 1) / 4;
+      result += table(in, i, j - 1) / 4;
-      if (table(out, i, j) >= 4)
+      table(out, i, j) = result;
+      if (result >= 4)
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