#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "omp.h"
#define Nx 3840
#define Ny 2160
#define M 32
#define LEN 121
#define MAPSIZE 65536
char a[Ny][Nx] = {0};
char tmp[Ny][Nx] = {0};

void readPattern(char *filename)
{
  FILE *file = fopen(filename, "r");
  char line[256];
  int x = -1, y = -1;

  int row, col;
  for (int i = 0; fgets(line, 256, file) != NULL; ++i)
  {
    if (i == 0)
    {
      continue;
    }
    else if (i == 1)
    {
      sscanf(line, "x = %d, y = %d", &x, &y);
      row = (Ny - y) / 2;
      col = (Nx - x) / 2;
    }
    else
    {
      // printf("%s", line);
      int len = strlen(line);
      int run_count = 1, tag = 0, cur = 0;
      while (cur < len)
      {
        if ('0' <= line[cur] && line[cur] <= '9')
        {
          sscanf(line + cur, "%d", &run_count);
          // printf("run count %d\n", run_count);
          while ((cur < len) && ('0' <= line[cur] && line[cur] <= '9'))
          {
            // printf("cur %d %c \n", cur, line[cur]);
            cur++;
          }
        }
        else if (line[cur] == 'b' || line[cur] == 'o')
        {
          tag = line[cur];
          for (int k = 0; k < run_count; ++k)
          {
            if (tag == 'o')
              a[row][col + k] = 1;
          }
          col += run_count;
          run_count = 1;
          cur++;
        }
        else if (line[cur] == '$')
        {
          for (int k = 0; k < run_count; ++k)
          {
            row++;
          }
          col = (Nx - x) / 2;
          run_count = 1;
          cur++;
        }
        else if (line[cur] == '\n')
        {
          cur++;
        }
        else if (line[cur] == '!')
        {
          return;
        }
      }
    }
  }
}

void runConwayLifeGame(int max_iter, unsigned int **bita, unsigned int **bittmp, char *Table)
{
  // 上下更新限制
  int up = 1;
  int down = Ny + 1;
  for (int iter = 0; iter < max_iter; ++iter)
  {
    if (up > 1)
      up--;
    if (down < Ny + 1)
      down++;
    unsigned int kernel = 0;
#pragma omp parallel for schedule(dynamic) private(kernel) num_threads(48)
    for (int i = up; i < down; i += 2)
    {
      kernel = ((bita[i - 1][0] & 0xC0000000) >> 16) |
               ((bita[i][0] & 0xC0000000) >> 20) |
               ((bita[i + 1][0] & 0xC0000000) >> 24) |
               ((bita[i + 2][0] & 0xC0000000) >> 28);
      for (int j = 2; j < Nx + 2; j += 2)
      {
        char pos = j / M;
        char rest = j % M;
        kernel &= 0xCCCC;
        kernel |= (((bita[i - 1][pos] >> (30 - rest)) & 0x3)) << 12;
        kernel |= (((bita[i][pos] >> (30 - rest)) & 0x3)) << 8;
        kernel |= (((bita[i + 1][pos] >> (30 - rest)) & 0x3)) << 4;
        kernel |= (((bita[i + 2][pos] >> (30 - rest)) & 0x3));
        if (kernel == 0)
          continue;
        char cnt = Table[kernel];
        if (cnt == 0)
        {
          kernel <<= 2;
          continue;
        }
        if (rest)
        {
          bittmp[i][pos] ^= ((cnt & 0xC) >> 2) << (31 - rest);
          bittmp[i + 1][pos] ^= (cnt & 0x3) << (31 - rest);
        }
        else
        {
          bittmp[i][pos - 1] ^= (cnt & 0x8) >> 3;
          bittmp[i][pos] ^= (cnt & 0x4) << 29;
          bittmp[i + 1][pos - 1] ^= (cnt & 0x2) >> 1;
          bittmp[i + 1][pos] ^= (cnt & 1) << 31;
        }
        kernel <<= 2;
      }
    }
    for (int i = 1; i < Ny + 1; i++)
    {
      memcpy(bita[i], bittmp[i], sizeof(unsigned int) * LEN);
    }
    int flag = 1;
    for (int i = up; i < down; i++)
    {
      for (int j = 0; j < LEN; j++)
        if (bita[i][j] != 0)
        {
          flag = 0;
          break;
        }
      if (flag)
        up++;
      else
      {
        break;
      }
    }
    flag = 1;
    for (int i = down - 1; i > up; i--)
    {
      for (int j = 0; j < LEN; j++)
        if (bita[i][j] != 0)
        {
          flag = 0;
          break;
        }
      if (flag)
        down--;
      else
      {
        break;
      }
    }

#ifdef DEBUG
    char filename[10];
    sprintf(filename, "./output/iter%d", iter);
    FILE *f = fopen(filename, "w");
    for (int i = 0; i < Ny; ++i)
    {
      for (int j = 0; j < Nx; ++j)
      {
        fprintf(f, "%d ", tmp[i][j]);
      }
      fprintf(f, "\n");
    }
    fclose(f);
#endif
  }
#pragma omp parallel for
  for (int i = 1; i < Ny + 1; i++)
  {
    int j;
    unsigned int t;
    t = bita[i][0];
    for (int k = 30; k >= 0; k--)
    {
      a[i - 1][k] = t & 1;
      t >>= 1;
    }

    for (j = 1; j < LEN - 1; j++)
    {
      t = bita[i][j];
      for (int k = 31; k >= 0; k--)
      {
        a[i - 1][j * M + k - 1] = t & 1;
        t >>= 1;
      }
    }
    for (int k = 0; k < Nx - j * M + 1; k++)
    {
      a[i - 1][j * M + k - 1] = (bita[i][j] >> (31 - k)) & 1;
    }
  }
  FILE *f = fopen("test_output", "wb");
#ifdef BINARY_OUT
  fwrite(a, Nx * Ny, sizeof(int), f);
#else
  char *str = (char *)malloc((Nx * (Ny + 1) + 1) * sizeof(char));
  for (int i = 0; i < Ny; ++i)
  {
    int j;
    for (j = 0; j < Nx; ++j)
    {
      str[i * (Ny + 1) + j] = a[i][j] + '0';
    }
    str[i * (Ny + 1) + j] = '\n';
  }
  str[Nx * (Ny + 1)] = 0;
  fwrite(str, Nx * (Ny + 1), sizeof(char), f);
  free(str);
#endif
  fclose(f);
}
void initGame(unsigned int **bita, unsigned int **bittmp, char *Table)
{
// 初始化
#pragma omp parallel for
  for (int i = 0; i < Ny + 2; i++)
  {
    bita[i] = (unsigned int *)malloc(sizeof(unsigned int) * LEN);
    bittmp[i] = (unsigned int *)malloc(sizeof(unsigned int) * LEN);
    for (int j = 0; j < LEN; j++)
    {
      bita[i][j] = 0;
    }
  }
// a转化为bita
#pragma omp parallel for
  for (int i = 0; i < Ny; i++)
  {
    int j;
    for (int k = 0; k < M - 1; k++)
      bita[i + 1][0] |= a[i][k] << (30 - k);
    for (j = 1; j < LEN - 1; j++)
    {
      for (int k = 0; k < M; k++)
        bita[i + 1][j] |= a[i][j * M + k - 1] << (31 - k);
    }
    for (int k = 0; k < Nx - j * M + 1; k++)
      bita[i + 1][j] |= a[i][j * M + k - 1] << (31 - k);
  }
#pragma omp parallel for
  for (int i = 0; i < Ny + 2; i++)
  {
    memcpy(bittmp[i], bita[i], sizeof(unsigned int) * LEN);
  }

  // 计算右下角3*3矩阵改变量
  // for (int i = 0; i < 0x1000; i = ((i | 0x888) + 1) & 0x1777) // 0x1777 = 0111 0111 0111 0x888=1000 1000 1000
  // {
  //   char cnt = 0;
  //   unsigned int bitcnt = i;
  //   char flag = bitcnt & 0x20 ? 1 : 0;
  //   while (bitcnt)
  //   {
  //     cnt++;
  //     bitcnt &= (bitcnt - 1);
  //   }
  //   cnt = ((!(cnt ^ 3)) || (flag && (!(cnt ^ 4))));
  //   Table[i] = (cnt != flag);
  // }
  int Init[] = {7, 19, 21, 22, 32, 33, 34, 36, 48, 55, 67, 69, 70, 81, 82, 84, 96, 103, 115, 117, 118, 119, 259, 261, 262, 273, 274, 276, 288, 295, 307, 309, 310, 311, 321, 322, 324, 336, 355, 357, 358, 359, 369, 370, 371, 372, 373, 374, 375, 515, 517, 518, 529, 530, 532, 544, 551, 563, 565, 566, 567, 577, 578, 580, 592, 611, 613, 614, 615, 625, 626, 627, 628, 629, 630, 631, 769, 770, 772, 784, 803, 805, 806, 807, 817, 818, 819, 820, 821, 822, 823, 832, 865, 866, 867, 868, 869, 870, 871, 880, 881, 882, 883, 884, 885, 886, 887, 1027, 1029, 1030, 1041, 1042, 1044, 1056, 1063, 1075, 1077, 1078, 1079, 1089, 1090, 1092, 1104, 1123, 1125, 1126, 1127, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1281, 1282, 1284, 1296, 1315, 1317, 1318, 1319, 1329, 1330, 1331, 1332, 1333, 1334, 1335, 1344, 1377, 1378, 1379, 1380, 1381, 1382, 1383, 1392, 1393, 1394, 1395, 1396, 1397, 1398, 1399, 1537, 1538, 1540, 1552, 1571, 1573, 1574, 1575, 1585, 1586, 1587, 1588, 1589, 1590, 1591, 1600, 1633, 1634, 1635, 1636, 1637, 1638, 1639, 1648, 1649, 1650, 1651, 1652, 1653, 1654, 1655, 1792, 1825, 1826, 1827, 1828, 1829, 1830, 1831, 1840, 1841, 1842, 1843, 1844, 1845, 1846, 1847, 1888, 1889, 1890, 1891, 1892, 1893, 1894, 1895, 1904, 1905, 1906, 1907, 1908, 1909, 1910, 1911};
#pragma omp parallel for
  for (int i = 0; i < 228; i++)
  {
    Table[Init[i]] = 1;
  }
// 计算4*4矩阵改变量
#pragma omp parallel for
  for (int i = 0; i < MAPSIZE; i++)
  {
    Table[i] = (1 & Table[i & 0x777]) +
               ((1 & Table[(i >> 1) & 0x777]) << 1) +
               ((1 & Table[(i >> 4) & 0x777]) << 2) +
               ((1 & Table[(i >> 5) & 0x777]) << 3); // ↖↗↙↘
  }
}
int main(int argc, char **argv)
{
  // for (int i = 0; i < N; ++i) {
  // for (int j = 0; j < N; ++j) {
  // if ((i - N/2) * (i - N/2) + (j - N/2) * (j - N/2) <= N/4 * N/4 && rand()%10 >4)
  // a[i][j] = 1;
  // }
  // }
  readPattern("test_pattern");

  unsigned int *bita[Ny + 2];
  unsigned int *bittmp[Ny + 2];
  char Table[MAPSIZE];
  initGame(bita, bittmp, Table);
  struct timeval start;
  struct timeval end;
  gettimeofday(&start, NULL);

  runConwayLifeGame(1000, bita, bittmp, Table);

  gettimeofday(&end, NULL);
  printf("%lf ms\n", ((end.tv_sec - start.tv_sec) * 1000 + (end.tv_usec - start.tv_usec) / 1000.0));
}