#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

int main(int argc, char* argv[])
{
    uint64_t cache, cache2;
    uint64_t n;
    uint64_t m;
    uint64_t nnz;

    FILE *fread;
    fread = fopen(argv[1], "r");
    printf("Codefile :: %s\n", argv[1]);

    FILE *fp;
    fp = fopen("codefile_mindist.txt", "w");


    fscanf(fread, "nc: %lu\n", &n);
    fscanf(fread, "mc: %lu\n", &m);
    fscanf(fread, "nct: %lu\n", &cache);
    fscanf(fread, "mct: %lu\n", &cache);
    fscanf(fread, "nnz: %lu\n", &nnz);
    fscanf(fread, "puncture [%lu]: ", &cache);
    fscanf(fread, "shorten [%lu]: ", &cache);

    fprintf(fp, "%lu\n", n);
    fprintf(fp, "%lu\n", m);

    size_t* column_weight = calloc(m, sizeof(size_t));
    size_t** column = calloc(m, sizeof(size_t*));
    for (size_t i=0; i<m; ++i)
        column[i] = calloc(0, sizeof(size_t));

    size_t tmp = 0;
    for (size_t i=0; i<nnz; ++i)
    {
        fscanf(fread, "%lu %lu\n", &cache, &cache2);
        ++column_weight[cache];
        column[cache] = realloc(column[cache], column_weight[cache] * sizeof(size_t));
        column[cache][column_weight[cache]-1] = cache2;

        if (column_weight[cache] > tmp)
            tmp = column_weight[cache];
    }

    fprintf(fp, "%lu\n", tmp);

    for (size_t j=0; j<m; ++j)
    {
        for (size_t i=0; i<tmp; ++i)
        {
            if (i < column_weight[j])
                fprintf(fp, "%lu ", column[j][i]);
            else
                fprintf(fp, "0 ");
        }
        fprintf(fp, "\n");
    }

    free(column_weight);
    for (size_t i=0; i<m; ++i)
        free(column[i]);
    free(column);

    fclose(fread);
    fclose(fp);

    return 0;
}
