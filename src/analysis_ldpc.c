#include <unistd.h>
#include <string.h>

#include "function/ldpc_types.h"
#include "function/ldpc_functions.h"
#include "LDPC_SSET/ldpc_cycles.h"
#include "decoder/ldpc_decoder.h"
#include "LDPC_SSET/ldpc_stoppingsets.h"

int main(int argc, char* argv[])
{
    uint8_t abort = 0;

    char codeName[128];
    char stFile[128];
    char stcountFile[128];
    uint64_t maxSize = 0;
    uint64_t ImaxBP = 10;
    int ImaxE = 10;
    double InitLLR = -10.0;

    const int numArgs = (argc-1)/2;

    if (argc == 7)
        for (int i = 0; i < numArgs; ++i)
        {
            if (strcmp(argv[2*i+1], "-code") == 0)
                strcpy(codeName, argv[2*i+2]);
            else if (strcmp(argv[2*i+1], "-out") == 0)
            {
                strcpy(stFile, argv[2*i+2]);
                char* tmpstr = strtok(argv[2*i+2], ".");
                strcpy(stcountFile, tmpstr);
                strcat(stcountFile, "_count.txt");
            }
            else if (strcmp(argv[2*i+1], "-thresh") == 0)
                maxSize = (uint64_t) atoi(argv[2*i+2]);
            else
                abort = 1;
        }
    else if (argc == 9 || argc == 11 || argc == 13)
    {
        for (int i = 0; i < numArgs; ++i)
        {
            if (strcmp(argv[2*i+1], "-code") == 0)
                strcpy(codeName, argv[2*i+2]);
            else if (strcmp(argv[2*i+1], "-out") == 0)
            {
                strcpy(stFile, argv[2*i+2]);
                char* tmpstr = strtok(argv[2*i+2], ".");
                strcpy(stcountFile, tmpstr);
                strcat(stcountFile, "_count.txt");
            }
            else if (strcmp(argv[2*i+1], "-thresh") == 0)
                maxSize = (uint64_t) atoi(argv[2*i+2]);
            else if (strcmp(argv[2*i+1], "-iterBP") == 0)
                ImaxBP = (uint64_t) atoi(argv[2*i+2]);
            else if (strcmp(argv[2*i+1], "-iterE") == 0)
                ImaxE = (int) atoi(argv[2*i+2]);
            else if (strcmp(argv[2*i+1], "-llr") == 0)
                InitLLR = (double) atof(argv[2*i+2]);
            else
                abort = 1;
        }
    }
    else
        abort = 1;

    if (abort)
    {
        printf("====================== LDPC Stopping Sets ======================\n");
        printf("                         Usage Reminder:                        \n");
        printf("         Main -code CodeName -out OutName -thresh maxSize       \n");
        printf("                 optional: -iterBP Imax                         \n");
        printf("                           -iterE Imax                          \n");
        printf("                           -llr LLR                             \n");
        printf("                                                                \n");
        printf("                 CodeName: Name of the code file                \n");
        printf("       OutName: Name of the result file, eg. results.txt        \n");
        printf("       maxSize: Maximum size of Stopping Sets to be count       \n");
        printf("       Imax: [int] Maximum iterations of the belief propagation \n");
        printf("             and erasure decoder (default: Imax=10)             \n");
        printf("       LLR: [double] init LLRs for BP algorithm (see ReadMe.txt)\n");
        printf("            (LLR < 0.0) default: LLR=-10.0                      \n");
        printf("================================================================\n");
        exit(EXIT_FAILURE);
    }

    ldpc_code_t *code = calloc(1, sizeof(ldpc_code_t));

    read_ldpc_file(code, codeName);

    if (!maxSize || maxSize > code->nc)
        maxSize = code->nc;

    if (InitLLR >= 0.0)
        InitLLR = -10.0;

    lpdc_code_t_stopping_sets(code, stFile, stcountFile, maxSize, ImaxBP, ImaxE, InitLLR);

    destroy_ldpc_code_t(code);
    free(code);

    return 0;
}
