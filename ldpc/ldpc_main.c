#include "unistd.h"
#include "string.h"

#include "ldpc_types.h"
#include "functions.h"
#include "ldpc_cycles.h"
#include "ldpc_decoder.h"
#include "ldpc_stoppingsets.h"

int main(int argc, char* argv[])
{
    uint8_t abort = 0;

    char codeName[100];
    char stFile[100];
    char stcountFile[100];
    uint64_t maxSize = 0;

    const int numArgs = (argc-1)/2;

    if (argc == 5)
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
            else
                abort = 1;
        }
    else if (argc == 7)
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
        printf("                Main -code CodeName -out OutName                \n");
        printf("                 Option -thresh maxSize                         \n");
        printf("                                                                \n");
        printf("                 CodeName: Name of the code file                \n");
        printf("       OutName: Name of the result file, eg. results.txt        \n");
        printf("       maxSize: Maximum size of Stopping Sets to be count       \n");
        printf("================================================================\n");
        exit(EXIT_FAILURE);
    }

    ldpc_code_t *code = calloc(1, sizeof(ldpc_code_t));

    read_ldpc_file(code, codeName);

    if (!maxSize)
        maxSize = code->nc;

    lpdc_code_t_stopping_sets(code, stFile, stcountFile, maxSize);

    destroy_ldpc_code_t(code);
    free(code);

    return 0;
}
