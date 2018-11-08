#include "unistd.h"

#include "ldpc_types.h"
#include "functions.h"
#include "ldpc_cycles.h"

int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        printf("Code file needed!\n");
        exit(EXIT_FAILURE);
    }

    char fname[100];
    sprintf(fname, "code_%s.txt", argv[1]);

    if (access(fname, F_OK) == -1)
    {
        printf("Cannot find code file %s!\n", fname);
        exit(EXIT_FAILURE);
    }

    ldpc_code_t *code = calloc(1, sizeof(ldpc_code_t));

    read_ldpc_file(code, fname);
	
    /*
    print_ldpc_code_t(code);

    printf("Enter girth: (Enter 0 to calculate girth) ");
    scanf("%lu", &code.girth);

    if (code.girth == 0)
        girth_ldpc_code_t(&code);
    */

	girth_ldpc_code_t(code);

	/* fill in the new code here */
	cycle_ldpc_code_t(code, argv[1]);

    destroy_ldpc_code_t(code);
	free(code);
	
    return 0;
}
