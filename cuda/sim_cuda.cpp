#include "ldpc/ldpc.h"

int main()
{
    LDPC* sim = new LDPC("code12x24.txt");

    sim->print_ldpc_code();

    delete sim;

    return 0;
}
