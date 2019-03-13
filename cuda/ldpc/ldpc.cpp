#include "ldpc.h"

LDPC::LDPC(const char *filename)
{
    read_file(filename);
}

LDPC::~LDPC()
{
    for(size_t i = 0; i < n_c; i++)
        free(vn_c[i]);
    for(size_t i = 0; i < m_c; i++)
        free(cn_c[i]);

    free(vn_c);
    free(cn_c);
    free(vw_c);
    free(cw_c);
    free(r_c);
    free(c_c);
    free(puncture);
    free(shorten);
}

bool LDPC::read_file(const char* filename)
{
    FILE *fp;

    fp = fopen(filename, "r");
    if(!fp) {
        printf("can not open codefile %s for reading.\n", filename);
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "nc: %lu\n", &n_c);
    fscanf(fp, "mc: %lu\n", &m_c);
    fscanf(fp, "nct: %lu\n", &nct);
    fscanf(fp, "mct: %lu\n", &mct);
    fscanf(fp,  "nnz: %lu\n", &nnz_c);
    k_c = n_c-m_c;
    kct = nct-mct;

    fscanf(fp, "puncture [%lu]: ", &(num_puncture));
    num_puncture_sys = 0;
    num_puncture_par = 0;
    if(num_puncture != 0)
    {
        puncture = new size_t[num_puncture];
        for(size_t i = 0; i < num_puncture; i++)
        {
            fscanf(fp, " %lu ", &(puncture[i]));
            if(puncture[i] < k_c)
                num_puncture_sys++;
            else
                num_puncture_par++;

        }
    }
    else
    {
        puncture = nullptr;
    }
    fscanf(fp, "shorten [%lu]: ", &num_shorten);
    if(num_shorten != 0)
    {
        shorten = new size_t[num_shorten];
        for(size_t i = 0; i < num_shorten; i++)
            fscanf(fp, " %lu ", &(shorten[i]));
    }
    else
    {
        shorten = nullptr;
    }

    size_t* cw_tmp;
    size_t* vw_tmp;
    cw_c = new uint64_t[m_c] ();
    cw_tmp = new uint64_t[m_c] ();
    vw_c = new uint64_t[n_c] ();
    vw_tmp = new uint64_t[n_c] ();
    r_c = new uint64_t[nnz_c] ();
    c_c = new uint64_t[nnz_c] ();


    for(size_t i = 0; i < nnz_c; i++)
    {
        fscanf(fp, "%lu %lu\n", &(r_c[i]), &(c_c[i]));
        cw_c[r_c[i]]++;
        vw_c[c_c[i]]++;
    }

    cn_c = new size_t*[m_c] ();
    for(size_t i = 0; i < m_c; i++)
        cn_c[i] = new size_t[cw_c[i]] ();

    vn_c = new size_t*[n_c] ();
    for(size_t i = 0; i < n_c; i++)
        vn_c[i] = new size_t[vw_c[i]] ();

    for(size_t i = 0; i < nnz_c; i++)
    {
        cn_c[r_c[i]][cw_tmp[r_c[i]]++] = i;
        vn_c[c_c[i]][vw_tmp[c_c[i]]++] = i;
    }

    delete[] cw_tmp;
    delete[] vw_tmp;

    max_dc = 0;
    for(size_t i = 0; i < m_c; i++)
    {
        if(cw_c[i] > max_dc)
            max_dc = cw_c[i];
    }

    fclose(fp);

    return true;
}

void LDPC::print_ldpc_code()
{
    cout << "=========== LDPC ===========" << endl;
    cout << "nc : " << n_c << endl;
    cout << "mc : " << m_c << endl;
    cout << "kc : " << k_c << endl;
    cout << "nnz : " << nnz_c << endl;
    cout << "nct :" << nct << endl;
    cout << "mct : " << mct << endl;
    cout << "kct : " << kct << endl;
    cout << "max dc : " << max_dc << endl;
    cout << "num puncture: " << num_puncture << endl;
    cout << "num puncture sys: " << num_puncture_sys << endl;
    cout << "num puncture par: " << num_puncture_par << endl;
    cout << "num shorten: " << num_shorten << endl;
    cout << "=========== LDPC: END ===========" << endl;
}

uint64_t LDPC::n() const
{
    return n_c;
}

uint64_t LDPC::k() const
{
    return k_c;
}

uint64_t LDPC::m() const
{
    return m_c;
}

uint64_t LDPC::nnz() const
{
    return nnz_c;
}

size_t *LDPC::cw() const
{
    return cw_c;
}

size_t *LDPC::vw() const
{
    return vw_c;
}

size_t **LDPC::cn() const
{
    return cn_c;
}

size_t **LDPC::vn() const
{
    return vn_c;
}

size_t *LDPC::r() const
{
    return r_c;
}

size_t *LDPC::c() const
{
    return c_c;
}


