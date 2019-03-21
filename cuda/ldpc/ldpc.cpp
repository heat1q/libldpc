#include "ldpc.h"
#include <exception>

using namespace std;
using namespace ldpc;

Ldpc_Code_cl::Ldpc_Code_cl(const char* filename)
{
    //init
    puncture_c = nullptr;
    shorten_c = nullptr;
    cw_c = nullptr;
    cn_c = nullptr;
    vw_c = nullptr;
    vn_c = nullptr;
    r_c = nullptr;
    c_c = nullptr;

    try
    {
        FILE *fp;

        fp = fopen(filename, "r");
        if(!fp)
            throw runtime_error("can not open codefile for reading.");

        fscanf(fp, "nc: %lu\n", &n_c);
        fscanf(fp, "mc: %lu\n", &m_c);
        fscanf(fp, "nct: %lu\n", &nct_c);
        fscanf(fp, "mct: %lu\n", &mct_c);
        fscanf(fp,  "nnz: %lu\n", &nnz_c);
        k_c = n_c-m_c;
        kct_c = nct_c-mct_c;

        fscanf(fp, "puncture [%lu]: ", &(num_puncture_c));
        num_puncture_sys_c = 0;
        num_puncture_par_c = 0;
        if(num_puncture_c != 0)
        {
            puncture_c = new size_t[num_puncture_c];
            for(size_t i = 0; i < num_puncture_c; i++)
            {
                fscanf(fp, " %lu ", &(puncture_c[i]));
                if(puncture_c[i] < k_c)
                    num_puncture_sys_c++;
                else
                    num_puncture_par_c++;
            }
        }

        fscanf(fp, "shorten [%lu]: ", &num_shorten_c);
        if(num_shorten_c != 0)
        {
            shorten_c = new size_t[num_shorten_c];
            for(size_t i = 0; i < num_shorten_c; i++)
                fscanf(fp, " %lu ", &(shorten_c[i]));
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

        max_dc_c = 0;
        for(size_t i = 0; i < m_c; i++)
        {
            if(cw_c[i] > max_dc_c)
                max_dc_c = cw_c[i];
        }

        fclose(fp);
    }
    catch(exception &e)
    {
        cout << "Error: " << e.what() << endl;
        destroy_ldpc_code();

        exit(EXIT_FAILURE);
    }
}

Ldpc_Code_cl::~Ldpc_Code_cl() { destroy_ldpc_code(); }

#ifdef QC_LYR_DEC
Ldpc_Code_cl::Ldpc_Code_cl(const char* filename, const char* clfile)
{
    //init
    puncture_c = nullptr;
    shorten_c = nullptr;
    cw_c = nullptr;
    cn_c = nullptr;
    vw_c = nullptr;
    vn_c = nullptr;
    r_c = nullptr;
    c_c = nullptr;
    #ifdef QC_LYR_DEC
    lw_c = nullptr;
    layers_c = nullptr;
    #endif

    try
    {
        FILE *fp;

        fp = fopen(filename, "r");
        if(!fp)
            throw runtime_error("can not open codefile for reading.");

        fscanf(fp, "nc: %lu\n", &n_c);
        fscanf(fp, "mc: %lu\n", &m_c);
        fscanf(fp, "nct: %lu\n", &nct_c);
        fscanf(fp, "mct: %lu\n", &mct_c);
        fscanf(fp,  "nnz: %lu\n", &nnz_c);
        k_c = n_c-m_c;
        kct_c = nct_c-mct_c;

        fscanf(fp, "puncture [%lu]: ", &(num_puncture_c));
        num_puncture_sys_c = 0;
        num_puncture_par_c = 0;
        if(num_puncture_c != 0)
        {
            puncture_c = new size_t[num_puncture_c];
            for(size_t i = 0; i < num_puncture_c; i++)
            {
                fscanf(fp, " %lu ", &(puncture_c[i]));
                if(puncture_c[i] < k_c)
                    num_puncture_sys_c++;
                else
                    num_puncture_par_c++;
            }
        }

        fscanf(fp, "shorten [%lu]: ", &num_shorten_c);
        if(num_shorten_c != 0)
        {
            shorten_c = new size_t[num_shorten_c];
            for(size_t i = 0; i < num_shorten_c; i++)
                fscanf(fp, " %lu ", &(shorten_c[i]));
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

        max_dc_c = 0;
        for(size_t i = 0; i < m_c; i++)
        {
            if(cw_c[i] > max_dc_c)
                max_dc_c = cw_c[i];
        }

        fclose(fp);

        //layers
        setup_layers(clfile);
    }
    catch(exception &e)
    {
        cout << "Error: " << e.what() << endl;
        destroy_ldpc_code();

        exit(EXIT_FAILURE);
    }
}

void Ldpc_Code_cl::setup_layers(const char *clfile)
{
    FILE *fp = fopen(clfile, "r");
    if(!fp)
        throw runtime_error("Can not open layer file");

    fscanf(fp, "nl: %lu\n", &nl_c);

    lw_c = new uint64_t[nl_c];
    layers_c = new uint64_t*[nl_c];

    for (size_t i = 0; i < nl_c; ++i)
    {
        fscanf(fp, "cn[i]: %lu\n", &(lw_c[i]));
        layers_c[i] = new uint64_t[lw_c[i]];
        for (size_t j = 0; j < lw_c[i]; ++j)
            fscanf(fp, "%lu\n", &(layers_c[i][j]));
    }

    printf("=========== LDPC LAYERS ===========\n");
    printf("nl: %lu\n", nl_c);
    for (size_t i = 0; i < nl_c; ++i)
    {
        printf("cn[%lu]: %lu\n", i, lw_c[i]);
        printVector<uint64_t>(layers_c[i], lw_c[i]);
        printf("\n");
    }
    printf("========= LDPC LAYERS: END ========\n");

    fclose(fp);
}
#endif

void Ldpc_Code_cl::destroy_ldpc_code()
{
    if (vn_c != nullptr)
    {
        for(size_t i = 0; i < n_c; i++)
            delete[] vn_c[i];
        delete[] vn_c;
    }

    if (vn_c != nullptr)
    {
        for(size_t i = 0; i < m_c; i++)
            delete[] cn_c[i];
        delete[] cn_c;
    }

    if (vw_c != nullptr)
        delete[] vw_c;

    if (cw_c != nullptr)
        delete[] cw_c;

    if (r_c != nullptr)
        delete[] r_c;

    if (c_c != nullptr)
        delete[] c_c;

    if (puncture_c != nullptr)
        delete[] puncture_c;

    if (shorten_c != nullptr)
        delete[] shorten_c;

    #ifdef QC_LYR_DEC
    if (layers_c != nullptr)
    {
        for(size_t i = 0; i < nl_c; i++)
            delete[] layers_c[i];
        delete[] layers_c;
    }
    if (lw_c != nullptr)
        delete[] lw_c;
    #endif
}


void Ldpc_Code_cl::print_ldpc_code()
{
    cout << "=========== LDPC ===========" << endl;
    cout << "nc : " << n_c << endl;
    cout << "mc : " << m_c << endl;
    cout << "kc : " << k_c << endl;
    cout << "nnz : " << nnz_c << endl;
    cout << "nct :" << nct_c << endl;
    cout << "mct : " << mct_c << endl;
    cout << "kct : " << kct_c << endl;
    cout << "max dc : " << max_dc_c << endl;
    cout << "num puncture: " << num_puncture_c << endl;
    cout << "num puncture sys: " << num_puncture_sys_c << endl;
    cout << "num puncture par: " << num_puncture_par_c << endl;
    cout << "num shorten: " << num_shorten_c << endl;
    cout << "=========== LDPC: END ===========" << endl;
}

uint64_t Ldpc_Code_cl::nc() const { return n_c; }
uint64_t Ldpc_Code_cl::kc() const { return k_c; }
uint64_t Ldpc_Code_cl::mc() const { return m_c; }
uint64_t Ldpc_Code_cl::nnz() const { return nnz_c; }
size_t *Ldpc_Code_cl::cw() const { return cw_c; }
size_t *Ldpc_Code_cl::vw() const { return vw_c; }
size_t **Ldpc_Code_cl::cn() const { return cn_c; }
size_t **Ldpc_Code_cl::vn() const { return vn_c; }
size_t *Ldpc_Code_cl::r() const { return r_c; }
size_t *Ldpc_Code_cl::c() const { return c_c; }
uint64_t Ldpc_Code_cl::nct() const { return nct_c; }
uint64_t Ldpc_Code_cl::mct() const { return mct_c; }
size_t *Ldpc_Code_cl::puncture() const { return puncture_c; }
size_t Ldpc_Code_cl::num_puncture() const { return num_puncture_c; }
size_t *Ldpc_Code_cl::shorten() const { return shorten_c; }
size_t Ldpc_Code_cl::num_shorten() const { return num_shorten_c; }
uint64_t Ldpc_Code_cl::kct() const { return kct_c; }
size_t Ldpc_Code_cl::max_dc() const { return max_dc_c; }

#ifdef QC_LYR_DEC
uint64_t Ldpc_Code_cl::nl() const { return nl_c; }
uint64_t *Ldpc_Code_cl::lw() const { return lw_c; }
uint64_t **Ldpc_Code_cl::layers() const { return layers_c; }
#endif


void ldpc::dec2bin(uint64_t val, uint8_t m)
{
    for(size_t i = 0; i < m; i++)
        printf("%lu", (val>>(m-i-1) & 0x01));
}

template<typename T> void ldpc::printVector(T *x, const size_t &l)
{
    cout << "[";
    for (size_t i = 0; i < l-1; ++i)
        cout << x[i] << " ";
    cout << x[l-1] << "]";
}
