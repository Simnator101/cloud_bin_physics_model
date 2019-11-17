#include "../include/mathfuncs.h"
#include "../include/extio.h"
#include "../include/environment.h"

void free_sounding(sounding* snd)
{
    free_vec(snd->pe);
    free_vec(snd->Tde);
    free_vec(snd->Te);
    free_vec(snd->ze);
    memset(snd, 0, sizeof(sounding));
}

sounding* __alloc_sounding(sounding* snd)
{
    snd->Te = allocate_vec(0);
    snd->Tde = allocate_vec(0);
    snd->pe = allocate_vec(0);
    snd->ze = allocate_vec(0);
    return snd;
}

vec* __snd_data_ref(sounding* snd, unsigned char code)
{
    switch (code)
    {
    case SOUNDING_T:
    case SOUNDING_TH:
        return snd->Te;
    case SOUNDING_TD:
        return snd->Tde;
    case SOUNDING_P:
        return snd->pe;
    case SOUNDING_Z:
        return snd->ze;
    case SOUNDING_Q:
        return snd->Tde;
    default:
        perror("Error Unidentified fmt code passed to \'read_sounding_txt(...)\'");
        break;
    }
    return NULL;
}

int read_sounding_txt(sounding* snd, const char* file_nm, unsigned long long fmt)
{
    // Decode fmt tag
    unsigned int nentries = 0;
    while (((fmt >> (8 * nentries)) & 0xff) != SOUNDING_END) ++nentries;

    // Clear Old Data and make pointers available
    free_sounding(snd);
    __alloc_sounding(snd);

    FILE* pf = fopen(file_nm, "r");
    if (pf == NULL)
    {
        char __err[128];
        sprintf(__err, "Failed to open file \"%s\"\n", file_nm);
        perror(__err);
        return 1;
    }


    unsigned long i;
    double value;
    long fsz = file_size(pf);
    char* __buff = calloc(fsz, 1);
    while(fgets(__buff, fsz, pf) != NULL)
    {
        if (__buff[0] == '#') continue; // Ignore line
        strbuffer sbuf = str_split(__buff, ',');
        size_t nsbuf = MIN(strbuffer_len(sbuf), nentries);
        for (i = 0; i < nsbuf; ++i)
        {
            if (str_eval_func(sbuf[i], isalpha))
                value = NAN;
            else
                value = strtod(sbuf[i], NULL);

            vec* pvec = __snd_data_ref(snd, (fmt >> (8 * i)) & 0xff);
            push_back(pvec, value);
        }

        free_strbuffer(sbuf);
        memset(__buff, 0, fsz);
    }

    // Flags Checks for conversion
    int TTde_fnd = 0;
    int conv_Th = 0;
    int conv_q = 0;
    int conv_z = 1, conv_p = 1;
    for (i = 0; i < 8; ++i)
    {
        unsigned char flg = (unsigned char)((fmt >> (8 * i)) & 0xff);
        if (flg == SOUNDING_P) {conv_p = 0; scl_vec(snd->pe, 1e2); }
        else if (flg == SOUNDING_T) { ++TTde_fnd; conv_Th = 0; add_scalar_vec(snd->Te, 273.15); }
        else if (flg == SOUNDING_TH) {++TTde_fnd; conv_Th = 1;}
        else if (flg == SOUNDING_Z) {conv_z = 0; scl_vec(snd->ze, 1e3);}
        else if (flg == SOUNDING_TD) {conv_q = 0; ++TTde_fnd; add_scalar_vec(snd->Tde, 273.15); }
        else if (flg == SOUNDING_Q) {conv_q = 1; ++TTde_fnd; scl_vec(snd->Tde, 1e-3); }
    }

    // Check that vertical coordinate is defined
    // Check that we have temperature/humidity data
    assert(conv_z != 1 || conv_p != 1);
    assert(TTde_fnd >= 2);

    // Conversion
    if (conv_z)
    {
        free_vec(snd->ze);
        snd->ze = convert_p_to_z(snd->pe, snd->Te, MODEL_SETTINGS.z0);
    }
    if (conv_p)
    {
        free_vec(snd->pe);
        snd->pe = convert_z_to_p(snd->ze, snd->Te, MODEL_SETTINGS.p0);
    }
    if (conv_q)
    {
        for (i = 0; i < snd->Tde->N; ++i)
            snd->Tde->data[i] = dewpoint_temp(snd->pe->data[i], snd->Tde->data[i]);
    }
    if (conv_Th)
    {
        vec* Th_conv = copy_vec(snd->pe);
        for (i = 0; i < LEN(Th_conv); ++i) Th_conv->data[i] = pow(Th_conv->data[i] / snd->pe->data[0], RAIR / CPA);
        mlt_vec(snd->Te, 1, Th_conv);
        free_vec(Th_conv);
    }

    fclose(pf);
    free(__buff);
    return 0;
}

double sample_snd_qe(double z, sounding snd)
{
    double Td = sample_snd_Tde(z, snd);
    double p = interp(z, snd.ze, snd.pe, LINEAR);
    return sat_mixr_vapour(p, Td);
}