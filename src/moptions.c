#include "../include/extio.h"

#define STR_CMP(X,Y) (strcmp(X,Y) == 0)

// Default settings
__model_settings MODEL_SETTINGS =
{
    // Run time
    720UL,
    0.5, 10,

    // Dimensionality
    0.0, 3e3,
    0.0, 9e3,
    50., 50.,

    // Cloud Bin Settings
    69, 0,
    0.25, 0.055, 

    // Bottoms
    1e5, 0.0,

    // LH
    FALSE, -5.0,

    // Stream Function
    270.0, 200.0, 1e3, 1e3,

    // CCN
    1.2e8, 0.4,

    "./data/shallow_cumulus_profile.txt",
    0x100104UL,
    //"./data/stratocumulus_profile.txt",
    //0x102008UL,

    "./data/default_out.nc",
    2,
    NC_OUTPUT_CLOUD_LIQUID,

    NULL
};

char* model_settings_add_str(char* str)
{
    if (str == NULL) return NULL;
    unsigned long nlen = 2;
    if (MODEL_SETTINGS.__str_int != NULL)
        nlen += strbuffer_len(MODEL_SETTINGS.__str_int);

    MODEL_SETTINGS.__str_int = realloc(MODEL_SETTINGS.__str_int, sizeof(char*) * nlen);
    MODEL_SETTINGS.__str_int[nlen - 1] = NULL;
    MODEL_SETTINGS.__str_int[nlen - 2] = calloc(strlen(str) + 1, 1);
    strcpy(MODEL_SETTINGS.__str_int[nlen - 2], str);
    return MODEL_SETTINGS.__str_int[nlen - 2];
}

int read_job_settings(const char *file_name)
{
    FILE *pf = fopen(file_name, "r");
    long fsz = file_size(pf);
    if (fsz <= 0) return -1;

    char* file_data = calloc(fsz + 1, 1);
    fread(file_data, 1, fsz, pf);
    fclose(pf);

    // Split Into New strings
    unsigned long i;
    long comment_start;
    long cmd_start, cmd_end;
    strbuffer buf =  str_split(file_data, '\n');
    for (i = 0; i < strbuffer_len(buf);  ++i)
    {
        if (strlen(buf[i]) == 0) continue;

        // Find where the comments start
        comment_start = str_find_char(buf[i], '#', STR_FIND_FIRST);
        //if (comment_start == -1) comment_start = strlen(buf[i]);

        // Find where the command statement starts
        cmd_start = str_find_char(buf[i], '[', STR_FIND_FIRST);
        cmd_end = str_find_char(buf[i], ']', STR_FIND_LAST);
        if (cmd_end == cmd_start && cmd_start == -1) continue;

        // Check that a command can be interpreted
        assert(cmd_end > cmd_start && (comment_start < cmd_start || comment_start > cmd_end));
        if (comment_start != -1 && comment_start < cmd_start) continue;

        // Construct string option
        long  str_l = cmd_end - cmd_start;
        char* str_opt = calloc(str_l, 1);
        memcpy(str_opt, buf[i] + cmd_start + 1, str_l - 1);

        strbuffer opt_ln = str_split(str_opt, ',');
        free(str_opt);

        // Decode string
        if (STR_CMP(opt_ln[0], "TIMESTEPS") && strbuffer_len(opt_ln) == 2)
        {
            MODEL_SETTINGS.NT = strtoull(opt_ln[1], NULL, 10);
        }
        else if (STR_CMP(opt_ln[0], "DT") && strbuffer_len(opt_ln) == 2)
        {
            MODEL_SETTINGS.dt = strtod(opt_ln[1], NULL);
        }
        else if (STR_CMP(opt_ln[0], "FT") && strbuffer_len(opt_ln) == 2)
        {
            MODEL_SETTINGS.fNT = strtoull(opt_ln[1], NULL, 10);
        }
        else if (STR_CMP(opt_ln[0], "XDIM") && strbuffer_len(opt_ln) == 4)
        {
            MODEL_SETTINGS.xMin = strtod(opt_ln[1], NULL);
            MODEL_SETTINGS.xMax = strtod(opt_ln[2], NULL);
            MODEL_SETTINGS.dx = strtod(opt_ln[3], NULL);
        }
        else if (STR_CMP(opt_ln[0], "ZDIM") && strbuffer_len(opt_ln) == 4)
        {
            MODEL_SETTINGS.zMin = strtod(opt_ln[1], NULL);
            MODEL_SETTINGS.zMax = strtod(opt_ln[2], NULL);
            MODEL_SETTINGS.dz = strtod(opt_ln[3], NULL);
        }
        else if (STR_CMP(opt_ln[0], "BINOPTS") && strbuffer_len(opt_ln) == 5)
        {
            MODEL_SETTINGS.nbins = strtoul(opt_ln[1], NULL, 10);
            if (STR_CMP(opt_ln[2], "linexp"))
            {
                MODEL_SETTINGS.bgrid = LINEXP_RGRID;
            }
            else if (STR_CMP(opt_ln[2], "masmult"))
            {
                MODEL_SETTINGS.bgrid = LINMASS_RGRID;
            }
            else
            {
                fprintf(stderr, "Invalid bin grid option %s passed.\r\n", opt_ln[2]);
                return -1;
            }
            
            MODEL_SETTINGS.dr_lin = strtod(opt_ln[3], NULL);
            MODEL_SETTINGS.dr_ex = strtod(opt_ln[4], NULL);
        }
        else if (STR_CMP(opt_ln[0], "PROFILE") && strbuffer_len(opt_ln) > 1)
        {
            // Print Error The Specified profile format is lacking
            if (strbuffer_len(opt_ln) == 2)
            {
                fprintf(stderr, "Error: profile file \"%s\" has not been provided with a interpreter\r\n", opt_ln[1]);
                unsigned long j;
                char* stmp = str_trim(opt_ln[1], TRUE, TRUE);
                MODEL_SETTINGS.snd_file_nm = model_settings_add_str(stmp);
                MODEL_SETTINGS.snd_fmt = 0UL;
                free(stmp);

                for (j = 2; j < strbuffer_len(opt_ln); ++j)
                {
                    unsigned short sflag = 0U;
                    stmp = str_trim(opt_ln[j], TRUE, TRUE);

                    if (STR_CMP(stmp, "T")) sflag = 0x01;
                    else if (STR_CMP(stmp, "TD")) sflag = 0x02;
                    else if (STR_CMP(stmp, "P")) sflag = 0x04;
                    else if (STR_CMP(stmp, "Z")) sflag = 0x08;
                    else if (STR_CMP(stmp, "Q")) sflag = 0x10;
                    else if (STR_CMP(stmp, "TH")) sflag = 0x20;

                    free(stmp);
                    MODEL_SETTINGS.snd_fmt |= (sflag << (2 * (j - 2)));

                }
                return -2;
            }
        }
        else if (STR_CMP(opt_ln[0], "CCN") && strbuffer_len(opt_ln) == 3)
        {
            MODEL_SETTINGS.CCN_C0 = strtod(opt_ln[1], NULL);
            MODEL_SETTINGS.CCN_k = strtod(opt_ln[2], NULL);
        }
        else if (STR_CMP(opt_ln[0], "NCOUT") && strbuffer_len(opt_ln) >= 3)
        {
            char* tmpstr = str_trim(opt_ln[1], TRUE, TRUE);
            MODEL_SETTINGS.nc_output_nm = model_settings_add_str(tmpstr);
            free(tmpstr);
            if (strbuffer_len(opt_ln) >= 3) MODEL_SETTINGS.nc_write_freq = strtoul(opt_ln[2], NULL, 10);
            
        }
        else if (STR_CMP(opt_ln[0], "NCFLAGS") && strbuffer_len(opt_ln) >= 1)
        {
            unsigned long j;
            MODEL_SETTINGS.nc_flags = 0x0U;
            for (j = 1; j < strbuffer_len(opt_ln); ++j)
            {
                char* tmpstr = str_trim(opt_ln[j], TRUE, TRUE);
                if (STR_CMP(tmpstr, "L_FIELD")) MODEL_SETTINGS.nc_flags |= NC_OUTPUT_CLOUD_LIQUID;
                else if (STR_CMP(tmpstr, "RHOW_FIELD")) MODEL_SETTINGS.nc_flags |= NC_OUTPUT_MOMENTUM_W;
                else if (STR_CMP(tmpstr, "RHOU_FIELD")) MODEL_SETTINGS.nc_flags |= NC_OUTPUT_MOMENTUM_U;
                else if (STR_CMP(tmpstr, "NONE")) {MODEL_SETTINGS.nc_flags = 0; free(tmpstr); break;}
                free(tmpstr);
            }
        }
#ifdef OPTIONS_PRINT_UNKNOWN
        else
        {
            // Print out unknown Option
            fprintf(stderr, "Unknown or Incorrect Option: %s\r\n", opt_ln[0]);

        }
#endif
        


        free_strbuffer(opt_ln);
    }


    // Free Mem
    free(file_data);
    free_strbuffer(buf);
    return 0;
}

void fprint_opts(FILE* pf, __model_settings sets)
{
    const char* btype = (sets.bgrid == LINEXP_RGRID) ? "linexp" : "massmult";
    fprintf(pf, "-------------------------------------------------\r\n");
    fprintf(pf, "Time steps: %lu\r\n", sets.NT);
    fprintf(pf, "Time step: %.1f\r\n", sets.dt);
    fprintf(pf, "Fractional time steps: %u\r\n", sets.fNT);
    fprintf(pf, "Z-dimensions [%.1f, %.1f] m with step size %.1f\r\n", sets.zMin, sets.zMax, sets.dz);
    fprintf(pf, "X-dimensions [%.1f, %.1f] m with step size %.1f\r\n", sets.xMin, sets.xMax, sets.dx);
    fprintf(pf, "%lu Bins on a %s grid:\r\n\tLin. spacing %.3f\r\n\tExp. spacing %.3f\r\n", sets.nbins, btype, sets.dr_lin, sets.dr_ex);
    fprintf(pf, "Initial values\r\n\tp0: %.0f hPa\r\n\tz0: %.1f m\r\n", sets.p0 * 1e-2, sets.z0);
    // TODO add stream function 
    // TODO add flux parametrisation
    fprintf(pf, "CCN Options\r\n\tC0: %.1f cm^-3\r\n\tk: %.1f\r\n", sets.CCN_C0 * 1e-6, sets.CCN_k);
    fprintf(pf, "Environment background: \"%s\"\r\n", sets.snd_file_nm);
    fprintf(pf, "\tEnvironment read flags: 0x%llx\r\n", sets.snd_fmt);
    fprintf(pf, "NetCDF4 Info:\r\n");
    fprintf(pf, "\tOutput File: \"%s\"\r\n", sets.nc_output_nm);
    fprintf(pf, "\tWrite Frequency: %i\r\n", (int)sets.nc_write_freq);
    fprintf(pf, "\tExtra Flags: 0x%lx\r\n", sets.nc_flags);
    fprintf(pf, "-------------------------------------------------\r\n");
}