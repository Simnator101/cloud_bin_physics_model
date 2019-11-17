#include "../include/extio.h"

// Default settings
__model_settings MODEL_SETTINGS =
{
    // Run time
    7200UL,
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

    "./data/shallow_cumulus.nc",
    8,
    0U,

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
