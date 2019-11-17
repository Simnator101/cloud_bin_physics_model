#include "../include/extio.h"

long str_find_char(const char* str, const char delim, int mode)
{
    assert(mode == STR_FIND_FIRST || mode == STR_FIND_LAST);
    long fp = -1;
    unsigned long i;
    for (i = 0; i < strlen(str); ++i)
        if (str[i] == delim) 
        {
            if (mode == STR_FIND_LAST)
                fp = (long)i;
            else
                fp = (fp < 0) ? (long)i : fp;
        }

    return fp;
}

char* str_trim(char* str, int left, int right)
{
    unsigned long li = 0;
    unsigned long ri = strlen(str);
    if (left) while (isspace(str[li])) ++li;
    if (right) while (ri > 0 && isspace(str[ri])) --ri;

    char* ndata = calloc(ri - li + 1, 1);
    memcpy(ndata, str + li, ri - li);
    return ndata;
}

int str_eval_func(char* str, int(*func)(int))
{
    unsigned long i;
    for (i = 0; i < strlen(str); ++i) if (func(str[i]) != 0) return 1;
    return 0;
}

unsigned long strbuffer_len(strbuffer buffer)
{
    unsigned long ln = 0;
    while (buffer[ln] != NULL) ++ln;
    return ln;
}

void free_strbuffer(strbuffer strb)
{
    if (strb == NULL) return;
    unsigned long ln = strbuffer_len(strb);
    unsigned long i;
    for (i = 0; i < ln; ++i)
        free(strb[i]);
    free(strb);
}

strbuffer str_split(const char* str, const char delim)
{
    size_t sstrl = 1;
    size_t strl = strlen(str);
    unsigned long i, ri, li;
    for (i = 0; i < strl; ++i)
        if (str[i] == delim) ++sstrl;

    strbuffer buffer = allocate_strbuffer(sstrl + 1);
    for (i = 0, li = 0, ri = 0; ri < strl; ++ri)
    {
        if (str[ri] == delim)
        {
            // The length of characters to copy
            // Buffer size is c + 1
            size_t c = ri - li;
            buffer[i] = calloc(c + 1, 1);
            memcpy(buffer[i], str + li, c);

            li = ri + 1;
            ++i;
        }
    }

    // Add Final String
    buffer[sstrl-1] = calloc(ri - li + 1, 1);
    memcpy(buffer[sstrl-1], str + li, ri - li);
    return buffer;
}