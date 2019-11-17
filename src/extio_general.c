#include "../include/extio.h"

/*General Tools*/
long file_size(FILE* pf)
{
    long sz;
    long cpos = ftell(pf);
    fseek(pf, 0, SEEK_END);
    sz = ftell(pf);
    fseek(pf, cpos, SEEK_SET);
    return sz;
}