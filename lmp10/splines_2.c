#include "splines_2.h"                                                            

#include <stdlib.h>
#include <math.h>

#define MALLOC_FAILED( P, SIZE ) (((P)=malloc( (SIZE)*sizeof( *(P))))==NULL)
#define PI 3.14159265

    int
alloc_spl (spline_t * spl, int n)
{
    spl->n = n;
    return MALLOC_FAILED (spl->x, spl->n)
        || MALLOC_FAILED (spl->f, spl->n)
        || MALLOC_FAILED (spl->a, spl->n)
        || MALLOC_FAILED (spl->b, spl->n);
}

    int
read_spl (FILE * inf, spline_t * spl)
{
    int i;
    if (fscanf (inf, "%d", &(spl->n)) != 1 || spl->n < 0)
        return 1;

    if (alloc_spl (spl, spl->n))
        return 1;

    for (i = 0; i < spl->n; i++)
        if (fscanf (inf, "%lf %lf %lf %lf", spl->x + i, spl->f + i, spl->a + i, spl->b +i) !=4)
            return 1;
    return 0;
}

    void
write_spl (spline_t * spl, FILE * ouf)
{
    int i;
    fprintf (ouf, "%d\n", spl->n);
    for (i = 0; i < spl->n; i++)
        fprintf (ouf, "%g %g %g %g\n", spl->x[i], spl->f[i], spl->a[i], spl->b[i]);
}

    double
value_spl (spline_t * spl, double x)
{
    int i;
    double f = 0.0;  
    int m;

    if (spl->n % 2 == 0)
        m = (spl->n) / 2 - 1;
    else
        m = (spl->n - 1) / 2 - 1;

    for (i = 1; i <= m; i++)
        f += (spl->a[i]) * cos(2 *PI *i *x / spl->n) + (spl->b[i]) * sin(2 *PI *i *x / spl->n);
    f += spl->a[0];

    return f;
}                      
