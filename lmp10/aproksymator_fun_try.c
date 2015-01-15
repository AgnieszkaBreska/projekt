#include "makespl.h"                                                            
/*#include "piv_ge_solver.h"*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#define PI 3.14159265
    void
make_spl(points_t * pts, spline_t * spl)
{

    double         *x = pts->x;
    double         *y = pts->y;
    int             i, j, m;
    char *mEnv = getenv( "APPROX_BASE_SIZE" );

    if (mEnv != NULL && atoi(mEnv) > 0 && atoi(mEnv) <= ((pts->n - 1) / 2 - 1)) 
        m = atoi(mEnv);
    else {
        if (pts->n >= 4) {
            if (pts->n % 2 == 0)
                m = (pts->n) / 2 - 1;
            else
                m = (pts->n - 1) / 2 - 1;
        } else
            m = 1;
    }   

    if (alloc_spl (spl,pts->n) == 0) {
        spl->x = x;
        spl->f = y;
        spl->a[0] = 0.0;
        for (i = 1; i <= m; i++) {
            spl->a[i] = 0.0;
            spl->b[i] = 0.0;

            for (j = 1; j <= spl->n; j++) {
                spl->a[i] += (pts->y[j]) *( cos(2 *PI *i *j / pts->n));
                spl->b[i] += (pts->y[j]) *( sin(2 *PI *i *j / pts->n));
            }   
            spl->a[i] *= 2;
            spl->b[i] *= 2;
            spl->a[i] /= spl->n;
            spl->b[i] /= spl->n;
        }   
        for(j = 1; j < spl->n; j++){
            spl->a[0] += pts->y[j];
        }   
        spl->a[0] /= spl->n;
    }   
}

