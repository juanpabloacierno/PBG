//#include <stdio.h>
//#include <ctype.h>
//#include <string.h>
//#include <stdlib.h>
//#include <math.h>
//#include "sing_val.h"



BOOL SingValDecomp(double **, int, int, double *, double **);

#define MAX_STEP 30

/* Algorithm adapted from Numerical Recipes. */

/* computes sqrt(a^2 + b^2) without destructive overflow */

static double at, bt, ct;
#define PYTHAG(a, b) ((at = fabsf(a)) > (bt = fabsf(b)) ? \
    (ct = bt / at, at * sqrtf(1.0f + ct * ct)) : \
    (bt ? (ct = at / bt, bt * sqrtf(1.0f + ct * ct)) : 0.0f))



BOOL SingValDecomp(double **a, int m, int n, double *w, double **v)
/****************************************************************************/
/* Given a matrix a[0..m-1][0..n-1], this routine computes its singular     */
/* value decomposition, A = U W V'. The Matrix U replaces a on output. The  */
/* diagonal matrix of singular values W is output as a vector w[0..n-1].    */
/* The matrix V is output as v[0..n-1][0..n-1]. m must be greater or equal  */
/* to n; if it is smaller, then a should be filled up to square with zero   */
/* rows.                                                                    */
/****************************************************************************/
{
  int i, its, j, jj, k, l, nm;
  BOOL flag;
  double c, f, h, s, x, y, z;
  double anorm = 0.0f, g = 0.0f, scale = 0.0f;
  double *rv1;

  if (m < n)
    return FALSE;

 rv1 = (double *) malloc(n * sizeof(*rv1));



  /* Housholder reduction to bidiagonal form. */
  for (i = 0; i < n; i++) {
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0f;

    if (i < m) {
      for (k = i; k < m; k++)
        scale += fabsf(a[k][i]);

      if (scale) {
        for (k = i; k < m; k++) {
          a[k][i] /= scale;
          s += a[k][i] * a[k][i];
        }

        f = a[i][i];
        g = - SIGN(sqrtf(s), f);
        h = f * g - s;
        a[i][i] = f - g;

        if (i != n - 1) {
          for (j = l; j < n; j++) {
                for (s = 0.0f, k = i; k < m; k++)
              s += a[k][i] * a[k][j];
            f = s / h;
            for (k = i; k < m; k++)
              a[k][j] += f * a[k][i];
          }
        }

        for (k = i; k < m; k++)
          a[k][i] *= scale;
      }
    }

    w[i] = scale * g;
    g = s = scale = 0.0f;

    if (i < m && i != n - 1) {
      for (k = l; k < n; k++)
        scale += fabsf(a[i][k]);

      if (scale) {
        for (k = l; k < n; k++) {
          a[i][k] /= scale;
          s += a[i][k] * a[i][k];
        }

        f = a[i][l];
        g = - SIGN(sqrtf(s), f);
        h = f * g - s;
        a[i][l] = f - g;

        for (k = l; k < n; k++)
          rv1[k] = a[i][k] / h;

        if (i != m - 1)
          for (j = l; j < m; j++) {
            for (s = 0.0f, k = l; k < n; k++)
              s += a[j][k] * a[i][k];
            for (k = l; k < n; k++)
              a[j][k] += s * rv1[k];
          }

        for (k = l;k < n; k++)
          a[i][k] *= scale;
      }
    }

    anorm = FMAX(anorm, (fabsf(w[i]) + fabsf(rv1[i])));
  }

   /* Accumulation of right hand transformations. */
  for (i = n - 1; i >= 0; i--) {
    if (i < n - 1) {
      if (g != 0.0f) {
        for (j = l; j < n; j++)
          /* Double division to avoid possible underflow. */
          v[j][i] = (a[i][j] / a[i][l]) / g;

        for (j = l; j < n; j++) {
          for (s = 0.0f, k = l; k < n; k++)
            s += a[i][k] * v[k][j];
          for (k = l; k < n; k++)
            v[k][j] += s * v[k][i];
        }
      }

      for (j = l; j < n; j++)
        v[i][j] = v[j][i] = 0.0f;
    }

    v[i][i] = 1.0f;
    g = rv1[i];
    l = i;
  }

  /* Accumulation of left hand transformations. */
  for (i = n - 1; i >= 0; i--) {
    l = i + 1;
    g = w[i];

    if (i < n)
      for (j = l; j < n; j++)
        a[i][j] = 0.0;

    if (g != 0.0f) {
      g = 1.0f / g;

      if (i != n - 1) {
        for (j = l; j < n; j++) {
          for (s = 0.0f, k = l; k < m; k++)
            s += a[k][i] * a[k][j];
          f = (s / a[i][i]) * g;
          for (k = i; k < m; k++)
            a[k][j] += f * a[k][i];
        }
      }

      for (j = i; j < m; j++)
        a[j][i] *= g;
    } else {
      for (j = i; j < m; j++)
        a[j][i] = 0.0f;
    }

    a[i][i] += 1.0f;
  }

  /* Diagonalization of the bidiagonal form */
  for (k = n - 1; k >= 0; k--) {           /* Loop over singular values.    */
    for (its = 0; its < MAX_STEP; its++) { /* Loop over allowed iterations. */
      flag = TRUE;

      for (l = k; l >= 0; l--) {           /* Test for splitting.           */
        nm = l - 1;
        if ((fabsf(rv1[l]) + anorm) == anorm) {
          flag = FALSE;
          break;
        }

        if ((fabsf(w[nm]) + anorm) == anorm)
          break;
      }

      if (flag) {
        c = 0.0f;
        s = 1.0f;
        for (i = l; i <= k; i++) {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];
          if ((fabsf(f) + anorm) == anorm)
            continue;

          g = w[i];
          h = PYTHAG(f, g);
          w[i] = h;
          h = 1.0f / h;
          c = g * h;
          s = - f * h;

          for (j = 0; j < m; j++) {
            y = a[j][nm];
            z = a[j][i];
            a[j][nm] = y * c + z * s;
            a[j][i] = z * c - y * s;
          }
        }
      }

      z = w[k];

      if (l == k) {                        /* Convergence.                  */
        if (z < 0.0f) {             /* Singular value is made non negative. */
          w[k] = - z;
          for (j = 0; j < n; j++)
            v[j][k] = - v[j][k];
        }

        break;
      }

      if (its == MAX_STEP) {
        free(rv1);
        return FALSE;
      }

      /* Shift from bottom 2-by-2 minor. */
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0f * h * y);
      g = PYTHAG(f, 1.0f);
      f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

      /* Next QR transformation. */
      c = s = 1.0f;
      for (j = l; j <= nm; j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = PYTHAG(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;

        for (jj = 0;jj < n; jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z * s;
          v[jj][i] = z * c - x * s;
        }

        z = PYTHAG(f, h);
        w[j] = z;
        if (z != 0.0f) {            /* Rotation can be arbitrary if z == 0 */
          z = 1.0f / z;
          c = f * z;
          s = h * z;
        }

        f = c * g + s * y;
        x = c * y - s * g;

        for (jj = 0; jj < m; jj++) {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = y * c + z * s;
          a[jj][i] = z * c - y * s;
        }
      }

      rv1[l] = 0.0f;
      rv1[k] = f;
      w[k] = x;
    }
  }

  free(rv1);

  return TRUE;
}
