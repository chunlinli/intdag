/**********
    C++ Routines for Linear Algebra Operations.

    Copyright (C) 2020-2021  Chunlin Li

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.
**********/

#ifndef _UTILS_H_
#define _UTILS_H_

#ifdef __arm64__

#include <Accelerate/Accelerate.h> // Apple arm

inline double inner_product_simd(const double *array1,
                                 const double *array2,
                                 const int length,
                                 double init)
{
    return init + cblas_ddot(length, array1, 1, array2, 1);
}

inline double inner_product_simd(const double *array1,
                                 const double *array2,
                                 const double *array3,
                                 const int length,
                                 double init)
{
    for (int i = 0; i < length; ++i)
    {
        init += array1[i] * array2[i] * array3[i];
    }
    return init;
}

inline void vec_add_simd(const double *array1,
                         const double scalar1,
                         const double scalar2,
                         const int length,
                         double *dest)
{
    
    for (int i = 0; i < length; ++i)
    {
        dest[i] = scalar1 * array1[i] + scalar2;
    }
}

inline void vec_add_simd(const double *array1,
                         const double *array2,
                         const double scalar1,
                         const double scalar2,
                         const int length,
                         double *dest)
{
    
    for (int i = 0; i < length; ++i)
    {
        dest[i] = scalar1 * array1[i] + scalar2 * array2[i];
    }
}

#endif

#ifdef __x86_64__

#include <immintrin.h>

inline double inner_product_simd(const double *array1,
                                 const double *array2,
                                 const int length,
                                 double init)
{
    __m256d vsum = _mm256_setzero_pd();
    int i = 0;
    if (length >= 8 * 4)
    {
        __m256d vsum1 = _mm256_setzero_pd();
        __m256d vsum2 = _mm256_setzero_pd();
        __m256d vsum3 = _mm256_setzero_pd();
        __m256d vsum4 = _mm256_setzero_pd();
        __m256d vsum5 = _mm256_setzero_pd();
        __m256d vsum6 = _mm256_setzero_pd();
        __m256d vsum7 = _mm256_setzero_pd();
        __m256d vsum8 = _mm256_setzero_pd();
        for (; i + 8 * 4 - 1 < length; i += 8 * 4)
        { // could unroll further
            vsum1 = _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i),
                                    _mm256_loadu_pd(array2 + i), vsum1);
            vsum2 = _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 4),
                                    _mm256_loadu_pd(array2 + i + 4), vsum2);
            vsum3 = _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 8),
                                    _mm256_loadu_pd(array2 + i + 8), vsum3);
            vsum4 = _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 12),
                                    _mm256_loadu_pd(array2 + i + 12), vsum4);
            vsum5 = _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 16),
                                    _mm256_loadu_pd(array2 + i + 16), vsum5);
            vsum6 = _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 20),
                                    _mm256_loadu_pd(array2 + i + 20), vsum6);
            vsum7 = _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 24),
                                    _mm256_loadu_pd(array2 + i + 24), vsum7);
            vsum8 = _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 28),
                                    _mm256_loadu_pd(array2 + i + 28), vsum8);
        }
        vsum1 = _mm256_add_pd(vsum1, vsum2);
        vsum3 = _mm256_add_pd(vsum3, vsum4);
        vsum5 = _mm256_add_pd(vsum5, vsum6);
        vsum7 = _mm256_add_pd(vsum7, vsum8);
        vsum1 = _mm256_add_pd(vsum1, vsum3);
        vsum5 = _mm256_add_pd(vsum5, vsum7);
        vsum = _mm256_add_pd(vsum1, vsum5);
    }
    for (; i + 3 < length; i += 4)
    { // could unroll further
        vsum = _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i), _mm256_loadu_pd(array2 + i), vsum);
    }
    double buffer[4];
    _mm256_storeu_pd(buffer, vsum);
    init += buffer[0] + buffer[1] + buffer[2] + buffer[3];
    for (; i < length; ++i)
    {
        init += array1[i] * array2[i];
    }
    return init;
}

inline double inner_product_simd(const double *array1,
                                 const double *array2,
                                 const double *array3,
                                 const int length,
                                 double init)
{
    __m256d vsum = _mm256_setzero_pd();
    int i = 0;
    if (length >= 8 * 4)
    {
        __m256d vsum1 = _mm256_setzero_pd();
        __m256d vsum2 = _mm256_setzero_pd();
        __m256d vsum3 = _mm256_setzero_pd();
        __m256d vsum4 = _mm256_setzero_pd();
        __m256d vsum5 = _mm256_setzero_pd();
        __m256d vsum6 = _mm256_setzero_pd();
        __m256d vsum7 = _mm256_setzero_pd();
        __m256d vsum8 = _mm256_setzero_pd();
        for (; i + 8 * 4 - 1 < length; i += 8 * 4)
        { // could unroll further

            vsum1 = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(array1 + i),
                                                  _mm256_loadu_pd(array2 + i)),
                                    _mm256_loadu_pd(array3 + i), vsum1);
            vsum2 = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(array1 + i + 4),
                                                  _mm256_loadu_pd(array2 + i + 4)),
                                    _mm256_loadu_pd(array3 + i + 4), vsum2);
            vsum3 = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(array1 + i + 8),
                                                  _mm256_loadu_pd(array2 + i + 8)),
                                    _mm256_loadu_pd(array3 + i + 8), vsum3);
            vsum4 = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(array1 + i + 12),
                                                  _mm256_loadu_pd(array2 + i + 12)),
                                    _mm256_loadu_pd(array3 + i + 12), vsum4);
            vsum5 = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(array1 + i + 16),
                                                  _mm256_loadu_pd(array2 + i + 16)),
                                    _mm256_loadu_pd(array3 + i + 16), vsum5);
            vsum6 = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(array1 + i + 20),
                                                  _mm256_loadu_pd(array2 + i + 20)),
                                    _mm256_loadu_pd(array3 + i + 20), vsum6);
            vsum7 = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(array1 + i + 24),
                                                  _mm256_loadu_pd(array2 + i + 24)),
                                    _mm256_loadu_pd(array3 + i + 24), vsum7);
            vsum8 = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(array1 + i + 28),
                                                  _mm256_loadu_pd(array2 + i + 28)),
                                    _mm256_loadu_pd(array3 + i + 28), vsum8);
        }
        vsum1 = _mm256_add_pd(vsum1, vsum2);
        vsum3 = _mm256_add_pd(vsum3, vsum4);
        vsum5 = _mm256_add_pd(vsum5, vsum6);
        vsum7 = _mm256_add_pd(vsum7, vsum8);
        vsum1 = _mm256_add_pd(vsum1, vsum3);
        vsum5 = _mm256_add_pd(vsum5, vsum7);
        vsum = _mm256_add_pd(vsum1, vsum5);
    }
    for (; i + 3 < length; i += 4)
    { // could unroll further
        vsum = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_loadu_pd(array1 + i),
                                             _mm256_loadu_pd(array2 + i)),
                               _mm256_loadu_pd(array3 + i), vsum);
    }
    double buffer[4];
    _mm256_storeu_pd(buffer, vsum);
    init += buffer[0] + buffer[1] + buffer[2] + buffer[3];
    for (; i < length; ++i)
    {
        init += array1[i] * array2[i] * array3[i];
    }
    return init;
}

inline void vec_add_simd(const double *array1,
                         const double scalar1,
                         const double scalar2,
                         const int length,
                         double *dest)
{
    __m256d a = _mm256_set_pd(scalar1, scalar1, scalar1, scalar1);
    __m256d b = _mm256_set_pd(scalar2, scalar2, scalar2, scalar2);

    int i = 0;
    if (length >= 8 * 4)
    {
        for (; i + 8 * 4 - 1 < length; i += 8 * 4)
        {
            _mm256_storeu_pd(dest + i, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i), a, b));
            _mm256_storeu_pd(dest + i + 4, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 4), a, b));
            _mm256_storeu_pd(dest + i + 8, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 8), a, b));
            _mm256_storeu_pd(dest + i + 12, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 12), a, b));
            _mm256_storeu_pd(dest + i + 16, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 16), a, b));
            _mm256_storeu_pd(dest + i + 20, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 20), a, b));
            _mm256_storeu_pd(dest + i + 24, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 24), a, b));
            _mm256_storeu_pd(dest + i + 28, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 28), a, b));
        }
    }
    for (; i + 3 < length; i += 4)
    {
        _mm256_storeu_pd(dest + i, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i), a, b));
    }
    for (; i < length; ++i)
    {
        dest[i] = scalar1 * array1[i] + scalar2;
    }
}

inline void vec_add_simd(const double *array1,
                         const double *array2,
                         const double scalar1,
                         const double scalar2,
                         const int length,
                         double *dest)
{
    __m256d a = _mm256_set_pd(scalar1, scalar1, scalar1, scalar1);
    __m256d b = _mm256_set_pd(scalar2, scalar2, scalar2, scalar2);

    int i = 0;
    if (length >= 8 * 4)
    {
        for (; i + 8 * 4 - 1 < length; i += 8 * 4)
        {
            _mm256_storeu_pd(dest + i, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i), a,
                                                       _mm256_mul_pd(_mm256_loadu_pd(array2 + i), b)));
            _mm256_storeu_pd(dest + i + 4, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 4), a,
                                                           _mm256_mul_pd(_mm256_loadu_pd(array2 + i + 4), b)));
            _mm256_storeu_pd(dest + i + 8, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 8), a,
                                                           _mm256_mul_pd(_mm256_loadu_pd(array2 + i + 8), b)));
            _mm256_storeu_pd(dest + i + 12, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 12), a,
                                                            _mm256_mul_pd(_mm256_loadu_pd(array2 + i + 12), b)));
            _mm256_storeu_pd(dest + i + 16, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 16), a,
                                                            _mm256_mul_pd(_mm256_loadu_pd(array2 + i + 16), b)));
            _mm256_storeu_pd(dest + i + 20, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 20), a,
                                                            _mm256_mul_pd(_mm256_loadu_pd(array2 + i + 20), b)));
            _mm256_storeu_pd(dest + i + 24, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 24), a,
                                                            _mm256_mul_pd(_mm256_loadu_pd(array2 + i + 24), b)));
            _mm256_storeu_pd(dest + i + 28, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i + 28), a,
                                                            _mm256_mul_pd(_mm256_loadu_pd(array2 + i + 28), b)));
        }
    }
    for (; i + 3 < length; i += 4)
    {
        _mm256_storeu_pd(dest + i, _mm256_fmadd_pd(_mm256_loadu_pd(array1 + i), a,
                                                   _mm256_mul_pd(_mm256_loadu_pd(array2 + i), b)));
    }
    for (; i < length; ++i)
    {
        dest[i] = scalar1 * array1[i] + scalar2 * array2[i];
    }
}

#endif

inline double accumulate(const double *array,
                         const int length,
                         double init)
{
    for (int i = 0; i < length; ++i)
    {
        init += array[i];
    }
    return init;
}


// compute soft-thresholding function
inline double soft_thresh(double init, double thresh)
{
    if (init > thresh)
        init -= thresh;
    else if (init < -thresh)
        init += thresh;
    else
        init = 0.0;
    return init;
}

#endif
