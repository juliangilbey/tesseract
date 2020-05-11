///////////////////////////////////////////////////////////////////////
// File:        lores.cpp
// Description: Class to hold information about a single low-resolution
//              image and its scaled versions.
//              Structure based on imagedata.h, with algorithms from
//              Pillow and leptonica.
// Author:      Julian Gilbey
// Created:     Mon Apr 20 14:34:17 BST 2020
//
// (C) Copyright 2020, Julian Gilbey
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// The code based on Pillow is 
//       Copyright © 2010-2020 by Alex Clark and contributors
//
// Pillow is licensed under the PIL Software License:
// 
// By obtaining, using, and/or copying this software and/or its associated
// documentation, you agree that you have read, understood, and will comply
// with the following terms and conditions:
// 
// Permission to use, copy, modify, and distribute this software and its
// associated documentation for any purpose and without fee is hereby granted,
// provided that the above copyright notice appears in all copies, and that
// both that copyright notice and this permission notice appear in supporting
// documentation, and that the name of Secret Labs AB or the author not be
// used in advertising or publicity pertaining to distribution of the software
// without specific, written prior permission.
// 
// SECRET LABS AB AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
// SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
// IN NO EVENT SHALL SECRET LABS AB OR THE AUTHOR BE LIABLE FOR ANY SPECIAL,
// INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
// LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
// OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
// PERFORMANCE OF THIS SOFTWARE.
// 
// The relevant code is indicated below.
// 
///////////////////////////////////////////////////////////////////////

// Include automatically generated configuration file if running autoconf.
#ifdef HAVE_CONFIG_H
#include "config_auto.h"
#endif

#include "lores.h"

// TODO: go through these - most are probably unnecessary

#if defined(__MINGW32__)
#include <unistd.h>
#else
#include <thread>
#endif

#include "allheaders.h"  // for pixDestroy, pixGetHeight, pixGetWidth, lept_...
#include "helpers.h"     // for IntCastRounded, TRand, ClipToRange, Modulo
#include "rect.h"        // for TBOX, ICOORD, FCOORD
#include <cmath>         // for ceil
#include <limits>        // for maximum int limit checking
#include "tprintf.h"     // for tprintf

// Maximum scale factor to allow; beyond this probably means that there
// are no letters present, perhaps only some midline punctuation.
const float kMaxScale = 10.0;

namespace tesseract {

namespace {

// Bilinear and bicubic interpolation, for 8 bit greyscale only.
// All code in this anonymous namespace is from Pillow, version 7.0.0,
// source file libImaging/Resample.c, adapted to work with Leptonica
// data structures.
// It is translated from C code, hence the numerous reinterpret_casts.

struct filter {
  double (*filter)(double x);
  double support;
};

inline double box_filter(double x)
{
  if (x > -0.5 && x <= 0.5)
    return 1.0;
  return 0.0;
}

inline double bilinear_filter(double x)
{
  if (x < 0.0)
    x = -x;
  if (x < 1.0)
    return 1.0-x;
  return 0.0;
}

inline double bicubic_filter(double x)
{
  /* https://en.wikipedia.org/wiki/Bicubic_interpolation#Bicubic_convolution_algorithm */
  const double a = -0.5;

  if (x < 0.0)
    x = -x;
  if (x < 1.0)
    return ((a + 2.0) * x - (a + 3.0)) * x*x + 1;
  if (x < 2.0)
    return (((x - 5) * x + 8) * x - 4) * a;
  return 0.0;
}

struct filter BOXF = { box_filter, 0.5 };
struct filter BILINEAR = { bilinear_filter, 1.0 };
struct filter BICUBIC = { bicubic_filter, 2.0 };

/* 8 bits for result. Filter can have negative areas.
   In one cases the sum of the coefficients will be negative,
   in the other it will be more than 1.0. That is why we need
   two extra bits for overflow and int type. */
const int PRECISION_BITS (32 - 8 - 2);


/* Handles values form -640 to 639. */
l_uint8 _clip8_lookups[1280] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
    32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
    48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
    64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
    80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
    96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
    112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
    128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
    144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
    160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
    176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
    192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
    208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
    224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
    240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
};


l_uint8 *clip8_lookups = &_clip8_lookups[640];

inline l_uint8 clip8(int in)
{
  return clip8_lookups[in >> PRECISION_BITS];
}


int precompute_coeffs(int inSize, float in0, float in1, int outSize,
                      struct filter *filterp, int **boundsp, double **kkp)
{
  double support, scale, filterscale;
  double center, w, ww, ss;
  int xx, x, ksize, xmin, xmax;
  int *bounds;
  double *kk, *k;

  /* prepare for horizontal stretch */
  filterscale = scale = double(in1 - in0) / outSize;
  if (filterscale < 1.0) {
    filterscale = 1.0;
  }

  /* determine support size (length of resampling filter) */
  support = filterp->support * filterscale;

  /* maximum number of coeffs */
  ksize = int(ceil(support)) * 2 + 1;

  // check for overflow
  if (outSize > std::numeric_limits<int>::max() / (ksize * sizeof(double))) {
    return 0;
  }

  /* coefficient buffer */
  /* malloc check ok, overflow checked above */
  kk = static_cast<double *> (malloc(outSize * ksize * sizeof(double)));
  if (! kk) {
    return 0;
  }

  /* malloc check ok, ksize*sizeof(double) > 2*sizeof(int) */
  bounds = static_cast<int *> (malloc(outSize * 2 * sizeof(int)));
  if (! bounds) {
    free(kk);
    return 0;
  }

  for (xx = 0; xx < outSize; xx++) {
    center = in0 + (xx + 0.5) * scale;
    ww = 0.0;
    ss = 1.0 / filterscale;
    // Round the value
    xmin = int(center - support + 0.5);
    if (xmin < 0)
      xmin = 0;
    // Round the value
    xmax = int(center + support + 0.5);
    if (xmax > inSize)
      xmax = inSize;
    xmax -= xmin;
    k = &kk[xx * ksize];
    for (x = 0; x < xmax; x++) {
      w = filterp->filter((x + xmin - center + 0.5) * ss);
      k[x] = w;
      ww += w;
    }
    for (x = 0; x < xmax; x++) {
      if (ww != 0.0)
        k[x] /= ww;
    }
    // Remaining values should stay empty if they are used despite of xmax.
    for (; x < ksize; x++) {
      k[x] = 0;
    }
    bounds[xx * 2 + 0] = xmin;
    bounds[xx * 2 + 1] = xmax;
  }
  *boundsp = bounds;
  *kkp = kk;
  return ksize;
}


void normalize_coeffs(int outSize, int ksize, double *prekk)
{
  int x;
  l_int32 *kk;

  // use the same buffer for normalized coefficients
  // this is somewhat of a hack!
  kk = reinterpret_cast <l_int32 *>(prekk);

  for (x = 0; x < outSize * ksize; x++) {
    if (prekk[x] < 0) {
      kk[x] = int(-0.5 + prekk[x] * (1 << PRECISION_BITS));
    } else {
      kk[x] = int(0.5 + prekk[x] * (1 << PRECISION_BITS));
    }
  }
}


void resample_horizontal(Pix* pixOut, Pix* pixIn, int offset,
                         int ksize, int *bounds, double *prekk)
{
  int ss0;
  int xx, yy, x, xmin, xmax;
  l_int32 *k, *kk, wOut, hOut;
  l_uint8 **lineptrsIn, **lineptrsOut, *lineIn, *lineOut;

  /* We use the Leptonica usage pattern for 8-bit direct image access,
   * as described in the Leptonica source file pix1.c.  */

  lineptrsIn = reinterpret_cast<l_uint8 **>(pixGetLinePtrs(pixIn, nullptr));
  lineptrsOut = reinterpret_cast<l_uint8 **>(pixGetLinePtrs(pixOut, nullptr));
  wOut = pixGetWidth(pixOut);
  hOut = pixGetHeight(pixOut);
  
  // use the same buffer for normalized coefficients
  kk = reinterpret_cast<l_int32 *>(prekk);
  normalize_coeffs(wOut, ksize, prekk);

  /* The next part is much simpler than in the original, as we are
   * only dealing with 8-bit images. */

  for (yy = 0; yy < hOut; yy++) {
    lineIn = lineptrsIn[yy + offset];
    lineOut = lineptrsOut[yy];
    for (xx = 0; xx < wOut; xx++) {
      xmin = bounds[xx * 2 + 0];
      xmax = bounds[xx * 2 + 1];
      k = &kk[xx * ksize];
      ss0 = 1 << (PRECISION_BITS -1);
      for (x = 0; x < xmax; x++)
        ss0 += GET_DATA_BYTE(lineIn, x + xmin) * k[x];
      SET_DATA_BYTE(lineOut, xx, clip8(ss0));
    }
  }

  LEPT_FREE(lineptrsIn);
  LEPT_FREE(lineptrsOut);
}


void resample_vertical(Pix* pixOut, Pix* pixIn, int offset,
                       int ksize, int *bounds, double *prekk)
{
  int ss0;
  int xx, yy, y, ymin, ymax;
  l_int32 *k, *kk, wOut, hOut;
  l_uint8 **lineptrsIn, **lineptrsOut, *lineIn, *lineOut;

  /* We use the Leptonica usage pattern for 8-bit direct image access,
   * as described in the Leptonica source file pix1.c. */

  lineptrsIn = reinterpret_cast<l_uint8 **>(pixGetLinePtrs(pixIn, nullptr));
  lineptrsOut = reinterpret_cast<l_uint8 **>(pixGetLinePtrs(pixOut, nullptr));
  wOut = pixGetWidth(pixOut);
  hOut = pixGetHeight(pixOut);
  
  // use the same buffer for normalized coefficients
  kk = reinterpret_cast<l_int32 *>(prekk);
  normalize_coeffs(hOut, ksize, prekk);

  for (yy = 0; yy < hOut; yy++) {
    lineOut = lineptrsOut[yy];
    k = &kk[yy * ksize];
    ymin = bounds[yy * 2 + 0];
    ymax = bounds[yy * 2 + 1];
    for (xx = 0; xx < wOut; xx++) {
      ss0 = 1 << (PRECISION_BITS -1);
      for (y = 0; y < ymax; y++)
        ss0 += GET_DATA_BYTE(lineptrsIn[y + ymin], xx) * k[y];
      SET_DATA_BYTE(lineOut, xx, clip8(ss0));
    }
  }

  LEPT_FREE(lineptrsIn);
  LEPT_FREE(lineptrsOut);
}


/* Here, xsize and ysize are the requested size for the rescaled image,
 * filter specifies the filter to be used, and box_ul, box_lr
 * specify the portion of the original image to scale.  Note that
 * the y-coordinate is 0 at the top in leptonica (though it is at the
 * bottom in tesseract); hopefully this won't matter here!
 * (Note that we cannot use a TBOX to specify the box, as that takes only
 * integer coordinates.)
*/

Pix* rescale(Pix* pix, int xsize, int ysize, LoresScalingMethod filter,
             float left, float top, float right, float bottom)
{
  struct filter *filterp;
  Pix* pixTemp = nullptr;
  Pix* pixOut = nullptr;
  int i, need_horizontal, need_vertical;
  int ybox_first, ybox_last;
  int ksize_horiz, ksize_vert;
  int *bounds_horiz, *bounds_vert;
  double *kk_horiz, *kk_vert;

  /* check filter; we only implement bilinear and bicubic here,
     though PILLOW provides other filters too */
  switch (filter) {
  case LSM_BILINEAR:
    filterp = &BILINEAR;
    break;
  case LSM_BICUBIC:
    filterp = &BICUBIC;
    break;
  case LSM_BOX:
    filterp = &BOXF;
    break;
  default:
    return nullptr;
  }

  need_horizontal = xsize != pixGetWidth(pix) || left != 0 || right != xsize;
  need_vertical = ysize != pixGetHeight(pix) || top != 0 || bottom != ysize;

  ksize_horiz = precompute_coeffs(pixGetWidth(pix), left, right, xsize,
                                  filterp, &bounds_horiz, &kk_horiz);
  if (! ksize_horiz) {
    return nullptr;
  }

  ksize_vert = precompute_coeffs(pixGetHeight(pix), top, bottom, ysize,
                                 filterp, &bounds_vert, &kk_vert);
  if (! ksize_vert) {
    free(bounds_horiz);
    free(kk_horiz);
    return nullptr;
  }

  // First used row in the source image
  ybox_first = bounds_vert[0];
  // Last used row in the source image
  ybox_last = bounds_vert[ysize*2 - 2] + bounds_vert[ysize*2 - 1];

  /* two-pass resize, horizontal pass */
  if (need_horizontal) {
    // Shift bounds for vertical pass
    for (i = 0; i < ysize; i++) {
      bounds_vert[i * 2] -= ybox_first;
    }

    // The '8' here is because we are working solely with 8-bit images
    pixTemp = pixCreateNoInit(xsize, ybox_last - ybox_first, 8);
    if (pixTemp) {
      resample_horizontal(pixTemp, pix, ybox_first,
                          ksize_horiz, bounds_horiz, kk_horiz);
    }
    free(bounds_horiz);
    free(kk_horiz);
    if (! pixTemp) {
      free(bounds_vert);
      free(kk_vert);
      return nullptr;
    }
    pixOut = pix = pixTemp;
  } else {
    // Free in any case
    free(bounds_horiz);
    free(kk_horiz);
  }

  /* vertical pass */
  if (need_vertical) {
    // The '8' here is because we are working solely with 8-bit images
    pixOut = pixCreateNoInit(pixGetWidth(pix), ysize, 8);
    if (pixOut) {
      /* pix can be the original image or horizontally resampled one */
      resample_vertical(pixOut, pix, 0,
                        ksize_vert, bounds_vert, kk_vert);
    }
    /* it's safe to call pixDestroy with empty value
       if previous step was not performed. */
    pixDestroy(&pixTemp);
    free(bounds_vert);
    free(kk_vert);
    if (! pixOut) {
      return nullptr;
    }
  } else {
    // Free in any case
    free(bounds_vert);
    free(kk_vert);
  }
  
  /* none of the previous steps are performed, copying */
  if (! pixOut) {
    pixOut = pixCopy(nullptr, pix);
  }

  return pixOut;
}

} // end anonymous namespace

// The destructor must pixDestroy the created copy of the image.
LoresImage::LoresImage(Pix* image, int resolution, int target_resolution,
                       LoresScalingMethod scaling_method, double blur)
  : resolution_(resolution),
    target_resolution_(target_resolution),
    scaling_method_(scaling_method),
    blur_amount_(blur),
    image_(nullptr),
    scaled_image_(nullptr)
{
  Pix* pixtemp;

  // This forces the image to be 8-bit greyscale
  if (pixGetColormap(image))
    image_ = pixRemoveColormap(image, REMOVE_CMAP_TO_GRAYSCALE);
  else
    image_ = pixConvertTo8(image, false);

  if (image_ == nullptr)
    return;

  worig_ = pixGetWidth(image_);
  horig_ = pixGetHeight(image_);

  scale_factor_ = (double) target_resolution / resolution;
  if (scale_factor_ > kMaxScale)
    return;

  wtarget_ = IntCastRounded(worig_ * scale_factor_);
  htarget_ = IntCastRounded(horig_ * scale_factor_);

  pixtemp = rescale(image_, wtarget_, htarget_, scaling_method_,
                    0, 0, worig_, horig_);
  if (pixtemp == nullptr)
    return;

  if (blur_amount_ > 0) {
    // This is perfectly adequate for the accuracy we need
    kernel_halfsize_ = IntCastRounded(4 * blur);
    gauss_kernel_ = makeGaussianKernel(kernel_halfsize_, kernel_halfsize_,
                                       blur, 1.0);
    scaled_image_ = pixConvolve(pixtemp, gauss_kernel_, 8, 1);
  } else {
    kernel_halfsize_ = 0;
    gauss_kernel_ = nullptr;
    scaled_image_ = pixClone(pixtemp);
  }

  pixDestroy(&pixtemp);
}


LoresImage::~LoresImage()
{
  pixDestroy(&image_);
  pixDestroy(&scaled_image_);
}


Pix* LoresImage::GetImage() const
{
  return image_ ? pixClone(image_) : nullptr;
}


Pix* LoresImage::GetFullScaledImage() const
{
  return scaled_image_ ? pixClone(scaled_image_) : nullptr;
}


// Gets a scaled image of a portion of the original image, scaled using the
// stored scaling and blurring methods to achieve the specified
// target_height.  The specified bounding box is relative to the scaled
// whole image, when the lores image has been scaled to the originally
// requested resolution; it is given in Tesseract coordinates (y=0 at the
// bottom of the image).
// The return value is the scaled Pix, which must be pixDestroyed after use.
// If the resulting scale factor is too great or other errors occur,
// nullptr will be returned.
Pix* LoresImage::GetScaledImageBox(int target_height, const TBOX& box) const
{
  // We rescale a portion of the lores image large enough to give the
  // requested box plus enough overshoot to account for the scaling
  // method and blurring, to avoid edge effects.
  // Note that overshoot is measured in lores image pixels.

  float overshoot;
  switch (scaling_method_) {
  case LSM_BILINEAR:
    overshoot = 1;
    break;
  case LSM_BICUBIC:
    overshoot = 2;
    break;
  case LSM_BOX:
    overshoot = 0;
    break;
  }

  // The kernel_halfsize_ is in target pixels.  The effective
  // resolution will be somewhat different from the target_resolution,
  // but it is unlikely to be very different.  We therefore allow some
  // leeway in this part of the overshoot.
  overshoot += 1.5 * kernel_halfsize_ / scale_factor_;

  // We convert into Pillow/Leptonica coordinates at this point;
  // compare this with the Tesseract::GetRectImage function in linerec.cpp;
  // note, though, that our rescale function above requires the coordinates
  // of the bounding box rather than its width and height.
  float ulx0 = box.left() / scale_factor_ - overshoot;
  float uly0 = (htarget_ - box.top()) / scale_factor_ - overshoot;
  float lrx0 = box.right() / scale_factor_ + overshoot;
  float lry0 = (htarget_ - box.bottom()) / scale_factor_ + overshoot;
  float ulx = std::max(ulx0, 0.0f);
  float uly = std::max(uly0, 0.0f);
  float lrx = std::min(lrx0, static_cast<float>(worig_));
  float lry = std::min(lry0, static_cast<float>(horig_));
  
  float scaling = scale_factor_ * target_height / box.height();
  if (scaling > kMaxScale) {
    tprintf("scaling factor %g too large\n", scaling);
    return nullptr;
  }

  float image_ratio = float(box.right() - box.left()) /
                        (box.top() - box.bottom());
  int target_width = IntCastRounded(target_height * image_ratio);
  int xsize_over = IntCastRounded(scaling * (lrx - ulx));
  int ysize_over = IntCastRounded(scaling * (lry - uly));
  
  Pix *pixtemp = rescale(image_, xsize_over, ysize_over, scaling_method_,
                         ulx, uly, lrx, lry);
  if (!pixtemp) {
    tprintf("Could not rescale image to requested size "
            "(ul = (%.2f, %.2f), lr = (%.2f, %.2f))\n",
            ulx, uly, lrx, lry);
    return nullptr;
  }

  if (blur_amount_ > 0) {
    Pix *pixblur = pixConvolve(pixtemp, gauss_kernel_, 8, 1);
    pixDestroy(&pixtemp);
    pixtemp = pixblur;
  }

  // Now the pixtemp is the rescaled image with some overshoot.
  // The intended overshoot was "overshoot", but we had to trim the
  // overshoot by (ulx - ulx0) on the left and by (uly - uly0) on the top.
  // We therefore remove overshoot - (ulx - ulx0) from the left and
  // similarly from the top, rounded down to avoid removing more pixels
  // than we have available.  (One pixel difference should have little
  // impact on the recognition results.)

  int left = int(overshoot * scaling - (ulx - ulx0));
  int top = int(overshoot * scaling - (uly - uly0));

  Box *boxtemp = boxCreateValid(left, top, target_width, target_height);
  if (!boxtemp) {
    pixDestroy(&pixtemp);
    tprintf("Could not create valid box "
            "(ul = (%d, %d), (w,h) = (%d, %d))\n",
            left, top, target_width, target_height);
    return nullptr;
  }
  
  Pix *pixd = pixClipRectangle(pixtemp, boxtemp, nullptr);
  if (!pixd) {
    tprintf("Could not run pixClipRectangle: "
            "pixtemp (w,h) = (%d, %d), "
            "clip box = (ul = (%d, %d), (w,h) = (%d, %d))\n",
            pixGetWidth(pixtemp), pixGetHeight(pixtemp),
            left, top, target_width, target_height);
  }
  boxDestroy(&boxtemp);
  pixDestroy(&pixtemp);

  return pixd;
}


}  // namespace tesseract.
