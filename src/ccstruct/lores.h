///////////////////////////////////////////////////////////////////////
// File:        lores.h
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
// The code based on Leptonica is
// 
//   Copyright (C) 2001 Leptonica.  All rights reserved.
// 
// It is licensed under the following conditions:
// 
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//   2. Redistributions in binary form must reproduce the above
//      copyright notice, this list of conditions and the following
//      disclaimer in the documentation and/or other materials
//      provided with the distribution.
// 
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL ANY
//   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
//   OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The relevant code is indicated below.
// 
///////////////////////////////////////////////////////////////////////

#ifndef TESSERACT_IMAGE_LORES_H_
#define TESSERACT_IMAGE_LORES_H_

#include "genericvector.h"      // for GenericVector, PointerVector, FileReader

class TBOX;
#include "allheaders.h"       // for Pix and other Leptonica structures

namespace tesseract {

class TFile;

// Enum to store the method used to scale a LoresImage
enum LoresScalingMethod {
  LSM_BOX = 0,
  LSM_BILINEAR = 1,
  LSM_BICUBIC = 2
};

// Class to hold information on a single low-resolution greyscale image:
// Caches low-resolution image as a Pix*, original resolution,
// scaling method.
// It performs scaling to 300 dpi according to the specified scaling method.
class LoresImage {
 public:
  LoresImage();
  LoresImage(Pix* image, int32_t resolution, int32_t target_resolution,
             LoresScalingMethod scaling_method, double blur);
  ~LoresImage();
  LoresImage(const LoresImage&);

  static void InitializeScaling(LoresScalingMethod scaling_method,
                                double blur);

  // Accessor functions
  LoresScalingMethod scaling_method() const {
    return scaling_method_;
  }
  double blur_amount() const {
    return blur_amount_;
  }
  
  // Writes to the given file. Returns false in case of error.
  bool Serialize(TFile* fp) const;
  // Reads from the given file. Returns false in case of error.
  bool DeSerialize(TFile* fp);
  // As DeSerialize, but only seeks past the data - hence a static method.
  static bool SkipDeSerialize(TFile* fp);

  // Returns a copy of the full scaled-up Pix image.
  // Must be pixDestroyed after use.
  Pix* GetFullScaledImage() const;
  
  // Gets a scaled image of a portion of the original image, scaled using the
  // stored scaling and blurring methods to achieve the specified
  // target_height.  The specified bounding box is relative to the scaled
  // whole image, when the lores image has been scaled to the originally
  // requested resolution.
  // The return value is the scaled Pix, which must be pixDestroyed after use.
  // If the resulting scale factor is too great or other errors occur,
  // nullptr will be returned.
  Pix* GetScaledImageBox(int32_t target_height, const TBOX& box) const;

 private:
  // Make the full scaled image.
  // This method assumes that all private variables are suitably initialised
  // This internal function is marked const, as it only touches mutable
  // members, and is called by const functions.
  void MakeFullScaledImage() const;

  // These are (almost) identical to the static methods
  // ImageData::{Set,Get}PixInternal, except that the Set function does not
  // pixDestroy the given Pix.  If the ImageData methods were not private,
  // we would use those here.

  // Saves the given Pix as a PNG-encoded string and destroys it.
  // In case of missing PNG support in Leptonica use PNM format,
  // which requires more memory.
  static l_int32 PixToData(Pix* pix, GenericVector<char>* image_data);
  // Returns the Pix image for the image_data. Must be pixDestroyed after use.
  static Pix* DataToPix(const GenericVector<char>& image_data);

 private:
  Pix* image_;                          // Lores greyscale image.
  // We could probably suffice with int16_t for many of the following,
  // but then we would have to be more careful with implicit type conversions
  int32_t resolution_;                  // Stated resolution of image.
  int32_t worig_;                       // Width of lores image.
  int32_t horig_;                       // Height of lores image.
  int32_t target_resolution_;           // Stated resolution of image.
  int32_t wtarget_;                     // Width of full target image.
  int32_t htarget_;                     // Height of full target image.
  double scale_factor_;                 // target_resolution_ / resolution_.
  static bool scaling_initialized_;     // Has the scaling been initialized?
  static bool scaling_warning_;       // Have we warned about a scaling change?
  static LoresScalingMethod scaling_method_;   // Method used for scaling.
  static double blur_amount_;           // Gaussian blur s.d. to apply.
  static int32_t kernel_halfsize_;      // Kernel size is 2 * this + 1.
  static L_KERNEL *gauss_kernel_;       // Gaussian blurring kernel.
};

}  // namespace tesseract


#endif  // TESSERACT_IMAGE_LORES_H_
