/* gransimg.c - GranSim TIFF image output.
 *
 * Copyright (C) 2024 Kirschner Laboratory, Pariksheet Nanda <pnanda@umich.edu>
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <string.h>

#include <tiffio.h>
#include <omp.h>


/** Change this to your GranSim grid size */
const int SZ = 300;

typedef enum {
  AGENT,
  CHEMOKINE
} im_t;


/**
 * @internal
 * @brief Common code for functions write_im_agent() and write_im_chemokine().
 *
 * Opens a writable single-image TIFF file and writes the TIFF header.
 *
 * @param tif Output pointer to the TIFF file
 * @param name Path to the TIFF file to create or overwrite
 * @param im_type Data to store in image of integer or floating point types
 *
 * @return 0 on success
 *
 * @see im_end
 */
int
im_begin(TIFF** tif, char* name, im_t im_type)
{
  *tif = TIFFOpen(name, "w");
#ifdef DEBUG
  printf("TIFFOpen(\"...\", \"w\") = %p\n", tif);
#endif

  /* Match ImageJ2's TIFF tags. */
  TIFFSetField(*tif, TIFFTAG_IMAGEWIDTH, SZ);
  TIFFSetField(*tif, TIFFTAG_IMAGELENGTH, SZ);
  TIFFSetField(*tif, TIFFTAG_BITSPERSAMPLE, 32);
  TIFFSetField(*tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  TIFFSetField(*tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(*tif, TIFFTAG_SAMPLESPERPIXEL, 1);
  TIFFSetField(*tif, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(*tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  switch (im_type) {
  case AGENT:
    TIFFSetField(*tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    break;
  case CHEMOKINE:
    TIFFSetField(*tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
    break;
  }

  return 0;
}


/**
 * @internal
 * @brief Common code for functions write_im_agent() and write_im_chemokine().
 *
 * Flushes and closes a TIFF file.
 *
 * @param tif Output pointer to the TIFF file
 *
 * @return 0 on success
 *
 * @see im_begin
 */
int
im_end(TIFF** tif)
{
  TIFFClose(*tif);

  return 0;
}

/**
 * @brief Write agent data as a single-image TIFF.
 *
 * @param name Path to the TIFF file to create or overwrite
 * @param im Dense single-dimensional array of data to store
 *
 * @return 0 on success
 *
 * @see write_im_chemokine()
 */
int
write_im_agent(char* name, uint32_t* im)
{
  TIFF* tif;
  im_begin(&tif, name, AGENT);
  
  if (tif) {
    tmsize_t bufsz = TIFFStripSize(tif);
    uint32_t *buf = (uint32_t*) _TIFFmalloc(bufsz);
    for (size_t i = 0; i < SZ; ++i) {
      for (int j = 0; j < SZ; ++j) {
	buf[j] = im[i * SZ + j];
      }
      TIFFWriteEncodedStrip(tif, i, buf, bufsz);
    }
    _TIFFfree(buf);
  }

  im_end(&tif);

  return 0;
}


/**
 * @brief Write chemokine data as a single-image TIFF.
 *
 * @param name Path to the TIFF file to create or overwrite
 * @param im Dense single-dimensional array of data to store
 *
 * @return 0 on success
 *
 * @see write_im_agent()
 */
int
write_im_chemokine(char* name, float* im) {
  TIFF* tif;
  im_begin(&tif, name, CHEMOKINE);
  
  if (tif) {
    tmsize_t bufsz = TIFFStripSize(tif);
    float *buf = (float*) _TIFFmalloc(bufsz);
    for (size_t i = 0; i < SZ; ++i) {
      for (int j = 0; j < SZ; ++j) {
	buf[j] = im[i * SZ + j];
      }
      TIFFWriteEncodedStrip(tif, i, buf, bufsz);
    }
    _TIFFfree(buf);
  }

  im_end(&tif);

  return 0;
}


/**
 * @brief Export model directory agents and grid dumps to TIFF.
 *
 * @param name Path to the TIFF file to create or overwrite
 * @param im Dense single-dimensional array of data to store
 *
 * @return 0 on success
 */
int
main(int argc, char* argv[])
{
  if (argc == 2 &&
      strncmp(argv[1], "--test", 6)) {
    uint32_t im_a[SZ * SZ];
    float im_c[SZ * SZ];
    for (uint32_t i = 0; i < SZ * SZ; ++i) {
      im_a[i] = i;
      im_c[i] = ((float) i) / (SZ * SZ);
    }
    write_im_agent("a.tif", im_a);
    write_im_chemokine("c.tif", im_c);
  }
  
  return 0;
}
