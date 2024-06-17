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
#include <getopt.h>
#include <unistd.h>
#include <errno.h>

#include <tiffio.h>
#include <omp.h>


/** Change this to your GranSim grid size. */
const int SZ = 300;

/** Change this to your GranSim timesteps-per-day. */
const int TS_PER_DAY = 144;

typedef enum {
  AGENT,
  CHEMOKINE
} im_t;

/** Store program options. */
struct opts_s {
  char* i;
  char* o;
  char* e;
  char* d;
} opts_default = {".", "img", "0", "0"};
typedef struct opts_s opts_t;

/** Program usage. */
const char USAGE[] = "Usage: gransimg [--help] [--test] [--input DIR_INPUT] "
  "[--output DIR_OUTPUT]\n"
  "                [--exps LIST_EXPS] [--days LIST_DAYS]\n\n"
  "Convert GranSim agent and grid dumps to TIFF images.  Images from each\n"
  "experiment are created in parallel using OpenMP.\n\n"
  "options:\n"
  "--help    Show this message and exit.\n"
  "--test    Only generate two test images: agent.tif, chemokine.tif without\n"
  "          using any parallelism.\n"
  "--input   Path to the model run directory (default: .)\n"
  "--output  Path to the TIFF image output directory (default: img).\n"
  "--exps    String of integers separated by non-numeric characters without\n"
  "          corresponding replicate integers (default: 0 for all).  All\n"
  "          replicates are run and individual replicates are not currently\n"
             "supported.\n"
  "--days    String of integers separated by non-numeric characters\n"
  "          (default: 0 for the last day).\n";

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
 * @internal
 * @brief Parse CLI arguments.
 *
 * @param argc Passed from main().
 * @param argv Passed from main().
 */
void
parse_args(opts_t* opts, int argc, char* argv[])
{
  /* Populate the option structs described in `man 3 getopt_long`. */
  struct option longopts[] = {
    {"help", no_argument, NULL, 'h'},
    {"test-args", no_argument, NULL, 'a'},
    {"test-imgs", no_argument, NULL, 't'},
    {"input", required_argument, NULL, 'i'},
    {"output", required_argument, NULL, 'o'},
    {"exps", required_argument, NULL, 'e'},
    {"days", required_argument, NULL, 'd'},
    {NULL, 0, NULL, 0}
  };

  int test_args = 0;
  int ch;
  while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
    switch(ch) {
    case 'h':
      printf(USAGE);
      exit(0);
    case 't':
      {
	uint32_t im_a[SZ * SZ];
	float im_c[SZ * SZ];
	for (uint32_t i = 0; i < SZ * SZ; ++i) {
	  im_a[i] = i;
	  im_c[i] = ((float) i) / (SZ * SZ);
	}
	write_im_agent("agent.tif", im_a);
	write_im_chemokine("chemokine.tif", im_c);
	exit(0);
      }
    case 'a':
      test_args = 1;
      break;
    case 'i':
      opts->i = optarg;
      break;
    case 'o':
      opts->o = optarg;
      break;
    case 'e':
      opts->e = optarg;
      break;
    case 'd':
      opts->d = optarg;
      break;
    default:
      printf(USAGE);
      exit(1);
    }
  }

  if (test_args) {
    printf("input: %s\noutput: %s\nexps: %s\ndays: %s\n",
	   opts->i, opts->o, opts->e, opts->d);
    exit(0);
  }
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
  opts_t opts = opts_default;
  parse_args(&opts, argc, argv);

  return 0;
}
