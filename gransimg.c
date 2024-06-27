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

/* C standard library */
#include <ctype.h>
#include <errno.h>
#include <regex.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

/* Unix library */
#include <dirent.h>
#include <getopt.h>
#include <sys/syslimits.h>	/* PATH_MAX */
#include <unistd.h>

/* External libraries */
#include <omp.h>
#include <tiffio.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>


/* BEGIN USER CONFIGURATION */

/** Change this to your GranSim grid size. */
const int SZ = 300;

/** Change this to your GranSim timesteps-per-day. */
const int TS_PER_DAY = 144;

/** The value of `--agent-dump-file` in "lung-model-options-short.sh". */
const char FILE_AGENTS[] = "agents";

/** The value(s) of `--grids-to-dump` in "lung-model-options-short.sh". */
const size_t NCHARS_FILES_GRIDS = 10;
const char FILES_GRIDS[][NCHARS_FILES_GRIDS] = {
  "TNF",
  "nKillings"
};
const size_t N_FILES_GRIDS = sizeof(FILES_GRIDS) / NCHARS_FILES_GRIDS;

/** Mapping of agent states to groups. */
typedef enum {
  MAC,
  T_GAM,
  T_CYT,
  T_REG,
  FIB,
  MYOFIB,
  N_AGENTS,
  IGNORE
} agent_group_t;
const size_t NCHARS_AGENT_GROUP = 7;
const char AGENT_GROUP_NAMES[][NCHARS_AGENT_GROUP] = {
  "mac",
  "t_gam",
  "t_cyt",
  "t_reg",
  "fib",
  "myofib"
};
/* O(1) lookup of agent group from agent TYPE and STATE. */
const int N_AGENT_TYPES = 7;
const int N_AGENT_STATES = 4;
const int AGENT_GROUPS[N_AGENT_TYPES][N_AGENT_STATES] = {
  /*    0,      1,      2,      3 : STATES */
  {   MAC,    MAC,    MAC,    MAC}, /* 0 : MAC */
  {IGNORE, IGNORE, IGNORE, IGNORE}, /* 1 : DCELL_LUNG */
  {IGNORE,  T_GAM,  T_GAM, IGNORE}, /* 2 : TGAM */
  {IGNORE,  T_CYT,  T_CYT, IGNORE}, /* 3 : TCYT */
  {IGNORE,  T_REG,  T_REG, IGNORE}, /* 4 : TREG */
  {IGNORE, IGNORE, IGNORE, IGNORE}, /* 5 : MTB */
  {   FIB, MYOFIB, IGNORE, IGNORE}, /* 6 : FIBROBLAST */
};

/* END USER CONFIGURATION */


/* BEGIN CONSTANTS, ENUMS, AND STRUCTS */

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

/** Program usage matching LONGOPTS below. */
const char USAGE[] = "Usage: gransimg [--help] [--test] [--input DIR_INPUT] "
  "[--output DIR_OUTPUT]\n"
  "                [--exps LIST_EXPS] [--days LIST_DAYS]\n\n"
  "Convert GranSim agent and grid dumps to TIFF images.  Images from each\n"
  "experiment are created in parallel using OpenMP.\n\n"
  "options:\n"
  "--help       Show this message and exit.\n"
  "--test-imgs  Only generate two test images: agent.tif, chemokine.tif without\n"
  "             using any parallelism.\n"
  "--test-args  Show input arguments.\n"
  "--input      Path to the model run directory (default: .)\n"
  "--output     Path to the TIFF image output directory (default: img).\n"
  "--exps       String of integers separated by non-numeric characters without\n"
  "             corresponding replicate integers (default: 0 for all).  All\n"
  "             replicates are run and individual replicates are not currently\n"
  "             supported.\n"
  "--days       String of integers separated by non-numeric characters\n"
  "             (default: 0 for the last day).\n";

/** Populate the option structs described in `man 3 getopt_long`. */
const struct option LONGOPTS[] = {
  {"help", no_argument, NULL, 'h'},
  {"test-args", no_argument, NULL, 'a'},
  {"test-imgs", no_argument, NULL, 't'},
  {"input", required_argument, NULL, 'i'},
  {"output", required_argument, NULL, 'o'},
  {"exps", required_argument, NULL, 'e'},
  {"days", required_argument, NULL, 'd'},
  {NULL, 0, NULL, 0}
};

/** Metadata to query images from the model run directory. */
struct run_s {
  char* root;
  int* exps;
  int* times;
  int n_exps;
  int n_times;
} run_default = {NULL, NULL, NULL, 0, 0};
typedef struct run_s run_t;
const int EXPS_MAX = 65535;
const int TIMES_MAX = 65535;

/* END CONSTANTS, ENUMS, AND STRUCTS */


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
  fprintf(stderr, "TIFFOpen(\"...\", \"w\") = %p\n", tif);
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
 * @brief Print semi-parsed CLI arguments.
 *
 * @param opts Parsed options.
 */
void
print_args(opts_t* opts)
{
  printf("input: %s\noutput: %s\nexps: %s\ndays: %s\n",
	 opts->i, opts->o, opts->e, opts->d);
}


/**
 * @internal
 * @brief Semi-parse CLI arguments.
 *
 * @param opts Stores parsed options.
 * @param argc Passed from main().
 * @param argv Passed from main().
 */
void
parse_args(opts_t* opts, int argc, char* argv[])
{
  int test_args = 0;
  int ch;
  while ((ch = getopt_long(argc, argv, "", LONGOPTS, NULL)) != -1) {
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
    print_args(opts);
    exit(0);
  }
}


/**
 * @internal
 * @brief Read the model run directory seed.
 *
 * @param seed Output seed integer value.
 * @param path Output path to the seed file.
 * @param exp Experiment to read.
 * @param run Simulation run object containing the root directory.
 * @param buffer Temporary storage of seed file contents.
 */
void
run_read_seed(int* seed, char* path, int exp, run_t* run, char* buffer,
	      size_t buffer_sz)
{
  sprintf(path, "%s/exp%d/exp%d-1/seed", run->root, exp, exp);
  FILE* stream = fopen(path, "r");
  if (stream) {
    fgets(buffer, buffer_sz, stream);
    fclose(stream);
  } else {
    fprintf(stderr, "Error: Cannot read file %s\n", path);
    exit(1);
  }
  *seed = atoi(buffer);
}


/**
 * @brief Initialize metadata for the model directory.
 *
 * @param run Stores new run object.
 * @param opts Semi-parsed options from parse_opts().
 */
void
run_init(run_t* run, opts_t* opts)
{
  char* root = (char*) malloc(PATH_MAX * sizeof(char));
  run->root = realpath(opts->i, root);

  /* If no experiments are provided, process all experiment directories. */
  if (strncmp(opts->e, "0", 1) == 0) {
    /* Allocate maximum size and realloc after finding the paths. */
    regex_t preg;
    const char pat[] = "^exp([0-9]+)$";
    int ret = regcomp(&preg, pat, REG_EXTENDED);
    if (ret != 0) {
      fprintf(stderr,
	      "Error: Failed to compile regular expression: %s\n",
	      run->root);
      exit(1);
    }
    DIR* dir = opendir(run->root);
    if (dir == NULL) {
      fprintf(stderr,
	      "Error: Failed to open directory: %s\n",
	      run->root);
      exit(1);
    }
    int* exps = (int*) malloc(EXPS_MAX * sizeof(int));
    struct dirent* ent;
    while ((ent = readdir(dir)) != NULL) {
      if (ent->d_type == DT_DIR) {
#ifdef DEBUG
	fprintf(stderr,
		"Checking if directory matches %s pattern: %s ...",
		pat, ent->d_name);
#endif
	size_t nmatch = 2;
	regmatch_t pmatch[nmatch];
	ret = regexec(&preg, ent->d_name, nmatch, pmatch, 0);
	if (ret == 0) {
	  char match_c[5];
	  strncpy(match_c,
		  ent->d_name + pmatch[1].rm_so,
		  pmatch[1].rm_eo - pmatch[1].rm_so);
	  exps[run->n_exps] = atoi(match_c);
#ifdef DEBUG
	  fprintf(stderr, " %d\n", exps[run->n_exps]);
#endif
	  run->n_exps++;
	} else if (ret == REG_NOMATCH) {
#ifdef DEBUG
	  fprintf(stderr, " no\n");
#endif
	} else {
	  char errbuf[1024];
	  regerror(ret, &preg, errbuf, 1024);
	  fprintf(stderr, "Error: Regex match failed: %s\n", errbuf);
	  exit(1);
	}
      }
    }
    regfree(&preg);
    exps = (int*) realloc(exps, sizeof(int) * run->n_exps);
    run->exps = exps;
  } else {
    /* TODO. */
  }

  /* If no days are provided, process last timepoint. */
  int* times = (int*) malloc(TIMES_MAX * sizeof(int));
  if (strncmp(opts->d, "0", 1) == 0) {
    run->n_times = 1;
    /* The smallest file is the stats file. */
    int exp = run->exps[0];
    char file_seed[PATH_MAX];
    int seed = -1;
    char buffer[65535];
    run_read_seed(&seed, file_seed, exp, run, buffer, 65535);
    sprintf(file_seed, "%s/exp%d/exp%d-1/seed%d.csv",
	    run->root, exp, exp, seed);
    FILE* stream = fopen(file_seed, "r");
    if (stream) {
      /* Read the last line. */
      while (! feof(stream)) {
	fscanf(stream, "%[^\n]\n", buffer);
      }
      fclose(stream);
#ifdef DEBUG
      fprintf(stderr, "Last line: %s\n", buffer);
#endif
      for (size_t i = 0; i < sizeof(buffer); ++i) {
	if (buffer[i] == ',') {
	  buffer[i] = '\0';
	  break;
	}
      }
      times[0] = atoi(buffer);
    } else {
      fprintf(stderr, "Error: Cannot read file %s\n", file_seed);
      exit(1);
    }
    times = (int*) realloc(times, sizeof(int) * run->n_times);
  } else {
    /* TODO. */
  }
  run->times = times;
}


/**
 * @brief Print metadata for the model directory.
 *
 * @param run Run object to inspect.
 */
void
run_print(run_t* run)
{
  printf("root: %s\nexps: <%d> ", run->root, run->n_exps);
  printf("%d", run->exps[0]);
  for (int i = 1; i < run->n_exps; ++i) {
    printf(", %d", run->exps[i]);
  }
  printf("\ntimes: <%d> ", run->n_times);
  printf("%d", run->times[0]);
  for (int i = 1; i < run->n_times; ++i) {
    printf(", %d", run->times[i]);
  }
  printf("\n");
}


/**
 * @brief Free allocated memory of model directory metadata.
 *
 * @param run Run object to free.
 */
void
run_destroy(run_t* run)
{
  free(run->root);
  free(run->exps);
  free(run->times);
}


/**
 * @brief Read GranSim agents dump.
 *
 * @param im Output to fill from path.
 * @param path Agent dump CSV file to read.
 */
int
read_agents(uint32_t* im, char* path, int time, char* buffer, size_t buffer_sz)
{
  FILE* stream = fopen(path, "r");
  char time_c[128];
  if (stream) {
    /* Skip the header. */
    fscanf(stream, "%[^\n]\n", buffer);
    /* Find the line containing the timepoint. */
    bool found = false;
    while (! feof(stream) && ! found) {
      fscanf(stream, "%[^\n]\n", buffer);
      for (size_t i = 0; i < buffer_sz; ++i) {
	if (buffer[i] == ',') {
	  strncpy(time_c, buffer, i);
	  if (time == atoi(time_c)) {
	    found = true;
	  }
	  break;		/* break the for loop. */
	}
      }
    }
#ifdef DEBUG
    fprintf(stderr, "Image beginning line: %s\n", buffer);
    int counter[N_AGENTS] = {0};
#endif
    /* Scan the image values. */
    int line_time = time;
    while (! feof(stream) && time == line_time) {
      int id, type, state, x, y;
      fscanf(stream, "%d,%d,%d,%d,%d,%d,%*[^\n]\n",
	     &line_time, &id, &type, &state, &x, &y);
      if (time == line_time &&
	type < N_AGENT_TYPES &&
	  state < N_AGENT_STATES) {
	int agent = AGENT_GROUPS[type][state];
	if (agent != IGNORE) {
#ifdef DEBUG
	  ++counter[agent];
#endif
	  int offset = agent * SZ * SZ;
	  int row = y;
	  int col = x;
	  im[(row * SZ + col) + offset] = id;
	}
      }
    }
#ifdef DEBUG
    /* Image agent counts. */
    fprintf(stderr, "  Agents: %s = %d", AGENT_GROUP_NAMES[0], counter[0]);
    for (int i = 1; i < N_AGENTS; ++i) {
      fprintf(stderr, ", %s = %d", AGENT_GROUP_NAMES[i], counter[i]);
    }
    fprintf(stderr, "\n");
#endif
  }
  fclose(stream);
  return 0;
}


/**
 * @brief Read GranSim grid dump.
 *
 * @param im Output to fill from path.
 * @param path Agent dump CSV file to read.
 */
int
read_grid(float* im, char* path, int time, char* buffer, size_t buffer_sz)
{
  FILE* stream = fopen(path, "r");
  char time_c[128];
  size_t offset = -1;
  if (stream) {
    /* Find the line containing the timepoint. */
    bool found = false;
    while (! feof(stream) && ! found) {
      fscanf(stream, "%[^\n]\n", buffer);
      for (size_t i = 0; i < buffer_sz; ++i) {
	if (buffer[i] == ',') {
	  strncpy(time_c, buffer, i);
	  if (time == atoi(time_c)) {
	    found = true;
	    offset = i;
#ifdef DEBUG
	    fprintf(stderr, "Image offset: %lu\n", i);
#endif
	  }
	  break;		/* break the for loop. */
	}
      }
    }
#ifdef DEBUG
    /* Printing the entire buffer line is massive, so only print a preview. */
    strncpy(time_c, buffer, 127);
    fprintf(stderr, "Image line: %s...\n", time_c);
#endif
    /* Scan the image values. */
    char* ptr_begin = &buffer[offset + 1];
    char* ptr_end = NULL;
    for (int i = 0; i < SZ; ++i) {
      for (int j = 0; j < SZ; ++j) {
	ptr_end = strchr(ptr_begin, ',');
	im[i * SZ + j] = strtof(ptr_begin, &ptr_end);
	ptr_begin = ptr_end + 1;
      }
    }
#ifdef DEBUG
    /* Preview image. */
    fprintf(stderr, "Image: %f", im[0]);
    for (int i = 1; i < 63; ++i) {
      fprintf(stderr, ", %e", im[i]);
    }
    fprintf(stderr, ", ...");
    for (int i = SZ * SZ - 64; i < SZ * SZ; ++i) {
      fprintf(stderr, ", %e", im[i]);
    }
    fprintf(stderr, "\n");
#endif
  }
  fclose(stream);
  return 0;
}


/**
 * @brief Read XML XPath.
 *
 * @param value Output string corresponding to xpath.
 * @param sz_value Output string size.
 * @param xpath Input full XPath to retrieve.
 * @param path Path to XML file.
 *
 * @return 0 if XML document parsed successfully, otherwise 1.
 */
int
read_xml_xpath(char* value, size_t sz_value, char* xpath, char* path)
{
  /* Validate XPath argument. */
  char* ptr_end;
  ptr_end = strchr(xpath, '@');
  if (ptr_end == NULL) {
    fprintf(stderr, "Error: XPath missing '@' attribute: %s\n", xpath);
    return 1;
  }
  char xpath_node[strlen(xpath) + 1];
  strncpy(xpath_node, xpath, ptr_end - xpath);
  char* xpath_attr = (ptr_end + 1);
#ifdef DEBUG
  fprintf(stderr, "XPath %s split into %s node and %s attribute.\n",
	  xpath, xpath_node, xpath_attr);
#endif

  /* Retrieve XPath node attribute. */
  xmlDocPtr doc = xmlParseFile(path);
  if (doc == NULL) {
    fprintf(stderr, "Error: XML document not parsed successfully: %s\n", path);
    return 1;
  }
  xmlNodePtr cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr, "Error: XML document is empty: %s\n", path);
    xmlFreeDoc(doc);
    return 1;
  }
  if (xmlStrcmp(cur->name, (const xmlChar *) "GR")) {
    fprintf(stderr, "Error: XML document missing root node GR: %s\n", path);
    xmlFreeDoc(doc);
    return 1;
  }
  xmlXPathContextPtr context = xmlXPathNewContext(doc);
  xmlXPathObjectPtr result =
    xmlXPathEvalExpression((const xmlChar*) xpath_node, context);
  if (xmlXPathNodeSetIsEmpty(result->nodesetval)) {
    fprintf(stderr, "Error: XML document does not contain %s XPath node: %s\n",
	    xpath, xpath_node);
    xmlXPathFreeObject(result);
    xmlXPathFreeContext(context);
    xmlFreeDoc(doc);
    return 1;
  }
  strncpy(value,
	  (char*) xmlGetProp(*result->nodesetval->nodeTab,
			     (xmlChar*) xpath_attr),
	  sz_value);
  xmlXPathFreeObject(result);
  xmlXPathFreeContext(context);
  xmlFreeDoc(doc);
  return 0;
}


/**
 * @brief Convolve probability image and kernel ignoring edge values.
 *
 * @param res Output probability image.
 * @param im Image to convolve with the kernel.
 * @param kernel Kernel to convolve with im.
 * @param sz_kernel One of the dimensions of the square kernel.
 */
void
convolve_prob(float* res, uint32_t* im, float* kernel, size_t sz_kernel)
{
  /* Initialize result image. */
  for (int i = 0; i < SZ * SZ; ++i) {
    res[i] = 0.0;
  }
  int kern_c = sz_kernel / 2;
  for (int i = 0; i < SZ; ++i) {
    for (int j = 0; j < SZ; ++j) {
      /* Convolution. */
      for (size_t m = 0; m < sz_kernel; ++m) {
	int mm = sz_kernel - 1 - m;
	for (size_t n = 0; n < sz_kernel; ++n) {
	  int nn = sz_kernel - 1 - n;
	  int ii = i + (kern_c - mm);
	  int jj = j + (kern_c - nn);
	  if (ii >= 0 && ii <= SZ &&
	      jj >= 0 && jj <= SZ) {
	    res[i * SZ + j] += im[ii * SZ + jj] * kernel[mm * sz_kernel + nn];
	  }
	}
      }
      /* Maximum probability of 1.0. */
      if (res[i * SZ + j] > 1.0) {
	res[i * SZ + j] = 1.0;
      }
    }
  }
}


/**
 * @brief Write all images for a given experiment and time.
 *
 * @param exp Integer of experiment.
 * @param time Integer of time ticks.
 * @param run Initialized model run_t.
 * @param im_i Uninitialized 1-D uint32_t SZ * SZ * N_AGENTS image array.
 * @param im_f Uninitialized 1-D float SZ * SZ image array.
 * @param buffer Uninitialized char buffer large enough to read grid dump line.
 * @param buffer_sz Size of buffer.
 *
 * @return 0 on success, 1 if any image is missing.
 */
int
write_ims(int exp, int time, run_t* run, opts_t* opts, uint32_t* im_i,
	  float* im_f, char* buffer, size_t buffer_sz)
{
  int im_counter = 1;
#ifdef DEBUG
  printf("Image buffer size: %lu\n", buffer_sz);
#endif
  /* Read the seed. */
  char path[PATH_MAX];
  int seed;
  run_read_seed(&seed, path, exp, run, buffer, buffer_sz);

  /* Write the agents. */
  for (size_t i = 0; i < SZ * SZ * N_AGENTS; ++i) {
    im_i[i] = 0;
  }
  sprintf(path, "%s/exp%d/exp%d-1/%s%d.csv",
	  run->root, exp, exp, FILE_AGENTS, seed);
  fprintf(stderr, "Processing exp %d time %d agents from %s\n",
	  exp, time, path);
#ifdef DEBUG
  fprintf(stderr, "  read_agents(%p, \"%s\")\n",
	  im_i, path);
#endif
  read_agents(im_i, path, time, buffer, buffer_sz);
  char im_name[128];
  for (int i = 0; i < N_AGENTS; ++i) {
    strncpy(im_name, AGENT_GROUP_NAMES[i], NCHARS_AGENT_GROUP);
#ifdef DEBUG
    fprintf(stderr, "  write_im_agent(\"%s\")\n", im_name);
#endif
    sprintf(path, "%s/exp%d_time%d_%02d_%s.tif",
	    opts->o, exp, time, im_counter, im_name);
    write_im_agent(path, &im_i[i * SZ * SZ]);
    ++im_counter;
  }

  /* Write agent-derived grids. */
  char path_xml[PATH_MAX];
  sprintf(path_xml, "%s/%d.xml", run->root, exp);
  char xpath_prob[] = "/GR/Core/Tcell/Tgam@probIFNMooreExtend";
  char prob_c[10];
  int ret = read_xml_xpath(prob_c, 10, xpath_prob, path_xml);
  if (ret == 0) {
    float prob = atof(prob_c);
#ifdef DEBUG
    fprintf(stderr, "Setting exp %d IFN-g probability from %s = %e\n",
	    exp, xpath_prob, prob);
#endif
    float kernel[5 * 5] = {
      prob, prob, prob, prob, prob,
      prob,   1.,   1.,   1., prob,
      prob,   1.,   1.,   1., prob,
      prob,   1.,   1.,   1., prob,
      prob, prob, prob, prob, prob,
    };
    uint32_t* im_tgam = &im_i[T_GAM * SZ * SZ];
    convolve_prob(im_f, im_tgam, kernel, 5);
    strncpy(im_name, "ifng", sizeof("ifng"));
#ifdef DEBUG
    fprintf(stderr, "  write_im_grid(\"%s\")\n", im_name);
#endif
    sprintf(path, "%s/exp%d_time%d_%02d_%s.tif",
	    opts->o, exp, time, im_counter, im_name);
    write_im_chemokine(path, im_f);
  }
  ++im_counter;

  /* Write the grids. */
  for (size_t i = 0; i < N_FILES_GRIDS; ++i) {
    sprintf(path, "%s/exp%d/exp%d-1/%s%d.csv",
	    run->root, exp, exp, FILES_GRIDS[i], seed);
    fprintf(stderr, "Processing exp %d time %d grid %ld from %s\n",
	    exp, time, i, path);
    strncpy(im_name, FILES_GRIDS[i], sizeof(FILES_GRIDS[i]));
#ifdef DEBUG
    fprintf(stderr, "  read_grid(%p, \"%s\")\n",
	    im_f, path);
#endif
    read_grid(im_f, path, time, buffer, buffer_sz);
    if (strncmp(FILES_GRIDS[i], "nKillings", 9) == 0) {
      strncpy(im_name, "caseum", sizeof("caseum"));
      char threshold_c[3];
      char xpath_threshold[] = "/GR/Core@nrKillingsCaseation";
      ret = read_xml_xpath(threshold_c, 3, xpath_threshold, path_xml);
      if (ret == 0) {
	int threshold = atoi(threshold_c);
#ifdef DEBUG
	fprintf(stderr, "Setting exp %d caseum threshold from %s = %d\n",
		exp, xpath_threshold, threshold);
#endif
	for (int i = 0; i < SZ * SZ; ++i) {
	  im_f[i] = im_f[i] >= threshold;
	}
      }
    }
    for (unsigned long i = 0; i < strlen(im_name); ++i) {
      im_name[i] = tolower(im_name[i]);
    }
#ifdef DEBUG
    fprintf(stderr, "  write_im_chemokine(\"%s\")\n", im_name);
#endif
    sprintf(path, "%s/exp%d_time%d_%02d_%s.tif",
	    opts->o, exp, time, im_counter, im_name);
    write_im_chemokine(path, im_f);
    ++im_counter;
  }

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
  opts_t opts = opts_default;
  parse_args(&opts, argc, argv);
#ifdef DEBUG
  print_args(&opts);
#endif

  run_t run = run_default;
  run_init(&run, &opts);
  run_print(&run);

  /* Shared memory allocations. */
  const int max_sets = run.n_exps * run.n_times;
  const int ndigits_max_sets = 1 + max_sets / 10;
  atomic_int sum_sets = 0;
#pragma omp parallel
  {
    /* Thread-specific memory allocations:
       All agents. */
    uint32_t* im_i = NULL;
    im_i = (uint32_t*) malloc(sizeof(uint32_t) * SZ * SZ * N_AGENTS);
    if (im_i == NULL) {
      fprintf(stderr, "Error: Cannot allocate image agent buffer\n");
      exit(1);
    }
    /* Single grid dump. */
    float* im_f = NULL;
    im_f = (float*) malloc(sizeof(float) * SZ * SZ);
    if (im_f == NULL) {
      fprintf(stderr, "Error: Cannot allocate image grid buffer\n");
      exit(1);
    }
    /* Buffer large enough to hold an entire row of the grid dump CSV. */
    const size_t buffer_sz = 1<<24;
    char* buffer = NULL;
    buffer = (char*) malloc(sizeof(char) * buffer_sz);
    if (buffer == NULL) {
      fprintf(stderr, "Error: Cannot allocate CSV buffer\n");
      exit(1);
    }

    /* Main loop. */
#pragma omp for nowait
    for (int i = 0; i < max_sets; ++i) {
      int exp = run.exps[i - i / run.n_exps];
      int time = run.times[i / run.n_exps];
      write_ims(exp, time, &run, &opts, im_i, im_f, buffer, buffer_sz);
      ++sum_sets;
      fprintf(stderr,
	      "Progress: %*d / %*d (%5.1f%%) after finishing exp %d time %d\n",
	      ndigits_max_sets, sum_sets,
	      ndigits_max_sets, max_sets,
	      (100.0 * sum_sets) / max_sets,
	      exp, time);
    }

    free(im_i);
    free(im_f);
    free(buffer);
  }
  fprintf(stderr, "Finished\n");

  run_destroy(&run);

  return 0;
}
