/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "affinity.h"
#include "allocate.h"
#include "timing.h"

int main(int argc, char **argv) {
  if (argc == 1) {
    argv[1] = "";
  }

  char hostname[HOST_NAME_MAX];
  gethostname(hostname, sizeof(hostname));

  printf("hello, %s, on host %s\n", argv[1], hostname);

  return EXIT_SUCCESS;
}
