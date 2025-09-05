#!/bin/sh

set -euxo pipefail

gcc -O3 -Wall -Wextra -pedantic \
    -I. -I/opt/homebrew/opt/openblas/include -I/usr/local/lib \
    -L/opt/homebrew/opt/openblas/lib -L/usr/local/lib \
    -lopenblas -Wl,-rpath,/opt/homebrew/opt/openblas/lib \
    -lnlopt -lm -Wl,-rpath,/usr/local/lib \
    elt.c \
    -o elt

gcc -O3 -Wall -Wextra -pedantic \
    -I. -I/opt/homebrew/opt/openblas/include -I/usr/local/lib \
    -L/opt/homebrew/opt/openblas/lib -L/usr/local/lib \
    -lopenblas -Wl,-rpath,/opt/homebrew/opt/openblas/lib \
    -lnlopt -lm -Wl,-rpath,/usr/local/lib \
    esa.c \
    -o esa
