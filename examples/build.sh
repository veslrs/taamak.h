#!/usr/bin/env sh

set -xe

if [ "$1" = "-D" ]
then 
    CFLAGS="-Wall -Wextra -O0 -g -fsanitize=address,undefined"
elif [ $1 = "-R" ]
then
    CFLAGS="-Wall -Wextra -O2"
else
    echo "[ERROR] 1st argument: -D for debug build; -R for release build."
    exit -1
fi

gcc $CFLAGS \
    `pkg-config --cflags openblas` \
    -o $2 \
    $2.c \
    -L$HOME/opt/openblas/lib \
    -l:libopenblas.a -lm -lpthread
