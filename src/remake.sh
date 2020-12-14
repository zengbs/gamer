#!/bin/bash
make clean
make -f $1 -j $2
