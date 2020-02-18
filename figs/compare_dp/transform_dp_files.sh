#!/usr/bin/env bash

sed '/PARAMETER/d' $1 | sed '/nIterations/d' | sed 's/://g' > $1_out

mv $1_out $1
