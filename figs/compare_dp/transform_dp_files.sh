#!/usr/bin/env bash

sed '/Random/d' $1 | sed '/PARAMETER/d' | sed '/nIterations/d' | sed 's/://g' > $1_out

mv $1_out $1
