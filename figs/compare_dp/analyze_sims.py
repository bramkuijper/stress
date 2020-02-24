#!/usr/bin/env python3

# summarize the individual-based model
import summarizesims
import argparse


parser = argparse.ArgumentParser(description="Process simulation files")

parser.add_argument("--pattern"
                    ,help="Regular expression pattern to match files"
                    ,default="(sim.*\d+$|iter.*\d+)"
                    )

parser.add_argument("--path"
                    ,help="The path over which to iterate"
                    ,default=".")


parser.add_argument("--sep"
                    ,help="Column separator"
                    ,default=";")

args = vars(parser.parse_args())

summarizesims.SummarizeSims(
    path=args["path"]
    ,sep=args["sep"]
    ,pattern=args["pattern"]
    )
