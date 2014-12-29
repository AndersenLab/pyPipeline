#!/usr/bin/env python
from ast import literal_eval
import sys
from utils import remove_file

files_to_delete = literal_eval(sys.argv[1])

for file in files_to_delete:
	remove_file(file)