#!/bin/bash

cut -f 3-5 dasDiversity.tab | sort | uniq -c > dasDiversity.histo
