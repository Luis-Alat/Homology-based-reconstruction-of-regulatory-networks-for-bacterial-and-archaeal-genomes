#!/bin/bash

# This script is expecting four arguments, in order
# 1. Names of the columns to create a table
# 2. Col names to execute the GROUP BY instruction
# 3. Instruction to execute inside the SELECT command according to the colunm names
# 4. File containing the table values

COL_NAMES=$1
GROUP_BY=$2
SELECT_INSTRUCTION=$3
TABLE=$4

sqlite3 << EOF

CREATE TABLE data ($COL_NAMES);
.separator "\t"
.import $TABLE data
SELECT $SELECT_INSTRUCTION FROM data GROUP BY $GROUP_BY;

EOF
