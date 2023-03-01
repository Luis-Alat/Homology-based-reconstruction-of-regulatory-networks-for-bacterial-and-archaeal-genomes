#!/bin/sh

NUMBER_COLUMNS=$(head $1 | perl -F'\t' -nae "$col_number= scalar @F")

sqlite3 << EOF
create table data(a,b,c,d,e,f,g);
.separator "\t"
.import $1 data
select group_concat(a), b, c, group_concat(d), group_concat(e), group_concat(f), group_concat(g) from data group by b,c;
EOF
