#!/bin/sh
sqlite3 << EOF
create table data(a,b,c);
.separator "\t"
.import $1 data
select a, group_concat(b), group_concat(c) from data group by a;
EOF
