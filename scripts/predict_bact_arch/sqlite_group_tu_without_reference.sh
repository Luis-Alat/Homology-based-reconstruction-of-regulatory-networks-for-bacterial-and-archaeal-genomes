#!/bin/sh
sqlite3 << EOF
create table data(a,b,c,d,e,f,g);
.separator "\t"
.import $1 data
select a, b, group_concat(c), group_concat(d), group_concat(e), group_concat(f), group_concat(g) from data group by a,b;
EOF
