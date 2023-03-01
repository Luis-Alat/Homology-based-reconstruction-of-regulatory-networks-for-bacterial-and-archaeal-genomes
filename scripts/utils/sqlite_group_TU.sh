#!/bin/sh
sqlite3 << EOF
create table data(a,b,c,d,e,f,g,h,i);
.separator "\t"
.import $1 data
select group_concat(a), b, c, group_concat(d), group_concat(e), group_concat(f), group_concat(g), group_concat(h), group_concat(i) from data group by b,c;
EOF
