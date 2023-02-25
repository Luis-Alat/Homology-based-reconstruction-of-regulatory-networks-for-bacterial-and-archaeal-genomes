#!/bin/sh
sqlite3 << EOF
create table data(a,b,c,d);
.separator "\t"
.import $1 data
select c,d, group_concat(a) , group_concat(b) from data group by c,d;
EOF
