cd $1

sqlite3 << EOF
CREATE TABLE element (id, elements); 
CREATE TABLE strand  (id, strand);

.separator "\t"
.import operon_group_elements element
.import operon_strand_type strand

SELECT e.id, e.elements, s.strand FROM element as e
    LEFT JOIN strand as s ON s.id=e.id;

EOF