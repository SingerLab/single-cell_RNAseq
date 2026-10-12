BEGIN { OFS = "\t" }
/^#/ { print; next }
{ if ($3 == ".") { $3 = "hg38_" $1 "_" $2 }; print }
