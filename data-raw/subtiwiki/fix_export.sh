
# Problem: new lines at end of descriptive text
# Observation: '$' sign does not seem to appear
# Approach: for sed to work substitute newlines
#           then replace the problamatic char
#           and then convert back
cat export.csv                   | \
    tr '\n' '$'                  | \
    sed -e 's/\([^"]\)\$"/\1"/g' | \
    tr '$' '\n'                  > export_fixed.csv 


# save some space
gzip *.csv
