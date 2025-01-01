# Outcomes and Derived Variable Tables

Chirag Patel (chirag\@hms.harvard.edu)

December 30 2024

## This folder contains .rmds to define new variables in the nhanes database, such as new outcomes.

### General procedure the .rmd should follow:

1.  Connect to the database, e.g., `../db/nhanes_***.sqlite`

2.  Select directories per survey year to define new variables

3.  Define new variables in a new table

4.  Add new variable names to `variable_names_epcf`

5.  Add new table names to `table_names_epcf`
