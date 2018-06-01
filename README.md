# SLINK_CLINK
Implementation of SLINK(Optimized Single Linkage Clustering), CLINK(Complete Linkkage Clustering) algorithm


## SLINK
SLINK: an optimally efficient algorithm for the single-link cluster method
The Computer Journal 16 (1973), No. 1, p. 30-34


## CLINK
An Efficient Algorithm for a Complete Link Method
The Computer Journal, Volume 20, Issue 4, 1 January 1977, Pages 364–366

## Note
If removing the comments in the program, the result will show the detail value changes in each step. 


## Usage (setting)
numData: the number of numeric sample data 

numDimension: the number of dimension

cutThreshold: Double.POSITIVE_INFINITY for the whole hierarchy

directory: directory path

fileName: file name

Note: the data file should not contain the column name

## Output(example)
level: 1, height:0.5, list=[4, 6]

level: the order of merged 

height: the height of merged objects

list: the list of object index when they are merged


