# CHUI-PSO

This is a repository for CHUI-PSO, a _particle swarm optimization_ algorithm for _closed high-utility itemset mining_ in quantitative transactional data.

How to run:

* Select minimum utility threshold, population size, and number of iterations in CHUI_PSO.java.
* Set the "input" string in CHUI_PSO.java to the path of the dataset on file. The file must be in SPMF format, [here are some example datasets.](https://www.philippe-fournier-viger.com/spmf/index.php?link=datasets.php) 
* Set the "output" string to any .txt file path. The discovered patterns are written to this file during execution.
* Run main.java

The algorithm can also discover high-utility itemsets by changing the variable "closed" to False in CHUI_PSO.java
