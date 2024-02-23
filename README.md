# CATSequenceAnalysis

## Tests 
 Test ceases have 3 tetsCases sources, which can be edited and allows you to run the test with variety of sequence lengths. 
  - SequenceLengthProvider - is the source for  CanContributeForFormingTriangle and CollisionReportFor_1000
    - CanContributeForFormingTriangle - test if sequences with provided length can contribute for applying of trilateration method, i.e. calculated CAT profile can form a triangle which is essential for method of trilateration 
    - CollisionReportFor_1000 - generate report with number collisions which can occurs for all combination for a sequence with a particular length.
  - CompareCATwithNWProvider - is the source for CompareCATwithNW and CompareCATwithNWLong
    - CompareCATwithNW and CompareCATWithKMP generates report from CAT and Needleman–Wunsch or Knuth–Morris–Pratt comparison, where a subsequence with different length and from different places (from the first half, from the second half, from the middle ) is  taken from the input sequence. This tests also generates sequence with random length and makes the same comparison
  - CompareCATwithNWLongProvider - source for CompareCATwithNWLong, which is similar to CompareCATwithNW, but for each experiment (input source ) a new task is starter for a faster execution of the test

 - CombineCompareCATWithKMP - has to be ran after CompareCATWithKMP to combine all reports in a single csv for a easier processing of the results. 
 - CombineCompareCATWithNW - has to be ran after CompareCATwithNW to combine all reports in a single csv for a easier processing of the results.  

All reports are generated in the \SequenceAnaliseTests\bin\Debug\net6.0 folder 

 ## Reports' header legend 

- Dna1 Length - length of the 1-st DNA sequence
- Dna2 Length - length of the 1-st DNA sequence
- Kind - kind of teh second sequence, if it is a subsequence form the 1-st one and from which part is taken (from the first half, from the second half, from the middle) or it is random generated sequence 
- CAT Result - rate of similarity after comparison with CAT in [0-1]
- CAT ms - time in ms to perform the comparison with CAT
- NW Result	- rate of similarity after comparison with Needleman–Wunsch in [0-1]
- NW ms - time in ms to perform the comparison with Needleman–Wunsch
- KMP Result	- rate of similarity after comparison with Knuth–Morris–Pratt in [0-1]
- KMP ms - time in ms to perform the comparison with Knuth–Morris–Pratt

- dna string - raw data of the sequence 
- totalRecords - total number of all combination for a sequences with length given as input parameter 
- cat matches - total number of collisions 
- [0.1 - 1.0] - rate of similarity after comparison with Needleman–Wunsch for a sequence which is discovered as collision based on CAT

## How to run tests

 Open *sln file and run desired test from teh test explorer https://learn.microsoft.com/en-us/visualstudio/test/run-unit-tests-with-test-explorer?view=vs-2022. Keep in mind that some test takes longer (up to couple of hours) based on the machine and sequence length which is examined. This is the reason why test generates a reports in *.csv for storing the results. 

> **Warning:**
> For localizations like German where separator for csv is ';' and decimal separator is ',' instead of '.' opening the reports in Excel could appear weird  - https://answers.microsoft.com/en-us/msoffice/forum/all/power-query-how-to-import-a-csv-file-that-does-not/e19fe942-80ad-4ac1-a164-061e7cac9486      

## Block diagram  of CAT method for DNA sequence analysis

```plantuml
@startuml
skinparam classFontSize 25
skinparam classFontColor red

start
partition "Local Variables"{
:Benchmark.minPoint = 0.4d
Benchmark.maxPoint = 1d
baseDistance = new double[4] { 0d, 0.6d, Benchmark.minPoint, 0.6d } 
nearMatches = 0d
exactMatches = 0d 
prevMatches = 0d  
bonusTotal = 0d ;
}
partition "Loop"{
while (index < dnaStr.Length) is (yes)
-[#green]->
: @base = dnaStr[index];
:var nearMatch = 0d 
exactMatch = 0d;
partition "**function** Near Matches" {
    if(Benchmark[(index + 0) % Benchmark.Length] == @base) then (yes)
    -[#green]->
        :nearMatch += (baseDistance[0] + prevMatches * baseDistance[0]);
    endif 
    -[#green]->
    if(Benchmark[(index + 1) % Benchmark.Length] == @base) then (yes)
    -[#green]->
        :nearMatch += (baseDistance[1] + prevMatches * baseDistance[1]);
    endif
    -[#green]->
    if(Benchmark[(index + 2) % Benchmark.Length] == @base) then (yes)
    -[#green]->
        :nearMatch += (baseDistance[2] + prevMatches * baseDistance[2]);
    endif
    -[#green]->
    if(Benchmark[(index + 3) % Benchmark.Length] == @base) then (yes)
    -[#green]->
        :nearMatch += (baseDistance[3] + prevMatches * baseDistance[3]);
    endif
    -[#green]->
}

partition "**function** Exact Match" {

    if(Benchmark[index % Benchmark.Length] == @base) then (yes)
    -[#green]->
        :exactMatch = Benchmark.maxPoint;
    endif
}

partition "Dependency to previous matches and index" {
:nearMatches += nearMatch + prevMatches * nearMatch 
exactMatches += exactMatch + prevMatches * exactMatch 
prevMatches = sequencePrevLength / (i + 1) + exactMatch + nearMatch - Benchmark.minPoint
bonusTotal += i != dnaString.Length - 1 ? 0 : prevMatches
sequencePrevLength +=  prevMatches;
}

backward: index++;

endwhile (no)
-[#red]->
}
:**dnaDistance** = (nearMatches + exactMatches) / (double)(bonusTotal + dnaString.Length)
**cos** = (Math.Pow(dnaDistance, 2) + Math.Pow(benchmarkDistance, 2) - Math.Pow(profile.dnaDistance, 2)) / (2 * dnaDistance * benchmarkDistance)
**d** = dnaDistance * cos
**h** = Math.Sqrt(Math.Pow(dnaDistance, 2) - Math.Pow(dnaDistance * benchmarkCOS, 2));
stop
@enduml
```
