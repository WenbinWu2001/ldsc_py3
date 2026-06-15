**Use an example to illustrate each functionality module.**

1. What does this functionality do? -- use more plain language (what are the inputs and outputs, one concise sentence how do the estimated quantities facilitate our genetic understanding)

2. Required input formats (if input depends on the output of another functionality module, write out specifically what commands are used to create these input, or link to the wiki of that functionality, what steps are needed before running this module (e.g., munge-sumstats)).

3. minimal command to use (assume most users would use this)

   [Also give practical advice on memory usage and running time]

4. What are the outputs -- explain each artifact separately. Start from the main artifact storing the results (i.e.,. the one users practically use), describe what each row / column that is needed. Emphasize which quantities are most relevant to the current analysis (i.e., in standard, which columns or values should I use for interpretation and scienctific conclusions?)

5. Caveats (what to do and what not to do)

6. Advanced setups on command usage
   1. List all arguments in a table, like that in `io-argument-inventory.md`. Make sure you include concise description on the role of each argument and when each accepted value is needed, or what format is supported.
   2. If you want to ... [requirements], you can [change value of which argument], and the results will be [in what output format].

7. More details on output results
   1. Explain what each column means, how to interpret them, what scientific implication do they have, what interpretations we can and cannot obtain (especially for partitioned h2).

8. More details on the machinery under the hood
   1. explain more details in the pipeline -- e.g., how are ref panel handled (merging, matching, etc.), how can ld scores computed (concise). Can use math.
   2. Or consider relegate this to a separate design doc.