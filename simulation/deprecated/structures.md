## `simulation`
>input: 
>    - simulation functions
>        - otu population
>            - dimension, correlation
>        - mediation -> target population
>
>output:
>    - simulator object
>        - input:
>            - truth values (`beta`s)
>        - output:
>            - simulated otu populations abundance and metabolite abundance

|`variable`|dimension|
|------|------|
|otu population abundance (mediation)|`mediations` * `n`|
|otu population abundance (not mediation)| `p_otu`-`mediations` * `n`|
|target population abundance|`1` * `n`|
|metabolite abundance (mediation)|`mediations` * `n`|
|metabolite abundance (not mediation)|`p_metabolite`-`mediations` * `n`|

## `tester`
>input:
>   - simulator object
>   - number of simulation trials `t`
>   - method of prediction
>       - B&K steps
>
>output:
>   - statistic testing result table
>       - predicted `beta`s
>       - p-values
>
output table format:

For `t`=`some number of trials`:
|`beta`|truth value|`p_value`|
|:------:|:------:|:------:|
|`b11`|-|-|
|`b21`|-|-|
|`b32`|-|-|