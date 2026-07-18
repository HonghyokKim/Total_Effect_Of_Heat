**Update: July 17, 2026**

The previous code did not correctly implement the intended estimands. The updated code now supports two analyses:

1. **Total effect of high ambient temperature (HAT):** for example, comparing temperature ≥25°C with temperature fixed at 25°C or below. In this analysis, the direct-effect component of HAT is represented by the direct effect of heat waves (HW) like the paper.
2. **Total effect of heat waves:** comparing HW = 1 with HW = 0.

One simulation dataset is included. With the second and third scripts, users can estimate either effect. Users may also generate their own simulated datasets using the first script.

Please email me if you identify any errors or have questions.



This repository provides example R code and simulated data for the paper:
“Total Effect of Heat on Mortality Considering Heat-Mediated Air Pollution and Interaction Effects under Demographic, Mitigation, and Adaptation Scenarios,” 
(Environmental Research, 2026)

#1 Code creates a simulated dataset.
#2 Code fits a regression model
#3 Code estimates the total effect and its decomposition.

