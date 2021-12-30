## Internal Model Project. Actuarial statistics

In this Project, an insurance company internal model is constructed considering life and non-life insurances. The objective of the model is to properly evaluate the premium risk (the risk of having more claims than expected) so the economic capital can be estimated.

<p align="center">
   <img src="https://raw.githubusercontent.com/asiergs/csv_to_view/main/Economic_Capital.svg" alt="2400"/>
</p>

The insurance company to be modelled is composed by the following number of policies:
-	Car insurances: 25,234
-	Life insurances: 20,809

For the car claims frequency and severity, the data of last year 1,116 policyholders with 2,230 claims are available. As for the life insurances, the USA-2015 mortality table and the policyholders ages are provided.
The project has three main parts:
1.	Available data study
2.	Estimation methods validation
3.	Probability density function goodness of fit testing
4.	Aggregated models construction for life and non-life insurances
5.	Montecarlo simulation for the cost estimation

With this model, the expected cost will be evaluated together with the VaR<sub>99.5</sub> and TVaR<sub>99.5</sub> to quantify the worst-case scenarios.
See the document [Internal_Model_Project_Actuarial_Statistics.pdf](https://github.com/asiergs/Internal-Model-Project-Actuarial-Statistics/blob/main/Internal_Model_Project_Actuarial_Statistics.pdf) for the complete report about the project.

## How to run the code

To run the code, simply run the project.R script once. The code with the default simulation sizes takes about 8.5 hours to run (i5-7600k). The project_lite.R file is the same code with lower number of simulations so the time is reduced to about 15 minutes.

Note that for the data export, the **package _writexl_ needs to be installed**.
