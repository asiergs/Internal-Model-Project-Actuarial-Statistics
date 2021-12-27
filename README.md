## Internal Model Project. Actuarial statistics

In this Project, an insurance company internal model is constructed considering life and non-life insurances. The objective of the model is to properly evaluate the premium risk (the risk of having more claims than expected) so the economic capital can be estimated.
The insurance company to be modelled is composed by the following number of policies:
-	Car insurances: 25,234
-	Life insurances: 20,809

For the car claims frequency and severity, the data of last year 1.116 policyholders with 2.230 claims are available. As for the life insurances, the USA-2015 mortality table and the policyholders ages are provided.
The project has three main parts:
1.	Available data analysis
2.	Estimation methods validation
3.	Probability density function goodness of fit testing
4.	Aggregated models construction for life and non-life
5.	Montecarlo simulation

With this model, the expected cost will be evaluated together with the VaR (99.5%) and TVaR to quantify the worst-case scenarios.

## How to run the code

To run the code, simply run the project.R script once. The code with the default simulation sizes takes about 10 hours to run. The project_lite.R file is the same code with lower simulation sizes so the time is drastically reduced to about 30 minutes.

Note that for the data export, the package _writexl_ needs to be installed.
